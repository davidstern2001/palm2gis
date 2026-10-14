"""Output in a CRS other than PALM's own (crs.output, e.g. EPSG:4326).

A PALM domain is a regular grid in its PROJECTED CRS (crs.palm). An ArcGIS
scene in another CRS - a global scene in WGS 84, typically - needs layers
whose files are in that CRS, not merely displayable in it. Two kinds of
product need two different treatments:

  * Vector layers (surfaces) are transformed
    EXACTLY: every vertex goes through pyproj, z is left untouched.

  * Gridded netCDF (voxel, wind, canopy) cannot be: the PALM grid is not
    regular in the output CRS (at Prague, UTM 33N is rotated about 0.45
    degrees against the meridians and a cell is not square in degrees).
    ArcGIS accepts a voxel source only with 1-D dimension coordinates and
    NO auxiliary coordinate variables, so the data is RESAMPLED onto a
    regular grid in the output CRS, whose spacing matches the native cell
    at the centre of the cut (advanced.netcdf.resampling):

      conservative (default) - PALM is a finite-volume model: a value is
        the AVERAGE over its cell. Each output cell gets the average of
        the PALM cells it overlaps, weighted by the exact overlap area
        (polygon clipping below), over air cells only. Area means and
        totals (leaf area) are preserved; a cell less than half air is
        NoData, so walls stay within half a cell.
      nearest - each output cell copies the PALM cell under its centre:
        values exact, positions good to half a cell.

Horizontal wind components are rotated from PALM's grid axes into the
output CRS's local east/north (true east/north for a geographic CRS), so
u, v and dir_met keep their meaning. w, speeds and scalars are invariant.

With no crs.output (or one equal to crs.palm) everything here is an
identity and the outputs are byte-for-byte what they were before.
"""

import hashlib
import math

import numpy as np
import shapely
from pyproj import CRS, Transformer

from ..log import debug, info

FILL = -9999.0


# ------------------------------
# 1. REGULAR OUTPUT GRID FOR ONE COLUMN STRIDE
# ------------------------------
class OutputGrid:
    """Horizontal geometry of a gridded output.

    x, y   1-D output axes (cell centres, ascending), in the output CRS.
    jj, ii (ny_out, nx_out) indices into the CUT PALM grid of the column
           each output cell takes its value from.
    valid  (ny_out, nx_out) False where the output cell lies outside the
           cut - those cells are NoData.
    """

    def __init__(self, reproj, x, y, jj, ii, valid, stride,
                 method="nearest", weights=None):
        self._r = reproj
        self.x, self.y = x, y
        self.jj, self.ii, self.valid = jj, ii, valid
        self.stride = stride
        # "conservative" only when resampling; the identity grid and the
        # nearest method both work from jj/ii alone.
        self.method = method
        self._w = weights            # (src, w, area) for conservative
        self._cos = self._sin = None
        self._M = None               # scipy CSR of the weights, if any
        self._den_key = self._den = None

    @property
    def conservative(self):
        return self.method == "conservative"

    @property
    def dims(self):
        """Dimension names (y, x) for the variables."""
        return ("lat", "lon") if self._r.geographic else ("y", "x")

    @property
    def shape(self):
        return (self.y.size, self.x.size)

    def take(self, arr, fill=FILL, valid=None):
        """(..., ny_cut, nx_cut) -> (..., ny_out, nx_out).

        *valid* marks the source cells that hold data (air); it matters
        only to the conservative method, which averages over them alone.
        Without it, finite values other than *fill* count as valid.
        Identity at stride 1 returns the array itself, so the default
        (native CRS) path costs no copy of a large 3D record.
        """
        if self._r.identity and self.stride == 1:
            return arr
        if self.conservative:
            return self._remap(arr, fill, valid, mask_only=False)
        out = np.asarray(arr)[..., self.jj, self.ii]
        if not self.valid.all():
            out = np.array(out, copy=True)
            out[..., ~self.valid] = fill
        return out

    def mask(self, masked):
        """NoData mask (True = masked) on the output grid. Conservative: an
        output cell is masked when less than half of it is valid PALM air
        - walls then stay within half a cell of where they are."""
        if self._r.identity and self.stride == 1:
            return masked
        if self.conservative:
            return self._remap(None, None, ~np.asarray(masked, dtype=bool),
                               mask_only=True)
        return self.take(masked, fill=True)

    # Largest gather done at once, in elements (~200 MB of float64).
    _CHUNK = 25_000_000

    def _matrix(self):
        """The weights as a scipy CSR matrix (n_out x n_cut), or None
        without scipy. A sparse product is ~5x faster than the numpy
        gathers on a 1000 x 1000 domain; both give the same sums."""
        if self._M is None:
            try:
                import scipy.sparse as sp
            except ImportError:
                self._M = False
                debug("[proj] scipy not installed - conservative "
                      "resampling uses the slower numpy path.")
                return None
            src, w, _ = self._w
            rows = np.repeat(np.arange(src.shape[0]), src.shape[1])
            keep = w.ravel() > 0.0
            self._M = sp.csr_matrix(
                (w.ravel()[keep], (rows[keep], src.ravel()[keep])),
                shape=(src.shape[0], self._r.ny_cut * self._r.nx_cut))
        return self._M if self._M is not False else None

    def _denominator(self, mat, valid):
        """sum(w * valid) per output cell. The voxel stage passes the same
        air mask for every variable and record, so the last one is kept
        (float32) and reused when the mask has not changed."""
        key = (valid.shape, hashlib.blake2b(
            np.packbits(valid).tobytes(), digest_size=16).digest())
        if key != self._den_key:
            den = np.empty((valid.shape[0], mat.shape[0]), dtype="f4")
            step = max(1, self._CHUNK // max(1, valid.shape[1]))
            for a0 in range(0, valid.shape[0], step):
                a1 = min(valid.shape[0], a0 + step)
                v = np.ascontiguousarray(valid[a0:a1].T, dtype="f8")
                den[a0:a1] = (mat @ v).T
            self._den_key, self._den = key, den
        return self._den

    def _remap(self, arr, fill, valid, mask_only):
        """Area-weighted average over valid source cells.

        out = sum(w * v * a) / sum(w * v), w = overlap area, v = valid,
        a = value; NoData where sum(w * v) is below half the output cell's
        area. First-order conservative: each PALM cell is taken as
        constant, which is how a finite-volume model defines it, so area
        means and totals are kept exactly. Linear in `a`, so rotating
        remapped u/v afterwards is the same as remapping rotated
        components.
        """
        src, w, area = self._w
        ny_c = self._r.ny_cut
        nx_c = self._r.nx_cut
        ref = valid if arr is None else arr
        lead = np.shape(ref)[:-2]
        n_lead = int(np.prod(lead)) if lead else 1
        if valid is None:
            a = np.asarray(arr, dtype="f8")
            valid = np.isfinite(a) & (a != fill)
        valid = np.broadcast_to(np.asarray(valid, dtype=bool),
                                lead + (ny_c, nx_c)).reshape(n_lead, -1)
        vals = (None if mask_only else
                np.asarray(arr, dtype="f8").reshape(n_lead, -1))
        n_out = src.shape[0]
        out = np.empty((n_lead, n_out),
                       dtype=bool if mask_only else "f8")
        half = 0.5 * area
        mat = self._matrix()
        if mat is not None:
            den = self._denominator(mat, valid)
            enough = den >= half[None, :]
            if mask_only:
                out[:] = ~enough
            else:
                step = max(1, self._CHUNK // max(1, valid.shape[1]))
                for a0 in range(0, n_lead, step):
                    a1 = min(n_lead, a0 + step)
                    a = np.where(valid[a0:a1], vals[a0:a1], 0.0)
                    num = (mat @ np.ascontiguousarray(a.T)).T
                    with np.errstate(invalid="ignore", divide="ignore"):
                        out[a0:a1] = np.where(enough[a0:a1],
                                              num / den[a0:a1], fill)
            return out.reshape(lead + self.shape)

        # numpy fallback: the same sums as dense (n_out, K) gathers.
        w32 = w.astype("f4")
        step = max(1, self._CHUNK // max(1, src.size))
        for a0 in range(0, n_lead, step):
            a1 = min(n_lead, a0 + step)
            vf = valid[a0:a1].astype("f4")
            den = np.einsum("lnk,nk->ln", vf[:, src], w32, dtype="f8")
            enough = den >= half
            if mask_only:
                out[a0:a1] = ~enough
                continue
            a = np.where(valid[a0:a1], vals[a0:a1], 0.0).astype("f4")
            num = np.einsum("lnk,nk->ln", a[:, src], w32, dtype="f8")
            with np.errstate(invalid="ignore", divide="ignore"):
                out[a0:a1] = np.where(enough, num / den, fill)
        return out.reshape(lead + self.shape)

    def rotate(self, u, v):
        """Grid-axis components -> output-CRS east/north components."""
        if self._r.identity:
            return u, v
        if self._cos is None:
            c, s = self._r.rotation()
            self._cos, self._sin = c[self.jj, self.ii], s[self.jj, self.ii]
        return (u * self._cos - v * self._sin,
                u * self._sin + v * self._cos)

    def write_axes(self, fout):
        """Create the two horizontal dimensions and their coordinates."""
        yd, xd = self.dims
        fout.createDimension(xd, self.x.size)
        fout.createDimension(yd, self.y.size)
        xv = fout.createVariable(xd, "f8", (xd,))
        xv[:] = self.x
        yv = fout.createVariable(yd, "f8", (yd,))
        yv[:] = self.y
        if self._r.geographic:
            xv.standard_name, xv.long_name = "longitude", "longitude"
            xv.units = "degrees_east"
            yv.standard_name, yv.long_name = "latitude", "latitude"
            yv.units = "degrees_north"
        else:
            xv.standard_name = "projection_x_coordinate"
            yv.standard_name = "projection_y_coordinate"
            xv.units = yv.units = self._r.linear_units
        xv.axis, yv.axis = "X", "Y"
        if not self._r.identity:
            how = ("conservatively (area-weighted average over air cells)"
                   if self.conservative else "by nearest neighbour")
            for v in (xv, yv):
                v.comment = (
                    "Resampled {} from the PALM grid in {} ({:g} m cells); "
                    "spacing matches {} native cell(s) at the centre of "
                    "the cut.".format(how, self._r.src_name,
                                      self._r.native_dx, self.stride))
        return xv, yv


# ------------------------------
# 2. THE MAPPING FOR ONE RUN
# ------------------------------
class Reprojection:
    """PALM CRS -> output CRS for the cut domain of one run."""

    def __init__(self, grid, crs_out=None, method="conservative"):
        if method not in ("conservative", "nearest"):
            raise ValueError(f"resampling must be 'conservative' or "
                             f"'nearest', got {method!r}")
        self.method = method
        self.ny_cut = int(grid.nymax - grid.ny0)
        self.nx_cut = int(grid.nxmax - grid.nx0)
        self.src = grid.crs
        self.crs = CRS.from_user_input(crs_out) if crs_out else grid.crs
        self.identity = crs_out is None or self.crs == self.src
        self.geographic = bool(self.crs.is_geographic) and not self.identity
        # A CRS given as a PROJ string has no name ("unknown").
        self.src_name = (self.src.name if self.src.name != "unknown"
                         else str(grid.crs_raw).strip())
        self.native_dx = float(grid.dx)
        self._g = grid
        self._grids = {}
        self._rot = None
        if not self.identity:
            self._fwd = Transformer.from_crs(self.src, self.crs,
                                             always_xy=True)
            self._inv = Transformer.from_crs(self.crs, self.src,
                                             always_xy=True)

    @property
    def linear_units(self):
        try:
            name = self.crs.axis_info[0].unit_name
        except (IndexError, AttributeError):
            return "m"
        return "m" if name in ("metre", "meter", "m") else name

    # -- points and geometries -------------------------------------------
    def points(self, x, y):
        """Absolute PALM-CRS coordinates -> output CRS (arrays)."""
        if self.identity:
            return np.asarray(x), np.asarray(y)
        return self._fwd.transform(np.asarray(x, dtype="f8"),
                                   np.asarray(y, dtype="f8"))

    def geoms(self, geoms):
        """Transform shapely geometries built in PALM's CRS; z unchanged.

        2-D and 3-D geometries are transformed separately: asking shapely
        to include z for a 2-D point would give it a NaN z.
        """
        if self.identity:
            return geoms
        arr = np.asarray(geoms, dtype=object)
        out = np.empty_like(arr)

        def fwd(coords):
            x, y = self._fwd.transform(coords[:, 0], coords[:, 1])
            res = np.array(coords, dtype="f8", copy=True)
            res[:, 0], res[:, 1] = x, y
            return res

        has_z = shapely.has_z(arr)
        if has_z.any():
            out[has_z] = shapely.transform(arr[has_z], fwd, include_z=True)
        if (~has_z).any():
            out[~has_z] = shapely.transform(arr[~has_z], fwd,
                                            include_z=False)
        return list(out) if isinstance(geoms, list) else out

    # -- vector rotation --------------------------------------------------
    def rotation(self):
        """(cos a, sin a) on the cut grid: a is the angle of PALM's x axis
        against the output CRS's local east, positive anticlockwise.

        Found numerically - step 1 m along PALM's x axis and see where it
        lands - so it holds for any pair of CRSs, datum shifts included.
        For a geographic output CRS the step's direction is its geodesic
        azimuth on the CRS's own ellipsoid (a flat cos(latitude) scaling
        is off by ~1e-3 deg at 50 N). Both CRSs are conformal, so PALM's
        y axis is the x axis turned by 90 degrees.
        """
        if self._rot is None:
            g = self._g
            xx, yy = np.meshgrid(g.abs_x(g.x_stat), g.abs_y(g.y_stat))
            x0, y0 = self._fwd.transform(xx, yy)
            x1, y1 = self._fwd.transform(xx + 1.0, yy)
            if self.geographic:
                az, _, _ = self.crs.get_geod().inv(x0, y0, x1, y1)
                a = np.radians(90.0 - np.asarray(az))
            else:
                a = np.arctan2(y1 - y0, x1 - x0)
            self._rot = (np.cos(a), np.sin(a))
            info("[proj] wind components rotated into {} east/north: "
                 "{:+.3f} .. {:+.3f} deg", self.crs.name,
                 float(np.degrees(a.min())), float(np.degrees(a.max())))
        return self._rot

    def rotate_at(self, u, v, j, i):
        """Rotate components sampled at cut-grid columns (j, i)."""
        if self.identity:
            return u, v
        c, s = self.rotation()
        c, s = c[j, i], s[j, i]
        return u * c - v * s, u * s + v * c

    # -- gridded outputs --------------------------------------------------
    def grid(self, stride=1):
        """OutputGrid for columns thinned by *stride* (cached)."""
        stride = max(1, int(stride))
        if stride not in self._grids:
            self._grids[stride] = (self._identity_grid(stride)
                                   if self.identity
                                   else self._resampled_grid(stride))
        return self._grids[stride]

    def _identity_grid(self, stride):
        g = self._g
        jj = np.arange(0, g.nymax - g.ny0, stride)
        ii = np.arange(0, g.nxmax - g.nx0, stride)
        jg, ig = np.meshgrid(jj, ii, indexing="ij")
        return OutputGrid(self, g.abs_x(g.x_stat)[ii], g.abs_y(g.y_stat)[jj],
                          jg, ig, np.ones(jg.shape, dtype=bool), stride)

    def _resampled_grid(self, stride):
        g = self._g
        xc, yc = g.abs_x(g.x_stat), g.abs_y(g.y_stat)
        nx, ny = xc.size, yc.size
        dx, dy = g.dx, g.dy

        # Footprint of the cut: its outer cell edges, one point per cell,
        # so a curved edge in the output CRS is followed.
        xe = np.append(xc - dx / 2.0, xc[-1] + dx / 2.0)
        ye = np.append(yc - dy / 2.0, yc[-1] + dy / 2.0)
        bx = np.concatenate([xe, xe, np.full(ye.size, xe[0]),
                             np.full(ye.size, xe[-1])])
        by = np.concatenate([np.full(xe.size, ye[0]),
                             np.full(xe.size, ye[-1]), ye, ye])
        bX, bY = self._fwd.transform(bx, by)

        # Output spacing: one native cell, measured at the centre.
        cx, cy = float(xc[nx // 2]), float(yc[ny // 2])
        (xa, xb), _ = self._fwd.transform([cx - dx / 2, cx + dx / 2],
                                          [cy, cy])
        _, (ya, yb) = self._fwd.transform([cx, cx],
                                          [cy - dy / 2, cy + dy / 2])
        dX, dY = abs(xb - xa) * stride, abs(yb - ya) * stride

        x0, x1 = float(np.min(bX)), float(np.max(bX))
        y0, y1 = float(np.min(bY)), float(np.max(bY))
        nX = max(1, int(math.ceil((x1 - x0) / dX - 1e-9)))
        nY = max(1, int(math.ceil((y1 - y0) / dY - 1e-9)))
        X = x0 + (np.arange(nX) + 0.5) * dX
        Y = y0 + (np.arange(nY) + 0.5) * dY

        XX, YY = np.meshgrid(X, Y)
        sx, sy = self._inv.transform(XX, YY)
        fi = (sx - xc[0]) / dx
        fj = (sy - yc[0]) / dy
        valid = (fi >= -0.5) & (fi < nx - 0.5) & (fj >= -0.5) & (fj < ny - 0.5)
        ii = np.clip(np.rint(fi), 0, nx - 1).astype(int)
        jj = np.clip(np.rint(fj), 0, ny - 1).astype(int)

        weights = None
        if self.method == "conservative":
            # Output cell corners in PALM index space, where source cell
            # (j, i) is exactly the unit square [i, i+1) x [j, j+1).
            XE = x0 + np.arange(nX + 1) * dX
            YE = y0 + np.arange(nY + 1) * dY
            cx_, cy_ = self._inv.transform(*np.meshgrid(XE, YE))
            ci = (np.asarray(cx_) - (xc[0] - dx / 2.0)) / dx
            cj = (np.asarray(cy_) - (yc[0] - dy / 2.0)) / dy
            weights = _overlap_weights(ci, cj, nx, ny)
            covered = weights[1].sum(axis=1) >= 0.5 * weights[2]
            valid = valid | covered.reshape(valid.shape)

        info("[proj] output grid {} x {}, {:.0f} % outside the domain",
             nX, nY, 100.0 * (1.0 - valid.mean()))
        debug("[proj] {} grid: {} x {} cells of {:.3g} x {:.3g} {}, "
              "{} from {} x {} PALM columns{} "
              "({:.0f} % outside the domain, NoData)",
              self.crs.name, nX, nY, dX, dY,
              "deg" if self.geographic else self.linear_units,
              "conservative (area-weighted)" if weights is not None
              else "nearest neighbour",
              nx, ny, " (every {})".format(stride) if stride > 1 else "",
              100.0 * (1.0 - valid.mean()))
        return OutputGrid(self, X, Y, jj, ii, valid, stride,
                          method=self.method, weights=weights)


# ------------------------------
# 3. EXACT OVERLAP AREAS (CONSERVATIVE REMAPPING)
# ------------------------------
def _overlap_weights(ci, cj, nx, ny, chunk=500_000):
    """Overlap of every output cell with the PALM cells beneath it.

    ci, cj: (nY+1, nX+1) output cell corners in PALM index space. An
    output cell is the quadrilateral of its four corners - exact up to
    the curvature of its edges, which is ~1e-10 of a cell for a 1 m cell.

    Returns (src, w, area): src (n_out, K) flat PALM cell indices, w the
    overlap areas (0 for a candidate outside the domain or not touched),
    area (n_out,) the output cell's own area - all in PALM cell units.
    """
    # CCW corners of each output cell: (i, j) (i+1, j) (i+1, j+1) (i, j+1)
    qx = np.stack([ci[:-1, :-1], ci[:-1, 1:], ci[1:, 1:], ci[1:, :-1]],
                  axis=-1).reshape(-1, 4)
    qy = np.stack([cj[:-1, :-1], cj[:-1, 1:], cj[1:, 1:], cj[1:, :-1]],
                  axis=-1).reshape(-1, 4)
    area = np.abs(_shoelace(qx, qy, np.full(qx.shape[0], 4)))

    # Candidate PALM cells: the window of cells the quad's bounding box
    # touches, the same size for every output cell so it vectorises.
    i0 = np.floor(qx.min(axis=1)).astype(int)
    j0 = np.floor(qy.min(axis=1)).astype(int)
    wx = int((np.floor(qx.max(axis=1)).astype(int) - i0).max()) + 1
    wy = int((np.floor(qy.max(axis=1)).astype(int) - j0).max()) + 1
    di, dj = np.meshgrid(np.arange(wx), np.arange(wy))
    ci_ = i0[:, None] + di.ravel()[None, :]             # (n_out, K)
    cj_ = j0[:, None] + dj.ravel()[None, :]
    inside = (ci_ >= 0) & (ci_ < nx) & (cj_ >= 0) & (cj_ < ny)

    # Only candidates inside the quad's own bounding box and the domain
    # can overlap it; the rest of the common window is skipped unclipped.
    near = (inside & (ci_ <= np.floor(qx.max(axis=1))[:, None])
            & (cj_ <= np.floor(qy.max(axis=1))[:, None]))
    n_out, k = ci_.shape
    w = np.zeros((n_out, k))
    rows, cols = np.nonzero(near)
    # Clip in coordinates relative to each candidate cell, so the square
    # is always [0, 1]^2 and precision does not depend on the index.
    for a in range(0, rows.size, chunk):
        b = min(rows.size, a + chunk)
        r, c = rows[a:b], cols[a:b]
        px = qx[r] - ci_[r, c][:, None]
        py = qy[r] - cj_[r, c][:, None]
        w[r, c] = _clip_unit_square_area(px, py)
    w[w < 1e-12] = 0.0
    src = (np.clip(cj_, 0, ny - 1) * nx + np.clip(ci_, 0, nx - 1))

    # Pack the overlaps to the front of each row and drop the columns no
    # row uses: usually 4 of the 9 candidates remain, and the remap cost
    # is proportional to that width.
    order = np.argsort(w <= 0.0, axis=1, kind="stable")
    w = np.take_along_axis(w, order, axis=1)
    src = np.take_along_axis(src, order, axis=1)
    kmax = max(1, int((w > 0.0).sum(axis=1).max()))
    return (np.ascontiguousarray(src[:, :kmax]),
            np.ascontiguousarray(w[:, :kmax]), area)


def _shoelace(x, y, n):
    """Signed area of polygons with n[r] leading vertices per row."""
    k = np.arange(x.shape[1])[None, :]
    nxt = (k + 1) % np.maximum(n, 1)[:, None]
    xn = np.take_along_axis(x, nxt, axis=1)
    yn = np.take_along_axis(y, nxt, axis=1)
    term = np.where(k < n[:, None], x * yn - xn * y, 0.0)
    return 0.5 * term.sum(axis=1)


def _clip_unit_square_area(px, py):
    """Area of convex polygons (rows of px, py, CCW) inside [0, 1]^2.

    Sutherland-Hodgman against the four sides, vectorised over rows: the
    polygon lives in a fixed 8-vertex buffer (a quad clipped by a square
    has at most 8 vertices) with a per-row vertex count.
    """
    m = px.shape[0]
    vmax = 8
    x = np.zeros((m, vmax))
    y = np.zeros((m, vmax))
    x[:, :px.shape[1]], y[:, :py.shape[1]] = px, py
    n = np.full(m, px.shape[1])
    for axis, bound, keep_ge in ((0, 0.0, True), (0, 1.0, False),
                                 (1, 0.0, True), (1, 1.0, False)):
        x, y, n = _clip_halfplane(x, y, n, axis, bound, keep_ge, vmax)
    return np.abs(_shoelace(x, y, n))


def _clip_halfplane(x, y, n, axis, bound, keep_ge, vmax):
    k = np.arange(vmax)[None, :]
    active = k < n[:, None]
    nxt = (k + 1) % np.maximum(n, 1)[:, None]
    xn = np.take_along_axis(x, nxt, axis=1)
    yn = np.take_along_axis(y, nxt, axis=1)
    c, cn = (x, xn) if axis == 0 else (y, yn)
    if keep_ge:
        cur_in, nxt_in = c >= bound, cn >= bound
    else:
        cur_in, nxt_in = c <= bound, cn <= bound
    cross = active & (cur_in != nxt_in)
    denom = np.where(cn != c, cn - c, 1.0)
    t = np.clip((bound - c) / denom, 0.0, 1.0)
    ix, iy = x + t * (xn - x), y + t * (yn - y)
    # Each edge emits [intersection if it crosses, next vertex if inside].
    ox = np.stack([ix, xn], axis=-1).reshape(x.shape[0], -1)
    oy = np.stack([iy, yn], axis=-1).reshape(x.shape[0], -1)
    ok = np.stack([cross, active & nxt_in], axis=-1).reshape(x.shape[0], -1)
    # Compact the kept vertices to the front, in order (at most vmax).
    pos = np.cumsum(ok, axis=1) - 1
    keep = ok & (pos < vmax)
    r = np.nonzero(keep)[0]
    xo = np.zeros((x.shape[0], vmax))
    yo = np.zeros((x.shape[0], vmax))
    xo[r, pos[keep]] = ox[keep]
    yo[r, pos[keep]] = oy[keep]
    return xo, yo, np.minimum(ok.sum(axis=1), vmax)
