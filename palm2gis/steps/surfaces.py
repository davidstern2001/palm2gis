"""PALM surface output -> GIS layers.

GIS counterpart of paraview_vis/utils/surf_gen.py. Instead of building VTU
quads for ParaView, the unstructured surface elements (xs, ys, zs with
azimuth/zenith) are split by orientation:

  * UPWARD-facing elements (zenith == 0: ground + roofs). For each grid
    column the TOPMOST element is taken - i.e. what a nadir thermal camera
    would see - and rasterised into one multiband GeoTIFF per variable
    (band = timestep, band description = timestamp). Loads directly in
    ArcGIS Pro / QGIS with correct georeferencing.

  * WALL elements (zenith == 90). Each element becomes a vertical quad
    (PolygonZ) in absolute coordinates, written to a GeoPackage in long
    format (elem_id, azimuth, z_centre, time, one column per variable), so
    the ArcGIS time slider / QGIS temporal controller can animate it and
    the quads render correctly in 3D scenes.

Elements with other zenith values (sloped roofs, downward-facing) are
counted and reported but not converted.
"""

import numpy as np
from netCDF4 import Dataset

import rasterio
from rasterio.transform import from_origin
import geopandas as gpd
import pyogrio
from shapely.geometry import LineString, MultiPolygon, Polygon

from ..log import debug, get_logger, info, warning
from .georef import check_product, check_raw_input, select_times
from .thermo import convert_temperature, exner_simple, is_celsius_exempt, \
    is_temperature, temperature_scale, theta_to_ta
from .qgis_style import (embed_gpkg_style, raster_qml, vector_qml,
                         write_sidecar)

log = get_logger()

FILL = -9999.0


# ------------------------------
# 1. LOAD AND CLASSIFY ELEMENTS
# ------------------------------
def _load_geometry(nc):
    xs = np.asarray(nc.variables["xs"][:])
    ys = np.asarray(nc.variables["ys"][:])
    zs = np.asarray(nc.variables["zs"][:])
    zen = np.ma.filled(nc.variables["zenith"][:], FILL)
    az = np.ma.filled(nc.variables["azimuth"][:], FILL)
    return xs, ys, zs, zen, az


def _bounds_mask(cfg, xs, ys):
    bnd = cfg["bounds"]
    m = np.ones(xs.size, dtype=bool)
    m &= xs >= bnd.get("x0", -1.0e30)
    m &= xs <= bnd.get("xmax", 1.0e30)
    m &= ys >= bnd.get("y0", -1.0e30)
    m &= ys <= bnd.get("ymax", 1.0e30)
    return m


# ------------------------------
# 2. OUTPUT VARIABLE LIST
# ------------------------------
# PALM naming: `t_surf` is the actual surface (skin) temperature, while
# `theta_surface` / `theta_surf` is the surface POTENTIAL temperature.
# Both spellings of the potential form occur depending on PALM version and
# output block, so both are accepted.
SURF_ACTUAL = "t_surf"
SURF_POTENTIAL = ("theta_surface", "theta_surf")


def output_vars(sf, available=None, verbose=False):
    """(output name, source variable, derive?) for every surface product.

    With surfaces.derive_t_surf the actual SURFACE temperature `t_surf` is
    ENSURED: if the file already holds it, PALM's own value is used as-is;
    only when it is absent is it derived from the surface potential
    temperature. The same list drives the raster, the mesh NetCDF and the
    element layers, so a field is never present in one product and missing
    from another.

    `available` is the set of variable names in the source file. Listed
    variables that are missing but present under the other accepted
    spelling are silently substituted.
    """
    have = set(available or ())
    out, names = [], []

    for var in sf["variables"]:
        src = var
        if have and var not in have:
            alias = next((a for a in SURF_POTENTIAL
                          if a != var and a in have), None)
            if var in SURF_POTENTIAL and alias:
                if verbose:
                    warning("[surfaces] '{}' not in the file; using '{}'.",
                            var, alias)
                src = alias
        out.append((var, src, False))
        names.append(var)

    if not sf["derive_t_surf"]:
        return out

    if SURF_ACTUAL in names:
        return out                       # user asked for it explicitly
    if have and SURF_ACTUAL in have:
        if verbose:
            info("[surfaces] the file holds '{}' - using PALM's own "
                 "surface temperature instead of deriving it.",
                 SURF_ACTUAL)
        out.append((SURF_ACTUAL, SURF_ACTUAL, False))
        return out

    src = next((v for _n, v, _d in out if v in SURF_POTENTIAL), None)
    if src is None:
        if verbose:
            warning("[surfaces] derive_t_surf is on but neither '{}' nor a "
                    "surface potential temperature ({}) is available.",
                    SURF_ACTUAL, ", ".join(SURF_POTENTIAL))
        return out
    out.append((SURF_ACTUAL, src, True))
    return out


def _derive_values(cfg, vals, z_elem, derive, name=None, src_units=None):
    """theta -> T conversion (when *derive*) and the output scale
    conversion for temperature variables (cfg["celsius"]).

    The Exner-based derivation runs in kelvin; the scale conversion is
    applied last, only to valid entries - FILL sentinels and NaN stay
    untouched.

    `src_units` is the source variable's `units` attribute and decides
    the INPUT scale. Assuming kelvin here was wrong for PALM's own
    Celsius variables; see thermo.temperature_scale.

    NOTE `z_elem` is height above the DOMAIN BASE, not above sea level.
    exner_simple pairs it with cfg["surface_pressure"], defined as the
    pressure at the model's z = 0. Passing absolute (origin_z-shifted)
    heights here would put a second, larger error on every derived
    temperature.
    """
    out = vals
    if derive:
        good = np.isfinite(vals) & (vals != FILL)
        if good.any():
            pi = exner_simple(z_elem,
                              theta_ref=float(np.nanmean(vals[good])),
                              p_surface=cfg["surface_pressure"])
            out = np.array(vals, dtype="f8", copy=True)
            out[good] = theta_to_ta(vals[good], pi[good])
    if name is None or not is_temperature(name) or is_celsius_exempt(name):
        return out
    # A derived value came out of theta_to_ta, so it is kelvin whatever
    # the source variable was labelled.
    source = "K" if derive else temperature_scale(
        src_units, name, dict(cfg["temperature_units"] or {}))
    return convert_temperature(out, source, "C" if cfg["celsius"] else "K",
                               fill=FILL)


# ------------------------------
# 3. TOP-VIEW NETCDF (QGIS MESH LAYER)
# ------------------------------
def _topview_netcdf(cfg, grid, nc, its, times, top_s, valid, zs,
                    i0, i1, j0, j1, path):
    """Same top-view fields as the GeoTIFF, but as a CF NetCDF.

    QGIS opens this as a MESH layer, and MDAL reads the CF time axis
    natively - the Temporal Controller then works with no style sidecar
    and no band-to-time mapping to get wrong. All configured variables go
    into one file, appearing as separate dataset groups on the layer.
    """
    from netCDF4 import Dataset as NCDataset
    from .voxel import write_crs_var

    sf = cfg["surfaces"]
    nj, ni = top_s.shape
    h, w = j1 - j0, i1 - i0

    with NCDataset(path, "w", format="NETCDF4") as fout:
        fout.Conventions = "CF-1.8"
        fout.title = "PALM surface output (top view) for QGIS mesh layers"

        fout.createDimension("x", w)
        fout.createDimension("y", h)
        fout.createDimension("time", len(its))

        xv = fout.createVariable("x", "f8", ("x",))
        xv[:] = (grid.origin_x
                 + (np.arange(i0, i1) + 0.5) * grid.dx)
        xv.standard_name = "projection_x_coordinate"
        xv.units = "m"
        xv.axis = "X"

        yv = fout.createVariable("y", "f8", ("y",))
        yv[:] = (grid.origin_y
                 + (np.arange(j0, j1) + 0.5) * grid.dy)
        yv.standard_name = "projection_y_coordinate"
        yv.units = "m"
        yv.axis = "Y"

        tv = fout.createVariable("time", "f8", ("time",))
        t_sel = np.asarray(times)[its]
        if np.all(np.isfinite(t_sel)) and np.all(np.abs(t_sel) < 1.0e12):
            tv[:] = np.round(t_sel)
            tv.standard_name = "time"
            tv.units = "seconds since {}".format(
                grid.origin_time.strftime("%Y-%m-%d %H:%M:%S"))
            tv.calendar = "standard"
        else:
            tv[:] = np.arange(len(its), dtype="f8")
            tv.long_name = "record index (no valid time values in input)"
            tv.units = "1"
        tv.axis = "T"

        gm = write_crs_var(fout, grid.crs)

        for name, var, derive in output_vars(sf, nc.variables):
            ov = fout.createVariable(name, "f4", ("time", "y", "x"),
                                     fill_value=FILL, zlib=True,
                                     complevel=4)
            ov.grid_mapping = gm
            # See the note in voxel.py: no auxiliary coordinates.
            if derive:
                ov.units = "degC" if cfg["celsius"] else "K"
                ov.long_name = ("surface air temperature derived from {} "
                                "via the Exner function".format(var))
            else:
                src = nc.variables.get(var)
                if src is not None:
                    for att in ("units", "long_name"):
                        if hasattr(src, att):
                            ov.setncattr(att, getattr(src, att))
                if is_temperature(name) and not is_celsius_exempt(name):
                    ov.units = "degC" if cfg["celsius"] else "K"
            # RELATIVE height: this feeds the Exner function, which pairs
            # it with surface_pressure at the model's z = 0. Not abs_z.
            z_top = zs[top_s[valid]]
            for b, it in enumerate(its):
                data = np.ma.filled(nc.variables[var][it, :], FILL)
                band = np.full((nj, ni), FILL, dtype=np.float32)
                band[valid] = _derive_values(
                    cfg, data[top_s[valid]], z_top, derive, name=name,
                    src_units=getattr(nc.variables[var], "units", None))
                ov[b, ...] = band[j0:j1, i0:i1]

    info("[surfaces] wrote {}: {} vars x {} times (mesh layer)",
         cfg.disp(path), len(sf["variables"]), len(its))


# ------------------------------
# 4. TOP-VIEW RASTER (ROOFS + GROUND)
# ------------------------------
def _topview_raster(cfg, grid, nc, its, times, xs, ys, zs, zen, sel,
                    targets, nc_path=None):
    sf = cfg["surfaces"]
    dx, dy = grid.dx, grid.dy

    up = sel & (np.abs(zen) < 1.0e-3)
    idx_up = np.where(up)[0]
    if idx_up.size == 0:
        warning("[surfaces] no upward-facing elements in bounds - "
                "top-view raster skipped.")
        return

    # map each upward element to a grid cell, keep the topmost per cell.
    #
    # This is an argmax-by-group, which numpy does natively. The Python
    # loop it replaces ran once per upward surface element, which on a real
    # domain is millions of iterations and dominated the stage. Sorting by
    # (cell, height) and taking the LAST entry of each cell group leaves the
    # highest element per cell, exactly as the running comparison did.
    i = np.round((xs[idx_up] - dx / 2.0) / dx).astype(int)
    j = np.round((ys[idx_up] - dy / 2.0) / dy).astype(int)
    ni, nj = int(i.max()) + 1, int(j.max()) + 1
    top_z = np.full((nj, ni), -np.inf)
    top_s = np.full((nj, ni), -1, dtype=np.int64)

    cell = j.astype(np.int64) * ni + i.astype(np.int64)
    z_up = zs[idx_up]
    order = np.lexsort((z_up, cell))          # by cell, then by height
    cell_s, z_s, idx_s = cell[order], z_up[order], idx_up[order]
    last = np.empty(cell_s.size, dtype=bool)  # last of each cell run
    last[-1] = True
    last[:-1] = cell_s[1:] != cell_s[:-1]
    sel_cell, sel_z, sel_idx = cell_s[last], z_s[last], idx_s[last]
    top_z.reshape(-1)[sel_cell] = sel_z
    top_s.reshape(-1)[sel_cell] = sel_idx

    valid = top_s >= 0
    debug("[surfaces] {} upward elements -> {} raster cells",
          idx_up.size, int(valid.sum()))

    # georeferencing: raster origin at top-left, so flip y on write
    i0, i1 = int(i.min()), ni
    j0, j1 = int(j.min()), nj
    transform = from_origin(grid.origin_x + i0 * dx,
                            grid.origin_y + j1 * dy, dx, dy)

    if nc_path is not None:
        _topview_netcdf(cfg, grid, nc, its, times, top_s, valid, zs,
                        i0, i1, j0, j1, nc_path)

    dts = grid.times_to_datetime(times)
    z_top = zs[top_s[valid]]
    src_of = {n: (v, d) for n, v, d in output_vars(sf, nc.variables)}
    for name, path in targets.items():
        var, derive = src_of.get(name, (name, False))
        h, w = j1 - j0, i1 - i0
        with rasterio.open(path, "w", driver="GTiff", height=h, width=w,
                           count=len(its), dtype="float32",
                           crs=grid.rio_crs(),
                           transform=transform, nodata=FILL,
                           compress="deflate") as dst:
            for b, it in enumerate(its, start=1):
                data = np.ma.filled(nc.variables[var][it, :], FILL)
                band = np.full((nj, ni), FILL, dtype=np.float32)
                band[valid] = _derive_values(
                    cfg, data[top_s[valid]], z_top, derive, name=name,
                    src_units=getattr(nc.variables[var], "units", None))
                band = band[j0:j1, i0:i1]
                dst.write(band[::-1, :], b)          # north-up
                dst.set_band_description(b, grid.label(dts[it], it))
        if cfg["qgis_styles"]:
            write_sidecar(path, raster_qml([dts[it] for it in its]))
        _ensure_raster_crs(cfg, grid, path)
        info("[surfaces] wrote {}: {} bands", cfg.disp(path), len(its))


# ------------------------------
# 3. SURFACE ELEMENT POLYGONS
# ------------------------------
_AZ_FACING = {0: "N", 90: "E", 180: "S", 270: "W"}
_AZ_TOL = 1.0e-3


def _wall_facing(az):
    """Cardinal direction a wall faces, or None if it is off-axis.

    PALM writes wall azimuths as exact multiples of 90 degrees, but they
    arrive as float32 and pass through arithmetic, so an exact `==` test is
    not safe. This is the single place that decides which axis a wall is
    on: `_element_quad` and `_orientation` used to decide separately, the
    first with `az in (0.0, 180.0)` (exact equality) and the second with
    `round(az, 1)`. An azimuth of 179.99999 was therefore labelled wall_S
    while being given an x-normal quad - a polygon rotated 90 degrees from
    its own label, in the layer used to attribute wall temperatures.
    """
    a = float(az) % 360.0
    for k, v in _AZ_FACING.items():
        if abs(a - k) < _AZ_TOL or abs(a - k - 360.0) < _AZ_TOL:
            return v
    return None


def _element_quad(x, y, z, zen, az, grid):
    """Georeferenced quad (PolygonZ) for one surface element.

    Upward (zenith 0) and downward (zenith 180) faces become horizontal
    dx x dy quads at their own height; wall faces (zenith 90) become
    vertical dz-tall quads spanning the cell in the direction the wall
    runs. Absolute coordinates.
    """
    ox, oy = grid.origin_x, grid.origin_y
    dx, dy, dz = grid.dx, grid.dy, grid.dz
    # z is shifted to absolute here for the same reason x and y are: the
    # quad is output geometry. PALM's zs is height above the domain base,
    # so a quad built from it raw sat origin_z metres below the buildings
    # it describes.
    z = float(grid.abs_z(z))
    if abs(zen - 90.0) > 45.0:            # horizontal face (up or down)
        return Polygon([(ox + x - dx / 2.0, oy + y - dy / 2.0, z),
                        (ox + x + dx / 2.0, oy + y - dy / 2.0, z),
                        (ox + x + dx / 2.0, oy + y + dy / 2.0, z),
                        (ox + x - dx / 2.0, oy + y + dy / 2.0, z)])
    z0, z1 = z - dz / 2.0, z + dz / 2.0
    if _wall_facing(az) in ("N", "S"):    # y-normal wall, spans x
        return Polygon([(ox + x - dx / 2.0, oy + y, z0),
                        (ox + x + dx / 2.0, oy + y, z0),
                        (ox + x + dx / 2.0, oy + y, z1),
                        (ox + x - dx / 2.0, oy + y, z1)])
    return Polygon([(ox + x, oy + y - dy / 2.0, z0),   # x-normal wall
                    (ox + x, oy + y + dy / 2.0, z0),
                    (ox + x, oy + y + dy / 2.0, z1),
                    (ox + x, oy + y - dy / 2.0, z1)])


def _orientation(zen, az):
    """Human-readable face orientation for symbology and filtering."""
    if abs(zen) < 1.0e-3:
        return "up"
    if abs(zen - 180.0) < 1.0e-3:
        return "down"
    if abs(zen - 90.0) < 1.0e-3:
        return "wall_" + (_wall_facing(az) or "?")
    return "other"


def _wall_trace(poly):
    """Horizontal trace of a vertical quad, as a 2D LineString.

    The quad spans one cell horizontally and one dz vertically, so its
    footprint is a single segment; Z is dropped so the line renders in a
    plain 2D map view.
    """
    xy = sorted(set((round(x, 6), round(y, 6))
                    for x, y, *_ in poly.exterior.coords))
    if len(xy) < 2:                       # degenerate, should not happen
        return poly
    return LineString([xy[0], xy[-1]])


def _wall_quad(x, y, z, az, grid):
    """Backwards-compatible helper: a wall element is zenith 90."""
    return _element_quad(x, y, z, 90.0, az, grid)


# ------------------------------
# 3b. MULTIPATCH / FILE GEODATABASE
# ------------------------------
# Shapefile field names are capped at 10 characters by DBF. Truncating
# silently would collapse `theta_surface_derived` and `theta_surface_raw`
# onto the same column, so the map is explicit and collisions are
# reported rather than resolved by luck.
_DBF_NAME_MAX = 10


def _dbf_field_map(names):
    """Map attribute names to DBF-legal (<=10 char) names, uniquely."""
    out, taken = {}, set()
    for n in names:
        base = str(n)[:_DBF_NAME_MAX]
        cand, k = base, 1
        while cand in taken:
            suffix = str(k)
            cand = base[:_DBF_NAME_MAX - len(suffix)] + suffix
            k += 1
        if cand != n:
            warning("[surfaces] shapefile field '{}' -> '{}' (DBF caps "
                    "names at {} characters).", n, cand, _DBF_NAME_MAX)
        out[n] = cand
        taken.add(cand)
    return out


def _element_frames(cfg, grid, nc, its, times, xs, ys, zs, zen, az, idx,
                    multipatch=False):
    """Yield (timestep index, GeoDataFrame) for the selected elements.

    Shared by the GeoPackage, multipatch and File Geodatabase writers so
    the geometry and the attributes cannot drift apart between products.
    With `multipatch` the quads are wrapped in MultiPolygon, which is
    what GDAL's shapefile driver needs to emit shape type 31.
    """
    sf = cfg["surfaces"]
    quads = [_element_quad(xs[s], ys[s], zs[s], zen[s], az[s], grid)
             for s in idx]
    geoms = [MultiPolygon([q]) for q in quads] if multipatch else quads
    orient = [_orientation(zen[s], az[s]) for s in idx]

    i = np.clip(np.round((xs[idx] - grid.dx / 2.0) / grid.dx).astype(int)
                - grid.nx0, 0, grid.zt.shape[1] - 1)
    j = np.clip(np.round((ys[idx] - grid.dy / 2.0) / grid.dy).astype(int)
                - grid.ny0, 0, grid.zt.shape[0] - 1)
    z_agl = zs[idx] - grid.zt[j, i]

    dts_tz = grid.times_to_datetime(times, tz_aware=True)

    for n, it in enumerate(its):
        # int32 explicitly: OpenFileGDB has no 64-bit integer field, so
        # an int64 index silently arrives back as a float.
        rec = {"elem_id": idx.astype("i4"),
               "orient": orient,
               "zenith": zen[idx],
               "z": grid.abs_z(zs[idx]),
               "z_agl": z_agl,
               "rec": np.full(idx.size, int(it), dtype="i4")}
        if multipatch:
            # DBF has no datetime type: GDAL downgrades one to a string
            # with a warning. Writing the string deliberately means the
            # format is chosen rather than inherited, and ArcGIS Pro can
            # time-enable a text field on it.
            stamp = dts_tz[it]
            rec["time"] = [stamp.strftime("%Y-%m-%dT%H:%M:%S")
                           if stamp is not None else ""] * idx.size
            # `azimuth` is omitted on purpose: shapefiles have no numeric
            # NULL, so the -9999 sentinel for horizontal faces could not
            # be cleaned up. It is recoverable from `orient`.
        else:
            rec["time"] = [dts_tz[it]] * idx.size
            rec["azimuth"] = np.where(np.isclose(az[idx], FILL), np.nan,
                                      az[idx])
        for name, var, derive in output_vars(sf, nc.variables):
            vals = np.ma.filled(nc.variables[var][it, :], np.nan)[idx]
            rec[name] = _derive_values(
                cfg, vals, zs[idx], derive, name=name,
                src_units=getattr(nc.variables[var], "units", None))
        yield n, it, gpd.GeoDataFrame(rec, geometry=geoms, crs=grid.crs)


def _elements_multipatch(cfg, grid, nc, its, times, xs, ys, zs, zen, az,
                         idx, path):
    """Write surface elements as an ESRI MultiPatch shapefile.

    Multipatch (shape type 31) is ArcGIS's native 3D geometry type and
    the one a 3D Object scene layer is published from. GDAL reaches it
    only through the shapefile driver - OpenFileGDB writes MultiPolygon
    Z, not multipatch - which brings the DBF limits handled above.
    """
    sf = cfg["surfaces"]
    if idx.size == 0:
        warning("[surfaces] no elements in bounds - {} skipped.",
                cfg.disp(path))
        return

    split_at = int(sf["multipatch_split"] or 0)
    # A .shp and its .dbf are each capped at 2 GB, so a long run has to
    # be split. One file per timestep is the only split that keeps every
    # feature's geometry intact.
    per_time = bool(split_at) and idx.size * len(its) > split_at
    if per_time:
        info("[surfaces] {} elements x {} times exceeds "
             "surfaces.multipatch_split ({}) - writing one shapefile per "
             "timestep.", idx.size, len(its), split_at)

    written, frames = [], []
    for n, it, frame in _element_frames(cfg, grid, nc, its, times, xs, ys,
                                        zs, zen, az, idx, multipatch=True):
        if per_time:
            stamp = frame["time"].iloc[0].replace(":", "").replace("-", "")
            out = path.with_name("{}_{}{}".format(path.stem, stamp,
                                                  path.suffix))
            _write_multipatch(frame, out)
            written.append(out)
        else:
            frames.append(frame)

    if not per_time:
        frame = (frames[0] if len(frames) == 1
                 else gpd.pd.concat(frames, ignore_index=True))
        _write_multipatch(frame, path)
        written.append(path)

    info("[surfaces] wrote {}: {} multipatch features x {} times{}",
         cfg.disp(path), idx.size, len(its),
         " in {} files".format(len(written)) if len(written) > 1 else "")


def _write_multipatch(frame, path):
    """One shapefile, geometry forced to MultiPatch, fields DBF-safe."""
    if path.exists():
        path.unlink()
    cols = [c for c in frame.columns if c != "geometry"]
    frame = frame.rename(columns=_dbf_field_map(cols))
    pyogrio.write_dataframe(frame, path, driver="ESRI Shapefile",
                            SHPT="MULTIPATCH")
    # .prj comes from the frame's CRS; without it ArcGIS drops the layer
    # at the origin, which is a miserable thing to debug.
    if not path.with_suffix(".prj").exists():
        warning("[surfaces] {} has no .prj - the CRS was not written.",
                path.name)


def _elements_fgdb(cfg, grid, nc, its, times, xs, ys, zs, zen, az, idx,
                   path, layer):
    """Write the same elements as PolygonZ into a File Geodatabase.

    The GeoPackage keeps QGIS happy; this is the container ArcGIS opens
    natively. Unlike the multipatch shapefile it preserves real
    datetimes, full-length field names and NULLs, so nothing is lost
    relative to the GeoPackage.
    """
    if idx.size == 0:
        warning("[surfaces] no elements in bounds - {} skipped.",
                cfg.disp(path))
        return
    for n, it, frame in _element_frames(cfg, grid, nc, its, times, xs, ys,
                                        zs, zen, az, idx):
        pyogrio.write_dataframe(frame, path, driver="OpenFileGDB",
                                layer=layer, append=(n > 0))
    info("[surfaces] wrote {}: feature class '{}', {} elements x {} times",
         cfg.disp(path), layer, idx.size, len(its))


def _elements_gpkg(cfg, grid, nc, its, times, xs, ys, zs, zen, az, sel,
                   path, idx, kind):
    """Write one GeoPackage row per surface element per timestep.

    `idx` selects which elements go in (all of them, or just the walls).
    Every face of the static driver's geometry is represented: ground,
    roofs and each of the four wall orientations, each as a PolygonZ in
    absolute coordinates carrying the configured variables.
    """
    sf = cfg["surfaces"]
    if idx.size == 0:
        warning("[surfaces] no {} elements in bounds - {} skipped.",
                kind, cfg.disp(path))
        return

    n_rows = idx.size * len(its)
    if n_rows > 2_000_000:
        warning("[surfaces] {} will hold {} rows - this will be slow and "
                "large. Tighten bounds/n_times or disable it.",
                cfg.disp(path), n_rows)

    geoms = [_element_quad(xs[s], ys[s], zs[s], zen[s], az[s], grid)
             for s in idx]
    orient = [_orientation(zen[s], az[s]) for s in idx]

    # height above the terrain of that column: separates roofs from ground
    # without needing buildings_2d (which some static drivers leave empty)
    i = np.clip(np.round((xs[idx] - grid.dx / 2.0) / grid.dx).astype(int)
                - grid.nx0, 0, grid.zt.shape[1] - 1)
    j = np.clip(np.round((ys[idx] - grid.dy / 2.0) / grid.dy).astype(int)
                - grid.ny0, 0, grid.zt.shape[0] - 1)
    z_agl = zs[idx] - grid.zt[j, i]

    dts = grid.times_to_datetime(times)
    dts_tz = grid.times_to_datetime(times, tz_aware=True)
    layer = sf.get("layer_name") or kind

    # A GeoPackage is written one timestep at a time: the first with
    # mode="w", the rest appended. Building every frame in a list and
    # concatenating meant peak memory scaled with elements x timesteps,
    # so the row warning above fired only after the memory had already
    # been committed. Peak memory is now one timestep regardless of how
    # many are selected.
    is_wall = np.array([o.startswith("wall_") for o in orient])
    wall_idx = np.where(is_wall)[0]
    want_lines = bool(sf["wall_lines"]) and wall_idx.size > 0
    line_layer = layer + "_2d"
    wall_geoms = ([_wall_trace(geoms[k]) for k in wall_idx]
                  if want_lines else None)

    # An existing file must not be appended to, or a rerun doubles it.
    if path.exists():
        path.unlink()

    # The ramp is built over the FIRST output variable - with one variable
    # configured (the usual case) that is simply the variable of interest,
    # and with several it is a defined choice rather than an arbitrary one.
    _out = list(output_vars(sf, nc.variables))
    ramp_attr = _out[0][0] if _out else None
    v_lo = v_hi = None

    n_written = 0
    for n, it in enumerate(its):
        rec = {"elem_id": idx,
               "orient": orient,
               "zenith": zen[idx],
               # NULL, not a sentinel: a magic -9999 in a numeric attribute
               # has no NoData meaning to a GIS and skews every
               # classification and statistic computed on the field.
               "azimuth": np.where(np.isclose(az[idx], FILL), np.nan,
                                   az[idx]),
               "z": grid.abs_z(zs[idx]),
               "z_agl": z_agl,
               "rec": it,
               "time": [dts_tz[it]] * idx.size}
        for name, var, derive in output_vars(sf, nc.variables):
            vals = np.ma.filled(nc.variables[var][it, :], np.nan)[idx]
            rec[name] = _derive_values(
                cfg, vals, zs[idx], derive, name=name,
                src_units=getattr(nc.variables[var], "units", None))
        frame = gpd.GeoDataFrame(rec, geometry=geoms, crs=grid.crs)
        frame.to_file(path, layer=layer, driver="GPKG",
                      mode="w" if n == 0 else "a")

        # Range of the first output variable, accumulated across the
        # timesteps as they are written, so the graduated ramp in the
        # style spans the whole animation rather than one frame - and
        # without a second pass over the data.
        if ramp_attr is not None and ramp_attr in rec:
            col = np.asarray(rec[ramp_attr], dtype="f8")
            col = col[np.isfinite(col)]
            if col.size:
                v_lo = col.min() if v_lo is None else min(v_lo, col.min())
                v_hi = col.max() if v_hi is None else max(v_hi, col.max())

        # Vertical faces have a ZERO-AREA footprint in the map plane, so a
        # 2D view of the polygons shows nothing at all. Write their
        # horizontal trace as lines too, in the same GeoPackage, so wall
        # temperatures can be mapped and clicked in a normal plan view.
        if want_lines:
            lines = frame.iloc[wall_idx].copy()
            lines["geometry"] = wall_geoms
            lines.to_file(path, layer=line_layer, driver="GPKG",
                          mode="w" if n == 0 else "a")
        n_written += idx.size
        debug("[surfaces] {}: timestep {}/{} written ({} rows so far)",
              cfg.disp(path), n + 1, len(its), n_written)

    qml = None
    if cfg["qgis_styles"]:
        step = _median_step(dts, its)
        # The polygon layer carries true Z, so it gets the 3D renderer;
        # the 2D trace layer shares the SAME ramp over the same variable
        # and range, so plan view and scene agree.
        qml = vector_qml("time", step, geom="polygon", attr=ramp_attr,
                         vmin=v_lo, vmax=v_hi,
                         clamp_absolute=bool(cfg["qgis_styles_3d"]))
        write_sidecar(path, qml)
        embed_gpkg_style(path, layer, qml)

    if want_lines:
        if qml is not None:
            embed_gpkg_style(path, line_layer,
                             vector_qml("time", step, geom="line",
                                        attr=ramp_attr, vmin=v_lo,
                                        vmax=v_hi))
        info("[surfaces] wrote {}: layer '{}' with {} wall traces "
             "(2D-visible)", cfg.disp(path), line_layer, int(wall_idx.size))
    counts = gpd.pd.Series(orient).value_counts().to_dict()
    info("[surfaces] wrote {}: {} elements x {} times ({})",
         cfg.disp(path), idx.size, len(its),
         ", ".join(f"{k} {v}" for k, v in sorted(counts.items())))


def _walls_gpkg(cfg, grid, nc, its, times, xs, ys, zs, zen, az, sel, path):
    idx = np.where(sel & (np.abs(zen - 90.0) < 1.0e-3))[0]
    _elements_gpkg(cfg, grid, nc, its, times, xs, ys, zs, zen, az, sel,
                   path, idx, "walls")


def _all_elements_gpkg(cfg, grid, nc, its, times, xs, ys, zs, zen, az, sel,
                       path):
    idx = np.where(sel)[0]
    _elements_gpkg(cfg, grid, nc, its, times, xs, ys, zs, zen, az, sel,
                   path, idx, "surfaces")


# ------------------------------
# 4. STAGE ENTRY
# ------------------------------
def stage_surfaces(cfg, grid):
    sf = cfg["surfaces"]
    src_path = cfg.resolve(sf["file_from"])

    # Read only the header first: which variables exist decides the output
    # names (whether `t_surf` is copied or derived), and that has to be
    # settled before the overwrite check picks the files to write.
    with Dataset(src_path, "r") as nc:
        check_product(nc, "surfaces", sf["file_from"])
        check_raw_input(nc, grid, "surfaces", sf["file_from"])
        available = set(nc.variables)
    plan = output_vars(sf, available, verbose=True)
    derived = [n for n, _v, d in plan if d]
    if derived:
        info("[surfaces] deriving {} from the surface potential "
             "temperature.", ", ".join(derived))

    targets = {}
    if sf["write_topview_raster"]:
        for name, _var, _derive in plan:
            path = cfg.out_path("{}{}_topview.tif".format(sf["prefix"],
                                                          name))
            if cfg.prepare_output(path, "surfaces"):
                targets[name] = path
    nc_path = cfg.out_path("{}surface_2d.nc".format(sf["prefix"]))
    do_nc = (sf["write_topview_netcdf"]
             and cfg.prepare_output(nc_path, "surfaces"))
    walls_path = cfg.out_path("{}walls.gpkg".format(sf["prefix"]))
    do_walls = (sf["write_walls_gpkg"]
                and cfg.prepare_output(walls_path, "surfaces"))
    all_path = cfg.out_path("{}surfaces.gpkg".format(sf["prefix"]))
    do_all = (sf["write_all_elements_gpkg"]
              and cfg.prepare_output(all_path, "surfaces"))
    mp_path = cfg.out_path(sf["multipatch_to"]
                           or "{}surfaces.shp".format(sf["prefix"]))
    do_mp = (sf["write_multipatch"]
             and cfg.prepare_output(mp_path, "surfaces"))
    fgdb_path = cfg.out_path(sf["fgdb_to"]
                             or "{}surfaces.gdb".format(sf["prefix"]))
    do_fgdb = (sf["write_fgdb"]
               and cfg.prepare_output(fgdb_path, "surfaces"))
    if not (targets or do_nc or do_walls or do_all or do_mp or do_fgdb):
        info("[surfaces] nothing to write, skipping.")
        return

    with Dataset(cfg.resolve(sf["file_from"]), "r") as nc:
        its, times = select_times(cfg, nc, "surfaces", grid, cfg["surfaces"]["file_from"])
        xs, ys, zs, zen, az = _load_geometry(nc)
        sel = _bounds_mask(cfg, xs, ys)

        other = sel & (np.abs(zen) > 1.0e-3) & (np.abs(zen - 90.0) > 1.0e-3)
        if other.any():
            vals = sorted(set(np.round(zen[other], 1)))
            warning("[surfaces] {} elements with zenith not in (0, 90) - "
                    "not converted (values: {}).", int(other.sum()), vals)

        if targets or do_nc:
            _topview_raster(cfg, grid, nc, its, times, xs, ys, zs, zen, sel,
                            targets, nc_path if do_nc else None)
        if do_walls:
            _walls_gpkg(cfg, grid, nc, its, times, xs, ys, zs, zen, az, sel,
                        walls_path)
        if do_all:
            _all_elements_gpkg(cfg, grid, nc, its, times, xs, ys, zs, zen,
                               az, sel, all_path)
        if do_mp or do_fgdb:
            idx_all = np.where(sel)[0]
            if do_mp:
                _elements_multipatch(cfg, grid, nc, its, times, xs, ys, zs,
                                     zen, az, idx_all, mp_path)
            if do_fgdb:
                _elements_fgdb(cfg, grid, nc, its, times, xs, ys, zs, zen,
                               az, idx_all, fgdb_path, "surfaces")
                wall_idx = np.where(
                    sel & (np.abs(zen - 90.0) < 1.0e-3))[0]
                if wall_idx.size:
                    _elements_fgdb(cfg, grid, nc, its, times, xs, ys, zs,
                                   zen, az, wall_idx, fgdb_path, "walls")


# ------------------------------
# 5. TIMESTEP SPACING
# ------------------------------
def _median_step(dts, its, default=3600.0):
    """Typical spacing between the selected records, in seconds. Used as
    the feature duration in the QGIS temporal style so an hourly layer
    stays visible for the whole hour instead of a single instant."""
    ts = [dts[i] for i in its if dts[i] is not None]
    if len(ts) < 2:
        return default
    gaps = sorted((ts[i + 1] - ts[i]).total_seconds()
                  for i in range(len(ts) - 1))
    return gaps[len(gaps) // 2] or default


# ------------------------------
# 6. RASTER CRS REPAIR
# ------------------------------
def _ensure_raster_crs(cfg, grid, path):
    """Verify the written raster carries a CRS, and repair it if not.

    A CRS that constructs fine can still fail to be encoded into a GeoTIFF
    on installations with a stale or mismatched PROJ database, leaving the
    file unreferenced. Reopen it, and if the CRS is missing try the
    remaining candidates in place (WKT, the config string, the EPSG form)
    before giving up with a warning.
    """
    candidates = grid.rio_crs_candidates()
    with rasterio.open(path) as chk:
        if chk.crs is not None:
            return True
    for cand in candidates:
        try:
            with rasterio.open(path, "r+") as dst:
                dst.crs = cand
            with rasterio.open(path) as chk:
                if chk.crs is not None:
                    debug("[surfaces] CRS repaired in place on {}",
                          cfg.disp(path))
                    return True
        except Exception:
            continue
    warning("[surfaces] {} has NO CRS - GDAL could not encode '{}'. The "
            "layer will not line up on a map; check the PROJ installation "
            "(stale proj.db, or a PROJ_DATA/PROJ_LIB left over from "
            "another install).", cfg.disp(path), grid.crs_raw)
    return False
