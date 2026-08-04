"""Static driver processing and georeferencing.

GIS counterpart of paraview_vis/utils/static_driver.py. Reads grid geometry
(nx, ny, dx, dy and optionally nz, dz from the _p3d namelist), the combined
terrain + building elevation, the domain origin (origin_x, origin_y,
origin_time) and converts the metre-based `bounds` section into grid index
ranges. Everything downstream is written in ABSOLUTE projected coordinates
(crs.palm), so ArcGIS/QGIS place the layers correctly with no manual
shifting.
"""

import re
import datetime

import numpy as np
from netCDF4 import Dataset
from pyproj import CRS

from ..log import debug, get_logger, info, progress, warning
from .thermo import TEMP_PREFIXES as _TEMP_VAR_PREFIXES

log = get_logger()


# ------------------------------
# 1. GRID CONTAINER
# ------------------------------
class Grid:
    """Everything the conversion stages need to know about the domain:
    grid geometry, origin, origin_time, CRS, the terrain + building
    elevation of the cut subdomain, and the bounds index ranges."""

    __slots__ = ("nx", "ny", "nz", "dx", "dy", "dz",
                 "origin_x", "origin_y", "origin_z", "origin_time",
                 "crs", "crs_raw",
                 "zt", "elev", "x_stat", "y_stat", "p3d",
                 "nx0", "nxmax", "ny0", "nymax", "nz0", "nzmax",
                 "zt_min", "zt_max")

    def rio_crs_candidates(self):
        """CRS candidates for rasterio, best first.

        Handing rasterio a CRS that merely *constructs* is not enough: on
        some GDAL/PROJ installations a PROJ-string CRS builds fine but
        cannot be encoded into a GeoTIFF, and the file ends up with NO
        georeferencing. Explicit WKT is tried first because it needs no
        proj.db lookup, and each candidate is checked for a non-empty WKT
        before being offered.
        """
        import rasterio.crs

        out = []
        for src in (self.crs.to_wkt(), self.crs_raw,
                    "EPSG:{}".format(self.crs.to_epsg())
                    if self.crs.to_epsg() else None):
            if not src:
                continue
            try:
                c = rasterio.crs.CRS.from_user_input(src)
                if c and c.to_wkt():
                    out.append(c)
            except Exception:
                continue
        if not out:
            warning("[static] could not build a rasterio CRS from {} - "
                    "rasters may be written without georeferencing.",
                    self.crs_raw)
        return out

    def rio_crs(self):
        """Best CRS for rasterio, or None."""
        c = self.rio_crs_candidates()
        return c[0] if c else None

    def abs_x(self, x_local):
        return np.asarray(x_local) + self.origin_x

    def abs_y(self, y_local):
        return np.asarray(y_local) + self.origin_y

    def times_to_datetime(self, seconds, tz_aware=False):
        """Absolute timestamps for PALM's time-in-seconds. Entries that are
        NaN or absurdly large (an unwritten time variable left at its
        _FillValue) come back as None rather than overflowing datetime."""
        out = []
        for s in np.atleast_1d(seconds):
            s = float(s)
            if not np.isfinite(s) or abs(s) > 1.0e12:
                out.append(None)
            else:
                dt = self.origin_time + datetime.timedelta(seconds=s)
                # PALM origin_time is UTC. Tagging the datetimes as UTC
                # keeps GDAL from writing them as naive local times into
                # the GeoPackage, which shifts a temporal animation by the
                # local offset.
                out.append(dt.replace(tzinfo=datetime.timezone.utc)
                           if tz_aware else dt)
        return out

    @staticmethod
    def label(dt, it):
        """Timestamp string for an output label, or a record tag when the
        source file carried no usable time value."""
        return (dt.strftime("%Y-%m-%d %H:%M:%S") if dt is not None
                else "record {}".format(it))

    def z_levels(self):
        """PALM zu levels: 0, dz/2, 3dz/2, ..."""
        return np.append(0.0, np.arange(self.dz / 2.0,
                                        (self.nz - 1) * self.dz, self.dz))


# ------------------------------
# 2. P3D NAMELIST PARSER
# ------------------------------
_P3D_KEYS = ("nx", "ny", "nz", "dx", "dy", "dz")

# Namelist keys holding output timing, per product and per averaging mode.
# PALM keeps SEPARATE intervals for 3D and surface output, which is why one
# file can carry 10-minute records and another hourly ones in the same run.
TIME_KEYS = {
    "voxel": (("dt_data_output_av", "skip_time_data_output_av"),
              ("dt_data_output", "skip_time_data_output")),
    "wind": (("dt_data_output_av", "skip_time_data_output_av"),
             ("dt_data_output", "skip_time_data_output")),
    "surfaces": (("dt_dosurf_av", "skip_time_dosurf_av"),
                 ("dt_dosurf", "skip_time_dosurf")),
}


def read_p3d(path):
    """Parse every `key = value` numeric setting from a PALM namelist.

    Assignments are matched wherever they occur, not only when they are
    alone on their line: PALM namelists routinely pack several onto one
    line, e.g.

        nx = 511, ny = 607, nz = 512,
        dx = 1.0, dy = 1.0, dz = 1.0,

    An earlier version anchored the pattern to whole lines and so read
    nothing at all from those two.

    Comment lines (`!`) are stripped first, so a commented-out or
    documentary value cannot override the real one. Values that are not
    plain numbers - strings, arrays, logicals - are skipped; only the
    scalar numeric settings are of interest here.
    """
    out = {}
    for raw_line in open(path, "r"):
        line = raw_line.split("!", 1)[0]          # drop trailing comments
        for m in re.finditer(
                r"([A-Za-z_][A-Za-z0-9_]*)\s*=\s*"
                r"([-+]?[0-9]*\.?[0-9]+(?:[eEdD][-+]?[0-9]+)?)\s*(?=[,\s]|$)",
                line):
            key = m.group(1).lower()
            raw = m.group(2).replace("d", "e").replace("D", "e")
            try:
                val = float(raw)
            except ValueError:
                continue
            # First occurrence wins: PALM applies the first namelist block,
            # and a later &..._parameters block may reuse a name.
            if key not in out:
                out[key] = (int(val) if key.startswith("n")
                            and val.is_integer() else val)
    return out


# ------------------------------
# 3. STATIC DRIVER + BOUNDS
# ------------------------------
def process_static_driver(cfg):
    """Read the static driver (and optionally the p3d file); return a Grid."""
    progress("Reading the static driver")

    dom = cfg["domain"]
    bnd = cfg["bounds"]
    g = Grid()
    g.crs_raw = str(cfg["crs"]["palm"])
    g.crs = CRS.from_user_input(g.crs_raw)

    # nz, dz either from p3d or given explicitly (config wins)
    nz, dz = dom["nz"], dom["dz"]
    g.p3d = {}
    if dom["p3d"]:
        g.p3d = read_p3d(cfg.resolve(dom["p3d"]))
        log.debug("p3d: " + ", ".join(
            "{}={}".format(k, g.p3d[k]) for k in sorted(g.p3d)))
        nz = nz if nz is not None else g.p3d.get("nz")
        dz = dz if dz is not None else g.p3d.get("dz")
    if nz is None or dz is None:
        raise RuntimeError("nz/dz not found in p3d and not set in config.")
    g.nz, g.dz = int(nz), float(dz)

    with Dataset(cfg.resolve(dom["static_driver"]), "r") as nc:
        x = np.asarray(nc.variables["x"][:])
        y = np.asarray(nc.variables["y"][:])
        g.nx, g.ny = x.size, y.size
        g.dx = float(x[1] - x[0])
        g.dy = float(y[1] - y[0])

        # horizontal bounds (metres, local coords) -> index ranges
        ix = np.where((x >= bnd.get("x0", -1.0e30))
                      & (x <= bnd.get("xmax", 1.0e30)))[0]
        iy = np.where((y >= bnd.get("y0", -1.0e30))
                      & (y <= bnd.get("ymax", 1.0e30)))[0]
        if ix.size == 0 or iy.size == 0:
            raise RuntimeError("bounds cut leaves an empty domain.")
        g.nx0, g.nxmax = int(ix[0]), int(ix[-1]) + 1
        g.ny0, g.nymax = int(iy[0]), int(iy[-1]) + 1

        # terrain + buildings elevation of the cut subdomain
        zt = np.ma.filled(nc.variables["zt"][g.ny0:g.nymax, g.nx0:g.nxmax],
                          0.0)
        if "buildings_2d" in nc.variables:
            zb = np.ma.filled(
                nc.variables["buildings_2d"][g.ny0:g.nymax, g.nx0:g.nxmax],
                0.0)
            if not np.any(zb > 0.0):
                # Static drivers that describe buildings only in 3D
                # (building_id / building_surface_level) leave buildings_2d
                # entirely at its fill value. The local surface is then bare
                # terrain, so heights_agl reference the ground even on roofs.
                warning("[static] buildings_2d is empty in the cut - "
                        "heights above ground reference bare TERRAIN, not "
                        "rooftops.")
        else:
            warning("[static] buildings_2d not found in the static "
                    "driver - elevation is terrain only.")
            zb = np.zeros_like(zt)
        g.zt = zt
        g.elev = zt + zb
        g.x_stat = x[g.nx0:g.nxmax]
        g.y_stat = y[g.ny0:g.nymax]

        # vertical index range from zmax (height above the highest surface)
        g.zt_min = float(np.min(zt))
        g.zt_max = float(np.max(g.elev))
        z = g.z_levels()
        g.nz0 = int(max(np.searchsorted(z, g.zt_min) - 1, 0))
        zmax = bnd.get("zmax", 1.0e30)
        g.nzmax = int(min(np.searchsorted(z, g.zt_max + zmax) + 1, z.size))

        # domain origin
        g.origin_x = float(nc.origin_x)
        g.origin_y = float(nc.origin_y)
        g.origin_z = float(getattr(nc, "origin_z", 0.0))
        g.origin_time = _parse_origin_time(dom["origin_time"],
                                           getattr(nc, "origin_time", None))

    info("[static] {} x {} x {} pts, d = {}/{}/{} m, origin "
         "({:.2f}, {:.2f}) {}", g.nx, g.ny, g.nz, g.dx, g.dy, g.dz,
         g.origin_x, g.origin_y, cfg["crs"]["palm"])
    info("[static] cut x[{}:{}] y[{}:{}] z[{}:{}], t0 = {}",
         g.nx0, g.nxmax, g.ny0, g.nymax, g.nz0, g.nzmax, g.origin_time)
    return g


def _parse_origin_time(override, nc_attr):
    """origin_time from the config override or the static driver attribute.
    PALM-style '2023-08-23 17:00:00 +02' is accepted; the UTC-offset suffix
    is dropped."""
    s = override or nc_attr
    if not s:
        warning("[static] origin_time not found - using "
                "2000-01-01 00:00:00 (set domain.origin_time).")
        return datetime.datetime(2000, 1, 1)
    m = re.match(r"(\d{4}-\d{2}-\d{2}[ T]\d{2}:\d{2}:\d{2})", str(s).strip())
    if not m:
        raise RuntimeError(f"Cannot parse origin_time: {s}")
    return datetime.datetime.strptime(m.group(1).replace("T", " "),
                                      "%Y-%m-%d %H:%M:%S")


# ------------------------------
# 4. SHARED TIME SELECTION
# ------------------------------
def infer_time_params(grid, tag, source=""):
    """Output interval and first-record time for a stage, from the p3d.

    PALM writes averaged output at
        skip_time_* + n * dt_*
    with the FIRST record one full interval after the skip time, and it
    keeps separate settings for 3D output (`dt_data_output_av`) and
    surface output (`dt_dosurf_av`) - which is why one file in a run can
    hold 10-minute records and another hourly ones.

    Returns (step, start, key) or (None, None, None) when the namelist
    does not say. The averaged keys are preferred for a file whose name
    marks it as averaged.
    """
    p3d = getattr(grid, "p3d", None) or {}
    if not p3d:
        return None, None, None

    variants = TIME_KEYS.get(tag)
    if not variants:
        return None, None, None
    if "_av" not in str(source):          # instantaneous output: try it first
        variants = tuple(reversed(variants))

    for dt_key, skip_key in variants:
        step = p3d.get(dt_key)
        if step and step > 0:
            skip = float(p3d.get(skip_key, 0.0) or 0.0)
            return float(step), skip + float(step), dt_key
    return None, None, None


def select_times(cfg, nc, tag="", grid=None, source=""):
    """Timestep indices selected by bounds.skip_times / n_times, plus the
    time values. Shared by all stages reading PALM output.

    PALM pads unwritten records at the end of a file with the time variable's
    _FillValue (a large integer), so trailing fill entries are dropped. A
    file whose time variable is fill/masked THROUGHOUT still has usable data
    records - its raw values are used with a warning rather than leaving
    nothing to convert. Different output files of the same run routinely
    hold different numbers of records (3D and surface output use separate
    intervals), so the window is clamped per file and reported.
    """
    label = "[{}] ".format(tag) if tag else ""

    if "time" not in nc.variables:
        raise RuntimeError("No 'time' variable in the file.")
    tvar = nc.variables["time"][:]
    raw = np.atleast_1d(np.squeeze(np.ma.getdata(tvar).astype("f8")))
    mask = np.atleast_1d(np.squeeze(np.ma.getmaskarray(tvar)))
    if mask.shape != raw.shape:
        mask = np.zeros(raw.shape, dtype=bool)

    # fill/sentinel detection: masked, non-finite, or absurdly large
    bad = mask | ~np.isfinite(raw) | (np.abs(raw) > 1.0e15)
    fv = getattr(nc.variables["time"], "_FillValue", None)
    if fv is not None:
        bad |= raw == float(fv)

    n_all = raw.size
    step = cfg["bounds"].get("time_step")
    start = cfg["bounds"].get("time_start") or 0.0
    src_key = None
    if step is None and grid is not None:
        # Nothing set by hand: ask the namelist for this stage's interval.
        step, start_auto, src_key = infer_time_params(grid, tag, source)
        if step:
            start = (cfg["bounds"].get("time_start")
                     if cfg["bounds"].get("time_start") else start_auto)
    if bad.all():
        # Some PALM surface files are written with the time variable left at
        # its _FillValue (~1e36). Those numbers are not times - using them
        # would put every record ~3e29 years after origin_time (and overflow
        # datetime). Reconstruct from the configured interval if given,
        # otherwise carry NaN and let each stage fall back to record indices.
        if step:
            times = start + float(step) * np.arange(n_all)
            origin = ("the p3d ({} = {:g} s)".format(src_key, step)
                      if src_key else
                      "bounds.time_start/time_step")
            warning("{}time variable is all fill - times rebuilt from {}: "
                    "first record at {:g} s, every {:g} s.",
                    label, origin, start, step)
        else:
            times = np.full(n_all, np.nan)
            warning("{}time variable is all fill, and the p3d gives no "
                    "output interval for this product - labelling by "
                    "record index (set bounds.time_step to rebuild real "
                    "times).", label)
        n_valid = n_all
    else:
        # keep the leading run of good records (trailing fill is the norm)
        n_valid = int(np.argmax(bad)) if bad.any() else n_all
        if bad[:n_valid].any():          # fill in the middle: keep all good
            n_valid = int((~bad).sum())
            times = raw[~bad]
        else:
            times = raw
        if n_valid < n_all:
            debug("{}{} of {} records are fill/masked and were dropped.",
                  label, n_all - n_valid, n_all)

    skip = cfg["bounds"]["skip_times"]

    # `n_times` is a COUNT of records to keep, starting at skip_times. It
    # used to be passed straight to range() as the END index, so it behaved
    # as "stop at record n_times" instead: skip_times=10 with n_times=5
    # selected nothing at all, and the error message had to explain that
    # n_times must exceed skip_times. The name said count, the code meant
    # end index. `time_end` is available for anyone who wants the old
    # exclusive-end meaning explicitly.
    end = cfg["bounds"].get("time_end")
    if end is not None:
        stop = min(int(end), n_valid)
    else:
        want = cfg["bounds"]["n_times"]
        stop = n_valid if want >= 99999 else min(skip + int(want), n_valid)
    its = list(range(skip, stop))
    if not its:
        raise RuntimeError(
            "No timesteps selected: this file holds {} usable record(s) "
            "(of {} in the time dimension), while bounds.skip_times={} and "
            "bounds.{}={}. Lower skip_times, or raise the count.".format(
                n_valid, n_all, skip,
                "time_end" if end is not None else "n_times",
                end if end is not None else cfg["bounds"]["n_times"]))
    if len(its) == n_all:
        debug("{}all {} records selected.", label, n_all)
    else:
        info("{}records {}..{} of {} selected.", label, its[0], its[-1],
             n_all)
    return its, times


# ------------------------------
# 5. RAW-INPUT GUARD
# ------------------------------
# palm2gis reads RAW (or join-only) PALM output: local metre coordinates,
# temperatures in kelvin, original variable layout. Files processed by the
# palm_postproc `coord` step carry ABSOLUTE projected coordinates and
# degrees-C temperatures; feeding one in here would add origin_x/origin_y a
# second time (layers land ~440 km off for UTM 33N Prague) and break the
# ta-from-theta derivation. `coord` marks its work unambiguously, so this
# is detected and rejected with an explanation instead of failing silently.

# _TEMP_VAR_PREFIXES is imported from thermo.py (single source of truth).
_CELSIUS_UNITS = ("degrees_c", "degc", "celsius", "\N{DEGREE SIGN}c")


def check_raw_input(nc, grid, tag, source=""):
    """Raise RuntimeError if *nc* looks like palm_postproc `coord` output.

    Three independent signals, any one of which rejects the file:
      a) x/y coordinate variable with standard_name projection_*_coordinate
         (written by coord when it overwrites local coords with UTM);
      b) x[0] far outside the local domain extent (absolute easting);
      c) a temperature variable in degrees C (coord's K -> degC conversion).

    Surface-output files have no `x`/`y` dimension coordinates (elements
    are listed along `s` with xs/ys), so (a)/(b) are skipped when `x` is
    absent and (c) alone decides.
    """
    label = "[{}] ".format(tag) if tag else ""
    src = " ({})".format(source) if source else ""
    reject = None

    # a) + b) - dimension coordinates
    if "x" in nc.variables:
        xv = nc.variables["x"]
        std = str(getattr(xv, "standard_name", "")).lower()
        if std == "projection_x_coordinate":
            reject = ("x carries standard_name=projection_x_coordinate, "
                      "i.e. absolute projected coordinates")
        else:
            x0 = float(np.asarray(xv[:1])[0])
            extent = grid.nx * grid.dx
            if abs(x0) > max(10.0 * extent, 1.0e4):
                reject = ("x starts at {:.1f} m, far beyond the local "
                          "domain extent of {:.1f} m - these are absolute "
                          "coordinates".format(x0, extent))

    # c) - converted temperature units
    if reject is None:
        for name in nc.variables:
            if not name.lower().startswith(_TEMP_VAR_PREFIXES):
                continue
            units = str(getattr(nc.variables[name], "units", "")).lower()
            if units.strip() in _CELSIUS_UNITS:
                reject = ("variable '{}' has units '{}', i.e. already "
                          "converted from kelvin".format(
                              name, getattr(nc.variables[name], "units")))
                break

    if reject is not None:
        raise RuntimeError(
            "{}this input file{} appears to be palm_postproc `coord` "
            "output: {}. palm2gis must read raw PALM output (or "
            "palm_postproc `join` output) together with the static driver "
            "- it applies origin_x/origin_y and unit handling itself. "
            "Point file_from at the raw OUTPUT/ or OUTPUT_join/ file, not "
            "at OUTPUT_coord/.".format(label, src, reject))
    debug("{}raw-input guard passed.", label)
