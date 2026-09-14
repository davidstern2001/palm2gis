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
from .thermo import CELSIUS_UNITS as _THERMO_CELSIUS_UNITS
from .thermo import CELSIUS_PREFIXES as _THERMO_CELSIUS_PREFIXES

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
                 "zt_min", "zt_max", "zmax_req")

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

    def abs_z(self, z_local):
        """PALM height above the domain base -> height above sea level.

        PALM's z = 0 sits at origin_z, so every zu_3d / zs / terrain value
        is RELATIVE. x and y were already being shifted by origin_x /
        origin_y while z was written raw, which put voxels, LAD, surface
        quads and wind arrows origin_z metres too low - 175 m for the
        Holesovice domain.

        NOT for thermodynamics. exner_profile / exner_simple pair z with
        `surface_pressure`, which is defined as the pressure at the
        model's z = 0; feeding them absolute heights would put a second,
        larger error on every derived temperature. Keep a relative z for
        those and use this only for output geometry.
        """
        return np.asarray(z_local) + self.origin_z

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
    reset_time_windows()          # a fresh run: forget the previous windows

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

    cfg.resolve_surface_pressure(g.p3d)

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
        g.zmax_req = float(zmax)
        g.nzmax = int(min(np.searchsorted(z, g.zt_max + zmax) + 1, z.size))

        # domain origin
        g.origin_x = float(nc.origin_x)
        g.origin_y = float(nc.origin_y)
        g.origin_z = float(dom["origin_z"]) if dom.get("origin_z") is not None \
            else float(getattr(nc, "origin_z", 0.0))
        g.origin_time = _parse_origin_time(dom["origin_time"],
                                           getattr(nc, "origin_time", None))

    info("[static] {} x {} x {} pts, d = {}/{}/{} m, origin "
         "({:.2f}, {:.2f}) {}", g.nx, g.ny, g.nz, g.dx, g.dy, g.dz,
         g.origin_x, g.origin_y, cfg["crs"]["palm"])
    info("[static] cut x[{}:{}] y[{}:{}] z[{}:{}], t0 = {}",
         g.nx0, g.nxmax, g.ny0, g.nymax, g.nz0, g.nzmax, g.origin_time)
    if g.origin_z:
        info("[static] origin_z = {:.2f} m - z output is height above sea "
             "level.", g.origin_z)
    else:
        warning("[static] origin_z is 0 (absent or unset in the static "
                "driver): z output will be height above the DOMAIN BASE, "
                "not above sea level, and 3D layers will sit too low in "
                "ArcGIS Pro. Set domain.origin_z if the domain is above "
                "sea level.")
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


# ------------------------------
# 4b. TIME WINDOW PARSING
# ------------------------------
# `skip_times` / `n_times` count RECORDS, and a record means a different
# amount of wall-clock time in each file: PALM keeps separate output
# intervals for 3D and surface output, so record 18 of a 10-minute 3D file
# is 3 hours in while record 18 of an hourly surface file is 18 hours in.
# Selecting the same index range from both therefore converts two
# completely different periods. `bounds.time_from` / `time_to` state the
# window ONCE, in time, and every stage resolves it against its own file.
_UNIT_S = {"s": 1.0, "sec": 1.0, "secs": 1.0, "second": 1.0, "seconds": 1.0,
           "m": 60.0, "min": 60.0, "mins": 60.0, "minute": 60.0,
           "minutes": 60.0,
           "h": 3600.0, "hr": 3600.0, "hrs": 3600.0, "hour": 3600.0,
           "hours": 3600.0,
           "d": 86400.0, "day": 86400.0, "days": 86400.0}

_OFFSET_RE = re.compile(r"^([-+]?\d+(?:\.\d+)?)\s*([a-zA-Z]*)$")
_STAMP_RE = re.compile(r"^(\d{4})-(\d{1,2})-(\d{1,2})[ T](\d{1,2}):(\d{2})"
                       r"(?::(\d{2}))?$")
_CLOCK_RE = re.compile(r"^(\d{1,2}):(\d{2})(?::(\d{2}))?$")


def parse_time_bound(value, origin_time, key, after=None):
    """One end of the bounds.time_from / time_to window.

    Returns (seconds_since_origin_time, absolute_datetime_or_None), or
    (None, None) when the value is unset. Accepted forms:

        3600            seconds since origin_time (a bare number)
        "6h" / "90min"  the same, with a unit (s/min/h/d)
        "2023-08-24 13:00:00"   an absolute timestamp
        "13:00"         a clock time - the FIRST occurrence at or after
                        origin_time (or after `after`, so a window like
                        "22:00" -> "02:00" crosses midnight correctly)

    YAML quirks worth knowing: an unquoted `2023-08-24 13:00:00` is parsed
    by YAML itself into a datetime (handled here), while an unquoted
    `13:00:00` becomes the sexagesimal integer 46800 - which would be read
    as seconds. Clock times must be QUOTED.
    """
    if value is None or value == "":
        return None, None

    origin = getattr(origin_time, "year", None) and origin_time or None

    def _to_pair(dt):
        if dt.tzinfo is not None:                # compare like with like
            dt = dt.replace(tzinfo=None)
        if origin is None:
            raise RuntimeError(
                "{} is a timestamp, but the domain origin_time is unknown - "
                "give domain.origin_time, or use seconds instead."
                .format(key))
        return (dt - origin).total_seconds(), dt

    if isinstance(value, bool):
        raise RuntimeError("{}: expected a time, got {!r}.".format(key, value))
    if isinstance(value, datetime.datetime):     # YAML parsed it for us
        return _to_pair(value)
    if isinstance(value, datetime.date):
        return _to_pair(datetime.datetime(value.year, value.month, value.day))
    if isinstance(value, (int, float)):
        return float(value), None

    s = str(value).strip()

    m = _STAMP_RE.match(s)
    if m:
        y, mo, d, hh, mm, ss = m.groups()
        return _to_pair(datetime.datetime(int(y), int(mo), int(d), int(hh),
                                          int(mm), int(ss or 0)))

    m = _CLOCK_RE.match(s)
    if m:
        hh, mm, ss = m.groups()
        if origin is None:
            raise RuntimeError(
                "{}: a clock time needs the domain origin_time to be "
                "resolved.".format(key))
        ref = after if after is not None else origin
        cand = datetime.datetime.combine(
            ref.date(), datetime.time(int(hh), int(mm), int(ss or 0)))
        while cand < ref:                        # roll into the next day(s)
            cand += datetime.timedelta(days=1)
        return _to_pair(cand)

    m = _OFFSET_RE.match(s)
    if m:
        num, unit = m.groups()
        unit = (unit or "s").lower()
        if unit not in _UNIT_S:
            raise RuntimeError(
                "{}: unknown time unit {!r} in {!r} (use s, min, h or d)."
                .format(key, unit, s))
        return float(num) * _UNIT_S[unit], None

    raise RuntimeError(
        "{}: cannot read {!r} as a time. Use seconds (3600), a unit "
        "(\"6h\"), a timestamp (\"2023-08-24 13:00:00\") or a quoted clock "
        "time (\"13:00\").".format(key, s))


# Windows already resolved in this run, so a later stage can notice that it
# is converting a different period than an earlier one - the exact failure
# that index-based selection causes across files of unequal interval.
_RUN_WINDOWS = {}


def reset_time_windows():
    """Forget the per-stage windows of a previous run (same process)."""
    _RUN_WINDOWS.clear()


def _fmt_stamp(grid, seconds):
    """`seconds` as a timestamp when origin_time is known, else as s."""
    if grid is not None and getattr(grid, "origin_time", None) is not None \
            and np.isfinite(seconds):
        dt = grid.times_to_datetime([seconds])[0]
        if dt is not None:
            return dt.strftime("%Y-%m-%d %H:%M:%S")
    return "{:g} s".format(seconds)


def select_times(cfg, nc, tag="", grid=None, source=""):
    """Timestep indices to convert, plus the time value of every record.

    Two ways to pick the window:

    * bounds.time_from / time_to - a window in TIME. Each file is matched
      on its own time axis, so a 10-minute 3D file and an hourly surface
      file both come back covering the same period. Prefer this whenever
      more than one stage runs.
    * bounds.skip_times / n_times (or time_end) - a window in RECORD
      INDICES, clamped per file. Simple, but index n means a different
      moment in files written at different intervals.

    `times` is returned FULL LENGTH and aligned with the file's record
    dimension (NaN where a record carries no usable time), because callers
    index it with the returned record indices.

    PALM pads unwritten records at the end of a file with the time variable's
    _FillValue (a large integer), so trailing fill entries are dropped. A
    file whose time variable is fill/masked THROUGHOUT still has usable data
    records - its times are rebuilt from the p3d rather than leaving nothing
    to convert.
    """
    label = "[{}] ".format(tag) if tag else ""
    bnd = cfg["bounds"]

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
    step = bnd.get("time_step")
    start = bnd.get("time_start") or 0.0
    src_key = None
    if step is None and grid is not None:
        # Nothing set by hand: ask the namelist for this stage's interval.
        step, start_auto, src_key = infer_time_params(grid, tag, source)
        if step:
            start = (bnd.get("time_start")
                     if bnd.get("time_start") else start_auto)
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
        good = np.ones(n_all, dtype=bool)
        n_valid = n_all
    else:
        # Times stay FULL LENGTH with NaN in the fill slots: the indices
        # returned here are used to index the DATA arrays as well, so
        # compacting the time array would shift every label after a gap.
        times = np.where(bad, np.nan, raw)
        good = ~bad
        n_valid = int(np.argmax(bad)) if bad.any() else n_all
        if bad[:n_valid].any():          # fill in the middle: keep all good
            n_valid = int(good.sum())
        else:                            # trailing fill: nothing past n_valid
            good = good & (np.arange(n_all) < n_valid)
        if int(good.sum()) < n_all:
            debug("{}{} of {} records are fill/masked and were dropped.",
                  label, n_all - int(good.sum()), n_all)

    # ---- window in TIME (preferred: file-interval independent) ----------
    origin_time = getattr(grid, "origin_time", None) if grid is not None \
        else None
    t0, dt0 = parse_time_bound(bnd.get("time_from"), origin_time,
                               "bounds.time_from")
    t1, dt1 = parse_time_bound(bnd.get("time_to"), origin_time,
                               "bounds.time_to", after=dt0)
    by_time = t0 is not None or t1 is not None

    if by_time:
        if t0 is not None and t1 is not None and t1 < t0:
            raise RuntimeError(
                "bounds.time_to ({:g} s) is before bounds.time_from ({:g} s)."
                .format(t1, t0))
        for key in ("skip_times", "n_times", "time_end"):
            val = bnd.get(key)
            if val not in (None, 0, 99999):
                warning("{}bounds.{} = {} is ignored: the window is set by "
                        "bounds.time_from/time_to.", label, key, val)
        usable = good & np.isfinite(times)
        if not usable.any():
            raise RuntimeError(
                "bounds.time_from/time_to need real timestamps, but this "
                "file carries none and none could be rebuilt. Set "
                "bounds.time_step (and time_start), or select by "
                "skip_times/n_times instead.")
        # PALM writes output at the first model timestep AT OR AFTER the
        # target, so record times carry a residue: this run's hourly
        # surface file holds 3600.73, 7200.73, ... A record nominally at
        # 16:00 is really at 16:00:00.73, and an exact comparison against
        # `time_to: "16:00:00"` silently dropped it - the window looked
        # exclusive at the end when it is meant to be inclusive. Match with
        # a tolerance tied to the file's OWN spacing, so the residue cannot
        # cost a record while the next one stays safely out.
        tol = 1.0
        tv_all = times[usable]
        if tv_all.size > 1:
            step_file = float(np.median(np.diff(tv_all)))
            if np.isfinite(step_file) and step_file > 0:
                tol = float(np.clip(0.5 * step_file, 1.0, 60.0))
        keep = usable.copy()
        if t0 is not None:
            keep &= times >= t0 - tol
        if t1 is not None:
            keep &= times <= t1 + tol
        its = [int(i) for i in np.flatnonzero(keep)]
        if not its:
            tv = times[usable]
            raise RuntimeError(
                "No timesteps in the requested window: bounds.time_from/"
                "time_to ask for {} .. {}, while this file covers {} .. {} "
                "({} record(s)).".format(
                    _fmt_stamp(grid, t0) if t0 is not None else "start",
                    _fmt_stamp(grid, t1) if t1 is not None else "end",
                    _fmt_stamp(grid, float(tv.min())),
                    _fmt_stamp(grid, float(tv.max())), tv.size))
    # ---- window in RECORD INDICES (the original behaviour) --------------
    else:
        skip = bnd["skip_times"]

        # `n_times` is a COUNT of records to keep, starting at skip_times. It
        # used to be passed straight to range() as the END index, so it
        # behaved as "stop at record n_times" instead: skip_times=10 with
        # n_times=5 selected nothing at all, and the error message had to
        # explain that n_times must exceed skip_times. The name said count,
        # the code meant end index. `time_end` is available for anyone who
        # wants the old exclusive-end meaning explicitly.
        end = bnd.get("time_end")
        if end is not None:
            stop = min(int(end), n_valid)
        else:
            want = bnd["n_times"]
            stop = n_valid if want >= 99999 else min(skip + int(want), n_valid)
        its = [i for i in range(skip, stop) if good[i]]
        if not its:
            raise RuntimeError(
                "No timesteps selected: this file holds {} usable record(s) "
                "(of {} in the time dimension), while bounds.skip_times={} "
                "and bounds.{}={}. Lower skip_times, raise the count, or "
                "select by time with bounds.time_from/time_to.".format(
                    n_valid, n_all, skip,
                    "time_end" if end is not None else "n_times",
                    end if end is not None else bnd["n_times"]))

    # ---- report, and cross-check against the other stages ---------------
    sel_t = times[np.asarray(its)]
    span = ""
    if np.isfinite(sel_t).all():
        span = " ({} .. {})".format(_fmt_stamp(grid, float(sel_t[0])),
                                    _fmt_stamp(grid, float(sel_t[-1])))
    if len(its) == n_all:
        debug("{}all {} records selected{}.", label, n_all, span)
    else:
        info("{}records {}..{} of {} selected{}.", label, its[0], its[-1],
             n_all, span)

    if np.isfinite(sel_t).all() and tag:
        lo, hi = float(sel_t[0]), float(sel_t[-1])
        for other, (o_lo, o_hi, o_tag) in _RUN_WINDOWS.items():
            if other == tag:
                continue
            if hi < o_lo or lo > o_hi:    # no overlap at all
                warning("{}this window ({} .. {}) does not overlap the one "
                        "[{}] is converting ({} .. {}) - record indices mean "
                        "different times in files of unequal output "
                        "interval. Set bounds.time_from/time_to to pin the "
                        "same period in every stage.",
                        label, _fmt_stamp(grid, lo), _fmt_stamp(grid, hi),
                        o_tag, _fmt_stamp(grid, o_lo), _fmt_stamp(grid, o_hi))
                break
        _RUN_WINDOWS[tag] = (lo, hi, tag)

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

# _TEMP_VAR_PREFIXES and the unit lists are imported from thermo.py
# (single source of truth). The local copy this replaces listed
# "degrees_c" but not "degree_c", so it missed PALM's own truncated
# "degree_" entirely.
#
# `ta*` is deliberately EXCLUDED from the Celsius signal below. PALM
# writes `ta` and `ta_2m*` in degrees Celsius natively, so "this
# temperature is in degC" is not evidence of post-processing for those -
# it is what a raw file looks like. Only the variables PALM writes in
# kelvin (theta*, tsurf*, t_surf*) can testify here.
_COORD_EVIDENCE_PREFIXES = ("theta", "tsurf", "t_surf")


# Which PALM product each stage needs, and how to recognise it. `file_from`
# is one key name across four stages that want three DIFFERENT files, which
# is easy to get wrong and previously surfaced as a confusing KeyError deep
# inside a stage ("no variable 'theta_surface'") rather than "this is the
# 3D file, surfaces wants the SURFACE_DATA one".
PRODUCTS = {
    "3d": ("3D output", "DATA_3D_NETCDF*",
           "gridded volume: dimensions x, y, zu_3d/zw_3d"),
    "surface": ("surface output", "SURFACE_DATA_NETCDF*",
                "unstructured surface elements listed along the s "
                "dimension"),
    "static": ("static driver", "*_static*",
               "the domain description: zt, buildings, lad - no time axis"),
}

STAGE_PRODUCT = {"voxel": "3d", "wind": "3d", "surfaces": "surface",
                 "canopy": "static"}


def identify_product(nc):
    """Which PALM product a file is, from its dimensions. None if unclear."""
    dims = set(nc.dimensions)
    if "s" in dims and "x" not in dims:
        return "surface"
    if {"zu_3d", "zw_3d"} & dims:
        return "3d"
    if "zlad" in dims or ("time" not in dims and "x" in dims):
        return "static"
    return None


def check_product(nc, tag, source=""):
    """Raise RuntimeError if *nc* is not the PALM product *tag* needs.

    Recognition is by DIMENSIONS rather than by filename, so a renamed or
    joined file is judged on what it actually holds.
    """
    want = STAGE_PRODUCT.get(tag)
    if want is None:
        return
    got = identify_product(nc)
    if got is None or got == want:
        return                      # unrecognised: let the stage try anyway

    w_name, w_glob, w_how = PRODUCTS[want]
    g_name = PRODUCTS[got][0]
    raise RuntimeError(
        "[{tag}] this looks like the {got} file, but the {tag} stage needs "
        "the {want} ({glob}) - {how}.{src}\n"
        "  Set {tag}.file_from to your {want} file. The stages read "
        "different PALM products: voxel and wind take the 3D output, "
        "surfaces takes the surface output, canopy takes the static "
        "driver.".format(
            tag=tag, got=g_name, want=w_name, glob=w_glob, how=w_how,
            src="\n  Given: {}".format(source) if source else ""))


def resolve_z_range(nc, grid, tag, z_coord="zu_3d"):
    """Vertical index range for *nc*, against the FILE's own z axis.

    grid.nz0 / grid.nzmax come from `nz` and `dz`, i.e. from the FULL
    model grid. That is only the same axis the file has when nothing has
    already cut it. Feed palm2gis a palm_postproc splitz output and the
    same `bounds.zmax` silently produced a different result - 19 m from
    OUTPUT_splitvar, 11 m from OUTPUT_splitz, with no warning - because
    indices derived from the full axis were applied to a truncated one.

    Recomputing from the file's own levels makes `bounds.zmax` mean the
    same thing either way. Falls back to the configured range when the
    file has no usable z axis.
    """
    label = "[{}] ".format(tag) if tag else ""
    if z_coord not in nc.variables:
        return grid.nz0, grid.nzmax

    z_file = np.asarray(nc.variables[z_coord][:], dtype=float)
    nz0 = int(max(np.searchsorted(z_file, grid.zt_min) - 1, 0))
    nzmax = int(min(np.searchsorted(z_file, grid.zt_max + grid.zmax_req) + 1,
                    z_file.size))

    pre_cut = getattr(nc, "palm_postproc_z_max", None)
    if pre_cut is not None:
        warning("{}this file was already cut to z <= {} m by "
                "palm_postproc splitz, so bounds.zmax is resolved against "
                "the {} levels that survived, not the full model grid.",
                label, pre_cut, z_file.size)

    if (nz0, nzmax) != (grid.nz0, grid.nzmax):
        debug("{}z range from the file's own axis: [{}:{}] "
              "(nz/dz would have given [{}:{}])",
              label, nz0, nzmax, grid.nz0, grid.nzmax)
    return nz0, nzmax


def require_variables(nc, names, tag, source="", hint=""):
    """Raise a readable error when *nc* lacks any of *names*.

    Without this a palm_postproc splitvar file - one variable per file by
    design - failed with a bare `KeyError: 'u'` from somewhere inside the
    stage, which says nothing about the cause or the fix.
    """
    missing = [n for n in names if n not in nc.variables]
    if not missing:
        return
    have = ", ".join(sorted(v for v in nc.variables
                            if nc.variables[v].dimensions
                            and nc.variables[v].dimensions[0] == "time"))
    raise RuntimeError(
        "[{}] this file does not contain {}. It holds: {}.{}\n"
        "  A palm_postproc splitvar/splitz/splittime file holds ONE "
        "variable (or one group) per file by design, so a stage that "
        "needs several has to read the joined file instead - or be "
        "configured for what this one actually has.{}".format(
            tag, ", ".join("'{}'".format(m) for m in missing),
            have or "(no time-dependent variables)",
            "\n  Given: {}".format(source) if source else "",
            "\n  {}".format(hint) if hint else ""))


def check_raw_input(nc, grid, tag, source=""):
    """Raise RuntimeError if *nc* looks like palm_postproc `coord` output.

    Three independent signals, any one of which rejects the file:
      a) x/y coordinate variable with standard_name projection_*_coordinate
         (written by coord when it overwrites local coords with UTM);
      b) x[0] far outside the local domain extent (absolute easting);
      c) a KELVIN-defined temperature variable arriving in degrees C.

    Surface-output files have no `x`/`y` dimension coordinates (elements
    are listed along `s` with xs/ys), so (a)/(b) are skipped when `x` is
    absent and (c) alone decides.

    A file stamped `palm_postproc_stage = "join"` short-circuits all
    three. join output is the intended feed for palm2gis, and a positive
    marker beats inferring safety from the absence of coord's
    fingerprints - which is what signal (c) used to do badly enough to
    nearly reject raw files containing PALM's own Celsius `ta`.
    """
    label = "[{}] ".format(tag) if tag else ""

    stage = str(getattr(nc, "palm_postproc_stage", "")).strip().lower()
    if stage == "join":
        debug("{}input is palm_postproc join output - raw-input guard "
              "satisfied by its stage marker.", label)
        return
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

    # c) - kelvin-defined temperatures that arrived in Celsius
    if reject is None:
        for name in nc.variables:
            if not name.lower().startswith(_COORD_EVIDENCE_PREFIXES):
                continue
            units = str(getattr(nc.variables[name], "units", "")).strip().lower()
            units = units.replace(" ", "")
            if units in _THERMO_CELSIUS_UNITS or (
                    units and units.startswith(_THERMO_CELSIUS_PREFIXES)):
                reject = ("variable '{}' has units '{}' but PALM writes it "
                          "in kelvin, i.e. it has already been converted"
                          .format(name, getattr(nc.variables[name], "units")))
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
