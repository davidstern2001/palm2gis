"""PROJ database sanity check and repair.

A Python environment that mixes conda and pip packages can end up with
GDAL/rasterio linked against one PROJ but resolving its database from
another, which shows up as:

    PROJ: internal_proj_create_from_database: .../proj.db contains
    DATABASE.LAYOUT.VERSION.MINOR = 4 whereas a number >= 6 is expected.
    It comes from another PROJ installation.

The visible damage is that GDAL cannot encode a CRS, so rasters are
written with no georeferencing while everything else looks fine. This
module detects that at startup and repoints GDAL's PROJ search path at a
database that actually matches - rasterio's own bundled copy, or pyproj's
- for the duration of the run. It changes nothing on disk and nothing
outside this process.

It is a workaround, not a cure: the run says what it did and what to fix.
"""

import os

from .log import debug, info, warning

# EPSG code used to probe: resolving it requires a working proj.db.
_PROBE_EPSG = 32633


# ------------------------------
# 1. PROBE
# ------------------------------
def _proj_works():
    """True if GDAL can resolve a CRS from the EPSG database."""
    try:
        import rasterio.crs
        crs = rasterio.crs.CRS.from_epsg(_PROBE_EPSG)
        return bool(crs and crs.to_wkt())
    except Exception:
        return False


def _candidate_dirs():
    """Directories that might hold a usable proj.db, best first."""
    out = []
    try:
        import rasterio
        out.append(os.path.join(os.path.dirname(rasterio.__file__),
                                "proj_data"))
    except Exception:
        pass
    try:
        import pyproj.datadir
        out.append(pyproj.datadir.get_data_dir())
    except Exception:
        pass
    for var in ("PROJ_DATA", "PROJ_LIB"):
        if os.environ.get(var):
            out.append(os.environ[var])
    seen, uniq = set(), []
    for d in out:
        if d and d not in seen and os.path.isfile(os.path.join(d,
                                                               "proj.db")):
            seen.add(d)
            uniq.append(d)
    return uniq


# ------------------------------
# 2. REPAIR
# ------------------------------
def ensure_proj_data():
    """Check PROJ, and repoint it at a working database if needed.

    Returns True when CRS handling works (either it already did, or the
    repair succeeded).
    """
    if _proj_works():
        return True

    warning("[proj] PROJ cannot resolve EPSG codes (proj.db from another "
            "installation) - trying to repoint it.")

    try:
        from rasterio._env import (get_proj_data_search_paths,
                                   set_proj_data_search_path)
    except ImportError:
        warning("[proj] this rasterio cannot set PROJ search paths - "
                "rasters may be written without a CRS.")
        return False

    original = None
    try:
        original = get_proj_data_search_paths()
    except Exception:
        pass
    debug("[proj] search paths before repair: {}", original)

    for path in _candidate_dirs():
        try:
            set_proj_data_search_path(path)
        except Exception:
            continue
        if _proj_works():
            info("[proj] repointed at {} for this run.", path)
            warning("[proj] to fix it for good, install GDAL/rasterio/pyproj "
                    "from ONE source and clear stale PROJ_DATA / PROJ_LIB.")
            return True

    if original:
        try:
            set_proj_data_search_path(original[0])
        except Exception:
            pass
    warning("[proj] no usable proj.db - rasters will have NO CRS; install "
            "GDAL/rasterio/pyproj from one source.")
    return False
