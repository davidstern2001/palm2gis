"""PALM surface output -> GIS layers.

GIS counterpart of paraview_vis/utils/surf_gen.py. Instead of building VTU
quads for ParaView, the unstructured surface elements (xs, ys, zs with
azimuth/zenith) become 3D polygons:

  * Every ground, roof and wall element is a quad (PolygonZ) in absolute
    coordinates, one row per element per timestep (elem_id, orient,
    azimuth, z, z_agl, time, one column per variable), so the ArcGIS time
    slider / QGIS temporal controller can animate it and the quads render
    correctly in 3D scenes.

  * The same elements are written as a GeoPackage (QGIS) and a File
    Geodatabase (ArcGIS), each optional.

Elements with other zenith values (sloped roofs, downward-facing) are
counted and reported but not converted.
"""

import numpy as np
from netCDF4 import Dataset

import geopandas as gpd
from shapely.geometry import Polygon

from ..descriptions import write_gpkg_descriptions, write_with_aliases
from ..log import debug, get_logger, info, warning
from .georef import (check_product, check_raw_input, require_variables,
                     select_times)
from .thermo import convert_temperature, exner_basic_state, \
    is_celsius_exempt, is_temperature, temperature_scale, theta_to_ta
from .qgis_style import embed_gpkg_style, vector_qml, write_sidecar

log = get_logger()

FILL = -9999.0


# ------------------------------
# 1. LOAD AND CLASSIFY ELEMENTS
# ------------------------------
def _load_geometry(nc):
    # float64: PALM stores xs/ys/zs as float32, and float32 + origin stays
    # float32 - a UTM northing of 5.5e6 then resolves only 0.5 m, so every
    # quad corner was snapped to the nearest half metre.
    xs = np.asarray(nc.variables["xs"][:], dtype="f8")
    ys = np.asarray(nc.variables["ys"][:], dtype="f8")
    zs = np.asarray(nc.variables["zs"][:], dtype="f8")
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
# `t_surface` is the actual surface (skin) temperature, named to match
# `theta_surface`, the surface POTENTIAL temperature. Both spellings of the
# potential form occur depending on PALM version and output block, so both
# are accepted; a file's own actual surface temperature is recognised under
# either `t_surface` or the older `t_surf`, and is always WRITTEN as
# `t_surface` (palm2gis wrote `t_surf` up to 0.21).
SURF_ACTUAL = "t_surface"
SURF_ACTUAL_NATIVE = ("t_surface", "t_surf")
SURF_POTENTIAL = ("theta_surface", "theta_surf")


def output_vars(sf, available=None, verbose=False):
    """(output name, source variable, derive?) for every surface product.

    With surfaces.derive_t_surface the actual SURFACE temperature
    `t_surface` is ENSURED: if the file already holds it (as `t_surface` or
    `t_surf`), PALM's own value is used as-is; only when it is absent is it
    derived from the surface potential temperature. The same list drives
    every element layer, so a field is never present in one product and
    missing from another.

    `available` is the set of variable names in the source file. Listed
    variables that are missing but present under the other accepted
    spelling are silently substituted.
    """
    have = set(available or ())
    out, names = [], []

    for var in sf["variables"]:
        src = var
        if have and var not in have:
            for group in (SURF_POTENTIAL, SURF_ACTUAL_NATIVE):
                alias = next((a for a in group if a != var and a in have),
                             None)
                if var in group and alias:
                    if verbose:
                        warning("[surfaces] '{}' not in the file; using "
                                "'{}'.", var, alias)
                    src = alias
                    break
        out.append((var, src, False))
        names.append(var)

    if not sf["derive_t_surface"]:
        return out

    if SURF_ACTUAL in names:
        return out                       # user asked for it explicitly
    native = next((a for a in SURF_ACTUAL_NATIVE if a in have), None)
    if native:
        if verbose:
            info("[surfaces] the file holds '{}' - using PALM's own "
                 "surface temperature instead of deriving it.", native)
        out.append((SURF_ACTUAL, native, False))
        return out

    src = next((v for _n, v, _d in out if v in SURF_POTENTIAL), None)
    if src is None:
        if verbose:
            warning("[surfaces] derive_t_surface is on but neither '{}' nor "
                    "a surface potential temperature ({}) is available.",
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
    exner_basic_state pairs it with cfg["surface_pressure"], defined as the
    pressure at the model's z = 0. Passing absolute (origin_z-shifted)
    heights here would put a second, larger error on every derived
    temperature.
    """
    out = vals
    if derive:
        good = np.isfinite(vals) & (vals != FILL)
        if good.any():
            # PALM's own basic state: PALM made theta_surface from t_surface
            # with exactly this Exner function, so this inverts it.
            pi = exner_basic_state(z_elem, cfg["surface_pressure"],
                                   cfg["pt_surface"])
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


def _wall_quad(x, y, z, az, grid):
    """Backwards-compatible helper: a wall element is zenith 90."""
    return _element_quad(x, y, z, 90.0, az, grid)


# ------------------------------
# 3b. FILE GEODATABASE
# ------------------------------
def _element_frames(cfg, grid, nc, its, times, xs, ys, zs, zen, az, idx):
    """Yield (timestep index, GeoDataFrame) for the selected elements.

    Shared by the GeoPackage and File Geodatabase writers so the geometry
    and the attributes cannot drift apart between products.
    """
    sf = cfg["surfaces"]
    # Quads are built in PALM's CRS and transformed as a whole into the
    # output CRS (project.output_crs): exact per vertex, z untouched.
    geoms = list(grid.reproj.geoms(
        [_element_quad(xs[s], ys[s], zs[s], zen[s], az[s], grid)
         for s in idx]))
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
               "rec": np.full(idx.size, int(it), dtype="i4"),
               "time": [dts_tz[it]] * idx.size,
               "azimuth": np.where(np.isclose(az[idx], FILL), np.nan,
                                   az[idx])}
        for name, var, derive in output_vars(sf, nc.variables):
            vals = np.ma.filled(nc.variables[var][it, :], np.nan)[idx]
            rec[name] = _derive_values(
                cfg, vals, zs[idx], derive, name=name,
                src_units=getattr(nc.variables[var], "units", None))
        yield n, it, gpd.GeoDataFrame(rec, geometry=geoms,
                                      crs=grid.reproj.crs)


def _elements_fgdb(cfg, grid, nc, its, times, xs, ys, zs, zen, az, idx,
                   path, layer):
    """Write the same elements as PolygonZ into a File Geodatabase.

    The GeoPackage keeps QGIS happy; this is the container ArcGIS opens
    natively, with real datetimes, full-length field names and NULLs, so
    nothing is lost relative to the GeoPackage.
    """
    if idx.size == 0:
        warning("[surfaces] no elements in bounds - {} skipped.",
                cfg.disp(path))
        return
    for n, it, frame in _element_frames(cfg, grid, nc, its, times, xs, ys,
                                        zs, zen, az, idx):
        write_with_aliases(frame, path, "OpenFileGDB", layer,
                           append=(n > 0))
    info("[surfaces] wrote {}: feature class '{}', {} elements x {} times",
         cfg.disp(path), layer, idx.size, len(its))


def _elements_gpkg(cfg, grid, nc, its, times, xs, ys, zs, zen, az, sel,
                   path, idx, kind):
    """Write one GeoPackage row per surface element per timestep.

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
        warning("[surfaces] {} will hold {} rows - slow and large; tighten "
                "bounds/n_times or disable it.",
                cfg.disp(path), n_rows)

    geoms = list(grid.reproj.geoms(
        [_element_quad(xs[s], ys[s], zs[s], zen[s], az[s], grid)
         for s in idx]))
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
        frame = gpd.GeoDataFrame(rec, geometry=geoms, crs=grid.reproj.crs)
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
        n_written += idx.size
        debug("[surfaces] {}: timestep {}/{} written ({} rows so far)",
              cfg.disp(path), n + 1, len(its), n_written)

    write_gpkg_descriptions(path, layer)

    qml = None
    if cfg["qgis_styles"]:
        step = _median_step(dts, its)
        # The polygon layer carries true Z, so it gets the 3D renderer.
        qml = vector_qml("time", step, geom="polygon", attr=ramp_attr,
                         vmin=v_lo, vmax=v_hi,
                         clamp_absolute=bool(cfg["qgis_styles_3d"]))
        write_sidecar(path, qml)
        embed_gpkg_style(path, layer, qml)

    counts = gpd.pd.Series(orient).value_counts().to_dict()
    info("[surfaces] wrote {}: {} elements x {} times ({})",
         cfg.disp(path), idx.size, len(its),
         ", ".join(f"{k} {v}" for k, v in sorted(counts.items())))


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
    # names (whether `t_surface` is copied or derived), and that has to be
    # settled before the overwrite check picks the files to write.
    with Dataset(src_path, "r") as nc:
        check_product(nc, "surfaces", sf["file_from"])
        check_raw_input(nc, grid, "surfaces", sf["file_from"])
        available = set(nc.variables)
        plan = output_vars(sf, available, verbose=True)
        # Checked up front: a listed variable the file lacks used to fail
        # halfway through with a bare KeyError naming only the variable.
        require_variables(nc, sorted({v for _n, v, _d in plan}), "surfaces",
                          sf["file_from"],
                          hint="Adjust surfaces.variables to what this "
                               "file holds.")
    derived = [n for n, _v, d in plan if d]
    if derived:
        info("[surfaces] deriving {} from the surface potential "
             "temperature.", ", ".join(derived))

    all_path = cfg.out_path(sf["all_elements_to"]
                            or "{}surfaces.gpkg".format(sf["prefix"]))
    do_all = (sf["write_all_elements_gpkg"]
              and cfg.prepare_output(all_path, "surfaces"))
    fgdb_path = cfg.out_path(sf["fgdb_to"]
                             or "{}surfaces.gdb".format(sf["prefix"]))
    do_fgdb = (sf["write_fgdb"]
               and cfg.prepare_output(fgdb_path, "surfaces"))
    if not (do_all or do_fgdb):
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

        if do_all:
            _all_elements_gpkg(cfg, grid, nc, its, times, xs, ys, zs, zen,
                               az, sel, all_path)
        if do_fgdb:
            idx_all = np.where(sel)[0]
            _elements_fgdb(cfg, grid, nc, its, times, xs, ys, zs, zen,
                           az, idx_all, fgdb_path, "surfaces")
            wall_idx = np.where(sel & (np.abs(zen - 90.0) < 1.0e-3))[0]
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
