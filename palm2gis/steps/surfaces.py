"""PALM surface output -> GIS layers.

GIS counterpart of paraview_vis/utils/surf_gen.py. Instead of building VTU
quads for ParaView, the unstructured surface elements (xs, ys, zs with
azimuth/zenith) become 3D polygons:

  * Every ground, roof and wall element is a quad (PolygonZ) in absolute
    coordinates, one row per element per timestep (elem_id, orient,
    azimuth, z, z_agl, time, one column per variable), so the ArcGIS time
    slider / QGIS temporal controller can animate it and the quads render
    correctly in 3D scenes.

  * The same elements are written as a GeoPackage (QGIS), an ESRI
    MultiPatch shapefile and a File Geodatabase (ArcGIS), each optional.

Elements with other zenith values (sloped roofs, downward-facing) are
counted and reported but not converted.
"""

import numpy as np
from netCDF4 import Dataset

import geopandas as gpd
import pyogrio
from shapely.geometry import MultiPolygon, Polygon

from ..descriptions import write_gpkg_descriptions, write_with_aliases
from ..log import debug, get_logger, info, warning
from .georef import (check_product, check_raw_input, require_variables,
                     select_times)
from .thermo import convert_temperature, exner_simple, is_celsius_exempt, \
    is_temperature, temperature_scale, theta_to_ta
from .qgis_style import embed_gpkg_style, vector_qml, write_sidecar

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
    temperature. The same list drives every element layer, so a field is never present in one product and missing
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
    # names (whether `t_surf` is copied or derived), and that has to be
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
    mp_path = cfg.out_path(sf["multipatch_to"]
                           or "{}surfaces.shp".format(sf["prefix"]))
    do_mp = (sf["write_multipatch"]
             and cfg.prepare_output(mp_path, "surfaces"))
    fgdb_path = cfg.out_path(sf["fgdb_to"]
                             or "{}surfaces.gdb".format(sf["prefix"]))
    do_fgdb = (sf["write_fgdb"]
               and cfg.prepare_output(fgdb_path, "surfaces"))
    if not (do_all or do_mp or do_fgdb):
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
