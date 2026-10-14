"""3D output -> gridded wind netCDF at ABSOLUTE heights, for ArcGIS Pro's
flow / vector-field renderers.

No paraview_vis counterpart - this replaces the ParaView Glyph/Calculator
workflow (make_vector(u,v,w) + Slice + Glyph): one netCDF per height above
sea level, (time, z, y, x) with u, v, w, speed and meteorological
direction, all collected in one folder. The renderer draws its own
streamlines or arrows from u / v (or speed / dir_met) and the time slider
steps through the records.
"""

import re
from pathlib import Path

import numpy as np
from netCDF4 import Dataset

import geopandas as gpd
from shapely.geometry import LineString

from ..descriptions import long_name, write_gpkg_descriptions
from ..log import debug, get_logger, info, stopwatch, warning
from .georef import (check_product, check_raw_input, require_variables,
                     resolve_z_range, select_times)
from .qgis_style import embed_gpkg_style, vector_qml, write_sidecar
from .voxel import destagger, write_crs_var

log = get_logger()


# ------------------------------
# 1. GRIDDED NETCDF PER ABSOLUTE HEIGHT (ARCGIS FLOW RENDERER)
# ------------------------------
def level_dir(pattern):
    """Folder that collects the per-height files: the pattern's own name
    without the extension and the {z} placeholder ('x_wind.nc' ->
    'x_wind/', 'x_wind_{z}m.nc' -> 'x_wind/'), next to where the pattern
    points."""
    p = Path(pattern)
    stem = p.name[:-len(p.suffix)] if p.suffix else p.name
    stem = re.sub(r"[_\-. ]*z?\{z\}m?", "", stem).strip("_-. ") or "wind"
    return p.parent / stem


def level_path(pattern, z_abs):
    """Output file for one absolute height, inside level_dir(pattern): {z}
    in the pattern is replaced by the height (e.g. 200, 200.5), otherwise
    `_z<height>m` goes before the extension."""
    tag = "{:g}".format(round(float(z_abs), 3))
    name = Path(pattern).name
    if "{z}" in name:
        name = name.replace("{z}", tag)
    else:
        stem, dot, ext = name.rpartition(".")
        name = ("{}_z{}m.{}".format(stem, tag, ext) if dot else
                "{}_z{}m".format(name, tag))
    return level_dir(pattern) / name


def _bracket(z_rel, target):
    """Levels k0, k1 and weight t with target = (1-t) z[k0] + t z[k1], or
    None when the target lies outside the file's levels."""
    if target < z_rel[0] - 1e-6 or target > z_rel[-1] + 1e-6:
        return None
    k1 = int(np.searchsorted(z_rel, target))
    k1 = min(max(k1, 1), z_rel.size - 1) if z_rel.size > 1 else 0
    k0 = max(k1 - 1, 0)
    if z_rel.size == 1 or abs(target - z_rel[k1]) < 1e-6:
        return k1, k1, 0.0
    if abs(target - z_rel[k0]) < 1e-6:
        return k0, k0, 0.0
    return k0, k1, float((target - z_rel[k0]) / (z_rel[k1] - z_rel[k0]))


def _solid(nc, grid, it, k, mask_src, z_rel):
    """True where PALM level k is inside terrain or a building: the PALM
    mask of a scalar at that level (as the voxel stage uses), or, without
    one, the level lying below the local obstacle top."""
    if mask_src is not None:
        return np.ma.getmaskarray(
            nc.variables[mask_src][it, k, grid.ny0:grid.nymax,
                                   grid.nx0:grid.nxmax])
    return z_rel[k] < grid.elev


def _heights(cfg, grid, nc, z_rel):
    """The absolute heights [m a.s.l.] both wind outputs sample:
    wind.netcdf.heights, or every PALM level of the voxel's vertical
    range - the same region.height_above_surface cut and the same k = 0
    rule - so they sit exactly on the voxel file's z values."""
    wd = cfg["wind"]
    if wd["grid_heights"]:
        return [float(h) for h in wd["grid_heights"]]
    nz0, nzmax = resolve_z_range(nc, grid, "wind")
    k_off = 1 if (cfg["voxel"].get("drop_k0", True) and nz0 == 0) else 0
    heights = [float(h) for h in grid.abs_z(z_rel[nz0 + k_off:nzmax])]
    info("[wind] every voxel level: {} heights, {:.1f} .. {:.1f} m",
         len(heights), heights[0], heights[-1])
    return heights


def _sample(nc, grid, it, br, solid):
    """u, v, w on the cut grid at one height (linear between the
    bracketing levels), NaN where `solid` (terrain or building)."""
    k0, k1, t = br
    out = {}
    for var in ("u", "v", "w"):
        a = destagger(nc, var, it, k0, k0 + 1, grid.ny0, grid.nymax,
                      grid.nx0, grid.nxmax)[0]
        if k1 != k0:
            b = destagger(nc, var, it, k1, k1 + 1, grid.ny0, grid.nymax,
                          grid.nx0, grid.nxmax)[0]
            a = (1.0 - t) * a + t * b
        out[var] = np.where(solid, np.nan, a)
    return out


def _wind_levels(cfg, grid, nc, its, times, z_rel, heights):
    """One netCDF per ABSOLUTE height (m a.s.l.): (time, z, y, x) with u,
    v, w, speed and dir_met, z of length 1 holding that height.

    For ArcGIS Pro's flow / vector-field renderers: Add Multidimensional
    Raster Layer, tick u and v (Vector Field U-V) or speed and dir_met
    (Magnitude-Direction). A requested height between two PALM levels is
    interpolated linearly between them, so the file holds the wind at
    exactly that height; a cell is NoData where either level is inside
    terrain or a building, or where the height is below the ground.
    """
    wd = cfg["wind"]
    stride = max(1, int(wd["grid_stride"]))
    og = grid.reproj.grid(stride)
    yd, xd = og.dims
    mask_src = next((c for c in ("theta", "ta", "q", "p", "rh")
                     if c in nc.variables), None)
    written = []

    for z_abs in heights:
        path = cfg.out_path(level_path(wd["grid_to"], z_abs))
        if not cfg.prepare_output(path, "wind"):
            continue
        sw = stopwatch()
        target = z_abs - grid.origin_z          # PALM's relative height
        br = _bracket(z_rel, target)
        if br is None:
            warning("[wind] {:g} m is outside the file's levels "
                    "({:.1f} .. {:.1f} m) - skipped.", z_abs,
                    float(grid.abs_z(z_rel[0])), float(grid.abs_z(z_rel[-1])))
            continue
        k0, k1, t = br
        solid = _solid(nc, grid, its[0], k0, mask_src, z_rel)
        if k1 != k0:
            solid = solid | _solid(nc, grid, its[0], k1, mask_src, z_rel)
        air = ~solid
        mask_out = og.mask(solid)
        debug("[wind] {:g} m a.s.l.: levels {} / {} (t = {:.3f}), {} of {} "
              "columns in terrain or buildings", z_abs, k0, k1, t,
              int(solid.sum()), solid.size)

        with cfg.atomic_output(path, "wind") as tmp, \
                Dataset(tmp, "w", format="NETCDF4") as fout:
            fout.Conventions = "CF-1.8"
            fout.title = ("PALM wind at {:g} m above sea level"
                          .format(z_abs))
            fout.createDimension("time", len(its))
            fout.createDimension("z", 1)
            gm = write_crs_var(fout, grid.reproj.crs,
                               cfg["crs"]["netcdf_attributes"])
            og.write_axes(fout)

            zv = fout.createVariable("z", "f8", ("z",))
            zv[:] = [z_abs]
            zv.standard_name = "height_above_mean_sea_level"
            zv.long_name = "height above sea level"
            zv.units = "m"
            zv.axis = "Z"
            zv.positive = "up"
            zv.comment = (
                "Absolute height (PALM zu + origin_z {:g} m). Interpolated "
                "linearly between PALM levels {:.3f} and {:.3f} m a.s.l."
                .format(grid.origin_z, float(grid.abs_z(z_rel[k0])),
                        float(grid.abs_z(z_rel[k1]))))

            tv = fout.createVariable("time", "f8", ("time",))
            tv.standard_name = "time"
            tv.units = "seconds since {}".format(
                grid.origin_time.strftime("%Y-%m-%d %H:%M:%S"))
            tv.calendar = "standard"
            tv.axis = "T"
            tv[:] = [float(times[it]) for it in its]

            specs = (
                ("u", "eastward_wind", "m s-1"),
                ("v", "northward_wind", "m s-1"),
                ("w", "upward_air_velocity", "m s-1"),
                ("speed", "wind_speed", "m s-1"),
                ("dir_met", "wind_from_direction", "degree"),
            )
            # NoData is NaN, not -9999: ArcGIS Pro's flow renderer reads
            # a -9999 fill in u/v as a real (-9999, -9999) m/s vector, so
            # every building cell became a huge wind towards the south-west
            # that swamped the real field and streaked through the walls.
            # A NaN can never be advected and is NoData everywhere else.
            for name, std, units in specs:
                ov = fout.createVariable(name, "f4", ("time", "z", yd, xd),
                                         fill_value=np.nan, zlib=True,
                                         complevel=4)
                ov.grid_mapping = gm
                ov.standard_name = std
                ov.units = units
                ov.long_name = long_name(name)
            fout.variables["dir_met"].comment = (
                "Meteorological convention: the direction the wind blows "
                "FROM, clockwise from north. In ArcGIS Pro's vector field "
                "renderer set Reference System = Geographic and Angle "
                "Representation = From. Horizontal only, ignoring w.")

            for q, it in enumerate(its):
                comp = {}
                for var in ("u", "v", "w"):
                    a = destagger(nc, var, it, k0, k0 + 1, grid.ny0,
                                  grid.nymax, grid.nx0, grid.nxmax)[0]
                    if k1 != k0:
                        b = destagger(nc, var, it, k1, k1 + 1, grid.ny0,
                                      grid.nymax, grid.nx0, grid.nxmax)[0]
                        a = (1.0 - t) * a + t * b
                    a = np.where(air, a, np.nan)
                    # Onto the output grid (identity in PALM's CRS); the
                    # conservative remap averages over air cells only.
                    comp[var] = og.take(a, np.nan, valid=air)
                u, v = og.rotate(comp["u"], comp["v"])
                sp = np.hypot(u, v)
                dr = (270.0 - np.degrees(np.arctan2(v, u))) % 360.0
                for name, val in (("u", u), ("v", v), ("w", comp["w"]),
                                  ("speed", sp), ("dir_met", dr)):
                    arr = np.array(val, dtype="f4")
                    arr[mask_out | ~np.isfinite(arr)] = np.nan
                    fout.variables[name][q, 0, ...] = arr

        written.append(path)
        debug("[wind] {}: {} x {} x {} (t,{},{}) at {:g} m a.s.l., u, v, w, "
              "speed, dir_met", cfg.disp(path), len(its), og.shape[0],
              og.shape[1], yd, xd, z_abs)
        info("[wind] wrote {} in {}", cfg.disp(path), sw)
    if written:
        debug("[wind] ArcGIS Pro: Add Multidimensional Raster Layer as a "
              "Vector Field")
    return written


# ------------------------------
# 1b. 3D WIND ARROWS (GEOPACKAGE)
# ------------------------------
def _wind_arrows(cfg, grid, nc, its, times, z_rel, heights, path):
    """3D arrows (LineStringZ) at the same absolute heights as the netCDF.

    One arrow per kept grid column, height and record: it starts at the
    cell centre at that height and points along (u, v, w), scaled by
    wind.gpkg.arrow_length_per_ms. Built in PALM's CRS and transformed as
    a whole into the output CRS (exact, z untouched), so it keeps its
    metre length and points where the wind blows; the u / v attributes
    are rotated into the output CRS's east / north. Columns inside
    terrain or a building at that height get no arrow. Fields: time, z
    (m a.s.l.), u, v, w, speed_3d. Written one record at a time.
    """
    wd = cfg["wind"]
    sw = stopwatch()
    stride = max(1, int(wd["stride"]))
    scale = float(wd["arrow_scale"])
    mask_src = next((c for c in ("theta", "ta", "q", "p", "rh")
                     if c in nc.variables), None)
    jj = np.arange(0, grid.nymax - grid.ny0, stride)
    ii = np.arange(0, grid.nxmax - grid.nx0, stride)
    jg, ig = np.meshgrid(jj, ii, indexing="ij")
    jf, if_ = jg.ravel(), ig.ravel()
    x_abs = grid.abs_x(grid.x_stat)[if_]
    y_abs = grid.abs_y(grid.y_stat)[jf]

    levels = []
    for z_abs in heights:
        br = _bracket(z_rel, z_abs - grid.origin_z)
        if br is None:
            continue
        solid = _solid(nc, grid, its[0], br[0], mask_src, z_rel)
        if br[1] != br[0]:
            solid = solid | _solid(nc, grid, its[0], br[1], mask_src, z_rel)
        keep = ~solid[jf, if_]
        if keep.any():
            levels.append((z_abs, br, solid, keep))
    n_rows = sum(int(k.sum()) for *_x, k in levels) * len(its)
    if not levels:
        warning("[wind] no air cells at the wind heights - {} skipped.",
                cfg.disp(path))
        return
    if n_rows > 2_000_000:
        warning("[wind] {} will hold {} arrows - large; raise "
                "every_nth_column.", cfg.disp(path), n_rows)

    dts = grid.times_to_datetime(times)
    dts_tz = grid.times_to_datetime(times, tz_aware=True)
    if path.exists():
        path.unlink()           # never append onto a previous run
    for n, it in enumerate(its):
        parts = []
        for z_abs, br, solid, keep in levels:
            c = _sample(nc, grid, it, br, solid)
            j, i = jf[keep], if_[keep]
            ug, vg, wv = c["u"][j, i], c["v"][j, i], c["w"][j, i]
            u, v = grid.reproj.rotate_at(ug, vg, j, i)
            x0, y0 = x_abs[keep], y_abs[keep]
            z0 = np.full(x0.size, z_abs)
            segs = [LineString([(a, b, h), (a + du * scale, b + dv * scale,
                                            h + dw * scale)])
                    for a, b, h, du, dv, dw
                    in zip(x0, y0, z0, np.nan_to_num(ug), np.nan_to_num(vg),
                           np.nan_to_num(wv))]
            parts.append(gpd.GeoDataFrame(
                {"time": [dts_tz[it]] * x0.size,
                 "z": z0.astype("f4"),
                 "u": u.astype("f4"), "v": v.astype("f4"),
                 "w": wv.astype("f4"),
                 "speed_3d": np.sqrt(u * u + v * v + wv * wv).astype("f4")},
                geometry=grid.reproj.geoms(segs), crs=grid.reproj.crs))
        frame = (parts[0] if len(parts) == 1
                 else gpd.pd.concat(parts, ignore_index=True))
        frame.to_file(path, layer="wind", driver="GPKG",
                      mode="w" if n == 0 else "a")
    write_gpkg_descriptions(path, "wind")
    if cfg["qgis_styles"]:
        from .surfaces import _median_step
        qml = vector_qml("time", _median_step(dts, its), geom="line",
                         clamp_absolute=bool(cfg["qgis_styles_3d"]))
        write_sidecar(path, qml)
        embed_gpkg_style(path, "wind", qml)
    debug("[wind] {}: {} heights, every {} column(s), {} m per m/s",
          cfg.disp(path), len(levels), stride, scale)
    info("[wind] wrote {}: {} arrows in {}", cfg.disp(path), n_rows, sw)


# ------------------------------
# 2. STAGE ENTRY
# ------------------------------
def stage_wind(cfg, grid):
    wd = cfg["wind"]
    # One netCDF per absolute height; each is checked for existing output
    # on its own inside _wind_levels.
    arrows = cfg.out_path(wd["file_to"]) if wd["file_to"] else None
    want_arrows = arrows is not None and cfg.prepare_output(arrows, "wind")
    if not (wd["grid_to"] or want_arrows):
        return

    with Dataset(cfg.resolve(wd["file_from"]), "r") as nc:
        check_product(nc, "wind", wd["file_from"])
        check_raw_input(nc, grid, "wind", wd["file_from"])
        its, times = select_times(cfg, nc, "wind", grid, wd["file_from"])

        require_variables(
            nc, ["u", "v"], "wind", wd["file_from"],
            hint="wind needs u and v together, which among the splitvar "
                 "outputs only the grouped u/v/w file has "
                 "(steps.splitvar.group_vars).")
        z = np.asarray(nc.variables["zu_3d"][:], dtype="f8")
        heights = _heights(cfg, grid, nc, z)
        if wd["grid_to"]:
            _wind_levels(cfg, grid, nc, its, times, z, heights)
        if want_arrows:
            _wind_arrows(cfg, grid, nc, its, times, z, heights, arrows)
