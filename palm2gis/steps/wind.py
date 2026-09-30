"""3D output -> wind layers: points / 3D arrows (GeoPackage) at heights
above the local surface, and gridded netCDF at ABSOLUTE heights.

No paraview_vis counterpart - this replaces the ParaView Glyph/Calculator
workflow (make_vector(u,v,w) + Slice + Glyph) with a GIS-native product:
a point per grid column (optionally thinned by wind.stride) at one or more
heights ABOVE THE LOCAL SURFACE (terrain + buildings), carrying u, v,
speed and meteorological direction per timestep. In ArcGIS Pro / QGIS,
symbolise with arrow markers rotated by `dir_met` (geographic rotation)
and sized/coloured by `speed`; the `time` field drives the time slider.

Height-above-ground sampling follows terrain and rooftops, which is what
pedestrian-comfort maps actually need (a fixed-z ParaView slice cuts
through hillsides and buildings instead).
"""

import numpy as np
from netCDF4 import Dataset

import geopandas as gpd
from shapely.geometry import LineString, Point

from ..descriptions import long_name, write_gpkg_descriptions
from ..log import debug, get_logger, info, warning
from .georef import (check_product, check_raw_input, require_variables,
                     resolve_z_range, select_times)
from .voxel import FILL, destagger, write_crs_var
from .qgis_style import (embed_gpkg_style, vector_qml, write_sidecar)

log = get_logger()


# ------------------------------
# 1. GRIDDED NETCDF PER ABSOLUTE HEIGHT (ARCGIS FLOW RENDERER)
# ------------------------------
def _warn_above_top(inside, z):
    """Say how many samples fell above the file's top level (-> NoData).

    A vertically cut input (palm_postproc splitz) or a tall building with
    a large height above it puts surface + h above the last level. The
    index used to be clipped to that last level, so those samples carried
    wind from a different height under the requested label.
    """
    for h, ok in inside.items():
        n_out = int((~ok).sum())
        if n_out:
            warning("[wind] {} of {} columns: {:g} m AGL is above the top "
                    "level ({:g} m) - NoData.", n_out, ok.size, h,
                    float(z[-1]))



def level_path(pattern, z_abs):
    """Output file for one absolute height: {z} in the pattern is replaced
    by the height (e.g. 200, 200.5), otherwise `_z<height>m` goes before
    the extension."""
    tag = "{:g}".format(round(float(z_abs), 3))
    if "{z}" in pattern:
        return pattern.replace("{z}", tag)
    stem, dot, ext = pattern.rpartition(".")
    return "{}_z{}m.{}".format(stem, tag, ext) if dot else \
        "{}_z{}m".format(pattern, tag)


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


def _wind_levels(cfg, grid, nc, its, times, z_rel):
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

    if wd["grid_heights"]:
        heights = [float(h) for h in wd["grid_heights"]]
    else:
        # No list: every PALM level of the voxel's vertical range - the
        # same region.height_above_surface cut and the same k = 0 rule -
        # so these files sit exactly on the voxel file's z values.
        nz0, nzmax = resolve_z_range(nc, grid, "wind")
        k_off = 1 if (cfg["voxel"].get("drop_k0", True) and nz0 == 0) else 0
        heights = [float(h) for h in grid.abs_z(z_rel[nz0 + k_off:nzmax])]
        info("[wind] no wind.netcdf.heights - one file per PALM level of "
             "the voxel range: {} levels, {:.2f} .. {:.2f} m a.s.l.",
             len(heights), heights[0], heights[-1])

    for z_abs in heights:
        path = cfg.out_path(level_path(wd["grid_to"], z_abs))
        if not cfg.prepare_output(path, "wind"):
            continue
        target = z_abs - grid.origin_z          # PALM's relative height
        br = _bracket(z_rel, target)
        if br is None:
            warning("[wind] {:g} m a.s.l. is outside the file's levels "
                    "({:.2f} .. {:.2f} m a.s.l.) - {} skipped.", z_abs,
                    float(grid.abs_z(z_rel[0])), float(grid.abs_z(z_rel[-1])),
                    path.name)
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
            for name, std, units in specs:
                ov = fout.createVariable(name, "f4", ("time", "z", yd, xd),
                                         fill_value=FILL, zlib=True,
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
                    arr[mask_out | ~np.isfinite(arr)] = FILL
                    fout.variables[name][q, 0, ...] = arr

        written.append(path)
        info("[wind] wrote {}: {} x {} x {} (t,{},{}) at {:g} m a.s.l., "
             "u, v, w, speed, dir_met", cfg.disp(path), len(its),
             og.shape[0], og.shape[1], yd, xd, z_abs)
    if written:
        debug("[wind] ArcGIS Pro: Add Multidimensional Raster Layer as a "
              "Vector Field")
    return written


# ------------------------------
# 2. STAGE ENTRY
# ------------------------------
def stage_wind(cfg, grid):
    wd = cfg["wind"]
    file_to = cfg.out_path(wd["file_to"]) if wd["file_to"] else None
    want_gpkg = file_to is not None and cfg.prepare_output(file_to, "wind")
    # One netCDF per absolute height; each is checked for existing output
    # on its own inside _wind_levels.
    want_grid = bool(wd["grid_to"])
    if not (want_gpkg or want_grid):
        return

    with Dataset(cfg.resolve(wd["file_from"]), "r") as nc:
        check_product(nc, "wind", wd["file_from"])
        check_raw_input(nc, grid, "wind", wd["file_from"])
        its, times = select_times(cfg, nc, "wind", grid, cfg["wind"]["file_from"])

        require_variables(
            nc, ["u", "v"], "wind", wd["file_from"],
            hint="wind needs u and v together, which among the splitvar "
                 "outputs only the grouped u/v/w file has "
                 "(steps.splitvar.group_vars).")
        z = np.asarray(nc.variables["zu_3d"][:], dtype="f8")
        nz = z.size
        if want_grid:
            _wind_levels(cfg, grid, nc, its, times, z)
        if not want_gpkg:
            return
        stride = int(wd["stride"])
        heights = [float(h) for h in wd["heights_agl"]]

        # thinned column selection (indices into the cut subdomain)
        jj = np.arange(0, grid.nymax - grid.ny0, stride)
        ii = np.arange(0, grid.nxmax - grid.nx0, stride)
        jg, ig = np.meshgrid(jj, ii, indexing="ij")
        jf, if_ = jg.ravel(), ig.ravel()
        debug("[wind] {} columns (stride {}), heights AGL: {} m",
              len(jf), stride, ", ".join(str(h) for h in heights))


        lines_3d = bool(wd["write_lines_3d"])
        include_w = bool(wd["include_w"])
        scale = float(wd["arrow_scale"])

        # per column: k index of the first level >= local surface + h_agl
        k_idx = {}
        z_abs = {}
        inside = {}
        for h in heights:
            k = np.searchsorted(z, grid.elev[jf, if_] + h)
            inside[h] = k < nz
            k_idx[h] = np.clip(k, 0, nz - 1)
            # The true height of the sample, for the 3D geometry. The 2D
            # point layer carried this only as the `height_agl` attribute,
            # so every level collapsed onto the ground in a 3D scene and
            # the 2 m and 10 m samples sat exactly on top of each other.
            # grid.elev is terrain + buildings ABOVE THE DOMAIN BASE, so
            # the sample height needs origin_z to become a real elevation.
            # Without it the 3D arrows sat origin_z metres underground.
            z_abs[h] = grid.abs_z(grid.elev[jf, if_] + h)
        _warn_above_top(inside, z)

        # Geometry is built in PALM's CRS and transformed as a whole into
        # the output CRS (exact, z untouched), so the 3D arrows keep their
        # metre-based length and point where the wind actually blows.
        reproj = grid.reproj
        x_abs = grid.abs_x(grid.x_stat)
        y_abs = grid.abs_y(grid.y_stat)
        pts = reproj.geoms([Point(x_abs[i], y_abs[j])
                            for j, i in zip(jf, if_)])
        dts = grid.times_to_datetime(times)
        dts_tz = grid.times_to_datetime(times, tz_aware=True)

        # Written one timestep at a time (first with mode="w", the rest
        # appended) so peak memory is a single timestep instead of
        # points x heights x timesteps held at once.
        if file_to.exists():
            file_to.unlink()      # never append onto a previous run

        n_out = 0
        for n, it in enumerate(its):
            debug("[wind] time {} ({:.0f} s)", it, times[it])
            u3 = destagger(nc, "u", it, 0, nz, grid.ny0, grid.nymax,
                           grid.nx0, grid.nxmax)
            v3 = destagger(nc, "v", it, 0, nz, grid.ny0, grid.nymax,
                           grid.nx0, grid.nxmax)
            w3 = (destagger(nc, "w", it, 0, nz, grid.ny0, grid.nymax,
                            grid.nx0, grid.nxmax)
                  if (include_w and lines_3d) else None)
            parts = []
            parts3 = []
            for h in heights:
                k = k_idx[h]
                u_g = np.where(inside[h], u3[k, jf, if_], np.nan)
                v_g = np.where(inside[h], v3[k, jf, if_], np.nan)
                # Attributes in the output CRS's east/north; the arrow
                # geometry below is built from the grid components and
                # transformed, which gives the same direction.
                u, v = reproj.rotate_at(u_g, v_g, jf, if_)
                speed = np.hypot(u, v)
                # `speed` and `dir_met` stay HORIZONTAL even when the 3D
                # geometry tilts with w, so the two layers remain directly
                # comparable and dir_met keeps its meteorological meaning.
                dir_met = (270.0 - np.degrees(np.arctan2(v, u))) % 360.0
                attrs = {"time": [dts_tz[it]] * len(pts),
                         "rec": it,
                         "height_agl": h,
                         "u": u.astype(np.float32),
                         "v": v.astype(np.float32),
                         "speed": speed.astype(np.float32),
                         "dir_met": dir_met.astype(np.float32)}
                parts.append(gpd.GeoDataFrame(dict(attrs), geometry=pts,
                                              crs=reproj.crs))
                if not lines_3d:
                    continue
                w = (np.where(inside[h], w3[k, jf, if_], np.nan)
                     if w3 is not None else np.zeros_like(u))
                zb = z_abs[h]
                segs = [LineString([(x0, y0, z0),
                                    (x0 + ux * scale, y0 + vy * scale,
                                     z0 + wz * scale)])
                        for x0, y0, z0, ux, vy, wz
                        in zip(x_abs[if_], y_abs[jf], zb,
                               np.nan_to_num(u_g), np.nan_to_num(v_g),
                               np.nan_to_num(w))]
                segs = reproj.geoms(segs)
                a3 = dict(attrs)
                a3["z"] = zb.astype(np.float32)
                if w3 is not None:
                    a3["w"] = w.astype(np.float32)
                parts3.append(gpd.GeoDataFrame(a3, geometry=segs,
                                               crs=reproj.crs))
            frame = (parts[0] if len(parts) == 1
                     else gpd.pd.concat(parts, ignore_index=True))
            frame.to_file(file_to, layer="wind", driver="GPKG",
                          mode="w" if n == 0 else "a")
            n_out += len(frame)

            if lines_3d:
                f3 = (parts3[0] if len(parts3) == 1
                      else gpd.pd.concat(parts3, ignore_index=True))
                # Always appended: the point layer above has already
                # created (or replaced) the file for this timestep.
                f3.to_file(file_to, layer="wind_3d", driver="GPKG",
                           mode="a")

        write_gpkg_descriptions(file_to, "wind")
        if lines_3d:
            write_gpkg_descriptions(file_to, "wind_3d")

        if cfg["qgis_styles"]:
            from .surfaces import _median_step
            step = _median_step(dts, its)
            qml = vector_qml("time", step)
            write_sidecar(file_to, qml)
            embed_gpkg_style(file_to, "wind", qml)
            if lines_3d:
                embed_gpkg_style(file_to, "wind_3d",
                                 vector_qml("time", step, geom="line",
                                            clamp_absolute=True))
        info("[wind] wrote {}: {} pts x {} heights x {} times "
             "(stride {})", cfg.disp(file_to), len(pts), len(heights),
             len(its), stride)
        if lines_3d:
            info("[wind] wrote {}: layer 'wind_3d' with {} 3D arrows "
                 "({} m per m/s, {})", cfg.disp(file_to),
                 len(pts) * len(heights) * len(its), scale,
                 "u, v, w" if include_w else "u, v (horizontal)")
