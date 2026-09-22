"""3D output -> wind vector point layer (GeoPackage).

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
                     select_times)
from .voxel import FILL, destagger, write_crs_var
from .qgis_style import (embed_gpkg_style, vector_qml, write_sidecar)

log = get_logger()


# ------------------------------
# 1. GRIDDED NETCDF (ARCGIS FLOW RENDERER)
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



def _wind_netcdf(cfg, grid, nc, its, times, heights, z, path):
    """Regular (time, height, y, x) grid of u, v, w, speed and dir_met.

    This is the product ArcGIS Pro's flow / vector-field renderers want:
    Add Multidimensional Raster Layer, tick two variables, and choose
    Output Configuration "Vector Field (U-V)" or "(Magnitude-Direction)".
    Both pairings are written so neither requires regenerating the file.

    `height` is the sampling level ABOVE THE LOCAL SURFACE, terrain and
    rooftops included - the same terrain-following sampling the point
    layer uses, and what pedestrian-comfort work needs. It is therefore
    NOT an elevation, and deliberately carries no `standard_name`:
    claiming CF's `height` would assert a single elevation per level,
    which terrain-following sampling does not have.

    No `coordinates` attribute is written. That attribute declares
    AUXILIARY coordinates in CF, and Esri accepts only netCDF with no
    auxiliary variables as a multidimensional source.
    """
    wd = cfg["wind"]
    stride = max(1, int(wd["grid_stride"]))
    nz = z.size

    # Own column selection: `wind.stride` thins the GeoPackage point
    # layer, `wind.grid_stride` thins this raster, and the two have no
    # reason to agree.
    jj_g = np.arange(0, grid.nymax - grid.ny0, stride)
    ii_g = np.arange(0, grid.nxmax - grid.nx0, stride)
    jg, ig = np.meshgrid(jj_g, ii_g, indexing="ij")

    # grid.x_stat / y_stat are ALREADY the cut subdomain (x[nx0:nxmax]),
    # so they are indexed with cut-relative indices directly. Slicing
    # them by [nx0:nxmax] a second time shortens the array and then
    # indexes off the end of it.
    x_abs = grid.abs_x(np.asarray(grid.x_stat))[ii_g]
    y_abs = grid.abs_y(np.asarray(grid.y_stat))[jj_g]

    # k index of the first level at or above the local surface + h
    k_of, inside = {}, {}
    for h in heights:
        k = np.searchsorted(z, grid.elev[jg, ig] + float(h))
        inside[h] = k < nz
        k_of[h] = np.clip(k, 0, nz - 1)
    _warn_above_top(inside, z)

    with Dataset(path, "w", format="NETCDF4") as fout:
        fout.Conventions = "CF-1.8"
        fout.title = "PALM wind field sampled at fixed heights above ground"
        fout.createDimension("time", len(its))
        fout.createDimension("height", len(heights))
        fout.createDimension("y", y_abs.size)
        fout.createDimension("x", x_abs.size)

        gm = write_crs_var(fout, grid.crs,
                           cfg["crs"]["netcdf_attributes"])

        xv = fout.createVariable("x", "f8", ("x",))
        xv[:] = x_abs
        xv.standard_name = "projection_x_coordinate"
        xv.units = "m"
        xv.axis = "X"

        # y stays ASCENDING. GeoTIFF is a north-up raster format and needs
        # its rows flipped; netCDF is not, and flipping here would mirror
        # the field.
        yv = fout.createVariable("y", "f8", ("y",))
        yv[:] = y_abs
        yv.standard_name = "projection_y_coordinate"
        yv.units = "m"
        yv.axis = "Y"

        hv = fout.createVariable("height", "f8", ("height",))
        hv[:] = np.asarray(heights, dtype="f8")
        hv.long_name = long_name("height_agl")
        hv.units = "m"
        hv.axis = "Z"
        hv.positive = "up"

        tv = fout.createVariable("time", "f8", ("time",))
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
            ov = fout.createVariable(name, "f4",
                                     ("time", "height", "y", "x"),
                                     fill_value=FILL, zlib=True, complevel=4)
            ov.grid_mapping = gm
            ov.standard_name = std
            ov.units = units
            ov.long_name = long_name(name)
        fout.variables["dir_met"].comment = (
            "Meteorological convention: the direction the wind blows FROM, "
            "measured clockwise from north. In ArcGIS Pro's vector field "
            "renderer set Reference System = Geographic and Angle "
            "Representation = From. Horizontal only, ignoring w, so it "
            "stays comparable with the point layer.")

        for q, it in enumerate(its):
            u3 = destagger(nc, "u", it, 0, nz, grid.ny0, grid.nymax,
                           grid.nx0, grid.nxmax)
            v3 = destagger(nc, "v", it, 0, nz, grid.ny0, grid.nymax,
                           grid.nx0, grid.nxmax)
            w3 = destagger(nc, "w", it, 0, nz, grid.ny0, grid.nymax,
                           grid.nx0, grid.nxmax)
            for hidx, h in enumerate(heights):
                k = k_of[h]
                # A sample above the file's top level is NoData, not
                # the top level's value relabelled as the requested height.
                u = np.where(inside[h], u3[k, jg, ig], np.nan)
                v = np.where(inside[h], v3[k, jg, ig], np.nan)
                w = np.where(inside[h], w3[k, jg, ig], np.nan)
                sp = np.hypot(u, v)
                dr = (270.0 - np.degrees(np.arctan2(v, u))) % 360.0
                for name, val in (("u", u), ("v", v), ("w", w),
                                  ("speed", sp), ("dir_met", dr)):
                    arr = np.array(val, dtype="f4")
                    arr[~np.isfinite(arr)] = FILL
                    fout.variables[name][q, hidx, ...] = arr

    info("[wind] wrote {}: {} x {} x {} x {} (t,height,y,x), "
         "u, v, w, speed, dir_met", cfg.disp(path), len(its),
         len(heights), y_abs.size, x_abs.size)
    debug("[wind] ArcGIS Pro: Add Multidimensional Raster Layer as a "
          "Vector Field")


# ------------------------------
# 2. STAGE ENTRY
# ------------------------------
def stage_wind(cfg, grid):
    wd = cfg["wind"]
    file_to = cfg.out_path(wd["file_to"]) if wd["file_to"] else None
    grid_to = cfg.out_path(wd["grid_to"]) if wd["grid_to"] else None
    want_gpkg = file_to is not None and cfg.prepare_output(file_to, "wind")
    want_grid = grid_to is not None and cfg.prepare_output(grid_to, "wind")
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
        z = np.asarray(nc.variables["zu_3d"][:])
        nz = z.size
        stride = int(wd["stride"])
        heights = [float(h) for h in wd["heights_agl"]]

        # thinned column selection (indices into the cut subdomain)
        jj = np.arange(0, grid.nymax - grid.ny0, stride)
        ii = np.arange(0, grid.nxmax - grid.nx0, stride)
        jg, ig = np.meshgrid(jj, ii, indexing="ij")
        jf, if_ = jg.ravel(), ig.ravel()
        debug("[wind] {} columns (stride {}), heights AGL: {} m",
              len(jf), stride, ", ".join(str(h) for h in heights))

        if want_grid:
            _wind_netcdf(cfg, grid, nc, its, times, heights, z, grid_to)
        if not want_gpkg:
            return

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

        x_abs = grid.abs_x(grid.x_stat)
        y_abs = grid.abs_y(grid.y_stat)
        pts = [Point(x_abs[i], y_abs[j]) for j, i in zip(jf, if_)]
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
                u = np.where(inside[h], u3[k, jf, if_], np.nan)
                v = np.where(inside[h], v3[k, jf, if_], np.nan)
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
                                              crs=grid.crs))
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
                               np.nan_to_num(u), np.nan_to_num(v),
                               np.nan_to_num(w))]
                a3 = dict(attrs)
                a3["z"] = zb.astype(np.float32)
                if w3 is not None:
                    a3["w"] = w.astype(np.float32)
                parts3.append(gpd.GeoDataFrame(a3, geometry=segs,
                                               crs=grid.crs))
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
