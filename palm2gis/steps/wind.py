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

from ..log import debug, get_logger, info, warning
from .georef import check_product, check_raw_input, select_times
from .voxel import destagger
from .qgis_style import (embed_gpkg_style, vector_qml, write_sidecar)

log = get_logger()


# ------------------------------
# 1. STAGE ENTRY
# ------------------------------
def stage_wind(cfg, grid):
    wd = cfg["wind"]
    file_to = cfg.out_path(wd["file_to"])
    if not cfg.prepare_output(file_to, "wind"):
        return

    with Dataset(cfg.resolve(wd["file_from"]), "r") as nc:
        check_product(nc, "wind", wd["file_from"])
        check_raw_input(nc, grid, "wind", wd["file_from"])
        its, times = select_times(cfg, nc, "wind", grid, cfg["wind"]["file_from"])

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

        lines_3d = bool(wd["write_lines_3d"])
        include_w = bool(wd["include_w"])
        scale = float(wd["arrow_scale"])

        # per column: k index of the first level >= local surface + h_agl
        k_idx = {}
        z_abs = {}
        for h in heights:
            k = np.searchsorted(z, grid.elev[jf, if_] + h)
            k_idx[h] = np.clip(k, 0, nz - 1)
            # The true height of the sample, for the 3D geometry. The 2D
            # point layer carried this only as the `height_agl` attribute,
            # so every level collapsed onto the ground in a 3D scene and
            # the 2 m and 10 m samples sat exactly on top of each other.
            z_abs[h] = grid.elev[jf, if_] + h

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
                u = u3[k, jf, if_]
                v = v3[k, jf, if_]
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
                w = (w3[k, jf, if_] if w3 is not None
                     else np.zeros_like(u))
                zb = z_abs[h]
                segs = [LineString([(x0, y0, z0),
                                    (x0 + ux * scale, y0 + vy * scale,
                                     z0 + wz * scale)])
                        for x0, y0, z0, ux, vy, wz
                        in zip(x_abs[if_], y_abs[jf], zb, u, v, w)]
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
