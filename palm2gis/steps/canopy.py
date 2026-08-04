"""Plant canopy (LAD) from the static driver -> GIS layers.

GIS counterpart of paraview_vis/utils/tree_gen.py, the one stage of the
original tool that had no equivalent here.

The static driver holds `lad` (leaf area density) on the `zlad` axis,
which is height ABOVE THE LOCAL TERRAIN, not absolute height - a tree on a
hillside is stored at the same indices as one in a valley. The original
tool lifted that into absolute-z space by offsetting each column by its
terrain height, producing a 3D field for ParaView volume rendering. The
same lift is done here, and the result is served three ways:

  * LAI raster - LAD integrated over the column (sum * dz_lad), the
    standard 2D measure of leaf area above each square metre. The layer
    most urban-climate maps actually want, and the cheapest to produce.
  * canopy height raster - the highest level with LAD above a threshold,
    in metres above the LOCAL terrain (a CHM). Derived from LAD itself,
    so it does not depend on `buildings_2d`, which some static drivers
    leave empty.
  * absolute-z LAD volume as a CF NetCDF - the ArcGIS Pro voxel / QGIS
    mesh equivalent of what paraview_vis produced for ParaView.

NOTE on the original: tree_gen.py used `np.int`, removed in NumPy 1.24, so
it raises AttributeError on any current install; it also ignored the
`bounds` z-cut and always wrote the full column. Both are fixed here.
"""

import numpy as np
from netCDF4 import Dataset

import rasterio
from rasterio.transform import from_origin

from ..log import debug, get_logger, info, warning
from .voxel import write_crs_var

log = get_logger()

FILL = -9999.0


# ------------------------------
# 1. LIFT LAD TO ABSOLUTE HEIGHT
# ------------------------------
def lift_to_absolute(lad_rel, zt, dz, nz):
    """Move LAD from height-above-terrain onto the absolute model grid.

    lad_rel : (nzlad, ny, nx) as stored in the static driver.
    zt      : (ny, nx) terrain height of the same cut.
    Returns (nz, ny, nx) on the absolute z axis.

    Vectorised over the vertical: the original looped over every (j, i)
    column in Python, which on a 1000 x 1000 domain is a million
    iterations.
    """
    nzlad, ny, nx = lad_rel.shape
    k_terrain = np.round(np.asarray(zt) / dz).astype(int)
    out = np.zeros((nz, ny, nx), dtype="f4")
    kk = np.arange(nz)[:, None, None]
    for m in range(nzlad):
        target = k_terrain[None, :, :] + 1 + m       # absolute level index
        valid = (target >= 0) & (target < nz)
        np.add(out, np.where(valid & (kk == target), lad_rel[m], 0.0),
               out=out)
    return out


# ------------------------------
# 2. RASTER WRITER
# ------------------------------
def _write_raster(cfg, grid, data, path, description, band_name):
    """One single-band GeoTIFF, north-up, in absolute coordinates."""
    ny, nx = data.shape
    west = grid.origin_x + grid.x_stat[0] - grid.dx / 2.0
    north = grid.origin_y + grid.y_stat[-1] + grid.dy / 2.0
    transform = from_origin(west, north, grid.dx, grid.dy)
    with rasterio.open(path, "w", driver="GTiff", height=ny, width=nx,
                       count=1, dtype="float32", crs=grid.rio_crs(),
                       transform=transform, nodata=FILL,
                       compress="deflate") as dst:
        dst.write(data[::-1, :].astype("float32"), 1)   # north-up
        dst.set_band_description(1, band_name)
        dst.update_tags(description=description)
    info("[canopy] wrote {}: {}", cfg.disp(path), description)


# ------------------------------
# 3. STAGE ENTRY
# ------------------------------
def stage_canopy(cfg, grid):
    cn = cfg["canopy"]
    src = cfg.resolve(cn["file_from"] or cfg["domain"]["static_driver"])

    with Dataset(src, "r") as nc:
        if "lad" not in nc.variables:
            warning("[canopy] no 'lad' in {} - nothing to convert. Trees "
                    "may be described only by tree_id / vegetation_type in "
                    "this static driver.", cfg.disp(src))
            return
        zlad = (np.asarray(nc.variables["zlad"][:])
                if "zlad" in nc.variables else None)
        lad_rel = np.ma.filled(
            nc.variables["lad"][:, grid.ny0:grid.nymax,
                                grid.nx0:grid.nxmax], 0.0)

    lad_rel = np.where(np.isfinite(lad_rel) & (lad_rel > 0), lad_rel, 0.0)
    if not lad_rel.any():
        warning("[canopy] 'lad' is present but empty inside the bounds cut "
                "- no vegetation here.")
        return

    dz_lad = (float(zlad[1] - zlad[0])
              if zlad is not None and zlad.size > 1 else grid.dz)
    info("[canopy] {} canopy levels, dz_lad = {:g} m, LAD max = "
         "{:.3f} m2/m3", lad_rel.shape[0], dz_lad, float(lad_rel.max()))

    # -- LAI: vertical integral of LAD -------------------------------------
    if cn["write_lai"]:
        path = cfg.out_path("{}lai.tif".format(cn["prefix"]))
        if cfg.prepare_output(path, "canopy"):
            lai = lad_rel.sum(axis=0) * dz_lad
            _write_raster(cfg, grid, np.where(lai > 0, lai, FILL), path,
                          "leaf area index (LAD integrated over the "
                          "column, m2/m2)", "LAI")
            debug("[canopy] LAI max {:.2f} over {} vegetated cells",
                  float(lai.max()), int((lai > 0).sum()))

    # -- canopy height model ----------------------------------------------
    if cn["write_canopy_height"]:
        path = cfg.out_path("{}canopy_height.tif".format(cn["prefix"]))
        if cfg.prepare_output(path, "canopy"):
            thr = float(cn["lad_threshold"])
            above = lad_rel > thr
            has = above.any(axis=0)
            # index of the highest level above the threshold
            top_k = (lad_rel.shape[0] - 1
                     - np.argmax(above[::-1, :, :], axis=0))
            chm = np.where(has, (top_k + 1) * dz_lad, FILL)
            _write_raster(cfg, grid, chm, path,
                          "canopy height above local terrain (m), "
                          "LAD > {:g}".format(thr), "canopy_height")
            debug("[canopy] canopy height max {:.1f} m over {} cells",
                  float(chm[chm != FILL].max()) if has.any() else 0.0,
                  int(has.sum()))

    # -- absolute-z LAD volume --------------------------------------------
    if cn["write_volume"]:
        path = cfg.out_path("{}lad_3d.nc".format(cn["prefix"]))
        if not cfg.prepare_output(path, "canopy"):
            return
        z_all = grid.z_levels()
        lad_abs = lift_to_absolute(lad_rel, grid.zt, grid.dz, grid.nz)
        kz0, kzmax = grid.nz0, grid.nzmax
        lad_cut = lad_abs[kz0:kzmax]
        z_out = z_all[kz0:kzmax]
        if not lad_cut.any():
            warning("[canopy] the canopy falls entirely outside the "
                    "vertical cut z[{}:{}] - raise bounds.zmax.",
                    kz0, kzmax)

        with Dataset(path, "w", format="NETCDF4") as fout:
            fout.Conventions = "CF-1.8"
            fout.title = "PALM plant canopy (LAD) on absolute heights"
            fout.source = str(src)

            fout.createDimension("x", grid.nxmax - grid.nx0)
            fout.createDimension("y", grid.nymax - grid.ny0)
            fout.createDimension("z", z_out.size)

            xv = fout.createVariable("x", "f8", ("x",))
            xv[:] = grid.abs_x(grid.x_stat)
            xv.standard_name = "projection_x_coordinate"
            xv.units = "m"
            xv.axis = "X"

            yv = fout.createVariable("y", "f8", ("y",))
            yv[:] = grid.abs_y(grid.y_stat)
            yv.standard_name = "projection_y_coordinate"
            yv.units = "m"
            yv.axis = "Y"

            zv = fout.createVariable("z", "f8", ("z",))
            zv[:] = z_out
            zv.standard_name = "height"
            zv.long_name = "height above domain base"
            zv.units = "m"
            zv.axis = "Z"
            zv.positive = "up"

            gm = write_crs_var(fout, grid.crs)
            ov = fout.createVariable("lad", "f4", ("z", "y", "x"),
                                     fill_value=FILL, zlib=True,
                                     complevel=4)
            ov.grid_mapping = gm
            ov.coordinates = "z y x"
            ov.units = "m2 m-3"
            ov.long_name = "leaf area density"
            ov[:] = np.where(lad_cut > 0, lad_cut, FILL)

        info("[canopy] wrote {}: {} x {} x {} (z,y,x), {} canopy cells",
             cfg.disp(path), z_out.size, grid.nymax - grid.ny0,
             grid.nxmax - grid.nx0, int((lad_cut > 0).sum()))
