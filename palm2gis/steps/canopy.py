"""Plant canopy (LAD) from the static driver -> GIS layers.

GIS counterpart of paraview_vis/utils/tree_gen.py, the one stage of the
original tool that had no equivalent here.

The static driver holds `lad` (leaf area density) on the `zlad` axis,
which is height ABOVE THE LOCAL TERRAIN, not absolute height - a tree on a
hillside is stored at the same indices as one in a valley. The original
tool lifted that into absolute-z space by offsetting each column by its
terrain height, producing a 3D field for ParaView volume rendering. The
same lift is done here, and the result is written as an absolute-z LAD
volume in a CF NetCDF - the ArcGIS Pro voxel / QGIS mesh equivalent of what
paraview_vis produced for ParaView.

NOTE on the original: tree_gen.py used `np.int`, removed in NumPy 1.24, so
it raises AttributeError on any current install; it also ignored the
`bounds` z-cut and always wrote the full column. Both are fixed here.
"""

import numpy as np
from netCDF4 import Dataset

from ..log import debug, get_logger, info, warning
from .georef import parse_time_bound
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
# 2. STAGE ENTRY
# ------------------------------
def stage_canopy(cfg, grid):
    cn = cfg["canopy"]
    src = cfg.resolve(cn["file_from"] or cfg["domain"]["static_driver"])
    if not cn["write_volume"]:
        info("[canopy] write_volume is off - nothing to write, skipping.")
        return

    with Dataset(src, "r") as nc:
        if "lad" not in nc.variables:
            warning("[canopy] no 'lad' in {} - nothing to convert.",
                    cfg.disp(src))
            return
        zlad = (np.asarray(nc.variables["zlad"][:])
                if "zlad" in nc.variables else None)
        lad_rel = np.ma.filled(
            nc.variables["lad"][:, grid.ny0:grid.nymax,
                                grid.nx0:grid.nxmax], 0.0)
        # `bad` (basal area density: branch and trunk area per volume) is
        # the woody companion to `lad`, on the same zlad axis with the same
        # shape, so it needs no separate handling. It is OPTIONAL - PALM's
        # own use of it was still under development when the static-driver
        # specification was published - so its absence is normal and silent
        # beyond a debug line. Including it is what makes a 3D canopy read
        # as trees rather than floating crowns.
        bad_rel = None
        if "bad" in nc.variables:
            bad_rel = np.ma.filled(
                nc.variables["bad"][:, grid.ny0:grid.nymax,
                                    grid.nx0:grid.nxmax], 0.0)

    lad_rel = np.where(np.isfinite(lad_rel) & (lad_rel > 0), lad_rel, 0.0)
    if bad_rel is not None:
        bad_rel = np.where(np.isfinite(bad_rel) & (bad_rel > 0), bad_rel, 0.0)
        if bad_rel.shape != lad_rel.shape:
            warning("[canopy] 'bad' has shape {} but 'lad' has {} - the two "
                    "should share the zlad axis; skipping bad.",
                    bad_rel.shape, lad_rel.shape)
            bad_rel = None
        elif not bad_rel.any():
            debug("[canopy] 'bad' is present but empty inside the cut.")
            bad_rel = None
    if not lad_rel.any():
        warning("[canopy] 'lad' is present but empty inside the bounds cut "
                "- no vegetation here.")
        return

    dz_lad = (float(zlad[1] - zlad[0])
              if zlad is not None and zlad.size > 1 else grid.dz)
    info("[canopy] {} canopy levels, dz_lad = {:g} m, LAD max = "
         "{:.3f} m2/m3{}", lad_rel.shape[0], dz_lad, float(lad_rel.max()),
         ", BAD max = {:.3f} m2/m3".format(float(bad_rel.max()))
         if bad_rel is not None else "")

    # -- absolute-z LAD volume --------------------------------------------
    path = cfg.out_path(cn["volume_to"]
                        or "{}lad_3d.nc".format(cn["prefix"]))
    if not cfg.prepare_output(path, "canopy"):
        return
    z_all = grid.z_levels()
    kz0, kzmax = grid.nz0, grid.nzmax
    # Absolute height for the output axis. lift_to_absolute() below is
    # unrelated - it lifts the terrain-following LAD onto the model
    # grid and still works in RELATIVE levels, as grid.zt is relative.
    z_out = grid.abs_z(z_all[kz0:kzmax])

    vols = {"lad": ("leaf area density",
                    lift_to_absolute(lad_rel, grid.zt, grid.dz,
                                     grid.nz)[kz0:kzmax])}
    if bad_rel is not None:
        vols["bad"] = ("basal area density (branches and trunks)",
                       lift_to_absolute(bad_rel, grid.zt, grid.dz,
                                        grid.nz)[kz0:kzmax])
    lad_cut = vols["lad"][1]
    if not lad_cut.any():
        warning("[canopy] the canopy falls entirely outside the "
                "vertical cut z[{}:{}] - raise bounds.zmax.",
                kz0, kzmax)

    # Nominal timestamp for the static field: the start of the window
    # the rest of the run converts, so it lines up with the first
    # timestep of the voxel/surface layers in the same scene. Falls
    # back to origin_time when no window is set. Resolved WITHOUT
    # opening a 3D output file, so the canopy stage stays independent
    # of whether the voxel stage runs at all.
    t_static, t_dt = parse_time_bound(cfg["bounds"].get("time_from"),
                                      grid.origin_time,
                                      "bounds.time_from")
    if t_static is None:
        t_static = 0.0
    t_label = (grid.times_to_datetime([t_static])[0]
               or "{:g} s".format(t_static))

    with Dataset(path, "w", format="NETCDF4") as fout:
        fout.Conventions = "CF-1.8"
        fout.title = ("PALM plant canopy ({}) on absolute heights"
                      .format(", ".join(sorted(vols)).upper()))
        fout.source = str(src)

        fout.createDimension("time", 1)
        # Horizontal axes in the output CRS (PALM's own, or resampled).
        og = grid.reproj.grid()
        yd, xd = og.dims
        og.write_axes(fout)
        fout.createDimension("z", z_out.size)

        zv = fout.createVariable("z", "f8", ("z",))
        zv[:] = z_out
        zv.standard_name = "height"
        zv.long_name = ("height above sea level (PALM zu + origin_z)")
        zv.palm2gis_origin_z = float(grid.origin_z)
        zv.units = "m"
        zv.axis = "Z"
        zv.positive = "up"

        # The canopy is STATIC, but a scene that animates the voxel
        # output needs this volume to carry a time axis too, or it
        # drops out of view as soon as the time slider moves. One
        # timestep is enough for that and costs nothing in size; the
        # value is the start of the window the rest of the run is
        # converting, so the canopy sits at the beginning of the
        # scene's temporal range rather than at an unrelated instant.
        tv = fout.createVariable("time", "f8", ("time",))
        tv[:] = [t_static]
        tv.standard_name = "time"
        tv.units = "seconds since {}".format(
            grid.origin_time.strftime("%Y-%m-%d %H:%M:%S"))
        tv.calendar = "standard"
        tv.axis = "T"
        tv.long_name = ("nominal time of a static field (the canopy "
                        "does not vary in time)")

        gm = write_crs_var(fout, grid.reproj.crs,
                           cfg["crs"]["netcdf_attributes"])
        for name, (long_name, cut) in sorted(vols.items()):
            ov = fout.createVariable(name, "f4",
                                     ("time", "z", yd, xd),
                                     fill_value=FILL, zlib=True,
                                     complevel=4)
            ov.grid_mapping = gm
            # No `coordinates` attribute: it declares AUXILIARY
            # coordinates in CF, and Esri rejects a multidimensional
            # source that has any.
            ov.units = "m2 m-3"
            ov.long_name = long_name
            # Density, zeros included: a conservative remap then keeps
            # the total leaf area, and "no vegetation" becomes NoData only
            # after resampling.
            dens = og.take(np.asarray(cut, dtype="f8"), fill=0.0,
                           valid=np.ones(cut.shape, dtype=bool))
            ov[0] = np.where(dens > 0, dens, FILL)

    info("[canopy] wrote {}: 1 x {} x {} x {} (t,z,{},{}), {} at {}, "
         "{} canopy cells",
         cfg.disp(path), z_out.size, og.shape[0], og.shape[1], yd, xd, ", ".join(sorted(vols)),
         t_label, int((lad_cut > 0).sum()))
