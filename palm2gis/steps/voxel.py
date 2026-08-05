"""3D output -> CF-compliant NetCDF for ArcGIS Pro Voxel Layers.

GIS counterpart of paraview_vis/utils/operation_3d.py. Destaggers u, v, w
onto scalar mass centres, cuts the subdomain given by `bounds`, and writes
a NetCDF that ArcGIS Pro accepts directly as a Voxel Layer:

  * dims ordered (time, z, y, x), regular spacing in x and y
  * x/y in absolute projected coordinates with CF standard names
  * a grid_mapping variable carrying the CRS (WKT + Esri PE string)
  * time as CF "seconds since <origin_time>"
  * masked (terrain/building) cells as NoData via _FillValue

NOTE on the original code: operation_3d.py averaged staggered components as
`(var[..., 1:] + var[..., -1:]) / 2` - the `-1:` broadcasts the LAST slice
against everything, so every u value was averaged with the domain-edge
column instead of its neighbour. Fixed here (`:-1` + `1:`).

By default the k=0 level (z = 0) is dropped so vertical spacing is uniform
(dz/2, 3dz/2, ...), which voxel layers handle best (voxel.drop_k0).
"""

import numpy as np
from netCDF4 import Dataset

from ..log import debug, get_logger, info, warning
from .georef import check_raw_input, select_times
from .thermo import LAPSE_DRY, exner_profile, is_temperature, \
    kelvin_to_celsius, theta_to_ta

log = get_logger()

FILL = -9999.0


# ------------------------------
# 1. DESTAGGERING
# ------------------------------
def _pad_last(v, axis):
    """Append a copy of the last slice along `axis`: the staggered dim ended
    at the domain edge, so the outermost available face is reused."""
    idx = [slice(None)] * v.ndim
    idx[axis] = slice(-1, None)
    return np.concatenate([v, v[tuple(idx)]], axis=axis)


def destagger(nc, var, it, nz0, nzmax, ny0, nymax, nx0, nxmax):
    """Interpolate one staggered velocity component to mass centres.

    Handles both staggered-dimension conventions: xu/yv/zw one point LONGER
    than the scalar dims (idealised files), and the SAME length (real PALM
    3D output, where the outermost face is not stored). In the latter case,
    when a cut reaches the domain edge the slice comes back one face short
    and the last available face is reused - a zero-gradient assumption for
    one half-cell at the very edge.
    """
    if var == "u":
        v = np.ma.filled(nc.variables[var][it, nz0:nzmax, ny0:nymax,
                                           nx0:nxmax + 1], 0.0)
        if v.shape[2] < (nxmax - nx0) + 1:
            v = _pad_last(v, 2)
        return 0.5 * (v[:, :, :-1] + v[:, :, 1:])
    if var == "v":
        v = np.ma.filled(nc.variables[var][it, nz0:nzmax, ny0:nymax + 1,
                                           nx0:nxmax], 0.0)
        if v.shape[1] < (nymax - ny0) + 1:
            v = _pad_last(v, 1)
        return 0.5 * (v[:, :-1, :] + v[:, 1:, :])
    if var == "w":
        v = np.ma.filled(nc.variables[var][it, nz0 - 1 if nz0 > 0 else 0:
                                           nzmax, ny0:nymax, nx0:nxmax],
                         0.0)
        if nz0 > 0:
            if v.shape[0] < (nzmax - nz0) + 1:
                v = _pad_last(v, 0)
            return 0.5 * (v[:-1, :, :] + v[1:, :, :])
        # k=0: w on zw(0)=0 approximately equals the first mass level
        if v.shape[0] < nzmax - nz0:
            v = _pad_last(v, 0)
        out = np.empty((nzmax - nz0,) + v.shape[1:])
        out[1:] = 0.5 * (v[:-1, :, :] + v[1:, :, :])
        out[0] = v[0]
        return out
    raise ValueError(var)


# ------------------------------
# 2. CF / CRS METADATA
# ------------------------------
def write_crs_var(fout, crs):
    """Write a grid_mapping variable ArcGIS and GDAL both understand."""
    cv = fout.createVariable("crs", "i4")
    for k, v in crs.to_cf().items():
        if k == "crs_wkt":
            continue
        try:
            cv.setncattr(k, v)
        except TypeError:
            pass
    cv.spatial_ref = crs.to_wkt()
    cv.crs_wkt = crs.to_wkt()
    cv.esri_pe_string = crs.to_wkt(version="WKT1_ESRI")
    return "crs"


# ------------------------------
# 3. STAGE ENTRY
# ------------------------------
def stage_voxel(cfg, grid):
    vx = cfg["voxel"]
    file_from = cfg.resolve(vx["file_from"])
    file_to = cfg.out_path(vx["file_to"])
    if not cfg.prepare_output(file_to, "voxel"):
        return

    with Dataset(file_from, "r") as nc:
        check_raw_input(nc, grid, "voxel", vx["file_from"])
        its, times = select_times(cfg, nc, "voxel", grid, cfg["voxel"]["file_from"])

        z_all = np.asarray(nc.variables["zu_3d"][:])
        k_off = 1 if (vx["drop_k0"] and grid.nz0 == 0) else 0
        kz0, kzmax = grid.nz0 + k_off, grid.nzmax
        z_out = z_all[kz0:kzmax]

        scalars = [s for s in vx["scalars"] if s]
        # Air temperature is ENSURED rather than always derived: PALM's own
        # `ta` is exact, so it is copied when the file has it, and `ta` is
        # reconstructed from `theta` only when it does not.
        derive_ta = False
        if vx["derive_ta"] and "ta" not in scalars:
            if "ta" in nc.variables:
                info("[voxel] the file holds 'ta' - using PALM's own air "
                     "temperature instead of deriving it.")
                scalars = scalars + ["ta"]
            elif "theta" in nc.variables:
                derive_ta = True
                if "theta" not in scalars:
                    info("[voxel] no 'ta' in the file - deriving it from "
                         "'theta'.")
                    scalars = scalars + ["theta"]
            else:
                warning("[voxel] derive_ta is on but the file holds neither "
                        "'ta' nor 'theta'.")
        out_vars = (["u", "v", "w"] if vx["do_velocity"] else []) + scalars
        if derive_ta:
            out_vars = out_vars + ["ta"]

        with cfg.atomic_output(file_to, "voxel") as _tmp_out, \
                Dataset(_tmp_out, "w", format="NETCDF4") as fout:
            fout.Conventions = "CF-1.8"
            fout.title = "PALM 3D output for ArcGIS Pro voxel layer"
            fout.source = str(file_from)

            fout.createDimension("x", grid.nxmax - grid.nx0)
            fout.createDimension("y", grid.nymax - grid.ny0)
            fout.createDimension("z", z_out.size)
            fout.createDimension("time", len(its))

            xv = fout.createVariable("x", "f8", ("x",))
            xv[:] = grid.abs_x(
                np.asarray(nc.variables["x"][grid.nx0:grid.nxmax]))
            xv.standard_name = "projection_x_coordinate"
            xv.units = "m"
            xv.axis = "X"

            yv = fout.createVariable("y", "f8", ("y",))
            yv[:] = grid.abs_y(
                np.asarray(nc.variables["y"][grid.ny0:grid.nymax]))
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

            tv = fout.createVariable("time", "f8", ("time",))
            t_sel = np.asarray(times)[its]
            if np.all(np.isfinite(t_sel)) and np.all(np.abs(t_sel) < 1.0e12):
                tv[:] = np.round(t_sel)
                tv.standard_name = "time"
                tv.units = "seconds since {}".format(
                    grid.origin_time.strftime("%Y-%m-%d %H:%M:%S"))
                tv.calendar = "standard"
            else:
                warning("[voxel] no usable time values in the input - the "
                        "time axis carries record indices, not real times.")
                tv[:] = np.arange(len(its), dtype="f8")
                tv.long_name = ("record index (the source file had no valid "
                                "time values)")
                tv.units = "1"
            tv.axis = "T"

            gm = write_crs_var(fout, grid.crs)

            pi = None
            if derive_ta:
                # Exner profile from the horizontally averaged theta of the
                # first selected timestep, on absolute model heights.
                th0 = np.ma.masked_invalid(
                    nc.variables["theta"][its[0], kz0:kzmax,
                                          grid.ny0:grid.nymax,
                                          grid.nx0:grid.nxmax])
                theta_bar = np.ma.filled(th0.mean(axis=(1, 2)), np.nan)
                pi = exner_profile(z_out, theta_bar,
                                   p_surface=cfg["surface_pressure"])
                # What the conversion actually DOES to the numbers, per
                # level: ta - theta = theta * (Pi - 1).
                off = np.asarray(theta_bar) * (pi - 1.0)
                o0, o1 = float(off[0]), float(off[-1])
                z0, z1 = float(z_out[0]), float(z_out[-1])
                if o0 >= 0.0 and o1 >= 0.0:
                    info("[voxel] deriving ta from theta: ta comes out "
                         "{:.2f} K warmer at z = {:.0f} m, {:.2f} K warmer "
                         "at {:.0f} m.", o0, z0, o1, z1)
                elif o0 <= 0.0 and o1 <= 0.0:
                    info("[voxel] deriving ta from theta: ta comes out "
                         "{:.2f} K cooler at z = {:.0f} m, {:.2f} K cooler "
                         "at {:.0f} m.", -o0, z0, -o1, z1)
                else:
                    info("[voxel] deriving ta from theta: ta - theta = "
                         "{:+.2f} K at z = {:.0f} m, {:+.2f} K at {:.0f} m.",
                         o0, z0, o1, z1)

                # The offset must decay with height at the dry adiabatic
                # rate whatever the stratification, because theta * dPi/dz
                # = -g/cp. A departure means the profile, p_surface or the
                # height axis is wrong - worth saying so rather than
                # writing plausible-looking temperatures.
                lapse = (o0 - o1) / (z1 - z0) * 1000.0 if z1 > z0 else None
                if lapse is not None and abs(lapse - LAPSE_DRY) > 0.5:
                    warning("[voxel] the theta -> ta offset changes at "
                            "{:.1f} K/km, but a constant-theta column is a "
                            "dry adiabat and must give {:.1f} K/km - check "
                            "surface_pressure and the z axis.",
                            lapse, LAPSE_DRY)
                elif lapse is not None:
                    debug("[voxel] offset decays at {:.1f} K/km (dry "
                          "adiabat {:.1f}); Exner {:.5f}..{:.5f} from "
                          "p_s = {:.0f} Pa at z = 0.",
                          lapse, LAPSE_DRY, pi[0], pi[-1],
                          cfg["surface_pressure"])
                pi = pi[:, None, None]

            celsius = bool(cfg["celsius"])
            for var in out_vars:
                ov = fout.createVariable(var, "f4", ("time", "z", "y", "x"),
                                         fill_value=FILL, zlib=True,
                                         complevel=4)
                ov.grid_mapping = gm
                ov.coordinates = "time z y x"
                if var == "ta" and derive_ta:
                    ov.units = "degC" if celsius else "K"
                    ov.long_name = "air temperature"
                    ov.comment = ("derived from theta via a hydrostatic "
                                  "Exner profile (PALM did not output ta)")
                    ov.derived_from = "theta"
                    continue
                src = nc.variables.get(var)
                if src is not None:
                    for att in ("units", "long_name"):
                        if hasattr(src, att):
                            ov.setncattr(att, getattr(src, att))
                if celsius and is_temperature(var):
                    ov.units = "degC"

            # terrain mask from the first mass-grid variable at first time
            mask_src = next((c for c in scalars + ["w"]
                             if c in nc.variables), None)
            terrain_mask = np.ma.getmaskarray(
                nc.variables[mask_src][its[0], kz0:kzmax,
                                       grid.ny0:grid.nymax,
                                       grid.nx0:grid.nxmax])
            debug("[voxel] terrain mask from '{}': {} cells",
                  mask_src, int(terrain_mask.sum()))

            for q, it in enumerate(its):
                debug("[voxel] time {} ({:.0f} s)", it, times[it])
                if vx["do_velocity"]:
                    for var in ("u", "v", "w"):
                        val = destagger(nc, var, it, kz0, kzmax,
                                        grid.ny0, grid.nymax,
                                        grid.nx0, grid.nxmax)
                        val[terrain_mask] = FILL
                        fout.variables[var][q, ...] = val
                for var in scalars:
                    val = np.ma.filled(
                        nc.variables[var][it, kz0:kzmax,
                                          grid.ny0:grid.nymax,
                                          grid.nx0:grid.nxmax], FILL)
                    if celsius and is_temperature(var):
                        val = kelvin_to_celsius(val, fill=FILL)
                    val[terrain_mask] = FILL
                    fout.variables[var][q, ...] = val
                if derive_ta:
                    th = np.ma.masked_invalid(np.ma.masked_equal(
                        nc.variables["theta"][it, kz0:kzmax,
                                              grid.ny0:grid.nymax,
                                              grid.nx0:grid.nxmax], FILL))
                    ta = np.ma.filled(theta_to_ta(th, pi), FILL)
                    if celsius:
                        ta = kelvin_to_celsius(ta, fill=FILL)
                    ta[terrain_mask] = FILL
                    fout.variables["ta"][q, ...] = ta

    info("[voxel] wrote {}: {} x {} x {} x {} (t,z,y,x), {}",
         cfg.disp(file_to), len(its), z_out.size,
         grid.nymax - grid.ny0, grid.nxmax - grid.nx0, ", ".join(out_vars))
