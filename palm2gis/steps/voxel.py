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
from .georef import (check_product, check_raw_input, require_variables,
                     resolve_z_range, select_times)
from .thermo import LAPSE_DRY, convert_temperature, exner_profile, \
    fill_profile_gaps, is_celsius_exempt, is_temperature, \
    temperature_scale, theta_to_ta

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
    cf = crs.to_cf()
    for k, v in cf.items():
        if k == "crs_wkt":
            continue
        try:
            cv.setncattr(k, v)
        except TypeError:
            pass

    # Esri's netCDF reader builds the spatial reference from
    # grid_mapping_name and its parameters. It does NOT read crs_wkt or
    # epsg_code, so a CRS pyproj cannot express in CF (Krovak /
    # S-JTSK EPSG:5514, for one) leaves ArcGIS with nothing to work from.
    if "grid_mapping_name" not in cf:
        warning("[crs] {} has no CF grid mapping - pyproj cannot express "
                "it as one. GDAL and QGIS will read the WKT, but ArcGIS "
                "reads grid_mapping_name and will treat this file as "
                "having an unknown CRS. Reproject to a CRS with a CF "
                "grid mapping (e.g. UTM) if the output is for ArcGIS.",
                crs.name)

    # Esri also assumes the WGS 1984 datum for netCDF regardless of what
    # the file says, so a non-WGS84 datum silently shifts the data.
    # EPSG names the WGS 84 datum several ways ("WGS 84", "WGS_1984",
    # "World Geodetic System 1984 ensemble"); all of them are the datum
    # ArcGIS assumes, so none of them should warn.
    datum = (crs.datum.name if crs.datum is not None else "")
    _d = datum.lower().replace("_", " ")
    is_wgs84 = "wgs 84" in _d or "world geodetic system 1984" in _d
    if datum and not is_wgs84:
        warning("[crs] datum is '{}', not WGS 84. ArcGIS assumes WGS 1984 "
                "for netCDF whatever the file declares, so expect a datum "
                "shift when this is opened in ArcGIS Pro.", datum)
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
        check_product(nc, "voxel", vx["file_from"])
        check_raw_input(nc, grid, "voxel", vx["file_from"])
        its, times = select_times(cfg, nc, "voxel", grid, cfg["voxel"]["file_from"])

        z_all = np.asarray(nc.variables["zu_3d"][:])
        nz0_f, nzmax_f = resolve_z_range(nc, grid, "voxel")
        k_off = 1 if (vx["drop_k0"] and nz0_f == 0) else 0
        kz0, kzmax = nz0_f + k_off, nzmax_f
        # Two z arrays, deliberately.
        #   z_rel - height above the domain base, what PALM wrote. The
        #           Exner profile below pairs it with surface_pressure,
        #           which is the pressure at the model's z = 0.
        #   z_out - absolute height, origin_z added. Output geometry only.
        # Collapsing these would trade a geometry error for a temperature
        # one.
        z_rel = z_all[kz0:kzmax]
        z_out = grid.abs_z(z_rel)

        # An ArcGIS voxel layer requires a REGULARLY gridded z axis.
        # Two things break it: a PALM stretched grid (dz_stretch_level),
        # and a cut that keeps k=0, whose spacing runs dz/2 then dz.
        if z_rel.size > 2:
            dz_steps = np.diff(z_rel)
            spread = float(dz_steps.max() - dz_steps.min())
            if spread > 1.0e-3 * max(1.0, float(abs(dz_steps.mean()))):
                raise RuntimeError(
                    "[voxel] the z axis is not regularly spaced "
                    "(steps {:.3f}..{:.3f} m), which an ArcGIS voxel layer "
                    "does not accept. Either the PALM run uses a stretched "
                    "grid (dz_stretch_level), in which case a voxel layer "
                    "is not possible from these levels, or the cut keeps "
                    "k = 0 - set voxel.drop_k0: true, or raise bounds.z0 "
                    "above the first level.".format(
                        float(dz_steps.min()), float(dz_steps.max())))

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
        do_velocity = bool(vx["do_velocity"])
        if do_velocity:
            require_variables(
                nc, ["u", "v", "w"], "voxel", vx["file_from"],
                hint="Set voxel.do_velocity: false to write scalars only, "
                     "or point file_from at the joined file.")
        require_variables(nc, [v for v in scalars if v != "ta"], "voxel",
                          vx["file_from"],
                          hint="Adjust voxel.scalars to what this file "
                               "holds, or point file_from at the joined "
                               "file.")

        # P7.3: the staggered axes must have the same number of levels.
        # splitz cuts zu_3d and zw_3d independently, so a grouped u/v/w
        # file could arrive with zu_3d and zw_3d of different lengths and
        # the terrain mask then failed with a bare numpy boolean-index
        # error deep in the write loop.
        if do_velocity and "zw_3d" in nc.variables:
            n_zu = nc.variables["zu_3d"].size
            n_zw = nc.variables["zw_3d"].size
            if n_zu != n_zw:
                raise RuntimeError(
                    "[voxel] zu_3d has {} levels but zw_3d has {}. The "
                    "staggered axes have been cut to different lengths, "
                    "which happens when palm_postproc splitz slices only "
                    "one of them. Read the joined file, or re-run splitz "
                    "so both axes are cut alike.".format(n_zu, n_zw))

        out_vars = (["u", "v", "w"] if do_velocity else []) + scalars
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
            zv.long_name = ("height above sea level (PALM zu_3d + "
                            "origin_z)")
            zv.palm2gis_origin_z = float(grid.origin_z)
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
                pi = exner_profile(z_rel, theta_bar,
                                   p_surface=cfg["surface_pressure"])
                # Report from the GAP-FILLED profile, which is what the
                # integration used. The raw means are NaN on any level
                # lying entirely inside terrain - normal in a coarse
                # domain, where the lowest levels can sit below the terrain
                # across the whole cut - and reporting those gave a "+nan K"
                # that looked like a failure when nothing was wrong.
                th_fill, n_gap = fill_profile_gaps(theta_bar)
                if n_gap:
                    info("[voxel] {} of {} levels are entirely inside "
                         "terrain (lowest air at z = {:.0f} m); their theta "
                         "is interpolated for the hydrostatic integration, "
                         "and ta there is terrain fill either way.",
                         n_gap, len(z_rel),
                         float(z_out[int(np.argmax(np.isfinite(theta_bar)))]))
                off = th_fill * (pi - 1.0)
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
            target = "C" if celsius else "K"
            overrides = dict(cfg["temperature_units"] or {})
            # Resolved once, here, and reused by the data loop below so the
            # attribute and the number can never disagree.
            scales = {}
            for var in out_vars:
                if not is_temperature(var) or is_celsius_exempt(var):
                    continue
                if var == "ta" and derive_ta:
                    scales[var] = "K"       # derived from theta, so kelvin
                    continue
                src_v = nc.variables.get(var)
                scales[var] = temperature_scale(
                    getattr(src_v, "units", None) if src_v is not None else None,
                    var, overrides)
            for var in out_vars:
                ov = fout.createVariable(var, "f4", ("time", "z", "y", "x"),
                                         fill_value=FILL, zlib=True,
                                         complevel=4)
                ov.grid_mapping = gm
                # No `coordinates` attribute: it declares AUXILIARY
                # coordinates in CF, and Esri accepts only netCDF with no
                # auxiliary variables as a multidimensional/voxel source.
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
                if var in scales:
                    ov.units = "degC" if target == "C" else "K"
                    ov.palm2gis_source_units = str(
                        getattr(src, "units", "") if src is not None else "")
                    ov.palm2gis_output_units = ov.units
                elif is_celsius_exempt(var):
                    # Potential temperature stays kelvin whatever `celsius`
                    # says; it is a kelvin-defined quantity.
                    ov.palm2gis_output_units = str(getattr(ov, "units", "K"))

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
                    if var in scales:
                        val = convert_temperature(val, scales[var], target,
                                                  fill=FILL)
                    val[terrain_mask] = FILL
                    fout.variables[var][q, ...] = val
                if derive_ta:
                    th = np.ma.masked_invalid(np.ma.masked_equal(
                        nc.variables["theta"][it, kz0:kzmax,
                                              grid.ny0:grid.nymax,
                                              grid.nx0:grid.nxmax], FILL))
                    ta = np.ma.filled(theta_to_ta(th, pi), FILL)
                    ta = convert_temperature(ta, "K", target, fill=FILL)
                    ta[terrain_mask] = FILL
                    fout.variables["ta"][q, ...] = ta

    info("[voxel] wrote {}: {} x {} x {} x {} (t,z,y,x), {}",
         cfg.disp(file_to), len(its), z_out.size,
         grid.nymax - grid.ny0, grid.nxmax - grid.nx0, ", ".join(out_vars))
