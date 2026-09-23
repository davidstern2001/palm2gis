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

Optionally writes `speed_3d`, the full 3D wind speed sqrt(u^2 + v^2 + w^2)
of the destaggered components (voxel.derive_speed_3d).

By default the k=0 level (z = 0) is dropped so vertical spacing is uniform
(dz/2, 3dz/2, ...), which voxel layers handle best (voxel.drop_k0).
"""

import numpy as np
from netCDF4 import Dataset
from pyproj import CRS

from ..descriptions import long_name
from ..log import debug, get_logger, info, warning
from .georef import (check_product, check_raw_input, require_variables,
                     resolve_z_range, select_times)
from .thermo import absolute_humidity, convert_temperature, \
    exner_basic_state, humidity_scale, is_celsius_exempt, is_temperature, \
    pressure_from_exner, temperature_scale, theta_to_ta

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
def write_crs_var(fout, crs, mode="wkt"):
    """Write a grid_mapping variable ArcGIS and GDAL both understand.

    mode "wkt" writes the CRS as WKT only; "cf" also writes the CF
    grid_mapping parameters (see crs.netcdf_attributes in config.py).
    """
    cv = fout.createVariable("crs", "i4")
    # A CRS given as a PROJ string (or bare WKT) carries no EPSG code, so
    # its WKT has no AUTHORITY and ArcGIS cannot find the WKID. When pyproj
    # identifies it as an EPSG CRS, write that definition instead.
    epsg = crs.to_epsg()
    if epsg:
        crs = CRS.from_epsg(epsg)
    cf = crs.to_cf()

    if mode == "cf":
        # ArcGIS Pro ranks the CF parameters above every WKT attribute
        # (grid_mapping > srid > crs_wkt > EPSG_code > proj4_parms >
        # spatial_ref > esri_pe_string). CF carries no name, EPSG code or
        # datum, so ArcGIS builds a nameless custom CRS on an assumed
        # WGS 1984 from them.
        for k, v in cf.items():
            if k == "crs_wkt":
                continue
            try:
                cv.setncattr(k, v)
            except TypeError:
                pass
        # A CRS pyproj cannot express in CF (Krovak / S-JTSK EPSG:5514,
        # for one) has no CF parameters; ArcGIS falls back to the WKT.
        if "grid_mapping_name" not in cf:
            info("[crs] {} has no CF grid mapping - readers must use the "
                 "WKT (ArcGIS, GDAL and QGIS do).", crs.name)
        # The WGS 1984 assumption only bites when ArcGIS reads CF, so a
        # non-WGS84 datum silently shifts the data. EPSG names the WGS 84
        # datum several ways; none of them should warn.
        datum = (crs.datum.name if crs.datum is not None else "")
        _d = datum.lower().replace("_", " ")
        is_wgs84 = "wgs 84" in _d or "world geodetic system 1984" in _d
        if datum and not is_wgs84 and "grid_mapping_name" in cf:
            # ArcGIS assumes WGS 1984 for CF grid mappings.
            warning("[crs] datum '{}' is not WGS 84 - ArcGIS Pro will shift "
                    "it; set crs.netcdf_attributes: wkt.", datum)

    # Every WKT is written as ASCII WKT1, which netCDF stores as a classic
    # char attribute. pyproj's default WKT2 carries a USAGE block with a
    # non-ASCII degree sign, so netCDF4 stored it as NC_STRING; together
    # with the WKT2-only ENSEMBLE node that left ArcGIS Pro parsing
    # spatial_ref into an unknown CRS (WKID 0) and reporting it unsupported.
    # WKT1_GDAL keeps AUTHORITY["EPSG", ...], from which ArcGIS gets the
    # WKID.
    wkt_esri = _ascii_wkt(crs, "WKT1_ESRI")
    wkt_gdal = _ascii_wkt(crs, "WKT1_GDAL") or wkt_esri
    cv.crs_wkt = wkt_gdal
    cv.spatial_ref = wkt_gdal
    cv.esri_pe_string = wkt_esri
    if epsg:
        cv.EPSG_code = f"EPSG:{epsg}"
    # The voxel layer reads esri_pe_string as a GLOBAL attribute; on the
    # grid_mapping variable alone it is not seen.
    fout.setncattr("esri_pe_string", wkt_esri)
    return "crs"


def _ascii_wkt(crs, version):
    """WKT of *crs* in *version*, ASCII only (stored as NC_CHAR)."""
    try:
        wkt = crs.to_wkt(version=version)
    except Exception:
        return None
    if not wkt:
        return None
    return wkt.encode("ascii", "ignore").decode("ascii")


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
        #   z_rel - height above the domain base, what PALM wrote. PALM's
        #           basic-state Exner function below pairs it with
        #           surface_pressure, the pressure at the model's z = 0.
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
        # Absolute humidity, ensured the same way: PALM's own `ah` if the
        # file has one, otherwise derived from `q`. Its temperature comes
        # from PALM's `ta` when present (exact) and from theta otherwise.
        derive_ah = False
        ah_temp = None
        if vx["derive_ah"] and "ah" not in scalars:
            if "ah" in nc.variables:
                info("[voxel] the file holds 'ah' - using PALM's own "
                     "absolute humidity instead of deriving it.")
                scalars = scalars + ["ah"]
            elif "q" not in nc.variables:
                require_variables(
                    nc, ["q"], "voxel", vx["file_from"],
                    hint="Absolute humidity is derived from q. Add 'q' to "
                         "data_output in the p3d (needs humidity = .TRUE.), "
                         "or drop ah from voxel.variables.")
            elif "ta" in nc.variables or "theta" in nc.variables:
                derive_ah = True
                ah_temp = "ta" if "ta" in nc.variables else "theta"
                info("[voxel] no 'ah' in the file - deriving it from 'q', "
                     "with the air temperature from '{}'.", ah_temp)
            else:
                require_variables(
                    nc, ["ta"], "voxel", vx["file_from"],
                    hint="Absolute humidity needs an air temperature: "
                         "PALM's ta, or theta to derive it from.")
        do_velocity = bool(vx["do_velocity"])
        derive_speed = bool(vx["derive_speed_3d"])
        need_uvw = do_velocity or derive_speed
        if do_velocity:
            require_variables(
                nc, ["u", "v", "w"], "voxel", vx["file_from"],
                hint="Set voxel.do_velocity: false to write scalars only, "
                     "or point file_from at the joined file.")
        if derive_speed:
            require_variables(
                nc, ["u", "v", "w"], "voxel", vx["file_from"],
                hint="speed_3d is computed from u, v and w. Point file_from "
                     "at a file holding all three, or drop speed_3d.")
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
        if need_uvw and "zw_3d" in nc.variables:
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
        if derive_ah:
            out_vars = out_vars + ["ah"]
        if derive_speed:
            out_vars = out_vars + ["speed_3d"]

        with cfg.atomic_output(file_to, "voxel") as _tmp_out, \
                Dataset(_tmp_out, "w", format="NETCDF4") as fout:
            fout.Conventions = "CF-1.8"
            fout.title = "PALM 3D output for ArcGIS Pro voxel layer"
            fout.source = str(file_from)

            # Horizontal axes in the output CRS: PALM's own grid, or a
            # regular grid resampled from it (see reproject).
            og = grid.reproj.grid()
            yd, xd = og.dims
            og.write_axes(fout)
            fout.createDimension("z", z_out.size)
            fout.createDimension("time", len(its))

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

            gm = write_crs_var(fout, grid.reproj.crs,
                               cfg["crs"]["netcdf_attributes"])

            # PALM's basic-state Exner function on the output levels: linear
            # in z, fixed for the whole run (see thermo). One profile serves
            # both the ta and the ah derivation.
            pi = None
            if derive_ta or derive_ah:
                p_s, th_ref = cfg["surface_pressure"], cfg["pt_surface"]
                pi_col = exner_basic_state(z_rel, p_s, th_ref)
                pi = pi_col[:, None, None]
                debug("[voxel] Exner {:.5f} .. {:.5f} over {:.0f} .. {:.0f} "
                      "m asl ({:.2f} hPa, {:.2f} K)", float(pi_col[0]),
                      float(pi_col[-1]), float(z_out[0]), float(z_out[-1]),
                      p_s / 100.0, th_ref)

            p_col = q_factor = ta_scale = None
            if derive_ah:
                # Absolute pressure on the output levels, (nz, 1, 1), from
                # the same basic state ta is derived from.
                p_col = pressure_from_exner(pi)
                q_units = getattr(nc.variables["q"], "units", None)
                q_factor = humidity_scale(q_units)
                if ah_temp == "ta":
                    ta_scale = temperature_scale(
                        getattr(nc.variables["ta"], "units", None), "ta",
                        dict(cfg["temperature_units"] or {}))
                info("[voxel] deriving ah from q ({}): p {:.0f} .. {:.0f} Pa "
                     "over {:.0f} .. {:.0f} m",
                     q_units or "no units, taken as kg/kg",
                     float(p_col.ravel()[0]), float(p_col.ravel()[-1]),
                     float(z_out[0]), float(z_out[-1]))

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
                ov = fout.createVariable(var, "f4", ("time", "z", yd, xd),
                                         fill_value=FILL, zlib=True,
                                         complevel=4)
                ov.grid_mapping = gm
                # No `coordinates` attribute: it declares AUXILIARY
                # coordinates in CF, and Esri accepts only netCDF with no
                # auxiliary variables as a multidimensional/voxel source.
                if var == "ta" and derive_ta:
                    ov.units = "degC" if celsius else "K"
                    ov.long_name = long_name("ta")
                    ov.comment = ("derived from theta with PALM's basic-state "
                                  "Exner function (PALM did not output ta)")
                    ov.derived_from = "theta"
                    ov.palm2gis_surface_pressure = float(
                        cfg["surface_pressure"])
                    ov.palm2gis_pt_surface = float(cfg["pt_surface"])
                    continue
                if var == "speed_3d" and derive_speed:
                    # No CF standard_name: CF's wind_speed is horizontal.
                    ov.units = "m s-1"
                    ov.long_name = long_name("speed_3d")
                    ov.comment = ("sqrt(u^2 + v^2 + w^2) of u, v, w "
                                  "interpolated to cell centres")
                    ov.derived_from = "u v w"
                    continue
                if var == "ah" and derive_ah:
                    ov.units = "g m-3"
                    ov.long_name = long_name("ah")
                    ov.comment = (
                        "derived as q * p / (Rd * Tv), Tv = T (1 + 0.608 q), "
                        "with T from {} and p from PALM's basic-state "
                        "Exner function (PALM did not output ah). PALM's q "
                        "is the total water content, which equals specific "
                        "humidity when no cloud water is present."
                        .format("PALM's ta" if ah_temp == "ta"
                                else "theta"))
                    ov.derived_from = "q"
                    ov.palm2gis_surface_pressure = float(
                        cfg["surface_pressure"])
                    ov.palm2gis_pt_surface = float(cfg["pt_surface"])
                    continue
                src = nc.variables.get(var)
                if src is not None and hasattr(src, "units"):
                    ov.units = src.units
                # ArcGIS Pro shows long_name as the variable Description.
                # PALM writes only the bare name there ("theta"), so it is
                # replaced by a real description, falling back to PALM's
                # own text for variables the table does not list.
                ov.long_name = long_name(
                    var, getattr(src, "long_name", None)
                    if src is not None else None)
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
            mask_src = next((c for c in scalars + ["w", "theta", "q", "ta"]
                             if c in nc.variables), None)
            terrain_mask = np.ma.getmaskarray(
                nc.variables[mask_src][its[0], kz0:kzmax,
                                       grid.ny0:grid.nymax,
                                       grid.nx0:grid.nxmax])
            debug("[voxel] terrain mask from '{}': {} cells",
                  mask_src, int(terrain_mask.sum()))
            # The same mask on the output grid; cells outside the PALM
            # domain count as masked too. `air` tells the conservative
            # resampling which PALM cells to average over.
            mask_out = og.mask(terrain_mask)
            air = ~terrain_mask

            for q, it in enumerate(its):
                debug("[voxel] time {} ({:.0f} s)", it, times[it])
                if need_uvw:
                    # Destaggered once, shared by u/v/w and speed_3d.
                    comp = {var: destagger(nc, var, it, kz0, kzmax,
                                           grid.ny0, grid.nymax,
                                           grid.nx0, grid.nxmax)
                            for var in ("u", "v", "w")}
                    if derive_speed:
                        sp = np.sqrt(comp["u"] ** 2 + comp["v"] ** 2
                                     + comp["w"] ** 2)
                        sp = np.array(og.take(sp, valid=air), copy=True)
                        sp[mask_out] = FILL
                        fout.variables["speed_3d"][q, ...] = sp
                    if do_velocity:
                        # Taken onto the output grid and turned from PALM's
                        # grid axes into the output CRS's east/north before
                        # any NoData goes in, so a fill value is never
                        # mixed into a rotated component.
                        out = {var: og.take(comp[var], np.nan, valid=air)
                               for var in ("u", "v", "w")}
                        out["u"], out["v"] = og.rotate(out["u"], out["v"])
                        for var in ("u", "v", "w"):
                            val = np.array(out[var], copy=True)
                            val[mask_out | ~np.isfinite(val)] = FILL
                            fout.variables[var][q, ...] = val
                for var in scalars:
                    val = np.ma.filled(
                        nc.variables[var][it, kz0:kzmax,
                                          grid.ny0:grid.nymax,
                                          grid.nx0:grid.nxmax], FILL)
                    if var in scales:
                        val = convert_temperature(val, scales[var], target,
                                                  fill=FILL)
                    val = np.array(og.take(val, valid=air & (val != FILL)),
                                   copy=True)
                    val[mask_out] = FILL
                    fout.variables[var][q, ...] = val
                if derive_ta:
                    th = np.ma.masked_invalid(np.ma.masked_equal(
                        nc.variables["theta"][it, kz0:kzmax,
                                              grid.ny0:grid.nymax,
                                              grid.nx0:grid.nxmax], FILL))
                    ta = np.ma.filled(theta_to_ta(th, pi), FILL)
                    ta = convert_temperature(ta, "K", target, fill=FILL)
                    ta = np.array(og.take(ta, valid=air & (ta != FILL)),
                                  copy=True)
                    ta[mask_out] = FILL
                    fout.variables["ta"][q, ...] = ta
                if derive_ah:
                    fields = {}
                    for name in ("q", ah_temp):
                        fields[name] = np.ma.masked_invalid(
                            np.ma.masked_equal(
                                nc.variables[name][it, kz0:kzmax,
                                                   grid.ny0:grid.nymax,
                                                   grid.nx0:grid.nxmax],
                                FILL)).astype("f8")
                    qq = fields["q"] * q_factor
                    if ah_temp == "theta":
                        t_k = theta_to_ta(fields["theta"], pi)
                    else:
                        t_raw = np.ma.filled(fields["ta"], np.nan)
                        t_k = np.ma.masked_invalid(convert_temperature(
                            t_raw, ta_scale, "K"))
                    if q == 0:
                        q_max = float(qq.max()) if qq.count() else 0.0
                        if q_max > 0.05:
                            warning("[voxel] q reaches {:.3f} kg/kg, far "
                                    "above any real humidity - check its "
                                    "units ({}).", q_max,
                                    getattr(nc.variables["q"], "units",
                                            "none given"))
                    ah = np.ma.filled(absolute_humidity(qq, t_k, p_col),
                                      FILL)
                    ah = np.array(og.take(ah, valid=air & (ah != FILL)),
                                  copy=True)
                    ah[mask_out] = FILL
                    fout.variables["ah"][q, ...] = ah

    info("[voxel] wrote {}: {} x {} x {} x {} (t,z,{},{}), {}",
         cfg.disp(file_to), len(its), z_out.size, og.shape[0], og.shape[1],
         yd, xd, ", ".join(out_vars))
