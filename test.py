#!/usr/bin/env python3
"""palm2gis test.

Builds a tiny synthetic PALM-like dataset in a temporary directory (a
static driver, a staggered 3D output and an unstructured surface file),
runs the full pipeline over it, and checks that the expected outputs
appear with the right CRS, georeferencing and values. Nothing outside the
temporary directory is touched, and no real data is needed.

Run it after installing the dependencies, before pointing the pipeline at
real data:

    python test.py            # run, report, clean up
    python test.py --keep     # keep the temporary directory for inspection

Exit status is 0 on success, 1 on failure, so it also works in CI.
"""

# ------------------------------
# 1. TEST DOMAIN DEFINITION
# ------------------------------
import argparse
import shutil
import sys
import tempfile
import traceback
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

# Synthetic domain: 20 x 20 x 12 points, dx = dy = dz = 2 m, gently sloping
# terrain, one 8 m block building, three timesteps 10 minutes apart.
# UTM zone 33N given as an explicit PROJ definition rather than the bare
# code "EPSG:32633" on purpose: creating data from a bare EPSG code forces
# a lookup in PROJ's proj.db, which fails on environments with a stale or
# mismatched PROJ installation. The pipeline reads crs.palm from the
# config, where the test also uses the explicit definition.
CRS = "+proj=utm +zone=33 +datum=WGS84 +units=m +no_defs"
NX, NY, NZ = 20, 20, 12
DX = DY = DZ = 2.0
ORIGIN_X, ORIGIN_Y = 458000.0, 5546000.0
# Deliberately NOT zero: every z in every product is relative to this, and
# a 0 here is what let the vertical misalignment ship unnoticed.
ORIGIN_Z = 175.0
ORIGIN_TIME = "2026-07-15 10:00:00 +02"
NT = 3


# ------------------------------
# 2. DEPENDENCY CHECK
# ------------------------------
def check_dependencies():
    """Import every third-party dependency and report versions, so a broken
    environment fails here with a clear message instead of deep inside a
    stage."""
    missing, versions = [], {}
    for mod in ("numpy", "netCDF4", "yaml", "pyproj", "rasterio",
                "geopandas", "shapely"):
        try:
            m = __import__(mod)
            versions[mod] = getattr(m, "__version__", "?")
        except ImportError:
            missing.append(mod)
    if missing:
        print("Missing dependencies: {} - run "
              "`pip install -r requirements.txt`.".format(", ".join(missing)))
        return False
    print("Dependencies OK: " + ", ".join(
        "{} {}".format(k, v) for k, v in versions.items()))
    return True


# ------------------------------
# 3. SYNTHETIC DATASET
# ------------------------------
def build_dataset(root):
    import numpy as np
    from netCDF4 import Dataset

    x = np.arange(DX / 2.0, NX * DX, DX)
    y = np.arange(DY / 2.0, NY * DY, DY)
    zu = np.append(0.0, np.arange(DZ / 2.0, (NZ - 1) * DZ, DZ))
    zw = np.arange(0.0, NZ * DZ, DZ)[:NZ]

    zt = np.round(np.tile(np.linspace(0.0, 4.0, NX), (NY, 1)) / DZ) * DZ
    bld = np.zeros((NY, NX))
    bld[8:12, 8:12] = 8.0
    elev = zt + bld

    # -- static driver -----------------------------------------------------
    with Dataset(root / "test_static", "w", format="NETCDF4") as nc:
        nc.origin_x = ORIGIN_X
        nc.origin_y = ORIGIN_Y
        nc.origin_z = ORIGIN_Z
        nc.origin_time = ORIGIN_TIME
        nc.createDimension("x", NX)
        nc.createDimension("y", NY)
        nc.createVariable("x", "f4", ("x",))[:] = x
        nc.createVariable("y", "f4", ("y",))[:] = y
        nc.createVariable("zt", "f4", ("y", "x"))[:] = zt
        b = nc.createVariable("buildings_2d", "f4", ("y", "x"),
                              fill_value=-9999.0)
        b[:] = np.ma.masked_where(bld == 0.0, bld)

        # Plant canopy: a 3-level tree block on flat-ish ground, stored on
        # zlad = height ABOVE TERRAIN (as PALM does), away from the
        # building so the two do not overlap.
        nzlad = 3
        nc.createDimension("zlad", nzlad)
        nc.createVariable("zlad", "f4", ("zlad",))[:] = \
            np.arange(DZ / 2.0, nzlad * DZ, DZ)[:nzlad]
        lad = np.zeros((nzlad, NY, NX), dtype="f4")
        lad[:, 3:6, 3:6] = 0.5           # 3 x 3 columns, 3 levels deep
        lv = nc.createVariable("lad", "f4", ("zlad", "y", "x"),
                               fill_value=-9999.0)
        lv[:] = np.ma.masked_where(lad == 0.0, lad)

    # -- 3D output (staggered u/v/w + theta, terrain masked) ---------------
    kmask = np.zeros((NZ, NY, NX), dtype=bool)
    for j in range(NY):
        for i in range(NX):
            kmask[:, j, i] = zu < elev[j, i]

    rng = np.random.default_rng(1)
    with Dataset(root / "test_3d", "w", format="NETCDF4") as nc:
        # Real PALM 3D-output convention: the staggered dims xu/yv/zw_3d
        # have the SAME length as x/y/zu_3d (the outermost face is not
        # stored). This is exactly the case that clips the destagger slice
        # when the cut reaches the domain edge.
        nc.createDimension("time", NT)
        nc.createDimension("x", NX)
        nc.createDimension("xu", NX)
        nc.createDimension("y", NY)
        nc.createDimension("yv", NY)
        nc.createDimension("zu_3d", NZ)
        nc.createDimension("zw_3d", NZ)
        nc.createVariable("time", "f8", ("time",))[:] = [0.0, 600.0, 1200.0]
        nc.createVariable("x", "f4", ("x",))[:] = x
        nc.createVariable("xu", "f4", ("xu",))[:] = np.arange(
            0, NX * DX, DX)
        nc.createVariable("y", "f4", ("y",))[:] = y
        nc.createVariable("yv", "f4", ("yv",))[:] = np.arange(
            0, NY * DY, DY)
        nc.createVariable("zu_3d", "f4", ("zu_3d",))[:] = zu
        nc.createVariable("zw_3d", "f4", ("zw_3d",))[:] = zw

        def mk(name, dims, base, units=None):
            v = nc.createVariable(name, "f4", dims, fill_value=-9999.0)
            if units is not None:
                v.units = units
            for t in range(NT):
                shp = tuple(len(nc.dimensions[d]) for d in dims[1:])
                data = base + 0.1 * t + 0.05 * rng.standard_normal(shp)
                m = np.zeros(shp, dtype=bool)
                m[:NZ, :NY, :NX] = kmask[:shp[0], :shp[1], :shp[2]]
                v[t] = np.ma.masked_where(m, data)

        mk("u", ("time", "zu_3d", "y", "xu"), 2.0)
        mk("v", ("time", "zu_3d", "yv", "x"), 1.0)
        mk("w", ("time", "zw_3d", "y", "x"), 0.1)
        mk("theta", ("time", "zu_3d", "y", "x"), 300.0, units="K")

        # q in kg/kg: its own values, since mk's +0.1 per record and 0.05
        # noise would give negative or impossible humidities
        qv = nc.createVariable("q", "f4", ("time", "zu_3d", "y", "x"),
                               fill_value=-9999.0)
        qv.units = "kg/kg"
        for t in range(NT):
            shp = (NZ, NY, NX)
            data = 0.008 + 0.0005 * t + 0.0002 * rng.standard_normal(shp)
            qv[t] = np.ma.masked_where(kmask[:NZ, :NY, :NX], data)

    # -- surface output (unstructured: tops + two building faces) ----------
    xs, ys, zs, zen, az = [], [], [], [], []
    for j in range(NY):
        for i in range(NX):
            xs.append(x[i]); ys.append(y[j]); zs.append(elev[j, i])
            zen.append(0.0); az.append(-9999.0)
    for i in range(8, 12):
        for k in range(int(zt[8, i] / DZ), int(elev[8, i] / DZ)):
            zc = zw[k] + DZ / 2.0
            xs.append(x[i]); ys.append(8 * DY); zs.append(zc)
            zen.append(90.0); az.append(180.0)
            xs.append(x[i]); ys.append(12 * DY); zs.append(zc)
            zen.append(90.0); az.append(0.0)

    with Dataset(root / "test_surf", "w", format="NETCDF4") as nc:
        nc.createDimension("time", NT)
        nc.createDimension("s", len(xs))
        nc.createVariable("time", "f8", ("time",))[:] = [0.0, 600.0, 1200.0]
        nc.createVariable("xs", "f4", ("s",))[:] = xs
        nc.createVariable("ys", "f4", ("s",))[:] = ys
        nc.createVariable("zs", "f4", ("s",))[:] = zs
        nc.createVariable("zenith", "f4", ("s",))[:] = zen
        nc.createVariable("azimuth", "f4", ("s",),
                          fill_value=-9999.0)[:] = az
        import numpy as np
        v = nc.createVariable("theta_surface", "f4", ("time", "s"))
        for t in range(NT):
            v[t] = 295.0 + 5.0 * (np.asarray(zen) == 0.0) + 0.5 * t

    # -- config ------------------------------------------------------------
    (root / "config.yaml").write_text(f"""\
project:
  name: selftest
  output_dir: ./out
crs:
  palm: "{CRS}"
domain:
  static_driver: test_static
  nz: {NZ}
  dz: {DZ}
bounds:
  zmax: 10
voxel:
  enabled: true
  file_from: test_3d
  file_to: test_voxel.nc
  scalars: ["theta"]
  derive_ta: true
surfaces:
  enabled: true
  file_from: test_surf
  prefix: test_
  variables: ["theta_surface"]
  derive_t_surface: true
  write_all_elements_gpkg: true
canopy:
  enabled: true
  prefix: test_
  write_volume: true
wind:
  enabled: true
  file_from: test_3d
  file_to: test_wind.gpkg
  heights_agl: [2.0, 10.0]
  stride: 2
  grid_to: test_wind.nc
  grid_heights: [180.0, 183.0]
""")


# ------------------------------
# 4. OUTPUT VERIFICATION
# ------------------------------
def check(cond, msg, failures):
    tag = "ok  " if cond else "FAIL"
    print(f"  [{tag}] {msg}")
    if not cond:
        failures.append(msg)


_GRID_FOR_CHECKS = None


def verify_outputs(out):
    import numpy as _np_thermo
    import numpy as np
    import geopandas as gpd
    from netCDF4 import Dataset

    failures = []

    # -- voxel -------------------------------------------------------------
    with Dataset(out / "test_voxel.nc") as nc:
        nc.set_auto_mask(False)
        z = np.asarray(nc["z"][:])
        check(np.allclose(np.diff(z), DZ), "voxel: uniform z spacing "
              "(drop_k0)", failures)
        check(float(nc["x"][0]) == ORIGIN_X + DX / 2.0,
              "voxel: absolute x coordinates", failures)
        check("esri_pe_string" in nc["crs"].ncattrs(),
              "voxel: Esri PE string on the crs variable", failures)
        check("esri_pe_string" in nc.ncattrs(),
              "voxel: global Esri PE string (read by the voxel layer)",
              failures)
        check("grid_mapping_name" not in nc["crs"].ncattrs()
              and "AUTHORITY[\"EPSG\",\"32633\"]" in nc["crs"].crs_wkt
              and nc["crs"].EPSG_code == "EPSG:32633",
              "voxel: crs holds WKT with EPSG authority, no CF params "
              "(ArcGIS then reads the WKID)", failures)
        check(all(str(nc["crs"].getncattr(a)).isascii()
                  and "ENSEMBLE" not in str(nc["crs"].getncattr(a))
                  for a in ("spatial_ref", "crs_wkt", "esri_pe_string")),
              "voxel: CRS WKT attributes are ASCII WKT1", failures)
        check(nc["time"].units.startswith("seconds since 2026-07-15"),
              "voxel: CF time from origin_time", failures)
        u = nc["u"][0]
        check(int((u == -9999.0).sum()) > 0,
              "voxel: terrain/building cells are NoData", failures)
        check(abs(float(u[u != -9999.0].mean()) - 2.0) < 0.1,
              "voxel: destaggered u close to the prescribed 2.0 m/s",
              failures)
        # actual_range spans EVERY record, not the first one: an ArcGIS
        # voxel layer's own statistics come from the loaded time step.
        ok_rng = True
        for var in ("u", "theta"):
            a = np.asarray(nc[var][:])
            a = a[a != -9999.0]
            rng = np.asarray(nc[var].actual_range)
            ok_rng &= bool(np.allclose(rng, [a.min(), a.max()]))
        t0 = np.asarray(nc["theta"][0]); t0 = t0[t0 != -9999.0]
        # Display range = p2 .. p98 of every value over every record,
        # exact to the fine-histogram bin width.
        th = np.asarray(nc["theta"][:]); th = th[th != -9999.0]
        p2, p98 = np.percentile(th, [2, 98])
        dr = np.asarray(nc["theta"].palm2gis_display_range, dtype="f8")
        tol = 2 * (th.max() - th.min()) / 4000
        check(abs(dr[0] - p2) <= tol and abs(dr[1] - p98) <= tol,
              "voxel: palm2gis_display_range is p2 .. p98 over all records",
              failures)
        check(ok_rng and nc["time"].size > 1
              and float(nc["theta"].actual_range[1]) > float(t0.max()),
              "voxel: actual_range covers all records (wider than the "
              "first)", failures)

    import csv as _csv
    with open(out / "test_voxel_ranges.csv", newline="") as f:
        rows = {r["variable"]: r for r in _csv.DictReader(f)}
    check(set(rows) >= {"u", "v", "w", "theta"}
          and all(float(r["min"]) <= float(r["suggested_min"])
                  <= float(r["suggested_max"]) <= float(r["max"])
                  for r in rows.values())
          and (out / "test_voxel_histogram.csv").exists(),
          "voxel: ranges and histogram CSVs written, min <= p2 <= p98 <= max",
          failures)

    # -- derived air temperature -------------------------------------------
    # PALM's basic state: Pi linear in z, (p_s/p0)^kappa at z = 0, slope
    # -g/(cp theta_ref), so T of a constant-theta_ref column falls at g/cp.
    from palm2gis.steps.thermo import (CP, G, exner_basic_state,
                                       fit_basic_state, theta_to_ta)
    zz = _np_thermo.arange(0.0, 1001.0, 50.0)
    pi_bs = exner_basic_state(zz, 1.0e5, 300.0)
    T = theta_to_ta(_np_thermo.full_like(zz, 300.0), pi_bs)
    lapse = (T[0] - T[-1]) / (zz[-1] - zz[0])
    check(abs(lapse - G / CP) < 1e-9,
          f"thermo: dry adiabatic lapse rate {lapse * 1000:.3f} K/km "
          f"= g/cp", failures)
    check(abs(pi_bs[0] - 1.0) < 1e-12,
          "thermo: Exner = 1 at the reference pressure", failures)
    # anchored at z = 0: a cut starting at 100 m gives the full column's
    # values there, not the surface value again
    cut = _np_thermo.arange(100.0, 1001.0, 50.0)
    check(abs(exner_basic_state(cut, 1.0e5, 300.0)[0] - pi_bs[2]) < 1e-12,
          "thermo: Exner anchored at z = 0, not at the cut base", failures)
    # the fit recovers the basic state of the real holesovice parent run
    pi_real = exner_basic_state(zz, 99388.42, 285.0)
    fit = fit_basic_state(zz, pi_real)
    check(fit is not None and abs(fit[0] - 99388.42) < 0.01
          and abs(fit[1] - 285.0) < 1e-6 and fit[2] < 1e-9,
          "thermo: fit_basic_state recovers p_s and pt_surface", failures)
    check(fit_basic_state(zz, _np_thermo.full_like(zz, 1.0)) is None,
          "thermo: a flat Exner profile is not accepted as a basic state",
          failures)

    # derived temperature must reach the surface element layer
    allel_ta = gpd.read_file(out / "test_surfaces.gpkg", layer="surfaces")
    # t_surface is written in degrees C by default, theta_surface in K.
    check("t_surface" in allel_ta.columns and "t_surf" not in allel_ta.columns
          and (allel_ta["t_surface"] + 273.15
               > allel_ta["theta_surface"]).all(),
          "thermo: t_surface in the element layer (p_s > p0 -> T > theta)",
          failures)
    with Dataset(out / "test_voxel.nc") as nc:
        check("ta" in nc.variables
              and getattr(nc["ta"], "derived_from", "") == "theta",
              "thermo: derived ta in the voxel NetCDF, marked as derived",
              failures)

    # PALM's own variables must WIN over derivation
    from palm2gis.steps.surfaces import output_vars
    sf_native = {"variables": ["theta_surface"], "derive_t_surface": True}
    plan = output_vars(sf_native, {"theta_surface", "t_surface"})
    check(("t_surface", "t_surface", False) in plan,
          "thermo: native t_surface preferred over deriving it", failures)
    plan = output_vars(sf_native, {"theta_surface", "t_surf"})
    check(("t_surface", "t_surf", False) in plan,
          "thermo: a native t_surf is copied and written as t_surface",
          failures)
    plan = output_vars(sf_native, {"theta_surface"})
    check(("t_surface", "theta_surface", True) in plan,
          "thermo: t_surface derived when the file lacks it", failures)
    plan = output_vars({"variables": ["theta_surf"],
                        "derive_t_surface": True}, {"theta_surf"})
    check(("t_surface", "theta_surf", True) in plan,
          "thermo: the theta_surf spelling is accepted too", failures)

    # -- surfaces ----------------------------------------------------------
    # the removed products must stay removed
    gone = ["test_theta_surface_topview.tif", "test_t_surf_topview.tif",
            "test_surface_2d.nc", "test_walls.gpkg", "test_lai.tif",
            "test_canopy_height.tif"]
    check(not any((out / g).exists() for g in gone),
          "removed outputs: no top-view, walls-only, LAI or CHM files",
          failures)

    # every surface element: 400 ground/roof faces + 32 wall faces
    allel = gpd.read_file(out / "test_surfaces.gpkg", layer="surfaces")
    per_time = len(allel) // NT
    check(per_time == 432,
          f"surfaces: all elements present ({per_time} faces per timestep)",
          failures)
    check(len(allel) == 432 * NT,
          f"surfaces: exactly {432 * NT} rows after per-timestep append, "
          f"got {len(allel)}", failures)
    check(sorted(allel["rec"].unique()) == sorted(allel["rec"].unique())
          and allel.groupby("rec").size().nunique() == 1,
          "surfaces: no timestep duplicated or dropped by the append",
          failures)
    kinds = set(allel["orient"].unique())
    check(kinds == {"up", "wall_N", "wall_S"},
          f"surfaces: face orientations labelled ({sorted(kinds)})",
          failures)
    roof = allel[(allel.orient == "up") & (allel.z_agl > 1.0)]
    check(len(roof) // NT == 16 and abs(roof.z_agl.max() - 8.0) < 1e-6,
          "surfaces: z_agl separates roofs (8 m) from ground", failures)
    check(allel.geometry.iloc[0].has_z and allel["theta_surface"].notna().all(),
          "surfaces: all-element layer is PolygonZ with values", failures)
    up0 = allel[(allel.orient == "up") & (allel.rec == allel.rec.min())]
    check(np.allclose(up0["theta_surface"], 300.0),
          "surfaces: upward theta_surface = 300 at the first record",
          failures)
    walls = allel[allel.orient.str.startswith("wall_")]
    check(len(walls) == 32 * NT and walls["time"].nunique() == NT,
          f"surfaces: {32 * NT} wall rows over {NT} timesteps", failures)
    import pyogrio as _pyo
    check([l[0] for l in _pyo.list_layers(out / "test_surfaces.gpkg")
           if l[0] != "layer_styles"] == ["surfaces"],
          "surfaces: GeoPackage holds only the element layer (no _2d)",
          failures)

    # -- QGIS style sidecars ------------------------------------------------
    import sqlite3
    import xml.dom.minidom as minidom

    vqml = out / "test_wind.qml"
    check(vqml.exists(), "qgis: vector style sidecar written", failures)
    if vqml.exists():
        t = minidom.parse(str(vqml)).getElementsByTagName("temporal")[0]
        check(t.getAttribute("mode") == "1"
              and t.getAttribute("startField") == "time",
              "qgis: vector QML reads time from the time field", failures)

    con = sqlite3.connect(out / "test_wind.gpkg")
    rows = con.execute("SELECT f_table_name, useAsDefault FROM "
                       "layer_styles").fetchall()
    con.close()
    check(rows and rows[0][0] == "wind" and rows[0][1] == 1,
          "qgis: style embedded in the GeoPackage as default", failures)

    # -- p3d namelist parsing ----------------------------------------------
    # PALM packs several assignments onto one line, and the header comment
    # repeats the grid in prose - both must be handled.
    from palm2gis.steps.georef import infer_time_params, read_p3d
    nml = out / "_p3d_probe"
    nml.write_text(
        "!-- Grid:   512 x 608 x 512 cells  (nx=511, ny=607, nz=999)\n"
        "!-- Spacing: dx=dy=dz=9m\n"
        "&initialization_parameters\n"
        "         nx = 511, ny = 607, nz = 512,\n"
        "         dx = 1.0, dy = 1.0, dz = 1.0,\n"
        "/\n"
        "&runtime_parameters\n"
        "    dt_data_output_av = 600.0,\n"
        "    dt_dosurf_av      = 3600.0,\n"
        "    skip_time_dosurf_av = 0.0,\n"
        "/\n")
    p3d = read_p3d(nml)
    check(p3d.get("nz") == 512 and p3d.get("dz") == 1.0,
          "p3d: several assignments on one line are all read", failures)
    check(p3d.get("nx") == 511 and p3d.get("dx") == 1.0,
          "p3d: grid settings parsed from a packed line", failures)

    class _G:
        pass

    gg = _G()
    gg.p3d = p3d
    check(infer_time_params(gg, "voxel", "x_av_3d.nc")[0] == 600.0
          and infer_time_params(gg, "surfaces", "x_av_surf.nc")[0] == 3600.0,
          "p3d: 3D and surface intervals inferred separately", failures)
    nml.unlink()

    # -- canopy (LAD) ------------------------------------------------------
    with Dataset(out / "test_lad_3d.nc") as nc:
        nc.set_auto_mask(False)
        check("coordinates" not in nc["lad"].ncattrs(),
              "canopy: no auxiliary `coordinates` attribute (ArcGIS)",
              failures)
        lad = nc["lad"][:]
        check(int((lad > 0).sum()) == 27,
              "canopy: LAD volume holds every canopy cell", failures)
        zc = _np_thermo.asarray(nc["z"][:])
        check(_np_thermo.allclose(_np_thermo.diff(zc), DZ),
              "canopy: z axis regular (k = 0 dropped) - an ArcGIS voxel "
              "layer needs it", failures)
        # the test terrain slopes across the tree block, so the SAME
        # canopy must land on different absolute levels per column -
        # which is exactly what the height-above-terrain lift is for
        # lad is (time, z, y, x): the canopy volume carries a length-1
        # time axis so it stays visible while a scene's time slider
        # moves. Reduce over every axis except z rather than assuming a
        # 3D array, which is what this used to do - and it raised
        # IndexError, so the self-test could not certify anything.
        z_axis = nc["lad"].dimensions.index("z")
        other = tuple(a for a in range(lad.ndim) if a != z_axis)
        zocc = _np_thermo.asarray(nc["z"][:])[(lad > 0).any(axis=other)]
        check(zocc.size > 3,
              f"canopy: terrain lift spreads the canopy over {zocc.size} "
              f"absolute levels", failures)

    # -- wind --------------------------------------------------------------
    wind = gpd.read_file(out / "test_wind.gpkg", layer="wind")
    check(set(wind["height_agl"].unique()) == {2.0, 10.0},
          "wind: both AGL heights present", failures)
    check(wind.groupby("rec").size().nunique() == 1
          and wind["rec"].nunique() == NT,
          f"wind: {NT} timesteps appended cleanly, equal rows each",
          failures)
    # u = 2, v = 1 -> meteorological direction 270 - atan2(1, 2) = 243.4
    dmean = wind[wind.height_agl == 10.0].dir_met.mean()
    check(abs(dmean - 243.4) < 3.0,
          f"wind: mean direction {dmean:.1f} deg close to 243.4",
          failures)

    return failures


# ------------------------------
# 4b. CONFIG-CHANGE DETECTION
# ------------------------------
def verify_state(tmp, pipeline_main):
    """Editing a setting must invalidate the outputs it affects, and only
    those. Skipping on filename alone used to keep stale files silently."""
    import time as _time
    out = tmp / "out"
    cfgfile = tmp / "config.yaml"
    failures = []

    wind = out / "test_wind.gpkg"
    voxel = out / "test_voxel.nc"

    # settle, so mtime comparisons are meaningful
    pipeline_main(["-c", str(cfgfile)])
    _time.sleep(0.01)
    s_wind = wind.stat().st_mtime_ns
    s_voxel = voxel.stat().st_mtime_ns

    check((out / ".palm2gis_state.json").is_file(),
          "state: sidecar written to the output directory", failures)

    # a rerun with an UNCHANGED config still keeps everything
    pipeline_main(["-c", str(cfgfile)])
    check(wind.stat().st_mtime_ns == s_wind,
          "state: unchanged config still keeps existing outputs", failures)

    # change ONLY a wind setting
    text = cfgfile.read_text()
    assert "heights_agl: [2.0, 10.0]" in text, "test fixture changed"
    cfgfile.write_text(text.replace("heights_agl: [2.0, 10.0]",
                                    "heights_agl: [2.0, 10.0, 20.0]"))
    pipeline_main(["-c", str(cfgfile)])
    check(wind.stat().st_mtime_ns != s_wind,
          "state: editing wind.heights rewrites the wind output", failures)
    check(voxel.stat().st_mtime_ns == s_voxel,
          "state: a wind-only edit leaves the voxel output alone", failures)

    import geopandas as gpd
    w = gpd.read_file(wind, layer="wind")
    check(set(w["height_agl"].unique()) == {2.0, 10.0, 20.0},
          "state: the rewritten output reflects the NEW config", failures)

    # --ignore-state restores the old existence-only behaviour
    cfgfile.write_text(text)                       # back to two heights
    s_wind2 = wind.stat().st_mtime_ns
    pipeline_main(["-c", str(cfgfile), "--ignore-state"])
    check(wind.stat().st_mtime_ns == s_wind2,
          "state: --ignore-state skips on existence alone", failures)

    # and without it, the revert is picked up again
    pipeline_main(["-c", str(cfgfile)])
    w = gpd.read_file(wind, layer="wind")
    check(set(w["height_agl"].unique()) == {2.0, 10.0},
          "state: reverting the config rewrites the output back", failures)
    return failures


# ------------------------------
# 5. OVERWRITE POLICY
# ------------------------------
def verify_overwrite(tmp, pipeline_main):
    """A rerun must KEEP existing outputs, --overwrite must replace them,
    and a missing product must be refilled without touching the rest."""
    out = tmp / "out"
    cfgfile = str(tmp / "config.yaml")
    failures = []

    voxel = out / "test_voxel.nc"
    wind = out / "test_wind.gpkg"
    stamp_voxel = voxel.stat().st_mtime_ns
    stamp_wind = wind.stat().st_mtime_ns

    # plain rerun: everything kept
    pipeline_main(["-c", cfgfile])
    check(voxel.stat().st_mtime_ns == stamp_voxel,
          "overwrite: rerun keeps an existing output", failures)

    # one product deleted: refilled, the others still kept
    voxel.unlink()
    pipeline_main(["-c", cfgfile])
    check(voxel.exists(), "overwrite: a missing output is refilled",
          failures)
    check(wind.stat().st_mtime_ns == stamp_wind,
          "overwrite: refilling one product leaves the others untouched",
          failures)

    # --overwrite: everything regenerated
    stamp_wind2 = wind.stat().st_mtime_ns
    pipeline_main(["-c", cfgfile, "--overwrite"])
    check(wind.stat().st_mtime_ns != stamp_wind2,
          "overwrite: --overwrite replaces an existing output", failures)

    # GeoPackage layers must be replaced, not appended to
    import geopandas as gpd
    n_rows = len(gpd.read_file(wind, layer="wind"))
    pipeline_main(["-c", cfgfile, "--overwrite"])
    check(len(gpd.read_file(wind, layer="wind")) == n_rows,
          "overwrite: GeoPackage is replaced, not appended to", failures)

    return failures


# ------------------------------
# 6. FILL-ONLY TIME VARIABLE
# ------------------------------
def verify_regressions(tmp, out, pipeline_main):
    """Checks for bugs that shipped once and must not ship again.

    Every one of these passed silently before the fix: the values looked
    like values, the layers looked like layers.
    """
    import numpy as _np
    import pyogrio
    import struct
    from netCDF4 import Dataset
    import sys as _sys
    _sys.path.insert(0, str(Path(__file__).resolve().parent))
    from palm2gis.steps.thermo import (TemperatureUnitError,
                                       temperature_scale)

    failures = []

    # -- temperature scale -------------------------------------------------
    check(temperature_scale("degree_", "ta") == "C",
          "units: PALM's truncated 'degree_' reads as Celsius", failures)
    check(temperature_scale("K", "theta") == "K",
          "units: 'K' reads as kelvin", failures)
    try:
        temperature_scale("m/s", "ta")
        check(False, "units: an unresolvable unit raises", failures)
    except TemperatureUnitError:
        check(True, "units: an unresolvable unit raises", failures)

    # PALM's OWN `ta`, in the units PALM actually writes: Celsius,
    # truncated to "degree_" by its fixed-length character buffer.
    # Treating it as kelvin produced a voxel layer reading -253 degC.
    import shutil as _sh
    _sh.copy(tmp / "test_3d", tmp / "test_3d_ta")
    with Dataset(tmp / "test_3d_ta", "a") as nc:
        v = nc.createVariable("ta", "f4", ("time", "zu_3d", "y", "x"),
                              fill_value=-9999.0)
        v.units = "degree_"
        v[:] = 21.5
    cfg_ta = tmp / "config_ta.yaml"
    cfg_ta.write_text(
        (tmp / "config.yaml").read_text()
        .replace("file_from: test_3d", "file_from: test_3d_ta")
        .replace("file_to: test_voxel.nc", "file_to: test_voxel_ta.nc")
        # The base config never mentions `celsius`, so it takes the
        # default (false). Add it rather than substitute it.
        .replace("voxel:\n  enabled: true",
                 "celsius: true\nvoxel:\n  enabled: true"))
    pipeline_main(["-c", str(cfg_ta), "--stages", "voxel", "--overwrite"])
    with Dataset(out / "test_voxel_ta.nc") as nc:
        nc.set_auto_mask(False)
        ta = nc["ta"][:]
        ta = ta[ta != -9999.0]
        check(19.0 < float(ta.min()) < 25.0,
              f"units: PALM's own Celsius ta is not converted again "
              f"({float(ta.min()):.1f} degC)", failures)
        check(nc["ta"].palm2gis_source_units == "degree_",
              "units: the source unit is recorded in the output", failures)

    with Dataset(out / "test_voxel.nc") as nc:
        nc.set_auto_mask(False)
        th = nc["theta"][:]
        th = th[th != -9999.0]
        check(float(th.min()) > 250.0 and nc["theta"].units == "K",
              "units: theta stays in kelvin regardless of celsius",
              failures)

        # -- absolute height ---------------------------------------------
        z = _np.asarray(nc["z"][:])
        check(float(z.min()) > ORIGIN_Z,
              f"height: voxel z is above sea level "
              f"({float(z.min()):.1f} m, origin_z {ORIGIN_Z})", failures)
        check("coordinates" not in nc["ta"].ncattrs(),
              "arcgis: no auxiliary `coordinates` attribute", failures)
        check("esri_pe_string" in nc["crs"].ncattrs(),
              "arcgis: crs carries esri_pe_string", failures)

    with Dataset(out / "test_lad_3d.nc") as nc:
        # >=, not >: the lowest model level is z = 0, i.e. exactly
        # origin_z once shifted.
        check(float(_np.asarray(nc["z"][:]).min()) >= ORIGIN_Z,
              "height: LAD z is above sea level", failures)

    # -- the raw-input guard must not reject PALM's own Celsius ta ---------
    # The run above proves it end to end: a file whose `ta` is in Celsius
    # went through the voxel stage without being rejected as coord output.
    check((out / "test_voxel_ta.nc").exists(),
          "guard: a raw file holding Celsius `ta` is not rejected",
          failures)

    # -- multipatch --------------------------------------------------------
    shp = out / "test_surfaces.shp"
    st = struct.unpack("<i", shp.read_bytes()[32:36])[0]
    check(st == 31, f"multipatch: shapefile shape type is 31 (got {st})",
          failures)
    mp = pyogrio.read_dataframe(shp)
    check(mp.geometry.iloc[0].has_z,
          "multipatch: geometry carries z", failures)
    check(all(len(c) <= 10 for c in mp.columns if c != "geometry"),
          "multipatch: every field name fits DBF's 10 characters",
          failures)

    # -- file geodatabase --------------------------------------------------
    gdb = out / "test_surfaces.gdb"
    layers = {l[0] for l in pyogrio.list_layers(gdb)}
    check({"surfaces", "walls"} <= layers,
          "fgdb: both feature classes written", failures)
    fg = pyogrio.read_dataframe(gdb, layer="surfaces")
    check(str(fg["time"].dtype).startswith("datetime64"),
          "fgdb: time survives as a real datetime", failures)
    check(fg["azimuth"].isna().any(),
          "fgdb: horizontal faces carry NULL azimuth, not -9999", failures)
    check(str(fg["elem_id"].dtype) == "int32",
          "fgdb: index fields stay integers", failures)

    # -- removed config keys -----------------------------------------------
    from palm2gis.config import Config as _Config
    cfg_old = tmp / "config_removed_keys.yaml"
    cfg_old.write_text((tmp / "config.yaml").read_text()
                       .replace("  write_all_elements_gpkg: true",
                                "  write_all_elements_gpkg: true\n"
                                "  write_topview_raster: true\n"
                                "  write_walls_gpkg: false")
                       .replace("  write_volume: true",
                                "  write_volume: true\n  write_lai: true"))
    c_old = _Config(str(cfg_old))
    check("write_topview_raster" not in c_old["surfaces"]
          and "write_walls_gpkg" not in c_old["surfaces"]
          and "write_lai" not in c_old["canopy"],
          "config: removed keys load and are dropped", failures)

    # -- gridded wind: one netCDF per ABSOLUTE height ----------------------
    # 180 m a.s.l. is a PALM level (origin_z 175 + zu 5), 183 m lies
    # halfway between the levels at 182 and 184. The voxel file holds the
    # same destaggered wind on those levels, so it is the reference.
    with Dataset(out / "test_voxel.nc") as vx:
        vz = _np.asarray(vx["z"][:])
        vu = _np.ma.filled(vx["u"][0], _np.nan)
    k180, k182, k184 = (int(_np.argmin(abs(vz - h))) for h in (180, 182, 184))
    check(sorted(p.name for p in out.glob("test_wind_z*.nc"))
          == ["test_wind_z180m.nc", "test_wind_z183m.nc"],
          "wind grid: one file per absolute height, named by the height",
          failures)
    with Dataset(out / "test_wind_z180m.nc") as nc:
        check(nc["u"].dimensions == ("time", "z", "y", "x")
              and float(nc["z"][0]) == 180.0
              and nc["z"].standard_name == "height_above_mean_sea_level"
              and nc["time"].size > 1,
              "wind grid: (time, z, y, x), z = 180 m a.s.l., all records",
              failures)
        check({"u", "v", "w", "speed", "dir_met"} <= set(nc.variables),
              "wind grid: both renderer pairings present", failures)
        yv = _np.asarray(nc["y"][:])
        check(bool(_np.all(_np.diff(yv) > 0)),
              "wind grid: y ascending (not flipped like a GeoTIFF)",
              failures)
        g = _np.ma.filled(nc["u"][0, 0], _np.nan)
        ref = vu[k180]
        check(_np.array_equal(_np.isnan(g), _np.isnan(ref))
              and _np.allclose(g[~_np.isnan(g)], ref[~_np.isnan(ref)],
                               atol=1e-5),
              "wind grid: a height on a PALM level equals the voxel there "
              "(values and NoData)", failures)
    with Dataset(out / "test_wind_z183m.nc") as nc:
        g = _np.ma.filled(nc["u"][0, 0], _np.nan)
        ref = 0.5 * (vu[k182] + vu[k184])      # NaN where either is NoData
        check(float(nc["z"][0]) == 183.0
              and _np.array_equal(_np.isnan(g), _np.isnan(ref))
              and _np.allclose(g[~_np.isnan(g)], ref[~_np.isnan(ref)],
                               atol=1e-5),
              "wind grid: a height between levels is interpolated, NoData "
              "where either level is in terrain or a building", failures)
        check(int(_np.isnan(g).sum()) > 0 and not _np.isnan(g).all(),
              "wind grid: buildings are NoData, open air holds data",
              failures)

    # -- wind grid under a SUBDOMAIN CUT -----------------------------------
    # The main fixture has no x/y cut, so nx0 = 0 and nxmax = NX and every
    # coordinate-slicing mistake is a no-op. A cut must give the same
    # values at the same coordinates as the full grid.
    cfg_cut = tmp / "config_cut.yaml"
    cfg_cut.write_text(
        (tmp / "config.yaml").read_text()
        .replace("bounds:\n  zmax: 10",
                 "bounds:\n  zmax: 10\n  x0: 6\n  xmax: 30\n"
                 "  y0: 4\n  ymax: 32")
        .replace("file_to: test_wind.gpkg",
                 "file_to: test_wind_cut.gpkg")
        .replace("grid_to: test_wind.nc", "grid_to: test_wind_cut.nc"))
    rc = pipeline_main(["-c", str(cfg_cut), "--stages", "wind",
                        "--overwrite"])
    check(rc == 0, "wind grid: survives a subdomain cut", failures)
    with Dataset(out / "test_wind_z180m.nc") as full, \
            Dataset(out / "test_wind_cut_z180m.nc") as cut:
        fx, fy = _np.asarray(full["x"][:]), _np.asarray(full["y"][:])
        cx, cy = _np.asarray(cut["x"][:]), _np.asarray(cut["y"][:])
        i0 = int(_np.argmin(abs(fx - cx[0])))
        j0 = int(_np.argmin(abs(fy - cy[0])))
        fu = _np.ma.filled(full["u"][0, 0], _np.nan)[
            j0:j0 + cy.size, i0:i0 + cx.size]
        cu = _np.ma.filled(cut["u"][0, 0], _np.nan)
        ok = (_np.allclose(fx[i0:i0 + cx.size], cx)
              and _np.allclose(fy[j0:j0 + cy.size], cy)
              and _np.allclose(fu, cu, equal_nan=True))
    check(ok and cx.size < fx.size,
          "wind grid: geolocation and values match the full grid under a "
          "cut", failures)
    return failures


def verify_fill_times(tmp, pipeline_main):
    """Some PALM surface files leave the time variable at its _FillValue
    (~1e36). Those values must never be treated as seconds - they overflow
    datetime - and the stage must still produce usable layers."""
    import geopandas as gpd
    from netCDF4 import Dataset

    failures = []
    gpkg = tmp / "out" / "test_surfaces.gpkg"

    def second_time():
        el = gpd.read_file(gpkg, layer="surfaces")
        t = el[el.rec == sorted(el.rec.unique())[1]]["time"].iloc[0]
        return str(t)[:19]

    with Dataset(tmp / "test_surf", "a") as nc:
        nc["time"][:] = 9.969209968386869e36

    # (a) no reconstruction configured: record-index labels, empty times
    rc = pipeline_main(["-c", str(tmp / "config.yaml"), "--stages",
                        "surfaces", "--overwrite"])
    check(rc == 0, "fill times: the stage completes instead of overflowing",
          failures)
    el = gpd.read_file(gpkg, layer="surfaces")
    check(el["time"].isna().all() and "rec" in el.columns,
          "fill times: elements carry an empty time and a rec column",
          failures)

    # (b) the p3d supplies the interval - no manual settings at all
    (tmp / "p3d_probe").write_text(
        "&initialization_parameters\n nz = 12,\n dz = 2.0,\n/\n"
        "&runtime_parameters\n dt_dosurf_av = 600.0,\n"
        " skip_time_dosurf_av = 0.0,\n/\n")
    cfgp0 = tmp / "config.yaml"
    base = cfgp0.read_text()
    cfgp0.write_text(base.replace("  nz: 12\n  dz: 2.0", "  p3d: p3d_probe"))
    pipeline_main(["-c", str(cfgp0), "--stages", "surfaces", "--overwrite"])
    check(second_time() == "2026-07-15 10:20:00",
          "fill times: interval read from the p3d automatically", failures)
    cfgp0.write_text(base)

    # (c) bounds.time_step given: real timestamps reconstructed
    cfgp = tmp / "config.yaml"
    original = cfgp.read_text()
    cfgp.write_text(original.replace("bounds:\n  zmax: 10",
                                     "bounds:\n  zmax: 10\n"
                                     "  time_step: 600\n  time_start: 0.0"))
    pipeline_main(["-c", str(cfgp), "--stages", "surfaces", "--overwrite"])
    check(second_time() == "2026-07-15 10:10:00",
          "fill times: bounds.time_step reconstructs real timestamps",
          failures)
    cfgp.write_text(original)
    return failures


def verify_time_regularity():
    """PALM record times carry a per-record residue (output at the first
    timestep at or after the target). Rounded, they gave 601/599 s steps
    that ArcGIS Pro reads as an irregular voxel time axis."""
    import numpy as np
    from palm2gis.steps.georef import regularize_times

    failures = []
    raw = np.array([600.4, 1200.9, 1800.2, 2400.6, 3000.51])
    out = regularize_times(raw, range(raw.size), tag="test")
    check(np.allclose(np.diff(out), 600.0),
          "time: residue snapped to a regular 600 s interval", failures)
    gap = np.array([600.4, 1200.9, 2400.6])
    out = regularize_times(gap, range(gap.size), tag="test")
    check(np.allclose(out, [600.0, 1200.0, 2400.0]),
          "time: a missing record keeps the true 600 s interval", failures)
    odd = np.array([0.0, 700.0, 1200.0, 2000.0])
    out = regularize_times(odd, range(odd.size), tag="test")
    check(np.allclose(out, odd),
          "time: genuinely irregular times are left as written", failures)
    return failures


def verify_descriptions(out):
    """ArcGIS Pro shows a NetCDF variable's long_name as its Description;
    GeoPackage and File Geodatabase fields carry an alias."""
    import sqlite3
    from netCDF4 import Dataset
    from palm2gis.descriptions import DESCRIPTIONS

    failures = []
    with Dataset(out / "test_voxel.nc") as nc:
        bare = [v for v in nc.variables
                if nc[v].dimensions[:1] == ("time",) and len(nc[v].dimensions)
                == 4 and str(getattr(nc[v], "long_name", v)) == v]
        check(not bare, "descriptions: every voxel variable has a real "
              "long_name{}".format(" (missing: {})".format(bare)
                                   if bare else ""), failures)
        check(nc["theta"].long_name == DESCRIPTIONS["theta"][1],
              "descriptions: theta described as potential temperature",
              failures)

    con = sqlite3.connect(str(out / "test_surfaces.gpkg"))
    rows = dict(con.execute("SELECT column_name, description FROM "
                            "gpkg_data_columns WHERE table_name = "
                            "'surfaces'").fetchall())
    con.close()
    check(rows.get("t_surface") == DESCRIPTIONS["t_surface"][1]
          and "orient" in rows,
          "descriptions: surfaces GeoPackage fields carry descriptions",
          failures)

    con = sqlite3.connect(str(out / "test_wind.gpkg"))
    n = con.execute("SELECT count(*) FROM gpkg_data_columns WHERE "
                    "table_name = 'wind' AND column_name = 'dir_met'"
                    ).fetchone()[0]
    con.close()
    check(n == 1, "descriptions: wind GeoPackage fields carry descriptions",
          failures)

    try:
        import pyarrow  # noqa: F401
    except ImportError:
        print("  [skip] descriptions: pyarrow not installed, no FGDB aliases")
        return failures
    import pyogrio
    _meta, tbl = pyogrio.read_arrow(out / "test_surfaces.gdb",
                                    layer="surfaces")
    md = tbl.schema.field("t_surface").metadata or {}
    check(md.get(b"GDAL:OGR:alternative_name")
          == DESCRIPTIONS["t_surface"][0].encode(),
          "descriptions: File Geodatabase fields carry aliases", failures)
    gdb = pyogrio.read_dataframe(out / "test_surfaces.gdb", layer="surfaces")
    check(gdb["azimuth"].isna().any() and len(gdb) == 432 * NT,
          "descriptions: FGDB still holds every row, NULL azimuths intact",
          failures)
    return failures


def verify_basic_state(tmp):
    """PALM basic state (surface_pressure, pt_surface) per domain.

    Mirrors what real output showed (holesovice_stern_adapt): the root
    domain runs with the dynamic driver's first surface pressure, a nested
    child with the PARENT p3d's surface_pressure, both with pt_surface; and
    ta / theta, when present, is PALM's exact Exner function.
    """
    import shutil as _sh
    import numpy as _np
    from netCDF4 import Dataset
    from palm2gis.config import Config
    from palm2gis.steps.georef import process_static_driver
    from palm2gis.steps.thermo import exner_basic_state

    failures = []
    base = (tmp / "config.yaml").read_text()
    grid_lines = "  nz: {}\n  dz: {}".format(NZ, DZ)

    (tmp / "bs_p3d").write_text(
        " &initialization_parameters\n"
        "    nz = {}, dz = {},\n"
        "    surface_pressure = 1000.,\n"
        "    pt_surface = 285.,\n /\n".format(NZ, DZ))
    (tmp / "bs_p3d_N02").write_text(
        " &initialization_parameters\n"
        "    nz = {}, dz = {},\n /\n".format(NZ, DZ))

    def run(name, p3d, extra="", data_3d="test_3d"):
        cfgp = tmp / name
        cfgp.write_text(
            base.replace(grid_lines, "  p3d: " + p3d)
            .replace("file_from: test_3d", "file_from: " + data_3d) + extra)
        cfg = Config(cfgp)
        process_static_driver(cfg)
        return cfg["surface_pressure"], cfg["pt_surface"]

    # child without its own surface_pressure -> parent p3d, not 1013.25
    p_s, th = run("config_bs_child.yaml", "bs_p3d_N02")
    check(abs(p_s - 1.0e5) < 1e-6 and abs(th - 285.0) < 1e-9,
          f"basic state: child inherits the parent p3d ({p_s / 100:.2f} "
          f"hPa, {th:.1f} K)", failures)

    # root without a dynamic driver -> its own p3d
    p_s, th = run("config_bs_root.yaml", "bs_p3d")
    check(abs(p_s - 1.0e5) < 1e-6 and abs(th - 285.0) < 1e-9,
          "basic state: root domain reads its own p3d", failures)

    # root with <job>_dynamic next to the p3d -> first driver value wins
    with Dataset(tmp / "bs_dynamic", "w") as nc:
        nc.createDimension("time", 3)
        v = nc.createVariable("surface_forcing_surface_pressure", "f4",
                              ("time",), fill_value=-9999.0)
        v[:] = [99388.42, 99394.99, 98566.35]
    p_s, th = run("config_bs_dyn.yaml", "bs_p3d")
    check(abs(p_s - 99388.42) < 0.01 and abs(th - 285.0) < 1e-9,
          f"basic state: root domain takes the dynamic driver's first "
          f"surface pressure ({p_s / 100:.2f} hPa)", failures)
    # ... but a child never does
    p_s, _ = run("config_bs_child2.yaml", "bs_p3d_N02")
    check(abs(p_s - 1.0e5) < 1e-6,
          "basic state: a child ignores the root's dynamic driver", failures)
    (tmp / "bs_dynamic").unlink()

    # ta / theta consistent with 993.88 hPa while the p3d says 1000:
    # the fit replaces the namelist value
    _sh.copy(tmp / "test_3d", tmp / "test_3d_bs")
    with Dataset(tmp / "test_3d_bs", "a") as nc:
        z = _np.asarray(nc["zu_3d"][:], dtype="f8")
        th3 = nc["theta"][:]
        pi = exner_basic_state(z, 99388.42, 285.0)[None, :, None, None]
        v = nc.createVariable("ta", "f4", ("time", "zu_3d", "y", "x"),
                              fill_value=-9999.0)
        v.units = "degree_"
        v[:] = th3 * pi - 273.15
    p_s, th = run("config_bs_fit.yaml", "bs_p3d", data_3d="test_3d_bs")
    check(abs(p_s - 99388.42) < 5.0 and abs(th - 285.0) < 0.5,
          f"basic state: ta / theta overrides a wrong namelist value "
          f"({p_s / 100:.2f} hPa, {th:.2f} K)", failures)

    # an explicit config value is kept even when ta / theta disagrees
    p_s, _ = run("config_bs_explicit.yaml", "bs_p3d", data_3d="test_3d_bs",
                 extra="surface_pressure: 100000.0\n")
    check(abs(p_s - 1.0e5) < 1e-6,
          "basic state: an explicit config value wins over the fit",
          failures)
    return failures


def verify_t_surface_rename(tmp):
    """`t_surf` became `t_surface` in 0.22; old configs keep working."""
    from palm2gis.config import Config
    failures = []
    legacy = tmp / "config_legacy_tsurf.yaml"
    legacy.write_text((tmp / "config.yaml").read_text()
                      .replace("derive_t_surface: true",
                               "derive_t_surf: true"))
    cfg = Config(legacy)
    check(cfg["surfaces"]["derive_t_surface"] is True
          and "derive_t_surf" not in cfg["surfaces"],
          "rename: legacy surfaces.derive_t_surf -> derive_t_surface",
          failures)
    old_new = tmp / "config_new_tsurf.yaml"
    old_new.write_text((tmp / "config_new.yaml").read_text()
                       .replace("[theta_surface, t_surface]",
                                "[theta_surface, t_surf]"))
    cfg = Config(old_new)
    check(cfg["surfaces"]["derive_t_surface"] is True
          and cfg["surfaces"]["variables"] == ["theta_surface"],
          "rename: t_surf in surfaces.variables is read as t_surface",
          failures)
    return failures


def verify_layout(tmp, pipeline_main):
    """The newcomer layout must configure the same run as the legacy one,
    reject typos and legacy keys, and produce only the files it names."""
    import geopandas as gpd
    from netCDF4 import Dataset
    from palm2gis.config import Config, ConfigError

    failures = []
    new = tmp / "config_new.yaml"
    new.write_text(f"""\
project:
    name: test
    output_dir: ./out_new
input:
    static_driver: test_static
    data_3d: test_3d
    surface_data: test_surf
    crs: "{CRS}"
region:
    height_above_surface: 10
voxel:
    variables: [u, v, w, theta, ta]
    netcdf:
        file: '{{name}}_voxel.nc'
surfaces:
    variables: [theta_surface, t_surface]
    all_elements_gpkg:
        file: '{{name}}_surfaces.gpkg'
    all_elements_fgdb:
        file: '{{name}}_surfaces.gdb'
    all_elements_multipatch:
        file: '{{name}}_surfaces.shp'
wind:
    heights_above_surface: [2.0, 10.0]
    netcdf:
        file: '{{name}}_wind.nc'
        heights: [180, 183]
    gpkg:
        file: '{{name}}_wind.gpkg'
        every_nth_column: 2
canopy:
    lad_netcdf:
        file: '{{name}}_lad_3d.nc'
advanced:
    nz: {NZ}
    dz: {DZ}
""")
    old, cnew = Config(str(tmp / "config.yaml")), Config(str(new))
    check(cnew.layout == "new" and old.layout == "legacy",
          "layout: both layouts recognised", failures)
    same = all(
        {k: old[st][k] for k in keys} == {k: cnew[st][k] for k in keys}
        for st, keys in (
            ("domain", ("static_driver", "nz", "dz")),
            ("bounds", ("zmax", "time_from", "time_to")),
            ("voxel", ("file_from", "file_to", "do_velocity", "scalars",
                       "derive_ta")),
            ("surfaces", ("file_from", "variables", "derive_t_surface",
                          "write_all_elements_gpkg", "write_fgdb",
                          "write_multipatch")),
            ("wind", ("file_from", "file_to", "grid_to", "grid_heights",
                      "heights_agl", "stride")),
            ("canopy", ("write_volume",))))
    check(same, "layout: new config translates to the legacy settings",
          failures)

    rc = pipeline_main(["-c", str(new), "--overwrite"])
    out = tmp / "out_new"
    names = {p.name for p in out.iterdir()
             if p.suffix in (".nc", ".gpkg", ".gdb", ".shp")}
    want = {"test_voxel.nc", "test_surfaces.gpkg", "test_surfaces.gdb",
            "test_surfaces.shp", "test_wind_z180m.nc", "test_wind_z183m.nc",
            "test_wind.gpkg", "test_lad_3d.nc"}
    check(rc == 0 and names == want,
          f"layout: run writes exactly the named files ({sorted(names)})",
          failures)
    # voxel `speed` (horizontal) requested through the new layout
    sp_cfg = tmp / "config_speed_h.yaml"
    sp_cfg.write_text(new.read_text().replace("./out_new", "./out_speed_h")
                      .replace("[u, v, w, theta, ta]",
                               "[u, v, w, speed]"))
    rc_sp = pipeline_main(["-c", str(sp_cfg), "--stages", "voxel"])
    with Dataset(tmp / "out_speed_h" / "test_voxel.nc") as nc:
        import numpy as _npx
        su = _npx.ma.filled(nc["u"][:], _npx.nan)
        sv = _npx.ma.filled(nc["v"][:], _npx.nan)
        ss = _npx.ma.filled(nc["speed"][:], _npx.nan)
        ok_sp = (rc_sp == 0 and nc["speed"].standard_name == "wind_speed"
                 and _npx.array_equal(_npx.isnan(ss), _npx.isnan(su))
                 and _npx.allclose(ss[~_npx.isnan(ss)],
                                   _npx.hypot(su, sv)[~_npx.isnan(ss)],
                                   atol=1e-5))
    check(ok_sp, "voxel: `speed` in voxel.variables is sqrt(u^2 + v^2)",
          failures)

    a = gpd.read_file(tmp / "out" / "test_surfaces.gpkg", layer="surfaces")
    b = gpd.read_file(out / "test_surfaces.gpkg", layer="surfaces")
    check(len(a) == len(b) and list(a.columns) == list(b.columns),
          "layout: same surface layer as the legacy config", failures)

    # opt-in: only the wind NetCDF, no GeoPackage - legacy validation used
    # to demand file_to even though the stage handles this case
    only = tmp / "config_new_only.yaml"
    only.write_text(new.read_text()
                    .replace("./out_new", "./out_only")
                    .replace("file: '{name}_wind.gpkg'", "file: null")
                    .replace("file: '{name}_surfaces.gdb'", "file: null")
                    .replace("file: '{name}_surfaces.shp'", "file: null"))
    rc = pipeline_main(["-c", str(only), "--stages", "wind", "surfaces"])
    names = {p.name for p in (tmp / "out_only").iterdir()
             if p.suffix in (".nc", ".gpkg", ".gdb", ".shp")}
    check(rc == 0 and names == {"test_wind_z180m.nc", "test_wind_z183m.nc",
                                "test_surfaces.gpkg"},
          f"layout: null files are skipped, NetCDF-only wind runs "
          f"({sorted(names)})", failures)
    with Dataset(tmp / "out_only" / "test_wind_z183m.nc") as nc:
        check(nc["u"].dimensions == ("time", "z", "y", "x"),
              "layout: NetCDF-only wind grid is complete", failures)

    # -- absolute humidity ---------------------------------------------------
    from palm2gis.steps.thermo import KAPPA, P0, RD, humidity_scale
    ah_cfg = tmp / "config_ah.yaml"
    ah_cfg.write_text(new.read_text().replace("./out_new", "./out_ah")
                      .replace("[u, v, w, theta, ta]", "[theta, ta, ah]")
                      .replace("advanced:",
                               "advanced:\n    units:\n"
                               "        temperature: K", 1))
    rc = pipeline_main(["-c", str(ah_cfg), "--stages", "voxel"])
    import numpy as _np
    with Dataset(tmp / "out_ah" / "test_voxel.nc") as nc, \
            Dataset(tmp / "test_3d") as src:
        nc.set_auto_mask(False)
        ah = nc["ah"][0].astype("f8")
        ta = nc["ta"][0].astype("f8")
        th = nc["theta"][0].astype("f8")
        # drop_k0: output level k is source level k + 1
        qs = _np.ma.filled(src["q"][0, 1:, :, :], _np.nan).astype("f8")
        ok = (ah != -9999.0) & _np.isfinite(qs)
        # independent of palm2gis' pressure: ta / theta IS the Exner
        # function, so p = P0 * (ta / theta) ** (1 / kappa)
        p = P0 * (ta / th) ** (1.0 / KAPPA)
        want = 1.0e3 * qs * p / (RD * ta * (1.0 + 0.608 * qs))
        err = float(_np.max(_np.abs(ah[ok] - want[ok]) / want[ok]))
        check(rc == 0 and ok.sum() > 0 and err < 1.0e-4
              and 8.0 < float(_np.median(ah[ok])) < 11.0
              and nc["ah"].units == "g m-3"
              and nc["ah"].derived_from == "q",
              f"voxel: ah = q p / (Rd Tv), median "
              f"{float(_np.median(ah[ok])):.2f} g/m3, max rel. error "
              f"{err:.1e}", failures)
        check(bool((ah[~ok] == -9999.0).all()),
              "voxel: ah is NoData inside terrain", failures)
    check(humidity_scale("g/kg") == 1.0e-3 and humidity_scale("kg/kg") == 1.0
          and humidity_scale(None) == 1.0,
          "voxel: q units kg/kg and g/kg recognised", failures)

    # with PALM's own ta (degC, truncated "degree_") as the temperature
    via_ta = tmp / "config_ah_ta.yaml"
    via_ta.write_text(ah_cfg.read_text().replace("./out_ah", "./out_ah_ta")
                      .replace("data_3d: test_3d", "data_3d: test_3d_ta")
                      .replace("[theta, ta, ah]", "[ah]"))
    rc = pipeline_main(["-c", str(via_ta), "--stages", "voxel"])
    with Dataset(tmp / "out_ah_ta" / "test_voxel.nc") as nc:
        nc.set_auto_mask(False)
        a2 = nc["ah"][0].astype("f8")
        t_k = 21.5 + 273.15
        want2 = 1.0e3 * qs * p / (RD * t_k * (1.0 + 0.608 * qs))
        ok2 = (a2 != -9999.0) & _np.isfinite(qs)
        err2 = float(_np.max(_np.abs(a2[ok2] - want2[ok2]) / want2[ok2]))
        check(rc == 0 and ok2.sum() > 0 and err2 < 1.0e-4
              and "PALM's ta" in nc["ah"].comment,
              f"voxel: ah uses PALM's own Celsius ta when present "
              f"(max rel. error {err2:.1e})", failures)

    # PALM's own ah wins over deriving it
    import shutil as _sh
    _sh.copy(tmp / "test_3d", tmp / "test_3d_ah")
    with Dataset(tmp / "test_3d_ah", "a") as nc:
        v = nc.createVariable("ah", "f4", ("time", "zu_3d", "y", "x"),
                              fill_value=-9999.0)
        v.units = "g m-3"
        v[:] = 12.0
    own = tmp / "config_ah_own.yaml"
    own.write_text(ah_cfg.read_text().replace("./out_ah", "./out_ah_own")
                   .replace("data_3d: test_3d", "data_3d: test_3d_ah")
                   .replace("[theta, ta, ah]", "[ah]"))
    rc = pipeline_main(["-c", str(own), "--stages", "voxel"])
    with Dataset(tmp / "out_ah_own" / "test_voxel.nc") as nc:
        a = nc["ah"][0]
        check(rc == 0 and _np.allclose(a.compressed(), 12.0)
              and "derived_from" not in nc["ah"].ncattrs(),
              "voxel: PALM's own ah is copied, not derived", failures)

    def rejects(text, expect, msg):
        bad = tmp / "config_bad.yaml"
        bad.write_text(text)
        try:
            Config(str(bad))
            check(False, msg, failures)
        except ConfigError as exc:
            check(expect in str(exc), f"{msg} ({exc})", failures)

    from palm2gis.layout import to_new_names
    check(to_new_names("raise bounds.zmax, set wind.file_from") ==
          "raise region.height_above_surface, set input.data_3d"
          and to_new_names("mybounds.zmax") == "mybounds.zmax",
          "layout: messages use the new setting names", failures)

    base = new.read_text()

    # -- 0.25 key names, and the pre-0.25 ones they replaced ----------------
    def settings(text, name="config_names.yaml"):
        p = tmp / name
        p.write_text(text)
        return Config(str(p))

    c_new = settings(base)
    check(c_new["celsius"] is True and c_new["bounds"]["zmax"] == 10,
          "layout: temperature is Celsius by default, height_above_surface "
          "sets the cut", failures)
    c_old = settings(base.replace("height_above_surface",
                                  "height_above_highest_surface")
                     + "celsius: false\n", "config_oldnames.yaml")
    check(c_old["bounds"]["zmax"] == 10 and c_old["celsius"] is False,
          "layout: pre-0.25 height_above_highest_surface and celsius still "
          "work", failures)
    c_k = settings(base.replace("advanced:",
                                "advanced:\n    units:\n"
                                "        temperature: K", 1),
                   "config_kelvin.yaml")
    check(c_k["celsius"] is False,
          "layout: advanced.units.temperature: K selects kelvin", failures)
    c_t = settings(base + "time:\n    from: 100\n    to: 500\n",
                   "config_time.yaml")
    check(c_t["bounds"]["time_from"] == 100 and c_t["bounds"]["time_to"] == 500,
          "layout: the time section sets the window", failures)
    c_g = settings(base.replace("advanced:\n    nz:",
                                "advanced:\n    grid:\n        nz:")
                   .replace("\n    dz:", "\n        dz:"),
                   "config_grid.yaml")
    check(c_g["domain"]["nz"] == c_new["domain"]["nz"]
          and c_g["domain"]["dz"] == c_new["domain"]["dz"],
          "layout: advanced.grid.nz / dz are read like the flat names",
          failures)

    rejects(base.replace("advanced:",
                         "advanced:\n    units:\n"
                         "        potential_temperature: C", 1),
            "must be K", "layout: a Celsius potential temperature is "
                         "rejected")
    rejects(base.replace("advanced:",
                         "advanced:\n    units:\n"
                         "        temperature: F", 1),
            "must be K or C", "layout: an unknown temperature scale is "
                              "rejected")
    rejects(base.replace("every_nth_column: 2", "every_nth_colum: 2"),
            "Did you mean 'every_nth_column'",
            "layout: a misspelt key is rejected with a suggestion")
    rejects(base + "bounds:\n    zmax: 10\n", "old config layout",
            "layout: a legacy block in the new layout is explained")
    rejects(base.replace("[u, v, w, theta, ta]", "[u, v, theta]"),
            "u, v and w together",
            "layout: partial velocity list is rejected")
    rejects(base.replace("[theta_surface, t_surface]", "theta_surface"),
            "must be a list",
            "layout: a bare name instead of a list is rejected, not split "
            "into letters")
    rejects(base.replace("canopy:\n", "canopy:\n    enabled: 'false'\n"),
            "true or false", "layout: a quoted boolean is rejected")
    rejects(base.replace("every_nth_column: 2", "every_nth_column: 0"),
            "at least 1", "layout: every_nth_column 0 is rejected")
    rejects(base.replace("{name}_surfaces.gdb", "{name}_surfaces"),
            "must end in .gdb", "layout: a wrong output extension is "
                                "rejected before anything is written")
    rejects(base.replace("{name}_wind.nc", "{name}_voxel_z{z}.nc")
            .replace("heights: [180, 183]", "heights: [180]")
            .replace("{name}_voxel.nc", "{name}_voxel_z180.nc"),
            "would both write",
            "config: two outputs on one path are rejected")
    # no wind.netcdf.heights: one file per PALM level of the voxel range,
    # on exactly the voxel's z values
    import numpy as _np_thermo
    auto = tmp / "config_wind_auto.yaml"
    auto.write_text(base.replace("./out_new", "./out_wind_auto")
                    .replace("        heights: [180, 183]\n", ""))
    rc = pipeline_main(["-c", str(auto), "--stages", "voxel", "wind"])
    ad = tmp / "out_wind_auto"
    with Dataset(ad / "test_voxel.nc") as vx:
        vz = [float(z) for z in vx["z"][:]]
        vu0 = _np_thermo.ma.filled(vx["u"][0], _np_thermo.nan)
    files = sorted(ad.glob("test_wind_z*.nc"))
    zs = []
    for f in files:
        with Dataset(f) as nc:
            zs.append(float(nc["z"][0]))
    ok = rc == 0 and sorted(zs) == vz
    if ok:
        k = len(vz) // 2
        with Dataset(ad / "test_wind_z{:g}m.nc".format(vz[k])) as nc:
            g = _np_thermo.ma.filled(nc["u"][0, 0], _np_thermo.nan)
        ok = (_np_thermo.array_equal(_np_thermo.isnan(g),
                                     _np_thermo.isnan(vu0[k]))
              and _np_thermo.allclose(g[~_np_thermo.isnan(g)],
                                      vu0[k][~_np_thermo.isnan(g)],
                                      atol=1e-6))
    check(ok, "wind: without heights, one file per voxel level "
              f"({len(files)} files on the voxel's z, values equal)",
          failures)

    # air temperature alone is a valid voxel request
    ta_only = tmp / "config_ta_only.yaml"
    ta_only.write_text(base.replace("./out_new", "./out_ta_only")
                       .replace("[u, v, w, theta, ta]", "[ta]"))
    rc = pipeline_main(["-c", str(ta_only), "--stages", "voxel"])
    with Dataset(tmp / "out_ta_only" / "test_voxel.nc") as nc:
        check(rc == 0 and "ta" in nc.variables and "u" not in nc.variables,
              "voxel: ta alone is accepted and written", failures)

    # 3D wind speed: with the components, and on its own
    sp_cfg = tmp / "config_speed.yaml"
    sp_cfg.write_text(base.replace("./out_new", "./out_speed")
                      .replace("[u, v, w, theta, ta]", "[u, v, w, speed_3d]"))
    rc = pipeline_main(["-c", str(sp_cfg), "--stages", "voxel"])
    with Dataset(tmp / "out_speed" / "test_voxel.nc") as nc:
        import numpy as _np
        sp = nc["speed_3d"][:].astype("f8")
        uu, vv, ww = (nc[c][:].astype("f8") for c in ("u", "v", "w"))
        ok = ~_np.ma.getmaskarray(sp)
        want = _np.sqrt(uu ** 2 + vv ** 2 + ww ** 2)
        err = float(_np.max(_np.abs(sp[ok] - want[ok])))
        check(rc == 0 and ok.any() and err < 1.0e-4
              and float(_np.abs(ww[ok]).max()) > 0.0
              and nc["speed_3d"].units == "m s-1"
              and bool((_np.ma.getmaskarray(uu) == ~ok).all()),
              f"voxel: speed_3d = sqrt(u^2+v^2+w^2), max error {err:.2e}, "
              "same NoData as u", failures)
    sp_only = tmp / "config_speed_only.yaml"
    sp_only.write_text(sp_cfg.read_text().replace("./out_speed",
                                                  "./out_speed_only")
                       .replace("[u, v, w, speed_3d]", "[speed_3d]"))
    rc = pipeline_main(["-c", str(sp_only), "--stages", "voxel"])
    with Dataset(tmp / "out_speed_only" / "test_voxel.nc") as nc, \
            Dataset(tmp / "out_speed" / "test_voxel.nc") as ref:
        check(rc == 0 and "u" not in nc.variables
              and _np.ma.allclose(nc["speed_3d"][:], ref["speed_3d"][:]),
              "voxel: speed_3d alone is written without u, v, w", failures)
    cf_cfg = tmp / "config_cf.yaml"
    cf_cfg.write_text(sp_only.read_text()
                      .replace("./out_speed_only", "./out_cf")
                      .replace("advanced:", "advanced:\n    netcdf_crs: cf", 1))
    rc = pipeline_main(["-c", str(cf_cfg), "--stages", "voxel"])
    with Dataset(tmp / "out_cf" / "test_voxel.nc") as nc:
        check(rc == 0
              and nc["crs"].grid_mapping_name == "transverse_mercator"
              and "crs_wkt" in nc["crs"].ncattrs(),
              "voxel: advanced.netcdf_crs: cf adds the CF parameters",
              failures)
    old_sp = Config(str(sp_only))
    check(old_sp["voxel"]["derive_speed_3d"] is True
          and "speed_3d" not in old_sp["voxel"]["scalars"]
          and old_sp["voxel"]["do_velocity"] is False,
          "layout: speed_3d maps to derive_speed_3d", failures)

    # a sample above the file's top level is NoData, not the top value
    high = tmp / "config_high.yaml"
    high.write_text(base.replace("./out_new", "./out_high")
                    .replace("[2.0, 10.0]", "[2.0, 50.0]")
                    .replace("heights: [180, 183]", "heights: [180, 999]"))
    rc = pipeline_main(["-c", str(high), "--stages", "wind"])
    check(rc == 0 and (tmp / "out_high" / "test_wind_z180m.nc").exists()
          and not (tmp / "out_high" / "test_wind_z999m.nc").exists(),
          "wind: an absolute height above the file's levels is skipped",
          failures)
    pts = gpd.read_file(tmp / "out_high" / "test_wind.gpkg", layer="wind")
    check(pts[pts.height_agl == 50.0]["u"].isna().all(),
          "wind: the point layer carries NULL above the top level",
          failures)
    return failures


# ------------------------------
# 7. MAIN
# ------------------------------
# ------------------------------
# 9. OUTPUT CRS (project.output_crs: EPSG:4326)
# ------------------------------
def verify_output_crs(tmp, pipeline_main):
    """Every product in EPSG:4326 for an ArcGIS scene in WGS 84.

    Vector layers must be transformed exactly with z untouched; gridded
    netCDF resampled onto regular lon/lat axes with no auxiliary
    coordinates (Esri's voxel requirement); u/v turned from UTM grid axes
    to true east/north. The synthetic domain sits at ~14.4 E, west of
    zone 33's central meridian - where Prague is - so the rotation is a
    real +0.45 degrees, checked against pyproj's meridian convergence.
    """
    import geopandas as gpd
    import numpy as np
    from netCDF4 import Dataset
    from pyproj import Geod, Proj, Transformer

    failures = []
    cfg = tmp / "config_4326.yaml"
    cfg.write_text(f"""\
project:
    name: test
    output_dir: ./out_4326
    output_crs: 'EPSG:4326'
input:
    static_driver: test_static
    data_3d: test_3d
    surface_data: test_surf
    crs: "{CRS}"
region:
    height_above_surface: 10
voxel:
    variables: [u, v, w, theta]
    netcdf:
        file: '{{name}}_voxel.nc'
surfaces:
    variables: [theta_surface]
    all_elements_gpkg:
        file: '{{name}}_surfaces.gpkg'
    all_elements_multipatch:
        file: '{{name}}_surfaces.shp'
wind:
    heights_above_surface: [2.0, 10.0]
    netcdf:
        file: '{{name}}_wind.nc'
        heights: [180, 183]
        every_nth_column: 2
    gpkg:
        file: '{{name}}_wind.gpkg'
        every_nth_column: 2
        arrows_3d: true
canopy:
    lad_netcdf:
        file: '{{name}}_lad_3d.nc'
advanced:
    grid:
        nz: {NZ}
        dz: {DZ}
""")
    rc = pipeline_main(["-c", str(cfg), "--overwrite"])
    check(rc == 0, "output crs: EPSG:4326 run finishes", failures)
    if rc != 0:
        return failures
    out = tmp / "out_4326"
    to_ll = Transformer.from_crs(CRS, "EPSG:4326", always_xy=True)

    # -- voxel: regular lon/lat, WGS 84, no auxiliary coordinates --------
    with Dataset(out / "test_voxel.nc") as nc:
        lon, lat = nc["lon"][:], nc["lat"][:]
        check(nc["theta"].dimensions == ("time", "z", "lat", "lon")
              and nc["lon"].units == "degrees_east"
              and nc["lat"].standard_name == "latitude",
              "output crs: voxel on (time, z, lat, lon) with CF lon/lat",
              failures)
        check(nc["crs"].EPSG_code == "EPSG:4326"
              and 'AUTHORITY["EPSG","4326"]' in nc["crs"].crs_wkt,
              "output crs: voxel grid mapping is EPSG:4326", failures)
        check(not any("coordinates" in nc[v].ncattrs()
                      for v in nc.variables),
              "output crs: no auxiliary coordinates (voxel-layer source)",
              failures)
        check(np.ptp(np.diff(lon)) < 1e-12 and np.ptp(np.diff(lat)) < 1e-12
              and np.all(np.diff(lon) > 0) and np.all(np.diff(lat) > 0),
              "output crs: lon/lat axes regular and ascending", failures)
        # Native 2 m cells at ~50 N: ~2.8e-5 deg of lon, ~1.8e-5 of lat.
        check(abs(float(lon[1] - lon[0]) - 2.0 / 71_600) < 2e-6
              and abs(float(lat[1] - lat[0]) - 2.0 / 111_200) < 1e-6,
              "output crs: spacing matches the native 2 m cell", failures)
        th = nc["theta"][0].filled(np.nan)
        # Building block: local cells 8..11, top at terrain + 8 m. Its
        # centre (local x = y = 20 m) must be NoData at the ground and
        # air above the roof - the mask landed in the right place.
        blon, blat = to_ll.transform(ORIGIN_X + 20.0, ORIGIN_Y + 20.0)
        i = int(np.argmin(abs(lon - blon)))
        j = int(np.argmin(abs(lat - blat)))
        col = th[:, j, i]
        check(np.isnan(col[0]) and np.isfinite(col[-1]),
              "output crs: building mask lands on the building", failures)
        # Air cells only are averaged: every finite theta stays within the
        # model's own range (300 +- noise), none is a mix with NoData.
        fin = th[np.isfinite(th)]
        check(fin.size and 299.0 < fin.min() and fin.max() < 301.0,
              "output crs: no fill value blended into the average",
              failures)
        check(np.isnan(th[:, 0, 0]).all() or np.isnan(th[:, -1, -1]).all()
              or np.isnan(th[:, 0, -1]).all() or np.isnan(th[:, -1, 0]).all(),
              "output crs: a corner outside the rotated domain is NoData",
              failures)
        u = nc["u"][0].filled(np.nan)
        v = nc["v"][0].filled(np.nan)

    # -- rotation of u, v ------------------------------------------------
    # The angle palm2gis turns u/v by must be minus PALM's meridian
    # convergence at every column (pyproj computes that independently).
    from palm2gis.config import Config
    from palm2gis.steps.georef import process_static_driver
    pg = process_static_driver(Config(str(cfg)))
    cos_a, sin_a = pg.reproj.rotation()
    alpha = np.degrees(np.arctan2(sin_a, cos_a))
    xx, yy = np.meshgrid(pg.abs_x(pg.x_stat), pg.abs_y(pg.y_stat))
    lo, la = to_ll.transform(xx, yy)
    gamma = Proj(CRS).get_factors(lo, la).meridian_convergence
    check(np.max(abs(alpha + gamma)) < 1e-6 and abs(alpha.mean()) > 0.4,
          f"output crs: u/v turned by the meridian convergence "
          f"({alpha.mean():+.3f} deg, max error "
          f"{np.max(abs(alpha + gamma)):.1e})", failures)
    ug, vg = np.array([2.0]), np.array([0.0])
    ue, vn = pg.reproj.rotate_at(ug, vg, np.array([0]), np.array([0]))
    check(abs(np.hypot(ue, vn)[0] - 2.0) < 1e-12 and vn[0] > 0,
          "output crs: rotation keeps speed; grid-east wind gains a "
          "northward part west of the central meridian", failures)
    m = np.isfinite(u) & np.isfinite(v)
    check(m.any(), "output crs: voxel u/v hold data", failures)

    # -- wind: arrows exact, attributes in the same frame ----------------
    w3 = gpd.read_file(out / "test_wind.gpkg", layer="wind_3d")
    pts = gpd.read_file(out / "test_wind.gpkg", layer="wind")
    check(w3.crs.to_epsg() == 4326 and pts.crs.to_epsg() == 4326
          and w3.geometry.has_z.all() and not pts.geometry.has_z.any(),
          "output crs: wind points 2D and arrows 3D, both EPSG:4326",
          failures)
    g = w3[w3.speed > 0]
    c = np.array([np.array(ln.coords) for ln in g.geometry])
    az, _, _ = Geod(ellps="WGS84").inv(c[:, 0, 0], c[:, 0, 1],
                                       c[:, 1, 0], c[:, 1, 1])
    d = (90.0 - az) - np.degrees(np.arctan2(g.v, g.u))
    d = (d + 180.0) % 360.0 - 180.0
    check(np.max(abs(d)) < 1e-3,
          f"output crs: arrow bearing matches rotated u/v "
          f"(max {np.max(abs(d)):.3f} deg)", failures)
    check(np.allclose(c[:, 0, 2], g.z, atol=1e-3),
          "output crs: arrow z unchanged (m above sea level)", failures)
    with Dataset(out / "test_wind_z183m.nc") as nc:
        check(nc["u"].dimensions == ("time", "z", "lat", "lon")
              and nc["crs"].EPSG_code == "EPSG:4326",
              "output crs: wind netCDF on lat/lon in EPSG:4326", failures)

    # -- surfaces: exact transform, z untouched ---------------------------
    s = gpd.read_file(out / "test_surfaces.gpkg", layer="surfaces")
    ref = gpd.read_file(tmp / "out_new" / "test_surfaces.gpkg",
                        layer="surfaces")
    check(s.crs.to_epsg() == 4326 and len(s) == len(ref),
          "output crs: surfaces in EPSG:4326, same elements", failures)
    a = np.array(s.geometry.iloc[0].exterior.coords)
    b = np.array(ref.geometry.iloc[0].exterior.coords)
    bx, by = to_ll.transform(b[:, 0], b[:, 1])
    check(np.allclose(a[:, 0], bx, atol=1e-9)
          and np.allclose(a[:, 1], by, atol=1e-9)
          and np.allclose(a[:, 2], b[:, 2]),
          "output crs: surface vertices = UTM vertices transformed, z same",
          failures)
    shp = gpd.read_file(out / "test_surfaces.shp")
    check(shp.crs is not None and shp.crs.to_epsg() == 4326,
          "output crs: multipatch .prj is EPSG:4326", failures)

    # -- canopy: total leaf area conserved --------------------------------
    og = pg.reproj.grid()
    cell = og._w[2].reshape(og.shape)     # output cell areas, PALM cells
    with Dataset(out / "test_lad_3d.nc") as nc:
        lad = nc["lad"][:].filled(0.0).astype("f8")
        total = float((lad * cell).sum())
        check(nc["lad"].dimensions == ("time", "z", "lat", "lon")
              and abs(total - 27 * 0.5) < 1e-6 * 27 * 0.5,
              f"output crs: conservative remap keeps the total leaf area "
              f"(relative error {total / 13.5 - 1:.1e})", failures)

    # -- accuracy on analytic fields --------------------------------------
    # A constant must come back unchanged. A linear field comes back at each
    # output cell's centroid up to first-order remapping error: PALM cells
    # are treated as constant, which is exact only for axis-aligned grids
    # and at a 0.45 deg rotation leaves ~1e-4 of a cell.
    ny_c, nx_c = pg.nymax - pg.ny0, pg.nxmax - pg.nx0
    ones = np.ones((ny_c, nx_c), dtype=bool)
    c_out = og.take(np.full((ny_c, nx_c), 7.25), fill=np.nan, valid=ones)
    full = og._w[1].sum(axis=1).reshape(og.shape) >= \
        og._w[2].reshape(og.shape) * (1 - 1e-9)
    check(full.sum() > 200
          and np.abs(c_out[full] - 7.25).max() < 7.25 * 2e-7,
          "output crs: a constant field is reproduced (to the float32 "
          "precision the files store)", failures)
    jj_, ii_ = np.meshgrid(np.arange(ny_c) + 0.5, np.arange(nx_c) + 0.5,
                           indexing="ij")
    got = og.take(3.0 + 0.2 * ii_ + 0.1 * jj_, fill=np.nan, valid=ones)
    to_utm = Transformer.from_crs("EPSG:4326", CRS, always_xy=True)
    xe = og.x[0] + (og.x[1] - og.x[0]) * (np.arange(og.x.size + 1) - 0.5)
    ye = og.y[0] + (og.y[1] - og.y[0]) * (np.arange(og.y.size + 1) - 0.5)
    cx, cy = to_utm.transform(*np.meshgrid(xe, ye))
    ci = (np.asarray(cx) - ORIGIN_X) / DX
    cj = (np.asarray(cy) - ORIGIN_Y) / DY
    qi = 0.25 * (ci[:-1, :-1] + ci[:-1, 1:] + ci[1:, 1:] + ci[1:, :-1])
    qj = 0.25 * (cj[:-1, :-1] + cj[:-1, 1:] + cj[1:, 1:] + cj[1:, :-1])
    err = np.abs(got - (3.0 + 0.2 * qi + 0.1 * qj))[full] / 0.2
    check(err.max() < 1e-3,
          f"output crs: a linear field lands on the cell centroids "
          f"(max offset {err.max():.1e} of a cell)", failures)

    # -- scipy's sparse product and the numpy fallback agree -------------
    try:
        import scipy.sparse  # noqa: F401
        rnd = np.random.default_rng(3).random((4, ny_c, nx_c))
        air = np.random.default_rng(4).random((4, ny_c, nx_c)) > 0.2
        fast = og.take(rnd, fill=np.nan, valid=air)
        og._M, og._den_key = False, None
        slow = og.take(rnd, fill=np.nan, valid=air)
        og._M = None
        both = np.isfinite(fast) & np.isfinite(slow)
        check(np.array_equal(np.isfinite(fast), np.isfinite(slow))
              and np.abs(fast[both] - slow[both]).max() < 1e-6,
              "output crs: scipy and numpy remap paths agree", failures)
    except ImportError:
        print("  [skip] scipy not installed - sparse path not compared")

    # -- the nearest method is still there ---------------------------------
    cfg_nn = tmp / "config_4326_nn.yaml"
    cfg_nn.write_text(cfg.read_text()
                      .replace("./out_4326", "./out_4326_nn")
                      .replace("advanced:\n",
                               "advanced:\n    netcdf:\n"
                               "        resampling: nearest\n"))
    rc = pipeline_main(["-c", str(cfg_nn), "--overwrite",
                        "--stages", "canopy"])
    with Dataset(tmp / "out_4326_nn" / "test_lad_3d.nc") as nc:
        lad = nc["lad"][:].filled(0.0)
        vals = set(np.unique(lad[lad > 0]).round(6))
    check(rc == 0 and vals == {0.5},
          "output crs: resampling: nearest copies PALM's values unchanged",
          failures)
    return failures


def verify_log_last_run(tmp, pipeline_main):
    """<output_dir>/palm2gis.log holds the last run only, not every run."""
    failures = []
    cfg = str(tmp / "config.yaml")
    pipeline_main(["-c", cfg])
    pipeline_main(["-c", cfg])
    text = (tmp / "out" / "palm2gis.log").read_text(encoding="utf-8")
    check(text.count("palm2gis finished") == 1,
          "log: palm2gis.log holds only the last run", failures)
    return failures


def verify_float64_coords(tmp):
    """Absolute coordinates are float64. PALM stores x/y as float32 and
    float32 + origin stayed float32, rounding a UTM northing to 0.5 m."""
    from palm2gis.config import Config
    from palm2gis.steps.georef import process_static_driver
    failures = []
    g = process_static_driver(Config(str(tmp / "config.yaml")))
    y = g.abs_y(g.y_stat.astype("f4")) + 0.3
    check(y.dtype == "f8" and abs(float(y[0]) - (ORIGIN_Y + 1.3)) < 1e-6,
          "coords: absolute y keeps sub-0.5 m precision (float64)", failures)
    from netCDF4 import Dataset
    from palm2gis.steps.surfaces import _load_geometry
    with Dataset(tmp / "test_surf") as nc:
        xs, ys, zs = _load_geometry(nc)[:3]
    check(xs.dtype == ys.dtype == zs.dtype == "f8",
          "coords: surface xs/ys/zs read as float64", failures)
    return failures


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--keep", action="store_true",
                    help="keep the temporary directory for inspection")
    args = ap.parse_args()

    if not check_dependencies():
        return 1

    # A mismatched proj.db strips the CRS off every raster; repair it here
    # exactly as the pipeline does, so the test reflects a real run.
    from palm2gis.proj_env import ensure_proj_data
    from palm2gis.log import setup_logging
    setup_logging("info")
    if not ensure_proj_data():
        print("PROJ is unusable in this environment - the CRS checks below "
              "will fail. See the message above for the fix.")

    tmp = Path(tempfile.mkdtemp(prefix="palm2gis_test_"))
    print(f"Temporary directory: {tmp}")
    try:
        build_dataset(tmp)
        print("Synthetic dataset written; running the pipeline ...")

        from palm2gis.pipeline import main as pipeline_main
        rc = pipeline_main(["-c", str(tmp / "config.yaml")])
        if rc != 0:
            print("Pipeline returned a non-zero exit status.")
            return 1

        global _GRID_FOR_CHECKS
        from palm2gis.config import Config
        from palm2gis.steps.georef import process_static_driver
        _GRID_FOR_CHECKS = process_static_driver(Config(str(tmp /
                                                          "config.yaml")))

        print("Verifying the outputs ...")
        failures = verify_outputs(tmp / "out")

        print("Verifying variable descriptions ...")
        failures += verify_descriptions(tmp / "out")

        print("Verifying the overwrite policy ...")
        failures += verify_overwrite(tmp, pipeline_main)

        print("Verifying config-change detection ...")
        failures += verify_state(tmp, pipeline_main)

        print("Verifying handling of a fill-only time variable ...")
        failures += verify_regressions(tmp, tmp / "out", pipeline_main)

        print("Verifying the PALM basic state ...")
        failures += verify_basic_state(tmp)

        print("Verifying the config layout ...")
        failures += verify_layout(tmp, pipeline_main)
        failures += verify_t_surface_rename(tmp)
        failures += verify_fill_times(tmp, pipeline_main)

        print("Verifying the regular time axis ...")
        failures += verify_time_regularity()

        print("Verifying output in EPSG:4326 ...")
        failures += verify_output_crs(tmp, pipeline_main)
        failures += verify_float64_coords(tmp)
        failures += verify_log_last_run(tmp, pipeline_main)

        if failures:
            print(f"\nTEST FAILED: {len(failures)} check(s) failed.")
            return 1
        print("\nTEST OK: all checks passed.")
        return 0
    except Exception:
        traceback.print_exc()
        print("\nTEST FAILED with an exception (see traceback above).")
        return 1
    finally:
        if args.keep:
            print(f"Keeping {tmp}")
        else:
            shutil.rmtree(tmp, ignore_errors=True)


if __name__ == "__main__":
    sys.exit(main())
