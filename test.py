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

        def mk(name, dims, base):
            v = nc.createVariable(name, "f4", dims, fill_value=-9999.0)
            for t in range(NT):
                shp = tuple(len(nc.dimensions[d]) for d in dims[1:])
                data = base + 0.1 * t + 0.05 * rng.standard_normal(shp)
                m = np.zeros(shp, dtype=bool)
                m[:NZ, :NY, :NX] = kmask[:shp[0], :shp[1], :shp[2]]
                v[t] = np.ma.masked_where(m, data)

        mk("u", ("time", "zu_3d", "y", "xu"), 2.0)
        mk("v", ("time", "zu_3d", "yv", "x"), 1.0)
        mk("w", ("time", "zw_3d", "y", "x"), 0.1)
        mk("theta", ("time", "zu_3d", "y", "x"), 300.0)

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
  derive_ta: true
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
    import rasterio
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
        check(nc["time"].units.startswith("seconds since 2026-07-15"),
              "voxel: CF time from origin_time", failures)
        u = nc["u"][0]
        check(int((u == -9999.0).sum()) > 0,
              "voxel: terrain/building cells are NoData", failures)
        check(abs(float(u[u != -9999.0].mean()) - 2.0) < 0.1,
              "voxel: destaggered u close to the prescribed 2.0 m/s",
              failures)

    # -- derived air temperature -------------------------------------------
    # Exact analytic check: with theta constant, T must fall at g/cp.
    from palm2gis.steps.thermo import (CP, G, exner_profile, theta_to_ta)
    zz = _np_thermo.arange(0.0, 1001.0, 50.0)
    th_const = _np_thermo.full_like(zz, 300.0)
    T = theta_to_ta(th_const, exner_profile(zz, th_const, p_surface=1.0e5))
    lapse = (T[0] - T[-1]) / (zz[-1] - zz[0])
    check(abs(lapse - G / CP) < 1e-9,
          f"thermo: dry adiabatic lapse rate {lapse * 1000:.3f} K/km "
          f"= g/cp", failures)
    check(abs(exner_profile(zz, th_const, p_surface=1.0e5)[0] - 1.0) < 1e-12,
          "thermo: Exner = 1 at the reference pressure", failures)

    # p_surface applies at z = 0, not at the lowest OUTPUT level: a profile
    # over a cut starting at 100 m must agree with the full-column one there,
    # not repeat the surface value. Up to 0.12.0 the whole column was offset
    # by g/cp * z[0] - 1 K for this cut.
    cut = _np_thermo.arange(100.0, 1001.0, 50.0)
    pi_cut = exner_profile(cut, _np_thermo.full_like(cut, 300.0),
                           p_surface=1.0e5)
    pi_full = exner_profile(zz, th_const, p_surface=1.0e5)
    check(abs(pi_cut[0] - pi_full[2]) < 1e-12,
          "thermo: Exner profile is anchored at z = 0, not at the cut base",
          failures)

    # derived temperature must reach EVERY surface product, not just one
    with Dataset(out / "test_surface_2d.nc") as nc:
        check("t_surf" in nc.variables,
              "thermo: t_surf in the mesh NetCDF", failures)
    check((out / "test_t_surf_topview.tif").exists(),
          "thermo: t_surf raster written", failures)
    allel_ta = gpd.read_file(out / "test_surfaces.gpkg", layer="surfaces")
    check("t_surf" in allel_ta.columns
          and (allel_ta["t_surf"] > allel_ta["theta_surface"]).all(),
          "thermo: t_surf in the element layer (p_s > p0 -> T > theta)",
          failures)
    with Dataset(out / "test_voxel.nc") as nc:
        check("ta" in nc.variables
              and getattr(nc["ta"], "derived_from", "") == "theta",
              "thermo: derived ta in the voxel NetCDF, marked as derived",
              failures)

    # PALM's own variables must WIN over derivation
    from palm2gis.steps.surfaces import output_vars
    sf_native = {"variables": ["theta_surface"], "derive_ta": True}
    plan = output_vars(sf_native, {"theta_surface", "t_surf"})
    check(("t_surf", "t_surf", False) in plan,
          "thermo: native t_surf preferred over deriving it", failures)
    plan = output_vars(sf_native, {"theta_surface"})
    check(("t_surf", "theta_surface", True) in plan,
          "thermo: t_surf derived when the file lacks it", failures)
    plan = output_vars({"variables": ["theta_surf"], "derive_ta": True},
                       {"theta_surf"})
    check(("t_surf", "theta_surf", True) in plan,
          "thermo: the theta_surf spelling is accepted too", failures)

    # -- surfaces ----------------------------------------------------------
    with rasterio.open(out / "test_theta_surface_topview.tif") as r:
        check(r.count == NT, "surfaces: one raster band per timestep",
              failures)
        check(r.crs is not None, "surfaces: the raster carries a CRS",
              failures)
        check(r.bounds.left == ORIGIN_X and r.bounds.bottom == ORIGIN_Y,
              "surfaces: raster georeferencing at the domain origin",
              failures)
        b = r.read(1)
        check(np.allclose(b[b != r.nodata], 300.0),
              "surfaces: top-view values (upward theta_surface = 300)",
              failures)

    # CF NetCDF for QGIS mesh layers: time axis read natively by MDAL
    with Dataset(out / "test_surface_2d.nc") as nc:
        check(nc["time"].units.startswith("seconds since 2026-07-15"),
              "surfaces: mesh NetCDF carries a CF time axis", failures)
        check(nc["theta_surface"].shape == (NT, 20, 20),
              "surfaces: mesh NetCDF is (time, y, x)", failures)
        check(float(nc["x"][0]) == ORIGIN_X + DX / 2.0,
              "surfaces: mesh NetCDF x in absolute coordinates", failures)

    walls = gpd.read_file(out / "test_walls.gpkg", layer="walls")
    check(walls.geometry.iloc[0].has_z, "surfaces: walls are PolygonZ",
          failures)
    check(walls["time"].nunique() == NT,
          "surfaces: walls carry all timesteps", failures)
    # Exact, not floor-divided: the GeoPackage is written one timestep at a
    # time and appended, so a botched append would either duplicate rows or
    # drop them. 32 wall faces x NT timesteps, no more and no less.
    check(len(walls) == 32 * NT,
          f"surfaces: walls has exactly {32 * NT} rows, got {len(walls)}",
          failures)
    check(walls.groupby("rec").size().nunique() == 1,
          "surfaces: every timestep contributed the same row count",
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

    # -- QGIS style sidecars ------------------------------------------------
    import sqlite3
    import xml.dom.minidom as minidom

    qml = out / "test_theta_surface_topview.qml"
    check(qml.exists(), "qgis: raster style sidecar written", failures)
    if qml.exists():
        doc = minidom.parse(str(qml))          # also proves well-formedness
        temporal = doc.getElementsByTagName("temporal")[0]
        ranges = doc.getElementsByTagName("range")
        check(temporal.getAttribute("mode") == "3" and len(ranges) == NT,
              "qgis: raster QML maps every band to a time range", failures)
        check(ranges[0].getAttribute("begin") == "2026-07-15T10:00:00Z",
              "qgis: raster band times are explicit UTC", failures)

    vqml = out / "test_wind.qml"
    check(vqml.exists(), "qgis: vector style sidecar written", failures)
    if vqml.exists():
        t = minidom.parse(str(vqml)).getElementsByTagName("temporal")[0]
        check(t.getAttribute("mode") == "1"
              and t.getAttribute("startField") == "time",
              "qgis: vector QML reads time from the time field", failures)

    # the repair path must recover a raster that lands without a CRS
    import numpy as _np
    import rasterio as _rio
    from rasterio.transform import from_origin as _fo
    from palm2gis.steps.surfaces import _ensure_raster_crs

    class _Cfg:
        def disp(self, p):
            return Path(p).name

    probe = out / "_crs_probe.tif"
    with _rio.open(probe, "w", driver="GTiff", height=2, width=2, count=1,
                   dtype="float32", transform=_fo(0, 2, 1, 1)) as d:
        d.write(_np.zeros((2, 2), "float32"), 1)      # written with NO crs
    repaired = _ensure_raster_crs(_Cfg(), _GRID_FOR_CHECKS, probe)
    with _rio.open(probe) as r:
        check(repaired and r.crs is not None,
              "qgis: a raster written without a CRS is repaired in place",
              failures)
    probe.unlink()

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
    # 3 x 3 columns, 3 levels of LAD = 0.5 on dz_lad = 2 m.
    with rasterio.open(out / "test_lai.tif") as r:
        v = r.read(1)[r.read(1) != r.nodata]
        check(v.size == 9 and _np_thermo.allclose(v, 0.5 * 3 * DZ),
              f"canopy: LAI = LAD x depth ({v[0]:.2f} over {v.size} cells)",
              failures)
        check(r.crs is not None, "canopy: LAI raster carries a CRS",
              failures)
    with rasterio.open(out / "test_canopy_height.tif") as r:
        v = r.read(1)[r.read(1) != r.nodata]
        check(_np_thermo.allclose(v, 3 * DZ),
              f"canopy: canopy height = {v[0]:.1f} m above terrain",
              failures)
    with Dataset(out / "test_lad_3d.nc") as nc:
        nc.set_auto_mask(False)
        lad = nc["lad"][:]
        check(int((lad > 0).sum()) == 27,
              "canopy: LAD volume holds every canopy cell", failures)
        # the test terrain slopes across the tree block, so the SAME
        # canopy must land on different absolute levels per column -
        # which is exactly what the height-above-terrain lift is for
        zocc = _np_thermo.asarray(nc["z"][:])[
            (lad > 0).any(axis=(1, 2))]
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
def verify_fill_times(tmp, pipeline_main):
    """Some PALM surface files leave the time variable at its _FillValue
    (~1e36). Those values must never be treated as seconds - they overflow
    datetime - and the stage must still produce usable layers."""
    import geopandas as gpd
    import rasterio
    from netCDF4 import Dataset

    failures = []
    with Dataset(tmp / "test_surf", "a") as nc:
        nc["time"][:] = 9.969209968386869e36

    # (a) no reconstruction configured: record-index labels, empty times
    rc = pipeline_main(["-c", str(tmp / "config.yaml"), "--stages",
                        "surfaces", "--overwrite"])
    check(rc == 0, "fill times: the stage completes instead of overflowing",
          failures)
    with rasterio.open(tmp / "out" / "test_theta_surface_topview.tif") as r:
        check(r.descriptions[0] == "record 0",
              "fill times: raster bands fall back to record labels",
              failures)
    walls = gpd.read_file(tmp / "out" / "test_walls.gpkg", layer="walls")
    check(walls["time"].isna().all() and "rec" in walls.columns,
          "fill times: walls carry an empty time and a rec column", failures)

    # (b) the p3d supplies the interval - no manual settings at all
    (tmp / "p3d_probe").write_text(
        "&initialization_parameters\n nz = 12,\n dz = 2.0,\n/\n"
        "&runtime_parameters\n dt_dosurf_av = 600.0,\n"
        " skip_time_dosurf_av = 0.0,\n/\n")
    cfgp0 = tmp / "config.yaml"
    base = cfgp0.read_text()
    cfgp0.write_text(base.replace("  nz: 12\n  dz: 2.0", "  p3d: p3d_probe"))
    pipeline_main(["-c", str(cfgp0), "--stages", "surfaces", "--overwrite"])
    with rasterio.open(tmp / "out" / "test_theta_surface_topview.tif") as r:
        check(r.descriptions[1] == "2026-07-15 10:20:00",
              "fill times: interval read from the p3d automatically",
              failures)
    cfgp0.write_text(base)

    # (c) bounds.time_step given: real timestamps reconstructed
    cfgp = tmp / "config.yaml"
    original = cfgp.read_text()
    cfgp.write_text(original.replace("bounds:\n  zmax: 10",
                                     "bounds:\n  zmax: 10\n"
                                     "  time_step: 600\n  time_start: 0.0"))
    pipeline_main(["-c", str(cfgp), "--stages", "surfaces", "--overwrite"])
    with rasterio.open(tmp / "out" / "test_theta_surface_topview.tif") as r:
        check(r.descriptions[1] == "2026-07-15 10:10:00",
              "fill times: bounds.time_step reconstructs real timestamps",
              failures)
    cfgp.write_text(original)
    return failures


# ------------------------------
# 7. MAIN
# ------------------------------
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

        print("Verifying the overwrite policy ...")
        failures += verify_overwrite(tmp, pipeline_main)

        print("Verifying config-change detection ...")
        failures += verify_state(tmp, pipeline_main)

        print("Verifying handling of a fill-only time variable ...")
        failures += verify_fill_times(tmp, pipeline_main)

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
