#!/usr/bin/env python3
"""palm2gis diagnostics.

Prints a compact report on the inputs and the produced layers, so a bad
conversion can be diagnosed without moving large files around. Reads only;
nothing is written or modified.

    python diagnose.py -c config/<project>.yaml > diag.txt

Paste diag.txt (it is a few dozen lines) when reporting a problem.
"""

# ------------------------------
# 1. INPUT AND OUTPUT FILES
# ------------------------------
import argparse
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import numpy as np
from netCDF4 import Dataset

SEP = "-" * 62


def head(title):
    print("\n" + SEP)
    print(title)
    print(SEP)


def stat(name, a):
    """Robust min/mean/max of a possibly masked / fill-laden array."""
    a = np.ma.masked_invalid(np.ma.masked_equal(np.ma.asarray(a, dtype="f8"),
                                                -9999.0))
    if a.count() == 0:
        print(f"    {name:<22} all NoData ({a.size} cells)")
        return
    print(f"    {name:<22} min {a.min():12.4f}  mean {a.mean():12.4f}  "
          f"max {a.max():12.4f}   NoData {a.size - a.count()}/{a.size}")


# ------------------------------
# 2. CONFIG AND GRID
# ------------------------------
def report_config(cfg):
    head("CONFIG")
    print(f"  root         {cfg.root}")
    print(f"  output_dir   {cfg.output_dir}")
    print(f"  crs.palm     {cfg['crs']['palm']}")
    print(f"  overwrite    {cfg.overwrite}")
    print(f"  bounds       {dict(cfg['bounds']._settings)
                            if hasattr(cfg['bounds'], '_settings')
                            else cfg['bounds']}")
    for st in ("voxel", "surfaces", "wind"):
        s = cfg[st]
        if s["enabled"]:
            print(f"  {st:<12} from {s['file_from']}")
            for k in ("scalars", "variables", "grid_heights", "grid_stride",
                      "drop_k0"):
                if k in s:
                    print(f"               {k}: {s[k]}")


def report_static(cfg, grid):
    head("STATIC DRIVER / GRID")
    print(f"  nx x ny x nz   {grid.nx} x {grid.ny} x {grid.nz}")
    print(f"  dx/dy/dz       {grid.dx} / {grid.dy} / {grid.dz} m")
    print(f"  origin         ({grid.origin_x:.3f}, {grid.origin_y:.3f})")
    print(f"  origin_time    {grid.origin_time}")
    print(f"  cut indices    x[{grid.nx0}:{grid.nxmax}] "
          f"y[{grid.ny0}:{grid.nymax}] z[{grid.nz0}:{grid.nzmax}]")
    print(f"  cut extent     x {grid.x_stat[0]:.1f}..{grid.x_stat[-1]:.1f} m,"
          f"  y {grid.y_stat[0]:.1f}..{grid.y_stat[-1]:.1f} m (local)")
    print(f"  absolute       x {grid.abs_x(grid.x_stat)[0]:.1f}.."
          f"{grid.abs_x(grid.x_stat)[-1]:.1f},  "
          f"y {grid.abs_y(grid.y_stat)[0]:.1f}.."
          f"{grid.abs_y(grid.y_stat)[-1]:.1f}")
    stat("terrain zt", grid.zt)
    stat("terrain+buildings", grid.elev)
    print(f"  z levels (cut) {np.array2string(grid.z_levels()[grid.nz0:grid.nzmax][:6], precision=1)}"
          f" ... {grid.z_levels()[grid.nzmax - 1]:.1f}")

    with Dataset(cfg.resolve(cfg["domain"]["static_driver"])) as nc:
        print(f"  static vars    {', '.join(sorted(nc.variables))}")
        for a in ("origin_x", "origin_y", "origin_z", "origin_time",
                  "origin_lat", "origin_lon", "rotation_angle"):
            if hasattr(nc, a):
                print(f"    attr {a:<14} {getattr(nc, a)}")


# ------------------------------
# 3. PALM OUTPUT FILES
# ------------------------------
def report_input_file(cfg, path, label):
    head(f"INPUT: {label}  ({path})")
    p = cfg.resolve(path)
    if not Path(p).exists():
        print("  MISSING")
        return
    with Dataset(p) as nc:
        print("  dims: " + ", ".join(
            f"{d}={len(nc.dimensions[d])}" for d in nc.dimensions))
        t = nc.variables.get("time")
        if t is not None:
            raw = np.ma.getdata(t[:]).ravel()
            print(f"  time: {raw.size} records, first {raw[:3]}, "
                  f"last {raw[-1]}")
            print(f"        _FillValue "
                  f"{getattr(t, '_FillValue', None)}, masked "
                  f"{int(np.ma.getmaskarray(t[:]).sum())}")
        print("  variables:")
        for name, v in nc.variables.items():
            if name in ("time",):
                continue
            print(f"    {name:<18} {str(v.dimensions):<42} "
                  f"{getattr(v, 'units', '')}")
        # coordinate sanity
        for c in ("x", "y", "xu", "yv", "zu_3d", "zw_3d", "xs", "ys", "zs"):
            if c in nc.variables:
                a = np.asarray(nc.variables[c][:]).ravel()
                print(f"  coord {c:<6} n={a.size:<8} "
                      f"{a.min():.2f} .. {a.max():.2f}"
                      + (f"  step {a[1] - a[0]:.3f}" if a.size > 1 else ""))
        if "zenith" in nc.variables:
            z = np.ma.filled(nc.variables["zenith"][:], -9999.0)
            vals, cnt = np.unique(np.round(z, 1), return_counts=True)
            print("  zenith values: " + ", ".join(
                f"{v}({c})" for v, c in zip(vals, cnt)))
        if "azimuth" in nc.variables:
            a = np.ma.filled(nc.variables["azimuth"][:], -9999.0)
            vals, cnt = np.unique(np.round(a, 1), return_counts=True)
            print("  azimuth values: " + ", ".join(
                f"{v}({c})" for v, c in zip(vals, cnt)))
        # first-timestep value ranges for the configured variables
        for name in list(nc.variables):
            v = nc.variables[name]
            if v.ndim >= 3 and "time" in v.dimensions:
                try:
                    stat(name + " [t0]", v[0, ...])
                except Exception as exc:
                    print(f"    {name}: could not read ({exc})")


# ------------------------------
# 4. PRODUCED LAYERS
# ------------------------------
def report_outputs(cfg):
    head("OUTPUT LAYERS")
    import geopandas as gpd

    sf = cfg["surfaces"]
    checked = []

    if cfg["voxel"]["enabled"]:
        p = cfg.out_path(cfg["voxel"]["file_to"])
        checked.append(p)
        if p.exists():
            print(f"\n  {p.name}")
            with Dataset(p) as nc:
                print("    dims: " + ", ".join(
                    f"{d}={len(nc.dimensions[d])}" for d in nc.dimensions))
                print(f"    x {float(nc['x'][0]):.2f}..{float(nc['x'][-1]):.2f}"
                      f"  y {float(nc['y'][0]):.2f}..{float(nc['y'][-1]):.2f}"
                      f"  z {float(nc['z'][0]):.2f}..{float(nc['z'][-1]):.2f}")
                print(f"    time units: {getattr(nc['time'], 'units', '-')}, "
                      f"values {np.asarray(nc['time'][:])[:4]}")
                nc.set_auto_mask(False)
                for name in nc.variables:
                    if nc.variables[name].ndim == 4:
                        stat(name, nc.variables[name][0])

    if sf["enabled"] and sf["write_all_elements_gpkg"]:
        p = cfg.out_path(sf["all_elements_to"]
                         or f"{sf['prefix']}surfaces.gpkg")
        checked.append(p)
        if p.exists():
            g = gpd.read_file(p, layer=sf["layer_name"] or "surfaces")
            print(f"\n  {p.name}: {len(g)} rows, crs {g.crs.to_epsg()}, "
                  f"{g['time'].nunique()} times")
            print(f"    bounds {tuple(round(v, 1) for v in g.total_bounds)}")
            for c in g.columns:
                if (c != "geometry"
                        and gpd.pd.api.types.is_numeric_dtype(g[c])):
                    stat(c, g[c].values)

    missing = [p for p in checked if not p.exists()]
    if missing:
        print("\n  MISSING outputs: " + ", ".join(p.name for p in missing))


# ------------------------------
# 5. MAIN
# ------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("-c", "--config", required=True)
    args = ap.parse_args()

    from palm2gis import __version__
    from palm2gis.config import Config
    from palm2gis.log import setup_logging
    from palm2gis.steps.georef import process_static_driver

    setup_logging("warning")
    print(f"palm2gis {__version__} diagnostics")
    cfg = Config(args.config)

    report_config(cfg)
    grid = process_static_driver(cfg)
    report_static(cfg, grid)

    seen = set()
    for st in ("voxel", "surfaces", "wind"):
        if cfg[st]["enabled"] and cfg[st]["file_from"] not in seen:
            seen.add(cfg[st]["file_from"])
            report_input_file(cfg, cfg[st]["file_from"], st)

    report_outputs(cfg)
    print("\nDone.")


if __name__ == "__main__":
    main()
