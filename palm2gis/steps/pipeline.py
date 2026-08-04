"""palm2gis pipeline: stage orchestration and CLI.

Usage (via the root entry point):
    python run_palm2gis.py -c config/<project>.yaml [flags]

Stages:
    voxel     DATA_3D_NETCDF* -> CF-compliant NetCDF for the ArcGIS Pro
              Voxel Layer (destaggered u/v/w + scalars, terrain as NoData)
    surfaces  SURFACE_DATA_NETCDF* -> top-view GeoTIFF (topmost roof/ground
              element per column, band = timestep) + wall quads (PolygonZ)
              in a time-enabled GeoPackage
    wind      DATA_3D_NETCDF* -> wind vector point layer (u, v, speed,
              dir_met) at heights above the local surface, GeoPackage

The static driver is read first in every run - it provides the grid, the
georeferencing (origin_x/origin_y in crs.palm) and the bounds index ranges
shared by all stages.
"""

import argparse
import sys
from pathlib import Path

from .config import Config, ConfigError, STAGES
from .log import (emit_to_file_only, get_logger, progress,
                  setup_logging)
from .proj_env import ensure_proj_data
from .steps.georef import process_static_driver
from .steps.voxel import stage_voxel
from .steps.surfaces import stage_surfaces
from .steps.wind import stage_wind
from .steps.canopy import stage_canopy

log = get_logger()


# ------------------------------
# 1. STAGE TABLE
# ------------------------------
STAGE_LABELS = {
    "voxel": "Writing the ArcGIS voxel NetCDF",
    "surfaces": "Converting the surface output",
    "wind": "Sampling wind vectors",
    "canopy": "Converting the plant canopy",
}

STAGE_FUNCS = {"voxel": stage_voxel, "surfaces": stage_surfaces,
               "wind": stage_wind, "canopy": stage_canopy}


# ------------------------------
# 2. DRY-RUN PLAN
# ------------------------------
def stage_outputs(cfg, stage):
    """Output paths a stage would produce, in the order it writes them."""
    if stage == "voxel":
        return [cfg.out_path(cfg["voxel"]["file_to"])]
    if stage == "surfaces":
        sf = cfg["surfaces"]
        out = []
        if sf["write_topview_raster"]:
            from .steps.surfaces import output_vars
            # Plan only: without reading the file, assume the derived name.
            out += [cfg.out_path("{}{}_topview.tif".format(sf["prefix"], n))
                    for n, _v, _d in output_vars(sf)]
        if sf["write_topview_netcdf"]:
            out.append(cfg.out_path("{}surface_2d.nc".format(sf["prefix"])))
        if sf["write_walls_gpkg"]:
            out.append(cfg.out_path("{}walls.gpkg".format(sf["prefix"])))
        if sf["write_all_elements_gpkg"]:
            out.append(cfg.out_path("{}surfaces.gpkg".format(sf["prefix"])))
        return out
    if stage == "wind":
        return [cfg.out_path(cfg["wind"]["file_to"])]
    if stage == "canopy":
        cn = cfg["canopy"]
        out = []
        if cn["write_lai"]:
            out.append(cfg.out_path("{}lai.tif".format(cn["prefix"])))
        if cn["write_canopy_height"]:
            out.append(cfg.out_path(
                "{}canopy_height.tif".format(cn["prefix"])))
        if cn["write_volume"]:
            out.append(cfg.out_path("{}lad_3d.nc".format(cn["prefix"])))
        return out
    return []


def print_plan(cfg, stages):
    log.info("Plan (dry run - nothing written):")
    log.info(f"  static driver: {cfg.resolve(cfg['domain']['static_driver'])}")
    log.info(f"  output dir:    {cfg.output_dir}")
    log.info(f"  overwrite:     {cfg.overwrite}")
    for stage in stages:
        if not cfg[stage]["enabled"]:
            log.info(f"  {stage}: disabled")
            continue
        log.info(f"  {stage}: from "
                 f"{cfg[stage]['file_from'] or cfg['domain']['static_driver']}")
        for path in stage_outputs(cfg, stage):
            if path.exists() and not cfg.overwrite:
                action = "keep (exists)"
            elif path.exists():
                action = "REPLACE"
            else:
                action = "write"
            log.info(f"      {action}: {path}")


# ------------------------------
# 3. CLI ENTRY
# ------------------------------
def main(argv=None):
    ap = argparse.ArgumentParser(prog="palm2gis",
                                 description="PALM outputs -> ArcGIS Pro / "
                                             "QGIS layers.")
    ap.add_argument("-c", "--config", required=True, help="YAML config file")
    ap.add_argument("--stages", nargs="+", default=None,
                    help="subset of stages to run (default: from config)")
    ap.add_argument("--overwrite", action="store_true",
                    help="replace existing outputs (default: keep them and "
                         "write only what is missing)")
    ap.add_argument("-n", "--dry-run", action="store_true",
                    help="validate config and print plan, write nothing")
    ap.add_argument("-v", "--verbose", action="store_true",
                    help="DEBUG-level logging (per-timestep progress)")
    ap.add_argument("-q", "--quiet", action="store_true",
                    help="WARNING-level logging only")
    ap.add_argument("--log-file", default=None, metavar="PATH",
                    help="also append the log to this file (uncoloured, "
                         "always dated); default <output_dir>/palm2gis.log, "
                         "disable with --log-file ''")
    # argparse's own "version" action so it works WITHOUT -c: the flag is
    # useless if it still demands a config file.
    from . import __version__ as _ver
    ap.add_argument("--version", action="version",
                    version="palm2gis {}".format(_ver),
                    help="print the palm2gis version and exit")
    ap.add_argument("--ignore-state", action="store_true",
                    help="skip purely on file existence, ignoring the "
                         "recorded config hash (pre-0.10.0 behaviour)")
    ap.add_argument("--log-datetime", action="store_true",
                    help="full date in timestamps (time-only is the default)")
    args = ap.parse_args(argv)

    verbosity = "debug" if args.verbose else ("warning" if args.quiet
                                              else "info")
    from . import __version__

    setup_logging(verbosity, args.log_datetime, args.log_file or None)
    log.info(f"palm2gis {__version__}")

    try:
        cfg = Config(args.config)
    except (ConfigError, FileNotFoundError) as exc:
        log.error(f"Config error: {exc}")
        return 1

    # The output dir is only known once the config is read, so the default
    # log file is attached here. Re-initialising clears the handlers, hence
    # the version line is re-emitted into the freshly opened file.
    if args.log_file is None:
        setup_logging(verbosity, args.log_datetime,
                      cfg.output_dir / "palm2gis.log")
        # Into the FILE only: stderr already showed it above.
        emit_to_file_only(f"palm2gis {__version__}")

    # Records which config produced each output, so a settings change
    # invalidates the affected files instead of silently keeping them.
    from .state import State
    cfg._state = State(cfg.output_dir, enabled=not args.ignore_state)

    if args.overwrite:
        cfg["project"]["overwrite"] = True
    if cfg.overwrite:
        log.warning("project.overwrite: existing outputs will be REPLACED.")

    stages = args.stages or cfg["stages"]
    for stage in stages:
        if stage not in STAGES:
            log.error(f"Unknown stage: {stage} "
                      f"(known: {', '.join(STAGES)}).")
            return 1

    # A mismatched PROJ database silently strips the CRS from every
    # raster, so check (and repair) it before any work is done.
    ensure_proj_data()

    progress("Reading configuration")
    # Logged once, so every later "wrote ..." line can use a bare filename.
    log.info(f"output dir: {cfg.output_dir}")
    log.debug(f"Project: {cfg['project']['name']}")
    log.debug(f"Stages: {' -> '.join(stages)}")

    if args.dry_run:
        print_plan(cfg, stages)
        return 0

    try:
        grid = process_static_driver(cfg)
    except Exception as exc:
        log.error(f"Static driver failed: {exc}")
        return 1

    for stage in stages:
        if not cfg[stage]["enabled"]:
            log.info(f"Stage {stage}: disabled, skipping.")
            continue
        progress(STAGE_LABELS.get(stage, stage))
        try:
            STAGE_FUNCS[stage](cfg, grid)
        except Exception as exc:
            log.error(f"Stage {stage} failed: {exc}")
            return 1
        # Record the config behind whatever that stage actually produced,
        # and stamp each file so it can be identified on its own later.
        from .provenance import stamp_file
        # `canopy` reads the static driver when it has no file_from of its
        # own, so the source is resolved only once it is known to exist.
        src_cfg = (cfg[stage].get("file_from")
                   or (cfg["domain"]["static_driver"]
                       if stage == "canopy" else None))
        for path in stage_outputs(cfg, stage):
            if Path(path).exists() and cfg.was_written(path):
                cfg.note_written(path, stage)
                stamp_file(cfg, path, stage,
                           src=cfg.resolve(src_cfg) if src_cfg else None)
        cfg._state.save()

    progress("palm2gis finished OK")
    return 0


if __name__ == "__main__":
    sys.exit(main())
