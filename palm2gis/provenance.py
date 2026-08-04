"""
palm2gis.provenance
-------------------
Records HOW each output file was produced.

A GeoTIFF or GeoPackage carries no room for a run history the way a NetCDF
does, and six months later an output directory is anonymous: the filename
says which variable it holds and nothing about which bounds, which config
or which version produced it. Two runs that differ only in a setting are
indistinguishable on disk.

So every output gets a `<name>.palm2gis.json` sidecar, and formats that can
hold metadata natively get it there too:

  * GeoTIFF - written into the file's own metadata tags, so it survives
    being copied out of the directory and is visible in `gdalinfo` and in
    QGIS's layer properties.
  * NetCDF  - global attributes plus a CF-style `history` line, appended to
    rather than replacing any history already there.
  * GeoPackage - sidecar only (GDAL's metadata support for GPKG is patchy
    and writing it can rewrite the file).
"""

from __future__ import annotations

import datetime
import json
import shlex
import sys
from pathlib import Path

from .log import debug, warning


def provenance(cfg, stage: str, src=None) -> dict:
    """The provenance record for one output file."""
    from . import __version__
    from .state import stage_hash

    stamp = datetime.datetime.now().astimezone().isoformat(timespec="seconds")
    rec = {
        "version": __version__,
        "config": str(getattr(cfg, "_path", "") or ""),
        "config_sha": stage_hash(cfg, stage),
        "project": str(cfg["project"].get("name", "")),
        "stage": stage,
        "command": " ".join(shlex.quote(a) for a in sys.argv),
        "written": stamp,
    }
    if src:
        rec["source"] = str(src)
    try:
        rec["bounds"] = {k: v for k, v in dict(cfg["bounds"]).items()
                         if v is not None}
    except Exception:
        pass
    return rec


def _history_line(rec: dict) -> str:
    return f"{rec['written']}: palm2gis {rec['version']}: {rec['command']}"


def stamp_file(cfg, path, stage: str, src=None) -> None:
    """Attach provenance to *path*: sidecar always, native where possible."""
    path = Path(path)
    if not path.exists():
        return
    rec = provenance(cfg, stage, src)

    side = path.with_suffix(path.suffix + ".palm2gis.json")
    try:
        side.write_text(json.dumps(rec, indent=2, sort_keys=True))
        debug("[provenance] wrote {}", cfg.disp(side))
    except OSError as exc:
        warning("[provenance] could not write {} ({})", side, exc)

    suffix = path.suffix.lower()
    try:
        if suffix in (".tif", ".tiff"):
            import rasterio
            with rasterio.open(path, "r+") as dst:
                dst.update_tags(**{f"PALM2GIS_{k.upper()}": str(v)
                                   for k, v in rec.items()
                                   if not isinstance(v, dict)})
        elif suffix == ".nc":
            from netCDF4 import Dataset
            with Dataset(path, "a") as nc:
                for k, v in rec.items():
                    if isinstance(v, dict):
                        continue
                    nc.setncattr(f"palm2gis_{k}", str(v))
                prev = ""
                if "history" in nc.ncattrs():
                    prev = str(nc.getncattr("history")).rstrip()
                line = _history_line(rec)
                nc.setncattr("history",
                             f"{prev}\n{line}".strip() if prev else line)
    except Exception as exc:
        # Never fail a run over metadata: the sidecar is already written.
        debug("[provenance] native metadata not written for {} ({})",
              path.name, exc)
