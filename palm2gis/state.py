"""
palm2gis.state
--------------
Remembers WHICH CONFIG produced each output file.

Without this, the skip in `Config.prepare_output` is by file existence
alone: change `bounds`, `celsius`, `heights` or a variable list, rerun, and
every existing output is kept because the filename did not change. The
README used to have to warn about it and tell you to pass `--overwrite`
after editing bounds - a documented footgun that this removes.

The state is a small JSON sidecar in the output directory mapping each
output path to the hash of the settings that produced it. A file whose
recorded hash differs from the current one is treated as stale and
rewritten; a file with no record is left alone, so an existing output
directory from an older version is not invalidated wholesale on the first
run.

Only the settings a stage actually reads go into its hash, so editing the
wind heights does not invalidate the surface rasters.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

from .log import debug, warning

STATE_FILE = ".palm2gis_state.json"

# Which config sections each stage's output depends on. `domain` and
# `bounds` affect every product; the rest are per stage.
_STAGE_KEYS = {
    # surface_pressure and pt_surface feed the theta -> ta derivation, so a
    # change to either must invalidate the products that carry derived
    # temperatures.
    "voxel":    ("domain", "bounds", "crs", "celsius", "surface_pressure",
                 "pt_surface", "voxel"),
    # qgis_styles decides what goes INTO the GeoPackages (an embedded
    # layer_styles table), so toggling it must rewrite them.
    "surfaces": ("domain", "bounds", "crs", "celsius", "surface_pressure",
                 "pt_surface", "qgis_styles", "qgis_styles_3d", "surfaces"),
    "wind":     ("domain", "bounds", "crs", "qgis_styles", "qgis_styles_3d",
                 "wind"),
    "canopy":   ("domain", "bounds", "crs", "canopy"),
}


def _stable(obj):
    """JSON-serialisable form with stable ordering."""
    if isinstance(obj, dict):
        return {k: _stable(obj[k]) for k in sorted(obj)}
    if isinstance(obj, (list, tuple)):
        return [_stable(v) for v in obj]
    if isinstance(obj, Path):
        return str(obj)
    return obj


def stage_hash(cfg, stage: str) -> str:
    """Hash of the config settings that stage's outputs depend on."""
    keys = _STAGE_KEYS.get(stage, (stage,))
    payload = {}
    for k in keys:
        try:
            payload[k] = _stable(cfg[k])
        except Exception:
            payload[k] = None
    blob = json.dumps(payload, sort_keys=True, default=str)
    return hashlib.sha256(blob.encode("utf-8")).hexdigest()[:16]


class State:
    """Output-path -> config-hash record, persisted as JSON."""

    def __init__(self, output_dir, enabled: bool = True):
        self.path = Path(output_dir) / STATE_FILE
        self.enabled = bool(enabled)
        self._data = {}
        if self.enabled and self.path.is_file():
            try:
                self._data = json.loads(self.path.read_text())
            except (OSError, ValueError) as exc:
                warning("[state] could not read {} ({}) - treating all "
                        "existing outputs as unknown.", self.path, exc)
                self._data = {}

    def _key(self, path) -> str:
        try:
            return str(Path(path).resolve().relative_to(
                self.path.parent.resolve()))
        except ValueError:
            return str(Path(path).resolve())

    def is_stale(self, path, want_hash: str) -> bool:
        """True when *path* was produced by a DIFFERENT config.

        An unrecorded file is not stale: an output directory written before
        the state file existed should not be wiped on the first run.
        """
        if not self.enabled:
            return False
        rec = self._data.get(self._key(path))
        if rec is None:
            return False
        return rec.get("hash") != want_hash

    def record(self, path, want_hash: str) -> None:
        if not self.enabled:
            return
        self._data[self._key(path)] = {"hash": want_hash}

    def forget(self, path) -> None:
        self._data.pop(self._key(path), None)

    def save(self) -> None:
        if not self.enabled:
            return
        try:
            self.path.parent.mkdir(parents=True, exist_ok=True)
            self.path.write_text(json.dumps(self._data, indent=2,
                                            sort_keys=True))
            debug("[state] wrote {} ({} entries)", self.path,
                  len(self._data))
        except OSError as exc:
            warning("[state] could not write {} ({}) - the next run will "
                    "not detect config changes.", self.path, exc)
