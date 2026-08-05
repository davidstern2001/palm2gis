"""YAML config loading + defaults for palm2gis."""

import copy
import os
import re
from contextlib import contextmanager
from pathlib import Path

import yaml

from .log import debug, get_logger, info, warning

log = get_logger()

# ------------------------------
# 1. STAGE REGISTRY
# ------------------------------
STAGES = ("voxel", "surfaces", "wind", "canopy")

# ------------------------------
# 2. DEFAULTS
# ------------------------------
DEFAULTS = {
    "project": {
        "name": "palm2gis",
        # Base directory for ALL relative paths in this config. Default (null)
        # is the directory of the YAML file itself; set this when configs are
        # stored centrally (e.g. in config/) away from the simulation data.
        "root": None,
        "output_dir": ".",
        # false (default): an output file that already exists is KEPT and its
        # product is skipped, so a rerun only fills in what is missing.
        # true: regenerate everything, replacing existing outputs.
        "overwrite": False,
    },
    "crs": {
        # CRS of the PALM coordinate system. origin_x/origin_y from the
        # static driver are interpreted in this CRS and every output layer
        # is written in it (absolute coordinates - no manual shifting in
        # ArcGIS/QGIS).
        "palm": "EPSG:32633",
    },
    "domain": {
        "static_driver": None,     # [required] static driver nc file
        "p3d": None,               # optional: read nz/dz (and dx/dy) from
                                   # the _p3d namelist
        "nz": None,                # required when p3d is not given
        "dz": None,                # required when p3d is not given
        # Override the static driver's origin_time attribute,
        # format "YYYY-MM-DD HH:MM:SS" (any UTC-offset suffix is ignored).
        "origin_time": None,
    },
    "bounds": {
        # Subdomain cut in metres, LOCAL PALM coordinates. Unset = no cut.
        # x0 / xmax / y0 / ymax
        # zmax: height above the highest surface (terrain + buildings)
        #       inside the cut subdomain.
        # ---- window in TIME (preferred) ----------------------------------
        # The period to convert, stated ONCE and resolved against each
        # file's own time axis. PALM keeps separate output intervals for 3D
        # and surface output, so the same RECORD INDEX means a different
        # moment in each file - record 18 is 3 hours in for a 10-minute 3D
        # file and 18 hours in for an hourly surface file. Set these and
        # every stage converts the same period. Accepted:
        #   3600                  seconds since origin_time
        #   "6h" / "90min"        the same, with a unit (s/min/h/d)
        #   "2023-08-24 13:00:00" an absolute timestamp
        #   "13:00"               clock time, first occurrence at or after
        #                         origin_time (QUOTE it - YAML reads a bare
        #                         13:00:00 as the number 46800)
        # Both ends are INCLUSIVE. Either may be left null for "from the
        # start" / "to the end". When set, skip_times/n_times/time_end are
        # ignored.
        "time_from": None,
        "time_to": None,
        # ---- window in RECORD INDICES (simple, but per-file) -------------
        "skip_times": 0,
        # NUMBER of records to convert, counted from skip_times. The
        # default keeps everything. Up to 0.9.0 this was passed to range()
        # as the exclusive END index, so skip_times=10 with n_times=5
        # selected nothing; it is a genuine count from 0.10.0 on. Set
        # `time_end` instead for the old exclusive-end meaning.
        "n_times": 99999,
        "time_end": None,
        # Some PALM files (surface output in particular) are written with the
        # time variable left at its _FillValue, so they carry no timestamps.
        # Left as null, the output interval is read from the p3d namelist
        # (dt_dosurf_av / dt_data_output_av and the matching skip_time_*,
        # chosen per product), and the first record is placed one interval
        # in. Set these only to override that, or when there is no p3d.
        "time_step": None,
        "time_start": None,    # seconds since origin_time of the FIRST record
    },
    "voxel": {
        "enabled": False,
        "file_from": None,         # DATA_3D_NETCDF*
        "file_to": None,           # output .nc (relative -> output_dir)
        "do_velocity": True,       # destagger and write u, v, w
        "scalars": [],             # e.g. ["theta", "p"]
        # Drop the z = 0 level so vertical spacing is uniform
        # (dz/2, 3dz/2, ...) - voxel layers handle that best.
        "drop_k0": True,
        # Add `ta_derived`: air temperature from theta via a hydrostatic
        # Exner profile (T = theta * Pi, Pi integrated up the column from
        # surface_pressure). Requires theta, which is added to scalars if
        # missing. PALM's `p` output is the PERTURBATION pressure and is
        # deliberately NOT used. If PALM already writes `ta`, listing it in
        # scalars copies the model's own value, which is preferable.
        "derive_ta": False,
    },
    "surfaces": {
        "enabled": False,
        "file_from": None,         # SURFACE_DATA_NETCDF*
        "prefix": "palm_",         # output name prefix
        "variables": [],           # e.g. ["theta_surface"]
        # Topmost upward-facing element per column (ground + roofs) ->
        # one multiband GeoTIFF per variable, band = timestep.
        "write_topview_raster": True,
        # The same top-view fields as a CF NetCDF, which QGIS opens as a
        # MESH layer. MDAL reads the CF time axis natively, so the Temporal
        # Controller works with no style sidecar and no band-to-time
        # mapping - the reliable route for animation in QGIS. All
        # variables share one file, as separate dataset groups.
        "write_topview_netcdf": True,
        # Wall elements -> vertical PolygonZ quads in a GeoPackage,
        # long format with a `time` field (ArcGIS time slider ready).
        "write_walls_gpkg": True,
        # Add a ta_* column next to every theta_* variable on the element
        # layers: actual temperature at each face's own height.
        "derive_ta": False,
        # Vertical faces have a ZERO-AREA footprint, so the polygon layers
        # show nothing in a 2D map view. This adds a companion "<layer>_2d"
        # layer inside the same GeoPackage holding each wall's horizontal
        # trace as a LineString, which IS visible and clickable in plan
        # view. The polygons remain for 3D scenes.
        "wall_lines": True,
        # EVERY surface element of the static driver geometry - ground,
        # roofs and all four wall orientations - as PolygonZ in one
        # GeoPackage, with an `orient` column (up / down / wall_N / wall_E
        # / wall_S / wall_W) and `z_agl` (height above that column's
        # terrain). This is the complete surface-temperature picture; it is
        # also the largest output (rows = n_elements x n_times).
        "write_all_elements_gpkg": False,
        # Layer name inside the GeoPackage (null -> "surfaces" / "walls").
        "layer_name": None,
    },
    "wind": {
        "enabled": False,
        "file_from": None,         # DATA_3D_NETCDF* (needs u, v)
        "file_to": None,           # output .gpkg (relative -> output_dir)
        "heights_agl": [2.0],      # sampling heights above the LOCAL
                                   # surface (terrain + buildings) [m]
        "stride": 1,               # keep every n-th grid column
    },
    # Write QGIS style sidecars (<layer>.qml, plus an embedded style in
    # each GeoPackage) so layers open with their time settings already
    # configured: band-to-time mapping on the rasters, the `time` field
    # recognised on the vector layers. Harmless for ArcGIS, which ignores
    # them.
    "qgis_styles": True,
    # Write temperature variables (theta*, tsurf*/t_surf*, ta*) in degrees
    # Celsius instead of kelvin, in EVERY product: voxel NetCDF, top-view
    # GeoTIFF/mesh NetCDF, and the GeoPackage element layers. Applied at
    # write time only - all physics (the Exner/ta derivation) still runs in
    # kelvin, and NoData cells keep their fill value. The `units` attribute
    # in the NetCDF outputs is set to "degC" accordingly. Inputs must stay
    # raw kelvin either way (see the raw-input guard).
    "celsius": False,
    # Pressure at the model's z = 0 level (the domain base, i.e. at
    # origin_z - NOT sea level), used for the Exner function when deriving
    # air temperature. Left null (the default) it is read from the p3d:
    # PALM's own `surface_pressure`, which the namelist gives in hPa and
    # which is converted to Pa here. Falls back to 101325 Pa only when
    # there is no p3d and no value here - a standard-atmosphere placeholder
    # that is wrong for any domain much above sea level. This matters: T
    # scales as (p/p0)^0.286, so a 10 hPa error shifts the derived
    # temperature by about 0.85 K.
    "surface_pressure": None,
    "canopy": {
        "enabled": False,
        # Where `lad` lives - the static driver, unless overridden.
        "file_from": None,
        "prefix": "palm_",
        # LAD integrated over the column -> leaf area index [m2/m2].
        "write_lai": True,
        # Highest level with LAD above lad_threshold, in metres above the
        # LOCAL terrain (a canopy height model). Derived from LAD itself,
        # so it does not need buildings_2d.
        "write_canopy_height": True,
        "lad_threshold": 0.01,
        # LAD lifted from height-above-terrain onto absolute heights, as a
        # CF NetCDF: the ArcGIS voxel / QGIS mesh equivalent of what
        # paraview_vis produced for ParaView.
        "write_volume": False,
    },
    "stages": list(STAGES),
}


# ------------------------------
# 3. HELPERS
# ------------------------------
def _deep_merge(base, override):
    out = copy.deepcopy(base)
    for k, v in (override or {}).items():
        if isinstance(v, dict) and isinstance(out.get(k), dict):
            out[k] = _deep_merge(out[k], v)
        else:
            out[k] = copy.deepcopy(v)
    return out


class ConfigError(Exception):
    pass


# `time_from: 13:00` looks like a clock time but YAML 1.1 resolves it as a
# SEXAGESIMAL NUMBER - 780 here, 46800 for 13:00:00 - which then arrives as
# a plain seconds offset and silently converts the wrong period. There is no
# way to tell that apart from a genuine number once YAML has parsed it, so
# the raw text is checked instead.
_UNQUOTED_CLOCK = re.compile(
    r"^\s*(time_from|time_to)\s*:\s*(\d{1,2}:\d{2}(?::\d{2})?)\s*(?:#.*)?$",
    re.M)


def _check_unquoted_clock(text):
    m = _UNQUOTED_CLOCK.search(text)
    if m:
        key, val = m.group(1), m.group(2)
        raise ConfigError(
            "bounds.{key}: {val} must be QUOTED (\"{val}\"). Unquoted, YAML "
            "reads it as a sexagesimal number, not a clock time, and it "
            "would be taken as an offset in seconds from origin_time."
            .format(key=key, val=val))


# ------------------------------
# 4. CONFIG CLASS
# ------------------------------
class Config:
    """Loaded + validated palm2gis configuration.

    Access settings with item notation (cfg["voxel"]["file_to"]). Relative
    paths resolve against project.root (default: the config file's
    directory); output paths additionally resolve against
    project.output_dir.
    """

    def __init__(self, path):
        self._path = Path(path).resolve()
        if not self._path.is_file():
            raise FileNotFoundError(f"Config file not found: {path}")
        text = self._path.read_text()
        user = yaml.safe_load(text) or {}
        _check_unquoted_clock(text)

        self._cfg = _deep_merge(DEFAULTS, user)
        root = self._cfg["project"]["root"]
        self.root = Path(root).resolve() if root else self._path.parent
        self.output_dir = self.resolve(self._cfg["project"]["output_dir"])

        self._validate()

    # -- access ------------------------------------------------------------
    def __getitem__(self, key):
        return self._cfg[key]

    def as_plain_dict(self):
        return copy.deepcopy(self._cfg)

    # -- paths -------------------------------------------------------------
    def resolve(self, p):
        """Resolve a path against project.root."""
        p = Path(p)
        return p if p.is_absolute() else (self.root / p)

    def out_path(self, p):
        """Resolve an output path against project.output_dir."""
        p = Path(p)
        return p if p.is_absolute() else (self.output_dir / p)

    def disp(self, p):
        """Path as shown in logs: relative to output_dir where possible, so
        lines stay short (the output dir is logged once per run)."""
        p = Path(p)
        try:
            return str(p.relative_to(self.output_dir))
        except ValueError:
            return str(p)

    # -- namelist-derived settings -----------------------------------------
    P_STANDARD = 101325.0        # standard atmosphere at SEA LEVEL [Pa]

    def resolve_surface_pressure(self, p3d):
        """Fill in `surface_pressure` from the p3d when it was left null.

        PALM's namelist `surface_pressure` is in HECTOPASCALS (default
        1013.25); everything here works in Pa, so it is converted. A value
        set explicitly in the YAML is in Pa and wins untouched.

        Without this the config carried a hardcoded 101325 Pa - the
        standard atmosphere at SEA LEVEL. Prague sits near 200 m, where
        the real pressure is some 20 hPa lower, which biases every derived
        air temperature by about 1.7 K. Nothing downstream could detect
        that: the field simply comes out uniformly too warm.
        """
        if self._cfg["surface_pressure"] is not None:
            debug("[static] surface_pressure = {:.0f} Pa (from the config).",
                  float(self._cfg["surface_pressure"]))
            return

        raw = (p3d or {}).get("surface_pressure")
        if raw is None:
            self._cfg["surface_pressure"] = self.P_STANDARD
            info("[static] no surface_pressure in the p3d or the config - "
                 "using the standard {:.0f} Pa. That is SEA LEVEL: if the "
                 "domain base sits well above it, derived air temperatures "
                 "are biased warm by ~0.85 K per 10 hPa.", self.P_STANDARD)
            return

        raw = float(raw)
        # PALM reads this key in hPa. A four- or six-digit number means the
        # namelist was written in Pa instead - take it at face value rather
        # than turning 101325 into 10 132 500 Pa, but say so, because PALM
        # ITSELF would have read it as hPa and the run may be affected.
        if raw > 2000.0:
            warning("[static] p3d surface_pressure = {:g} looks like Pa, but "
                    "PALM reads that key in hPa (and would have used {:g} "
                    "Pa). Taking it as {:.0f} Pa here - check the namelist.",
                    raw, raw * 100.0, raw)
            p_pa = raw
        else:
            p_pa = raw * 100.0

        self._cfg["surface_pressure"] = p_pa
        info("[static] surface_pressure = {:.0f} Pa ({:g} hPa, from the "
             "p3d).", p_pa, p_pa / 100.0)

    # -- overwrite policy --------------------------------------------------
    @property
    def overwrite(self):
        return bool(self._cfg["project"]["overwrite"])

    def prepare_output(self, path, tag=""):
        """Decide whether `path` may be written, and get it ready.

        Returns False when the file already exists and project.overwrite is
        false - the caller then skips that product, so a rerun only fills in
        what is missing. With overwrite on, an existing file is REMOVED
        first: GDAL appends to an existing GeoPackage layer instead of
        replacing it, which would silently duplicate every feature.

        An existing file whose recorded config hash differs from the
        current one is rewritten even without --overwrite. Skipping on
        filename alone meant that editing bounds, celsius or a variable
        list and rerunning silently kept the old outputs.
        """
        path = Path(path)
        label = "[{}] ".format(tag) if tag else ""
        stale = False
        if path.exists():
            state = getattr(self, "_state", None)
            if state is not None and tag:
                from .state import stage_hash
                stale = state.is_stale(path, stage_hash(self, tag))
            if not self.overwrite and not stale:
                info("{}kept (exists): {}", label, self.disp(path))
                return False
            if stale and not self.overwrite:
                info("{}config changed since it was written - rewriting {}",
                     label, self.disp(path))
            path.unlink()
            debug("{}replacing {}", label, self.disp(path))
        path.parent.mkdir(parents=True, exist_ok=True)
        # Remember that this run is (re)writing the file. Only these get
        # re-stamped with provenance afterwards; stamping every existing
        # output would touch files a rerun is supposed to leave alone.
        try:
            self._written.add(Path(path).resolve())
        except AttributeError:
            self._written = {Path(path).resolve()}
        return True

    def was_written(self, path):
        """True when this run actually (re)wrote `path`."""
        return Path(path).resolve() in getattr(self, "_written", set())

    @contextmanager
    def atomic_output(self, path, tag=""):
        """Write to a temporary file, move it into place only on success.

        Yields the temporary path to write to. A stage that raises - or a
        run that is interrupted - therefore leaves no half-written output
        behind that a later rerun would mistake for a finished product
        (the skip logic goes by file existence, so a truncated NetCDF
        would be kept silently). The temporary sits in the destination
        directory so the final move stays on one filesystem and is atomic.
        """
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        tmp = path.with_name(path.name + ".partial")
        if tmp.exists():
            tmp.unlink()
        try:
            yield tmp
        except BaseException:
            if tmp.exists():
                tmp.unlink()
                debug("{}discarded the partial {}",
                      "[{}] ".format(tag) if tag else "", tmp.name)
            raise
        os.replace(tmp, path)

    def note_written(self, path, tag=""):
        """Record which config produced `path` (see palm2gis.state)."""
        state = getattr(self, "_state", None)
        if state is not None and tag:
            from .state import stage_hash
            state.record(path, stage_hash(self, tag))

    # -- validation --------------------------------------------------------
    def _validate(self):
        c = self._cfg
        if not c["domain"]["static_driver"]:
            raise ConfigError("domain.static_driver is required.")
        if not c["domain"]["p3d"] and (c["domain"]["nz"] is None
                                       or c["domain"]["dz"] is None):
            raise ConfigError("Give either domain.p3d or both domain.nz "
                              "and domain.dz.")

        # Catch a malformed time window at load time rather than halfway
        # through the first stage. A dummy origin is enough to check the
        # form; the real origin_time is only needed for the conversion.
        import datetime as _dt
        from .steps.georef import parse_time_bound
        ref = _dt.datetime(2000, 1, 1)
        try:
            _, dt0 = parse_time_bound(c["bounds"].get("time_from"), ref,
                                      "bounds.time_from")
            parse_time_bound(c["bounds"].get("time_to"), ref,
                             "bounds.time_to", after=dt0)
        except RuntimeError as exc:
            raise ConfigError(str(exc))

        for stage in c["stages"]:
            if stage not in STAGES:
                raise ConfigError(f"Unknown stage in stages: {stage} "
                                  f"(known: {', '.join(STAGES)}).")

        if c["voxel"]["enabled"]:
            for key in ("file_from", "file_to"):
                if not c["voxel"][key]:
                    raise ConfigError(f"voxel.{key} is required when "
                                      "voxel.enabled: true.")
            if not c["voxel"]["do_velocity"] and not c["voxel"]["scalars"]:
                raise ConfigError("voxel: nothing to write - enable "
                                  "do_velocity or list scalars.")
        if c["surfaces"]["enabled"]:
            if not c["surfaces"]["file_from"]:
                raise ConfigError("surfaces.file_from is required when "
                                  "surfaces.enabled: true.")
            if not c["surfaces"]["variables"]:
                raise ConfigError("surfaces.variables must list at least "
                                  "one variable.")
        if c["wind"]["enabled"]:
            for key in ("file_from", "file_to"):
                if not c["wind"][key]:
                    raise ConfigError(f"wind.{key} is required when "
                                      "wind.enabled: true.")
            if not c["wind"]["heights_agl"]:
                raise ConfigError("wind.heights_agl must list at least one "
                                  "height.")
