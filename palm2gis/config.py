"""YAML config loading + defaults for palm2gis."""

import copy
import os
import re
from contextlib import contextmanager
from pathlib import Path

import yaml

from .layout import LayoutError, is_new_layout, translate
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
        # Override the static driver's origin_z attribute [m above sea
        # level]. PALM's z = 0 sits at origin_z, so this is what turns
        # relative heights into absolute ones. Static drivers often carry
        # 0 because nobody set it; the run warns when that happens.
        "origin_z": None,
        # PALM dynamic driver, read for the basic-state surface pressure
        # (`surface_forcing_surface_pressure`, first value). Null = look
        # for <job>_dynamic next to the p3d, for the ROOT domain only: a
        # nested child runs with the namelist value instead.
        "dynamic_driver": None,
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
        # Ensure `ta` - AIR temperature. Not the surface skin temperature,
        # which is `t_surface` and comes from the surfaces stage; `ta` and
        # `t_surface` are different quantities and neither is ever named
        # ta_surf/ta_surface. Derived from theta with PALM's basic-state
        # Exner function (T = theta * Pi, Pi from surface_pressure and
        # pt_surface, see below). Requires theta, which is added to scalars
        # if missing. A derived value is tagged derived_from = "theta" so it
        # cannot be confused with PALM's own. PALM's `p` output is the
        # PERTURBATION pressure and is deliberately NOT used. If PALM
        # already writes `ta`, it is copied - the model's own value is
        # exact.
        "derive_ta": False,
        # Ensure `ah` - absolute humidity (water vapour density, g/m3).
        # PALM's own `ah` is copied if the file has one; otherwise it is
        # derived from `q` as q * p / (Rd * Tv), with T from PALM's `ta`
        # (or from theta, as for derive_ta) and p from the same basic
        # state. Requires q.
        "derive_ah": False,
        # Write `speed_3d` - the full 3D wind speed sqrt(u^2 + v^2 + w^2)
        # from the components destaggered to cell centres. Needs u, v and w
        # in the file; they are written only if do_velocity is also on.
        # Not the wind stage's `speed`, which is horizontal.
        "derive_speed_3d": False,
    },
    "surfaces": {
        "enabled": False,
        "file_from": None,         # SURFACE_DATA_NETCDF*
        "prefix": "palm_",         # output name prefix
        "variables": [],           # e.g. ["theta_surface"]
        # ESRI MultiPatch shapefile: ArcGIS's native 3D geometry type,
        # and what a 3D Object scene layer is published from. One file,
        # all element orientations (use `orient` to filter), all
        # timesteps. DBF limits apply - see the notes in surfaces.py.
        "multipatch_to": None,
        "write_multipatch": True,
        # Split into one shapefile per timestep once elements x times
        # exceeds this. A .shp and its .dbf are each capped at 2 GB.
        "multipatch_split": 2_000_000,
        # File Geodatabase holding the same elements as PolygonZ. The
        # GeoPackage is for QGIS; this is what ArcGIS opens natively,
        # and unlike the shapefile it keeps real datetimes, full-length
        # field names and NULLs.
        "fgdb_to": None,
        "write_fgdb": True,
        # Add `t_surface` - the actual SURFACE (skin) temperature - next
        # to every surface potential temperature (`theta_surface`),
        # converted at each face's own height. If the file already holds
        # it (`t_surface` or `t_surf`) its own value is used.
        #
        # NOT to be confused with voxel.derive_ta: `ta` is PALM's TRUE AIR
        # temperature, an atmospheric quantity, and nothing in this stage
        # produces one. This option was called `derive_ta` up to 0.17.0,
        # which was wrong twice over, and `derive_t_surf` up to 0.21 (the
        # output was then `t_surf`); both old names still work and warn.
        "derive_t_surface": False,
        # EVERY surface element of the static driver geometry - ground,
        # roofs and all four wall orientations - as PolygonZ in one
        # GeoPackage, with an `orient` column (up / down / wall_N / wall_E
        # / wall_S / wall_W) and `z_agl` (height above that column's
        # terrain). This is the complete surface-temperature picture; it is
        # also the largest output (rows = n_elements x n_times).
        "write_all_elements_gpkg": False,
        "all_elements_to": None,   # null -> <prefix>surfaces.gpkg
        # Layer name inside the GeoPackage (null -> "surfaces").
        "layer_name": None,
    },
    "wind": {
        "enabled": False,
        "file_from": None,         # DATA_3D_NETCDF* (needs u, v)
        "file_to": None,           # output .gpkg (relative -> output_dir)
        "heights_agl": [2.0],      # sampling heights above the LOCAL
                                   # surface (terrain + buildings) [m]
        "stride": 1,               # keep every n-th grid column
        # Gridded NetCDF for ArcGIS Pro's flow / vector-field renderers:
        # (time, height, y, x) with u, v, w, speed and dir_met. Written
        # alongside the GeoPackage, not instead of it.
        "grid_to": None,           # output .nc (relative -> output_dir)
        # Thinning for grid_to only. Default 1 (full resolution): the
        # renderer thins its own streamlines, and a 200x200x7x2 grid of
        # five float32 variables is about 11 MB.
        "grid_stride": 1,
        # A second layer, `wind_3d`, of LineStringZ arrows at the TRUE
        # sample height (local surface + height_agl). The point layer holds
        # 2D geometry with the height only as an attribute, so in a 3D
        # scene every level collapses onto the ground and the 2 m and 10 m
        # samples land on top of each other. Both layers are written: a
        # line layer reads badly from directly above, a flat point layer
        # reads badly in a scene.
        "write_lines_3d": False,
        # Arrow length in metres per m/s. Has to suit the domain: a 1 m
        # Holesovice cut and a 10 m Prague cut want very different values.
        "arrow_scale": 2.0,
        # Tilt the arrows with the vertical velocity. Off by default: it
        # shows updraughts in street canyons but is harder to read against
        # buildings. `speed` and `dir_met` stay horizontal either way, so
        # the two layers remain comparable; with this on, `w` is also
        # written as a field.
        "include_w": False,
    },
    # Write QGIS style sidecars (<layer>.qml, plus an embedded style in
    # each GeoPackage) so layers open with their time settings already
    # configured: the `time` field recognised on the vector layers. Harmless for ArcGIS, which ignores
    # them.
    "qgis_styles": True,
    # Add a 3D renderer block to the vector styles, so wall polygons and
    # wind arrows keep their true height in a QGIS 3D scene instead of
    # being draped on the terrain. The 3D symbol schema has moved between
    # QGIS releases more than the rest of the QML: if a layer refuses to
    # open, set this to false and the 2D style and temporal settings are
    # unaffected.
    "qgis_styles_3d": True,
    # Write temperature variables (theta*, tsurf*/t_surf*, t_surface, ta*)
    # in degrees Celsius instead of kelvin, in EVERY product: voxel NetCDF
    # and the surface element layers. Applied at
    # write time only - all physics (the Exner/ta derivation) still runs in
    # kelvin, and NoData cells keep their fill value. The `units` attribute
    # in the NetCDF outputs is set to "degC" accordingly. Inputs must stay
    # raw kelvin either way (see the raw-input guard).
    "celsius": False,
    # Explicit temperature scale for variables whose `units` attribute
    # cannot be read. Maps a variable name to "K" or "C", e.g.
    # {ta_2m*_xy: C}. Needed rarely: PALM's own unit strings are
    # recognised, including the truncated "degree_" it writes from a
    # fixed-length buffer. A stage RAISES on an unrecognised unit rather
    # than guessing, so this is the escape hatch when a new PALM version
    # writes something new.
    "temperature_units": {},
    # PALM's basic state, which fixes the Exner function used to derive ta,
    # t_surface and ah from potential temperature:
    #   Pi(z) = (surface_pressure / p0)^kappa - g z / (cp pt_surface)
    # Both are what PALM itself ran with, NOT the real atmosphere - PALM
    # converted T <-> theta with them, so only the same values invert that
    # exactly. Left null (the default) they are resolved per domain: see
    # georef.resolve_basic_state (dynamic driver for the root domain, the
    # p3d, the parent p3d for a nested child), and cross-checked against
    # ta / theta whenever a 3D file holds both. A value set here wins.
    # surface_pressure is in Pa (the namelist gives hPa), pt_surface in K.
    "surface_pressure": None,
    "pt_surface": None,
    "canopy": {
        "enabled": False,
        # Where `lad` lives - the static driver, unless overridden.
        "file_from": None,
        "prefix": "palm_",
        # LAD lifted from height-above-terrain onto absolute heights, as a
        # CF NetCDF: the ArcGIS voxel / QGIS mesh equivalent of what
        # paraview_vis produced for ParaView. The stage's only output.
        "write_volume": True,
        "volume_to": None,         # null -> <prefix>lad_3d.nc
    },
    "stages": list(STAGES),
}


# ------------------------------
# 3. HELPERS
# ------------------------------
# Keys renamed since a released version: {section: {old: new}}. Honoured
# on load with a warning, so an existing job config keeps working.
_RENAMED_KEYS = {
    "surfaces": {"derive_ta": "derive_t_surface",
                 "derive_t_surf": "derive_t_surface"},
}


def _migrate_renamed_keys(user):
    """Accept pre-0.18 key names, warning once per key.

    surfaces.derive_ta derived the SURFACE skin temperature - while `ta`
    in PALM is the true AIR temperature. The name claimed the wrong
    quantity and the wrong medium, and sat next to voxel.derive_ta, which
    does genuinely derive air temperature. surfaces.derive_t_surf became
    derive_t_surface in 0.22 with the output renamed to `t_surface`, to
    match `theta_surface`. A legacy `t_surf` in surfaces.variables is
    renamed too.
    """
    sf = user.get("surfaces")
    if isinstance(sf, dict) and isinstance(sf.get("variables"), list) \
            and "t_surf" in sf["variables"]:
        warning("config: surfaces.variables lists 't_surf', which is now "
                "written as 't_surface' (matching theta_surface). The old "
                "name still works for now.")
        sf["variables"] = ["t_surface" if v == "t_surf" else v
                           for v in sf["variables"]]
    for section, renames in _RENAMED_KEYS.items():
        block = user.get(section)
        if not isinstance(block, dict):
            continue
        for old, new in renames.items():
            if old not in block:
                continue
            if new in block:
                warning("config: {}.{} is deprecated and {}.{} is also "
                        "set - using {}.{}.",
                        section, old, section, new, section, new)
            else:
                warning("config: {}.{} has been renamed to {}.{}. It "
                        "derives t_surface, a SURFACE temperature; `ta` is "
                        "PALM's true AIR temperature. The old name still "
                        "works for now.", section, old, section, new)
                block[new] = block[old]
            block.pop(old)


# Keys whose output was removed: {section: (key, ...)}. Ignored on load with
# a warning, so an existing job config still runs instead of silently
# producing less than it asks for.
_REMOVED_KEYS = {
    "surfaces": ("write_topview_raster", "write_topview_netcdf",
                 "write_walls_gpkg", "wall_lines"),
    "canopy": ("write_lai", "write_canopy_height", "lad_threshold"),
}


def _drop_removed_keys(user):
    """Warn about and discard keys for outputs removed in 0.19.0."""
    for section, keys in _REMOVED_KEYS.items():
        block = user.get(section)
        if not isinstance(block, dict):
            continue
        for key in keys:
            if key not in block:
                continue
            if block.pop(key):
                warning("config: {}.{} was removed in 0.19.0 and is "
                        "ignored - that output is no longer produced.",
                        section, key)


def _deep_merge(base, override):
    out = copy.deepcopy(base)
    for k, v in (override or {}).items():
        if isinstance(v, dict) and isinstance(out.get(k), dict):
            out[k] = _deep_merge(out[k], v)
        elif v is None and isinstance(out.get(k), dict):
            # A section header whose every option is commented out parses as
            # None - `bounds:` with nothing under it. That means "leave this
            # section at its defaults", not "replace it with null", and
            # taking it literally made everything downstream fail on
            # NoneType. Nulling a whole section has no useful meaning, so
            # the defaults are kept.
            continue
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
            "{key}: {val} must be QUOTED (\"{val}\"). Unquoted, YAML "
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

        # Two layouts load: the newcomer layout (input / region / one
        # block per output, see palm2gis.layout) is translated onto the
        # internal settings; the legacy layout already IS those settings.
        if is_new_layout(user):
            try:
                user, notes = translate(user)
            except LayoutError as exc:
                raise ConfigError(str(exc))
            from .layout import to_new_names
            from .log import set_message_names
            set_message_names(to_new_names)
            for note in notes:
                info("config: {}", note)
            self.layout = "new"
        else:
            from .log import set_message_names
            set_message_names(None)
            _migrate_renamed_keys(user)
            _drop_removed_keys(user)
            self.layout = "legacy"
        self._cfg = _deep_merge(DEFAULTS, user)
        self._explicit_basic_state = {
            k: self._cfg[k] is not None
            for k in ("surface_pressure", "pt_surface")}
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

    # -- run-derived settings ----------------------------------------------
    def basic_state_explicit(self):
        """Which basic-state values the user set in the config (as loaded,
        before resolve_basic_state filled the rest in)."""
        return dict(self._explicit_basic_state)

    def set_basic_state(self, surface_pressure, pt_surface):
        """Store the resolved basic state [Pa, K]; see
        georef.resolve_basic_state."""
        self._cfg["surface_pressure"] = float(surface_pressure)
        self._cfg["pt_surface"] = float(pt_surface)

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
            # A File Geodatabase is a DIRECTORY, not a file, so unlink
            # raises on it. Everything else here is a plain file.
            if path.is_dir():
                import shutil
                shutil.rmtree(path, ignore_errors=True)
            else:
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
            if not (c["voxel"]["do_velocity"] or c["voxel"]["scalars"]
                    or c["voxel"]["derive_ta"] or c["voxel"]["derive_ah"]
                    or c["voxel"]["derive_speed_3d"]):
                raise ConfigError("voxel: nothing to write - enable "
                                  "do_velocity, derive_ta, derive_ah or "
                                  "derive_speed_3d, "
                                  "or list scalars.")
        if c["surfaces"]["enabled"]:
            if not c["surfaces"]["file_from"]:
                raise ConfigError("surfaces.file_from is required when "
                                  "surfaces.enabled: true.")
            if not c["surfaces"]["variables"]:
                raise ConfigError("surfaces.variables must list at least "
                                  "one variable.")
        if c["wind"]["enabled"]:
            if not c["wind"]["file_from"]:
                raise ConfigError("wind.file_from is required when "
                                  "wind.enabled: true.")
            if not (c["wind"]["file_to"] or c["wind"]["grid_to"]):
                raise ConfigError("wind: set file_to (GeoPackage), grid_to "
                                  "(NetCDF) or both when wind.enabled: "
                                  "true.")
            if (not isinstance(c["wind"]["heights_agl"], (list, tuple))
                    or not c["wind"]["heights_agl"]):
                raise ConfigError("wind.heights_agl must list at least one "
                                  "height.")
            for key in ("stride", "grid_stride"):
                v = c["wind"][key]
                if not isinstance(v, int) or isinstance(v, bool) or v < 1:
                    raise ConfigError(f"wind.{key} must be a whole number "
                                      f"of at least 1, got {v!r}.")
        for key, val in (("project.overwrite", c["project"]["overwrite"]),
                         ("celsius", c["celsius"])):
            if not isinstance(val, bool):
                raise ConfigError(f"{key} must be true or false (without "
                                  f"quotes), got {val!r}.")

        # Two outputs on one path: the later stage silently replaces the
        # earlier one's file (a wind grid named like the voxel NetCDF left
        # no voxel output at all, and the run still reported success).
        from .pipeline import stage_outputs
        seen = {}
        for stage in STAGES:
            if not c[stage]["enabled"]:
                continue
            for path in stage_outputs(self, stage):
                key = str(Path(path).resolve())
                if key in seen:
                    raise ConfigError(
                        "{} and {} would both write {} - give them "
                        "different filenames.".format(seen[key], stage,
                                                      Path(path).name))
                seen[key] = stage
