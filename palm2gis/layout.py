"""The newcomer config layout, translated onto the internal settings.

A config in this layout says what to READ once (`input`), what part of the
run to convert (`region`), and lists every output as its own block with a
`file:` line - a filename produces it, null or a missing block skips it.
Rarely needed settings live in `advanced`. See template.yaml.

Everything downstream still works on the internal structure (`domain`,
`bounds`, `voxel.file_from`, ...), which is also what the legacy layout
writes directly. This module only maps one onto the other, and checks the
new layout strictly: an unknown key is an error with a suggestion, because
a silently ignored typo is the commonest way a config does not do what it
says.
"""

import difflib
import re

# ------------------------------
# 1. SCHEMA
# ------------------------------
# Known keys per block. A dict value is a nested block with its own keys.
_OUT_FILE = {"file": None}

SCHEMA = {
    "project": {"name": None, "root": None, "output_dir": None,
                "output_crs": None, "overwrite": None},
    "input": {"static_driver": None, "p3d": None, "dynamic_driver": None,
              "data_3d": None, "surface_data": None, "crs": None},
    "region": {"x": None, "y": None, "height_above_surface": None,
               # pre-0.25 names, still accepted:
               "height_above_highest_surface": None,
               "time_from": None, "time_to": None},
    "time": {"from": None, "to": None},
    "celsius": None,             # pre-0.25; advanced.units now
    "voxel": {"enabled": None, "variables": None, "netcdf": _OUT_FILE},
    "surfaces": {"enabled": None, "variables": None,
                 "all_elements_gpkg": _OUT_FILE,
                 "all_elements_fgdb": _OUT_FILE,
                 "all_elements_multipatch": _OUT_FILE},
    "wind": {"enabled": None, "heights_above_surface": None,
             "netcdf": {"file": None, "heights": None,
                        "every_nth_column": None},
             "gpkg": {"file": None, "every_nth_column": None,
                      "arrows_3d": None, "arrow_length_per_ms": None,
                      "tilt_with_w": None}},
    "canopy": {"enabled": None, "lad_netcdf": _OUT_FILE},
    "advanced": {
        # Subsections (0.25 onwards).
        "units": {"temperature": None, "potential_temperature": None,
                  "overrides": None},
        "grid": {"nz": None, "dz": None, "origin_time": None,
                 "origin_z": None},
        "basic_state": {"surface_pressure": None, "pt_surface": None},
        "netcdf": {"crs": None, "voxel_drop_k0": None, "resampling": None},
        "qgis": {"styles": None, "styles_3d": None},
        "shapefile": {"multipatch_split": None},
        "records": {"skip": None, "count": None, "end": None},
        "time_axis": {"step": None, "start": None},
        "canopy_file_from": None, "stages": None,
        # Flat pre-0.25 names, still accepted.
        "nz": None, "dz": None, "origin_time": None, "origin_z": None,
        "surface_pressure": None, "pt_surface": None,
        "temperature_units": None, "qgis_styles": None,
        "qgis_styles_3d": None, "voxel_drop_k0": None,
        "multipatch_split": None, "netcdf_crs": None,
        "skip_times": None, "n_times": None, "time_end": None,
        "time_step": None, "time_start": None,
    },
}

# advanced: subsection key -> (internal name, flat pre-0.25 name).
_ADV_GROUPS = {
    "units": {"temperature": ("units_temperature", None),
              "potential_temperature": ("units_potential_temperature", None),
              "overrides": ("temperature_units", "temperature_units")},
    "grid": {"nz": ("nz", "nz"), "dz": ("dz", "dz"),
             "origin_time": ("origin_time", "origin_time"),
             "origin_z": ("origin_z", "origin_z")},
    "basic_state": {"surface_pressure": ("surface_pressure",
                                         "surface_pressure"),
                    "pt_surface": ("pt_surface", "pt_surface")},
    "netcdf": {"crs": ("netcdf_crs", "netcdf_crs"),
               "voxel_drop_k0": ("voxel_drop_k0", "voxel_drop_k0"),
               "resampling": ("resampling", None)},
    "qgis": {"styles": ("qgis_styles", "qgis_styles"),
             "styles_3d": ("qgis_styles_3d", "qgis_styles_3d")},
    "shapefile": {"multipatch_split": ("multipatch_split",
                                       "multipatch_split")},
    "records": {"skip": ("skip_times", "skip_times"),
                "count": ("n_times", "n_times"),
                "end": ("time_end", "time_end")},
    "time_axis": {"step": ("time_step", "time_step"),
                  "start": ("time_start", "time_start")},
}

# Top-level keys of the legacy layout, and where each went.
_LEGACY_HINT = {
    "domain": "input (static_driver, p3d) and advanced.grid (nz, dz, "
              "origin_*)",
    "bounds": "region (x, y, height_above_surface), time (from, to) and "
              "advanced.records",
    "crs": "input.crs",
    "stages": "advanced.stages",
    "surface_pressure": "advanced.basic_state.surface_pressure",
    "pt_surface": "advanced.basic_state.pt_surface",
    "temperature_units": "advanced.units.overrides",
    "qgis_styles": "advanced.qgis.styles",
    "qgis_styles_3d": "advanced.qgis.styles_3d",
}

VELOCITY = ("u", "v", "w")


class LayoutError(Exception):
    pass


# ------------------------------
# 2. DETECTION AND KEY CHECKING
# ------------------------------
def is_new_layout(user):
    """The new layout is recognised by its own top-level blocks."""
    return isinstance(user, dict) and any(
        k in user for k in ("input", "region", "time", "advanced"))


def _check_keys(block, schema, where):
    if block is None:
        return
    if not isinstance(block, dict):
        if schema is _OUT_FILE:
            raise LayoutError(f"{where} must be a block with a file line:\n"
                              f"    {where.rsplit('.', 1)[-1]}:\n"
                              f"        file: {block!r}")
        raise LayoutError(f"{where} must be a block of settings, got "
                          f"{block!r}.")
    for key, val in block.items():
        path = f"{where}.{key}" if where else key
        if key not in schema:
            if not where and key in _LEGACY_HINT:
                raise LayoutError(
                    f"'{key}' belongs to the old config layout. In this "
                    f"layout it is {_LEGACY_HINT[key]}.")
            close = difflib.get_close_matches(key, list(schema), n=1)
            hint = f" Did you mean '{close[0]}'?" if close else ""
            raise LayoutError(f"Unknown setting '{path}'.{hint} Known "
                              f"here: {', '.join(schema)}.")
        if isinstance(schema[key], dict):
            _check_keys(val, schema[key], path)


# ------------------------------
# 3. HELPERS
# ------------------------------
def _bool(val, where):
    """true/false, or None when unset. YAML reads "false" in quotes as a
    non-empty string, which Python would take as True."""
    if val is None or isinstance(val, bool):
        return val
    raise LayoutError(f"{where} must be true or false (without quotes), "
                      f"got {val!r}.")


def _names(val, where):
    """A list of variable names, without duplicates."""
    if val is None:
        return []
    if not isinstance(val, (list, tuple)):
        raise LayoutError(f"{where} must be a list, e.g. [{val}], got "
                          f"{val!r}.")
    names = [str(v) for v in val]
    dup = sorted({v for v in names if names.count(v) > 1})
    if dup:
        raise LayoutError(f"{where} lists {', '.join(dup)} more than once.")
    return names


def _number(val, where, minimum=None, integer=False):
    if val is None:
        return None
    ok = (isinstance(val, int) if integer else isinstance(val, (int, float)))
    if not ok or isinstance(val, bool):
        kind = "a whole number" if integer else "a number"
        raise LayoutError(f"{where} must be {kind}, got {val!r}.")
    if minimum is not None and val < minimum:
        raise LayoutError(f"{where} must be at least {minimum}, got {val}.")
    return val


# Extension each output writer needs. A GeoPackage written as `.nc`
# opens nowhere, and GDAL refuses a File Geodatabase without `.gdb` only
# AFTER the other surface layers have been written.
_EXTENSIONS = {"netcdf": ".nc", "lad_netcdf": ".nc", "gpkg": ".gpkg",
               "all_elements_gpkg": ".gpkg", "all_elements_fgdb": ".gdb",
               "all_elements_multipatch": ".shp"}


def _set(out, section, key, value):
    """Copy a setting only when the user gave it, so defaults stay put."""
    if value is not None:
        out.setdefault(section, {})[key] = value


def _file(block, name, where):
    """Filename of an output block, with {name} expanded, or None."""
    if block is not None and not isinstance(block, dict):
        raise LayoutError(f"{where} must be a block with a 'file:' line, "
                          f"got {block!r}.")
    f = (block or {}).get("file")
    if f is None:
        return None
    if not isinstance(f, str) or not f.strip():
        raise LayoutError(f"{where}.file must be a filename or null, "
                          f"got {f!r}.")
    f = f.replace("{name}", name)
    # The wind netCDF is one file per height: {z} is filled in per file.
    rest = f.replace("{z}", "") if where == "wind.netcdf" else f
    if "{" in rest or "}" in rest:
        allowed = ("{name} and {z}" if where == "wind.netcdf"
                   else "{name}")
        raise LayoutError(f"{where}.file: only {allowed} can be filled in, "
                          f"got {f!r}.")
    ext = _EXTENSIONS[where.rsplit(".", 1)[-1]]
    if not f.lower().endswith(ext):
        raise LayoutError(f"{where}.file must end in {ext}, got {f!r}.")
    return f


def _span(region, axis):
    val = region.get(axis)
    if val is None:
        return None, None
    if (not isinstance(val, (list, tuple)) or len(val) != 2
            or not all(isinstance(v, (int, float)) for v in val)):
        raise LayoutError(f"region.{axis} must be [from, to] in metres, "
                          f"e.g. [100, 300]; got {val!r}.")
    if val[1] <= val[0]:
        raise LayoutError(f"region.{axis}: 'to' must be larger than 'from' "
                          f"({val!r}).")
    return val[0], val[1]


def _advanced(user, notes):
    """advanced: -> {internal name: value}, reading the subsections and the
    flat pre-0.25 keys alike. A flat key is reported once and the
    subsection wins if both are given."""
    adv = user.get("advanced") or {}
    if not isinstance(adv, dict):
        raise LayoutError(f"advanced must be a block of settings, got "
                          f"{adv!r}.")
    out, where = {}, {}
    for group, keys in _ADV_GROUPS.items():
        for key, (internal, flat) in keys.items():
            if flat is not None and adv.get(flat) is not None:
                out[internal] = adv[flat]
                where[internal] = f"advanced.{flat}"
                notes.append(f"advanced.{flat} is now "
                             f"advanced.{group}.{key} - the old name still "
                             f"works for now.")
            val = (adv.get(group) or {}).get(key)
            if val is not None:
                out[internal] = val
                where[internal] = f"advanced.{group}.{key}"
    for key in ("canopy_file_from", "stages"):
        if adv.get(key) is not None:
            out[key] = adv[key]
            where[key] = f"advanced.{key}"
    return out, where


def _units(user, adv, where, notes):
    """advanced.units -> the internal `celsius` flag."""
    if user.get("celsius") is not None:
        notes.append("celsius is now advanced.units.temperature (C or K) "
                     "- the old name still works for now.")
        return _bool(user["celsius"], "celsius")
    scale = adv.get("units_temperature")
    pot = adv.get("units_potential_temperature")
    if pot is not None and str(pot).upper() != "K":
        raise LayoutError(
            "advanced.units.potential_temperature must be K: theta and "
            "theta_surface are kelvin-defined, so a Celsius value would "
            "not be a potential temperature. List ta or t_surface in the "
            "stage's variables for an actual temperature in degrees C.")
    if scale is None:
        return None
    if str(scale).upper() not in ("K", "C"):
        raise LayoutError(f"{where['units_temperature']} must be K or C, "
                          f"got {scale!r}.")
    return str(scale).upper() == "C"


def _time_window(user, reg, notes):
    """The time section, or the pre-0.25 region.time_from / time_to."""
    blk = user.get("time") or {}
    if not isinstance(blk, dict):
        raise LayoutError(f"time must be a block with from: and to: lines, "
                          f"got {blk!r}.")
    out = {"time_from": blk.get("from"), "time_to": blk.get("to")}
    for old, new in (("time_from", "from"), ("time_to", "to")):
        if reg.get(old) is not None:
            notes.append(f"region.{old} is now time.{new} - the old name "
                         f"still works for now.")
            if out[old] is None:
                out[old] = reg[old]
    return out


def _height(reg, notes):
    """The vertical cut, metres above the surface below."""
    val = reg.get("height_above_surface")
    old = reg.get("height_above_highest_surface")
    if old is not None:
        notes.append("region.height_above_highest_surface is now "
                     "region.height_above_surface - the old name still "
                     "works for now.")
        if val is None:
            val = old
    return _number(val, "region.height_above_surface", 0)


def _enabled(block, where):
    """A stage runs when its block exists and enabled is not false."""
    if block is None:
        return False
    if not isinstance(block, dict):
        raise LayoutError(f"{where} must be a block of settings, got "
                          f"{block!r}.")
    return _bool(block.get("enabled"), f"{where}.enabled") is not False


# ------------------------------
# 4. TRANSLATION
# ------------------------------
def translate(user):
    """New-layout dict -> internal settings dict (still merged with the
    defaults afterwards). Also returns notes for the log."""
    _check_keys(user, SCHEMA, "")
    notes = []
    out = {}

    project = dict(user.get("project") or {})
    _bool(project.get("overwrite"), "project.overwrite")
    name = str(project.get("name") or "palm2gis")
    out["project"] = {k: v for k, v in project.items()
                      if v is not None and k != "output_crs"}

    inp = user.get("input") or {}
    reg = user.get("region") or {}
    adv, adv_where = _advanced(user, notes)

    # -- input -------------------------------------------------------------
    if not inp.get("static_driver"):
        raise LayoutError("input.static_driver is required.")
    _set(out, "domain", "static_driver", inp.get("static_driver"))
    _set(out, "domain", "p3d", inp.get("p3d"))
    _set(out, "domain", "dynamic_driver", inp.get("dynamic_driver"))
    _set(out, "crs", "palm", inp.get("crs"))
    _set(out, "crs", "output", project.get("output_crs"))
    _set(out, "crs", "resampling", adv.get("resampling"))
    _set(out, "crs", "netcdf_attributes", adv.get("netcdf_crs"))
    _number(adv.get("nz"), adv_where.get("nz", "advanced.grid.nz"), 1,
            integer=True)
    _number(adv.get("dz"), adv_where.get("dz", "advanced.grid.dz"), 0)
    for key in ("nz", "dz", "origin_time", "origin_z"):
        _set(out, "domain", key, adv.get(key))
    if not inp.get("p3d") and (adv.get("nz") is None
                               or adv.get("dz") is None):
        raise LayoutError("Give input.p3d, or both advanced.grid.nz and "
                          "advanced.grid.dz.")

    # -- region and time -----------------------------------------------------
    x0, xmax = _span(reg, "x")
    y0, ymax = _span(reg, "y")
    window = _time_window(user, reg, notes)
    for key, val in (("x0", x0), ("xmax", xmax), ("y0", y0),
                     ("ymax", ymax),
                     ("zmax", _height(reg, notes)),
                     ("time_from", window["time_from"]),
                     ("time_to", window["time_to"])):
        _set(out, "bounds", key, val)
    for key in ("skip_times", "n_times", "time_end", "time_step",
                "time_start"):
        _set(out, "bounds", key, adv.get(key))

    # -- global options ------------------------------------------------------
    celsius = _units(user, adv, adv_where, notes)
    if celsius is not None:
        out["celsius"] = celsius
    for key in ("qgis_styles", "qgis_styles_3d"):
        _bool(adv.get(key), adv_where.get(key, f"advanced.qgis.{key}"))
    _number(adv.get("surface_pressure"),
            adv_where.get("surface_pressure",
                          "advanced.basic_state.surface_pressure"), 1)
    _number(adv.get("pt_surface"),
            adv_where.get("pt_surface",
                          "advanced.basic_state.pt_surface"), 1)
    for key in ("surface_pressure", "pt_surface", "temperature_units",
                "qgis_styles",
                "qgis_styles_3d"):
        if adv.get(key) is not None:
            out[key] = adv[key]
    if adv.get("stages") is not None:
        out["stages"] = _names(adv["stages"], "advanced.stages")

    data_3d = inp.get("data_3d")
    prefix = f"{name}_"

    # -- voxel ---------------------------------------------------------------
    vx = user.get("voxel")
    vx_on = _enabled(vx, "voxel")
    vx_file = _file((vx or {}).get("netcdf"), name, "voxel.netcdf")
    if vx_on and vx_file:
        variables = _names(vx.get("variables"), "voxel.variables")
        if not variables:
            raise LayoutError("voxel.variables must list at least one "
                              "variable, e.g. [u, v, w, theta].")
        vel = [v for v in VELOCITY if v in variables]
        if vel and len(vel) != 3:
            raise LayoutError("voxel.variables: list u, v and w together "
                              f"(got only {', '.join(vel)}).")
        if not data_3d:
            raise LayoutError("voxel needs input.data_3d.")
        out["voxel"] = {
            "enabled": True, "file_from": data_3d, "file_to": vx_file,
            "do_velocity": bool(vel), "derive_ta": "ta" in variables,
            "derive_ah": "ah" in variables,
            "derive_speed_3d": "speed_3d" in variables,
            "derive_speed": "speed" in variables,
            "scalars": [v for v in variables
                        if v not in VELOCITY
                        and v not in ("ta", "ah", "speed", "speed_3d")],
        }
        _set(out, "voxel", "drop_k0",
             _bool(adv.get("voxel_drop_k0"),
                   adv_where.get("voxel_drop_k0",
                                 "advanced.netcdf.voxel_drop_k0")))
    else:
        out["voxel"] = {"enabled": False}
        if vx_on:
            notes.append("voxel: no netcdf.file set - stage skipped.")

    # -- surfaces ------------------------------------------------------------
    sf = user.get("surfaces")
    sf_on = _enabled(sf, "surfaces")
    sf_files = {k: _file((sf or {}).get(k), name, f"surfaces.{k}")
                for k in ("all_elements_gpkg", "all_elements_fgdb",
                          "all_elements_multipatch")}
    if sf_on and any(sf_files.values()):
        variables = _names(sf.get("variables"), "surfaces.variables")
        if not variables:
            raise LayoutError("surfaces.variables must list at least one "
                              "variable, e.g. [theta_surface, t_surface].")
        if not inp.get("surface_data"):
            raise LayoutError("surfaces needs input.surface_data.")
        if "t_surf" in variables:
            notes.append("surfaces.variables: 't_surf' is now 't_surface' "
                         "- rename it in the config.")
            variables = ["t_surface" if v == "t_surf" else v
                         for v in variables]
            variables = list(dict.fromkeys(variables))
        others = [v for v in variables if v != "t_surface"]
        # t_surface next to another variable means "ensure it": PALM's own
        # if the file has it, derived from the potential temperature if
        # not. Listed alone it can only be PALM's own.
        derive = "t_surface" in variables and bool(others)
        out["surfaces"] = {
            "enabled": True, "file_from": inp["surface_data"],
            "prefix": prefix,
            "variables": others if derive else variables,
            "derive_t_surface": derive,
            "write_all_elements_gpkg": bool(sf_files["all_elements_gpkg"]),
            "all_elements_to": sf_files["all_elements_gpkg"],
            "write_fgdb": bool(sf_files["all_elements_fgdb"]),
            "fgdb_to": sf_files["all_elements_fgdb"],
            "write_multipatch": bool(sf_files["all_elements_multipatch"]),
            "multipatch_to": sf_files["all_elements_multipatch"],
        }
        _set(out, "surfaces", "multipatch_split",
             _number(adv.get("multipatch_split"),
                     adv_where.get("multipatch_split",
                                   "advanced.shapefile.multipatch_split"),
                     1, integer=True))
    else:
        out["surfaces"] = {"enabled": False}
        if sf_on:
            notes.append("surfaces: no output file set - stage skipped.")

    # -- wind ----------------------------------------------------------------
    wd = user.get("wind")
    wd_on = _enabled(wd, "wind")
    wd_nc = _file((wd or {}).get("netcdf"), name, "wind.netcdf")
    wd_gp = _file((wd or {}).get("gpkg"), name, "wind.gpkg")
    nc_blk = (wd or {}).get("netcdf") or {}
    gp_blk = (wd or {}).get("gpkg") or {}
    if wd_on and (wd_nc or wd_gp):
        heights = wd.get("heights_above_surface")
        if wd_gp:
            if not isinstance(heights, (list, tuple)) or not heights:
                raise LayoutError("wind.heights_above_surface must list at "
                                  "least one height for wind.gpkg, e.g. "
                                  f"[2.0]; got {heights!r}.")
            for h in heights:
                _number(h, "wind.heights_above_surface", 0)
            if len(set(heights)) != len(heights):
                raise LayoutError("wind.heights_above_surface lists a "
                                  "height more than once.")
        grid_heights = nc_blk.get("heights")
        if wd_nc and grid_heights is not None:
            # Left out: every PALM level of the voxel's vertical range.
            if not isinstance(grid_heights, (list, tuple)) or not grid_heights:
                raise LayoutError(
                    "wind.netcdf.heights must list heights ABOVE SEA LEVEL "
                    "[m], e.g. [200, 210], or be left out for one file per "
                    f"level of the voxel range; got {grid_heights!r}.")
            for h in grid_heights:
                if h is None:
                    raise LayoutError("wind.netcdf.heights holds an empty "
                                      "entry.")
                _number(h, "wind.netcdf.heights", None)
            if len(set(grid_heights)) != len(grid_heights):
                raise LayoutError("wind.netcdf.heights lists a height more "
                                  "than once.")
        _number(nc_blk.get("every_nth_column"),
                "wind.netcdf.every_nth_column", 1, integer=True)
        _number(gp_blk.get("every_nth_column"),
                "wind.gpkg.every_nth_column", 1, integer=True)
        _number(gp_blk.get("arrow_length_per_ms"),
                "wind.gpkg.arrow_length_per_ms", 0)
        _bool(gp_blk.get("arrows_3d"), "wind.gpkg.arrows_3d")
        _bool(gp_blk.get("tilt_with_w"), "wind.gpkg.tilt_with_w")
        if not data_3d:
            raise LayoutError("wind needs input.data_3d.")
        out["wind"] = {"enabled": True, "file_from": data_3d,
                       "grid_to": wd_nc, "file_to": wd_gp}
        if wd_gp:
            out["wind"]["heights_agl"] = [float(h) for h in heights]
        if wd_nc and grid_heights is not None:
            out["wind"]["grid_heights"] = [float(h) for h in grid_heights]
        _set(out, "wind", "grid_stride", nc_blk.get("every_nth_column"))
        _set(out, "wind", "stride", gp_blk.get("every_nth_column"))
        _set(out, "wind", "write_lines_3d", gp_blk.get("arrows_3d"))
        _set(out, "wind", "arrow_scale", gp_blk.get("arrow_length_per_ms"))
        _set(out, "wind", "include_w", gp_blk.get("tilt_with_w"))
    else:
        out["wind"] = {"enabled": False}
        if wd_on:
            notes.append("wind: no output file set - stage skipped.")

    # -- canopy --------------------------------------------------------------
    cn = user.get("canopy")
    cn_on = _enabled(cn, "canopy")
    cn_file = _file((cn or {}).get("lad_netcdf"), name, "canopy.lad_netcdf")
    if cn_on and cn_file:
        out["canopy"] = {"enabled": True, "prefix": prefix,
                         "write_volume": True, "volume_to": cn_file}
        _set(out, "canopy", "file_from", adv.get("canopy_file_from"))
    else:
        out["canopy"] = {"enabled": False}
        if cn_on:
            notes.append("canopy: no lad_netcdf.file set - stage skipped.")

    return out, notes


# ------------------------------
# 5. MESSAGES IN THE USER'S VOCABULARY
# ------------------------------
# Stage code and its messages speak the internal names. For a config in
# the new layout, the log rewrites them into the names that config uses,
# so "raise bounds.zmax" reads "raise region.height_above_highest_surface".
MESSAGE_NAMES = {
    "bounds.time_from": "time.from",
    "bounds.time_to": "time.to",
    "bounds.zmax": "region.height_above_surface",
    "bounds.x0": "region.x", "bounds.xmax": "region.x",
    "bounds.y0": "region.y", "bounds.ymax": "region.y",
    "bounds.skip_times": "advanced.records.skip",
    "bounds.n_times": "advanced.records.count",
    "bounds.time_end": "advanced.records.end",
    "bounds.time_step": "advanced.time_axis.step",
    "bounds.time_start": "advanced.time_axis.start",
    "bounds cut": "region window",
    "domain.static_driver": "input.static_driver",
    "domain.p3d": "input.p3d",
    "domain.dynamic_driver": "input.dynamic_driver",
    "domain.nz": "advanced.grid.nz", "domain.dz": "advanced.grid.dz",
    "domain.origin_time": "advanced.grid.origin_time",
    "domain.origin_z": "advanced.grid.origin_z",
    "surface_pressure": "advanced.basic_state.surface_pressure",
    "pt_surface": "advanced.basic_state.pt_surface",
    "crs.palm": "input.crs",
    "crs.output": "project.output_crs",
    "crs.resampling": "advanced.netcdf.resampling",
    "voxel.drop_k0": "advanced.netcdf.voxel_drop_k0",
    "Set voxel.do_velocity: false to write scalars only":
        "Remove u, v and w from voxel.variables to write scalars only",
    "voxel.scalars": "voxel.variables",
    "derive_ta is on": "ta is listed",
    "derive_t_surface is on": "t_surface is listed",
    "derive_ah is on": "ah is listed",
    "derive_speed_3d is on": "speed_3d is listed",
    "voxel.file_from": "input.data_3d",
    "wind.file_from": "input.data_3d",
    "surfaces.file_from": "input.surface_data",
    "canopy.file_from": "advanced.canopy_file_from",
    "point file_from at": "point the input file at",
    "Point file_from at": "Point the input file at",
    "wind.stride": "wind.gpkg.every_nth_column",
    "wind.grid_stride": "wind.netcdf.every_nth_column",
    "wind.grid_heights": "wind.netcdf.heights",
    "wind.heights_agl": "wind.heights_above_surface",
}
_MESSAGE_RE = re.compile(
    r"(?<![\w.])(" + "|".join(re.escape(k) for k in
                              sorted(MESSAGE_NAMES, key=len, reverse=True))
    + r")(?!\w)")


def to_new_names(text):
    """Rewrite internal setting names in a message into the new layout."""
    return _MESSAGE_RE.sub(lambda m: MESSAGE_NAMES[m.group(1)], text)
