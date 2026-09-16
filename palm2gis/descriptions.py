"""Human-readable descriptions of every variable and attribute field.

One table, used by every output:

  * NetCDF (voxel, wind grid, canopy) - written as `long_name`, which
    ArcGIS Pro shows as the variable's Description when a voxel or
    multidimensional layer is added.
  * GeoPackage (surfaces, wind) - written to `gpkg_data_columns`: `name`
    is the field alias, `description` the longer text.
  * File Geodatabase (surfaces) - written as the field alias.

A shapefile has no place for either, so the multipatch output carries
none.
"""

import sqlite3
from pathlib import Path

import pyogrio

from .log import debug, warning


# ------------------------------
# 1. DESCRIPTION TABLE
# ------------------------------
# name -> (alias, description). The alias is short enough for a table
# header; the description says what the quantity is and where it comes
# from. Units are deliberately absent: they are carried by the `units`
# attribute, and a temperature's units depend on `celsius`.
DESCRIPTIONS = {
    # -- 3D atmosphere ------------------------------------------------------
    "u": ("Eastward wind",
          "eastward wind component (u), interpolated to cell centres"),
    "v": ("Northward wind",
          "northward wind component (v), interpolated to cell centres"),
    "w": ("Upward wind",
          "upward wind component (w), interpolated to cell centres"),
    "theta": ("Potential temperature", "potential temperature"),
    "pt": ("Potential temperature", "potential temperature"),
    "thetal": ("Liquid water pot. temperature",
               "liquid water potential temperature"),
    "thetav": ("Virtual pot. temperature", "virtual potential temperature"),
    "ta": ("Air temperature", "air temperature"),
    "q": ("Total water content",
          "total water content (mixing ratio); equals specific humidity "
          "when no cloud water is present"),
    "qv": ("Water vapour content", "water vapour content (mixing ratio)"),
    "ql": ("Liquid water content", "liquid water content (mixing ratio)"),
    "qc": ("Cloud water content", "cloud water content (mixing ratio)"),
    "qr": ("Rain water content", "rain water content (mixing ratio)"),
    "ah": ("Absolute humidity", "absolute humidity (water vapour density)"),
    "rh": ("Relative humidity", "relative humidity"),
    "p": ("Perturbation pressure",
          "perturbation pressure (deviation from the hydrostatic "
          "reference)"),
    "e": ("Subgrid TKE", "subgrid-scale turbulent kinetic energy"),
    "s": ("Passive scalar", "passive scalar concentration"),
    "wspeed": ("Wind speed", "horizontal wind speed"),
    "wdir": ("Wind direction", "wind direction"),
    "lad": ("Leaf area density", "leaf area density of the plant canopy"),
    "bad": ("Basal area density",
            "basal area density (branches and trunks)"),
    "bio_mrt": ("Mean radiant temperature", "mean radiant temperature"),
    "bio_pet": ("PET", "physiological equivalent temperature"),
    "bio_utci": ("UTCI", "universal thermal climate index"),
    "bio_perct": ("Perceived temperature", "perceived temperature"),

    # -- wind samples (wind stage) -----------------------------------------
    "speed": ("Wind speed", "horizontal wind speed, sqrt(u^2 + v^2)"),
    "dir_met": ("Wind direction",
                "meteorological wind direction: where the wind blows FROM, "
                "clockwise from north"),
    "height_agl": ("Height above surface",
                   "sampling height above the local surface (terrain and "
                   "buildings)"),

    # -- surface output (surfaces stage) -----------------------------------
    "t_surf": ("Surface temperature", "surface (skin) temperature"),
    "theta_surface": ("Surface pot. temperature",
                      "surface potential temperature"),
    "theta_surf": ("Surface pot. temperature",
                   "surface potential temperature"),
    "q_surface": ("Surface humidity", "specific humidity at the surface"),
    "theta1": ("Pot. temperature (1st level)",
               "potential temperature at the first grid level above the "
               "surface"),
    "qv1": ("Water vapour (1st level)",
            "water vapour content at the first grid level above the "
            "surface"),
    "uvw1": ("Wind speed (1st level)",
             "wind speed at the first grid level above the surface"),
    "us": ("Friction velocity", "friction velocity"),
    "ts": ("Temperature scale", "Monin-Obukhov temperature scale"),
    "qs": ("Humidity scale", "Monin-Obukhov humidity scale"),
    "ol": ("Obukhov length", "Obukhov length"),
    "z0": ("Roughness length", "roughness length for momentum"),
    "z0h": ("Roughness length (heat)", "roughness length for heat"),
    "z0q": ("Roughness length (moisture)", "roughness length for moisture"),
    "shf": ("Sensible heat flux", "surface sensible heat flux"),
    "qsws": ("Latent heat flux", "surface latent heat flux"),
    "ghf": ("Ground heat flux", "ground (storage) heat flux"),
    "usws": ("u momentum flux", "surface momentum flux, u component"),
    "vsws": ("v momentum flux", "surface momentum flux, v component"),
    "rad_net": ("Net radiation", "net radiation flux at the surface"),
    "rad_sw_in": ("Incoming shortwave", "incoming shortwave radiation flux"),
    "rad_sw_out": ("Outgoing shortwave",
                   "outgoing (reflected) shortwave radiation flux"),
    "rad_lw_in": ("Incoming longwave", "incoming longwave radiation flux"),
    "rad_lw_out": ("Outgoing longwave",
                   "outgoing (emitted and reflected) longwave radiation "
                   "flux"),
    "waste_heat": ("Waste heat", "anthropogenic waste heat flux"),
    "r_a": ("Aerodynamic resistance", "aerodynamic resistance"),
    "r_s": ("Surface resistance", "total surface resistance"),

    # -- geometry / bookkeeping fields -------------------------------------
    "elem_id": ("Element index", "PALM surface element index"),
    "orient": ("Orientation",
               "face orientation: up, down, wall_N, wall_E, wall_S, wall_W"),
    "zenith": ("Zenith angle",
               "zenith angle of the face normal (0 = up, 90 = wall, "
               "180 = down)"),
    "azimuth": ("Azimuth",
                "azimuth of the wall normal, clockwise from north (empty "
                "for horizontal faces)"),
    "z": ("Height above sea level",
          "height above sea level of the element (PALM height + "
          "origin_z)"),
    "z_agl": ("Height above terrain",
              "height above the local terrain"),
    "rec": ("Record index", "PALM output record index in the source file"),
    "time": ("Time", "time of the output record (UTC)"),
}

# Families recognised by prefix when the exact name is not listed.
_PREFIXES = (
    ("kc_", "Concentration: {rest}", "concentration of {rest}"),
    ("em_", "Emission: {rest}", "emission of {rest}"),
)


# ------------------------------
# 2. LOOKUP
# ------------------------------
def describe(name, source_long_name=None):
    """(alias, description) for *name*.

    Unlisted names fall back to a prefix family, then to the source
    file's own `long_name` when it says more than the name itself (PALM
    usually writes the bare variable name there), and finally to the
    name.
    """
    key = str(name).rstrip("*")
    if key in DESCRIPTIONS:
        return DESCRIPTIONS[key]
    for prefix, alias, text in _PREFIXES:
        if key.startswith(prefix) and len(key) > len(prefix):
            rest = key[len(prefix):]
            return alias.format(rest=rest), text.format(rest=rest)
    src = str(source_long_name or "").strip()
    if src and src.lower() != key.lower():
        return src, src
    return key, key


def long_name(name, source_long_name=None):
    """NetCDF `long_name` (ArcGIS Pro's Description) for *name*."""
    return describe(name, source_long_name)[1]


def arrow_field_metadata(name):
    """Arrow field metadata GDAL turns into a field alias and comment."""
    alias, text = describe(name)
    return {b"GDAL:OGR:alternative_name": alias.encode("utf-8"),
            b"GDAL:OGR:comment": text.encode("utf-8")}


# ------------------------------
# 3. GEOPACKAGE FIELD ALIASES
# ------------------------------
_DATA_COLUMNS = (
    "CREATE TABLE IF NOT EXISTS gpkg_data_columns ("
    "table_name TEXT NOT NULL, column_name TEXT NOT NULL, name TEXT, "
    "title TEXT, description TEXT, mime_type TEXT, constraint_name TEXT, "
    "CONSTRAINT pk_gdc PRIMARY KEY (table_name, column_name), "
    "CONSTRAINT gdc_tn UNIQUE (table_name, name))")
_EXTENSION = ("gpkg_data_columns", None, "gpkg_schema",
              "http://www.geopackage.org/spec121/#extension_schema",
              "read-write")


def write_gpkg_descriptions(path, layer):
    """Alias and description for every attribute field of *layer*.

    Written to the GeoPackage schema extension table, exactly as GDAL
    stores an OGR alternative name and comment. Never fails the run.
    """
    try:
        con = sqlite3.connect(str(path))
        cur = con.cursor()
        cols = [r[1] for r in cur.execute(
            "PRAGMA table_info(\"{}\")".format(layer.replace('"', '""')))]
        geom = {r[0] for r in cur.execute(
            "SELECT column_name FROM gpkg_geometry_columns "
            "WHERE table_name = ?", (layer,))}
        cur.execute(_DATA_COLUMNS)
        cur.execute("CREATE TABLE IF NOT EXISTS gpkg_extensions ("
                    "table_name TEXT, column_name TEXT, "
                    "extension_name TEXT NOT NULL, definition TEXT NOT NULL, "
                    "scope TEXT NOT NULL, CONSTRAINT ge_tce UNIQUE "
                    "(table_name, column_name, extension_name))")
        if not cur.execute("SELECT 1 FROM gpkg_extensions WHERE table_name "
                           "= ? AND extension_name = ?",
                           _EXTENSION[:1] + _EXTENSION[2:3]).fetchone():
            cur.execute("INSERT INTO gpkg_extensions VALUES (?, ?, ?, ?, ?)",
                        _EXTENSION)
        cur.execute("DELETE FROM gpkg_data_columns WHERE table_name = ?",
                    (layer,))
        used = set()
        for col in cols:
            if col == "fid" or col in geom:
                continue
            alias, text = describe(col)
            if alias in used:                 # (table_name, name) is UNIQUE
                alias = "{} ({})".format(alias, col)
            used.add(alias)
            cur.execute("INSERT INTO gpkg_data_columns (table_name, "
                        "column_name, name, title, description) "
                        "VALUES (?, ?, ?, ?, ?)",
                        (layer, col, alias, alias, text))
        con.commit()
        con.close()
        debug("[fields] descriptions written for {} layer '{}'",
              Path(path).name, layer)
    except Exception as exc:
        warning("[fields] could not write field descriptions to {} ({}).",
                Path(path).name, exc)


# ------------------------------
# 4. FILE GEODATABASE FIELD ALIASES
# ------------------------------
def write_with_aliases(frame, path, driver, layer, append=False):
    """Write *frame* with pyogrio, carrying field aliases where possible.

    pyogrio's DataFrame writer has no way to pass an alias, but GDAL
    reads one from Arrow field metadata. The layer is therefore CREATED
    through the Arrow writer (the first frame), after which later frames
    are appended the usual way and inherit the field definitions. Without
    pyarrow the frame is written plainly, with no aliases.
    """
    if append:
        pyogrio.write_dataframe(frame, path, driver=driver, layer=layer,
                                append=True)
        return
    try:
        import pyarrow as pa
        import pyarrow.compute as pc
        ok = (hasattr(pyogrio, "write_arrow")
              and tuple(pyogrio.__gdal_version__) >= (3, 8, 0))
    except ImportError:
        ok = False
    if not ok:
        debug("[fields] pyarrow, pyogrio.write_arrow or GDAL >= 3.8 missing "
              "- {} written without field aliases.", Path(path).name)
        pyogrio.write_dataframe(frame, path, driver=driver, layer=layer)
        return

    geom = frame.geometry.name
    tbl = pa.table(frame.to_arrow(index=False, geometry_encoding="WKB"))
    cols, fields = [], []
    for field, col in zip(tbl.schema, tbl.columns):
        if field.name != geom:
            if pa.types.is_null(field.type):
                # an all-empty column (e.g. time when the file has none):
                # GDAL cannot create a field of Arrow's null type, and the
                # DataFrame writer would have made it a string field
                field = pa.field(field.name, pa.string(), True,
                                 field.metadata)
                col = col.cast(pa.string())
            elif pa.types.is_floating(field.type):
                # NaN -> NULL, as write_dataframe(nan_as_null=True) does
                col = pc.if_else(pc.is_nan(col), pa.scalar(None, field.type),
                                 col)
            meta = dict(field.metadata or {})
            meta.update(arrow_field_metadata(field.name))
            field = field.with_metadata(meta)
        cols.append(col)
        fields.append(field)
    tbl = pa.Table.from_arrays(
        cols, schema=pa.schema(fields, metadata=tbl.schema.metadata))

    kinds = frame.geom_type.dropna().unique()
    gtype = str(kinds[0]) if len(kinds) == 1 else "Unknown"
    if gtype != "Unknown" and bool(frame.has_z.any()):
        gtype += " Z"
    pyogrio.write_arrow(tbl, path, driver=driver, layer=layer,
                        geometry_name=geom, geometry_type=gtype,
                        crs=frame.crs.to_wkt() if frame.crs else None)
