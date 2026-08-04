"""QGIS style sidecars, so the layers open ready to animate.

QGIS applies a `<name>.qml` file sitting next to a data source when the
layer is added, and temporal settings are part of what a QML style
carries. This module writes those sidecars, so a multiband GeoTIFF arrives
with its band-to-time mapping already filled in and a GeoPackage arrives
with its `time` field already recognised - no manual Temporal tab, no
expression typing.

The XML follows QgsRasterLayerTemporalProperties::writeXml and
QgsVectorLayerTemporalProperties::writeXml:

  raster  <temporal enabled mode="3" bandNumber fetchMode>
              <ranges><range band begin end .../></ranges>
          mode 3 = FixedRangePerBand, added in QGIS 3.38 - older QGIS
          ignores the block and simply opens the raster untimed.

  vector  <temporal enabled mode="1" startField="time"
                    durationUnit fixedDuration accumulate limitMode>
          mode 1 = FeatureDateTimeInstantFromField.

Rasters also get a sensible single-band pseudocolour ramp so the first
view is not a grey blur; GeoPackages additionally get the style stored in
the file's own `layer_styles` table, which QGIS reads in preference to a
sidecar.
"""

import datetime
import sqlite3
from pathlib import Path
from xml.sax.saxutils import escape

from ..log import debug, warning

QGIS_VERSION = "3.38.0"          # minimum that understands FixedRangePerBand


# ------------------------------
# 1. RASTER: TIME RANGE PER BAND
# ------------------------------
def _iso(dt):
    """ISO timestamp with an explicit UTC marker.

    PALM's origin_time is UTC, so every timestamp palm2gis derives is UTC.
    QGIS writes per-band ranges WITHOUT a timezone (unlike its other
    temporal modes, which use an explicit offset), and a naive string is
    read back as LOCAL time while the map's temporal range is UTC - which
    shifts every band by the local offset and can make the layer render
    nothing at all. The trailing Z removes the ambiguity.
    """
    return dt.strftime("%Y-%m-%dT%H:%M:%SZ")


def raster_qml(band_times, default_step_s=3600.0):
    """QML for a multiband raster whose bands are timesteps.

    band_times : list of datetime (or None) in band order.
    Each band covers from its own time to the next one; the last band gets
    the same span as the previous gap so it stays visible at the end of an
    animation.

    Only temporal settings are written. An incomplete renderer here can
    leave the layer blank, and QGIS's own default (min/max stretch on the
    data) is a better starting point than anything guessed in advance.
    """
    known = [(i, t) for i, t in enumerate(band_times) if t is not None]
    if not known:
        return None

    ranges = []
    for n, (i, t) in enumerate(known):
        if n + 1 < len(known):
            end = known[n + 1][1]
        else:
            span = ((known[n][1] - known[n - 1][1]) if n > 0
                    else datetime.timedelta(seconds=default_step_s))
            end = t + span
        ranges.append(
            '      <range band="{}" begin="{}" end="{}" '
            'includeBeginning="1" includeEnd="0"/>'.format(
                i + 1, _iso(t), _iso(end)))

    return (
        '<!DOCTYPE qgis PUBLIC \'http://mrcc.com/qgis.dtd\' \'SYSTEM\'>\n'
        '<qgis version="{v}">\n'
        '  <temporal enabled="1" mode="3" bandNumber="1" fetchMode="0">\n'
        '    <ranges>\n{r}\n    </ranges>\n'
        '  </temporal>\n'
        '</qgis>\n'
    ).format(v=QGIS_VERSION, r="\n".join(ranges))


# ------------------------------
# 2. VECTOR: TIME FROM A FIELD
# ------------------------------
def vector_qml(field="time", duration_s=3600.0):
    """QML for a layer whose features carry an instant in `field`.

    Each feature stays visible for `duration_s` after its timestamp, so an
    hourly layer is continuous rather than flashing on exact instants.
    """
    if duration_s % 3600 == 0:
        unit, value = "h", duration_s / 3600.0
    elif duration_s % 60 == 0:
        unit, value = "min", duration_s / 60.0
    else:
        unit, value = "s", duration_s
    return (
        '<!DOCTYPE qgis PUBLIC \'http://mrcc.com/qgis.dtd\' \'SYSTEM\'>\n'
        '<qgis version="{v}">\n'
        '  <temporal enabled="1" mode="1" limitMode="0" '
        'startField="{f}" endField="" startExpression="" endExpression="" '
        'durationField="" durationUnit="{u}" fixedDuration="{d}" '
        'accumulate="0">\n'
        '    <fixedRange><start></start><end></end></fixedRange>\n'
        '  </temporal>\n'
        '</qgis>\n'
    ).format(v=QGIS_VERSION, f=escape(field), u=unit, d=value)


# ------------------------------
# 3. WRITING SIDECARS
# ------------------------------
def write_sidecar(path, qml):
    """Write <path>.qml next to the data file (QGIS picks it up on add)."""
    if not qml:
        return None
    out = Path(path).with_suffix(".qml")
    out.write_text(qml)
    debug("[qgis] wrote {}", out.name)
    return out


_STYLE_TABLE = """
CREATE TABLE IF NOT EXISTS layer_styles (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    f_table_catalog TEXT,
    f_table_schema TEXT,
    f_table_name TEXT,
    f_geometry_column TEXT,
    styleName TEXT,
    styleQML TEXT,
    styleSLD TEXT,
    useAsDefault BOOLEAN,
    description TEXT,
    owner TEXT,
    ui TEXT,
    update_time DATETIME DEFAULT (strftime('%Y-%m-%dT%H:%M:%fZ','now'))
)
"""


def embed_gpkg_style(path, layer, qml, style_name="palm2gis"):
    """Store the style inside the GeoPackage's own layer_styles table.

    QGIS reads that table when the layer is added and prefers it over a
    sidecar, so the styled layer travels inside the single .gpkg file.
    """
    if not qml:
        return
    try:
        con = sqlite3.connect(str(path))
        cur = con.cursor()
        cur.execute(_STYLE_TABLE)
        cur.execute("DELETE FROM layer_styles WHERE f_table_name = ? "
                    "AND styleName = ?", (layer, style_name))
        cur.execute(
            "INSERT INTO layer_styles (f_table_catalog, f_table_schema, "
            "f_table_name, f_geometry_column, styleName, styleQML, "
            "styleSLD, useAsDefault, description, owner, ui) "
            "VALUES ('', '', ?, 'geom', ?, ?, '', 1, "
            "'Generated by palm2gis', '', NULL)",
            (layer, style_name, qml))
        con.commit()
        con.close()
        debug("[qgis] embedded style in {}", Path(path).name)
    except Exception as exc:                      # never fail the run on this
        warning("[qgis] could not embed the style in {} ({}); the .qml "
                "sidecar still applies.", Path(path).name, exc)
