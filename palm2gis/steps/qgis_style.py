"""QGIS style sidecars, so the layers open ready to animate.

QGIS applies a `<name>.qml` file sitting next to a data source when the
layer is added, and temporal settings are part of what a QML style
carries. This module writes those sidecars, so a GeoPackage arrives with
its `time` field already recognised - no manual Temporal tab, no
expression typing.

The XML follows QgsVectorLayerTemporalProperties::writeXml:

  vector  <temporal enabled mode="1" startField="time"
                    durationUnit fixedDuration accumulate limitMode>
          mode 1 = FeatureDateTimeInstantFromField.

GeoPackages additionally get the style stored in the file's own
`layer_styles` table, which QGIS reads in preference to a sidecar.
"""

import sqlite3
from pathlib import Path
from xml.sax.saxutils import escape

from ..log import debug, warning

QGIS_VERSION = "3.38.0"


# ------------------------------
# 1. VECTOR: TIME FROM A FIELD
# ------------------------------
# ------------------------------
# 1b. COLOUR RAMP + 3D RENDERER
# ------------------------------
# One ramp definition, referenced by BOTH the 2D renderer and the 3D
# symbol, so a layer looks the same in plan view and in a scene and
# recolouring means editing one place.
RAMP_STOPS = ("0.25;43,131,186,255:"
              "0.5;171,221,164,255:"
              "0.75;253,174,97,255")
RAMP_C1 = "44,123,182,255"
RAMP_C2 = "215,25,28,255"


def _ramp_xml(name="palm2gis"):
    return (
        '    <colorramp name="{n}" type="gradient">\n'
        '      <Option type="Map">\n'
        '        <Option name="color1" type="QString" value="{c1}"/>\n'
        '        <Option name="color2" type="QString" value="{c2}"/>\n'
        '        <Option name="discrete" type="QString" value="0"/>\n'
        '        <Option name="rampType" type="QString" value="gradient"/>\n'
        '        <Option name="stops" type="QString" value="{s}"/>\n'
        '      </Option>\n'
        '    </colorramp>\n'
    ).format(n=name, c1=RAMP_C1, c2=RAMP_C2, s=RAMP_STOPS)


def _graduated_renderer(attr, vmin, vmax, n=5, geom="polygon"):
    """Graduated colour renderer over `attr`, equal intervals."""
    span = (vmax - vmin) or 1.0
    edges = [vmin + span * i / n for i in range(n + 1)]
    ranges, symbols = [], []
    for i in range(n):
        ranges.append(
            '      <range lower="{lo:.6f}" upper="{up:.6f}" symbol="{i}" '
            'label="{lo:.3g} - {up:.3g}" render="true"/>\n'.format(
                lo=edges[i], up=edges[i + 1], i=i))
        f = i / max(n - 1, 1)
        rgb = [int(round(a + (b - a) * f)) for a, b in
               zip([int(x) for x in RAMP_C1.split(",")],
                   [int(x) for x in RAMP_C2.split(",")])]
        col = ",".join(str(c) for c in rgb)
        symbols.append(
            '      <symbol name="{i}" type="{g}" alpha="1">\n'
            '        <layer class="{cls}" enabled="1" pass="0" locked="0">\n'
            '          <Option type="Map">\n'
            '            <Option name="color" type="QString" value="{c}"/>\n'
            '            <Option name="outline_style" type="QString" '
            'value="no"/>\n'
            '          </Option>\n'
            '        </layer>\n'
            '      </symbol>\n'.format(
                i=i, g=geom, c=col,
                cls="SimpleFill" if geom == "polygon" else "SimpleLine"))
    return (
        '  <renderer-v2 type="graduatedSymbol" attr="{a}" '
        'graduatedMethod="GraduatedColor" forceraster="0" '
        'enableorderby="0" symbollevels="0">\n'
        '    <ranges>\n{r}    </ranges>\n'
        '    <symbols>\n{s}    </symbols>\n'
        '{ramp}'
        '  </renderer-v2>\n'
    ).format(a=escape(attr), r="".join(ranges), s="".join(symbols),
             ramp=_ramp_xml())


def _renderer_3d(geom="polygon", attr=None, vmin=0.0, vmax=1.0):
    """A <layer3DRenderer> block clamping geometry to ABSOLUTE height.

    The wall polygons and the wind arrows already carry true Z, so the
    3D symbol only has to stop QGIS draping them on the terrain
    (alt-clamping="absolute", alt-binding="vertex"). Without this a scene
    flattens every wall onto the ground.

    NOTE: the 3D symbol schema has changed across QGIS releases more than
    the rest of the QML. If a layer refuses to open, set qgis_styles_3d
    to false - the 2D style and the temporal settings are unaffected.
    """
    data = ('        <data alt-clamping="absolute" alt-binding="vertex" '
            'height="0" extrusion-height="0" add-back-faces="0" '
            'invert-normals="0"/>\n')
    if geom == "line":
        data = ('        <data alt-clamping="absolute" alt-binding="vertex" '
                'height="0" extrusion-height="0" width="0.4" '
                'render-as-simple-lines="0"/>\n')
    return (
        '  <layer3DRenderer type="vector" layer="" '
        'max-screen-error="3" tiling-max-level="3" show-bounding-boxes="0">\n'
        '    <symbol type="{g}">\n'
        '{d}'
        '        <material shininess="0" ambient="25,25,25,255" '
        'diffuse="{c}" specular="10,10,10,255"/>\n'
        '    </symbol>\n'
        '  </layer3DRenderer>\n'
    ).format(g=geom, d=data, c=RAMP_C1)


def vector_qml(field="time", duration_s=3600.0, geom="polygon",
               attr=None, vmin=None, vmax=None, clamp_absolute=False):
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
    body = ''
    if attr is not None and vmin is not None and vmax is not None:
        body += _graduated_renderer(attr, float(vmin), float(vmax),
                                    geom=geom)
    if clamp_absolute:
        body += _renderer_3d(geom=geom, attr=attr)
    return (
        '<!DOCTYPE qgis PUBLIC \'http://mrcc.com/qgis.dtd\' \'SYSTEM\'>\n'
        '<qgis version="{v}">\n'
        '  <temporal enabled="1" mode="1" limitMode="0" '
        'startField="{f}" endField="" startExpression="" endExpression="" '
        'durationField="" durationUnit="{u}" fixedDuration="{d}" '
        'accumulate="0">\n'
        '    <fixedRange><start></start><end></end></fixedRange>\n'
        '  </temporal>\n'
        '{b}'
        '</qgis>\n'
    ).format(v=QGIS_VERSION, f=escape(field), u=unit, d=value, b=body)


# ------------------------------
# 2. WRITING SIDECARS
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
