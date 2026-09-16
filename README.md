# palm2gis

**PALM outputs → ArcGIS Pro / QGIS layers** — companion to `palm_preproc` and `palm_postproc`.

From a PALM run's static driver and NetCDF outputs, `palm2gis` produces georeferenced GIS layers: a CF-compliant NetCDF that ArcGIS Pro opens directly as a 3D Voxel Layer (destaggered velocities and scalars, terrain as NoData), surface temperatures on every ground, roof and wall element as time-enabled 3D polygons, and pedestrian-level wind vectors sampled at heights above the local surface. Every layer is written in absolute projected coordinates taken from the static driver's `origin_x`/`origin_y`, so it lands in the right place on the map with no manual shifting. It began as an adaptation of `paraview_vis`, a PALM visualisation script by [Martin Bureš](https://www.cs.cas.cz/staff/buresm/en) (Institute of Computer Science, Czech Academy of Sciences), and keeps its configuration concepts (domain, subdomain cut, timestep window) while targeting GIS software rather than ParaView. See [Origins and attribution](#origins-and-attribution).

---

## Table of contents

- [How it works (the three stages)](#how-it-works-the-three-stages)
- [Installing and running](#installing-and-running)
- [Files you get back, and how to load them](#files-you-get-back-and-how-to-load-them)
- [Configuration](#configuration)
- [Stage details](#stage-details)
  - [voxel — the 3D volume in ArcGIS Pro](#voxel--the-3d-volume-in-arcgis-pro)
  - [surfaces — roofs, ground, and walls](#surfaces--roofs-ground-and-walls)
  - [wind — vectors at heights above ground](#wind--vectors-at-heights-above-ground)
- [Command-line reference](#command-line-reference)
- [Things worth knowing](#things-worth-knowing)
- [Origins and attribution](#origins-and-attribution)
- [Authors](#authors)

---

## How it works (the three stages)

Every run first reads the **static driver** (and optionally the `_p3d` namelist for `nz`/`dz`): it provides the grid geometry, the terrain + building elevation, the domain origin and `origin_time`, and converts the metre-based `region` section into grid index ranges shared by all stages. Then the enabled stages run in order:

| # | Stage | What it produces |
|---|-------|------------------|
| 1 | `voxel` | CF-compliant NetCDF for the ArcGIS Pro Voxel Layer — destaggered `u`/`v`/`w` plus chosen scalars on `(time, z, y, x)`, terrain and building cells as NoData |
| 2 | `surfaces` | Every ground, roof and wall element as a time-enabled PolygonZ quad: GeoPackage (QGIS), MultiPatch shapefile and File Geodatabase (ArcGIS) |
| 4 | `canopy` | Plant canopy from the static driver: the LAD volume lifted onto absolute heights, as a voxel NetCDF |
| 3 | `wind` | A GeoPackage point layer with `u`, `v`, `speed`, `dir_met` per timestep, sampled at configurable heights above the local surface |

Each stage is enabled/disabled in its own config section; the order (and the default for `--stages`) is the `stages:` list.

---

## Installing and running

```bash
# 1. Dependencies
pip install -r requirements.txt

# 2. Check that the environment works. This builds a small synthetic dataset
#    in a temporary directory, runs the whole pipeline over it and verifies
#    the outputs — no real data needed, a few seconds.
python test.py

# 3. Make a config from the annotated template, then edit it
#    NOTE: config/ holds project configs. If you ever re-extract this
#    package over an existing checkout, extract to a scratch directory and
#    copy files across, so your own configs are not replaced.
cp template.yaml config/my_project.yaml
#    → set at least project.root, the input files, and a `file:` for
#      each output you want

# 4. Run
python run_palm2gis.py -c config/my_project.yaml
```

If `test.py` fails, the message says whether a dependency is missing or a stage broke — fix that before pointing the pipeline at real data.

**A note on PROJ.** If GDAL and its `proj.db` come from different installations — the usual result of mixing conda and pip packages in one environment — GDAL cannot encode a CRS, and rasters get written with no georeferencing while everything else looks fine. It announces itself as:

```
PROJ: internal_proj_create_from_database: .../proj.db contains
DATABASE.LAYOUT.VERSION.MINOR = 4 whereas a number >= 6 is expected.
```

palm2gis detects this at startup and repoints GDAL at a matching database (rasterio's bundled copy, or pyproj's) for the duration of the run, so the outputs are correct anyway. It is a workaround and says so. The real fix is to install `gdal`, `rasterio` and `pyproj` from a **single** source — all conda-forge, or all pip — and to clear any stale `PROJ_DATA` / `PROJ_LIB` variable:

```bash
conda list | grep -E "gdal|rasterio|pyproj|proj"   # look for mixed channels ("pypi")
conda install -c conda-forge --force-reinstall rasterio pyproj libgdal proj
echo $PROJ_DATA $PROJ_LIB                          # should be empty or match
```

Useful variations:

```bash
# Validate the config and print the plan without writing anything
python run_palm2gis.py -c config/my_project.yaml --dry-run

# Re-run only certain stages (e.g. after adding a wind height)
python run_palm2gis.py -c config/my_project.yaml --stages wind

# Per-timestep progress lines
python run_palm2gis.py -c config/my_project.yaml -v
```

`python -m palm2gis -c ...` works as an alternative to `run_palm2gis.py`.

---

## Files you get back, and how to load them

Everything is written under `project.output_dir` (relative to `project.root`), with the filename given in each output block:

| File | Stage | Loading it |
|------|-------|-----------|
| `voxel.netcdf` | voxel | ArcGIS Pro: **Local Scene → Add Data → Multidimensional Voxel Layer**. Pick a variable, get a true 3D volume with isosurfaces and section slicing — the GIS counterpart of ParaView volume rendering |
| `surfaces.all_elements_gpkg` | surfaces | Every surface element — ground, roofs and all four wall orientations — as PolygonZ, one row per element per timestep, with `orient`, `z` and `z_agl` columns. The largest output |
| `surfaces.all_elements_multipatch` (`.shp`) | surfaces | **ESRI MultiPatch** (shape type 31), ArcGIS's native 3D geometry type and what a 3D Object scene layer is published from. One file, every orientation — filter on `orient`. DBF limits apply: field names are capped at 10 characters and `time` is an ISO **string**, not a datetime, so time-enable it as text. `azimuth` is omitted (no numeric NULL in a shapefile); it is recoverable from `orient` |
| `surfaces.all_elements_fgdb` (`.gdb`) | surfaces | **File Geodatabase**, feature classes `surfaces` and `walls`, the same elements as PolygonZ. What ArcGIS opens natively, and unlike the shapefile it keeps real datetimes, full-length field names and NULLs — nothing is lost relative to the GeoPackage |
| `wind.netcdf` | wind | Gridded `(time, height, y, x)` with `u`, `v`, `w`, `speed`, `dir_met`, for ArcGIS Pro's **flow / vector-field** renderers. See the recipe below |
| `wind.gpkg` | wind | Layer `wind`. Symbolise with an arrow marker rotated by `dir_met` (geographic rotation — `dir_met` is the meteorological direction the wind blows FROM) and coloured or sized by `speed`; filter or animate on `time` and `height_agl` |

Because the style is stored inside each GeoPackage, those files contain a second, non-spatial `layer_styles` table. QGIS hides it automatically, but when reading a GeoPackage in your own scripts, name the layer explicitly — `gpd.read_file(path, layer="wind")` — or GeoPandas/pyogrio warns about the ambiguity. Set `advanced.qgis_styles: false` if you would rather have single-layer files.

**Opening them in QGIS.** With `advanced.qgis_styles: true` (the default) every layer is written with a `.qml` style sidecar, and GeoPackages additionally carry the style in their own `layer_styles` table. QGIS applies these automatically when the layer is added, so the time setup is already done: the vector layers have their `time` field registered with a feature duration equal to the output interval. Add the layer, open the Temporal Controller, press play. Older QGIS versions ignore the temporal block and simply open the layer untimed.

**The ArcGIS flow renderer.** Add Data → **Add Multidimensional Raster Layer**, pick the `wind.netcdf` file, tick **two** variables, and set Output Configuration to either *Vector Field (U-V)* — tick `u` and `v` — or *Vector Field (Magnitude-Direction)* — tick `speed` and `dir_met`, then set Reference System to **Geographic** and Angle Representation to **From**, which is what `dir_met` means. Both pairings are written so neither needs a regenerated file. Verify the sign convention once against `dir_met` on a case with known flow: this is the classic place for a 180° flip.

**Vertical datum.** Every z in every product is **height above sea level**: PALM's z = 0 sits at `origin_z`, which is added on output exactly as `origin_x` / `origin_y` are added to the horizontal coordinates. If your static driver leaves `origin_z` at 0 — many do — the run warns and you can set `advanced.origin_z` yourself. In an ArcGIS scene, set the layer's elevation to *absolute height from the feature's z-values*.

**Temperature units.** PALM writes `ta` and `ta_2m*` in **degrees Celsius** and `theta`, `tsurf*`, `t_surf*` in **kelvin**. `celsius` selects the output scale and the input scale is read from each variable's `units`, so an already-Celsius variable is never shifted again. `theta` is exempt: potential temperature stays kelvin either way. An unrecognised unit string raises rather than being guessed — use the `temperature_units` map to state the scale explicitly.

All layers share the CRS from `input.crs` (default `EPSG:32633`) and absolute coordinates, so they overlay each other — and any basemap — correctly. The `gis` module of `paraview_vis` (draping an orthophoto exported from QGIS) has no counterpart here for exactly that reason: in GIS the basemap is already underneath.

---

## Configuration

A config says what to read once (`input`), which part of the run to convert (`region`), and lists every output as its own block with a `file:` line. A filename produces that output; `null` or a missing block skips it, so nothing is written that the config does not name. Settings you leave out take their defaults, and an unknown or misspelt setting stops the run with a suggestion instead of being silently ignored. `template.yaml` lists every setting; `config/example_minimal.yaml` shows how short a real config can be:

```yaml
project:
    name: my_project
    root: /storage/data/palm/my_project   # base for the input paths
    output_dir: ./gis

input:
    static_driver: INPUT/my_project_static_N02
    p3d:           INPUT/my_project_p3d_N02
    data_3d:       OUTPUT/DATA_3D_NETCDF_N02
    surface_data:  OUTPUT/SURFACE_DATA_NETCDF_N02

region:
    x: [100, 300]                         # metres, local PALM coordinates
    y: [100, 300]
    height_above_highest_surface: 20
    time_from: '2023-08-24 13:00'
    time_to:   '2023-08-24 16:00'

voxel:
    variables: [u, v, w, theta, ta]
    netcdf:
        file: '{name}_voxel.nc'           # {name} = project.name

surfaces:
    variables: [theta_surface, t_surf]
    all_elements_gpkg:
        file: '{name}_surfaces.gpkg'
    all_elements_fgdb:
        file: null                        # skipped

wind:
    heights_above_surface: [2.0, 10.0]
    netcdf:
        file: '{name}_wind.nc'
```

A stage runs when its section is present; `enabled: false` switches it off without deleting its settings. Rarely needed settings (`nz`/`dz` without a p3d, `origin_z`, `surface_pressure`, `qgis_styles`, record-index time selection, ...) sit in an `advanced` block.

**Configs in the earlier layout still load.** Up to 0.19.0 a config named the internal settings directly — `domain`, `bounds`, `crs`, and per stage `file_from` / `file_to` / `write_*` flags. Such a file runs unchanged; palm2gis recognises the newer layout by its `input`, `region` or `advanced` block, and mixing the two in one file is rejected with a pointer to where each key went. The rest of this README uses the new names.

**Existing outputs are never replaced by default.** A file that is already there is kept and its product skipped, so a rerun costs almost nothing and fills in only what is missing — delete a single output and rerun to regenerate just that one. To regenerate everything, set `project.overwrite: true` or pass `--overwrite`; the run then logs a warning and replaces each output. Use `--dry-run` to see, per file, whether it would be kept, written, or replaced.

**Select the timestep window in time, not in record indices.** PALM keeps separate output intervals for 3D and surface output — `dt_data_output_av` and `dt_dosurf_av` — so a single run routinely produces a 10-minute 3D file next to an hourly surface file. A record index is then meaningful for one file only: `advanced.skip_times: 18` starts 3 hours into the 3D file and 18 hours into the surface file, and the stages quietly convert two different periods. `region.time_from` / `region.time_to` state the window once and each stage resolves it against its own time axis. Accepted forms are an absolute time (`'2023-08-24 13:00'`), a clock time (`'13:00'`, first occurrence after `origin_time`), an offset (`'18h'`, `'90min'`) or seconds since `origin_time`. Both ends are inclusive, either may be left `null` for "from the first record" / "to the last", and the run reports what each stage picked:

```
[surfaces] records 17..20 of 24 selected (2023-08-24 13:00:00 .. 2023-08-24 16:00:00).
[voxel]    records 107..125 of 144 selected (2023-08-24 13:00:00 .. 2023-08-24 16:00:00).
```

Quote clock times. Unquoted, YAML 1.1 reads `13:00` as the sexagesimal number 780, which would arrive as an offset in seconds; the config is rejected with an explanation rather than converting the wrong hours.

`advanced.skip_times` and `advanced.n_times` still select by record index (`n_times` is a count, `time_end` the exclusive end index) and are clamped per file, but they are ignored while `time_from`/`time_to` are set. When they are used and two stages end up on windows that do not even overlap, the run warns and points at `time_from`/`time_to`.

---

## Stage details

### voxel — the 3D volume in ArcGIS Pro

The 3D output is destaggered onto scalar mass centres (`u` averaged over the two x-faces, `v` over the y-faces, `w` over the z-faces), cut by `region`, and written with everything a voxel layer needs: `(time, z, y, x)` dimension order, absolute `x`/`y` with CF standard names, a `crs` grid-mapping variable carrying ASCII WKT1 (`spatial_ref`, `crs_wkt`) and an Esri PE string, plus the Esri PE string as a global attribute (where the voxel layer reads it), and CF time (`seconds since <origin_time>`). Cells inside terrain or buildings carry the `_FillValue` and render as empty space. With `advanced.voxel_drop_k0: true` (default) the z = 0 level is dropped, making the vertical spacing uniform (`dz/2, 3dz/2, …`), which voxel layers handle best.

### Which PALM file goes where

palm2gis reads three different PALM files, each named once in `input`:

| setting | PALM file | used by | typical name |
| --- | --- | --- | --- |
| `input.static_driver` | static driver (model input) | all stages; `canopy` reads LAD from it | `*_static*` |
| `input.p3d` | `_p3d` namelist (model input) | grid, output intervals, `surface_pressure` | `*_p3d*` |
| `input.data_3d` | 3D output | `voxel`, `wind` | `DATA_3D_NETCDF*` |
| `input.surface_data` | surface output | `surfaces` | `SURFACE_DATA_NETCDF*` |

If the names were changed, the contents tell them apart: the 3D output is a gridded volume (dimensions `x`, `y`, `zu_3d`) holding `u`, `v`, `w`, `theta`; the surface output lists individual surface elements along a dimension `s` with `xs`/`ys`/`zs` and `azimuth`/`zenith`. Stages check this at start and stop immediately if they were handed the wrong product, naming both what they got and what they need, rather than failing partway through with a missing-variable error.

### surfaces — roofs, ground, and walls

The unstructured surface elements are classified by orientation, the same `zenith`/`azimuth` logic `paraview_vis` used to pick quad geometry, and reconstructed as the same quads the VTU writer produced, but as georeferenced PolygonZ features. Elements with any other zenith than 0 or 90 (sloped roofs) are counted and reported in a warning rather than silently dropped.

The element layers hold **every** face of the static driver's geometry in a single layer, horizontal faces as dx × dy quads at their own height and walls as dz-tall vertical quads, each carrying the configured variables. The `orient` column (`up`, `down`, `wall_N`, `wall_E`, `wall_S`, `wall_W`) drives symbology and filtering, and `z_agl` (face height minus that column's terrain) separates roofs from ground **without** needing `buildings_2d`, which some static drivers leave empty. In a 3D scene this is the layer that shows a whole city block's surface temperature — roofs and every wall facing — at once. It is also by far the largest output: rows = elements × timesteps, and the run warns above two million.

**Actual temperature from potential temperature.** PALM's prognostic variables are potential temperatures: `theta` in the volume, `theta_surface` (or `theta_surf`) at the surface. The actual temperatures are `ta` and `t_surf`, related by `T = theta * Pi` with the Exner function integrated hydrostatically (`dPi/dz = -g/(cp*theta)`) from `surface_pressure`.

`surface_pressure` is read from the `_p3d` namelist unless it is set in the config, so it cannot drift from what the run actually used. **PALM gives that key in hPa** and palm2gis converts to Pa; a value set as `advanced.surface_pressure` is taken as Pa and wins untouched. Only when there is neither does it fall back to 101325 Pa — the standard atmosphere at *sea level*, which is wrong for any elevated domain and biases every derived temperature warm by roughly 0.85 K per 10 hPa. Prague's ~200 m is worth about 1.7 K, a uniform offset that no pattern in the output would reveal.

Listing `ta` in `voxel.variables` or `t_surf` in `surfaces.variables` *ensures* those variables rather than always computing them: **if the file already holds PALM's own `ta` or `t_surf`, that value is used**, because it is exact; only when it is absent is it reconstructed, and the reconstruction is tagged `derived_from = "theta"` so the two can never be confused. Surface temperature is evaluated at each face's own height and appears in every element layer. Both `theta_surface` and `theta_surf` are accepted, and a listed variable missing under one spelling is substituted from the other. PALM's `p` output is the *perturbation* pressure and is deliberately never used as absolute pressure.

**Absolute humidity from `q`.** Listing `ah` in `voxel.variables` adds the water vapour density in g/m³, ensured like `ta`: PALM's own `ah` is copied if the file has one, otherwise it is derived as `ah = q · p / (Rd · Tv)` with `Tv = T (1 + 0.608 q)`. `T` is PALM's `ta` when the file holds it (exact) and the theta-derived air temperature otherwise; `p` is the same hydrostatic pressure the `ta` derivation uses, so a wrong `surface_pressure` biases `ah` by about 1 % per 10 hPa. `q` is read in kg/kg, or g/kg when its `units` say so. PALM's `q` is the total water content, which equals specific humidity when no cloud water is present. The result is tagged `derived_from = "q"`.

The two names are different quantities, and only one of them is an air temperature. `ta` is PALM's **true air temperature**, an atmospheric field. `t_surf` is a **surface skin temperature**. In the earlier config layout they were `voxel.derive_ta` and `surfaces.derive_t_surf` (called `surfaces.derive_ta` up to 0.17.0, which claimed both the wrong quantity and the wrong medium).

### wind — vectors at heights above ground

For every kept grid column (`every_nth_column` thins the grid), the destaggered horizontal wind is sampled at the first model level at or above *local surface + h* for each height in `heights_above_surface`, where the local surface is terrain **plus buildings**. The 2 m level therefore follows streets and rooftops — which is what pedestrian-comfort maps want, and what a fixed-z ParaView slice cannot do (it cuts through hillsides and buildings instead). `dir_met` is the meteorological direction (degrees, direction the wind blows from), ready for arrow-marker rotation in both ArcGIS Pro and QGIS.

---

## Command-line reference

```
python run_palm2gis.py -c CONFIG [options]

  -c, --config      YAML config file (required)
  --stages S [S..]  subset of stages to run (default: stages list from config)
  --overwrite       replace existing outputs (same as project.overwrite: true)
  -n, --dry-run     validate config and print the plan, write nothing
  -v, --verbose     DEBUG-level logging (per-timestep progress)
  -q, --quiet       WARNING-level logging only
  --log-datetime    full date in timestamps (time-only is the default)
  --log-file PATH   also append the log to a file (uncoloured, dated);
                    defaults to <output_dir>/palm2gis.log
  --ignore-state    skip on file existence alone, ignoring the recorded
                    config hash (pre-0.10.0 behaviour)
  --version         print the version and exit
```

---

## Things worth knowing

**The original destaggering had a bug.** The `paraview_vis` `utils/operation_3d.py` averaged staggered components as `(var[..., 1:] + var[..., -1:]) / 2` — the `-1:` broadcasts the *last* slice against everything, so every interpolated value was averaged with the domain-edge column instead of its neighbouring face. `palm2gis` uses the correct `0.5 * (v[..., :-1] + v[..., 1:])`. Interpolated NetCDFs produced by the old tool carry the wrong velocities.

**`nz`/`dz` are not in the static driver** — that is why `input.p3d` (or `advanced.nz` + `advanced.dz`) is required, exactly as in `paraview_vis`. The `_p3d` parser matches `key = value` pairs with a regex, so the usual namelist formatting variations are fine.

**Timestamps** in every output are `origin_time` + PALM's time-in-seconds. `origin_time` comes from the static driver attribute (a PALM-style UTC-offset suffix such as `+02` is accepted and dropped); override it with `advanced.origin_time` when the attribute is missing or wrong.

**Some PALM files carry no usable times.** Surface output in particular is sometimes written with the `time` variable left at its `_FillValue` (~1e36) — the same quirk the ParaView guide flags for the last timestep, but occasionally for every record. Those values are not seconds and are never treated as such. The run warns, and the affected outputs fall back to record labels: the GeoPackages get an empty `time` field alongside a `rec` column you can still animate or filter on. In that case palm2gis reads the output interval from the `_p3d` namelist and rebuilds the timestamps itself: `dt_dosurf_av` / `dt_dosurf` for surface output, `dt_data_output_av` / `dt_data_output` for 3D, each with its matching `skip_time_*`, and the first record placed one interval in. PALM keeps those settings separate per product, so a run whose 3D file holds 10-minute records and whose surface file holds hourly ones is handled correctly without being told. The run says where the numbers came from (`times rebuilt from the p3d (dt_dosurf_av = 3600 s)`). `advanced.time_step` / `advanced.time_start` override it, and are only needed when there is no `p3d` or the namelist is wrong. Trailing fill records, the more common case, are simply dropped and need no configuration at all.

**Reruns are cheap, and config changes are noticed.** Existing outputs are kept, so a rerun after adding a wind height or a surface variable only writes the new products. `--stages` narrows it further, and `--overwrite` forces a full regeneration.

Since 0.10.0 the skip is no longer by file existence alone. A `.palm2gis_state.json` sidecar in the output directory records which config produced each file, and any file whose settings have changed is rewritten automatically — so editing `region`, `celsius`, `heights_above_surface` or a variable list and rerunning no longer silently keeps stale outputs. The hash is per stage and covers only the sections that stage reads, so changing the wind heights does not invalidate the surface layers. A file with no recorded hash is left alone, which means an output directory written by an older version is not wiped on the first run after upgrading. `--ignore-state` restores the old existence-only behaviour.

**Every output says where it came from.** Each file gets a `<name>.palm2gis.json` sidecar with the version, config path and hash, stage, source file, bounds and command line. NetCDF outputs get global attributes plus a CF `history` line.

**Voxel z is height above the domain base** (the PALM `zu` axis), not metres above sea level. The volume overlays 2D layers correctly in map view; in a Local Scene, place it with the voxel layer's vertical offset if you combine it with a real terrain surface.

---

## Origins and attribution

`palm2gis` is derived from `paraview_vis`, a PALM visualisation script by
[Martin Bureš](https://www.cs.cas.cz/staff/buresm/en) at the Institute of
Computer Science of the Czech Academy of Sciences. The domain-georeferencing
approach and the general shape of the 3D conversion come from that work, as do
the configuration concepts this tool kept (the domain, the subdomain
cut, the timestep window). The logging style follows PALM-GeM's (Martin Bureš,
Jaroslav Rešler).

The code has since been restructured into a configurable, stage-based pipeline
writing GIS-native formats rather than ParaView-oriented ones. One behavioural
difference matters for anyone holding older results: the destaggering has been
corrected, and interpolated NetCDFs produced by the earlier script carry
incorrect velocities — see [Things worth knowing](#things-worth-knowing) for
the detail.

## Authors

David Stern (IPR Praha / KFA MFF UK) — design, PALM workflow, this
implementation.
