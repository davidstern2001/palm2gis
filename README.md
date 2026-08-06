# palm2gis

**PALM outputs → ArcGIS Pro / QGIS layers** — companion to `palm_preproc` and `palm_postproc`.

From a PALM run's static driver and NetCDF outputs, `palm2gis` produces georeferenced GIS layers: a CF-compliant NetCDF that ArcGIS Pro opens directly as a 3D Voxel Layer (destaggered velocities and scalars, terrain as NoData), surface temperatures as a top-view raster time series plus 3D wall polygons, and pedestrian-level wind vectors sampled at heights above the local surface. Every layer is written in absolute projected coordinates taken from the static driver's `origin_x`/`origin_y`, so it lands in the right place on the map with no manual shifting. It began as an adaptation of `paraview_vis`, a PALM visualisation script by [Martin Bureš](https://www.cs.cas.cz/staff/buresm/en) (Institute of Computer Science, Czech Academy of Sciences), and keeps its configuration concepts (`domain`, `bounds` subdomain cut, timestep window) while targeting GIS software rather than ParaView. See [Origins and attribution](#origins-and-attribution).

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

Every run first reads the **static driver** (and optionally the `_p3d` namelist for `nz`/`dz`): it provides the grid geometry, the terrain + building elevation, the domain origin and `origin_time`, and converts the metre-based `bounds` section into grid index ranges shared by all stages. Then the enabled stages run in order:

| # | Stage | What it produces |
|---|-------|------------------|
| 1 | `voxel` | CF-compliant NetCDF for the ArcGIS Pro Voxel Layer — destaggered `u`/`v`/`w` plus chosen scalars on `(time, z, y, x)`, terrain and building cells as NoData |
| 2 | `surfaces` | One multiband GeoTIFF per surface variable (topmost roof/ground element per column, band = timestep) and the wall elements as vertical PolygonZ quads in a time-enabled GeoPackage |
| 4 | `canopy` | Plant canopy from the static driver: leaf area index and canopy height rasters, and optionally the LAD volume lifted onto absolute heights |
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
#    → set at least project.root, domain.static_driver and the
#      file_from/file_to paths of the stages you enable

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

Everything is written under `project.output_dir` (relative to `project.root`). For a `surfaces.prefix` of `<prefix>`:

| File | Stage | Loading it |
|------|-------|-----------|
| `<voxel.file_to>.nc` | voxel | ArcGIS Pro: **Local Scene → Add Data → Multidimensional Voxel Layer**. Pick a variable, get a true 3D volume with isosurfaces and section slicing — the GIS counterpart of ParaView volume rendering |
| `<prefix>surface_2d.nc` | surfaces | The top-view fields as a CF NetCDF. Open in QGIS with **Add Mesh Layer**: MDAL reads the CF time axis natively, so the Temporal Controller works immediately — the most reliable way to animate surface fields in QGIS |
| `<prefix><var>_topview.tif` | surfaces | Drag into ArcGIS Pro or QGIS; each band is one timestep (the timestamp is in the band description). What a nadir thermal camera would see: roofs where there are buildings, ground elsewhere |
| `<prefix>surfaces.gpkg` | surfaces | Every surface element — ground, roofs and all four wall orientations — as PolygonZ, one row per element per timestep, with `orient`, `z` and `z_agl` columns. Optional (`write_all_elements_gpkg`), and the largest output |
| `<prefix>walls.gpkg` | surfaces | Layer `walls`, one row per wall element per timestep. Enable time on the `time` field for the ArcGIS time slider / QGIS temporal controller; add to a Local Scene / QGIS 3D view for wall temperatures in 3D |
| `<wind.file_to>.gpkg` | wind | Layer `wind`. Symbolise with an arrow marker rotated by `dir_met` (geographic rotation — `dir_met` is the meteorological direction the wind blows FROM) and coloured or sized by `speed`; filter or animate on `time` and `height_agl` |

Because the style is stored inside each GeoPackage, those files contain a second, non-spatial `layer_styles` table. QGIS hides it automatically, but when reading a GeoPackage in your own scripts, name the layer explicitly — `gpd.read_file(path, layer="wind")` — or GeoPandas/pyogrio warns about the ambiguity. Set `qgis_styles: false` if you would rather have single-layer files.

**Animating in QGIS: prefer the mesh layer.** GeoTIFF bands carry timestamps only as band metadata, and getting QGIS to map them to times depends on the `Fixed Time Range Per Band` mode behaving as expected — it can leave the layer blank if the map has no active temporal range, and the per-band ranges are written without a timezone by QGIS itself. The NetCDF route avoids all of that: **Layer → Add Layer → Add Mesh Layer**, pick `<prefix>surface_2d.nc`, and the time axis is read from the file. Use the GeoTIFFs for static maps, print layouts and ArcGIS.

**Opening them in QGIS.** With `qgis_styles: true` (the default) every layer is written with a `.qml` style sidecar, and GeoPackages additionally carry the style in their own `layer_styles` table. QGIS applies these automatically when the layer is added, so the time setup is already done: each raster band is mapped to its own time range (the `Fixed Time Range Per Band` mode, QGIS 3.38+), and the vector layers have their `time` field registered with a feature duration equal to the output interval. Add the layer, open the Temporal Controller, press play. Older QGIS versions ignore the temporal block and simply open the layer untimed.

All layers share the CRS from `crs.palm` (default `EPSG:32633`) and absolute coordinates, so they overlay each other — and any basemap — correctly. The `gis` module of `paraview_vis` (draping an orthophoto exported from QGIS) has no counterpart here for exactly that reason: in GIS the basemap is already underneath.

---

## Configuration

Configuration is layered: hardcoded defaults are applied first, then your YAML file overrides them, and the merge is **deep** — you can override a single nested key (say `voxel.drop_k0`) without restating the rest of the block. The fully-commented `template.yaml` documents every setting; `config/example_minimal.yaml` shows how short a real project config can be:

```yaml
project:
  name: my_project
  root: /storage/data/palm/my_project   # base for all relative paths below
  output_dir: ./gis

domain:
  static_driver: INPUT/my_project_static_N02
  p3d: INPUT/my_project_p3d_N02         # or: nz + dz explicitly

voxel:
  enabled: true
  file_from: OUTPUT/DATA_3D_NETCDF_N02
  file_to: my_project_voxel.nc
  scalars: ["theta"]

surfaces:
  enabled: true
  file_from: OUTPUT/SURFACE_DATA_NETCDF_N02
  prefix: my_project_
  variables: ["theta_surface"]

wind:
  enabled: true
  file_from: OUTPUT/DATA_3D_NETCDF_N02
  file_to: my_project_wind.gpkg
  heights_agl: [2.0, 10.0]
  stride: 4
```

**Existing outputs are never replaced by default.** A file that is already there is kept and its product skipped, so a rerun costs almost nothing and fills in only what is missing — delete a single GeoTIFF and rerun to regenerate just that one. To regenerate everything, set `project.overwrite: true` or pass `--overwrite`; the run then logs a warning and replaces each output. Use `--dry-run` to see, per file, whether it would be kept, written, or replaced.

The `bounds` section cuts a subdomain in metres of the local PALM coordinate system (`x0`/`xmax`/`y0`/`ymax`, plus `zmax` as height above the highest surface inside the cut) and selects the timestep window.

**Select the timestep window in time, not in record indices.** PALM keeps separate output intervals for 3D and surface output — `dt_data_output_av` and `dt_dosurf_av` — so a single run routinely produces a 10-minute 3D file next to an hourly surface file. A record index is then meaningful for one file only: `skip_times: 18` starts 3 hours into the 3D file and 18 hours into the surface file, and the stages quietly convert two different periods. Since 0.12.0, `bounds.time_from` / `bounds.time_to` state the window once and each stage resolves it against its own time axis:

```yaml
bounds:
  time_from: "2023-08-24 13:00:00"   # absolute timestamp
  time_to:   "2023-08-24 16:00:00"
  # also accepted: 64800 (seconds since origin_time), "18h", "90min",
  #                "13:00" (clock time, first occurrence after origin_time)
```

Both ends are inclusive, either may be left `null` for "from the first record" / "to the last", and the run reports what each stage picked:

```
[surfaces] records 17..20 of 24 selected (2023-08-24 13:00:00 .. 2023-08-24 16:00:00).
[voxel]    records 107..125 of 144 selected (2023-08-24 13:00:00 .. 2023-08-24 16:00:00).
```

Quote clock times. Unquoted, YAML 1.1 reads `13:00` as the sexagesimal number 780, which would arrive as an offset in seconds; the config is rejected with an explanation rather than converting the wrong hours.

`skip_times` and `n_times` still work as before (`n_times` is a count, `time_end` the exclusive end index) and are clamped per file, but they are ignored while `time_from`/`time_to` are set. When they are used and two stages end up on windows that do not even overlap, the run now warns and points at `time_from`/`time_to`.

---

## Stage details

### voxel — the 3D volume in ArcGIS Pro

The 3D output is destaggered onto scalar mass centres (`u` averaged over the two x-faces, `v` over the y-faces, `w` over the z-faces), cut by `bounds`, and written with everything a voxel layer needs: `(time, z, y, x)` dimension order, absolute `x`/`y` with CF standard names, a `crs` grid-mapping variable carrying both the OGC WKT and an Esri PE string, and CF time (`seconds since <origin_time>`). Cells inside terrain or buildings carry the `_FillValue` and render as empty space. With `drop_k0: true` (default) the z = 0 level is dropped, making the vertical spacing uniform (`dz/2, 3dz/2, …`), which voxel layers handle best.

### Which PALM file goes where

palm2gis reads three different PALM files. `file_from` appears in four stages but does not mean the same file each time:

| setting | PALM file | typical name |
| --- | --- | --- |
| `domain.static_driver` | static driver (model input) | `*_static*` |
| `domain.p3d` | `_p3d` namelist (model input) | `*_p3d*` |
| `voxel.file_from` | 3D output | `DATA_3D_NETCDF*` |
| `wind.file_from` | 3D output — the same file | `DATA_3D_NETCDF*` |
| `surfaces.file_from` | surface output | `SURFACE_DATA_NETCDF*` |
| `canopy.file_from` | leave `null` → the static driver | — |

If the names were changed, the contents tell them apart: the 3D output is a gridded volume (dimensions `x`, `y`, `zu_3d`) holding `u`, `v`, `w`, `theta`; the surface output lists individual surface elements along a dimension `s` with `xs`/`ys`/`zs` and `azimuth`/`zenith`. Stages check this at start and stop immediately if they were handed the wrong product, naming both what they got and what they need, rather than failing partway through with a missing-variable error.

### surfaces — roofs, ground, and walls

The unstructured surface elements are classified by orientation, the same `zenith`/`azimuth` logic `paraview_vis` used to pick quad geometry. Upward-facing elements (zenith 0) are mapped to their grid column and the **topmost** one per column wins — ground in the streets, the roof on buildings — giving a clean 2D raster per variable. Wall elements (zenith 90) are reconstructed as the same vertical quads the VTU writer produced, but as georeferenced PolygonZ features. Elements with any other zenith (sloped roofs) are counted and reported in a warning rather than silently dropped.

Both of those are selective: the raster keeps only the topmost face per column, the walls layer only the vertical ones. `write_all_elements_gpkg: true` adds the complete picture instead — **every** face of the static driver's geometry in a single layer, horizontal faces as dx × dy quads at their own height and walls as dz-tall vertical quads, each carrying the configured variables. The `orient` column (`up`, `down`, `wall_N`, `wall_E`, `wall_S`, `wall_W`) drives symbology and filtering, and `z_agl` (face height minus that column's terrain) separates roofs from ground **without** needing `buildings_2d`, which some static drivers leave empty. In a 3D scene this is the layer that shows a whole city block's surface temperature — roofs and every wall facing — at once. It is also by far the largest output: rows = elements × timesteps, and the run warns above two million.

**Actual temperature from potential temperature.** PALM's prognostic variables are potential temperatures: `theta` in the volume, `theta_surface` (or `theta_surf`) at the surface. The actual temperatures are `ta` and `t_surf`, related by `T = theta * Pi` with the Exner function integrated hydrostatically (`dPi/dz = -g/(cp*theta)`) from `surface_pressure`.

`surface_pressure` is read from the `_p3d` namelist unless it is set in the config, so it cannot drift from what the run actually used. **PALM gives that key in hPa** and palm2gis converts to Pa; a value set in the YAML is taken as Pa and wins untouched. Only when there is neither does it fall back to 101325 Pa — the standard atmosphere at *sea level*, which is wrong for any elevated domain and biases every derived temperature warm by roughly 0.85 K per 10 hPa. Prague's ~200 m is worth about 1.7 K, a uniform offset that no pattern in the output would reveal.

`voxel.derive_ta` and `surfaces.derive_ta` *ensure* those variables rather than always computing them: **if the file already holds PALM's own `ta` or `t_surf`, that value is used**, because it is exact; only when it is absent is it reconstructed, and the reconstruction is tagged `derived_from = "theta"` so the two can never be confused. Surface temperature is evaluated at each face's own height and appears in the raster, the mesh NetCDF and the element layers alike. Both `theta_surface` and `theta_surf` are accepted, and a listed variable missing under one spelling is substituted from the other. PALM's `p` output is the *perturbation* pressure and is deliberately never used as absolute pressure.

### wind — vectors at heights above ground

For every kept grid column (`stride` thins the grid), the destaggered horizontal wind is sampled at the first model level at or above *local surface + h* for each height in `heights_agl`, where the local surface is terrain **plus buildings**. The 2 m level therefore follows streets and rooftops — which is what pedestrian-comfort maps want, and what a fixed-z ParaView slice cannot do (it cuts through hillsides and buildings instead). `dir_met` is the meteorological direction (degrees, direction the wind blows from), ready for arrow-marker rotation in both ArcGIS Pro and QGIS.

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

**The walls GeoPackage grows as n_walls × n_times rows.** For a large domain with many timesteps that can get heavy; trim `bounds` and `n_times`, or disable it with `surfaces.write_walls_gpkg: false`. A per-building aggregation (mean wall value per footprint) is a possible future lighter alternative.

**`nz`/`dz` are not in the static driver** — that is why `domain.p3d` (or explicit `nz` + `dz`) is required, exactly as in `paraview_vis`. The `_p3d` parser matches `key = value` pairs with a regex, so the usual namelist formatting variations are fine.

**Timestamps** in every output are `origin_time` + PALM's time-in-seconds. `origin_time` comes from the static driver attribute (a PALM-style UTC-offset suffix such as `+02` is accepted and dropped); override it with `domain.origin_time` when the attribute is missing or wrong.

**Some PALM files carry no usable times.** Surface output in particular is sometimes written with the `time` variable left at its `_FillValue` (~1e36) — the same quirk the ParaView guide flags for the last timestep, but occasionally for every record. Those values are not seconds and are never treated as such. The run warns, and the affected outputs fall back to record labels: raster bands become `record 0`, `record 1`, …, and the GeoPackages get an empty `time` field alongside a `rec` column you can still animate or filter on. In that case palm2gis reads the output interval from the `_p3d` namelist and rebuilds the timestamps itself: `dt_dosurf_av` / `dt_dosurf` for surface output, `dt_data_output_av` / `dt_data_output` for 3D, each with its matching `skip_time_*`, and the first record placed one interval in. PALM keeps those settings separate per product, so a run whose 3D file holds 10-minute records and whose surface file holds hourly ones is handled correctly without being told. The run says where the numbers came from (`times rebuilt from the p3d (dt_dosurf_av = 3600 s)`). `bounds.time_step` / `bounds.time_start` override it, and are only needed when there is no `p3d` or the namelist is wrong. Trailing fill records, the more common case, are simply dropped and need no configuration at all.

**Reruns are cheap, and config changes are noticed.** Existing outputs are kept, so a rerun after adding a wind height or a surface variable only writes the new products. `--stages` narrows it further, and `--overwrite` forces a full regeneration.

Since 0.10.0 the skip is no longer by file existence alone. A `.palm2gis_state.json` sidecar in the output directory records which config produced each file, and any file whose settings have changed is rewritten automatically — so editing `bounds`, `celsius`, `heights_agl` or a variable list and rerunning no longer silently keeps stale outputs. The hash is per stage and covers only the sections that stage reads, so changing the wind heights does not invalidate the surface rasters. A file with no recorded hash is left alone, which means an output directory written by an older version is not wiped on the first run after upgrading. `--ignore-state` restores the old existence-only behaviour.

**Every output says where it came from.** Each file gets a `<name>.palm2gis.json` sidecar with the version, config path and hash, stage, source file, bounds and command line. GeoTIFFs carry the same fields in their own metadata tags (so `gdalinfo` and QGIS layer properties show them, and they survive the file being copied elsewhere), and NetCDF outputs get global attributes plus a CF `history` line.

**Voxel z is height above the domain base** (the PALM `zu` axis), not metres above sea level. The volume overlays 2D layers correctly in map view; in a Local Scene, place it with the voxel layer's vertical offset if you combine it with a real terrain surface.

---

## Origins and attribution

`palm2gis` is derived from `paraview_vis`, a PALM visualisation script by
[Martin Bureš](https://www.cs.cas.cz/staff/buresm/en) at the Institute of
Computer Science of the Czech Academy of Sciences. The domain-georeferencing
approach and the general shape of the 3D conversion come from that work, as do
the configuration concepts this tool kept (`domain`, the `bounds` subdomain
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
