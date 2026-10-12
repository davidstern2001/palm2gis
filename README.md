# palm2gis

**PALM outputs → ArcGIS Pro / QGIS layers** — companion to `palm_preproc` and `palm_postproc`.

From a PALM run's static driver and NetCDF outputs, `palm2gis` produces georeferenced GIS layers: a CF-compliant NetCDF that ArcGIS Pro opens directly as a 3D Voxel Layer (destaggered velocities and scalars, terrain as NoData), surface temperatures on every ground, roof and wall element as time-enabled 3D polygons, gridded wind at absolute heights for ArcGIS Pro's flow renderer, and the plant canopy's leaf area density as a voxel NetCDF. Every layer is written in absolute projected coordinates taken from the static driver's `origin_x`/`origin_y`, so it lands in the right place on the map with no manual shifting. It began as an adaptation of `paraview_vis`, a PALM visualisation script by [Martin Bureš](https://www.cs.cas.cz/staff/buresm/en) (Institute of Computer Science, Czech Academy of Sciences), and keeps its configuration concepts (domain, subdomain cut, timestep window) while targeting GIS software rather than ParaView. See [Origins and attribution](#origins-and-attribution).

---

## Table of contents

- [How it works (the four stages)](#how-it-works-the-four-stages)
- [Installing and running](#installing-and-running)
- [Files you get back, and how to load them](#files-you-get-back-and-how-to-load-them)
- [Configuration](#configuration)
- [Stage details](#stage-details)
  - [voxel — the 3D volume in ArcGIS Pro](#voxel--the-3d-volume-in-arcgis-pro)
  - [surfaces — roofs, ground, and walls](#surfaces--roofs-ground-and-walls)
  - [wind — gridded at absolute heights](#wind--gridded-at-absolute-heights)
- [Command-line reference](#command-line-reference)
- [Things worth knowing](#things-worth-knowing)
- [Origins and attribution](#origins-and-attribution)
- [Authors](#authors)

---

## How it works (the four stages)

Every run first reads the **static driver** (and optionally the `_p3d` namelist for `nz`/`dz`): it provides the grid geometry, the terrain + building elevation, the domain origin and `origin_time`, and converts the metre-based `region` section into grid index ranges shared by all stages. Then the enabled stages run in order:

| # | Stage | What it produces |
|---|-------|------------------|
| 1 | `voxel` | CF-compliant NetCDF for the ArcGIS Pro Voxel Layer — destaggered `u`/`v`/`w` plus chosen scalars on `(time, z, y, x)`, terrain and building cells as NoData |
| 2 | `surfaces` | Every ground, roof and wall element as a time-enabled PolygonZ quad: GeoPackage (mainly QGIS) and File Geodatabase (local ArcGIS Pro); the walls alone as MultiPatch (web ArcGIS) |
| 3 | `wind` | Gridded netCDF per absolute height with `u`, `v`, `w`, `speed`, `dir_met` per timestep, for ArcGIS Pro's flow / vector-field renderers — all heights in one folder |
| 4 | `canopy` | Plant canopy from the static driver: the LAD volume lifted onto absolute heights, as a voxel NetCDF |

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
| `surfaces.all_elements_gpkg` (`.gpkg`) | surfaces | **Mainly for QGIS.** Every surface element — ground, roofs and all four wall orientations — as PolygonZ, one row per element per timestep, with `orient`, `z` and `z_agl` columns and the time setup embedded as a QGIS style. The largest output |
| `surfaces.all_elements_fgdb` (`.gdb`) | surfaces | **For local ArcGIS Pro.** File Geodatabase, feature classes `surfaces` and `walls`, the same elements as PolygonZ, one row per element per timestep. What ArcGIS opens natively, with real datetimes, full-length field names and NULLs — nothing is lost relative to the GeoPackage, and the time slider works on `time` |
| `surfaces.walls_multipatch` (`.shp`) | surfaces | **For web ArcGIS.** The walls alone as ESRI MultiPatch — what a 3D object scene layer is published from (ArcGIS Pro: Share → Web Layer, or Create 3D Object Scene Layer Content). One feature per wall instead of one per wall per timestep; each variable and record is its own field. See below |
| `wind.netcdf` | wind | **One file per absolute height** (m a.s.l.): every level of the voxel's range by default, or those in `wind.netcdf.heights`, `(time, z, y, x)` with `u`, `v`, `w`, `speed`, `dir_met`, for ArcGIS Pro's **flow / vector-field** renderers; collected in one folder named after the file, `<file>/<file>_z<height>m.nc` (or put `{z}` in the name). A height between two PALM levels is interpolated linearly; cells where either level is inside terrain or a building are NoData. See the recipe below |

Because the style is stored inside each GeoPackage, those files contain a second, non-spatial `layer_styles` table. QGIS hides it automatically, but when reading a GeoPackage in your own scripts, name the layer explicitly — `gpd.read_file(path, layer="surfaces")` — or GeoPandas/pyogrio warns about the ambiguity. Set `advanced.qgis.styles: false` if you would rather have single-layer files.

**Opening them in QGIS.** With `advanced.qgis.styles: true` (the default) every layer is written with a `.qml` style sidecar, and GeoPackages additionally carry the style in their own `layer_styles` table. QGIS applies these automatically when the layer is added, so the time setup is already done: the vector layers have their `time` field registered with a feature duration equal to the output interval. Add the layer, open the Temporal Controller, press play. Older QGIS versions ignore the temporal block and simply open the layer untimed.

**The ArcGIS flow renderer.** Add Data → **Add Multidimensional Raster Layer**, pick one of the `wind.netcdf` files (one per absolute height), tick **two** variables, and set Output Configuration to either *Vector Field (U-V)* — tick `u` and `v` — or *Vector Field (Magnitude-Direction)* — tick `speed` and `dir_met`, then set Reference System to **Geographic** and Angle Representation to **From**, which is what `dir_met` means. Both pairings are written so neither needs a regenerated file. Verify the sign convention once against `dir_met` on a case with known flow: this is the classic place for a 180° flip.

**Horizontal CRS in the netCDF outputs.** The `crs` grid-mapping variable carries the CRS as WKT only (`crs_wkt`, `spatial_ref`, `esri_pe_string`, `EPSG_code`), plus a global `esri_pe_string`. ArcGIS Pro reads CF grid-mapping parameters *before* any WKT. CF has no name, EPSG code or datum, so with those parameters present ArcGIS shows a nameless "Custom" Transverse_Mercator on an assumed WGS 1984 instead of *WGS 1984 UTM zone 33N (32633)*. Without them it reads `crs_wkt` and recognises the WKID, and a non-WGS84 datum (S-JTSK) is honoured. GDAL/QGIS read `spatial_ref` either way. `advanced.netcdf.crs: cf` adds the CF parameters back for CF-only readers (xarray/cf-xarray), at the cost of the custom CRS in ArcGIS. Strict CF checkers flag the default `crs` variable for lacking `grid_mapping_name`.

**Vertical datum.** Every z in every product is **height above sea level**: PALM's z = 0 sits at `origin_z`, which is added on output exactly as `origin_x` / `origin_y` are added to the horizontal coordinates. If your static driver leaves `origin_z` at 0 — many do — the run warns and you can set `advanced.grid.origin_z` yourself. In an ArcGIS scene, set the layer's elevation to *absolute height from the feature's z-values*.

**Temperature units.** PALM writes `ta` and `ta_2m*` in **degrees Celsius** and `theta`, `tsurf*`, `t_surf*` (and palm2gis' `t_surface`) in **kelvin**. `advanced.units.temperature` selects the output scale for actual temperatures and is **`C` by default** (it was kelvin before 0.25.0, under the `celsius` flag, which still works). The input scale is read from each variable's `units`, so an already-Celsius variable is never shifted again. Potential temperature is exempt: `theta` and `theta_surface` are kelvin-defined, `advanced.units.potential_temperature` accepts only `K`, and a Celsius value is rejected rather than silently producing something that is not a potential temperature. An unrecognised unit string raises rather than being guessed — use `advanced.units.overrides` to state the scale explicitly.

**Output in another CRS (`project.output_crs`).** An ArcGIS scene in WGS 84 wants files that *are* in EPSG:4326, not files it can reproject on the fly. `project.output_crs: 'EPSG:4326'` (any CRS pyproj knows works) writes every product in that CRS, and the two kinds of product are handled differently:

- *Vector layers* (surfaces, in both formats) are transformed **exactly**: every vertex goes through pyproj, z stays metres above sea level.
- *Gridded netCDF* (voxel, wind, canopy) **cannot** be transformed exactly: PALM's grid is regular in UTM and not in lon/lat (at Prague UTM 33N is turned about 0.45° against the meridians, and a cell is not square in degrees). ArcGIS takes a voxel source only with 1-D coordinates and no auxiliary variables, so the data is **resampled** onto a regular `(lat, lon)` grid whose spacing is one native cell at the centre of the cut. The corners of the lon/lat box outside the (rotated) PALM domain are NoData. `advanced.netcdf.resampling` chooses how:
  - `conservative` (default): PALM is a finite-volume model, so each value is the *average* over its cell. Each output cell gets the average of the PALM cells it overlaps, weighted by the exact overlap area and taken over air cells only. Area means and totals are kept (the canopy's total leaf area to 10⁻⁸, the float32 storage precision), and a cell less than half air is NoData, so walls stay within half a cell of their true position. It is first-order: a smooth field is reproduced to ~10⁻⁴ of a cell at Prague's rotation. With `scipy` installed this takes about 2 s per variable and record on a 1000 × 1000 × 60 cut, roughly 4× nearest neighbour; without it, the numpy path is about 4× slower again.
  - `nearest`: each output cell copies the PALM cell under its centre. Values are exactly PALM's, positions good to half a cell.
- *Wind components* `u`, `v` are turned from PALM's grid axes to true east/north (the local meridian convergence, found per column), and `speed` and `dir_met` follow; `w`, `speed_3d` and scalars are unchanged.

The voxel file then carries `lon` / `lat` (`degrees_east` / `degrees_north`) and a `crs` variable with the EPSG:4326 WKT, which ArcGIS Pro reads as WKID 4326. The default, `null`, writes everything in `input.crs` exactly as before.

Without `output_crs`, `u`, `v` and `dir_met` refer to the UTM **grid** north, which at Prague is 0.45° from true north.

All layers share the CRS from `project.output_crs`, or `input.crs` (default `EPSG:32633`) when it is not set, and absolute coordinates, so they overlay each other — and any basemap — correctly. The `gis` module of `paraview_vis` (draping an orthophoto exported from QGIS) has no counterpart here for exactly that reason: in GIS the basemap is already underneath.

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
    height_above_surface: 20              # above terrain and rooftops

time:
    from: '2023-08-24 13:00'
    to:   '2023-08-24 16:00'

voxel:
    variables: [u, v, w, theta, ta]
    netcdf:
        file: '{name}_voxel.nc'           # {name} = project.name

surfaces:
    variables: [theta_surface, t_surface]
    all_elements_gpkg:
        file: '{name}_surfaces.gpkg'
    all_elements_fgdb:
        file: null                        # skipped

wind:
    netcdf:
        file: '{name}_wind.nc'            # -> {name}_wind/{name}_wind_z200m.nc, ...
        #heights: [200.0, 210.0]          # m a.s.l.; default: voxel levels
```

A stage runs when its section is present; `enabled: false` switches it off without deleting its settings. Rarely needed settings (`nz`/`dz` without a p3d, `origin_z`, `surface_pressure`, `qgis_styles`, record-index time selection, ...) sit in an `advanced` block.

**Configs in the earlier layout still load.** Up to 0.19.0 a config named the internal settings directly — `domain`, `bounds`, `crs`, and per stage `file_from` / `file_to` / `write_*` flags. Such a file runs unchanged; palm2gis recognises the newer layout by its `input`, `region` or `advanced` block, and mixing the two in one file is rejected with a pointer to where each key went. The rest of this README uses the new names.

**Existing outputs are never replaced by default.** A file that is already there is kept and its product skipped, so a rerun costs almost nothing and fills in only what is missing — delete a single output and rerun to regenerate just that one. To regenerate everything, set `project.overwrite: true` or pass `--overwrite`; the run then logs a warning and replaces each output. Use `--dry-run` to see, per file, whether it would be kept, written, or replaced.

**Select the timestep window in time, not in record indices.** PALM keeps separate output intervals for 3D and surface output — `dt_data_output_av` and `dt_dosurf_av` — so a single run routinely produces a 10-minute 3D file next to an hourly surface file. A record index is then meaningful for one file only: `advanced.records.skip: 18` starts 3 hours into the 3D file and 18 hours into the surface file, and the stages quietly convert two different periods. `time.from` / `time.to` state the window once and each stage resolves it against its own time axis. Accepted forms are an absolute time (`'2023-08-24 13:00'`), a clock time (`'13:00'`, first occurrence after `origin_time`), an offset (`'18h'`, `'90min'`) or seconds since `origin_time`. Both ends are inclusive, either may be left `null` for "from the first record" / "to the last", and the run reports what each stage picked:

```
[surfaces] records 17..20 of 24 selected (2023-08-24 13:00:00 .. 2023-08-24 16:00:00).
[voxel]    records 107..125 of 144 selected (2023-08-24 13:00:00 .. 2023-08-24 16:00:00).
```

Quote clock times. Unquoted, YAML 1.1 reads `13:00` as the sexagesimal number 780, which would arrive as an offset in seconds; the config is rejected with an explanation rather than converting the wrong hours.

`advanced.records.skip` and `advanced.records.count` still select by record index (`count` is a count, `end` the exclusive end index) and are clamped per file, but they are ignored while `time.from`/`time.to` are set. When they are used and two stages end up on windows that do not even overlap, the run warns and points at `time.from`/`time.to`.

---

## Stage details

### voxel — the 3D volume in ArcGIS Pro

The 3D output is destaggered onto scalar mass centres (`u` averaged over the two x-faces, `v` over the y-faces, `w` over the z-faces), cut by `region`, and written with everything a voxel layer needs: `(time, z, y, x)` dimension order, absolute `x`/`y` with CF standard names, a `crs` grid-mapping variable carrying ASCII WKT1 (`spatial_ref`, `crs_wkt`) and an Esri PE string, plus the Esri PE string as a global attribute (where the voxel layer reads it), and CF time (`seconds since <origin_time>`). Cells inside terrain or buildings carry the `_FillValue` and render as empty space. With `advanced.netcdf.voxel_drop_k0: true` (default) the z = 0 level is dropped, making the vertical spacing uniform (`dz/2, 3dz/2, …`), which voxel layers handle best. Each variable also carries CF `actual_range`, its minimum and maximum over **all** written records (the log prints them too): ArcGIS Pro computes a voxel layer's statistics from the loaded time step, so its default stretch fits one record and hides the change over time. That full range is set by a few extreme cells, though, so palm2gis also bins every value of every record and writes `<name>_voxel_ranges.csv` (min, p1 … p99, max and a **suggested Min/Max = p2 … p98** per variable), `<name>_voxel_histogram.csv` and, with `matplotlib` installed, `<name>_voxel_histogram.png` (one histogram per variable with p2/p98 marked). Type the suggested Min/Max into the voxel layer's symbology; each variable also carries it as `palm2gis_display_range`.

### Which PALM file goes where

palm2gis reads three different PALM files, each named once in `input`:

| setting | PALM file | used by | typical name |
| --- | --- | --- | --- |
| `input.static_driver` | static driver (model input) | all stages; `canopy` reads LAD from it | `*_static*` |
| `input.p3d` | `_p3d` namelist (model input) | grid, output intervals, `surface_pressure`, `pt_surface` | `*_p3d*` |
| `input.dynamic_driver` | dynamic driver (model input, optional) | root-domain `surface_pressure`; found automatically next to the p3d | `*_dynamic` |
| `input.data_3d` | 3D output | `voxel`, `wind` | `DATA_3D_NETCDF*` |
| `input.surface_data` | surface output | `surfaces` | `SURFACE_DATA_NETCDF*` |

If the names were changed, the contents tell them apart: the 3D output is a gridded volume (dimensions `x`, `y`, `zu_3d`) holding `u`, `v`, `w`, `theta`; the surface output lists individual surface elements along a dimension `s` with `xs`/`ys`/`zs` and `azimuth`/`zenith`. Stages check this at start and stop immediately if they were handed the wrong product, naming both what they got and what they need, rather than failing partway through with a missing-variable error.

### surfaces — roofs, ground, and walls

The unstructured surface elements are classified by orientation, the same `zenith`/`azimuth` logic `paraview_vis` used to pick quad geometry, and reconstructed as the same quads the VTU writer produced, but as georeferenced PolygonZ features. Elements with any other zenith than 0 or 90 (sloped roofs) are counted and reported in a warning rather than silently dropped.

The element layers hold **every** face of the static driver's geometry in a single layer, horizontal faces as dx × dy quads at their own height and walls as dz-tall vertical quads, each carrying the configured variables. The `orient` column (`up`, `down`, `wall_N`, `wall_E`, `wall_S`, `wall_W`) drives symbology and filtering, and `z_agl` (face height minus that column's terrain) separates roofs from ground **without** needing `buildings_2d`, which some static drivers leave empty. In a 3D scene this is the layer that shows a whole city block's surface temperature — roofs and every wall facing — at once. It is also by far the largest output: rows = elements × timesteps, and the run warns above two million.

**Which format for what.**

| Output | For | Why |
| --- | --- | --- |
| `all_elements_gpkg` (`.gpkg`) | QGIS, mainly | Open format, temporal style embedded, every element |
| `all_elements_fgdb` (`.gdb`) | local ArcGIS Pro | ArcGIS's native container: datetimes for the time slider, full field names, NULLs, every element |
| `walls_multipatch` (`.shp`) | web ArcGIS (scene viewer, ArcGIS Maps SDK for JavaScript) | Multipatch is what a 3D object scene layer is built from, and a scene layer is tiled with levels of detail, so a browser can draw tens of thousands of walls |

**Walls for the web (`walls_multipatch`).** A browser cannot draw hundreds of thousands of polygons as a feature layer, and an elevation layer cannot hold a vertical face, so the web output is the walls alone — the added value over PALM's own top-view output of ground and roofs — as a multipatch, published as a 3D object scene layer. Three things keep it small and publishable:

- **Walls only.** Ground and roofs come from PALM's top-view output.
- **One feature per wall.** Each variable at each selected record is its own field, `<code>_<nnn>`: `<code>` is the variable shortened to fit DBF's 10 characters (`t_surface` → `t_surf`, `theta_surface` → `theta`) and `nnn` counts the selected records from `000`. `<name>_fields.csv` beside the shapefile maps every field to its variable and time. Animate in the web app by switching the field the renderer reads; scene layers do not use the time slider. A shapefile holds at most 255 fields, so variables × records can be at most 250 — narrow `time` or `variables` otherwise (the run says so).
- **Fixed fields** `elem_id`, `orient` (`wall_N` / `wall_E` / `wall_S` / `wall_W`), `azimuth`, `z` (m a.s.l.) and `z_agl`. NoData is **−9999**: a shapefile has no NULL, so exclude it in the renderer's range.

**Actual temperature from potential temperature.** PALM's prognostic variables are potential temperatures: `theta` in the volume, `theta_surface` (or `theta_surf`) at the surface. The actual temperatures are `ta` and `t_surface`, related by `T = theta * Pi`. The Exner function is **PALM's own basic state**, because PALM converted between the two with it: a constant-theta column fixed at initialisation, `Pi(z) = (surface_pressure / p0)^kappa - g z / (cp pt_surface)`, linear in z and unchanged for the whole run. (On real output `ta / theta` fits that line to under 0.0001 K; integrating the current theta profile instead, as up to 0.21 did, was off by about 0.3 K per km.)

`surface_pressure` and `pt_surface` must be what the run used, not the real atmosphere, and PALM takes them from different places per domain:

| domain | `surface_pressure` | `pt_surface` |
| --- | --- | --- |
| root | dynamic driver `surface_forcing_surface_pressure`, **first value only** (`input.dynamic_driver`, or `<job>_dynamic` next to the p3d); else its p3d | its p3d |
| nested child (`*_p3d_N02`) | its p3d, else the **parent** p3d (`*_p3d`) | same order |

**The namelist gives `surface_pressure` in hPa**; palm2gis converts to Pa. `advanced.basic_state.surface_pressure` (Pa) and `advanced.basic_state.pt_surface` (K) override everything. When nothing is found, PALM's defaults (1013.25 hPa, 300 K) are used with a warning — roughly 0.85 K per 10 hPa of error. Whenever the 3D input holds both `ta` and `theta`, their ratio is PALM's exact Exner function: palm2gis fits it, confirms the resolved values, and replaces them if they disagree by more than 0.05 K (a config value is kept, with a warning). The values used are written as `palm2gis_surface_pressure` / `palm2gis_pt_surface` attributes on derived variables.

Listing `ta` in `voxel.variables` or `t_surface` in `surfaces.variables` *ensures* those variables rather than always computing them: **if the file already holds PALM's own `ta` or surface temperature (`t_surface` or `t_surf`), that value is used**, because it is exact; only when it is absent is it reconstructed, and the reconstruction is tagged `derived_from = "theta"` so the two can never be confused. Surface temperature is evaluated at each face's own height and appears in every element layer. Both `theta_surface` and `theta_surf` are accepted, and a listed variable missing under one spelling is substituted from the other. PALM's `p` output is the *perturbation* pressure and is deliberately never used as absolute pressure.

**3D wind speed.** Listing `speed_3d` in `voxel.variables` adds `sqrt(u² + v² + w²)` in m/s, computed from the components after destaggering to cell centres. The file must hold `u`, `v` and `w`; they are written only if also listed. It is the full 3D magnitude and is tagged `derived_from = "u v w"`. Listing `speed` instead (or as well) adds the **horizontal** speed `sqrt(u² + v²)`, the same quantity as `speed` in the wind files (CF `wind_speed`).

**Absolute humidity from `q`.** Listing `ah` in `voxel.variables` adds the water vapour density in g/m³, ensured like `ta`: PALM's own `ah` is copied if the file has one, otherwise it is derived as `ah = q · p / (Rd · Tv)` with `Tv = T (1 + 0.608 q)`. `T` is PALM's `ta` when the file holds it (exact) and the theta-derived air temperature otherwise; `p` comes from the same basic state the `ta` derivation uses, so a wrong `surface_pressure` biases `ah` by about 1 % per 10 hPa. `q` is read in kg/kg, or g/kg when its `units` say so. PALM's `q` is the total water content, which equals specific humidity when no cloud water is present. The result is tagged `derived_from = "q"`.

The two names are different quantities, and only one of them is an air temperature. `ta` is PALM's **true air temperature**, an atmospheric field. `t_surface` is a **surface skin temperature**, named to match `theta_surface`. In the earlier config layout they were `voxel.derive_ta` and `surfaces.derive_t_surface` (called `surfaces.derive_ta` up to 0.17.0, which claimed both the wrong quantity and the wrong medium, and `derive_t_surf` up to 0.21, when the output was `t_surf`; the old names, and `t_surf` in `surfaces.variables`, still work with a warning).

### wind — gridded at absolute heights

**Gridded netCDF (`wind.netcdf`).** One file per **absolute height** — metres above sea level, the same z as the voxel file, not heights above the ground. By default the heights are the voxel's own levels: every PALM level from the bottom of the cut up to `region.height_above_surface` above the highest surface (the voxel's k = 0 rule included), so each wind file lies exactly on one voxel z. `wind.netcdf.heights` lists specific heights instead. Each file is `(time, z, y, x)` with `z` of length 1 holding that height, all selected records on the time axis, and `u`, `v`, `w`, `speed` (horizontal) and `dir_met`. A listed height that falls between two PALM levels is interpolated linearly between them, so the file holds the wind at exactly that height; a cell is NoData where either of the two levels is inside terrain or a building. A height outside the file's levels is skipped with a warning. All files go into one folder named after `file` without its extension — `{name}_wind.nc` writes `{name}_wind/{name}_wind_z<height>m.nc` — so a run with several dozen levels leaves the output directory tidy. Put `{z}` in `file` to place the height in the name yourself (`{name}_wind_{z}m.nc` → `{name}_wind/{name}_wind_200m.nc`).

The point / 3D-arrow GeoPackage (`wind.gpkg`) was removed in 0.29.0: the flow renderer draws its own arrows and streamlines from the netCDF. A config that still has `wind.gpkg` or `wind.heights_above_surface` runs, and the run notes that they are ignored.

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
  --log-file PATH   also write this run's log to a file, replacing the
                    previous run's (uncoloured, dated);
                    defaults to <output_dir>/palm2gis.log
  --ignore-state    skip on file existence alone, ignoring the recorded
                    config hash (pre-0.10.0 behaviour)
  --version         print the version and exit
```

---

## Things worth knowing

**The original destaggering had a bug.** The `paraview_vis` `utils/operation_3d.py` averaged staggered components as `(var[..., 1:] + var[..., -1:]) / 2` — the `-1:` broadcasts the *last* slice against everything, so every interpolated value was averaged with the domain-edge column instead of its neighbouring face. `palm2gis` uses the correct `0.5 * (v[..., :-1] + v[..., 1:])`. Interpolated NetCDFs produced by the old tool carry the wrong velocities.

**`nz`/`dz` are not in the static driver** — that is why `input.p3d` (or `advanced.grid.nz` + `advanced.grid.dz`) is required, exactly as in `paraview_vis`. The `_p3d` parser matches `key = value` pairs with a regex, so the usual namelist formatting variations are fine.

**Timestamps** in every output are `origin_time` + PALM's time-in-seconds. `origin_time` comes from the static driver attribute (a PALM-style UTC-offset suffix such as `+02` is accepted and dropped); override it with `advanced.grid.origin_time` when the attribute is missing or wrong.

**Some PALM files carry no usable times.** Surface output in particular is sometimes written with the `time` variable left at its `_FillValue` (~1e36) — the same quirk the ParaView guide flags for the last timestep, but occasionally for every record. Those values are not seconds and are never treated as such. The run warns, and the affected outputs fall back to record labels: the GeoPackages get an empty `time` field alongside a `rec` column you can still animate or filter on. In that case palm2gis reads the output interval from the `_p3d` namelist and rebuilds the timestamps itself: `dt_dosurf_av` / `dt_dosurf` for surface output, `dt_data_output_av` / `dt_data_output` for 3D, each with its matching `skip_time_*`, and the first record placed one interval in. PALM keeps those settings separate per product, so a run whose 3D file holds 10-minute records and whose surface file holds hourly ones is handled correctly without being told. The run says where the numbers came from (`times rebuilt from the p3d (dt_dosurf_av = 3600 s)`). `advanced.time_axis.step` / `advanced.time_axis.start` override it, and are only needed when there is no `p3d` or the namelist is wrong. Trailing fill records, the more common case, are simply dropped and need no configuration at all.

**Reruns are cheap, and config changes are noticed.** Existing outputs are kept, so a rerun after adding a wind height or a surface variable only writes the new products. `--stages` narrows it further, and `--overwrite` forces a full regeneration.

Since 0.10.0 the skip is no longer by file existence alone. A `.palm2gis_state.json` sidecar in the output directory records which config produced each file, and any file whose settings have changed is rewritten automatically — so editing `region`, `advanced.units`, `wind.netcdf.heights` or a variable list and rerunning no longer silently keeps stale outputs. The hash is per stage and covers only the sections that stage reads, so changing the wind heights does not invalidate the surface layers. A file with no recorded hash is left alone, which means an output directory written by an older version is not wiped on the first run after upgrading. `--ignore-state` restores the old existence-only behaviour.

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

David Stern ([IPR Praha](https://iprpraha.cz/) /
[KFA MFF UK](https://www.mff.cuni.cz/en/faculty/organizational-structure/department?code=115))
— design, PALM workflow, this implementation.
