# Changelog

## 0.16.0

**Wind in 3D**

- New `wind_3d` layer in the wind GeoPackage: LineStringZ arrows at the
  TRUE sample height (local surface + `height_agl`). The point layer holds
  2D geometry with the height as an attribute only, so in a 3D scene every
  level collapsed onto the ground and the 2 m and 10 m samples sat exactly
  on top of each other. Both layers are written - a line layer reads badly
  from directly above, a flat point layer reads badly in a scene.
- `wind.write_lines_3d` (default false), `wind.arrow_scale` (metres per
  m/s) and `wind.include_w` (tilt with vertical velocity, default false).
  `speed` and `dir_met` stay HORIZONTAL whether or not w is included, so
  the two layers stay comparable; with `include_w` the `w` value is also
  written as a field.

**Canopy volume**

- `lad_3d.nc` gained a length-1 time dimension, so it stays visible while
  a scene animates the voxel layer instead of dropping out as soon as the
  time slider moves. The timestamp is the start of the converted window
  (`bounds.time_from`, else `origin_time`), resolved without opening a 3D
  output file, so the canopy stage stays independent of the voxel stage.
- `bad` (basal area density - branches and trunks) is written as a second
  variable in the SAME file when the static driver has it, since ArcGIS
  voxel layers carry several variables. Optional and often absent, so its
  absence is silent. Adding the woody part is what makes a 3D canopy read
  as trees rather than floating crowns.
- `lad` is therefore now `(time, z, y, x)`, not `(z, y, x)` - external
  scripts indexing the file directly need the extra leading axis.

**QGIS styles**

- Vector styles gained a graduated colour ramp over the first output
  variable, spanning the whole animation (the range is accumulated as the
  timesteps are written, with no second pass). Until now GeoPackage layers
  carried only temporal settings and opened in a flat default colour.
- The polygon layer and its 2D trace layer share ONE ramp definition over
  the same variable and range, so plan view and 3D scene agree.
- A `<layer3DRenderer>` block clamps geometry to absolute height, so walls
  and wind arrows keep their Z in a scene instead of being draped on the
  terrain. New `qgis_styles_3d` (default true) turns just that block off:
  the 3D symbol schema has moved between QGIS releases more than the rest
  of the QML, and the 2D style and temporal settings are unaffected by it.

**Documentation**

- `template.yaml` now says which of the four surface products duplicate
  each other - `write_walls_gpkg` is a strict subset of
  `write_all_elements_gpkg` - so the redundant ones can be turned off
  deliberately. No defaults changed.

## 0.15.0

**Usability**

- Rewrote `template.yaml` in the palmmeteo style: every option is listed
  but COMMENTED OUT at its default, so what is uncommented is exactly what
  a project has to decide. 314 lines down to 176. The long explanations
  move to a separate guide; the template keeps only the file map, the time
  window forms, and the traps that cost a run (the quoted clock time, the
  hPa/Pa split on surface_pressure).

**Correctness**

- A config section whose every option is commented out parses as None -
  `bounds:` with nothing under it - and `_deep_merge` took that literally,
  replacing the whole section with null and failing downstream on
  NoneType. It now means "leave this section at its defaults", which is
  the only sensible reading and the one the new template relies on.

## 0.14.0

**Usability**

- Made it clear which PALM file each setting wants. `file_from` appears in
  four stages that read THREE different products - the 3D output for voxel
  and wind, the surface output for surfaces, the static driver for canopy -
  and one key name gave no hint of that. The template opens with a table
  mapping every setting to its file, each `file_from` says which product it
  needs and how to recognise it, and the README carries the same table.
- Stages now verify the product at start, by DIMENSIONS rather than by
  filename, so a renamed or joined file is judged on what it holds. Pointing
  surfaces at the 3D output used to fail deep inside the stage with a
  missing-variable error; it now stops immediately:

      [surfaces] this looks like the 3D output file, but the surfaces stage
      needs the surface output (SURFACE_DATA_NETCDF*) - unstructured surface
      elements listed along the s dimension.
        Given: OUTPUT/DATA_3D_NETCDF_N02
        Set surfaces.file_from to your surface output file. ...

  An unrecognised file is passed through to the stage as before, so nothing
  unusual is blocked.

## 0.13.1

**Correctness**

- `bounds.time_to` is now genuinely inclusive. PALM writes output at the
  first model timestep AT OR AFTER the target, so record times carry a
  residue - the prague_center files hold 600.73, 1200.73, ... A record
  nominally at 16:00 is really at 16:00:00.73, and the exact comparison
  dropped it, costing the last record of every window in every stage.
  Records are matched with a tolerance of half the file's own spacing
  (capped at 60 s), so the residue cannot cost a record while the next one
  stays out.

**Diagnostics**

- The voxel ta line reported `+nan K` for the lowest level whenever that
  level lay entirely inside terrain, which is normal for a coarse domain
  where the first output levels sit below the terrain across the whole
  cut. The offset is now reported from the gap-filled profile the
  integration actually used, and the buried levels are counted and named
  with the height of the lowest level holding air. Nothing was wrong with
  the data: `exner_profile` already interpolated those gaps.
- An all-terrain z-cut now says so instead of reporting that no valid
  potential temperature was found.

## 0.13.0

**Correctness**

- `surface_pressure` is read from the `_p3d` namelist instead of defaulting
  to a hardcoded 101325 Pa. That constant is the standard atmosphere at SEA
  LEVEL; the Prague domains sit near 200 m, some 20 hPa lower, so every
  derived `ta` / `t_surf` came out about 1.7 K too warm - a uniform offset
  that nothing in the output could reveal. PALM gives the key in
  HECTOPASCALS and it is converted to Pa here. A value set in the config is
  in Pa and still wins; the standard-atmosphere fallback applies only when
  neither source has it, and now says so. A namelist value above 2000 is
  taken as Pa with a warning, since PALM itself would have read it as hPa.
- `surface_pressure` joined the state hash for the voxel and surfaces
  stages, so changing it invalidates the products carrying derived
  temperatures. Previously a rerun kept them.

## 0.12.1

**Correctness**

- `exner_profile` anchored the Exner function at the LOWEST OUTPUT LEVEL
  rather than at z = 0, where `surface_pressure` actually applies. A
  bounds z-cut therefore shifted the whole derived-ta column by
  g/cp * z[0]: 0.04 K for a cut starting at 4 m, but 1 K for one starting
  at 100 m. The profile now descends from z = 0 to the first output level.

**Diagnostics**

- The voxel ta derivation says what it does to the numbers instead of
  reporting the Exner function:

      before: [voxel] Exner 1.00377..1.00227 over z 4..48 m
              (p_s = 101325 Pa): ta - theta = +1.11..+0.67 K
      after:  [voxel] deriving ta from theta: ta comes out 1.08 K warmer
              at z = 4 m, 0.65 K warmer at 48 m.

  The Exner values, the decay rate and p_surface moved to debug.
- Warn when the theta -> ta offset does not decay at the dry adiabatic
  rate. theta * dPi/dz = -g/cp holds whatever the stratification, so a
  departure means `surface_pressure` or the height axis is wrong - the
  kind of error that otherwise produces plausible-looking temperatures.
- `test.py` gained a check that the Exner profile is anchored at z = 0.

## 0.12.0

**Correctness**

- The timestep window can now be given in TIME
  (`bounds.time_from` / `bounds.time_to`) instead of record indices.
  `skip_times` / `n_times` count records, and a record is not the same
  amount of time in every file: PALM keeps separate output intervals for
  3D and surface output, so a run whose 3D file holds 10-minute records
  and whose surface file holds hourly ones gave

      [surfaces] records 18..21 of 24 selected     -> 18-21 h into the run
      [voxel]    records 18..21 of 144 selected    ->  3-3.5 h into the run

  from one config - two stages converting completely different periods,
  with nothing in the log to say so. Each stage now resolves the same
  window against its own time axis. Both ends are inclusive; either may be
  left null. Accepted: seconds since `origin_time` (`64800`), an offset
  with a unit (`"18h"`, `"90min"`), an absolute timestamp
  (`"2023-08-24 13:00:00"`), or a quoted clock time (`"13:00"`, resolved
  to its first occurrence at or after `origin_time`, so a window crossing
  midnight works).
- `select_times` returned a COMPACTED time array when fill records sat in
  the middle of a file, while the record indices it returned alongside
  index the data arrays. Every label after such a gap was therefore taken
  from the wrong record. Times are now returned full length and aligned
  with the file's record dimension, with NaN in the fill slots, and those
  records are excluded from the selection.

**Diagnostics**

- The selection line now reports the period, not just the indices:
  `[voxel] records 107..125 of 144 selected (2023-08-24 13:00:00 .. 2023-08-24 16:00:00).`
- When two stages are selected by record index and land on windows that do
  not overlap at all, the run warns and points at `time_from`/`time_to` -
  the case above would have been caught immediately.
- An unquoted clock time (`time_from: 13:00`) is rejected at config load.
  YAML 1.1 resolves it as the sexagesimal number 780, which would have
  been read as an offset in seconds; nothing downstream could tell it from
  a genuine number.
- A window that falls outside a file now says what the file actually
  covers instead of only reporting that nothing was selected.
- Malformed `time_from`/`time_to` values fail at config load rather than
  part-way through the first stage.

**Housekeeping**

- Removed stale duplicates of `config.py`, `pipeline.py` and `state.py`
  that were sitting unused (and unimportable) under `palm2gis/steps/`.

## 0.11.1

**Correctness**

- The `_p3d` parser read nothing from lines holding more than one
  assignment. PALM namelists routinely pack them, e.g.

      nx = 511, ny = 607, nz = 512,
      dx = 1.0, dy = 1.0, dz = 1.0,

  and 0.7.0 anchored the pattern to whole lines when it was generalised
  from six fixed keys to every setting - so those two lines matched
  nothing and a run with a valid namelist failed with "nz/dz not found in
  p3d and not set in config". Assignments are now matched wherever they
  occur.
- Comments are stripped before parsing, so the descriptive header PALM
  writes (`!-- Grid: ... (nx=511, ny=607, nz=512)`) cannot override the
  real settings, and the first occurrence of a key wins.

## 0.11.0

**Correctness**

- `Config.atomic_output()` was called by the voxel stage but never
  existed, so `stage voxel failed: 'Config' object has no attribute
  'atomic_output'` on every run - the stage could not complete at all, and
  `test.py` failed at the first check. Implemented as intended: the output
  is written to a `.partial` file in the destination directory and moved
  into place only on success, so an interrupted or failing run leaves no
  truncated file behind for the existence-based skip logic to keep.
- `--version` required `-c/--config` to be given as well, so it could
  never be used on its own. It is now argparse's own version action and
  exits immediately.

**Features**

- New `canopy` stage, the last part of `paraview_vis` that had no
  counterpart here. PALM stores `lad` on `zlad` - height ABOVE THE LOCAL
  TERRAIN - and the stage serves it three ways:
  - `<prefix>lai.tif`, LAD integrated over the column (leaf area index);
  - `<prefix>canopy_height.tif`, the highest level above `lad_threshold`
    in metres above the local terrain. Derived from LAD itself, so unlike
    a buildings-based height it works on static drivers that leave
    `buildings_2d` empty;
  - `<prefix>lad_3d.nc` (optional), LAD lifted onto absolute heights - the
    ArcGIS voxel / QGIS mesh equivalent of what paraview_vis rendered in
    ParaView. Respects the `bounds` z-cut, which the original ignored.
- The original `tree_gen.py` used `np.int`, removed in NumPy 1.24, so it
  raises AttributeError on any current install; the per-column Python loop
  over the whole domain is replaced by a vectorised lift.

**Tests**

- The synthetic dataset gains a canopy block on sloping terrain, and
  checks that LAI equals LAD x depth, that the canopy height matches, that
  the volume holds every canopy cell, and that the terrain lift spreads
  one canopy over several absolute levels - which is the whole point of
  the lift and what a naive copy would get wrong.

## 0.10.0

**Logging**

- The startup version line is no longer printed twice on the terminal when
  the default log file is attached (it is emitted to the file handler alone
  via `emit_to_file_only()`).

**Correctness**

- Wall geometry and wall labels now agree. `_element_quad` decided a wall's
  axis with `az in (0.0, 180.0)` (exact float equality) while
  `_orientation` used `round(az, 1)`, so an azimuth of 179.99999 was
  labelled `wall_S` but given an x-normal quad - a polygon rotated 90
  degrees from its own label, in the layer used to attribute wall
  temperatures. Both now go through one tolerance-based `_wall_facing()`.

**Breaking**

- `bounds.n_times` is a COUNT of records, not the exclusive END index. It
  was passed straight to `range(skip_times, n_times)`, so `skip_times: 10`
  with `n_times: 5` selected nothing at all and the error message had to
  explain that n_times must exceed skip_times. The name said count, the
  code meant end index. Behaviour is unchanged when `skip_times` is 0 or
  `n_times` is left at its default; it differs only in the case that used
  to misbehave. Set `bounds.time_end` for the old meaning explicitly.

**Performance and memory**

- GeoPackages are written one timestep at a time (first `mode="w"`, then
  appended) instead of building every frame in a list and concatenating.
  Peak memory was elements x timesteps, so the "this will hold N rows"
  warning only fired after the memory had already been committed; it is
  now one timestep regardless of how many are selected. Applies to the
  surface element, wall and wind layers.

- `_topview_raster` selects the topmost element per cell with a
  sort-and-take-last instead of a Python loop over every upward surface
  element - millions of iterations on a real domain, and the dominant cost
  of the stage. Verified to give identical results, including ties.

**Features**

- Config-change detection. The skip used to be by file existence alone, so
  editing `bounds`, `celsius`, a variable list or the wind heights and
  rerunning silently kept the old outputs - the README had to warn about
  it and tell you to pass `--overwrite`. A `.palm2gis_state.json` sidecar
  in the output directory now records which config produced each file, and
  a file whose recorded hash differs is rewritten. Only the settings a
  stage actually reads go into its hash, so a wind edit does not
  invalidate the surface rasters. A file with no record is left alone, so
  an existing output directory is not wiped on the first run after
  upgrading. `--ignore-state` restores the old behaviour.

- Provenance on every output. A `<name>.palm2gis.json` sidecar records the
  version, config path, config hash, project, stage, source file, bounds
  and command line. GeoTIFFs also carry it in their own metadata tags
  (visible in `gdalinfo` and QGIS layer properties, and it survives being
  copied elsewhere); NetCDF outputs get global attributes plus a CF-style
  `history` line appended to any existing history. Only files this run
  actually wrote are stamped, so a rerun still leaves untouched outputs
  untouched.

- `--log-file PATH` appends the full log to a file, uncoloured and always
  dated. Defaults to `<output_dir>/palm2gis.log`; `--log-file ''` opts out.

- `--version` prints the version and exits. The version is now logged in
  the startup banner instead of only at DEBUG level.

- Added `LICENSE` (MIT).

**Tests**

- `test.py` gains checks for config-change detection (a wind edit rewrites
  the wind layer, leaves the voxel output alone, and the rewritten file
  reflects the new setting; `--ignore-state` suppresses it) and exact
  GeoPackage row counts, which is what would catch a botched incremental
  append duplicating or dropping a timestep.

## 0.9.0

**Features**

- New top-level option `celsius` (default `false`): temperature variables
  (`theta*`, `tsurf*`/`t_surf*`, `ta*`) are written in degrees Celsius in
  every product - voxel NetCDF, top-view GeoTIFF and mesh NetCDF, and the
  GeoPackage element layers. Conversion is applied at write time only: the
  Exner/`ta` derivation still runs in kelvin, NoData cells keep their
  -9999 fill value, and the NetCDF `units` attributes are set to `degC`.
- Raw-input guard: each stage now rejects, at start and with an
  explanation, an input file that appears to be `palm_postproc` `coord`
  output (absolute projected coordinates or Celsius temperatures).
  palm2gis and palm_postproc are parallel branches after `join`, not a
  chain - feeding `coord` output in would apply `origin_x`/`origin_y` a
  second time and break the kelvin-based `ta` derivation. Detection uses
  three independent signals: `standard_name = projection_x_coordinate` on
  `x`, an `x[0]` far beyond the local domain extent, or a temperature
  variable already in degrees C.

## 0.8.0

**Naming (breaking)**

- The derived surface temperature is `t_surf`, not `ta_surface`, matching
  PALM's own naming: `t_surf` is the actual surface (skin) temperature,
  `theta_surface` / `theta_surf` the surface POTENTIAL temperature. Both
  spellings of the potential form are accepted, and a listed variable
  missing under one spelling is substituted from the other with a warning.
- The derived volume variable is `ta`, not `ta_derived`. A derived field
  carries `derived_from = "theta"` and an explanatory `comment`, so it is
  never mistaken for the model's own output.

**Features**

- `derive_ta` now ENSURES the actual temperature rather than always
  computing it: if the source file already holds `ta` (volume) or `t_surf`
  (surface), PALM's own - exact - value is used, and the reconstruction
  runs only when it is absent. The run says which happened.
- The surface stage reads the file header before choosing output names, so
  the decision is made against what the file actually contains.

## 0.7.1

**Correctness**

- `surfaces.derive_ta` produced `ta_surface` only in the GeoPackage
  element layers; the top-view GeoTIFF and the mesh NetCDF still held the
  potential temperature alone. All three products now share one variable
  list (`output_vars`), so a derived field can no longer be present in one
  and missing from another. `<prefix>ta_surface_topview.tif` appears
  alongside the theta raster, and `ta_surface` alongside `theta_surface`
  in the mesh NetCDF.
- Derived values are computed only where the source is valid, so NoData
  cells stay NoData instead of being scaled.

## 0.7.0

**Features**

- Time parameters are read from the files themselves. When a file's `time`
  variable is unusable, the output interval is taken from the `_p3d`
  namelist - `dt_dosurf_av` / `dt_dosurf` for surface output,
  `dt_data_output_av` / `dt_data_output` for 3D, each with its matching
  `skip_time_*` - and the first record is placed one interval in. PALM
  keeps these settings SEPARATE per product, so a run with a 10-minute 3D
  file and an hourly surface file is now handled correctly with no manual
  configuration. The averaged keys are preferred for files whose name
  marks them as averaged. The run reports which key it used.
- `bounds.time_step` / `bounds.time_start` become overrides rather than
  requirements, needed only when there is no p3d or the namelist is wrong.
- The p3d parser now reads every scalar numeric setting in the namelist,
  not just nx/ny/nz/dx/dy/dz.

## 0.6.0

**Features**

- New `palm2gis.steps.thermo`, and `voxel.derive_ta` / `surfaces.derive_ta`:
  actual air temperature from PALM's potential temperature, `T = theta *
  Pi`. The Exner function is integrated hydrostatically up the column
  (`dPi/dz = -g/(cp*theta)`, exact for dry air) from the new
  `surface_pressure` setting, using the horizontally averaged theta
  profile from the data. The voxel stage adds a `ta_derived` variable; the
  surface element layers gain a `ta_*` column beside every `theta_*` one,
  evaluated at each face's own height.
- PALM's `p` output is the PERTURBATION pressure (order tens of Pa) and is
  deliberately NOT used as absolute pressure; its true effect on T is
  about 0.02 K.
- The self-test checks the thermodynamics analytically: with theta
  constant the derived temperature must fall at exactly g/cp
  (9.761 K/km), and the Exner function must be 1 at the reference
  pressure.

## 0.5.0

**Features**

- `surfaces.write_topview_netcdf` (default true): the top-view fields are
  also written as a CF NetCDF, `<prefix>surface_2d.nc`, which QGIS opens
  as a MESH layer. MDAL reads the CF time axis natively, so the Temporal
  Controller works with no style sidecar and no band-to-time mapping -
  the reliable route for animation in QGIS, where the GeoTIFF per-band
  temporal mode proved fragile. All variables share the file as separate
  dataset groups.

**Correctness**

- Raster QML band timestamps now carry an explicit UTC marker. QGIS writes
  per-band ranges without a timezone (unlike its other temporal modes), so
  a naive string is read back as LOCAL time while the map range is UTC -
  shifting every band by the local offset and potentially rendering
  nothing.
- The raster QML no longer contains a renderer. The pseudocolour block it
  used to write was incomplete (empty colour ramp, NaN classification
  bounds) and could leave the layer blank; QGIS's own default stretch is a
  better starting point.

## 0.4.2

**Correctness**

- New `palm2gis.proj_env`: environments that mix conda and pip can end up
  with GDAL linked against one PROJ and resolving its database from
  another ("proj.db contains DATABASE.LAYOUT.VERSION.MINOR = 4 whereas a
  number >= 6 is expected"). The only visible symptom is that every raster
  is written WITHOUT a CRS. palm2gis now probes PROJ at startup and, when
  it is broken, repoints GDAL's search path at a matching database
  (rasterio's bundled copy, or pyproj's) for the duration of the run.
  Nothing on disk or outside the process is changed, and the run reports
  both the workaround and the permanent fix.
- `test.py` performs the same check before its own assertions, so a broken
  environment no longer produces a confusing CRS failure.

## 0.4.1

**Correctness**

- The raster CRS fallback never triggered: a CRS that merely *constructs*
  was accepted, but on installations with a stale or mismatched PROJ
  database a PROJ-string CRS builds fine and still cannot be encoded into
  a GeoTIFF, so the file was written unreferenced and the fallback was
  never reached. Candidates are now validated (explicit WKT first, which
  needs no proj.db lookup), and a raster that still lands without a CRS is
  REPAIRED in place by reopening it and trying the remaining candidates.
  Only if all of them fail does the run warn.
- `test.py` reads GeoPackages with an explicit `layer=`, since the
  embedded QGIS style makes them multi-layer files (pyogrio warned about
  the ambiguity). Same applies to your own scripts - see the README.

## 0.4.0

**Features**

- `qgis_styles` (default true): every layer is written with a QGIS `.qml`
  style sidecar, and GeoPackages additionally carry the style in their own
  `layer_styles` table. QGIS applies both automatically when a layer is
  added, so no Temporal tab work is needed:
  - multiband rasters arrive with every band mapped to its own time range
    (`FixedRangePerBand`, QGIS 3.38+), derived from the band timestamps;
  - vector layers arrive with the `time` field registered and a feature
    duration equal to the median output interval, so an hourly layer stays
    visible for the whole hour instead of a single instant;
  - rasters also get a pseudocolour renderer so the first view is not grey.
  Older QGIS versions ignore the temporal block and open the layer untimed;
  ArcGIS ignores the sidecars entirely.
- The XML follows QgsRasterLayerTemporalProperties::writeXml and
  QgsVectorLayerTemporalProperties::writeXml from the QGIS sources.

## 0.3.0

**Features**

- `surfaces.write_all_elements_gpkg`: EVERY surface element of the static
  driver geometry in one GeoPackage - ground, roofs and all four wall
  orientations - as PolygonZ carrying the configured variables. The
  existing products are selective (the raster keeps only the topmost face
  per column, the walls layer only the vertical faces); this is the
  complete surface picture, meant for 3D scenes where a whole block's
  roofs and wall facings are shown at once.
- New columns on the element layers: `orient` (`up`, `down`, `wall_N`,
  `wall_E`, `wall_S`, `wall_W`) for symbology and filtering, `z` (absolute
  face height) and `z_agl` (height above that column's terrain), which
  separates roofs from ground even when `buildings_2d` is empty.
- `surfaces.layer_name` overrides the layer name inside the GeoPackage.
- The run reports the orientation breakdown of what it wrote, and warns
  when a layer would exceed two million rows.

## 0.2.3

**Correctness**

- The top-view GeoTIFF could be written WITHOUT a CRS on installations
  where GDAL silently rejects a pyproj CRS object (typically a stale
  `proj.db`). The CRS is now handed to rasterio as the literal config
  string, falling back to explicit WKT, and the file is reopened after
  writing to confirm - a missing CRS is reported instead of shipping an
  unreferenced raster.
- A static driver that describes buildings only in 3D
  (`building_id` / `building_surface_level`) leaves `buildings_2d` entirely
  at its fill value. The local surface is then bare terrain, so
  `wind.heights_agl` references the ground even on rooftops; this is now
  warned about rather than passing silently.
- `origin_z` is read from the static driver (the voxel z axis is height
  above the domain base; add `origin_z` for metres above sea level).

## 0.2.2

**Logging**

- Output paths are logged relative to `project.output_dir`, which is itself
  logged once per run - the full absolute path no longer repeats on every
  line.
- One line per product instead of two: `[voxel] wrote x_voxel.nc: 29 x 46 x
  200 x 200 (t,z,y,x), u, v, w, ta`, `[wind] wrote x_wind.gpkg: 2500 pts x 2
  heights x 29 times (stride 4)`.
- The static-driver summary is two lines instead of three, and the
  fill-value warnings are one line each.
- The record-selection line appears only when the window is actually a
  subset; selecting every record is a debug detail. The version banner
  moved to debug as well.
- Everything removed from the default output is still available under `-v`.

## 0.2.1

**Correctness**

- A time variable left at its `_FillValue` (~1e36) throughout - as some PALM
  surface files are written - is no longer passed on as seconds. Doing so
  crashed the stage with "Python int too large to convert to C int" when
  building timestamps (1e36 s is ~3e29 years).
- Such files now degrade gracefully: raster bands are labelled `record N`,
  the GeoPackages carry an empty `time` field plus a new `rec` (record
  index) column that can still drive animation and filtering, and the voxel
  time axis falls back to record indices with a warning.
- New `bounds.time_step` / `bounds.time_start` reconstruct real timestamps
  for such a file from its known output interval
  (`time_start + n * time_step`).
- `rec` is written to the walls and wind GeoPackages in all cases, so a
  layer can be stepped through by record regardless of the time situation.

## 0.2.0

**Features**

- `project.overwrite` (and the `--overwrite` flag), following the
  `palm_preproc` convention: by default an output file that already exists
  is KEPT and its product skipped, so a rerun only fills in what is
  missing - delete one GeoTIFF and rerun to regenerate just that one. With
  overwrite on, every output is regenerated and the run logs a warning.
  The check is per product, not per stage, so the `surfaces` stage refills
  a single missing raster without rewriting the others, and skips opening
  the surface file at all when there is nothing to write.
- Existing files are REMOVED before being rewritten. GDAL appends to an
  existing GeoPackage layer rather than replacing it, which would silently
  duplicate every feature on each rerun.
- `--dry-run` now reports the overwrite setting and, for every output file,
  whether it would be kept, written, or replaced.

## 0.1.2

**Correctness**

- Time selection is now per file and fill-aware. PALM pads unwritten records
  at the end of a file with the time variable's `_FillValue` (a large
  integer - the same quirk the ParaView guide warns about for the last
  timestep); trailing fill records are dropped, and a file whose time
  variable is fill/masked throughout falls back to its raw values with a
  warning instead of selecting nothing.
- Different output files of the same run hold different numbers of records
  (3D and surface output use separate intervals), so `bounds.n_times` is
  clamped per file and the selected window is logged per stage
  (`[surfaces] 3 of 3 records selected (indices 0..2)`).
- When the window really is empty, the error now reports what the file
  contains and which setting to change, instead of "check bounds.skip_times
  / n_times".

## 0.1.1

**Correctness**

- `voxel`/`wind`: destaggering now handles real PALM 3D-output files, where
  the staggered dimensions `xu`/`yv`/`zw_3d` have the **same length** as
  `x`/`y`/`zu_3d` (the outermost face is not stored). Previously, when the
  `bounds` cut reached the domain edge, the staggered slice came back one
  face short and the stage failed with a boolean-index shape mismatch
  ("size of axis is N but size of corresponding boolean axis is N+1"). The
  last available face is now reused for the edge half-cell (zero-gradient
  assumption), and the idealised one-point-longer convention keeps working
  unchanged.
- `test.py`: the synthetic 3D file now uses the real PALM same-length
  staggered-dimension convention, so the self-test exercises exactly the
  case above.

## 0.1.0

Initial release — GIS-oriented re-implementation of the `paraview_vis`
pipeline for ArcGIS Pro / QGIS targets.

**Stages**

- `voxel`: destaggered 3D output as a CF-compliant NetCDF for the ArcGIS Pro
  Voxel Layer — `(time, z, y, x)`, absolute projected coordinates, `crs`
  grid-mapping variable with an Esri PE string, CF time from `origin_time`,
  terrain/building cells as NoData, optional `drop_k0` for uniform vertical
  spacing.
- `surfaces`: topmost upward-facing element per column rasterised into one
  multiband GeoTIFF per variable (band = timestep, timestamp in the band
  description); wall elements as vertical PolygonZ quads in a time-enabled
  GeoPackage.
- `wind`: wind vector point layer (`u`, `v`, `speed`, `dir_met`) sampled at
  heights above the local surface (terrain + buildings), thinnable by
  `stride`.

**Correctness**

- Destaggering fixed relative to the original `operation_3d.py`, which
  averaged each staggered component with the **domain-edge** slice
  (`var[..., -1:]` broadcast) instead of the neighbouring face
  (`var[..., :-1]`).
- The `_p3d` parser matches `key = value` pairs with a regex instead of
  whitespace splitting, so namelist formatting variations parse correctly.
- Elements with zenith outside (0, 90) are counted and reported instead of
  being silently misclassified.

**Infrastructure**

- Configuration, logging, CLI, self-test and repository layout follow
  `palm_preproc`: deep-merged defaults + annotated `template.yaml` +
  `config/example_minimal.yaml`, the PALM-GeM-style logger
  (`progress`/info/debug, `-v`/`-q`, `--log-datetime`), `--stages` and
  `--dry-run` flags, and `test.py` building a synthetic dataset in a
  temporary directory (`--keep` to inspect it).
