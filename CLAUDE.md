# [CLAUDE.md](http://CLAUDE.md)

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Commands

### Environment setup

```bash
uv sync                      # install/update the environment from uv.lock
uv run python -m build
uv pip install -e .          # editable install of the `oge` package
```

(`conda env create -f environment.yml` + `pip install --editable .` is the conda alternative — see README.md for full steps.)

### Lint / format (CI runs both via `.github/workflows/lint.yml`)

```bash
ruff check .                 # lint
ruff format .                # format (Black-compatible, line-length 88, double quotes)
```

`.flake8` also exists but `ruff` is the enforced tool in CI.

### Tests

```bash
pytest                       # run the full suite from repo root
pytest test/test_emissions.py                 # single file
pytest test/test_emissions.py::test_co2e_AR6_100   # single test
```

Note: some files under `test/` (e.g. `test_eia.py`, `test_logging.py`) still import from a legacy `src.*` module layout rather than the current `oge.*` layout and may not run as-is — check imports before trusting a green/red result from those specific files.

### Running the pipeline

```bash
cd src/oge
python data_pipeline.py --year 2024 [--flat] [--skip_outputs]
```

Set `PUDL_DATA_STORE=s3` (or `1`) to load PUDL tables from S3 in memory instead of downloading locally (`local`/`2`/unset = local disk). `--flat` uses flat hourly profiles instead of imputed ones; `--skip_outputs` skips writing CSVs for faster iteration.

## Architecture



### Pipeline flow

`src/oge/data_pipeline.py` is the single entry point and orchestrates every other module in sequence for a given `--year`. It is a straight-line script (not a DAG framework): each stage calls into one module, and later stages consume the DataFrames produced by earlier ones. Reading `main()` top-to-bottom is the fastest way to understand how the modules below fit together:

1. `download_data` — fetches source data (EIA-930, EPA CEMS, etc.) into `$HOME/open_grid_emissions_data/downloads`.
2. `load_data` — reads downloaded files and PUDL tables into DataFrames.
3. `data_cleaning` — cleans/standardizes raw EIA/CEMS data.
4. `subplant_identification` / `gross_to_net_generation` — groups generators into subplants and derives gross-to-net generation conversion factors.
5. `eia930` — cleans and shapes EIA-930 balancing-authority-level data.
6. `impute_hourly_profiles` / `anomaly_screening` — assigns hourly shapes to monthly data and flags/screens timeseries anomalies (per Ruggles et al. 2020 methodology).
7. `emissions` — calculates emissions (CO2/CH4/N2O/NOx/SO2, CO2e) from fuel consumption and emission factors.
8. `consumed` — consumption-based (as opposed to production-based) emissions calculations.
9. `validation` / `column_checks` — validate output correctness and enforce expected column names/dtypes on every output table.
10. `output_data` — writes intermediate outputs to `$HOME/open_grid_emissions_data/outputs` and final results to `.../results`.

Cross-cutting modules used throughout: `helpers` (shared utilities, e.g. plant/BA tables, geocoding), `constants` (conversion factors, year bounds like `earliest_data_year`/`latest_validated_year`/`current_early_release_year`), `filepaths` (all on-disk path resolution), `logging_util` (shared logger setup).

### Data sources & external dependencies

- **PUDL** (`catalystcoop-pudl`, installed from the `singularity-energy/pudl` fork) is the primary source for EIA/EPA data, exposed as a SQLite DB (`PUDL_ENGINE` in `src/oge/__init__.py`) or loaded from S3 depending on `PUDL_DATA_STORE`.
- `gridemissions` (Singularity Energy fork) supplies EIA-930 grid balancing data.
- The `OGE_DATA_STORE` env var (`s3` vs `local`) controls whether *OGE's own* published outputs are fetched from S3 rather than regenerated — this is for consumers importing `oge` as a package, distinct from `PUDL_DATA_STORE`.
- `constants.current_early_release_year` / `latest_validated_year` gate whether a given `--year` runs against stable vs. early-release EIA data (see `validation.validate_year`); early-release runs additionally require `PUDL_BUILD=nightly`.



### Year-bound constants (`constants.py`)

A handful of module-level ints gate pipeline behavior for a given `--year` and are updated once per release cycle (tracked with `# TODO:` comments naming the next update date): `earliest_data_year`/`earliest_validated_year` (2005 — earlier PUDL data exists but isn't considered clean), `earliest_hourly_data_year` (2019 — earliest year hourly profiles can be produced), `latest_validated_year`, and `current_early_release_year`. `ConversionFactors` is a `float` subclass used purely as a namespace of class attributes (`ConversionFactors.lb_to_kg`, etc.) — it's never instantiated, just referenced for its constants.

### Data store resolution (`filepaths.py`)

All path logic funnels through `get_oge_data_store()`/`get_pudl_data_store()`, which branch on the `OGE_DATA_STORE`/`PUDL_DATA_STORE`/`PUDL_BUILD` env vars and return either a local `~/open_grid_emissions_data/...` path or an `s3://...` URI — callers (`downloads_folder`, `outputs_folder`, `results_folder`, `pudl_folder`, etc.) never branch on these env vars themselves, they just build on top of these two functions. Because pandas/pyarrow (via `s3fs`) accept both forms transparently, the rest of the codebase reads/writes through these helper functions without caring which backend is active. `find_downloaded_file()` exists specifically to tolerate EIA's inconsistent filenames across releases (underscore/spacing drift, optional date suffixes) by globbing and preferring the shortest match.

### Loading & the canonical schema (`load_data.py`, `data_cleaning.py`)

`load_data` is where raw source column names/units get translated into OGE's internal naming convention — e.g. CEMS's `gross_load_mw` → `gross_generation_mwh`, `heat_content_mmbtu` → `fuel_consumed_mmbtu`, `co2_mass_tons` converted to `co2_mass_lb` via `ConversionFactors`. Once data leaves `load_data`, downstream modules assume OGE's naming/units, not the source system's. `data_cleaning.clean_eia923` doesn't reimplement EIA-923 allocation logic — it calls into PUDL's own `pudl.analysis.allocate_gen_fuel` to allocate generation/fuel down to the generator-energy-source level, then layers OGE-specific filtering (e.g. dropping cancelled-status generators via `constants.cancelled_status_codes`) on top.

### Schema registry (`column_checks.py`)

This module is the single source of truth for what columns every named intermediate/output table is allowed to have, and what dtype every column name should be, project-wide:

- `COLUMNS`: a dict keyed by table name (e.g. `"eia923_allocated"`, `"cems_cleaned"`, `"hourly_profiles"`) → the set of columns that table must contain. `check_columns(df, file_name)` looks up this set and treats the two directions asymmetrically: extra columns not in the expected set only `logger.warning` (the registry is allowed to lag reality), but any *missing* expected column is a hard `ValueError` — a table can gain columns quietly but can't silently lose one.
- `get_dtypes()`: returns one merged `{column_name: dtype}` dict (`pudl_dtypes` extended/overridden by `oge_dtypes`) covering every column name used anywhere in the project — dtype is treated as a property of the *column name*, not of a specific table, so the same column name must mean the same dtype everywhere it appears. Prefers nullable dtypes (`"Int32"`, `"boolean"`, `"string"`) over numpy's non-nullable equivalents.
- `apply_dtypes(df)`: applies whatever subset of that dict matches `df.columns`, `logger.warning`s (doesn't raise) about any column missing from the registry, and handles datetime columns separately via a hardcoded `datetime_cols` list — datetime needs special-case handling for tz-aware vs. tz-naive input (`try`/`except TypeError` on `.astype("datetime64[s]")`) and re-localizes to UTC when the column name contains `"_utc"`.
- `DATA_COLUMNS`: a flat (non-table-scoped) list of the core pollutant-mass/generation column names, reused across modules (e.g. `emissions.py`, `helpers.py`) to iterate generically over "every emissions/generation data column" rather than hardcoding the list again at each use site.
- **Maintenance order**: the module's own docstring specifies the workflow for schema changes — grep the whole project for every use of a column/file name, update all of them, *then* update `column_checks.py` last, and re-run `data_pipeline.py` to regenerate outputs and re-validate. Treat `column_checks.py` as the final step of a rename/add/remove, not the first.



### Recurring pattern: fallback method hierarchies

Both `gross_to_net_generation.convert_gross_to_net_generation` and `impute_hourly_profiles.calculate_hourly_profiles` follow the same shape: try the best-quality method first, `fillna()` progressively weaker fallback methods (e.g. subplant ratio → plant ratio → subplant shift → plant shift → fleet ratio → EIA default → assumed constant), and record *which* method was used per row in a companion `_method`/`gtn_method` categorical column (numbered strings like `"1_annual_subplant_ratio"`) rather than silently overwriting. When adding a new fallback tier to either of these, follow this same numbered-tier-plus-tracking-column convention.

### Release packaging (`output_data.py`)

Beyond writing per-year CSVs, `output_data` also contains one-off release-management functions (`prepare_files_for_upload`, `zip_results_for_s3`, `zip_data_for_zenodo`) that zip/reorganize published results for upload to S3 and Zenodo. These are only run manually when cutting a new versioned data release, not as part of a normal `data_pipeline.py` run.

### On-disk layout produced by the pipeline

All generated data lives under `$HOME/open_grid_emissions_data/` (not in the repo): `downloads/` (raw source files), `outputs/` (intermediate, per-year), `results/` (final published outputs, wiped and regenerated each run unless `--skip_outputs`). Manual/static reference tables live in-repo at `src/oge/reference_tables`.

### Notebooks

`notebooks/` is organized by purpose, not by module: `explore_data`, `explore_methods`, `manual_data` (used to regenerate files in reference tables), `validation`, `visualization`, `work_in_progress` (scratch/branch-specific, not maintained). Clear all cell outputs before committing a notebook.

### Reference tables vs. `constants.py`

Manual changes to reported categorical data, emission/conversion factor tables, or ID crosswalks belong in a `.csv` file under `src/oge/reference_tables` (e.g. `default_gross_to_net_ratios.csv`, `emission_factors_for_co2_ch4_n2o.csv`, `epa_eia_crosswalk_manual.csv`, `ba_reference.csv`) — not hardcoded as a dict/list in Python. `notebooks/manual_data` is where these CSVs get regenerated/updated. `constants.py` is reserved for true code-level constants that aren't "data" in this sense: unit conversions (`ConversionFactors`), year bounds, and small fixed classification lists (`CLEAN_FUELS`, `BIOMASS_FUELS`).

## Code style

Enforced by `ruff` (see `pyproject.toml`): line length 88, double quotes, Google-convention docstrings. Beyond what's enforced, the existing codebase consistently follows these conventions — match them in new/edited code. Column/table naming generally follows [PUDL's naming conventions](https://catalystcoop-pudl.readthedocs.io/en/latest/dev/naming_conventions.html), since OGE builds directly on PUDL tables — match PUDL's naming when adding a column derived from or parallel to a PUDL field, rather than inventing a new convention.

- **Docstrings**: Google style with `Args:`/`Returns:`/`Raises:` sections and full type annotations in both the signature and the docstring, e.g. `df (pd.DataFrame): data frame to operate on.` Every public function has one.
- **Imports**: stdlib first, then third-party (pandas/numpy typically first), blank line, then local `oge` imports. Local imports mix `import oge.<module> as <module>` (for modules used many times) and `from oge.<module> import <thing>` (for specific helpers/constants); multi-name imports from `oge.constants` use a parenthesized multi-line list.
- **Logging**: every module gets `logger = get_logger(__name__)` right after imports (never a bare `logging.getLogger`).
- **Pandas style**: heavy use of chained-but-reassigned transformations (`df = df.merge(...)` rather than deep one-liners), always passing `validate=` (`"1:1"`/`"m:1"`) on merges, and `dropna=False` on `groupby` when NaN groups matter. Short lowercase comments precede each logical block explaining *why* a step is needed, not what the pandas call does.
- `merge(..., validate=...)` **is mandatory**: every merge specifies `"1:1"`, `"m:1"`, `"1:m"`, or `"m:m"` so an unexpected cardinality change (row duplication/loss) raises immediately instead of silently corrupting downstream data. Pick the tightest cardinality that's actually true for the join, not `"m:m"` as a default.
- `groupby(..., dropna=...)` **is mandatory**: always pass `dropna` explicitly rather than relying on the pandas default (`True`, which silently drops NaN-key groups). Use `dropna=False` whenever a NaN group key is meaningful (e.g. missing `plant_id_eia`/`generator_id`) so those rows aren't dropped without comment.
- **Function naming**: verb-first snake_case describing the transformation, e.g. `calculate_x`, `add_x`, `create_x`, `assign_x`, `identify_x`.
- **Prefer one flexible function over several near-duplicates**: reuse a single function across similar call sites via parameters rather than forking near-identical copies. The common shape is boolean `include_x`/`requires_x` flags that toggle behavior — e.g. `calculate_ghg_emissions_from_fuel_consumption`'s `include_co2`/`include_ch4`/`include_n2o`, `add_geothermal_emission_factors`'s `include_co2`/`include_nox`/`include_so2`, or `download_helper`'s `requires_unzip`/`requires_untar`/`requires_gzip`/`should_clean` covering every download-and-extract case in `download_data.py`. When a new case is "almost the same as an existing function," add a parameter to the existing one instead of copy-pasting a variant.
- **Errors**: `ValueError` for bad/invalid data, `UserWarning` (raised, not warned) for pipeline-level validation failures the user must act on (e.g. unsupported `--year`, missing emission factors), `FileNotFoundError`/`OSError` for filesystem/config problems. Plain `assert` shows up in adapted third-party algorithmic code (e.g. `consumed.py`) rather than in original OGE logic — don't add new bare `assert`s outside that context.
- **Tests**: `pytest`, function-scoped fixtures for sample DataFrames, one `test_<function>_<scenario>` per case, assertions checking output columns/values rather than mocking internals.
- **"Coordinating function" naming convention**: a function whose job is to orchestrate several helpers (little logic of its own) says so in the first line of its docstring — e.g. `clean_eia923` ("This is the coordinating function for..."), `calculate_hourly_profiles`, `convert_gross_to_net_generation`, the builders in `create_nox_so2_factors.py`. Use this phrase for new top-level orchestration functions so a reader knows to look at what it calls rather than its own body.
- **File layout: top-down, not bottom-up**: within a module, the main/coordinating function is defined first, followed by the helper functions it calls, in roughly the order they're invoked (e.g. `gross_to_net_generation.py` defines `convert_gross_to_net_generation` first, then `calculate_gross_to_net_conversion_factors` and `filter_gtn_conversion_factors`, which it calls, right after; same in `data_cleaning.py`, `emissions.py`, `impute_hourly_profiles.py`). Place new helper functions after the function that calls them, not before.
- **Module-level UPPER_CASE list constants for column groups**: group a column list once at module scope and reuse it for selecting/renaming/iterating, rather than retyping the same list at each call site — e.g. `GENERATED_EMISSION_RATE_COLS`/`CONSUMED_EMISSION_RATE_COLS` (`output_data.py`), `POLLUTANTS`/`ADJUSTMENTS`/`EMISSION_COLS` (`consumed.py`), `DATA_COLUMNS` (`column_checks.py`), `CLEAN_FUELS`/`BIOMASS_FUELS` (`constants.py`).
- `.copy()[...]` **slicing**: when narrowing a DataFrame to a column subset for isolated use, chain `.copy()` immediately followed by the column selection on one line (e.g. `eia923_allocated.copy()[["plant_id_eia", "plant_primary_fuel"]]`) rather than slicing then copying separately, to avoid `SettingWithCopyWarning` at the point of narrowing.
- `apply_dtypes(...)` **right before** `return`: most data-producing functions end by coercing to expected dtypes as one of the last steps before returning, rather than deferring dtype enforcement to a later stage.
- **Explicit column reordering**: larger table-building functions build a `new_column_ordering`/`column_order` list and reindex (`df = df[new_column_ordering]`) near the end, so output column order is deliberate rather than whatever order merges happened to produce.
- **ALL-CAPS section-divider comments**: a `# SECTION NAME` line followed by a full-width `#####...` line marks major sub-sections within the largest, oldest modules (`validation.py`, `eia930.py`). Not used in every file, but the convention to reach for once a module accumulates multiple distinct groups of functions.
- **Upstream workarounds**: temporary fixes for bugs in PUDL or source data are marked with a `# NOTE:` comment citing the specific upstream issue/PR link and the date the workaround was added (and, once fixed upstream, a follow-up dated comment on why the workaround was removed/kept). Follow this pattern rather than silently patching around external bugs.



## Data quality practices

New pipeline code should be checked against these dimensions — each maps to concrete checks already used in `validation.py`/`column_checks.py` and should be extended the same way for new tables/columns rather than skipped:

- **Completeness** — no unexpected missing values. On input: check each column for nulls, and for timeseries check that the full expected set of timestamps exists at the expected interval over the date range. Mid-pipeline: confirm row counts/dimensions don't change unexpectedly, intermediate helper columns get dropped before output, and (for allocation steps) that output totals reconcile against input totals. On output: reject unexpected extra columns (see `column_checks.check_columns`).
- **Validity** — data has the expected shape. Declare expected dtypes per column and apply them on output (`column_checks.apply_dtypes`/`get_dtypes`); prefer pandas' nullable dtypes (`"Int64"`, `"boolean"`) over `"int"`/`"bool"` so a stray NA doesn't silently upcast a column to `float`. Check values fall in the expected domain and that timestamps are in the expected timezone (watch for DST shifts).
- **Uniqueness** — no duplicate values where duplicates aren't expected: unique timestamps where a table should be one-row-per-timestamp, unique primary keys, and (per the pandas rules above) `validate=` on every merge to catch accidental fan-out.
- **Accuracy** — values are plausible against a source of truth: define acceptable min/max ranges per numeric column (e.g. no negative fuel consumption, context-specific bounds) and flag violations rather than passing them through silently.
- **Consistency** — values agree across time/sources: screen for global extremes and for large percent-changes step-to-step (local extremes) — this is what `anomaly_screening.py`'s Ruggles et al.-based methodology is for; extend it rather than adding a parallel ad hoc check.
- **Timeliness & versioning** — for slow-moving reference inputs (e.g. emission factors), flag when the data is older than expected for the given `--year`. Intermediate static/reference files should either always be regenerated/overwritten by the pipeline or be tied to the `oge` package version (see `filepaths.get_oge_data_store`'s `vX.X.0` versioning) so a stale file from a previous version is never silently reused.

## Documentation

`docs/docs/` contains the published methodology documentation (https://docs.singularity.energy/docs/open-grid-emissions), written in the same level of narrative detail as this file's "Recurring pattern: fallback method hierarchies" note above — e.g. `Methodology/Converting Gross to Net Generation.md` walks through the exact same numbered fallback hierarchy implemented in `gross_to_net_generation.convert_gross_to_net_generation`, including *why* each method exists and its tradeoffs, not just that it exists.

**Whenever a code change alters a methodology already documented at this level of detail** — adding/removing/reordering a fallback tier, changing a calculation approach, switching which data source or column feeds a calculation, changing a filtering/validation threshold that's called out by name in the docs, etc. — update the corresponding file(s) under `docs/docs/` as part of the same change, not as a separate follow-up. The directory layout mirrors the pipeline stages, so the right file is usually easy to locate:
- `Methodology/Data Cleaning/` — EIA-923/EIA-930/CEMS cleaning steps
- `Methodology/Data Aggregation/` — subplant aggregation, primary fuel assignment, BA aggregation
- `Methodology/Converting Gross to Net Generation.md` — gross-to-net conversion hierarchy
- `Methodology/Assigning Hourly Profiles to Monthly Data/` — hourly profile imputation/shaping
- `Methodology/Emissions Calculations/` — GHG/NOx/SO2 calculations, CHP/biomass adjustments, consumption-based emissions
- `Data Validation/` — data quality metrics, eGRID comparison, known issues
- `Overview/` — input data sources, general project description

If a change doesn't have an obvious matching doc file, say so and ask rather than silently skipping the update — it may mean a new doc file is needed, or that the change isn't actually methodological (e.g. a refactor with no behavior change, which doesn't need a doc update).

## Contribution workflow

`main` is protected by convention, not just tooling — work happens on a descriptively-named feature branch and lands via PR, never a direct commit to `main`. The PR template (`.github/pull_request_template.md`) asks for purpose, what the code does, how it was tested, where to look, and a checklist that echoes conventions already noted above: update `docs/docs/` for methodology changes, clear notebook outputs, add docstrings/type hints to new functions.

