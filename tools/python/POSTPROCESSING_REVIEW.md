# uDALES Python postprocessing library — open product-quality tasks

**Scope:** `tools/python`: `udbase.py`, `udvis/`, `udgeom/`, `udprep/`, tests, packaging, and developer workflow.

**Purpose:** this is an open-task backlog for turning the Python postprocessing tools into a dependable internal software product. It should not track how the review evolved or which earlier comments were corrected. Keep it focused on what still needs doing.

---

## Target outcome

The Python postprocessing tools should be usable as a stable internal product:

- Users can set up an environment, import the library, and run examples without knowing repository internals.
- CI exercises the same entry points users rely on.
- Optional features fail with actionable dependency messages.
- Library code raises catchable exceptions rather than exiting the host process.
- Public plotting APIs are stable and documented.
- Numerical kernels can be tested independently from file IO, subprocesses, and mutable case state.
- Refactors are protected by characterization tests for scientific behaviour.

---

## Open tasks

### 1. Settle the environment and import model

- Choose one virtual-environment convention and apply it consistently across setup scripts, CI, docs, examples, and local developer guidance.
- Make the test runner prefer the active interpreter (`sys.executable`) unless an explicit override is provided.
- Decide whether `tools/python` should be an internal editable package.
- If using an editable package, add the minimum package metadata needed for `pip install -e tools/python`.
- If not using packaging, enforce one launcher/path convention across every entry point.
- Remove contradictory environment guidance from docs.
- Ensure notebooks and scripts can import the library without ad-hoc `sys.path` edits.

### 2. Finish dependency-boundary cleanup

- Decide whether `plotly`, `numba`, and `pvlib` are standard-environment dependencies or true optional extras.
- If they are optional, split them into feature-specific requirements files or extras.
- Make optional-import failures actionable: error messages should name the missing feature and the install command.
- Keep heavyweight visualization dependencies separate from the base environment.
- Ensure CI has coverage for the base environment and targeted optional-feature environments.

### 3. Define one library-wide error policy

- Add a small shared exception hierarchy, for example:
  - `UDALESError`
  - `ConfigurationError`
  - `DataFormatError`
  - `DependencyError`
  - `GeometryError`
  - `RadiationError`
- Replace `sys.exit(...)` in library code with typed exceptions.
- Use `warnings.warn` only for degraded-but-continuable behaviour.
- Use module loggers for progress and diagnostics.
- Remove print-based error reporting from library internals.
- Narrow broad `except Exception` handlers.
- Document fallback behaviour where catch-and-continue is genuinely intentional.

### 4. Reduce `UDBase` responsibility

- Keep `UDBase` as the user-facing case object, but extract cohesive modules for:
  - statistics: `merge_stat`, `time_average`, `coarsegrain_field`;
  - configuration and namelist read/write;
  - grid metadata;
  - facet data and facet analysis;
  - NetCDF loading;
  - facet-to-field conversion;
  - visualization delegation.
- Deprecate private visualization helper shims exposed through `UDBase`.
- Keep public compatibility wrappers thin and documented.
- Avoid moving code into new modules without also reducing responsibility or duplication.

### 5. Clean up visualization after the backend refactor

- Remove or deprecate remaining private render shims such as `_render_scene`, `_render_plotly`, `_render_trimesh`, and outline helper wrappers.
- Audit plotting docstrings against actual return types.
- Document backend dependency requirements in user-facing docs.
- Document the `plot_quiver=False` default as an explicit API/product decision.
- Keep backend-specific behaviour behind the common `backend=` API.
- Avoid reintroducing parallel public methods for individual renderers.

### 6. Reduce `UDGeom` duplication

- Extract shared mesh utility primitives.
- Consolidate connected-face component detection into one deterministic primitive.
- Preserve caller-specific ordering policies outside the shared primitive.
- Replace repeated cache-invalidation blocks with one named method.
- Document the method-vs-free-function convention for geometry helpers.
- Consolidate shared machinery in ground-generation code without collapsing active compatibility paths prematurely.
- Remove confirmed dead geometry-generation helpers only after characterization tests protect the active paths.

### 7. Separate numerics from orchestration

- Split radiation, shortwave, and SEB computation from:
  - file reads/writes;
  - subprocess execution;
  - cache management;
  - parameter-file mutation;
  - mutable `UDBase` / `UDPrep` state.
- Extract pure compute functions with arrays/config in and arrays/results out.
- Share solar-state construction.
- Share ray setup and DDA stepping primitives.
- Share energy-reconciliation logic.
- Batch repeated namelist or parameter writes rather than rewriting whole files repeatedly.
- Add characterization tests before changing high-risk scientific paths.

### 8. Strengthen behavioural coverage

- Add shortwave/SEB behavioural assertions, not just orchestration checks.
- Broaden loader tests for `UDBase`.
- Add golden-file tests for current tree/vegetation formats.
- Add tests for optional dependency error paths.
- Add example smoke tests for documented user workflows.
- Keep tests on the existing `unittest` runner unless there is a deliberate project-wide decision to change runner.

### 9. Fix data-format drift in tree and vegetation loading

- Reconcile `UDBase._load_tree_data` with the current `trees.inp` format.
- Pin parser expectations with a golden-file test:
  - skipped header rows;
  - column count;
  - 0-based vs 1-based indexing;
  - units;
  - missing or empty file behaviour.
- Verify the fields tutorial path that calls `plot_veg()`.
- Keep legacy-format compatibility only if it is intentional and tested.

### 10. Make resource ownership explicit

- Fix NetCDF handle lifetime in `UDBase._load_ncdata`.
- Prefer load-and-close for requested arrays.
- If returning lazy datasets, document ownership and provide a clear close path or context manager.
- Avoid hidden persistent file handles.
- Make cache invalidation explicit and centralized.

### 11. Fix small concrete defects

- Fix `read_matrix`: it is called as an instance attribute but defined without `self` or `@staticmethod`.
- Decide whether the `f*dump` filename attributes are used; either use them or remove them.
- Consider `f90nml` for namelist read/write instead of hand-rolled string surgery.
- Audit docstrings for stale backend or return-type claims.
- Replace library `print(...)` diagnostics with warnings or logging as appropriate.

---

## Working rules

- Prefer small, coherent commits on this improvement branch.
- Keep each commit reviewable: one defect fix, one refactor seam, one cleanup, or one documented API decision.
- It is acceptable to make API changes while the Python tools are still beta, but make them deliberate and document the new behaviour.
- Do not mix unrelated behaviour changes with mechanical cleanup in the same commit.
- Add characterization tests before changing high-risk scientific paths.
- Delete code when consolidation makes it redundant.
- Avoid adding a second abstraction beside the old one indefinitely.
- Preserve compatibility only where it is genuinely useful; otherwise prefer the cleaner beta API and update tests/docs accordingly.
- Keep optional features optional. If a module imports an optional dependency eagerly, it is no longer optional.
- Prefer small pure functions under large orchestrators.
