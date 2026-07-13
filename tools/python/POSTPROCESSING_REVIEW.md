# uDALES Python postprocessing library — architecture review & plan

**Date:** 2026-07-13
**Reviewed at:** branch `fix/issue-322-pyvista-support` (PR #324, commit `60bfef3`)
**Trigger:** software-engineering review of the PyVista integration in `udbase`/`udvis`, widened
into a full audit of the `tools/python` postprocessing library (~23k LOC).
**Scope:** `udbase.py`, `udvis/`, `udgeom/`, `udprep/`, `tests/`, packaging.

> This document is intended to be picked up and worked on in a **separate branch** from the
> PyVista feature branch. It captures the review so the PyVista work (PR #324) does not have to
> carry the whole cleanup. Findings cite `file:line` at the reviewed commit — verify line numbers
> before editing, as they will drift.

---

## Verdict

The library is **not yet "as strong as possible", but it does not need a rewrite.** It rests on
good bones — correct domain code, documentation well above research-code norms, a behaviourally
real test suite, and two clearly-right instincts (the `sim.vis` facade; `check() → fix() →
recheck`). The weaknesses are **structural debt and inconsistency**, and the remedy is dominated by
*deletion, consolidation, and two well-chosen abstractions* rather than new code.

**Honest grade:** domain correctness and documentation are top-decile for scientific code;
packaging, structural consistency, and error discipline are the weak axes. A solid **B−** that can
reach **A** mostly by removing code and making two abstractions explicit.

---

## CORRECTION (2026-07-13) — supersedes the earlier "CI / packaging" framing

An independent review (and direct verification against the code) corrected several load-bearing
claims in the first draft of this document. The corrected facts:

- **The Python tests DO run in CI — and mostly gate PRs.** CI runs `tests/run_tests.py supported`
  (`.github/workflows/ci.yml:90`); the `supported` group in `tests/test_suites.yml` (lines 82–268)
  runs the `tools/python/tests` suites via **`unittest`**. So "the suite never gates a PR" was
  **wrong**. The real gaps are narrower and specific:
    - `supported` **manually enumerates** the test modules and **omits `test_udprep_init.py`**.
      Full discovery exposes a genuine **Python-3.9 incompatibility** in that omitted file
      (`str | None` at runtime with no `from __future__ import annotations`, `test_udprep_init.py:28`)
      — invisible precisely because CI never ran it.
    - There are **two parallel test inventories**: the enumerated `supported` group and an orphan
      `python-library` group that uses `unittest discover` (`test_suites.yml:22`) but is not included
      by `supported` or `all`. Maintaining both is how modules get silently dropped.
- **The repo standard is `unittest`, not `pytest`.** A `pytest.ini`/`conftest.py` + `pythonpath` is
  therefore the **wrong** import-root fix (it would not help notebooks, scripts, `run_tests.py`, or
  plain `unittest`). Do **not** introduce pytest as a second runner. (The earlier claim that
  `test_solar.py` / direct-shortwave tests break under plain `unittest` was also overstated — they
  import `_common`, which puts `tools/python` on `sys.path`; untidy but functional.)
- **`pyproject.toml` was rejected for the wrong reason.** Repo coupling to Fortran/data does *not*
  preclude an **editable/internal** package — `pip install -e` needs no wheel/publish. The import
  model is a real decision (below), not a foregone "no packaging".
- **Environment selection is more broken than import style**, and should be Tier-0 #1 (below).
- **Dependency classification in the first draft was wrong**: `pandas` is a hard, **undeclared**
  dependency (`udbase.py:18`; absent from `requirements.txt` and `environment.yml`); `trimesh` is
  effectively required (UDBase refuses to operate without it); `shapely` is imported eagerly by
  `udgeom`. Only `pvlib` is genuinely backend-optional; `pyvista` is optional **and undeclared**.

## Address NOW (revised)

Ordered; the first item is the true prerequisite. All use the repo's existing `unittest` +
`run_tests.py` machinery — **no new runner, no packaging churn.**

- [ ] **1. Unify environment selection.** Three conventions currently contradict each other:
      `run_tests.py` hardcodes `../.venv` (`tests/run_tests.py:16`) and prefers it over the active
      interpreter; `setup_venv.sh` creates `tools/python/.venv` (`:227`); `AGENTS.md` calls `../.venv`
      canonical (`:139`). During review the dispatcher selected a stale shared env missing `shapely`.
      Fix: prefer the **active `sys.executable`** with an explicit override flag, and align the setup
      scripts, docs, `run_tests.py`, the f2py output location, and cluster guidance on one story.
- [ ] **2. Close the CI coverage gap without duplicating inventories.** Make the merge gate use the
      **discovery** suite already defined by `python-library`, or have `supported` and `python-library`
      both `include` a shared discovery group — so new/omitted test files (like `test_udprep_init.py`)
      can't silently escape the gate. Avoid maintaining two hand-written lists.
- [x] **3. Fix the Python-3.9 incompatibility** in `test_udprep_init.py` (add
      `from __future__ import annotations`). *(Done on this branch.)*
- [ ] **4. Declare actual dependencies (classify by real import boundaries).** Add **`pandas`** (a
      *hard*, currently undeclared dep) to `requirements.txt` + `environment.yml`. `pyvista` is
      optional **and** must not pollute the default install — put it in its own
      `requirements-pyvista.txt` and add it **with the backend PR, not the foundation branch**
      (uncommenting it in `requirements.txt` would wrongly make it a default install). Also note in
      the dep docs: **Plotly is feature-optional** and **numba is method-dependent** (both currently
      declared as if mandatory); `pvlib` is backend-optional; `trimesh` is effectively required. Do
      **not** add `pytest`.
- [ ] **5. Decide the import model** (explicit choice, not "no packaging"): either (a) an
      internal/editable package (`pip install -e`, optionally a `pyproject.toml` — no wheel publish),
      or (b) a launcher-managed convention enforced by *every* entry point and test env. Note entry
      points like `tools/write_inputs.py` already self-insert the path; `setup_venv.sh` does **not**
      permanently arrange it. Pick one and apply it uniformly. (Lower urgency than 1–4.)
- [ ] **6. Delete dead/committed artefacts.** Remove `udgeom/geometry_generation_backup_pre_triangle.py`
      (1040 lines, ~90% identical to the live file, imported nowhere — git history is the backup).
      `.gitignore` **already** covers `**/__pycache__/` and `*.pyc` (`.gitignore:40-41`) — only add the
      specific `tools/python/tests/.tmp/` rule.
- [ ] **7. Decide PR #324's fate** (see "PR #324"): salvage the Plotly display-once fixes separately;
      return PyVista behind a small, tested backend seam (with `requirements-pyvista.txt`); separately
      justify the `plot_quiver` default change bundled into it.

---

## Cross-cutting themes (the "why", condensed)

1. **Fragmented environment/import conventions.** Three contradictory venv conventions (see Address
   NOW #1), an import root that each entry point/test arranges ad-hoc (`sys.path` inserts in tests, a
   dual `try/except` fallback in `udbase.py`, self-inserting launchers), and undeclared deps
   (`pandas`, `pyvista`). NB: the tests themselves *do* run in CI via the `unittest` manifest — the
   problem is *how* the environment and import root are chosen, plus a manifest that omits one file
   and duplicates the inventory. (Not a call to publish a package; the import model is an open
   decision — Address NOW #5.)
2. **Five error-handling philosophies coexist** where there should be one: `raise`,
   `warnings.warn`, `print`-banner-to-stderr, silent `except Exception: pass`, and — worst —
   `sys.exit(1)` **inside library code** (`udprep/udprep_init.py:167`), which kills the host
   interpreter from a notebook. Identical preconditions are handled differently within one file.
3. **Two god-objects.** `UDBase` spans ~8 responsibilities and re-exposes the vis facade's *private*
   helpers back onto itself. `UDGeom` pastes the same cache-invalidation block ~9× and has an
   accidental method-vs-free-function split.
4. **No rendering-backend abstraction** (this is PR #324 generalised): three backends chosen by
   `IPython`-sniffing + a `pyvista=True` boolean with parallel `*_pyvista` methods; triplicated
   axis code; `plot_fac_type` reimplements the render pipeline inline.
5. **The `Section` proxy metaprogramming is a net-negative abstraction**
   (`udprep/udprep.py:52`): `__getattr__`/`__setattr__` write-through-to-`sim` buys terseness at the
   cost of type-checker/IDE blindness, silent write-typos, and one attribute access having three
   possible homes.
6. **Heavy duplication / dead code:** the 1040-line backup; connected-components reimplemented 5+
   times (orderings differ *by design* — dedup the primitive, keep per-caller ordering); three
   facet→field converters; three DDA ray kernels; solar-state assembly triplicated. (The three ground
   generators are *active* paths, not dead — consolidate shared machinery only, with golden tests.)
7. **Numerics tangled with IO → hard to test, plus leaked handles.** `udprep_radiation` interleaves
   physics with `subprocess`/netCDF/caching/writeback, which makes the *orchestration* methods hard to
   unit-test even though the sections do have coverage (IBM unit tests, a View3D integration test at
   `tests/integration/udprep/test_view3d.py`, BCS defaults/override, SEB orchestration — the gap is
   full *shortwave/SEB behavioural* coverage, not "zero"). `UDBase._load_ncdata` opens a NetCDF
   dataset (~922) and never reliably closes it.

## Genuine strengths (keep; build around these)

`check_mesh.check()` as a structured-report backbone · `solar.py` as a pure, fully-typed, IO-free
compute core · behavioural tests (analytic flux checks, MATLAB golden references, pvlib
cross-validation) · numpydoc docstrings throughout · `udbase.py:_load_sparse_file` · the `sim.vis`
facade concept · fail-fast validation in the numeric code · load-time STL normal-winding alignment ·
mtime-keyed view-factor caching.

---

## Sequenced plan

### Tier 0 — Foundation → **this is the "Address NOW" list above.**

### Tier 1 — Settle the two big design questions

- [ ] **Replace the `Section` proxy** (`udprep/udprep.py:52-81`) with an explicit, typed `sim` store.
      *Directionally right, but NOT a small change:* sections access **25+ named `sim` members** and
      implicitly proxy many more fields — not "~15". A `Protocol` alone also does **not** prevent
      runtime typos. Sequence: define field ownership + section dataclasses first, then migrate
      incrementally. *Recommendation: replace, but staged.*
- [x] **Backend-neutral 3D scene description + `backend="plotly" | "pyvista"`** — **DONE on this
      branch.** New `udvis/scene.py` defines `Scene` (meshes/lines/points/glyphs + title, colour bar,
      axis labels, bounds) and `render_scene(scene, backend, show)` with Plotly and PyVista renderers.
      `show_geometry`, `show_geometry_outline`, and `plot_fac` now build one `Scene` and dispatch by
      `backend=`; the three parallel `*_pyvista` methods are deleted and `pyvista=True` is a
      back-compat alias. `udbase_vis.py` 2063→1574 lines. Verified: 235 tests (incl. an off-screen
      PyVista build test) — only the 2 pre-existing env failures. **Update:** all 3-D plots now
      build a `Scene` and support `backend=` — the overlays (`plot_veg`, `plot_scalar_source`,
      `plot_solid`, `plot_fluid_boundary`, `plot_fac_type`, `plot_independent_surfaces`) too; the
      backend is set once on the `UDVis`/`UDGeom`/`UDBase` constructor (`pyvista=True` alias removed).
      `Scene` grew colour-bar, legend and point-label support consumed by both renderers. A Linux CI
      job exercises the PyVista backend off-screen. *Remaining:* retire the legacy
      `_render_scene`/`_render_plotly` helpers (kept only for a couple of tests now).

### Tier 2 — Decompose & consolidate

- [ ] **Split `UDBase`** (`udbase.py`): extract stateless stats (`merge_stat`, `time_average`,
      `coarsegrain_field`) into a `stats` module; a `CaseConfig`/namelist reader; a `Grid` object; a
      `FacetData`/facet-analysis module. Keep `UDBase` a thin orchestrator.
- [ ] **Delete the facade's private-helper shims re-exported onto `UDBase`** (`udbase.py:1478-1511`);
      keep only deprecated public wrappers pointing at `sim.vis.*`. Remove the `pyvista: bool` hoisted
      onto `UDBase.plot_fac` (`udbase.py:1466`).
- [ ] **Consolidate udgeom duplication:** one `connected_face_components(mesh)` primitive (currently
      reimplemented in `udgeom._split_buildings`, `split_buildings.py`,
      `check_mesh.calculate_independent_surfaces`, `fix_mesh._face_components`, `identify_ground_faces`).
      **The differing orderings are intentional** (geographic / size / face order): the shared
      primitive should return components in a **deterministic neutral order**, and each caller applies
      its documented ID policy afterward — do not collapse the orderings. Also: one `_invalidate_cache()`
      (block pasted ~9× in `udgeom.py`); a small `_meshutil` for
      `_iter_polygons`/`_project_vertical_face`/`_copy_mesh`/weld-by-coord (each duplicated 2-3×). Pick
      and document a method-vs-free-function convention.
- [ ] **Consolidate the ground generators — carefully.** The three functions in
      `geometry_generation.py` are **active** compatibility / constrained-triangulation / fallback
      paths, *not* three dead alternatives. Share their common machinery, but do **not** reduce to one
      implementation without golden/reference validation.
- [ ] **Remove dead code inside the live `geometry_generation.py`** (`_triangulate_polygon_piece`,
      `_triangulate_simple_polygon`, `_refine_skewed_triangles` chain, `_rings_to_triangle_pslg`,
      `_ordered_rectangle_boundary_points`, `_sample_interval`, `_triangles_to_mesh`, …) — likely
      ~400+ lines.
- [ ] **Separate numerics from IO** in `udprep_radiation.py` (`calc_view_factors`,
      `run_short_wave*`) and `directshortwave.py`: a pure compute layer (arrays in/out) vs an
      IO/persistence/toolchain layer (`subprocess`, netCDF, caching, `save_param` writeback,
      `gfortran` builds). Factor the shared ray-setup/energy-reconciliation and single DDA stepping
      primitive (3 near-duplicate kernels today). Route solar-state through `solar_state` (3 copies).
- [ ] **Consolidate the three facet→field converters** (`udbase.py:1523-1715`) into one
      `np.add.at`-based scatter; the others become thin wrappers. Same for the duplicated
      facet-section parsing (`udbase.py:525-558` vs `1717-1775`).

### Tier 3 — Consistency & coverage

- [ ] **One error policy library-wide:** `raise` typed exceptions for programmer/precondition errors
      (kill `sys.exit(1)` in `udprep_init.py`), `warnings.warn` for degradation, a module `logging`
      logger for progress/diagnostics (replace `print`-banners and `_missing_plot_data`/`_warn_load`).
      Narrow every bare `except Exception` (`udprep.py:257`, `udgeom.py:533/663`,
      `udbase_vis.py` outline overlay, `write_changed_params`).
- [ ] **Fix the NetCDF handle lifetime** (`udbase.py:_load_ncdata`, opens at ~922, never reliably
      closed). Prefer **load-and-close** for requested arrays, with explicit ownership for any returned
      lazy `Dataset` — **not** an automatic persistent dataset cache (safer; avoids stale/locked
      handles).
- [ ] **Fill the *remaining* coverage holes** (much exists already — IBM/radiation/BCS/SEB have unit
      +/or integration tests; `examples/101` already drives `calc_view_factors`): add **full
      shortwave/SEB behavioural** assertions (not just orchestration), and broaden `solar.py`
      (currently one loose pvlib check) and the `udbase` loaders. Batch `save_param` writes (currently
      O(vars × filesize) — 13 full parse+rewrite cycles in `IBMSection._update_counts`).
- [ ] **Smaller cleanups:** **fix** `read_matrix` (`udbase.py:264`) — it is *not* dead (called by
      `udprep_forcing.py:120`) but *broken*: defined without `self`/`@staticmethod`, so
      `self.sim.read_matrix(path, 1)` mis-binds args. Add `@staticmethod` (or `self`) — do not delete.
      Use the `f*dump` filename attributes or delete them (`udbase.py:142-161` set but ignored);
      consider `f90nml` for namelist read/write instead of hand-rolled string surgery. Renaming
      `UDGeom.stl` → `mesh` is **low value vs. compatibility cost**; if done, add `mesh` and deprecate
      `stl` gradually rather than renaming outright.

---

## PR #324 (PyVista integration) — specific guidance

Reviewed as part of this audit. **Do not merge #324 as-is.**

- **Salvage now (own small PR):** the Plotly *display-once* fixes and their `TestDisplayOnceSemantics`
  regression tests. These are correct, well-tested, and valuable independently of PyVista — the
  `if show: fig.show(); return None` additions and the `_render_scene` no-early-display change.
- **Separately justify the bundled API change:** #324 also flips `plot_quiver`'s default from `True`
  to `False` (in `UDGeom.show`, `UDVis.show_geometry`, `show_geometry_pyvista`) and adds tests
  pinning that new default. That is an unrelated behaviour change riding inside a "backend support" PR
  — it should stand on its own rationale, not be smuggled in.
- **Defer the PyVista backend** until the backend seam exists (Tier 1 — a scene description +
  `backend=` selector, per that item's qualifications), then reintroduce it behind that seam with
  real docstrings and off-screen (`pv.Plotter(off_screen=True)`) render tests + a declared `pyvista`
  dep.
- **If kept short-term anyway**, the minimum fixes are: extract the axis-drawing helper (triplicated
  verbatim at `udbase_vis.py:1680/1840/2005`); delete the unused `z_range`
  (`1678/1838/2003`); replace the silent `except Exception: pass` around the outline overlay
  (`~1995`); correct the ImportError strings naming nonexistent methods (`1611` "show_pyvista",
  `1756` "show_outline_pyvista"); document return-type is a `pv.Plotter` when `show=False` (the
  `udgeom`/`udbase` docstrings currently claim a Plotly figure); bring the one-line `*_pyvista`
  docstrings up to the module standard.

---

## Suggested first branches of work (split — don't bundle behaviour + cleanup)

1. **`fix/python-test-environment-and-discovery`** — interpreter selection (prefer active
   `sys.executable`; align `run_tests.py`, `setup_venv.sh`, `AGENTS.md`, docs), CI gate via
   **discovery** (collapse the `supported` vs `python-library` double inventory), and the Python-3.9
   `test_udprep_init.py` fix. Behaviour/topology only.
2. **`chore/python-dependencies-and-dead-files`** — declare `pandas`; add the `tools/python/tests/.tmp/`
   ignore; delete the 1040-line backup file. Pure declaration/deletion.

Both use the repo's existing **unittest** runner — no pytest, no distributable packaging. The
Plotly-fix salvage and the PyVista backend (with `requirements-pyvista.txt`) are separate PRs.

**Sequencing adjustment:** add **characterization tests before** the `Section` and radiation
refactors (Tier 1/2), not after (don't leave coverage to Tier 3) — they're the safety net for those
high-risk changes. The **NetCDF handle fix is a concrete resource bug** and can move earlier than
Tier 3.
