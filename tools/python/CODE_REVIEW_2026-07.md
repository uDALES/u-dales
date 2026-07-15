# uDALES Python library — code review findings (July 2026)

**Scope:** `tools/python`: core (`udbase.py` + helpers), `udgeom/`, `udprep/`, `udvis/`, examples, docs, tests, CI, and structure.

**Method:** five parallel in-depth review passes (core, udgeom, udprep, udvis/examples/docs, tests/packaging), each claim verified against the actual code before being reported. Findings reference `file:line` at the state of branch `fix/issue-322-pyvista-support`.

**Relation to `POSTPROCESSING_REVIEW.md`:** the one task still open there (settle the environment and import model) is confirmed by this review as the top structural priority — see "Structure and packaging" below. Everything else in this document is new.

**Two maintainer decisions recorded during review:**

1. *Packaging vs Fortran:* making `tools/python` an editable package does **not** interfere with the f2py workflow. `pip install -e` registers the source dir on `sys.path` without copying or building; the `.pyd`/`.so` extensions built by `setup_venv.*` continue to be found in-place inside `udprep/`. The `pyproject.toml` declares only the pure-Python modules; the Fortran build stays in `requirements-build.txt` + setup scripts.
2. *Committed notebook outputs:* keeping executed outputs in the tutorial notebooks is deliberate — end users do not have the experiment folders needed to regenerate them. The residual concerns are only (a) avoid re-committing re-executed notebooks when content hasn't meaningfully changed (each execution rewrites the multi-MB embedded scenes and adds a full copy to git history; `.git` is already ~310 MB), and (b) optionally use static screenshots where interactivity adds little. Notebooks requiring external experiment data is by design, not a defect.

---

## Status — updated 2026-07-15

The "same day" and "this week" priority tiers were implemented on this branch (commits `d279ebf`..`9d9908d`), each fix TDD-first and independently reviewed, with a final whole-branch review. Full unit suite green (317 tests).

**Fixed (closed):**

| Finding | Fix commit | Notes |
| --- | --- | --- |
| T1 (slow-test CI gate) | `d279ebf` + follow-up | Gated direct-shortwave tests now run in the experimental suites — and pass. Follow-up (post-review): the skip is now raised at module level BEFORE importing the solver, since the numba import cost hit ordinary discovery even when the class decorator skipped the tests |
| G1 (self-destructing lazy wrappers) | `386968c` | Eager imports + 3 regression tests |
| C1, V1, C4, C7 | `79a440f` | geom=None init; plot_trees→plot_veg; veg cache keyed on zero_based; ndmin on facet loaders |
| C2, C3 | `313a705` | DataFormatError on corrupt factypes width / prof.inp; missing-file fallback preserved. The code uses the correct width `6 + 4*nfaclyrs + 1`; `cb699a2` only corrected the formula in this document's C2 text (docs-only) |
| P1, P3 | `3484b6a` | s_veg array contract (no-veg timedep no longer crashes); per-call maxD honoured; caches unified |
| P2a, P2b | `1f36b23` | sveg files now true W/m3 (single scaling owner); Möller attenuates before crediting facets. Budget closes ~1e-12; Möller == facsec. **Regenerated cases will see a step-change in sveg values — intended** |
| P6, P8, P9 | `a3f24f1` | IBM outputs reload into sim; writeback errors propagate; grid stretch-fit raises instead of NaN |
| E1, E2 | `2c9d791` | Radiation tutorial exit()/backend fixed; vegetation tutorial rewritten to current API |
| Final-review regression + P28 annotation | `9d9908d` | Veg sections set cache-provenance flag; S_veg/load_stl docstrings honest; empty-factypes message |
| S1 packaging + env/import model | `e2c5ebf`, `3c1d049` | `tools/python` is an editable package (`pyproject.toml`); one venv convention (`tools/python/.venv`); runner uses active interpreter; per-test `sys.path` hacks + internal dual-import fallbacks deleted (see POSTPROCESSING_REVIEW #1) |
| D2 | (prior branch) | `setup_venv` import check requires core `pyvista`, not the now-optional `plotly` |
| P12, P25 | `909915f` | `_parse_shell_config` catches `OSError` (Windows no-bash fallback reachable); `write_scalarsources` skips only the existing file, not the whole loop |
| P4, P5 | `c3bd49c` | facsec raises `RadiationError` on missing `Sc`; Moller `ktot` includes the mesh z-extent; `calc_short_wave` keys Sdir/netsw reuse on an input signature so a sun sweep recomputes |
| T2, T3 (partial) | `c9e9d76` | Characterization tests added for the `load_*`/analysis API and `udprep_seb` (safety net for the remaining refactors) |
| #5 viz cleanup, V9 | `98d0abe` | Legacy render shims already gone with the Scene refactor; removed orphaned test scaffolding; fixed stale plotly/trimesh backend docstrings |
| C23, C13 | `9de8584` | `udfacet.facsec_to_field` vectorised (`np.add.at`) + shared by `convert_facflx_to_field`; dead `_file_has_data` removed |
| G25, G26, G27 | `9e8f4cd` | 443 lines of dead `geometry_generation.py` code removed; byte-identical mesh utils extracted to `udgeom/_meshutil.py`; `ndarray.ptp()`→`np.ptp()` (NumPy-2.0) |

**Still open:** the two biggest structural refactors — reduce-`UDBase` further and fully separate-numerics-from-IO in radiation (POSTPROCESSING_REVIEW #4/#7; note udstats/udgrid/udfacet/udconfig/udnetcdf/_radiation_compute are already extracted, so this is now incremental) — plus the remaining unmarked minors (C-list/G-list/P-list). Deferred sub-items noted in code: `run_short_wave`/`run_short_wave_timedep` file-existence skips are not yet signature-validated (use `force=True` after namelist/time changes); the P7 SEB gate discrepancy (characterized, not resolved).

**Newly discovered during fixing:** the committed `ibm_preproc_f2py` binary can drift from the caller signature (confirmed on one Windows machine: zero-arg extension vs ~20-arg call → `run_all` cannot complete until rebuilt). The unit tests mock one level above the f2py call, so this drift is invisible to CI — worth a smoke test that imports the extension and checks its signature, and a rebuild note in `fortran/README.md` (extends S2).

---

## Verdict

Architecturally sound: clean section/orchestrator split in `udprep`, check/repair design in `udgeom`, Scene abstraction in `udvis`, extracted stateless helper modules, and an unusually strong test suite (MATLAB golden files at 1e-12, analytic physics checks, disciplined negative-path coverage). The verified problems cluster in four places: (1) correctness bugs that silently corrupt science, mostly in vegetation/radiation coupling and `udbase` loaders; (2) an error-handling philosophy of broad `except` + warnings that converts corrupt inputs into plausible-but-wrong results; (3) internal duplication that lets paired implementations drift apart; (4) packaging debt whose symptoms appear in every layer.

---

## Priority actions

*(The "same day" and "this week" tiers below are **done** — see Status above. The "structural" tier remains open.)*

**Same day (small, high impact):**
- Enable `UDALES_RUN_SLOW_TESTS` in CI (one line of YAML) — see T1.
- Fix `s_veg` shape crash (P1), `maxD` no-op (P3), `plot_trees` forwarding (V1), `geom = None` default (C1), the udgeom lazy-wrapper bug (G1), the stray `exit()` in the radiation tutorial (E1).

**This week:**
- sveg irradiance factor (P2a) and Möller veg/facet double-count (P2b), with characterization tests first.
- Replace silent fallbacks in `udbase.py` with `DataFormatError` (C2, C3, C4).
- `ndmin=2` on the `np.loadtxt` calls (C7).
- Reload `Sc`/`facsec` after the IBM step in `run_all` (P6).

**Structural (next):**
- `pyproject.toml` + `pip install -e tools/python`; delete all dual-import machinery (S1). Everything else gets easier after this.
- udgeom consolidation pass: one building definition, one outline, one ground estimator, one tolerance policy (G-theme).
- Unit tests for the `load_*`/averaging analysis API and `udprep_seb.py` before the issue-#313 SEB redesign begins (T2, T3).

---

## Core: `udbase.py` and helper modules

### Critical

- **C1. `udbase.py:197, 391, 709` — `UDBase(load_geometry=False)` produces an object that crashes on display.** `self.geom` is only assigned inside `_load_geometry()`; with `load_geometry=False` the attribute never exists. Plain `repr(sim)` and all `udvis` plot methods (`udbase_vis.py:291, 351, 427, 493`) hit unguarded `self.sim.geom` → `AttributeError`. Fix: set `self.geom = None` unconditionally in `__init__`.
- **C2. `udbase.py:465-476` — silent guess of `nfaclyrs = 3` mis-slices `factypes` columns.** If the actual file has a different layer count, `data[:, 6:9]` etc. return the wrong physical columns without error; `load_seb` then computes G = −λ·dT/dz with the wrong conductivity. Validate the width and raise on mismatch. *(Correction during fixing: the true width is `6 + 4*nfaclyrs + 1` — four per-layer blocks d, C, l, k with k having nfaclyrs+1 entries — verified against the udprep writer and example files; fixed on this branch.)*
- **C3. `udbase.py:345-354` — any parse error in `prof.inp` silently degrades a stretched grid to uniform.** Bare `except Exception` around `_load_grid` substitutes a uniform z-grid with only a warning. Every z-coordinate, `dzt`, flux normalisation and frontal-area result becomes silently wrong. Should raise `DataFormatError`.

### Major

- **C4. `udbase.py:606-609` — `load_veg` cache ignores `zero_based`**; a later `load_veg(zero_based=False)` returns cached 0-based data. Key the cache on the flag (or remove the flag).
- **C5. `udbase.py:529-531, 1295-1304, 1383-1385` — dead guards for missing `facsec['c']`.** `self.facsec` is always a dict, never `None`, so `if ... self.facsec is None` never fires and users get a bare `KeyError: 'c'` instead of the written error message. Check `if 'c' not in self.facsec`.
- **C6. `udbase.py:506-546` vs `1422-1480` — `_load_facet_sections` and `load_facsec` duplicate parsing with divergent behaviour.** Public `load_facsec('u')` on a valid empty file crashes with `IndexError`; the private loader handles it. Consolidate in `udfacet.py`.
- **C7. `udbase.py:437, 449, 463` — `np.loadtxt` without `ndmin` breaks single-row facet files.** A one-facet `facets.inp` gives shape `(4,)` → IndexError → swallowed into a warning → facet data silently "missing". Add `ndmin=2` (`ndmin=1` for areas).
- **C8. `udnetcdf.py:66-71` — browse mode loads the entire NetCDF into RAM.** `sim.load_field()` with no argument executes `.load()` on the whole fielddump (routinely tens of GB) just to list variables. Return metadata only. Related: `display_ncinfo` (udnetcdf.py:85) lists only `data_vars`, hiding the `t` coordinate `load_seb` depends on.
- **C9. `udbase.py:25-52` — try/except-ImportError dual-import pattern is fragile.** (a) Masks real ImportErrors raised *inside* a module; (b) puts generic names (`exceptions`, `udconfig`) on the flat namespace where they can shadow/be shadowed; (c) if both branches run in one process, two copies of every module exist and `except UDALESError`/`isinstance` checks fail across them. Resolved by packaging (S1).
- **C10. `udconfig.py` + `namelists.json` — no schema validation; the JSON's two maps have diverged.** `parse_namoptions` accepts any key (typos become dead attributes); `udconfig` never reads `namelists.json` (only `udprep/_section.py:146` uses half of it); nine keys (`ltdump`, `tsample`, …) appear in both the `NAMSTATSDUMP` and `OUTPUT` arrays (lines 138-147 vs 173-183) while `"variables"` maps them to `NAMSTATSDUMP` only — `save_param` risks duplicate keys across blocks. Keep one map, derive the other.
- **C11. `udbase.py:247-248` — namoptions keys are `setattr` directly onto the instance.** Any key colliding with UDBase state (`path`, `trees`, `geom`, method names…) silently clobbers it. A `self.nml` dict with compatibility `__getattr__` preserves the MATLAB-style ergonomics safely.

### Minor

- C12. `udbase.py:1-14` — module "docstring" sits after `from __future__ import annotations`, so it is a discarded expression (`ast.get_docstring` → None). Swap.
- C13. `udbase.py:55-66` — `_file_has_data` is dead code (no callers).
- C14. `udbase.py:163-175` — most "configurable filename" attributes (`fnamoptions`, `fsolid`, `ffacEB`, `ffacT`, `ffac`, `ffacets`, `ffactypes`, `ffacetarea`, `ftrees`) are never used to build paths; loaders hardcode names. Delete or honour them.
- C15. `udbase.py:187` — `_lfgeom` is set True and never updated; dead flag.
- C16. `udbase.py:157` — `if callable(path): path = path()` is undocumented, contradicts the type hint, and exists only to paper over a notebook typo (`.resolve` without parens — see E3). Fix the notebooks, delete the hack.
- C17. `udconfig.py:22-33` — `parse_value` returns strings for Fortran d-exponents (`1.5d2`) and array values (`1, 5, 9`), both legal namelist syntax; `udconfig.py:43` opens without `encoding=` (cp1252 on Windows).
- C18. `udbase.py:340` vs `752` — `_load_grid` skips 1 header row of `prof.inp`, `load_prof` skips 2; both work only because example headers start with `#`.
- C19. `udbase.py:1150-1180` — `area_average_seb` drops `Tsurf` although its docstring says "all surface energy balance terms".
- C20. `udbase.py:1127-1144` — `area_average_fac` axis auto-detection is ambiguous when `n_facets == n_times`; add a warning or explicit `axis=`.
- C21. `udbase.py:1070-1080` — `assign_prop_to_fac('name')` crashes (`factypes['name']` is a list, no `.ndim`); unmatched `typeid`s silently become NaN.
- C22. `udbase.py:1297-1299, 1544-1546` — error messages name files that don't exist (`facet_sections_(u,v,w,c).{expnr}` vs actual `facet_sections_u.txt`); comment at 1413 ("convert from 1-indexed") is false — data is already 0-based; "fixing" the code to match it would introduce an off-by-one.
- C23. `udbase.py:1347-1420` — `convert_facflx_to_field` re-implements `udfacet.facsec_to_field` inline with a pure-Python loop; `np.add.at` would be ~100× faster. Mixes float32 accumulation (1388) with float64 elsewhere.
- C24. `udbase.py:696-727` — `__repr__` dumps grid summary *plus* every scalar attribute *plus* all namoptions params (hundreds of lines in a notebook). Move the full dump to a `describe()` method.
- C25. `udstats.py:36-45` — `time_average` returns a trailing singleton window dim, undocumented; `merge_stat` (87-102) silently interprets a two-positional call as `(XpXp, n)`.
- C26. `udbase.py:287-291` — `read_matrix`: commented-out dead code; "Unsupported file type" raised for any `loadtxt` failure.
- C27. Plot wrapper `show` defaults are inconsistent: `plot_veg`/`plot_scalar_source`/`plot_trees` default False; `plot_fac`/`plot_fac_type`/`plot_building_ids`/`plot_2dmap` default True.
- C28. `udbase.py:989-1004` — `load_seb` opens/re-transposes `facEB.nc` five times, `facT.nc` twice; only the `t` load is wrapped for FileNotFoundError.
- C29. `udbase.py:1422` — `load_facsec(var)` uses `var` for a grid designator (`'u'|'v'|'w'|'c'`) and doesn't validate it (unlike `load_slice`'s `plane`).

**Good:** the extraction of `udstats`/`udgrid`/`udfacet`/`udconfig`/`udnetcdf` out of the god-class is the right architecture — genuinely stateless and tested; `udnetcdf.load_ncdata` closes handles correctly in both modes; the exception hierarchy is clean (`DependencyError(UDALESError, ImportError)` is a thoughtful compatibility move); `_load_tree_data` is a model loader.

---

## `udgeom` package

### Critical

- **G1. `udgeom/__init__.py:56-65` — the lazy `truncate_below_ground`/`extrude_to_ground` wrappers self-destruct.** Running the wrapper body's `from .truncate_below_ground import ...` rebinds the package attribute to the *submodule*, replacing the wrapper (verified empirically): second call raises `TypeError: 'module' object is not callable`; one `extrude_to_ground` call clobbers `truncate_below_ground` too. Import the submodule under a different name, or import eagerly like every other symbol.
- **G2. `geometry_generation.py:55-100, 1550` — `_divide_faces` is 4-way, MATLAB's `divideFaces` is 2-way, but `divisions = round(Hx/edgelength)` was kept.** Face count grows as 4^d: `Hx=8, edgelength=0.5` → ~4.3×10⁹ triangles per cube (OOM); even `Hx=4, edgelength=1` yields 16× more faces than requested at the wrong facet size. Use `ceil(log2(Hx/edgelength))` and document the divergence.

### Major

- **G3. `truncate_below_ground.py:194-206` — unconstrained Delaunay on concave clip pieces.** `shapely.ops.triangulate` ignores concave boundaries: ground can develop holes near the footprint or spill into the building base. Fires routinely for L-shaped/notched footprints; also reused by `extrude_to_ground.py:149`. The package owns the correct tool (constrained triangulation via `triangle` with domain validation) — use it here, or add a post-hoc coverage check.
- **G4. `udgeom.py:514-515` vs `split_buildings.py` — two inconsistent building splitters.** `get_buildings()` hard-codes ground at z>0 and uses manifold-only `face_adjacency`; `split_buildings()` is topography-aware and non-manifold-tolerant. They disagree on terrain meshes and towers meeting at non-manifold edges; their building IDs cannot be reconciled. Share one definition.
- **G5. `delete_ground.py:56` — `identify_ground_faces` always nominates something as ground**, so `delete_ground` on a building-only STL silently deletes the lowest roof/slab. Add a confidence gate (minimum area fraction, or require the component to underlie other geometry).
- **G6. `fix_mesh.py:440-449` — duplicate-face removal keeps arbitrary winding.** Back-to-back internal walls collapse to a single wall oriented by STL file order; half can end up inward-facing, which the IBM misreads, and `check` doesn't re-validate winding. Delete the pair (interior wall) or resolve orientation against the enclosing volume.
- **G7. `check_mesh.py:206, 281-283` / `fix_mesh.py:181` — tolerance binning by `int(round(v/tol))` misses defects at bin boundaries and breaks at UTM-scale coordinates** (the `1e-8`-binned line *moment* exceeds bin width for coordinates ~10⁵–10⁷ m, so `weld_touching_boundaries` quietly does nothing). Use neighbour-bin lookup or sort-and-cluster.
- **G8. `udgeom.py:996-1213` — inconsistent in-place/return conventions:** `fix`/`resolve_vertical_coplanar_overlaps`/`weld_touching_boundaries` return `(Trimesh, report)` while `repair_adjacent_buildings`/`truncate_below_ground`/`extrude_to_ground` return `(self, report)` and expose `return_trimesh`. `UDGeom.fix` also silently drops the `weld_touching_boundaries` option that `fix_mesh.fix` supports. One convention everywhere.
- **G9. `geometry_generation.py:880-914` — `_matlab_domain_boundary` emits a dangling constraint when `nedges_y == 1`**, crashing ground generation with an inscrutable `KeyError` for domains one edge-length deep.
- **G10. `__init__.py:12-29` + `check_mesh.py:32` + `fix_mesh.py:39-41` — the soft-degradation design is dead.** `trimesh` is carefully guarded (`TRIMESH_AVAILABLE`) but `__init__` eagerly imports modules with unguarded `shapely`/`triangle` imports, so `import udgeom` fails hard anyway and the guards are unreachable. Declare hard dependencies and delete the guards, or guard consistently.
- **G11. Silent fallbacks that change results:** `calculate_outline2d` maps any internal bug to "empty polygon" (`udgeom.py:676-678`); `_generate_ground` silently falls back to the MATLAB-style generator (`geometry_generation.py:1329-1334`); `existing.contains(...)` wrapped in `except Exception: pass` (1357-1362); failed submesh extraction drops the building with only a warning, leaving its faces labelled "ground" (`udgeom.py:546-547`, `split_buildings.py:150-152`).

### Minor

- G12. `udgeom.py:209-210, 300-301` — `load`/`save` convert every exception to `ValueError` without `from e`; the docstring example `geom.load('geometry.001')` itself fails (no loader for `.001`).
- G13. `udgeom.py:195-199` — Scene→mesh concatenation ignores scene-graph transforms.
- G14. `udgeom.py:29-33` — `ConvexHull`/`SCIPY_AVAILABLE` imported, never used; `calculate_outline2d` docstring claims "convex hull" but returns the free-boundary loop.
- G15. Dead code in `udgeom.py`: unused `split_buildings` import (56); redundant `_outline3d is None` re-check (723-724); unused counters (1309-1319).
- G16. `udgeom.py:909-913` — `add_ground(inplace=True, return_debug=True)` resets cache fields by hand instead of `_invalidate_cache()`.
- G17. Inconsistent empty sentinels: `bounds` returns a fake zero box; `face_centers`/`face_normals`/`face_areas` return shape `(0,)` not `(0,3)`; `face_incenters` returns `(0,3)`.
- G18. `calculate_outline.py:141` — empty result has shape `(0,)`, not the documented `(K, 2)`.
- G19. `truncate_below_ground.py:239-242, 327` — `edgelength` parameter accepted, documented, then ignored; report field contains the *estimated* spacing instead.
- G20. `fix_mesh.py:497` — `merge_tolerance` → digits conversion rounds the wrong way (0.005 → 1e-3, *less* merging than requested).
- G21. `fix_mesh.py:319` — T-junction weld skips faces with junctions on two edges and is not iterated to a fixed point; three different effective tolerances in one routine (12-decimal keying, 8-digit merge, 1e-8 point tolerance).
- G22. `fix_mesh.py:503-505` — `removed_unused_vertices` measures only the final pass; nearly always 0.
- G23. `check_mesh.py:707-712` — degenerate-face mask computed twice; `fix` runs full `check` twice; for 10⁶-face meshes this is minutes of pure-Python looping. `n_connected_components`/`n_independent_surfaces` are duplicate keys.
- G24. `check_mesh.py:752-778` — below-ground detection suppressed entirely when a building is flush with the domain edge (`extends_beyond_non_ground` gate), undocumented.
- G25. `geometry_generation.py:1296-1298` — SciPy-less fallback calls `ndarray.ptp()`, removed in NumPy 2.0 (verified crash under 2.3.5); the advertised fallback is dead code.
- G26. ~500 lines (~28%) of `geometry_generation.py` are verified dead: `_triangles_to_mesh`, `_signed_area_2d`/`_point_in_triangle`/`_triangulate_simple_polygon`, `_sample_linestring_points`, `_refine_skewed_triangles`, the `_polygon_to_triangle_pslg` chain, `_rings_to_triangle_pslg`, `_ordered_rectangle_boundary_points`, `_sample_interval`, unused `reduce` import, discarded `footprint_polys` in `create_canyons`, the `tol = 0.0` block in `create_cubes`.
- G27. Duplication across check/fix: `_iter_polygons` and `_project_vertical_face` verbatim copies; plane-grouping loop copied; `_copy_mesh` duplicated; `_estimate_planar_ground_level` duplicated with raise-vs-None semantics; edge→face-count map hand-rolled ≥5 times; "clipped grid lines" debug block ×3; `UDGeom._calculate_outline_edges` is a third outline implementation kept alive only for udvis.
- G28. Magic absolute tolerances throughout (1e-12 areas, 1e-8 planes, 1e-6 elevations, `1.35 * target_spacing`, digits 8 vs 12), none derived from mesh scale — behaviour shifts for mm- or km-unit models. Adopt one tolerance policy (scaled by bbox diagonal or exposed constants).
- G29. `view3d.py` not exported in `__init__.py`; `resolve_view3d_exe` returns a non-existent path silently; `read_view3d_output` outformat 0 unconditionally drops the last row; `write_svf` uses `fmt="%4f"` (width-4, not 4-decimals).
- G30. Inherited MATLAB quirks worth documenting: `create_canyons` requires `ysize % (B+W) == 0`; `rotate90` translation only correct for square domains; staggered-cube y-clip only applies when `divisions < 2`; `_unit_cube` comment claims "open bottom" while including bottom faces.

**Good:** check/repair split with before/after reports; shared `_meshgraph` primitives; constrained-Triangle ground generation with explicit boundary/domain validation; copy-before-repair discipline; no `sys.exit`, no stray `print`.

---

## `udprep` package

### Critical

- **P1. `udprep_radiation.py:606-608, 664-690` — `run_short_wave_timedep` crashes or mis-shapes on vegetation.** `s_veg_all = np.zeros((albedo.size, nt))` is allocated with nfcts rows, but the solver returns `veg_absorb` of length nveg — and `np.zeros(0)`, never `None`, for no-veg/scanline runs (`directshortwave.py:1321, 1453, 1478`). First daytime step of any no-veg case raises `ValueError`; veg cases mis-shape whenever nveg ≠ nfcts. The unit test passes only because it mocks `s_veg=None` (`tests/test_udprep_core.py:690-691`).
- **P2. Vegetation radiation is physically wrong twice.** (a) `directshortwave.py:461-467, 657-665` vs `udprep_radiation.py:713-739`: `veg_absorb` is a transmittance-based quantity (the budget code multiplies it by irradiance) yet is written verbatim to `sveg.inp`/`timedepsveg.inp` under a "[W/m3]" header — off by the irradiance factor, time-varying in the timedep case. (b) `directshortwave.py:648-672`: the Möller kernel credits facet hits with *pre*-attenuation energy while vegetation in the same cell also absorbs from it — double counting in the tree-next-to-facade case (the facsec kernel gets the order right).
- **P3. `udprep_radiation.py:275, 297, 327, 350` — `calc_view_factors(maxD=...)` is validated then ignored**; the cache key and the View3D export both use `self.maxD`. Signature default (250) also contradicts the section default (1000). Bonus: the cache key here is a 6-tuple, `calc_short_wave` builds a 5-tuple — the two can never share the cache.

### Major

- **P4. `directshortwave.py:94-112, 1126-1129` — solver extent/blocking depend silently on `sim.Sc`.** Missing solid mask ⇒ Möller truncates the traced volume to ktot=2 (zero direct SW above the bottom two layers, no warning); facsec sets `has_solid=False` so rays are never blocked (sdir≈0). Reachable because `udbase._load_solid_masks` warns-and-sets-None. Raise (facsec); include mesh z-extent in kmax (Möller).
- **P5. `udprep_radiation.py:465-471, 499-503` — `calc_short_wave(nsun, irradiance, ...)` silently returns stale file contents when outputs exist**, ignoring its own arguments — a sun-position sweep returns the first run's answer every time. `run_short_wave`'s early return (522) has the same staleness for namoptions changes.
- **P6. `udprep.py:177-197` + `udprep_ibm.py:284-299` — `run_all()` on a fresh case feeds radiation stale (empty) `Sc`/`facsec`.** IBM writes the files but nothing reloads them into `sim`; the pipeline only works if IBM outputs existed at construction. Hidden ordering constraint: grid → ibm → *reload masks* → radiation.
- **P7. `udprep.py:186-197` vs `udprep_seb.py:30` — `run_all` gates SEB on `libm and lEB`**, but `SEBSection.run_all` also handles `iwallmom==2`/`iwalltemp==2`/`iwallmoist==2` without lEB — a non-EB case with `iwalltemp=2` never gets `Tfacinit.inp` from the pipeline. Remaining `**kwargs` are silently dropped when radiation doesn't run.
- **P8. `_section.py:259-262` — `write_changed_params` swallows all exceptions at INFO level.** Permissions errors or bugs silently leave namoptions un-updated while the user believes parameters persisted. Catch narrowly; warn or raise.
- **P9. `udprep_grid.py:167-187` — stretch-fit loop has no lower bound on `gf`.** With the shipped default `stretchconst=0.01`, `gf -= 0.01` reaches exactly 0.0 → 0/0 → NaN grid → literal "nan" text in `prof.inp`/`lscale.inp`, no error. Guard `gf > 0`; raise on failure to fit.
- **P10. `directshortwave.py:1406-1436` — facsec energy redistribution is two pure-Python loops over all sections, per time step.** For 10⁵–10⁷ sections this dominates all of preprocessing; vectorizes directly with `np.add.at`. Also `_build_cell_facet_lookup` (115-156) is a Python triple loop over triangles×cells, and `energy_in` (447, 608, 1244, 1336) is accumulated in the hot loop and never read.
- **P11. `directshortwave.py:1152-1153` vs `udprep_radiation.py:648, 673` — sun just above the horizon crashes the run.** `compute()` raises below |dir_z| < 1e-2 (zenith ≳ 89.43°) but the timedep loop gates only at zenith < 90 — one sunrise step aborts a multi-hour computation. Skip or floor such steps.
- **P12. `udprep_init.py:37-54` — `_parse_shell_config` catches only `CalledProcessError`;** on Windows without bash, `subprocess.run(["bash", ...])` raises `FileNotFoundError`, killing `UDPrep(...)` construction whenever `config.sh` exists. The text-parsing fallback is unreachable exactly where it's needed. Catch `(OSError, CalledProcessError)`.
- **P13. `udprep.py:120-130` — `addvar` rebuilds the section without `defaults`**, so change tracking (`show_changed_params`/`write_changed_params`) silently stops working for that section; previously computed outputs get promoted into `_fields`.

### f2py fallback behaviour

- There is no silent pure-Python fallback anywhere — missing extensions raise `RuntimeError` with build instructions (`udprep_ibm.py:301-308`, `directshortwave.py:1069-1076`). This fail-loud policy is correct; the messages should cross-reference the numba-based alternatives (`facsec`/`moller` for scanline; `backend="legacy"` for IBM).
- `directshortwave.py:1071` uses an absolute import (`import udprep.directshortwave_f2py`) while `udprep_ibm.py:303` uses relative — the absolute form breaks under vendoring; pick relative.
- The user cannot tell which solver path ran: method/backend logged only at INFO (invisible by default), the returned energy-budget dict is discarded, and budget closure (in ≈ fac+veg+sol+out) is never checked. Surface the chosen path and warn when the budget doesn't close.
- Equivalence effort is genuine (`_as_fortran_real_input` reproduces the legacy text→real pipeline bit-for-bit; pyf signatures match call sites). Committed binaries are cp312-win_amd64 only.
- Verified correct: the numba periodic-tiling scheme (wrap via `i % itot`, credit only `base_inside`, upstream pad) is self-consistent. Caveat: `max_ray_length = 10·max(xlen,ylen)` (1286, 1385) can truncate grazing-sun rays on small domains with tall geometry in periodic mode.

### Minor

- P14. `solar.py:1-9`, `directshortwave.py:1-9` — docstring after `from __future__` (module `__doc__` is None).
- P15. `udprep_grid.py:55-58, 71-72` — float-step `np.arange` endpoint pitfall; prefer `x0 + dx*np.arange(n)`.
- P16. `_section.py:45-54` — all sections write into the single shared `sim` namespace; cross-section coupling is implicit (e.g. IBM reads `sim.nfaclyrs` owned by seb); duplicate keys would be last-writer-wins.
- P17. `_section.py:131-139` — the `"num/den"` string-division hack misfires on any string default containing `/`.
- P18. `_section.py:93-94`, `udprep_ibm.py:66-73` — monkeypatching global `warnings.formatwarning`, duplicated in two places, process-global, not thread-safe.
- P19. `udprep.py:66, 90-91` — `DEFAULTS_JSON` set as *instance* attribute (defaults.json reparsed per UDPrep); unused imports; `_DefaultContext` defined mid-import-block.
- P20. `udprep_radiation.py:124` — `veg_key: id(veg_data)`; ids are reused after GC — a recreated dict at the same address silently reuses a stale solver.
- P21. `udprep_radiation.py:633-637, 789` — isolar=3 assumes exactly hourly rows starting 00:00; `_solar_state_time` truncates to whole hours — silently wrong for sub-hourly weather data.
- P22. `udprep_radiation.py:861-862` — `_read_weather_table` silently drops malformed rows.
- P23. `udprep_ibm.py:296-299` — `run_ibm` unconditionally deletes generically named files (`faces.txt`, …) from the case dir even on the f2py path that never created them.
- P24. `udprep_ibm.py:277-281` — `copy_geom_outputs` takes `sources[0]` of an unsorted glob; with multiple `facets.inp.*` variants an arbitrary one is copied.
- P25. `udprep_scalars.py:148-166` — `return` on an existing `scalarsourcep.inp.1.*` aborts writing sources for scalars 2..nsv (should be `continue`).
- P26. `udprep_forcing.py:136` — hardcoded surface anchor `[0,0,0,293.0,0]` injects 293 K regardless of `thl0`; `CubicSpline` extrapolates silently above the profile top.
- P27. `udprep_seb.py:53-55, 86, 90` — `fmt="%4f"` (width-4); "tempereatures" typo; `Tfac[:, :, -1]` bakes in an undocumented dimension order — add an explicit dimension check.
- P28. `udprep_vegetation.py:187, 274-275` — `load_stl` annotated `-> Dict[str, np.ndarray]` but returns a matplotlib figure and spawns a plot as a side effect (headless CI hazard).
- P29. `directshortwave.py:1547` — scanline `facet_points` come from `sim.geom.face_incenters` even when a `surface_mesh` override was supplied.
- P30. `solar.py:1113-1115` — SPA failures surface as `ValueError("SPA error code: 9")`; map codes to messages.

**Verified correct, for the record:** `interp_makima` matches MATLAB's makima exactly; `net_shortwave_reflections` matches the MATLAB iteration with a sensible denominator guard; the facsec kernel's solid-hit attribution matches the 0-based `locs` convention.

---

## `udvis`, examples, and documentation

### Major

- **V1. `udbase.py:1220-1222` — `sim.plot_trees()` delegates to `UDVis.plot_trees`, which does not exist** (only `plot_veg` does) — instant `AttributeError`, and `docs/udales-fields-tutorial.md` advertises the method. One-line fix: forward to `plot_veg`.
- **V2. `udvis/udbase_vis.py:695-769` — `plot_fac_type(building_ids=...)` accepts and silently ignores its argument** (forwarded by `udbase.py:1236-1246`, never referenced in the body).
- **V3. `udvis/udbase_vis.py:51-52` — a `UDVis` built from a `UDGeom` sets `self.sim = None`;** all sim-requiring methods (`plot_veg`, `plot_solid`, `plot_fac_type`, `plot_2dmap`, …) then crash with a bare `AttributeError` instead of a clear "requires a UDBase" error.
- **V4. Backend-parity gaps in `scene.py`:** `MeshPrimitive.opacity` honoured by plotly only; plotly colour-scale mapping supports only `{viridis, greys, greys_r}` and silently shows a Viridis colour*bar* for other cmaps while PyVista honours the real one; `_render_pyvista` adds a scalar bar unconditionally when `scene.colorbar` is set (latent crash plotly guards against); default cameras differ per backend (plotly az 225°/el 20° vs PyVista isometric+azimuth 180/el −10), undocumented.
- **V5. `udvis/udbase_vis.py:786-794` — `plot_building_ids` adds a second colorbar after `plt.show()` has already run**, and the `cmap="hsv"` colorbar misrepresents patches drawn with viridis.
- **E1. `examples/udprep_radiation_tutorial.py:162-163` — stray debugging `exit()` after `run_all`;** the three tutorial cases below it (~90 lines) never execute. Independently, `plot_shortwave` calls `fig.add_trace(...)` on the return of `plot_fac(show=False)`, which is now a PyVista `Plotter` by default → `AttributeError` even without the `exit()`.
- **E2. `examples/udprep_vegetation_tutorial.py:32-42` — written against a removed API:** `load_block` returns `None` and `load_stl` returns a dict (or a matplotlib figure, see P28), but the tutorial calls `.show(renderer="browser")` on both.
- **E3. Three notebooks hardcode the author's personal path** (`udbase_tutorial`, `facets_tutorial`, `fields_tutorial`, cell 2: `C:/Users/mvr/...`) and use `.resolve` without parens (passing a bound method as the path — the origin of the `callable(path)` hack, C16). *Note:* needing external experiment data (065/110/525/807) is by design for the data tutorials; the defects are the personal paths, the typo, and Run-All breakage — use a clearly-marked placeholder path instead. `fields_tutorial` cell 41's tree section says "uncomment if you have a tree simulation" but is live code, and cell 44 indexes experiment 525's tree data with experiment 110's grid.
- **D1. The docs site is a generation behind:** `mkdocs.yml:74-77` navigation exposes only the MATLAB tutorials; all four tutorial pages instruct `addpath(...tools/matlab)` and reference `docs/tutorial_mlx`, which does not exist (the files are in `docs/tutorial_udbase/`); `udales-post-processing.md` never mentions the Python tools; the Python API is invisible on the website.
- **D2. `setup_venv.sh:138` — the documented Quick Start fails on a fresh machine:** the import check requires `plotly` (now optional, moved to `requirements-plotly.txt`) and under `set -euo pipefail` aborts before building View3D/f2py; `pyvista` — the mandatory default backend — is *not* checked, so a broken install could pass.
- **D3. `setup_venv.sh:229-235` — the legacy root-`.venv` fallback documented in `README_VENV.md:90` is not implemented in the sh script** (dead conditional; only the ps1 implements it) — Linux users with a legacy venv get a second fresh one.
- **D4. Stale environment docs:** `udales-installation.md:25` says Python ≥ 3.6 (actual ≥ 3.9); `udales-preprocessing-windows.md:49` lists plotly as part of `requirements.txt` and omits `pyvista`/`trimesh`/`shapely`/`h5netcdf`.
- **D5. `udbase_tutorial.ipynb` uses `%matplotlib widget` but `ipympl` is not in `requirements.txt`** — the first cell errors in the venv the README tells users to register as the kernel.

### Minor

- V6. `udvis/udbase_vis.py:307-317` — `plot_veg` reports the 50 000-point subsample count as if it were the total; only disclosure is an invisible `logger.info` (`plot_solid`/`plot_fluid_boundary` do it right: "N of M shown").
- V7. `udvis/udbase_vis.py:83-90, 552-559` — per-face Python loops in `_collect_mesh_edges`, and O(edges×faces) scans in `_outline_segments(building_ids=...)` — minutes-long stalls on 10⁵–10⁶-face city meshes.
- V8. `scene.py:25-28` — `from exceptions import ...` tried *before* the relative fallback; a stray top-level `exceptions.py` on the path silently wins. `udgeom/udgeom.py:61` uses flat `from udvis import UDVis` with no relative fallback at all.
- V9. Docstring inaccuracies: `udbase_vis.py:1-6` says "matplotlib and pyvista" (omits plotly); `UDGeom.show()` claims trimesh is the plotly backend requirement.
- E4. Every notebook preamble promises an export example "at the end"; none has one (the sole `screenshot` mention is a comment).
- E5. `facets_tutorial.ipynb` cells 49 and 51 are byte-identical; cell 51's markdown intro suggests a different example was intended.
- E6. Unused `from unittest.mock import patch` in `geometry_tutorial`/`geometry_QA_tutorial` cell 2.
- E7. Untracked local leftovers in `examples/` (`aligned_cubes.stl`, `canyon.stl`, …) referenced by nothing.

**Hygiene (verified clean):** no `.pyc`/`__pycache__`/`.ipynb_checkpoints`/HTML artifacts are git-tracked; the STLs the notebooks need are tracked and whitelisted correctly.

---

## Tests and CI

### Critical

- **T1. The direct-shortwave physics tests never execute in any automated run.** `tests/_common.py:20-31` gates them on `UDALES_RUN_SLOW_TESTS` with a comment claiming CI sets it; the variable is set nowhere (`.github/workflows/ci.yml`, `tests/run_tests.py`, `tests/test_suites.yml` all checked — the `env_<NAME>` per-suite mechanism exists but no suite uses it). `test_directshortwave.py` and `test_directshortwave_periodic.py` silently skip everywhere, *including* the `experimental: python directshortwave unit` suite (`test_suites.yml:150-179`) that exists specifically to run them. These are the only tests validating the three ray tracers against the analytic flat-terrain flux. Fix: `env_UDALES_RUN_SLOW_TESTS: "1"` on those suites; correct the `_common.py` comment.

### Major

- **T2. The primary user-facing analysis API is untested at unit level.** `udbase.py` coverage stops at the constructor, scalar sources, trees, `read_matrix`, `_load_ncdata`, `load_facsec`, `convert_fac_to_field`. Untested: `load_field`, `load_stat_xyt/t/tree`, `load_slice`, `load_fac_momentum/eb/temperature`, `load_seb`, `assign_prop_to_fac`, `area_average_fac/seb`, `convert_facvar_to_field`, `convert_facflx_to_field`, `calculate_frontal_properties`. The repo-level MATLAB integration suite is a coarse net, not per-function coverage. The synthetic-NetCDF fixture pattern in `test_udbase_core.py:406-443` makes extending this cheap.
- **T3. Whole modules with zero tests:** `udprep/udprep_seb.py` (no characterization net under the active #313 redesign); `udvis/scene.py` composition (only error paths tested outside the GL-gated job); the compiled `ibm_preproc_f2py` extension is never invoked by unit tests (`test_udprep_core.py:177` mocks `run_ibm`) — covered only by the repo-level integration suite.

### Minor

- T4. `test_solar.py:17-60` hard-requires pvlib with no skip guard (fine in CI, errors on a trimmed install; add `skipUnless`).
- T5. 64 orphaned `case_*` dirs in `tests/.tmp` with fresh timestamps — something local resolves scratch there, and `_CaseDir.cleanup`'s `rmtree(ignore_errors=True)` (`_common.py:42-43`) swallows OneDrive/Windows lock failures; leaks accumulate silently (git-ignored, but OneDrive syncs them). Warn when rmtree leaves residue.
- T6. Unit tests depend on fixtures outside `tests/`: `examples/uDALES.stl` (`test_udgeom_api.py:44, 1446-1466`), repo-root `tests/cases/525` (`test_namelist.py:53`), repo-root `examples/101` (`test_directshortwave_periodic.py:15`). Reorganizing those breaks unit tests with no local signal.
- T7. The full Python unit suite runs identically in all four `build-and-supported-tests` matrix legs (ubuntu/macos × Debug/Release) — 4× cost, identical coverage.
- T8. `test_udbase_vis.py:585-602` hand-rolls monkeypatching of module globals with try/finally; `mock.patch` on both targets is safer against import refactors.
- T9. `tests/test_directshortwave.py:9-12` uses a divergent import preamble (`try/except ModuleNotFoundError` → `tools.python.tests._common`) — the first symptom of preamble drift (see S1).

**Good (genuinely):** assertions are scientific, not smoke — MATLAB golden files at atol=1e-12 across four fixtures including Xie–Castro (`test_udgeom_against_matlab.py:69-145`); analytic flux budgets for the ray tracers; pinned interpolation values and exact stretched-grid face positions; CubicSpline cross-checks. Negative paths and regressions are first-class (wrong-count rejection, "must not clobber authored `facets.inp`", Windows NetCDF handle-leak test, warning-hygiene assertions). Mocking is proportionate; numerics are never mocked. Fully deterministic; optional deps skip with reasons; platform quirks handled explicitly. CI runs the full suite plus a dedicated dual-backend viz job under xvfb.

---

## Structure and packaging

- **S1. No `pyproject.toml`/`setup.py`; five path-setup mechanisms coexist** (confirmed absent: pyproject, setup.py, pytest.ini, tox.ini, conftest.py):
  1. ~25 `try: import X / except ImportError: from . import X` blocks across source (see C9 for why this pattern is actively harmful);
  2. a 6-line `TESTS_DIR`/`sys.path.insert` preamble repeated in every test file plus `_common.py:33-34`;
  3. `sys.path.insert` cells in the notebooks;
  4. a different variant in `examples/udprep_tutorial.py:21-22`;
  5. `PYTHONPATH` injection in `setup_venv.sh:121, 176`.

  One divergent variant already slipped in (T9). Latent hazard: flat top-level modules named `exceptions`, `udconfig`, `udgrid`, `udstats` occupy the global namespace when the sys.path hack is active.

  **Recommendation (minimal, f2py-compatible):** one `pyproject.toml` (setuptools; `py-modules` for the flat files; `packages = udgeom, udprep, udvis`; dependencies from `requirements.txt`; `[project.optional-dependencies]` mapping `plotly` and `build` onto the other two requirements files), rename `exceptions.py` → `uderrors.py` (or fold everything under one `udales/` package), then `pip install -e tools/python` in `environment.yml`/CI/`setup_venv.*`. This deletes every dual-import block, every test preamble (a 2-line `_common.py` survives for `copy_case`/slow-gate), every notebook path cell, and the `callable(path)` hack. The editable install does not copy or build anything, so the in-place f2py `.pyd`/`.so` workflow is unaffected (extensions keep living inside `udprep/` and keep being built by `setup_venv.*`). The three-way requirements split is defensible and maps 1:1 onto extras.

- **S2. The committed f2py binaries are cp312-win_amd64 only** — any other interpreter/OS hits the RuntimeError path by design; fine, but worth stating in `fortran/README.md`.

---

## Systemic themes

1. **Error-handling philosophy is the number-one scientific risk.** Broad `except Exception` + `warnings.warn`/INFO logs repeatedly convert corrupt inputs and internal bugs into silently degraded results: uniform-grid substitution (C3), dropped facet data (C7), mis-sliced factypes (C2), stale physics (P5), unclosed energy budgets, swallowed writebacks (P8), empty outlines (G11). Adopt a library-wide rule: raise `DataFormatError` unless degradation is explicitly requested.
2. **Drift through duplication.** Two facsec loaders (C6), two building splitters (G4), three outline implementations (G27), two ground estimators, duplicated overlap-detection loops — each pair can disagree without either being wrong in isolation. Consolidation raises trustworthiness more than any single bug fix.
3. **Packaging debt is the root cause of much surface damage** — the dual-import blocks, the `callable(path)` hack, the notebook path cells, the preamble drift, the `exceptions.py` shadowing hazard. S1 retires all of it at once and does not disturb the Fortran workflow.
4. **Invisibility of diagnostics.** Critical signals (solver path chosen, budget closure, subsample counts, writeback failures) are logged at INFO, which no default configuration displays. Anything a user would act on belongs at WARNING or in a returned report.
