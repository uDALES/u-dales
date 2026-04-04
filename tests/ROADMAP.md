# Test Roadmap

This document sketches a phased roadmap for evolving the uDALES test program
toward a stronger, more maintainable standard.

It uses the repository's current contract terminology:

- `supported`: stable merge-gating coverage
- `experimental`: important automated coverage for developing functionality
- `heavy`: resource-intensive or environment-specific coverage, often run on a cluster

It also keeps the distinction between test purpose and test contract:

- `unit`: isolated checks of one routine or API
- `integration`: end-to-end or cross-component checks of interacting workflows
- `regression`: checks that accepted behaviour has not drifted relative to a reference

## Aspiration

The long-term aspiration is for the uDALES test suite to reach a quality level
comparable to mature open CFD codes such as SU2, while also matching the needs
of a turbulence-resolving urban LES codebase more directly in the style of
PALM.

In practical terms, that means aiming for:

- a trusted merge-gating test selection
- a clear separation between unit, integration, and regression purposes
- a larger bank of committed reference cases
- routine regression comparisons against accepted behaviour
- stronger visibility into what is tested, what is experimental, and what is
  too heavy for routine CI

## External Benchmarks

Two public projects are especially relevant reference points:

- SU2: the public repository separates `UnitTests` and `TestCases`, and the
  project advertises dedicated regression testing in CI. It also maintains a
  separate public `TestCases` repository containing a broad library of
  reference problems used for regular regression runs.
- PALM Model System: the official public repository is on GitLab rather than
  GitHub. The public tree exposes dedicated test areas under
  `packages/palm/model/tests/`, including `builds/` and `cases/`, and also has
  additional package-specific test areas such as
  `packages/dynamic_driver/inifor/tests/`. Public indexed paths also show that
  PALM keeps committed case inputs and monitoring/reference outputs inside the
  test-case tree, and ships dedicated shared testing utilities under
  `packages/palm/model/share/palmtest/`. This is a particularly relevant
  comparison for uDALES because PALM is also a turbulence-resolving LES model
  with urban and environmental modelling use cases.

The target is not to copy either project mechanically. The target is to reach a
similar level of discipline:

- SU2-like clarity in unit-vs-case-vs-regression structure
- PALM-like relevance in solver and workflow test coverage for an LES code

Relevant public references:

- SU2 repository: https://github.com/su2code/SU2
- SU2 TestCases repository: https://github.com/su2code/TestCases
- PALM repository: https://gitlab.palm-model.org/releases/palm_model_system
- PALM model tests tree: https://gitlab.palm-model.org/releases/palm_model_system/-/tree/v25.04-rc.1/packages/palm/model/tests
- PALM model case tests: https://gitlab.palm-model.org/releases/palm_model_system/-/tree/v25.04-rc.1/packages/palm/model/tests/cases
- PALM shared test utilities: https://gitlab.palm-model.org/releases/palm_model_system/-/tree/v25.04-rc.1/packages/palm/model/share/palmtest
- PALM package-specific tests: https://gitlab.palm-model.org/releases/palm_model_system/-/tree/v21.10/packages/dynamic_driver/inifor/tests

## Mapping To uDALES

uDALES already follows a similar overall direction in several important ways.
The remaining work is mainly to deepen and harden the structure rather than to
invent it from scratch.

Current rough mapping:

- SU2 `UnitTests` and PALM package-specific test areas map to
  `tools/python/tests/` and, over time, to any future isolated solver-side unit
  tests under `tests/unit/`.
- SU2 `TestCases` and PALM `packages/palm/model/tests/cases/` map to
  `tests/cases/`, `tests/integration/`, and the future expansion of
  `tests/regression/` case assets.
- PALM shared test tooling under `packages/palm/model/share/palmtest/` maps to
  the role currently played by `tests/run_tests.py`, `tests/test_suites.yml`,
  and helper scripts under `tests/regression/david_tests/scripts/`.
- PALM's separation between core model tests and package-specific tests maps to
  the existing uDALES split between repo-level solver tests in `tests/` and
  Python-tooling tests in `tools/python/tests/`.
- SU2's explicit regression identity maps to the intended long-term role of
  `tests/regression/`, even though the current uDALES regression harness is not
  yet fully realized as a case-output comparison program.

This means uDALES is already aligned with the benchmark projects in these ways:

- it already distinguishes test purpose from implementation language
- it already keeps committed shared case fixtures
- it already separates Python/tooling checks from solver-facing checks
- it already has a curated dispatcher for supported and experimental contracts
- it already documents the distinction between integration, regression, and
  heavier test classes

The main remaining differences are:

- benchmark projects have a larger and more mature bank of routinely exercised
  reference cases
- benchmark projects expose a more complete regression workflow around those
  cases
- benchmark projects make the merge-gating versus non-gating contract more
  operationally visible through their testing infrastructure

## Current Status

The repository has now completed an initial pass of phases 1 and 2.

Current state:

- `tests/test_suites.yml` is the source of truth for the current `supported`
  and `experimental` memberships.
- stale test metadata has been corrected, including outdated fixture paths in
  the lightweight Python unit layer
- `supported` now includes lightweight Python unit coverage in addition to the
  existing solver-facing checks
- the regression harness is still only partially mature: it remains in the
  supported path, but its full run-and-compare case-output flow is still a
  later roadmap item

This means:

- Phase 1 is complete at an initial cleanup level.
- Phase 2 is complete at an initial merge-gate level.
- Phases 3 onward remain active roadmap items.

## v3.0 Alignment

As of 2026-04-03, the GitHub repository has an open `3.0` milestone
([milestone 10](https://github.com/uDALES/u-dales/milestone/10)) with the
following milestone description:

- `V3.0 will make some substantial changes to input and output data, including`
  `intermediate files (after pre-processing the data) and post-processing the`
  `data (particularly not every field file will output their own data). There`
  `will be changes to filenames, variable names etc.`

At that date the milestone contains 18 tracked items, with 16 open and 2
closed. The issue set is consistent with the local repo direction and makes the
major-upgrade priorities much clearer than the local tree alone.

The main `v3.0` themes visible in the GitHub plan are:

- IO and statistics refactors:
  [#242](https://github.com/uDALES/u-dales/issues/242),
  [#276](https://github.com/uDALES/u-dales/issues/276),
  [#277](https://github.com/uDALES/u-dales/issues/277),
  [#278](https://github.com/uDALES/u-dales/issues/278),
  [#304](https://github.com/uDALES/u-dales/pull/304),
  [#305](https://github.com/uDALES/u-dales/pull/305)
- input-format evolution and schema/tooling work:
  [#279](https://github.com/uDALES/u-dales/issues/279),
  [#305](https://github.com/uDALES/u-dales/pull/305)
- preprocessing and postprocessing updates for new file contracts:
  [#196](https://github.com/uDALES/u-dales/issues/196),
  [#280](https://github.com/uDALES/u-dales/issues/280),
  [#281](https://github.com/uDALES/u-dales/issues/281)
- documentation catch-up for the new formats and workflows:
  [#258](https://github.com/uDALES/u-dales/issues/258),
  [#282](https://github.com/uDALES/u-dales/issues/282)
- solver cleanup that supports the release quality bar:
  [#68](https://github.com/uDALES/u-dales/issues/68),
  [#308](https://github.com/uDALES/u-dales/issues/308)
- geometry/facet robustness still considered part of the major release:
  [#265](https://github.com/uDALES/u-dales/issues/265)

For roadmap purposes, this means the test program should not be treated as an
isolated concern. It should be used to de-risk the actual GitHub-scoped `v3.0`
upgrade, especially where the milestone explicitly changes the user-facing data
contract.

## Action Plan From Current Review

The recent repo review identified a small set of architectural and workflow
gaps that matter disproportionately for a major release. The most useful way to
act on them is as a set of parallel workstreams rather than as isolated fixes.

### Workstream A: Build And Environment Reproducibility

Goal: make solver and preprocessing builds reproducible enough for a major
release and easier to reason about in CI, on local systems, and on HPC.

Priority actions:

- replace helper-project side effects in the top-level CMake configure path
  with a more standard dependency-discovery approach where practical
- reduce GNU-specific assumptions in preprocessing builds, especially in the
  f2py and OpenMP wrapper paths
- make Python-header and compiled-wrapper checks consistent across all
  preprocessing targets, not only the direct shortwave path
- tighten the contract between `tools/build_executable.sh`,
  `tools/build_preprocessing.sh`, CI, and the installation docs so they describe
  one coherent supported workflow
- treat [#196](https://github.com/uDALES/u-dales/issues/196) as a key `v3.0`
  build-system deliverable rather than just a preprocessing convenience

`v3.0` relevance:

- this is foundational work for any major release that wants broader adoption
- it directly supports the newer Python-first tooling direction
- it narrows one of the biggest current gaps versus SU2 and PALM

Suggested contract:

- start as documentation + build-system cleanup
- add or promote lightweight build-smoke tests to `supported` only after the
  resulting path is stable on standard CI runners

### Workstream B: Test Infrastructure Hardening

Goal: make the curated test layer trustworthy enough to carry a larger release.

Priority actions:

- replace the handwritten parser in `tests/run_tests.py` with a proper YAML
  loader so `tests/test_suites.yml` can evolve safely
- keep the manifest as the source of truth, but validate it structurally rather
  than relying on indentation-sensitive custom parsing
- improve per-suite reporting so failures leave behind enough context to debug
  without rerunning everything locally
- align manifest and runner evolution with the schema-driven input work in
  [#305](https://github.com/uDALES/u-dales/pull/305), so new input-contract
  tests do not depend on brittle runner behaviour

`v3.0` relevance:

- the `v3.0` upgrade will likely add more suites, more metadata, and more
  branch-specific exceptions if the runner remains fragile
- this is low-risk engineering work with a high leverage payoff

Suggested contract:

- keep the migration itself in `supported`
- treat any new richer reporting or artifact upload logic as `experimental`

### Workstream C: MPI Averaging Regression Lock-In

Goal: prove that the recent intrinsic-averaging and reduction refactor has not
changed accepted behaviour relative to `v2.2`.

Priority actions:

- add a dedicated regression test that compares `v2.2` and current-branch
  outputs for all actively used averaging/reduction paths:
  `av_intr`, `av_y_intr`, `sum_y_intr`, `sum_x_intr`, `reduce_xy_sum`, and
  `reduce_yz_sum`
- include both IBM and non-IBM coverage so intrinsic masking and all-fluid
  behaviour are both exercised
- explicitly compare staggered locations (`LOC_C`, `LOC_U`, `LOC_V`, `LOC_W`,
  `LOC_UV`, `LOC_WU`, `LOC_VW`) where they are scientifically meaningful
- include at least one case that would expose a ghost-cell/windowing error and
  one case that would expose an intrinsic-versus-comprehensive normalization
  error
- treat this test as a release-safety check for the averaging refactor before
  expanding the abstraction further

`v3.0` relevance:

- this protects the scientific meaning of slab/profile statistics during the
  broader IO/statistics changes already planned for `v3.0`
- it reduces the risk that refactored reduction APIs silently harden a changed
  definition into the release contract

Suggested contract:

- start as `experimental` while the `v2.2` comparison harness is being built
- promote to `supported` once the reference workflow is stable and cheap enough
  to run routinely

### Workstream D: MPI Abstraction Boundary

Goal: remove direct MPI calls from solver and workflow code outside the
parallel-architecture layer, following the `sparkle` pattern more closely.

Priority actions:

- treat [src/architecture.f90](/rds/general/user/mvr/home/udales/u-dales/src/architecture.f90)
  as the single home for collective abstractions and shared communicator state
- replace direct `MPI_BCAST`, `MPI_ALLREDUCE`, `MPI_ABORT`, and related calls
  in the rest of `src/` with named wrappers or generic interfaces exported from
  `architecture`
- migrate in phases, starting with repeated scalar/vector broadcasts in:
  `modfielddump`, `modstatsdump`, `modsubgrid`, `modtimedep`, `initfac`,
  `modEB`, `modibm`, and `modinlet`
- add focused wrappers for the operations that are already repeated many times,
  for example:
  broadcast of scalars and arrays, global sums, global maxima, and domain
  barriers
- keep low-level neighbour exchange kernels in the parallel layer, but stop
  calling raw MPI collectives directly from physics, IO, and setup modules

Current raw-MPI hotspots outside the parallel layer include:

- `tstep.f90`, `modchecksim.f90`, and `modforces.f90` for repeated
  `MPI_ALLREDUCE`
- `modfielddump.f90`, `modstatsdump.f90`, `modsubgrid.f90`,
  `modtimedep.f90`, `modinlet.f90`, `initfac.f90`, `modEB.f90`,
  `modibm.f90`, `modpurifiers.f90`, and `modglobal.f90` for repeated
  `MPI_BCAST`
- `modsave.f90` and some safety paths for direct `MPI_ABORT`

`v3.0` relevance:

- a major release is the right time to tighten architectural boundaries
- the same abstraction layer will make future test harnesses and portability
  work easier, because collective semantics will no longer be scattered through
  the solver

Suggested contract:

- do the wrapper introduction first without changing semantics
- migrate module-by-module with regression coverage after each cluster of
  replacements
- avoid mixing this refactor with unrelated physics changes

### Workstream C: Scientific Semantics And Reduction Contracts

Goal: make slab-averaging and masked-reduction semantics explicit before `v3.0`
locks in new IO, statistics, and postprocessing contracts.

Priority actions:

- treat `src/modmpi.f90` `avexy_ibm` as an intrinsic-average routine by
  definition, because it computes `sum(var * II) / IIs` rather than a
  comprehensive/superficial average over total slab area
- adopt a reduction naming rule in `src/modmpi.f90` and `src/definitions.f90`:
  reduction results are global by default, and only explicitly processor-local
  routines should use a `_local` suffix
- rename `avexy_ibm` to a semantics-explicit name such as `av_intr` or a
  longer intrinsic-average variant once the migration path is agreed
- review all current `avexy_ibm` call sites and classify them as:
  clearly intrinsic and valid, ambiguous, or scientifically invalid
- prioritize `src/modforces.f90` for that review, because `uvol`, `vvol`, and
  `fluidvolume` appear to use an intrinsic average where a masked sum or
  comprehensive quantity may be required instead
- add a focused regression test that distinguishes intrinsic averages from
  comprehensive/superficial reductions on an IBM-masked slab, so future helper
  refactors cannot silently change the scientific meaning
- document the distinction alongside the `xy`/`xyt` output contract and tie it
  to the existing user documentation on intrinsic averages

`v3.0` relevance:

- this is a scientific-correctness risk, not just a naming or cleanup issue
- `v3.0` is already changing filenames, variable names, and output contracts,
  so this is the right release to make the averaging semantics explicit
- if invalid call sites survive the IO/statistics refactor, downstream
  postprocessing could preserve a wrong definition in a much harder-to-detect
  form

Known hotspot to flag explicitly:

- `src/modforces.f90` currently uses `avexy_ibm` in places where the intended
  quantity may be a total fluid-area or fluid-volume integral rather than an
  intrinsic mean; this should be treated as a release-blocking review item for
  the `v3.0` statistics and solver cleanup track
  until it proves stable in CI

### Workstream C: Regression Coverage For Major-Upgrade Features

Goal: turn accepted `v3.0` behaviour into protected behaviour.

Priority actions:

- formalize canonical regression cases for facet, tree, and radiation-related
  workflows that are changing most in the major upgrade
- add explicit regression protection for the `v3.0` IO/statistics change set in
  [#242](https://github.com/uDALES/u-dales/issues/242),
  [#276](https://github.com/uDALES/u-dales/issues/276),
  [#277](https://github.com/uDALES/u-dales/issues/277), and
  [#278](https://github.com/uDALES/u-dales/issues/278)
- keep the existing `526` vegetation comparison, but broaden the reference set
  to cover the new facet/tree handling and Python preprocessing paths
- prioritize output comparisons that reflect actual user-visible contracts:
  generated inputs, field dumps, tree dumps, radiation outputs, and selected
  statistics or budgets
- make solver decomposition-parity checks part of this story rather than a
  separate diagnostic island

`v3.0` relevance:

- the Oct 2025 facet/tree upgrade note in `tools/matlab/udbase.m` is exactly
  the kind of change that should gain explicit regression protection before a
  major release
- this is the main mechanism for preventing silent behavioural drift while the
  upgrade is still moving

Suggested contract:

- start new major-upgrade coverage in `experimental`
- promote only the cheapest and most trusted reference cases to `supported`
- keep larger or cluster-only comparisons in `heavy`

### Workstream D: Solver Robustness And Failure Reporting

Goal: make solver failures easier to interpret during `v3.0` stabilization.

Priority actions:

- reduce reliance on scattered raw `stop 1` aborts in key startup, boundary,
  IBM, and Poisson-solver paths
- move toward more consistent fatal-error reporting that preserves context such
  as the subsystem, namelist setting, file path, or rank-local condition that
  triggered the stop
- add targeted tests for known invalid-input and unsupported-configuration
  branches where the solver should fail clearly rather than opaquely
- use [#68](https://github.com/uDALES/u-dales/issues/68) and
  [#308](https://github.com/uDALES/u-dales/issues/308) as the anchor issues
  for solver-side cleanup that should be test-backed before release

`v3.0` relevance:

- this matters most when interfaces and input contracts are changing, which is
  exactly what a major upgrade tends to do
- clearer failure modes reduce the support cost of bringing users onto the new
  release

Suggested contract:

- add narrow invalid-input tests as `experimental`
- only promote them to `supported` once they are deterministic across supported
  MPI environments

### Workstream E: Developer Workflow And Documentation

Goal: make the `v3.0` upgrade legible to contributors and users.

Priority actions:

- replace placeholder pages such as `docs/udales-workflow.md` and incomplete
  docs metadata such as `docs/udales-docs-software.md`
- document one recommended developer path for local Linux/WSL, macOS, and
  cluster workflows instead of leaving contributors to infer the supported
  combinations from scripts
- make the relationship between MATLAB legacy tooling and Python-first tooling
  explicit, especially for upgraded facet/tree workflows
- tie the documentation plan directly to the active `v3.0` docs issues
  [#258](https://github.com/uDALES/u-dales/issues/258) and
  [#282](https://github.com/uDALES/u-dales/issues/282), rather than treating
  docs as a generic follow-up

`v3.0` relevance:

- a major release needs a coherent contributor and user story, not just working
  code
- this is where PALM currently presents a more operationally complete model
  workflow than uDALES

Suggested contract:

- no special test contract is needed for prose itself
- but docs should point to the real `supported`, `experimental`, and `heavy`
  workflows and should stop describing obsolete paths

## Recommended Order For v3.0

To keep the major upgrade tractable, the recommended order is:

1. Harden the build and test infrastructure first.
2. Add regression protection around the facet/tree/Python-preprocessing upgrade.
3. Improve solver failure reporting in the highest-friction paths.
4. Clean up developer and installation documentation around the stabilized path.

This ordering is intentionally conservative. A major release becomes much less
risky if the project first improves its ability to build, test, compare, and
diagnose the code it already has.

## Core Principles

- New functionality should usually get `experimental` unit or integration tests first.
- Once behaviour stabilizes, add regression protection where appropriate.
- Promote mature, reliable coverage from `experimental` to `supported`.
- Keep `supported` small, trusted, and suitable for standard CI.
- Put expensive or environment-specific coverage in `heavy`.

## Phase 1: Make The Current Suite Honest

Goal: align documentation, manifests, and actual behaviour.

- Audit existing automated tests and assign both a purpose and a contract.
- Fix stale tests, missing fixtures, and outdated paths.
- Update documentation so it reflects what really runs today.
- Keep `tests/test_suites.yml` as the source of truth for contract membership.

This phase is mainly cleanup, but it is necessary before scaling the suite.

## Phase 2: Make `supported` A Real Merge Gate

Goal: establish a small, stable, trusted CI gate.

- Limit `supported` to tests that are reliable in standard GitHub Actions.
- Include only coverage that the project is willing to block merges on.
- Keep the runtime and environment requirements modest.
- Ensure solver-facing integration tests in `supported` are representative and maintainable.

This is the main PR gate and should be treated as the most trusted layer.

## Phase 3: Use `experimental` As The Incubation Lane

Goal: give new functionality automated coverage early without overcommitting the merge gate.

- Add new unit or integration tests to `experimental` first when interfaces are still changing.
- Use `experimental` for new physics paths, preprocessing/tooling, and tests with evolving tolerances.
- Run `experimental` regularly when practical, but do not require it to block merges until it is mature.
- Define clear expectations for promoting suites from `experimental` to `supported`.

`experimental` is not throwaway coverage. It is the place where new test coverage starts.

## Phase 4: Restore True Regression Testing

Goal: protect accepted behaviour, not just buildability.

- Finish the actual run-and-compare path in `tests/regression/david_tests/run_tests.py`.
- Start with a small set of canonical cases.
- Compare outputs that matter, such as field dumps, statistics, budgets, and restart outputs.
- Define explicit tolerances and document why they are acceptable.
- Put mature, cheap regression cases in `supported`.
- Keep unstable or still-tuning regression checks in `experimental`.

This phase is the largest current gap between the existing test program and a more mature CFD testing posture.

## Phase 5: Add `heavy` For Expensive Or Cluster-Oriented Coverage

Goal: give large or environment-specific tests a clear home.

- Use `heavy` for larger MPI-count runs, longer LES cases, and broader physics combinations.
- Keep `heavy` coverage out of routine GitHub Actions unless it proves cheap and stable enough.
- Document how heavy tests are expected to run: manually, on the cluster, or in a dedicated external workflow.
- Treat `heavy` as part of the overall test program, even if it is not part of the standard merge gate.

This avoids overloading `supported` with tests that are valuable but operationally unsuitable for routine CI.

## Phase 6: Improve CI Structure And Visibility

Goal: make the test program legible and actionable for contributors.

- Make CI reporting reflect the test contracts clearly.
- Treat `supported` as required coverage.
- Treat `experimental` as optional or non-blocking until suites are promoted.
- Keep `heavy` outside routine GitHub Actions unless a dedicated path exists.
- Upload logs, comparison artifacts, and case-level summaries on failure where useful.

This phase improves developer confidence and makes failures easier to interpret.

## Phase 7: Promote Mature Coverage Over Time

Goal: steadily improve coverage without bloating the required gate.

- Promote mature suites from `experimental` to `supported`.
- Add regression tests after new workflows have stabilized.
- Keep representative cheap cases in `supported`.
- Keep larger or more expensive variants in `heavy`.
- Revisit suite membership periodically so `supported` stays lean and trusted.

The long-term direction is:

- integration first for new workflows
- regression after stabilization for accepted behaviour
- `supported` for trusted gates
- `experimental` for growth
- `heavy` for scale

## Practical Interpretation

In this repository, the main decision rule should be:

- If the main risk is workflow breakage, add an integration test.
- If the main risk is behaviour drift in an accepted path, add a regression test.
- If both are true, start with integration and add regression once the feature stabilizes.

That keeps the distinction between test purpose and test contract clear while matching how uDALES development currently works.
