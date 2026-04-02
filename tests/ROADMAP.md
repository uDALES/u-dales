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
- `tree_sparse_compare` has been moved from `supported` to `experimental`
  because it is still better treated as developing coverage than as a trusted
  merge-gating path
- the regression harness is still only partially mature: it remains in the
  supported path, but its full run-and-compare case-output flow is still a
  later roadmap item

This means:

- Phase 1 is complete at an initial cleanup level.
- Phase 2 is complete at an initial merge-gate level.
- Phases 3 onward remain active roadmap items.

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
