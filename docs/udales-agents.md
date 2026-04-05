# Using Agents

This page describes how we use Codex agents for uDALES development work.
It is aimed at contributors working on this repository.

## When to use an agent

Agents are useful for:

- repeated operational workflows on clusters
- test execution and triage across multiple suites
- refactors with careful bookkeeping
- large documentation updates with cross-file consistency

Use an agent when the task benefits from:

- consistent, repeatable steps
- multi-file updates with traceability
- detailed run logs or comparisons

## Source of truth

Project-specific instructions for agents live in:

- `AGENTS.md` (repo-level guidance)
- `agent_skills/udales-exec/` (cluster execution runbook skill)

When in doubt, prefer the repo wrappers in `tools/` and the curated test
selectors under `tests/`.

## Skill usage

The `udales-exec` skill is intended for cluster execution and test workflows.
It detects the environment, loads the correct modules, and chooses the right
build/run/test commands.

## First run on a new system

When you first use a new machine or cluster, run the lightweight self-test to
confirm that wrappers, Python tooling, and the environment fingerprint are
correct before attempting heavier builds or MPI jobs. If it fails, update the
cluster profile in `agent_skills/udales-exec/references/clusters.md` and rerun.

To run the self-test:

```bash
UD_VENV=/rds/general/user/$USER/home/udales/.venv
UD_PYTHON=python3
UD_VENV="$UD_VENV" UD_PYTHON="$UD_PYTHON" agent_skills/udales-exec/scripts/skill_selftest.sh
```

Use `--full` for execution-only verification with a prebuilt solver binary:

```bash
UD_VENV=/rds/general/user/$USER/home/udales/.venv
UD_PYTHON=python3
UD_BUILD=/path/to/u-dales
UD_VENV="$UD_VENV" UD_PYTHON="$UD_PYTHON" UD_BUILD="$UD_BUILD" \
  agent_skills/udales-exec/scripts/skill_selftest.sh --full
```

## Do and don’t

- Do keep the top-level repo workflow solver-first.
- Do use `tests/run_tests.py` for curated selections.
- Do record new cluster profiles in `agent_skills/udales-exec/references/clusters.md`.
- Don’t invent ad hoc module stacks when repo wrappers exist.
- Don’t run heavy MPI tests inside sandboxed environments.
