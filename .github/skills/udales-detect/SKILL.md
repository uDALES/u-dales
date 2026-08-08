---
name: udales-detect
description: Use when users ask how to install, build, run, or test uDALES on this machine (Windows/WSL/Linux/HPC), what modules/toolchain to use, or how to detect the current environment before execution; gather environment fingerprints and route execution tasks to udales-exec.
---

# uDALES Detection Skill

This skill identifies the execution environment, gathers system fingerprints,
checks for site documentation, and determines which uDALES workflow is needed.

## When to use

Use this skill when a user asks:

- "How do I run this here?"
- "What modules do I need?"
- "Which MPI/Python stack should I use on this system?"
- "How do I install preprocessing on my machine?"
- "Should I use Windows, WSL, local Linux, or cluster modules for this?"
- "Can you detect this environment first and then give run/build steps?"

## Steps

1. Run the local probe to capture a baseline fingerprint.
2. Check known cluster profiles in `references/clusters.md`.
3. If no match, query online documentation for the host/cluster.
4. Ask what the user needs:
   - preprocessing and postprocessing
   - solver execution
   - both
5. Based on the answer, direct execution work to the `udales-exec` skill.
6. Record new details in `references/clusters.md`.

## Local probe

Run:

```bash
.github/skills/udales-detect/scripts/detect_env.sh
```

This prints hostname, OS, module availability, mpiexec flavor, and Python paths.

## Online lookup

If the cluster is unknown, search for:

- cluster name or hostname suffix (e.g. `*.imperial.ac.uk`)
- module environment documentation
- MPI launcher guidance
- Python module guidance

Record any relevant links in `references/sources.md`.

## Minimal checks

Only run lightweight checks during detection:

- `mpiexec --version`
- `python3 --version`

Leave build/run checks to `udales-exec` after the user confirms the intent.

## Files to consult

- `docs/cluster_workflows.md`
- `tools/build_executable.sh`
- `tools/build_preprocessing.sh`
- `tests/README.md`
- `.github/skills/udales-exec/references/clusters.md`
