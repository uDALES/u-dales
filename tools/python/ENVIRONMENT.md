# Python Environment Setup

This repository ships a small Python toolchain for post-processing, tutorials, and tests. Use the steps below to create a reproducible virtual environment in the repo root.

## 1) Create and activate a virtual environment
From the repository root (`u-dales/`):

- Windows (PowerShell)
  - `python -m venv .venv`
  - `./.venv/Scripts/Activate.ps1`
- macOS/Linux (bash/zsh)
  - `python3 -m venv .venv`
  - `source .venv/bin/activate`

## 2) Install required packages
Install the packages needed for the Python tools, notebooks, and test harness:

```
pip install --upgrade pip
pip install -r tools/python/requirements.txt
```

## 3) (Optional) Register a Jupyter kernel
If you plan to run the tutorial notebooks, register the environment as a kernel:

```
python -m ipykernel install --user --name udales-py --display-name "uDALES Python"
```

## Package list
Core: numpy, scipy, xarray, netCDF4, trimesh, matplotlib, f90nml

Notebook: ipykernel (installed by the requirements file)

Optional: pyvista (advanced 3D visualization)

## Notes
- Keep `.venv/` untracked (it is already ignored).
- Re-run `pip install -r tools/python/requirements.txt` after pulling changes to stay current.
