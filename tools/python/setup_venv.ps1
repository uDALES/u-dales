<#
.SYNOPSIS
    Setup script for Windows.
    Creates a Python virtual environment, installs dependencies, and builds
    the required preprocessing artifacts (View3D and f2py extension modules)
    via CMake and Ninja.

.DESCRIPTION
    This script is the Windows equivalent of tools/python/setup_venv.sh.
    It performs the same steps on native Windows PowerShell:
      - validates build_system and build_target arguments
      - selects or creates the virtual environment directory
      - checks that the chosen Python interpreter is available
      - creates the virtual environment
      - installs runtime dependencies from requirements.txt
      - installs build-time dependencies from requirements-build.txt
      - builds the preprocessing tools (View3D and/or f2py extension modules)
        via CMake/Ninja (since build_preprocessing.sh is bash-only)

    Run from the repository root.
    Requires: Python 3.9+, CMake, Ninja, GCC, GFortran on PATH.

    If the repository path contains spaces (e.g. OneDrive), map the parent
    directory to a drive letter first to avoid f2py failures:
        subst U: "C:\Users\<user>\...\uDALES"
        cd /d U:\u-dales
    Then use that Python as the interpreter:
        $env:PYTHON_BIN = "U:\.venv\Scripts\python.exe"

.PARAMETER BuildSystem
    Build environment to use.
    Allowed values : common
    Default        : common
      common - local Windows system (no module loading)
    Note: 'icl' (Imperial College London HPC) is Linux-only and not
    supported on Windows.

.PARAMETER BuildTarget
    CMake preprocessing target to build.
    Allowed values : view3d | preprocessing_tools
    Default        : preprocessing_tools
      view3d               - View3D executable only (requires a C++ compiler)
      preprocessing_tools  - View3D + directshortwave and IBM f2py modules
                             (requires GCC, GFortran, numpy, and Python headers)

.PARAMETER PythonBin
    Python interpreter to use for creating the venv.
    Default: value of $env:PYTHON_BIN, or 'python' if not set.

.PARAMETER VenvDir
    Explicit path for the virtual environment directory.
    Default: value of $env:VENV_DIR, or tools\python\.venv inside the repo.
    Falls back to .venv at the repo root if that legacy directory exists.

.EXAMPLE
    # Default local setup (from repository root)
    .\tools\python\setup_venv.ps1

.EXAMPLE
    # View3D only
    .\tools\python\setup_venv.ps1 -BuildTarget view3d

.EXAMPLE
    # Use a specific Python interpreter
    $env:PYTHON_BIN = "C:\Python312\python.exe"
    .\tools\python\setup_venv.ps1
#>

param(
    [string]$BuildSystem = "common",
    [string]$BuildTarget = "preprocessing_tools",
    [string]$PythonBin   = $(if ($env:PYTHON_BIN) { $env:PYTHON_BIN } else { "python" }),
    [string]$VenvDir     = ""
)

$ErrorActionPreference = "Stop"

$SCRIPT_DIR  = Split-Path -Parent $MyInvocation.MyCommand.Path
$UDALES_ROOT = (Get-Item (Join-Path $SCRIPT_DIR "..\..")).FullName

# --- Validate build_system ---
if ($BuildSystem -ne "common") {
    if ($BuildSystem -eq "icl") {
        Write-Host "Error: build system 'icl' is only available on the ICL Linux HPC cluster and cannot be used on Windows." -ForegroundColor Red
    } else {
        Write-Host "Error: invalid build system '$BuildSystem'." -ForegroundColor Red
    }
    Write-Host ""
    Write-Host "Usage: .\tools\python\setup_venv.ps1 [-BuildSystem <system>] [-BuildTarget <target>]"
    Write-Host ""
    Write-Host "  -BuildSystem   Build environment to use (default: common)"
    Write-Host "                   common  - local Windows system"
    Write-Host ""
    Write-Host "  -BuildTarget   CMake target to build (default: preprocessing_tools)"
    Write-Host "                   view3d               - View3D executable only"
    Write-Host "                   preprocessing_tools  - View3D + f2py extension modules"
    exit 1
}

# --- Validate build_target ---
if ($BuildTarget -ne "view3d" -and $BuildTarget -ne "preprocessing_tools") {
    Write-Host "Error: invalid build target '$BuildTarget'." -ForegroundColor Red
    Write-Host ""
    Write-Host "Usage: .\tools\python\setup_venv.ps1 [-BuildSystem <system>] [-BuildTarget <target>]"
    Write-Host ""
    Write-Host "  -BuildSystem   Build environment to use (default: common)"
    Write-Host "                   common  - local Windows system"
    Write-Host ""
    Write-Host "  -BuildTarget   CMake target to build (default: preprocessing_tools)"
    Write-Host "                   view3d               - View3D executable only"
    Write-Host "                   preprocessing_tools  - View3D + f2py extension modules"
    exit 1
}

# --- Resolve venv directory (same priority chain as setup_venv.sh) ---
$DEFAULT_VENV_DIR = Join-Path $UDALES_ROOT "tools\python\.venv"
$LEGACY_VENV_DIR  = Join-Path $UDALES_ROOT ".venv"

if ($VenvDir -ne "") {
    # explicit override via parameter or env var
} elseif ($env:VENV_DIR) {
    $VenvDir = $env:VENV_DIR
} elseif (Test-Path $DEFAULT_VENV_DIR) {
    $VenvDir = $DEFAULT_VENV_DIR
} elseif (Test-Path $LEGACY_VENV_DIR) {
    $VenvDir = $LEGACY_VENV_DIR
} else {
    $VenvDir = $DEFAULT_VENV_DIR
}

# --- Header ---
Write-Host "==========================================" -ForegroundColor Cyan
Write-Host "Setting up Python virtual environment" -ForegroundColor Cyan
Write-Host "==========================================" -ForegroundColor Cyan
Write-Host "u-dales repository: $UDALES_ROOT"
Write-Host "Python interpreter: $PythonBin"
Write-Host "Virtual environment: $VenvDir"
Write-Host "Build system:        $BuildSystem"
Write-Host "Build target:        $BuildTarget"
if ($VenvDir -eq $LEGACY_VENV_DIR) {
    Write-Host "Note: using legacy repo-root virtual environment. New setups default to $DEFAULT_VENV_DIR" -ForegroundColor Yellow
}
Write-Host ""

# --- Check Python ---
try {
    $pythonVersion = & $PythonBin --version 2>&1
    Write-Host "Python version: $pythonVersion"
} catch {
    Write-Host "Error: '$PythonBin' is not available. Install Python 3.9+ and ensure it is on PATH." -ForegroundColor Red
    exit 1
}

# --- Create or recreate venv ---
if (Test-Path $VenvDir) {
    Write-Host "Virtual environment already exists at: $VenvDir" -ForegroundColor Yellow
    $response = Read-Host "Do you want to recreate it? (y/N)"
    if ($response -eq 'y' -or $response -eq 'Y') {
        Write-Host "Removing existing virtual environment..." -ForegroundColor Yellow
        Remove-Item -Recurse -Force $VenvDir
    } else {
        Write-Host "Using existing virtual environment." -ForegroundColor Green
        Write-Host "To activate: $VenvDir\Scripts\Activate.ps1"
        exit 0
    }
}

Write-Host "Creating virtual environment at: $VenvDir"
& $PythonBin -m venv $VenvDir

Write-Host "Activating virtual environment..."
& "$VenvDir\Scripts\Activate.ps1"

Write-Host "Upgrading pip..."
python -m pip install --upgrade pip

# --- Install runtime dependencies ---
Write-Host "Installing dependencies from requirements.txt..."
$reqFile = Join-Path $SCRIPT_DIR "requirements.txt"
if (Test-Path $reqFile) {
    pip install -r $reqFile
} else {
    Write-Host "Warning: requirements.txt not found at $reqFile" -ForegroundColor Yellow
}

# --- Install build dependencies ---
Write-Host "Installing build dependencies from requirements-build.txt..."
$reqBuildFile = Join-Path $SCRIPT_DIR "requirements-build.txt"
if (Test-Path $reqBuildFile) {
    pip install -r $reqBuildFile
} else {
    Write-Host "Warning: requirements-build.txt not found at $reqBuildFile" -ForegroundColor Yellow
}

# --- Build preprocessing tools via CMake (build_preprocessing.sh is bash-only) ---
Write-Host "Building preprocessing tools..."

$activePython = (Get-Command python).Source
$buildDir = Join-Path $UDALES_ROOT "tools\preprocessing\build"

$cmakeArgs = @(
    "-G", "Ninja",
    "-S", (Join-Path $UDALES_ROOT "tools\preprocessing"),
    "-B", $buildDir,
    "-DCMAKE_POLICY_VERSION_MINIMUM=3.5",
    "-DPREPROCESSING_PYTHON_EXECUTABLE=$activePython",
    "-DCMAKE_C_COMPILER=gcc",
    "-DCMAKE_Fortran_COMPILER=gfortran"
)

# Check numpy and headers before enabling f2py targets
$hasNumpy  = (& python -c "import numpy" 2>$null; $LASTEXITCODE -eq 0)
$hasHeader = (& python -c "import pathlib, sysconfig; h = pathlib.Path(sysconfig.get_paths().get('include','')) / 'Python.h'; exit(0 if h.is_file() else 1)" 2>$null; $LASTEXITCODE -eq 0)

if ($BuildTarget -eq "preprocessing_tools") {
    if ($hasNumpy -and $hasHeader) {
        $cmakeArgs += "-DBUILD_PREPROCESSING_DIRECTSHORTWAVE_F2PY=ON"
        $cmakeArgs += "-DBUILD_PREPROCESSING_IBM_F2PY=ON"
    } else {
        if (-not $hasNumpy)  { Write-Host "WARNING: numpy not available; f2py targets will be disabled." -ForegroundColor Yellow }
        if (-not $hasHeader) { Write-Host "WARNING: Python.h not found; f2py targets will be disabled." -ForegroundColor Yellow }
    }
}

cmake @cmakeArgs
cmake --build $buildDir --target $BuildTarget

# --- Post-build summary ---
$view3dExe = Join-Path $buildDir "bin\view3d.exe"
if (Test-Path $view3dExe) {
    Write-Host "View3D executable available at $view3dExe"
    $legacyDir = Join-Path $UDALES_ROOT "tools\View3D\src"
    New-Item -ItemType Directory -Force -Path $legacyDir | Out-Null
    Copy-Item -Path $view3dExe -Destination (Join-Path $legacyDir "View3D.exe") -Force
    Write-Host "MATLAB compatibility copy available at $legacyDir\View3D.exe"
}
$dsSo = Get-ChildItem (Join-Path $UDALES_ROOT "tools\python\udprep") -Filter "directshortwave_f2py*.pyd" -ErrorAction SilentlyContinue
if ($dsSo) {
    Write-Host "directshortwave f2py module available at tools\python\udprep\"
} elseif ($BuildTarget -eq "preprocessing_tools") {
    Write-Host "WARNING: directshortwave f2py module missing. Check that GFortran and numpy are available." -ForegroundColor Yellow
}
$ibmSo = Get-ChildItem (Join-Path $UDALES_ROOT "tools\python\udprep") -Filter "ibm_preproc_f2py*.pyd" -ErrorAction SilentlyContinue
if ($ibmSo) {
    Write-Host "IBM preprocessing f2py module available at tools\python\udprep\"
} elseif ($BuildTarget -eq "preprocessing_tools") {
    Write-Host "WARNING: ibm_preproc f2py module missing. Check that GFortran and numpy are available." -ForegroundColor Yellow
}

Write-Host ""
Write-Host "==========================================" -ForegroundColor Green
Write-Host "Setup complete!" -ForegroundColor Green
Write-Host "==========================================" -ForegroundColor Green
Write-Host ""
Write-Host "To use the virtual environment:" -ForegroundColor Cyan
Write-Host "  1. Activate:   $VenvDir\Scripts\Activate.ps1"
Write-Host "  2. Set View3D threads: `$env:PREPROC_NCPU = `"8`""
Write-Host "  3. Run script: python tools\write_inputs.py [case_dir]"
Write-Host "  4. Deactivate: deactivate"
Write-Host ""
Write-Host "Example workflow:" -ForegroundColor Cyan
Write-Host "  cd $UDALES_ROOT"
Write-Host "  $VenvDir\Scripts\Activate.ps1"
Write-Host "  `$env:PREPROC_NCPU = `"8`""
Write-Host "  python tools\write_inputs.py examples\999"
Write-Host "  deactivate"
Write-Host ""
