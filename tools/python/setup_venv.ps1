# Setup script for Windows PowerShell
# Creates a Python virtual environment and installs all dependencies

$ErrorActionPreference = "Stop"

$SCRIPT_DIR = Split-Path -Parent $MyInvocation.MyCommand.Path
# Go to u-dales root (tools/python -> u-dales)
$UDALES_ROOT = (Get-Item (Join-Path $SCRIPT_DIR "..\..")).FullName
# Create venv in parent directory of u-dales
$VENV_DIR = Join-Path (Split-Path -Parent $UDALES_ROOT) "venv-udales"

Write-Host "=========================================" -ForegroundColor Cyan
Write-Host "Setting up Python virtual environment" -ForegroundColor Cyan
Write-Host "=========================================" -ForegroundColor Cyan
Write-Host "u-dales repository: $UDALES_ROOT" -ForegroundColor White
Write-Host "Virtual environment: $VENV_DIR" -ForegroundColor White
Write-Host ""

# Check if Python is available
try {
    $pythonVersion = python --version 2>&1
    Write-Host "Python version: $pythonVersion"
} catch {
    Write-Host "Error: Python is not installed or not in PATH" -ForegroundColor Red
    exit 1
}

# Create virtual environment if it doesn't exist
if (Test-Path $VENV_DIR) {
    Write-Host "Virtual environment already exists at: $VENV_DIR" -ForegroundColor Yellow
    $response = Read-Host "Do you want to recreate it? (y/N)"
    if ($response -eq 'y' -or $response -eq 'Y') {
        Write-Host "Removing existing virtual environment..." -ForegroundColor Yellow
        Remove-Item -Recurse -Force $VENV_DIR
    } else {
        Write-Host "Using existing virtual environment." -ForegroundColor Green
        Write-Host "To activate: ..\venv-udales\Scripts\Activate.ps1 (from u-dales root)" -ForegroundColor Green
        exit 0
    }
}

Write-Host "Creating virtual environment at: $VENV_DIR"
python -m venv $VENV_DIR

Write-Host "Activating virtual environment..."
& "$VENV_DIR\Scripts\Activate.ps1"

Write-Host "Upgrading pip..."
python -m pip install --upgrade pip

Write-Host "Installing dependencies from requirements.txt..."
$reqFile = Join-Path $SCRIPT_DIR "requirements.txt"
if (Test-Path $reqFile) {
    pip install -r $reqFile
} else {
    Write-Host "Warning: requirements.txt not found at $reqFile" -ForegroundColor Yellow
}

Write-Host ""
Write-Host "==========================================" -ForegroundColor Green
Write-Host "Setup complete!" -ForegroundColor Green
Write-Host "==========================================" -ForegroundColor Green
Write-Host ""
Write-Host "To use the virtual environment:" -ForegroundColor Cyan
Write-Host "  1. Activate:   ..\venv-udales\Scripts\Activate.ps1 (from u-dales root)" -ForegroundColor White
Write-Host "  2. Run script: python tools\write_inputs.py [config_dir]" -ForegroundColor White
Write-Host "  3. Deactivate: deactivate" -ForegroundColor White
Write-Host ""
Write-Host "Example workflow:" -ForegroundColor Cyan
Write-Host "  cd $UDALES_ROOT" -ForegroundColor White
Write-Host "  ..\venv-udales\Scripts\Activate.ps1" -ForegroundColor White
Write-Host "  python tools\write_inputs.py" -ForegroundColor White
Write-Host "  deactivate" -ForegroundColor White
Write-Host ""
