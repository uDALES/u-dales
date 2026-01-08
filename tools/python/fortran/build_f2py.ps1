param(
    [string]$ModuleName = "directshortwave_f2py",
    [string]$Source = "tools/python/fortran/directShortwave_f2py.f90",
    [string]$Opt = "-O3",
    [string]$TargetDir = "tools/python/udprep"
)

$ErrorActionPreference = "Stop"

Write-Host "Building $ModuleName from $Source..."
python -m numpy.f2py -c -m $ModuleName $Source --opt="$Opt"
$pyds = Get-ChildItem -Path . -Filter "$ModuleName*.pyd" -ErrorAction SilentlyContinue
if ($pyds.Count -gt 0) {
    if (!(Test-Path $TargetDir)) {
        New-Item -ItemType Directory -Path $TargetDir | Out-Null
    }
    foreach ($pyd in $pyds) {
        Copy-Item -Path $pyd.FullName -Destination $TargetDir -Force
        Write-Host "Copied $($pyd.Name) to $TargetDir"
    }
} else {
    Write-Host "No .pyd files found to copy."
}

Write-Host "Done."
