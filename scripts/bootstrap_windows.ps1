& {
    $ErrorActionPreference = 'Stop'
    $downloads = Join-Path $env:USERPROFILE 'Downloads'
    $package = Get-ChildItem -LiteralPath $downloads -File -Filter 'STRESS_COMBINATION_LOCAL_v1*.zip' |
        Sort-Object LastWriteTimeUtc -Descending | Select-Object -First 1
    if (-not $package) { throw 'Save STRESS_COMBINATION_LOCAL_v1.zip in Downloads first.' }
    if (-not (Test-Path -LiteralPath 'D:\' -PathType Container)) { throw 'D: drive is unavailable.' }
    $project = 'D:\projects\STRESS_COMBINATION_PREDICTION'
    $staging = Join-Path (Join-Path $project 'package_staging') ([guid]::NewGuid().ToString('N'))
    New-Item -ItemType Directory -Path $staging -Force | Out-Null
    Expand-Archive -LiteralPath $package.FullName -DestinationPath $staging
    $runner = Join-Path $staging 'stress-combination-prediction\run.ps1'
    if (-not (Test-Path -LiteralPath $runner -PathType Leaf)) { throw 'The extracted package has no run.ps1.' }
    & powershell.exe -NoProfile -ExecutionPolicy Bypass -File $runner -ProjectRoot $project
    if ($LASTEXITCODE -ne 0) {
        Write-Host 'Validation failed. Return the RESULT_STRESS_COMBINATION_VALIDATION ZIP printed above.'
    }
}
