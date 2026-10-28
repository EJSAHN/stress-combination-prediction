[CmdletBinding()]
param(
    [string]$ProjectRoot = 'D:\projects\STRESS_COMBINATION_PREDICTION',
    [string]$PythonExe = '',
    [string]$ResultDirectory = (Join-Path $env:USERPROFILE 'Downloads')
)

# Windows PowerShell 5.1. Reuse an existing interpreter; never install anything.
Set-StrictMode -Version 2.0
$ErrorActionPreference = 'Stop'
$packageRoot = $PSScriptRoot
$runId = (Get-Date -Format 'yyyyMMdd_HHmmss') + '_' + [guid]::NewGuid().ToString('N').Substring(0, 8)
$runRoot = $null
$logPath = $null
$resultZip = $null
$exitCode = 1
$chosenPython = ''
$savedEnvironment = @{}
$utf8 = New-Object System.Text.UTF8Encoding($false)
$savedOutputEncoding = $OutputEncoding
$savedConsoleOutputEncoding = [Console]::OutputEncoding

function Write-RunLog([string]$Message) {
    $line = '[' + (Get-Date -Format 'yyyy-MM-dd HH:mm:ss') + '] ' + $Message
    Write-Host $line
    if ($null -ne $script:logPath) {
        [System.IO.File]::AppendAllText($script:logPath, $line + [Environment]::NewLine, $script:utf8)
    }
}

function Set-RunEnvironment([string]$Name, [string]$Value) {
    if (-not $script:savedEnvironment.ContainsKey($Name)) {
        $script:savedEnvironment[$Name] = [Environment]::GetEnvironmentVariable($Name, 'Process')
    }
    [Environment]::SetEnvironmentVariable($Name, $Value, 'Process')
}

try {
    $OutputEncoding = $utf8
    [Console]::OutputEncoding = $utf8
    if (-not [System.IO.Path]::IsPathRooted($ProjectRoot)) {
        throw 'ProjectRoot must be an absolute path.'
    }
    $ProjectRoot = [System.IO.Path]::GetFullPath($ProjectRoot)
    $driveRoot = [System.IO.Path]::GetPathRoot($ProjectRoot)
    if (-not (Test-Path -LiteralPath $driveRoot -PathType Container)) {
        throw ('Project drive is unavailable: ' + $driveRoot)
    }
    $runRoot = Join-Path (Join-Path $ProjectRoot 'runs') $runId
    if (Test-Path -LiteralPath $runRoot) { throw 'Run directory already exists.' }
    $logDirectory = Join-Path $runRoot 'logs'
    $tempDirectory = Join-Path $runRoot 'tmp'
    $cacheDirectory = Join-Path $runRoot 'cache'
    $matplotlibDirectory = Join-Path $cacheDirectory 'matplotlib'
    foreach ($directory in @($logDirectory, $tempDirectory, $cacheDirectory, $matplotlibDirectory)) {
        New-Item -ItemType Directory -Path $directory -Force | Out-Null
    }
    $logPath = Join-Path $logDirectory 'launcher.log'
    $resultZip = Join-Path $ResultDirectory ('RESULT_STRESS_COMBINATION_VALIDATION_' + $runId + '.zip')
    Write-RunLog ('RUN_ID=' + $runId)
    Write-RunLog ('PROJECT_ROOT=' + $ProjectRoot)
    Write-RunLog ('PACKAGE_ROOT=' + $packageRoot)
    Write-RunLog ('RUN_ROOT=' + $runRoot)
    Write-RunLog 'MODE=EXISTING_PYTHON_ONLY; no environment creation, installation, or data download.'

    Set-RunEnvironment 'TEMP' $tempDirectory
    Set-RunEnvironment 'TMP' $tempDirectory
    Set-RunEnvironment 'MPLCONFIGDIR' $matplotlibDirectory
    Set-RunEnvironment 'XDG_CACHE_HOME' $cacheDirectory
    Set-RunEnvironment 'PYTHONDONTWRITEBYTECODE' '1'
    Set-RunEnvironment 'PYTHONUNBUFFERED' '1'
    Set-RunEnvironment 'PYTHONIOENCODING' 'utf-8'
    Set-RunEnvironment 'STRESS_PIPELINE_LAUNCHER' 'powershell'
    foreach ($threadVariable in @('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS')) {
        Set-RunEnvironment $threadVariable '1'
    }
    foreach ($requiredFile in @('package_manifest.json', 'scripts\run_pipeline.py', 'scripts\local_run.py')) {
        if (-not (Test-Path -LiteralPath (Join-Path $packageRoot $requiredFile) -PathType Leaf)) {
            throw ('Incomplete package: ' + $requiredFile)
        }
    }

    $probePath = Join-Path $logDirectory 'check_existing_environment.py'
    $probeCode = @'
import importlib
import json
import platform
import sys
names = ('numpy', 'pandas', 'matplotlib', 'openpyxl', 'PIL')
if sys.version_info < (3, 10):
    raise RuntimeError('An existing Python 3.10 or later is required.')
versions = {name: importlib.import_module(name).__version__ for name in names}
if sys.maxsize <= 2**32:
    raise RuntimeError('A 64-bit existing Python is required.')
print(json.dumps({'ok': True, 'python': platform.python_version(), 'executable': sys.executable, 'bits': 64, 'packages': versions}))
'@
    [System.IO.File]::WriteAllText($probePath, $probeCode, $utf8)
    $candidates = New-Object 'System.Collections.Generic.List[string]'
    if ($PythonExe) {
        $candidates.Add($PythonExe)
    } else {
        $candidates.Add((Join-Path $env:USERPROFILE 'anaconda3\python.exe'))
        $candidates.Add((Join-Path $env:USERPROFILE 'miniconda3\python.exe'))
        if ($env:CONDA_PREFIX) { $candidates.Add((Join-Path $env:CONDA_PREFIX 'python.exe')) }
        foreach ($name in @('python.exe', 'python3.exe')) {
            $application = Get-Command -Name $name -CommandType Application -ErrorAction SilentlyContinue | Select-Object -First 1
            if ($null -ne $application) { $candidates.Add($application.Source) }
        }
    }
    $seenCandidates = @{}
    foreach ($candidate in $candidates) {
        if ($seenCandidates.ContainsKey($candidate)) { continue }
        $seenCandidates[$candidate] = $true
        if ($candidate -match '\\Microsoft\\WindowsApps\\') { continue }
        if (-not (Test-Path -LiteralPath $candidate -PathType Leaf)) { continue }
        Write-RunLog ('CHECK_EXISTING=' + $candidate)
        $priorPreference = $ErrorActionPreference
        try {
            $ErrorActionPreference = 'Continue'
            $probeOutput = & $candidate -B $probePath 2>&1
            $probeExit = $LASTEXITCODE
        } finally {
            $ErrorActionPreference = $priorPreference
        }
        foreach ($line in $probeOutput) { Write-RunLog ('PROBE: ' + $line.ToString()) }
        if ($probeExit -eq 0) { $chosenPython = $candidate; break }
    }
    if (-not $chosenPython) {
        throw 'No existing 64-bit Python has the required packages. Nothing was installed. Use -PythonExe with an existing suitable interpreter.'
    }
    Write-RunLog ('REUSED_PYTHON=' + $chosenPython)
    $launcherEnvironment = [ordered]@{
        PowerShellVersion = $PSVersionTable.PSVersion.ToString()
        Python = $chosenPython
        ProjectRoot = $ProjectRoot
        PackageRoot = $packageRoot
        RunRoot = $runRoot
        ExistingEnvironmentOnly = $true
        NativeLauncher = 'Windows PowerShell'
        TemporaryDirectory = $tempDirectory
        MatplotlibCache = $matplotlibDirectory
    }
    [System.IO.File]::WriteAllText((Join-Path $logDirectory 'launcher_environment.json'), ($launcherEnvironment | ConvertTo-Json -Depth 5), $utf8)
    $helper = Join-Path $packageRoot 'scripts\local_run.py'
    $priorPreference = $ErrorActionPreference
    try {
        $ErrorActionPreference = 'Continue'
        & $chosenPython -B -u $helper --project-root $ProjectRoot --run-root $runRoot --result-directory $ResultDirectory 2>&1 |
            ForEach-Object { Write-RunLog ($_.ToString()) }
        $exitCode = $LASTEXITCODE
    } finally {
        $ErrorActionPreference = $priorPreference
    }
    if ($exitCode -ne 0) { throw ('Local pipeline validation failed, exit code ' + $exitCode) }
    if (-not (Test-Path -LiteralPath $resultZip -PathType Leaf)) { throw 'The result ZIP was not created.' }
    Write-RunLog 'LOCAL_VALIDATION=PASS'
} catch {
    $exitCode = 1
    Write-RunLog ('FAILED=' + $_.Exception.Message)
    if ($null -ne $runRoot -and (Test-Path -LiteralPath $runRoot -PathType Container)) {
        $failure = [ordered]@{ status = 'FAIL'; message = $_.Exception.Message; run_id = $runId; python = $chosenPython }
        [System.IO.File]::WriteAllText((Join-Path $runRoot 'launcher_failure.json'), ($failure | ConvertTo-Json -Depth 4), $utf8)
        if ($null -ne $resultZip -and -not (Test-Path -LiteralPath $resultZip)) {
            try {
                New-Item -ItemType Directory -Path $ResultDirectory -Force | Out-Null
                $fallbackItems = @((Join-Path $runRoot 'logs'), (Join-Path $runRoot 'launcher_failure.json'))
                Compress-Archive -LiteralPath $fallbackItems -DestinationPath $resultZip -CompressionLevel Optimal
                Write-RunLog 'A diagnostic ZIP was created without original input data.'
            } catch {
                Write-RunLog ('DIAGNOSTIC_BUNDLE_FAILED=' + $_.Exception.Message)
            }
        }
    }
} finally {
    foreach ($name in $savedEnvironment.Keys) {
        [Environment]::SetEnvironmentVariable($name, $savedEnvironment[$name], 'Process')
    }
    $OutputEncoding = $savedOutputEncoding
    [Console]::OutputEncoding = $savedConsoleOutputEncoding
}
if ($null -ne $resultZip -and (Test-Path -LiteralPath $resultZip -PathType Leaf)) {
    Write-Host ('UPLOAD_THIS_ZIP=' + $resultZip)
}
if ($null -ne $runRoot) { Write-Host ('RUN_ROOT=' + $runRoot) }
Write-Host ('EXIT_CODE=' + $exitCode)
exit $exitCode
