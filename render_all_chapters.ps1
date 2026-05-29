$env:PATH = "C:\Windows\System32;C:\Program Files\Quarto\bin;$env:PATH"
$logDir = "render_logs"
if (-not (Test-Path $logDir)) { New-Item -ItemType Directory -Path $logDir | Out-Null }

$summary = @()
foreach ($f in Get-ChildItem ch??_*.qmd | Where-Object Name -NotLike "*.bak" | Sort-Object Name) {
    $name = $f.Name
    $base = [System.IO.Path]::GetFileNameWithoutExtension($name)
    $log  = Join-Path $logDir "$base.log"
    $start = Get-Date
    Write-Host ">>> $name (started $start)"
    $output = & quarto render $name --to html 2>&1
    $exit = $LASTEXITCODE
    $duration = ((Get-Date) - $start).TotalSeconds
    $output | Out-File -FilePath $log -Encoding utf8
    $status = if ($exit -eq 0) { "PASS" } else { "FAIL" }
    $line = "{0,-50} {1,-4} {2,7:N1}s   log: {3}" -f $name, $status, $duration, $log
    Write-Host $line
    $summary += $line
}

Write-Host ""
Write-Host "=== Summary ==="
$summary | ForEach-Object { Write-Host $_ }
