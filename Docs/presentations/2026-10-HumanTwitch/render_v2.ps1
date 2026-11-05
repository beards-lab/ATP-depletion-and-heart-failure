# Render every slide of a .pptx to PNG using installed PowerPoint (COM). Usage:
#   powershell.exe -NoProfile -ExecutionPolicy Bypass -File render_v2.ps1 C:\path\deck.pptx [C:\out\dir]
param([string]$Pptx, [string]$OutDir = "")
$Pptx = (Resolve-Path $Pptx).Path
if (-not $OutDir) { $OutDir = Join-Path (Split-Path $Pptx) "render_v2" }
New-Item -ItemType Directory -Force -Path $OutDir | Out-Null
Get-ChildItem $OutDir -Filter "slide-*.png" | Remove-Item
$app = New-Object -ComObject PowerPoint.Application
try {
  # ReadOnly, Untitled, WithWindow=false
  $p = $app.Presentations.Open($Pptx, -1, 0, 0)
  $i = 1
  foreach ($s in $p.Slides) { $s.Export((Join-Path $OutDir ("slide-{0}.png" -f $i)), "PNG", 1600, 900); $i++ }
  $p.Close()
  Write-Output "render_v2ed $($i-1) slides to $OutDir"
} finally {
  if ($app.Presentations.Count -eq 0) { $app.Quit() }
  [System.Runtime.InteropServices.Marshal]::ReleaseComObject($app) | Out-Null
}
