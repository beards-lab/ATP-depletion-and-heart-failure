# build_status_deck.ps1 -- generate the 2026-09-16 status deck via PowerPoint COM.
#
# WHY COM. This box has no node and no python, and MATLAB Report Generator is
# listed in the toolbox folder but is not actually installed (open() throws
# "MATLAB Report Generator is not installed"), so pptxgenjs and mlreportgen.ppt
# are both unavailable. PowerPoint itself is installed, so COM it is.
#
# Content source: Docs/experiments/status-2026-09-16.md
#
#   powershell -ExecutionPolicy Bypass -File Docs/presentations/build_status_deck.ps1

$ErrorActionPreference = 'Stop'
$root = 'C:\home\git\ATP-wt-floored'
$out  = Join-Path $root 'Docs\presentations\status-2026-09-16.pptx'
$figR = Join-Path $root 'Analyses\RecoveryOscillation\results'
$figD = Join-Path $root 'Analyses\DataPooling\results'

# --- palette: cardiac. Cherry dominates, navy supports, warm grey recedes. ----
$CHERRY = 0x110099   # BGR! 990011
$NAVY   = 0x7E3C2F   # 2F3C7E
$INK    = 0x1A1A1A
$PAPER  = 0x00FFFFFF
$OFFW   = 0xF5F6FC   # FCF6F5
$GREY   = 0x6B6B6B
$LGREY  = 0xE8E8E8
$GREEN  = 0x2D5F2C   # 2C5F2D forest, for "good"
$AMBER  = 0x007AB8   # B87A00

function Hex2Ole([int]$r,[int]$g,[int]$b){ return $r + ($g -shl 8) + ($b -shl 16) }

$pp = New-Object -ComObject PowerPoint.Application
$pres = $pp.Presentations.Add(0)
# Set the canvas explicitly rather than via a PpSlideSizeType enum. The enum is
# a trap: ppSlideSizeOnScreen16x9 is 15, while 3 is ppSlideSizeA4Paper -- and
# every coordinate below assumes 960 x 540 pt, so the wrong enum silently pushes
# a third of each slide off the right edge instead of erroring.
$pres.PageSetup.SlideWidth  = 960
$pres.PageSetup.SlideHeight = 540
$W = 960.0; $H = 540.0; $M = 52.0
$BLANK = $pres.SlideMaster.CustomLayouts.Item(7)

function New-Slide([int]$bg = -1) {
    $s = $pres.Slides.AddSlide($pres.Slides.Count + 1, $BLANK)
    if ($bg -ge 0) {
        $s.FollowMasterBackground = 0
        $s.Background.Fill.Visible = -1
        $s.Background.Fill.Solid()
        $s.Background.Fill.ForeColor.RGB = $bg
    }
    return $s
}

function Add-Text($slide, [double]$x, [double]$y, [double]$w, [double]$h, [string]$text, [double]$size, [int]$color, [bool]$bold, [int]$align) {
    # align: 1 left, 2 centre, 3 right
    $tb = $slide.Shapes.AddTextbox(1, $x, $y, $w, $h)
    $tb.TextFrame.WordWrap = -1
    $tb.TextFrame.MarginLeft = 0; $tb.TextFrame.MarginRight = 0
    $tb.TextFrame.MarginTop = 0;  $tb.TextFrame.MarginBottom = 0
    $tb.TextFrame.TextRange.Text = $text
    $tb.TextFrame.TextRange.Font.Name = 'Calibri'
    $tb.TextFrame.TextRange.Font.Size = $size
    $tb.TextFrame.TextRange.Font.Color.RGB = $color
    if ($bold) { $tb.TextFrame.TextRange.Font.Bold = -1 } else { $tb.TextFrame.TextRange.Font.Bold = 0 }
    $tb.TextFrame.TextRange.ParagraphFormat.Alignment = $align
    return $tb
}

function Add-Title($slide, [string]$text, [int]$color) {
    Add-Text $slide $M 44 ($W - 2*$M) 60 $text 34 $color $true 1 | Out-Null
}

function Add-Card($slide, [double]$x, [double]$y, [double]$w, [double]$h, [int]$fill) {
    $c = $slide.Shapes.AddShape(5, $x, $y, $w, $h)   # 5 = rounded rectangle
    $c.Fill.Solid(); $c.Fill.ForeColor.RGB = $fill
    $c.Line.Visible = 0
    $c.Shadow.Visible = -1
    $c.Shadow.Blur = 8
    $c.Shadow.OffsetX = 0; $c.Shadow.OffsetY = 2
    $c.Shadow.Transparency = 0.82
    $c.Adjustments.Item(1) = 0.04
    return $c
}

function Add-Bullets($slide, [double]$x, [double]$y, [double]$w, [double]$h, [string[]]$items, [double]$size, [int]$color) {
    $tb = $slide.Shapes.AddTextbox(1, $x, $y, $w, $h)
    $tb.TextFrame.WordWrap = -1
    $tb.TextFrame.MarginLeft = 0; $tb.TextFrame.MarginRight = 0
    $tb.TextFrame.MarginTop = 0;  $tb.TextFrame.MarginBottom = 0
    $tb.TextFrame.TextRange.Text = ($items -join "`r")
    $tr = $tb.TextFrame.TextRange
    $tr.Font.Name = 'Calibri'; $tr.Font.Size = $size; $tr.Font.Color.RGB = $color
    $tr.ParagraphFormat.Bullet.Visible = -1
    $tr.ParagraphFormat.Bullet.Character = 8226
    $tr.ParagraphFormat.SpaceAfter = 7
    $tr.ParagraphFormat.SpaceWithin = 0.92
    return $tb
}

function Add-Pic($slide, [string]$path, [double]$x, [double]$y, [double]$maxW, [double]$maxH) {
    if (-not (Test-Path $path)) { Write-Warning "missing figure: $path"; return $null }
    $p = $slide.Shapes.AddPicture($path, 0, -1, $x, $y, -1, -1)
    $sc = [Math]::Min($maxW / $p.Width, $maxH / $p.Height)
    $p.LockAspectRatio = -1
    $p.Width = $p.Width * $sc
    $p.Left = $x + ($maxW - $p.Width) / 2
    $p.Top  = $y + ($maxH - $p.Height) / 2
    return $p
}

function Add-Stat($slide, [double]$x, [double]$y, [double]$w, [string]$big, [string]$label, [int]$color) {
    Add-Text $slide $x $y $w 62 $big 48 $color $true 2 | Out-Null
    Add-Text $slide $x ($y + 60) $w 44 $label 12 $GREY $false 2 | Out-Null
}

function Add-Note($slide, $text) {
    # blank layouts may not expose a notes placeholder; never fail the build for a note
    try { $slide.NotesPage.Shapes.Item(2).TextFrame.TextRange.Text = $text } catch { }
}

# =============================================================================
# 1  TITLE
# =============================================================================
$s = New-Slide $INK
Add-Text $s $M 170 ($W - 2*$M) 120 "Where the cardiac`ncross-bridge fit stands" 46 $PAPER $true 1 | Out-Null
Add-Text $s $M 320 ($W - 2*$M) 30 "A three-preparation pooled target, a deadband objective, and what is left to fix" 17 0x00C8C8C8 $false 1 | Out-Null
Add-Text $s $M 420 ($W - 2*$M) 24 "F. Jezek   -   16 September 2026   -   8 mM pooled target + 2 mM hypothesis ladder" 12 $GREY $false 1 | Out-Null
Add-Note $s "Every number in this deck is re-evaluated from disk, never read off a snapshot header. The project has previously ranked a campaign on a header that turned out to be fake."

# =============================================================================
# 2  HEADLINE NUMBERS
# =============================================================================
$s = New-Slide $PAPER
Add-Title $s "Where we are, in five numbers" $INK
$cw = ($W - 2*$M - 3*22) / 4
Add-Card $s $M 150 $cw 150 $OFFW | Out-Null
Add-Stat $s $M 172 $cw "0.512" "8 mM objective (floored)" $CHERRY
Add-Card $s ($M + $cw + 22) 150 $cw 150 $OFFW | Out-Null
Add-Stat $s ($M + $cw + 22) 172 $cw "-47.7%" "from the 0.980 seed, 15 h" $CHERRY
Add-Card $s ($M + 2*($cw + 22)) 150 $cw 150 $OFFW | Out-Null
Add-Stat $s ($M + 2*($cw + 22)) 172 $cw "14 / 21" "features inside 1 prep-width" $CHERRY
Add-Card $s ($M + 3*($cw + 22)) 150 $cw 150 $OFFW | Out-Null
Add-Stat $s ($M + 3*($cw + 22)) 172 $cw "0.0%" "of the gain went to F-V" $CHERRY

Add-Card $s $M 330 ($W - 2*$M) 150 $OFFW | Out-Null
Add-Text $s ($M+28) 352 ($W - 2*$M - 56) 24 "Two things moved that had resisted several campaigns" 17 $NAVY $true 1 | Out-Null
Add-Bullets $s ($M+28) 384 ($W - 2*$M - 56) 90 @(
 "The absolute force LEVEL is now on the pool: 71.5 kPa vs 69.3 (+3%, z = 0.20). Every previous generation sat at +15%, carrying 03/27's own scale.",
 "Timing came inside 2 sd: peak1_t, vall_t and peak1_dSL were all beyond 2 sd on the previous core."
) 13 $INK | Out-Null
Add-Note $s "The last six weeks were mostly about the objective, not the model. Five re-weightings, each cost-neutral by construction, so absolute numbers reset every generation. Only within a campaign is a drop a fit improvement."

# =============================================================================
# 3  THE OBJECTIVE LINEAGE
# =============================================================================
$s = New-Slide $PAPER
Add-Title $s "Five generations of objective -- why the numbers reset" $INK
Add-Text $s $M 106 ($W-2*$M) 24 "Each seed is the previous optimum reweighted, with lambda chosen for cost-neutrality." 14 $GREY $false 1 | Out-Null

$rows = @(
 @('spreadw','weights = measured between-prep spread','3.85 -> 0.98'),
 @('estiff125','structural exponent 0.965 -> 1.25','staircase 5.16 -> 1.60'),
 @('pooled8','target -> 3-preparation pooled mean','4.43 -> 2.70'),
 @('rigaware','rsK80 restricted to cycles 1-4','14 h, -10%'),
 @('floored','one-prep-width deadband on every term','2.41 -> 0.908'),
 @('lean5','cycle 5 back + lean pool + coolDownLS out','0.980 -> 0.512')
)
$y = 150
foreach ($r in $rows) {
    $isLast = ($r[0] -eq 'lean5')
    $fill = if ($isLast) { Hex2Ole 0xF2 0xE2 0xE4 } else { $OFFW }
    Add-Card $s $M $y ($W - 2*$M) 48 $fill | Out-Null
    $tc = if ($isLast) { $CHERRY } else { $NAVY }
    Add-Text $s ($M+22) ($y+14) 130 22 $r[0] 15 $tc $true 1 | Out-Null
    Add-Text $s ($M+165) ($y+14) 470 22 $r[1] 13 $INK $false 1 | Out-Null
    Add-Text $s ($M+645) ($y+14) 190 22 $r[2] 13 $tc $true 3 | Out-Null
    $y += 56
}
Add-Text $s $M 496 ($W-2*$M) 24 "Only within a campaign is a drop a fit improvement. Across generations the yardstick itself changed." 13 $CHERRY $true 1 | Out-Null
Add-Note $s "This is the single most common way to misread the project's history: treating 0.98 and 0.51 as comparable to 3.85. They are not."

# =============================================================================
# 4  POOLING: OLD vs NEW TARGET
# =============================================================================
$s = New-Slide $PAPER
Add-Title $s "The fit target: one preparation, or three?" $INK
$hw = ($W - 2*$M - 30) / 2
Add-Card $s $M 140 $hw 300 $OFFW | Out-Null
Add-Text $s ($M+24) 162 ($hw-48) 26 "BEFORE -- 03/27 alone" 17 $GREY $true 1 | Out-Null
Add-Bullets $s ($M+24) 200 ($hw-48) 230 @(
 "One preparation, in that fibre's raw kPa",
 "No error bar anywhere",
 "So every residual is unfalsifiable: a 15% miss and a 2% miss look the same to the objective, because nothing says how big 15% is",
 "Cross-section CV across three days: 22%"
) 13 $INK | Out-Null

Add-Card $s ($M+$hw+30) 140 $hw 300 (Hex2Ole 0xF2 0xE2 0xE4) | Out-Null
Add-Text $s ($M+$hw+54) 162 ($hw-48) 26 "NOW -- three preparations, pooled" 17 $CHERRY $true 1 | Out-Null
Add-Bullets $s ($M+$hw+54) 200 ($hw-48) 230 @(
 "03/27 M, 04/03 F, 04/10 M -- both ATP conditions",
 "A value AND a spread on every feature",
 "Every term becomes a chi-square: 1.0 = one preparation-width from the data",
 "Spread 29% -> 15%. Absolute tension goes from worst-determined quantity to among the best"
) 13 $INK | Out-Null
Add-Note $s "Tension is force divided by ESTIMATED cross-sectional area. Across the three days the pre-slack isometric reads 79.8 / 70.0 / 50.9 kPa. Almost none of that is biology; it is the area estimate."

# =============================================================================
# 5  POOLING: THE FOUR STEPS + THE TEST
# =============================================================================
$s = New-Slide $PAPER
Add-Title $s "How the pooling works -- and how we know it is allowed" $INK
$q = ($W - 2*$M - 3*18) / 4
$steps = @(
 @('1  SCALE','Divide each prep''s force features by its own Fss22. ONE divisor per preparation, shared by both ATP conditions -- so the 2/8 mM ratio stays a measurement, not a definition.'),
 @('2  RUNDOWN','Lift the run recorded SECOND in each session by phi = 0.452, on effective ACTIVATED time. Without it the ATP effect depends on running order.'),
 @('3  POOL','Mean of the per-preparation FEATURES -- never features of an averaged trace.'),
 @('4  SPREAD','s = between-prep relative SD, variance-moderated toward its class median (Smyth 2004), d0 = 2.58 -> 4.58 dof.')
)
for ($i=0; $i -lt 4; $i++) {
    $x = $M + $i*($q+18)
    Add-Card $s $x 140 $q 190 $OFFW | Out-Null
    Add-Text $s ($x+16) 158 ($q-32) 22 $steps[$i][0] 14 $CHERRY $true 1 | Out-Null
    Add-Text $s ($x+16) 186 ($q-32) 132 $steps[$i][1] 11 $INK $false 1 | Out-Null
}
Add-Card $s $M 352 ($W-2*$M) 132 (Hex2Ole 0xE9 0xEE 0xE9) | Out-Null
Add-Text $s ($M+28) 372 ($W-2*$M-56) 24 "One common scale is TESTED, not assumed" 16 $GREEN $true 1 | Out-Null
Add-Bullets $s ($M+28) 402 ($W-2*$M-56) 76 @(
 "PC1 holds 88.3% of the between-preparation log-variance, and its loading is flat (|cos| with all-ones = 0.955). A pure common scale is exactly a rank-1 matrix with a flat loading -- both conditions hold.",
 "Control group: ktr, t0, rsK, rsR2 are scale-free, so no normalisation can touch them. After normalising, the force features sit at 4-11% -- better than the best rate."
) 12 $INK | Out-Null
Add-Note $s "Which anchor barely matters: a paired bootstrap over 4000 draws cannot separate the top three (P(best) 64/18/15%). Normalising versus not is decisive; which divisor is not a decision n=3 can support."

# =============================================================================
# 6  WHERE 03/27 MISLED
# =============================================================================
$s = New-Slide $INK
Add-Text $s $M 90 ($W-2*$M) 40 "Pooling adds error bars, not a different target --" 26 0x00C8C8C8 $false 1 | Out-Null
Add-Text $s $M 132 ($W-2*$M) 44 "except in the one feature that mattered." 30 $PAPER $true 1 | Out-Null
$hw2 = ($W - 2*$M - 30)/2
Add-Card $s $M 216 $hw2 128 0x00332A2C | Out-Null
Add-Text $s ($M+24) 238 ($hw2-48) 46 "3% BELOW" 30 0x0088AAFF $true 2 | Out-Null
Add-Text $s ($M+24) 292 ($hw2-48) 40 "the pool's first re-stretch peak,`nrelative to its own isometric" 12 0x00C8C8C8 $false 2 | Out-Null
Add-Card $s ($M+$hw2+30) 216 $hw2 128 0x00332A2C | Out-Null
Add-Text $s ($M+$hw2+54) 238 ($hw2-48) 46 "12% ABOVE" 30 0x008888FF $true 2 | Out-Null
Add-Text $s ($M+$hw2+54) 292 ($hw2-48) 40 "03/27's first re-stretch peak,`nrelative to its own isometric" 12 0x00C8C8C8 $false 2 | Out-Null
Add-Text $s $M 378 ($W-2*$M) 60 "Fitting 03/27 alone was REWARDING the model's lengthening overshoot -- a known artefact." 19 $PAPER $true 1 | Out-Null
Add-Text $s $M 432 ($W-2*$M) 46 "Otherwise 03/27 is a fair stand-in: chi-2 = 0.39 against the pooled mean, where 0.67 is what three preparations scattering around their own mean must give. It is the most typical of the three." 13 0x00A0A0A0 $false 1 | Out-Null
Add-Note $s "The target moved by more than the model's own residual, so the switch was not cosmetic."

# =============================================================================
# 7  CHI-2 WEIGHTING
# =============================================================================
$s = New-Slide $PAPER
Add-Title $s "Step 1 -- make every term a chi-square" $INK
Add-Card $s $M 140 ($W-2*$M) 96 (Hex2Ole 0xF2 0xE2 0xE4) | Out-Null
Add-Text $s $M 166 ($W-2*$M) 46 "w  =  lambda - R / s^2" 36 $CHERRY $true 2 | Out-Null
Add-Bullets $s $M 262 ($W-2*$M) 180 @(
 "s = the MEASURED between-preparation spread of that feature; R = an ATP-relevance bonus; lambda = one global constant fixing cost-neutrality.",
 "The cost sums squared RELATIVE residuals, so dividing by s^2 turns every term into a chi-square: 1.0 = the model is one preparation-width from the data.",
 "Terms become directly comparable -- which they were not when the weights 0.2 ... 50 were set by hand.",
 "It immediately exposed two mis-weightings: peak1_dSL was over-weighted 5x (76% spread -- barely measurable), F-V over-weighted 25x."
) 14 $INK | Out-Null
Add-Text $s $M 462 ($W-2*$M) 24 "Diagnostic that falls out for free: how far is the model from the data, in units the data itself sets?" 13 $NAVY $true 1 | Out-Null
Add-Note $s "ATP relevance is measured and plotted but deliberately does not dominate the weights."

# =============================================================================
# 8  THE FAILURE THAT WEIGHTING COULD NOT FIX
# =============================================================================
$s = New-Slide $PAPER
Add-Title $s "Step 2 -- the failure weighting alone could not fix" $INK
Add-Card $s $M 140 ($W-2*$M) 126 (Hex2Ole 0xF2 0xE2 0xE4) | Out-Null
Add-Text $s ($M+28) 160 ($W-2*$M-56) 90 "The rigaware campaign: force-velocity held 40% of the objective at chi-2 = 0.10 -- a tenth of one preparation-width -- and absorbed 110% of a 14-hour gain, while every feature the campaign existed to fix went backwards." 17 $CHERRY $true 1 | Out-Null
Add-Bullets $s $M 292 ($W-2*$M) 160 @(
 "An objective is a fixed budget, and a search spends it where descent is cheapest -- reliably whichever term is ALREADY best determined, because it is the smoothest.",
 "Re-pinning the weight cannot fix this, and it was tried: FV_SHARE 0.70 -> 0.40 is what produced that run.",
 "Any nonzero weight on a smooth, nearly-solved term is still the cheapest descent available.",
 "Strip F-V from that campaign's result and the fit is 1.8% WORSE."
) 14 $INK | Out-Null
Add-Text $s $M 470 ($W-2*$M) 26 "The fix has to change whether a solved term pays at all -- not how much it is worth." 15 $NAVY $true 1 | Out-Null
Add-Note $s "This is the pivotal diagnosis in the whole objective story."

# =============================================================================
# 9  THE DEADBAND
# =============================================================================
$s = New-Slide $PAPER
Add-Title $s "Step 3 -- the deadband floor" $INK
$hw3 = ($W - 2*$M - 30)/2
Add-Card $s $M 138 $hw3 152 $OFFW | Out-Null
Add-Text $s ($M+22) 156 ($hw3-44) 22 "The mechanism -- one line" 14 $NAVY $true 1 | Out-Null
# Built as a variable first: PowerShell's command-argument parser trips over a
# bare '(' inside an inline argument string.
$codeTxt = 'c = max(0, c - floor_raw);' + [char]13 + 'c = c + nNaN * NAN_COST;'
$code = Add-Text $s ($M+22) 186 ($hw3-44) 90 $codeTxt 13 $INK $false 1
$code.TextFrame.TextRange.Font.Name = 'Consolas'
Add-Text $s ($M+22) 254 ($hw3-44) 30 "A DEADBAND, not a deletion. NaN penalty applied AFTER." 11 $GREY $false 1 | Out-Null

Add-Card $s ($M+$hw3+30) 138 $hw3 152 $OFFW | Out-Null
Add-Text $s ($M+$hw3+52) 156 ($hw3-44) 22 "How the value is set" 14 $NAVY $true 1 | Out-Null
Add-Text $s ($M+$hw3+52) 184 ($hw3-44) 30 "fl  =  FLOOR_CHI2 - nEl - s^2" 17 $CHERRY $true 1 | Out-Null
Add-Text $s ($M+$hw3+52) 220 ($hw3-44) 64 "ONE number was chosen: FLOOR_CHI2 = 1.0. Every per-term value is derived. steady: 1.0 x 5 x 0.157604^2 = 0.124196 -- exactly the fl: in the seed." 11 $GREY $false 1 | Out-Null

Add-Card $s $M 306 ($W-2*$M) 96 (Hex2Ole 0xE9 0xEE 0xE9) | Out-Null
Add-Text $s ($M+28) 324 ($W-2*$M-56) 24 "It reads as: being within one between-preparation SD of the pooled mean is FREE." 16 $GREEN $true 1 | Out-Null
Add-Text $s ($M+28) 352 ($W-2*$M-56) 40 "A floored term still DEFENDS itself -- it starts paying the moment it drifts back out. Effect: objective 2.4088 -> 0.9081, i.e. 62% of it was residual inside one preparation-width. Composition flipped: F-V 40% -> 0%, re-stretch block 28% -> 48%." 12 $INK $false 1 | Out-Null
Add-Note $s "FLOOR_CHI2 = 1.0 is itself a judgement call, not a measurement. And s is measured at n=3, so each weight is known only to a factor of about 18 -- the floor inherits that uncertainty."

# =============================================================================
# 10  THE FLOOR'S HONEST CAVEAT
# =============================================================================
$s = New-Slide $PAPER
Add-Title $s "The floor stops a solved term being bought -- not sold" $CHERRY
Add-Pic $s (Join-Path $figR 'lean5_fv.png') $M 130 ($W-2*$M) 230 | Out-Null
Add-Card $s $M 378 ($W-2*$M) 110 (Hex2Ole 0xF2 0xE2 0xE4) | Out-Null
Add-Text $s ($M+28) 394 ($W-2*$M-56) 22 "F-V raw residual grew 16x and its worst point went 0.5 -> 2.1 sd, while the scored cost stayed at exactly zero." 15 $CHERRY $true 1 | Out-Null
Add-Text $s ($M+28) 422 ($W-2*$M-56) 58 "Seed: raw 0.000535 = 1.0% of the deadband, worst point 0.48 sd.   Best: raw 0.008819 = 17.1% of the deadband, worst point 2.06 sd.`nWhy so wide: s = 11.3% there is the disagreement between two LABS -- the loosest uncertainty F-V has -- while the per-point weights use the between-CURVE variance. Fix: a per-term FLOOR_CHI2 (F-V ~ 0.25), or per-element floors. Not applied: it changes the objective." 11 $INK $false 1 | Out-Null
Add-Note $s "DataPooling section 6 found the in-preparation FV2 ramp lands 0.21 sd from this borrowed curve at 8 mM, so the target is trustworthy and 2 sd is a real drift."

# =============================================================================
# 11  CAMPAIGN OVERVIEW FIGURE
# =============================================================================
$s = New-Slide $PAPER
Add-Title $s "The 15-hour campaign: 0.980 -> 0.512" $INK
Add-Pic $s (Join-Path $figR 'lean5_overview.png') $M 118 ($W-2*$M) 330 | Out-Null
Add-Text $s $M 462 ($W-2*$M) 40 "Panel 4 is the control: F-V pays zero before and after, so the drop came from vall2_dy, peak1_y, peak1_dSL and rsK80 -- the re-stretch block the objective exists to fix. Panel 2 is the campaign's own question: slack cycle 5, back in the objective." 12 $GREY $false 1 | Out-Null
Add-Note $s "Profile a plateaued at round 6; profile c was still descending at round 13. Budget-matched at round 9, a is genuinely ahead, so its lead is the search and not extra evaluations."

# =============================================================================
# 12  GOOD / BAD / UGLY
# =============================================================================
$s = New-Slide $PAPER
Add-Title $s "How well do we fit 8 mM? The good, the bad, the ugly" $INK
$t = ($W - 2*$M - 2*20)/3
Add-Card $s $M 134 $t 330 (Hex2Ole 0xE9 0xEE 0xE9) | Out-Null
Add-Text $s ($M+20) 152 ($t-40) 24 "GOOD" 18 $GREEN $true 1 | Out-Null
Add-Text $s ($M+20) 178 ($t-40) 20 "at or inside the data's resolution" 11 $GREY $false 1 | Out-Null
Add-Bullets $s ($M+20) 206 ($t-40) 246 @(
 "Am, peak2, A, steady, rsA:  0.03-0.15 sd",
 "vall2_dy, t0, vall2_t, vall_y, rsK:  0.18-0.51 sd",
 "peak1_y 0.58 - ktr 0.82 - rsK80 0.98",
 "Force LEVEL finally on the pool: 71.5 vs 69.3 kPa (z = 0.20)",
 "Timing block came inside 2 sd"
) 11 $INK | Out-Null

Add-Card $s ($M+$t+20) 134 $t 330 $OFFW | Out-Null
Add-Text $s ($M+$t+40) 152 ($t-40) 24 "BAD" 18 $AMBER $true 1 | Out-Null
Add-Text $s ($M+$t+40) 178 ($t-40) 20 "coherent, not random" 11 $GREY $false 1 | Out-Null
Add-Bullets $s ($M+$t+40) 206 ($t-40) 246 @(
 "peak1_t 1.63 sd",
 "vall_t 1.50 sd",
 "peak1_dSL 1.37 sd",
 "rsT20 1.35 sd",
 "restretchSlopeStart 1.03 sd",
 "These are FOUR READINGS OF ONE THING: the shape and timing of the re-stretch transient. They move together -- which is why a 10-parameter draw struggles."
) 11 $INK | Out-Null

Add-Card $s ($M+2*($t+20)) 134 $t 330 (Hex2Ole 0xF2 0xE2 0xE4) | Out-Null
Add-Text $s ($M+2*($t+20)+20) 152 ($t-40) 24 "UGLY" 18 $CHERRY $true 1 | Out-Null
Add-Text $s ($M+2*($t+20)+20) 178 ($t-40) 20 "genuinely wrong" 11 $GREY $false 1 | Out-Null
Add-Bullets $s ($M+2*($t+20)+20) 206 ($t-40) 246 @(
 "rsT0 at 3.14 sd -- the ONLY feature beyond 2 sd, and it is UNSCORED, so an honest out-of-sample miss",
 "rsK80 per-cycle SHAPE: model climbs 51 -> 64 s-1 where the data is flat at ~47 then steps to 54",
 "Recovery too fast in segments 2-4: model plateaus by 360 ms, data still climbing at 420 ms",
 "F-V drifted to 2 sd -- invisible to the objective"
) 11 $INK | Out-Null
Add-Note $s "Putting cycle 5 back improved rsK80's LEVEL (chi-2 1.4 -> 1.15) but bought none of the shape."

# =============================================================================
# 13  TIMECOURSE FIGURE
# =============================================================================
$s = New-Slide $PAPER
Add-Title $s "Model vs pooled timecourse -- all five slack cycles" $INK
Add-Pic $s (Join-Path $figD 'target_timecourse.png') $M 116 ($W-2*$M) 336 | Out-Null
Add-Text $s $M 462 ($W-2*$M) 40 "Top: full segments, grey band = between-preparation SD. Bottom: the re-stretch recovery window. Segment 5 is the odd one -- it re-stretches to ML 1.000, not 1.101, and plateaus at 0.76 of reference. It is the only probe of how the recovery depends on the length it recovers at." 12 $GREY $false 1 | Out-Null
Add-Note $s "The visible defect: in segments 2-4 the red model curve reaches its plateau well before the data, which is still climbing at the right edge of the window."

# =============================================================================
# 14  FEATURE RESIDUALS
# =============================================================================
$s = New-Slide $PAPER
Add-Title $s "Every feature, in units of the data's own spread" $INK
Add-Pic $s (Join-Path $figD 'target_features.png') $M 122 ($W-2*$M) 312 | Out-Null
Add-Text $s $M 448 ($W-2*$M) 50 "Right panel: |model - target| / s. Fourteen of twenty-one features sit inside ONE preparation-width; only rsT0 is beyond two. On the previous core FOUR features were beyond 2 sd -- peak1_t, vall_t and peak1_dSL have all come inside." 12 $GREY $false 1 | Out-Null
Add-Note $s "rsT0 is not in the objective, so its 3.14 sd is an out-of-sample miss rather than a fitting failure."

# =============================================================================
# 15  IS THERE ROOM?
# =============================================================================
$s = New-Slide $PAPER
Add-Title $s "Is there room for improvement -- and is it realistic?" $INK
$hw4 = ($W - 2*$M - 26)/2
Add-Card $s $M 130 $hw4 250 $OFFW | Out-Null
Add-Text $s ($M+22) 148 ($hw4-44) 22 "We are near the limit" 16 $NAVY $true 1 | Out-Null
Add-Bullets $s ($M+22) 178 ($hw4-44) 190 @(
 "MechanismAblation: 35 mechanisms flipped one at a time, 39.6 h -- NOTHING, added or removed, beats the incumbent beyond noise (sigma = 0.055, measured)",
 "An 11-DOF identifiability ceiling: adding parameters does not add constraint",
 "14 of 21 features are inside one preparation-width -- you cannot fit below the noise",
 "6 lattice parameters carry only 2 identifiable combinations; 4 more measured inert"
) 11.5 $INK | Out-Null

Add-Card $s ($M+$hw4+26) 130 $hw4 250 (Hex2Ole 0xF2 0xE2 0xE4) | Out-Null
Add-Text $s ($M+$hw4+48) 148 ($hw4-44) 22 "Room remains" 16 $CHERRY $true 1 | Out-Null
Add-Bullets $s ($M+$hw4+48) 178 ($hw4-44) 190 @(
 "The residual is STRUCTURED, not random: one coherent defect across five features, plus a sign error in the rsK80 cycle shape. Structure means mechanism",
 "Two protocols are not scored AT ALL -- staircase and ktr. Free out-of-sample tests",
 "Two objective terms are unprincipled raw RMSEs -- 18% of the objective",
 "Searches STALL rather than converge: estiff moved 0.09% in four rounds and never again"
) 11.5 $INK | Out-Null

Add-Card $s $M 396 ($W-2*$M) 92 $INK | Out-Null
Add-Text $s ($M+28) 414 ($W-2*$M-56) 26 "The objective is close to its floor. The model has one structural defect left." 19 $PAPER $true 1 | Out-Null
Add-Text $s ($M+28) 446 ($W-2*$M-56) 30 "These are different statements, and conflating them is how a project concludes it is finished. More optimiser hours are NOT the lever -- that is the ablation's measured finding." 12 0x00B0B0B0 $false 1 | Out-Null
Add-Note $s "What would move things: more preparations, scoring the held-out protocols, fixing the unfloored terms, and a third ATP concentration."

# =============================================================================
# 16  2 mM PHENOTYPE
# =============================================================================
$s = New-Slide $PAPER
Add-Title $s "The 2 mM phenotype -- measured, not estimated" $INK
$q2 = ($W - 2*$M - 3*20)/4
$stats = @(@('x1.36','isometric force'), @('x0.55','ktr (CV 4%)'), @('x0.345','rsK'), @('x0.74','force x ktr'))
for ($i=0; $i -lt 4; $i++) {
    $x = $M + $i*($q2+20)
    Add-Card $s $x 132 $q2 128 $OFFW | Out-Null
    Add-Stat $s $x 152 $q2 $stats[$i][0] $stats[$i][1] $CHERRY
}
Add-Text $s $M 282 ($W-2*$M) 26 "Low ATP makes cardiac muscle STRONGER and SLOWER -- consistent across four preparations and both ATP running orders, once rundown is removed." 15 $INK $true 1 | Out-Null
Add-Card $s $M 322 ($W-2*$M) 132 (Hex2Ole 0xE9 0xEE 0xE9) | Out-Null
Add-Text $s ($M+28) 342 ($W-2*$M-56) 24 "The structural result matters more than the numbers" 16 $GREEN $true 1 | Out-Null
Add-Text $s ($M+28) 372 ($W-2*$M-56) 72 "Once each preparation is scaled by its own steady force, MOST FEATURES SHOW NO ATP EFFECT AT ALL -- A, peak1_y, vall_y, peak2, vall2_dy all sit at 1.0 +- 0.1.`n`nThe 2 mM phenotype is a force scale-up plus two rate slowdowns, NOT a change of waveform shape -- and the model reproduces that structure." 13 $INK $false 1 | Out-Null
Add-Note $s "This supersedes the earlier x1.18 figure, which was uncorrected 03/27 only."

# =============================================================================
# 17  THE GATE LADDER
# =============================================================================
$s = New-Slide $PAPER
Add-Title $s "Which rates are allowed to feel [ATP]?" $INK
Add-Text $s $M 106 ($W-2*$M) 22 "Every gate is normalised to 1 at 8 mM, so the high-ATP fit is invariant BY CONSTRUCTION -- asserted at |delta| = 0 before any fitting." 13 $GREY $false 1 | Out-Null
$hdr = @('rung','gates','DOF','2 mM cost','% of span')
$data = @(
 @('H0','no ATP effect at all','0','385.13','0%'),
 @('H1','k2, kah, kamh','3','8.20','99.5%'),
 @('H2','+ four SRX gates','7','8.20','99.5%'),
 @('H3','+ mechanics','10','7.64','99.7%'),
 @('H4','all 21','21','6.43','100%')
)
$colX = @(0, 110, 430, 530, 690); $colW = @(100, 310, 90, 150, 160)
for ($c=0; $c -lt 5; $c++) {
    Add-Text $s ($M+$colX[$c]) 146 $colW[$c] 20 $hdr[$c] 12 $GREY $true 1 | Out-Null
}
$y = 172
foreach ($r in $data) {
    $hl = ($r[0] -eq 'H1')
    if ($hl) { Add-Card $s ($M-10) ($y-6) ($W-2*$M+20) 34 (Hex2Ole 0xF2 0xE2 0xE4) | Out-Null }
    $tc = if ($hl) { $CHERRY } else { $INK }
    for ($c=0; $c -lt 5; $c++) {
        Add-Text $s ($M+$colX[$c]) $y $colW[$c] 22 $r[$c] 13 $tc $hl 1 | Out-Null
    }
    $y += 38
}
Add-Card $s $M 376 ($W-2*$M) 108 $OFFW | Out-Null
Add-Text $s ($M+28) 394 ($W-2*$M-56) 22 "Two gates do essentially everything -- and one of them is not defensible" 15 $NAVY $true 1 | Out-Null
Add-Text $s ($M+28) 422 ($W-2*$M-56) 56 "k2 (0.440) and kah (0.300) carry 99.5% of the whole span; the four SRX gates buy EXACTLY nothing, and everything beyond is compliance-bending, not pathway discovery.`nBut k2 is the canonical ATP step, while kah is hydrolysis of ALREADY-BOUND ATP -- free [ATP] does not appear in its rate law. The best-fitting rung is fitted, not mechanistic." 12 $INK $false 1 | Out-Null
Add-Note $s "H4 never converged: 21 gates x 4 steps needs about 84 evaluations per sweep and its budget allowed one. 6.43 is a floor, not an optimum."

# =============================================================================
# 18  DEFENSIBLE SETS  (recalibration)  -- placeholder filled by patch step
# =============================================================================
$s = New-Slide $PAPER
Add-Title $s "Defensible gate sets -- P1 and C4 recalibrated" $INK
Add-Text $s $M 104 ($W-2*$M) 22 "Recalibrating P1 and C4 onto the lean5 core did not reproduce the old ranking -- and finding out why exposed a hole in the low-ATP objective." 13 $GREY $false 1 | Out-Null

$hdr2 = @('id','gates','old core','lean5 core','fitted exponents')
$colX2 = @(0, 75, 360, 470, 600); $colW2 = @(70, 280, 105, 125, 250)
for ($c=0; $c -lt 5; $c++) { Add-Text $s ($M+$colX2[$c]) 140 $colW2[$c] 20 $hdr2[$c] 12 $GREY $true 1 | Out-Null }
$sets = @(
 @('C4','k2, ksrd2sr, kmsrd','8.26','13.50  (not conv)','k2 +0.300, kmsrd +0.300'),
 @('C0','k2 alone','11.62','13.95','k2 +0.300'),
 @('P1','k2, k_1','7.61','13.95','k2 +0.300, k_1 0.000'),
 @('H1','k2, kah, kamh','8.20','13.95','k2 +0.300, rest 0.000')
)
$y2 = 164
foreach ($r in $sets) {
    $hl = ($r[0] -eq 'P1' -or $r[0] -eq 'C4')
    if ($hl) { Add-Card $s ($M-10) ($y2-5) ($W-2*$M+20) 30 (Hex2Ole 0xF2 0xE2 0xE4) | Out-Null }
    $tc = if ($hl) { $CHERRY } else { $INK }
    for ($c=0; $c -lt 5; $c++) { Add-Text $s ($M+$colX2[$c]) $y2 $colW2[$c] 20 $r[$c] 12.5 $tc $hl 1 | Out-Null }
    $y2 += 34
}
Add-Text $s $M 306 ($W-2*$M) 20 "On the new core every set collapses to ONE lever -- k2 at +0.300 -- and every extra gate goes to exactly zero. C0 = P1 = H1 to four decimals." 13 $NAVY $true 1 | Out-Null

Add-Card $s $M 336 ($W-2*$M) 72 (Hex2Ole 0xF2 0xE2 0xE4) | Out-Null
Add-Text $s ($M+24) 348 ($W-2*$M-48) 20 "Do not read that as 'one gate is enough'. The fit had no rsK80." 14 $CHERRY $true 1 | Out-Null
Add-Text $s ($M+24) 372 ($W-2*$M-48) 30 "rsK80 is in NO curated .mat file -- not even the 8 mM one. On the data side it exists ONLY by injection from a pooled target. There was no pooled 2 mM target, so the low-ATP fit was scored WITHOUT its strongest ATP-carrying feature (z = 3.42), against one preparation, with 8 mM-sized deadbands." 10.5 $INK $false 1 | Out-Null

Add-Card $s $M 418 ($W-2*$M) 66 (Hex2Ole 0xE9 0xEE 0xE9) | Out-Null
Add-Text $s ($M+24) 430 ($W-2*$M-48) 20 "Fixed: a pooled 2 mM target now exists, and the re-run is in flight." 14 $GREEN $true 1 | Out-Null
Add-Text $s ($M+24) 454 ($W-2*$M-48) 26 "Same three preparations, same per-prep anchor, its OWN spread s2 -- and the ratio between the two frozen targets reproduces the measured effect: steady x1.331, ktr x0.544, rsK80 x0.382, rsK x0.338.  Port verified throughout: 8 mM bit-identical, |delta| = 0.000e+00." 10.5 $INK $false 1 | Out-Null
Add-Note $s ("P1 fits best on the OLD core and explains least: k_1 is the reverse stroke driven by Pi rebinding, " + `
             "so P1 says the 2 mM solution differs in Pi as well as ATP. Pi predicts LOWER force; the data show " + `
             "force x1.34 HIGHER. Port verified throughout: 8 mM bit-identical, |delta| = 0.000e+00.")

# =============================================================================
# 19  TWO HARD LIMITS
# =============================================================================
$s = New-Slide $PAPER
Add-Title $s "Two hard limits on the 2 mM side" $INK
$hw5 = ($W - 2*$M - 26)/2
Add-Card $s $M 138 $hw5 250 $OFFW | Out-Null
Add-Text $s ($M+22) 158 ($hw5-44) 22 "Km is not identifiable -- in principle" 15 $NAVY $true 1 | Out-Null
Add-Bullets $s ($M+22) 190 ($hw5-44) 186 @(
 "The power law is the K -> infinity limit of a Hill gate",
 "With data at only TWO [ATP], a gate is one number (its 2 mM factor), so Km follows in closed form. No refit is possible, and none would help",
 "The two forms diverge x0.4 on ktr/rsK by 0.25 mM and x0.8 at 20 mM",
 "=> A THIRD [ATP] discriminates them, and 20 mM is enough"
) 12 $INK | Out-Null

Add-Card $s ($M+$hw5+26) 138 $hw5 250 (Hex2Ole 0xF2 0xE2 0xE4) | Out-Null
Add-Text $s ($M+$hw5+48) 158 ($hw5-44) 22 "One gate is falsified on its own physiology" 15 $CHERRY $true 1 | Out-Null
Add-Bullets $s ($M+$hw5+48) 190 ($hw5-44) 186 @(
 "kmsrd -- the largest ATP effect in the H-series fit -- implies an ATP on-rate 100-1000x TOO SLOW to be ATP binding",
 "kah, ksr0 and ksr2srd come out textbook",
 "So the dominant fitted gate cannot be the mechanism it is named after"
) 12 $INK | Out-Null

Add-Text $s $M 406 ($W-2*$M) 24 "And what the 2 mM fit still fails" 15 $NAVY $true 1 | Out-Null
Add-Text $s $M 434 ($W-2*$M) 54 "rsK (model 0.627 vs data 0.345 as a ratio) and vall2_dy. But rsK is NOT an ATP-pathway question -- no rung of the ladder moves it below +8.8 sd, and it is 60% off at 8 mM too. It is the same structural defect as the 8 mM 'ugly'." 12.5 $INK $false 1 | Out-Null
Add-Note $s "The 8 mM core is fitted to a 3-preparation pooled target; the 2 mM leg is fitted to 03/27 alone. The two conditions are not on the same footing and the 2 mM numbers carry no error bar."

# =============================================================================
# 20  TASKS
# =============================================================================
$s = New-Slide $PAPER
Add-Title $s "Current challenges and tasks" $INK
$tasks = @(
 @('DONE','Pooled 2 mM target + condition-aware scoring','built today: own spread s2, own prior, and rsK80 finally scorable at 2 mM. Re-fit in flight'),
 @('1','Re-weight the 2 mM objective from s2 throughout','the values are pooled now; BuildWeightedObjective still derives every weight from the 8 mM spread'),
 @('2','Put the staircase in the objective','free out-of-sample test, known wrong in three specific ways, features and weights already derived'),
 @('3','Fix the two unfloored RMSE terms','18% of the objective, unprincipled, same species as the retired coolDownLS'),
 @('4','Narrow the F-V floor (per-term FLOOR_CHI2)','it drifted to 2 sd while paying exactly zero'),
 @('5','The re-stretch transient -- one defect behind five features','the only remaining structural error; needs MECHANISM, not parameters'),
 @('6','A third [ATP], ideally 20 mM','the only way to separate Hill from power gates -- currently unidentifiable in principle'),
 @('7','More preparations (n = 3 -> 6)','weights known only to a factor of ~18; protocol already designed (8-4-2-8 / 2-8-2)')
)
# Eight rows at 48 pt pitch from y = 118 ends at 502 -- inside the 540 pt slide.
# The first row is a completed item, so it is keyed on the LABEL, not on a cast:
# [int]"DONE" throws.
$y = 118
foreach ($t in $tasks) {
    $done = ($t[0] -eq 'DONE')
    if ($done)            { $fill = Hex2Ole 0xE9 0xEE 0xE9; $lc = $GREEN }
    elseif ($t[0] -le '2') { $fill = Hex2Ole 0xF2 0xE2 0xE4; $lc = $CHERRY }
    else                  { $fill = $OFFW;                  $lc = $CHERRY }
    Add-Card $s $M $y ($W-2*$M) 44 $fill | Out-Null
    $lw = if ($done) { 52 } else { 30 }
    $ls = if ($done) { 11 } else { 15 }
    Add-Text $s ($M+16) ($y+13) $lw 22 $t[0] $ls $lc $true 1 | Out-Null
    Add-Text $s ($M+76) ($y+4) 320 38 $t[1] 12.5 $INK $true 1 | Out-Null
    Add-Text $s ($M+408) ($y+4) ($W-2*$M-426) 38 $t[2] 11 $GREY $false 1 | Out-Null
    $y += 48
}
Add-Note $s "Also pending: delete the dead switches. UseSerialStiffness is listed in CLAUDE.md as a key behavioural switch and is read nowhere in the ODE."

# =============================================================================
# 21  CLOSING
# =============================================================================
$s = New-Slide $INK
Add-Text $s $M 130 ($W-2*$M) 40 "Is this the best we can do?" 24 0x00A0A0A0 $false 1 | Out-Null
Add-Text $s $M 186 ($W-2*$M) 110 "On this objective and this`ndataset -- close to it." 42 $PAPER $true 1 | Out-Null
Add-Text $s $M 330 ($W-2*$M) 90 "The remaining gain is not in the optimiser -- 39.6 hours of systematic ablation says so. It is in MORE DATA (preparations, a third [ATP], scoring the held-out protocols) and in ONE MECHANISM (the re-stretch transient)." 17 0x00D0D0D0 $false 1 | Out-Null
Add-Text $s $M 434 ($W-2*$M) 40 "Both are tractable. Neither is another fitting campaign." 20 $CHERRY $true 1 | Out-Null
Add-Note $s ""

# =============================================================================
$nSlides = $pres.Slides.Count
if (Test-Path $out) { Remove-Item $out -Force }
$pres.SaveAs($out)
$pres.Close()
$pp.Quit()
[System.Runtime.InteropServices.Marshal]::ReleaseComObject($pp) | Out-Null
"slides: $nSlides"
"written: $out"
