# The 12 h run: what it is for, and what would make me stop it

Started 2026-09-22 ~20:45. Six 2 h legs, `RESUME` between them, `TAG = refit6`.

## The bet

`ktr_shape` goes from 4 % to ~18 % of the objective, and the two dead force-gated SRX
branches (`kmsr`, `ksrd`, with `sigma1` / `sigma_srd2` behind them) become tunable for the
first time. refit5 proved the slope levers exist and are simply outbid; this is the
re-pricing plus the one mechanism the optimiser has never had.

## The constraint that outranks the objective

The deliverable is the **low-ATP fit**, at both concentrations, on the ratios and the
behaviour. 2 mM is **ktr x0.54** against **force x1.18**, and the earlier analysis
concluded the 2 mM fit is spoiled by high-ATP **re-stretch** error -- whatever is not
caught at 8 mM is worse at 2 mM. Two consequences:

* **`ktr` is the ATP-sensitive observable.** A model whose ktr has the wrong dependence on
  LENGTH has the wrong ktr mechanism, so its dependence on ATP is not to be trusted either.
  That is why the ladder slope is worth 18 % and not 4 %.
* **The re-stretch block is not currency.** `peak1_y, peak1_dSL, vall_y, vall_t, vall2_dy,
  peak2, restretchSlopeStart, rsK80`. A leg that lowers the total while this block rises is
  going the wrong way even though the scalar says otherwise.

## Check-in every 2 h -- the rules, decided in advance

`matlab -batch "TAG='refit6'; run('Analyses/KtrSegmentSlope/RunCheckin.m')"` reads
`params/refit6_state.mat` only, so it costs no cores.

| trigger | threshold | action |
|---|---|---|
| **selling re-stretch** | block > 1.20x its seed value | stop; re-floor the offending term and restart the leg from the last good snapshot |
| **ktr_shape not moving** | > 0.70x seed after leg 2 (4 h) | the price is still not enough OR the mechanism is unreachable -> go to the fallback below, do not keep paying for the same negative result |
| **F-V off zero** | any nonzero `FV_normpowerAvg` | the force scale is being traded; stop and look, this is the term that killed every structural variant |
| **hard stall** | 3+ consecutive no-improvement rounds with kicks failing | rotate the compulsory set rather than burn the remaining legs on the same draw |
| **held-out drift** | staircase or ktr-protocol rmse/noise worse at a 4 h checkpoint | overfitting; stop and report rather than bank it |

Otherwise: relaunch the next 2 h leg with `RESUME` and say nothing.

## Fallback if the ladder still will not move

Do not spend legs 5 and 6 repeating a failure. Switch to answering the question properly:
score `UsePassiveForSR = 1` **with a short refit** (it buys -0.63 of slope and its whole
cost is F-V, which a refit is exactly the thing that could absorb). If that also fails,
the honest result is that this parametrisation has no length-dependent RATE mechanism the
F-V curve will tolerate, and that is a finding worth writing down rather than a campaign
worth extending.

## Log

| leg | wall | best | re-stretch | ktr_shape | decision |
|---|---|---|---|---|---|
