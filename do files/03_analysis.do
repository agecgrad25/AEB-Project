*03_analysis

*******************************************************
* 03_analysis.do — Step 3: Benchmark correlations for AEB
* Project paths/globals defined in master.do:
*   $RAW  $PROC  $FIG  $TAB  $DOC  $CODE
* Requires: 00_setup_install.do already ran (packages)
* Stata: v19
*******************************************************

version 19.0
clear all
set more off

* ------------------ Safety: ensure output dirs exist ------------------
cap mkdir "$TAB"
cap mkdir "$FIG"
cap mkdir "$DOC"

* ------------------ Load source panel ------------------
* Use the compiled monthly merge from 01_clean_raw.do
capture confirm file "$PROC\First compilation.dta"
if _rc {
    di as err "Cannot find $PROC\First compilation.dta (from 01_clean_raw.do)."
    exit 601
}
use "$PROC\First compilation.dta", clear

* Expect a monthly date var mdate
capture confirm variable mdate
if _rc {
    di as err "mdate not found; please ensure earlier steps built a monthly index."
    exit 198
}
format mdate %tm
tsset mdate, monthly

* ------------------ Harmonize common variable names ------------------
* We standardize aliases into canonical names expected by downstream code.
* (Only created if an original exists; otherwise skipped)

* MCSI / UMCSENT
capture confirm variable MCSI
if _rc {
    capture confirm variable UMCSENT
    if !_rc gen MCSI = UMCSENT
}

* Trade policy uncertainty: TradeEPU / TPU
capture confirm variable TradeEPU
if _rc {
    capture confirm variable TPU
    if !_rc gen TradeEPU = TPU
}

* Headline EPU variants: USEPU / EPU
capture confirm variable EPU
if _rc {
    capture confirm variable USEPU
    if !_rc gen EPU = USEPU
}

* Consumer Confidence: CCI / CCI_Value
capture confirm variable CCI
if _rc {
    capture confirm variable CCI_Value
    if !_rc gen CCI = CCI_Value
}

* Partisanship proxies: partis_r / partis_d — try to synthesize if you only have party splits
* (Optional) If you have SI_Rep / SI_Dem style vars, you can uncomment to build rollups:
* capture confirm variable partis_r
* if _rc {
*     capture confirm variable SI_Rep
*     if !_rc gen partis_r = SI_Rep
* }
* capture confirm variable partis_d
* if _rc {
*     capture confirm variable SI_Dem
*     if !_rc gen partis_d = SI_Dem
* }

* ------------------ Package checks (defensive) ------------------
cap which estpost
if _rc ssc install estout, replace
cap which esttab
if _rc ssc install estout, replace
cap which rangestat
if _rc ssc install rangestat, replace

* ------------------ 3.1 Stationarity checks & transforms ------------------
local series AEB ICC IFE MCSI CCI partis_r partis_d EPU TradeEPU vix rv_corn rv_soy machinery_total

file close _all
capture erase "$DOC\stationarity_decisions.md"
file open SM using "$DOC\stationarity_decisions.md", write replace text
file write SM "# Stationarity decisions (automated log)" _n
file write SM "Generated: " c(current_date) " " c(current_time) _n _n
file write SM "| Variable | DF (lags=12) p | PP p | Decision | Transform |" _n
file write SM "|---|---:|---:|---|---|" _n

tempname dfpp
tempfile testdump
postfile `dfpp' str32 varname double p_df p_pp using "`testdump'", replace

foreach v of local series {
    capture confirm variable `v'
    if _rc continue
    quietly dfuller `v', lags(12)
    scalar p_df = r(pvalue)
    quietly pperron `v'
    scalar p_pp = r(pvalue)
    post `dfpp' ("`v'") (p_df) (p_pp)
}
postclose `dfpp'

preserve
use "`testdump'", clear
gen byte I1 = (p_df>0.05 & p_pp>0.05)
gen str12 decision  = cond(I1,"I(1)","I(0)")
gen str12 transform = cond(I1,"difference","level")

quietly {
    forvalues i=1/`=_N' {
        local v   = varname[`i']
        local pdf = string(p_df[`i'],"%9.4f")
        local ppp = string(p_pp[`i'],"%9.4f")
        local dcs = decision[`i']
        local trn = transform[`i']
        file write SM "| `v' | `pdf' | `ppp' | `dcs' | `trn' |" _n
    }
}
levelsof varname if I1, local(I1vars)
file close SM
restore

* Differences for I(1) variables (only if they exist)
foreach v of local I1vars {
    capture drop D_`v'
    gen D_`v' = D.`v'
}

* Z-score standardization (levels and differences)
foreach v of local series {
    capture confirm variable `v'
    if !_rc {
        capture drop z_`v'
        quietly egen z_`v' = std(`v')
    }
}
foreach v of local I1vars {
    capture confirm variable D_`v'
    if !_rc {
        capture drop z_D_`v'
        quietly egen z_D_`v' = std(D_`v')
    }
}

* Save transformed panel for reuse
save "$PROC\monthly_panel_transformed.dta", replace

* ------------------ 3.2 Static correlations (levels & changes) ------------------
use "$PROC\monthly_panel_transformed.dta", clear

* Candidate sets (existence-filtered)
local Lcands z_AEB z_MCSI z_CCI z_partis_r z_partis_d z_EPU z_TradeEPU z_vix z_rv_corn z_rv_soy
local Lexist
foreach v of local Lcands {
    capture confirm variable `v'
    if !_rc local Lexist `Lexist' `v'
}

local Dcands z_D_AEB z_D_MCSI z_D_CCI z_D_partis_r z_D_partis_d z_D_EPU z_D_TradeEPU
local Dexist
foreach v of local Dcands {
    capture confirm variable `v'
    if !_rc local Dexist `Dexist' `v'
}

* Levels
if "`Lexist'" != "" {
    estpost corr `Lexist', listwise
    esttab using "$TAB\T3_static_correlations.csv", ///
        cells("rho(fmt(3))") unstack replace nonote noobs
}
else di as err "No level variables found for static correlations."

* Changes
if "`Dexist'" != "" {
    estpost corr `Dexist', listwise
    esttab using "$TAB\T3_static_correlations_changes.csv", ///
        cells("rho(fmt(3))") unstack replace nonote noobs
}
else di as txt "No differenced variables present for Δ correlations (skip)."

* ------------------ 3.3 Rolling 24-month correlations ------------------
use "$PROC\monthly_panel_transformed.dta", clear
tsset mdate, monthly

* We'll attempt these pairs if present
local pairs "EPU TradeEPU partis_r rv_corn"
foreach p of local pairs {
    capture confirm variable z_`p'
    if _rc {
        di as txt "Skip rolling corr: z_`p' not found."
        continue
    }
    capture drop r_AEB_`p'
    rangestat (corr) r_AEB_`p' = z_AEB z_`p', interval(mdate -23 0)
}

* Example single plot: AEB vs EPU
capture confirm variable r_AEB_EPU
if !_rc {
    twoway (tsline r_AEB_EPU), ///
        title("Rolling 24m Corr: AEB vs EPU") ///
        ytitle("corr") xtitle("") ///
        name(G_roll_AEB_EPU, replace)
    graph export "$FIG\F3_rollcorr_AEB_EPU.png", replace width(1600)
}

* Multi overlay (TradeEPU/TPU, Partisan-R, Corn RV)
local overlay ""
foreach p in TradeEPU partis_r rv_corn {
    capture confirm variable r_AEB_`p'
    if !_rc local overlay `overlay' r_AEB_`p'
}
if "`overlay'" != "" {
    twoway (tsline `overlay'), ///
        title("Rolling 24m Corr with AEB") ///
        legend(order(1 "TPU" 2 "Partisan-R" 3 "Corn RV")) ///
        ytitle("corr") xtitle("") ///
        name(G_roll_panel, replace)
    graph export "$FIG\F4_rollcorr_panel.png", replace width(1600)
}
else di as txt "No rolling-correlation overlay could be drawn (missing series)."

* ------------------ 3.4 Optional: Cointegration & VECM ------------------
* Only if both exist and were likely I(1)
local cpair "AEB MCSI"
local ok = 1
foreach v of local cpair {
    capture confirm variable `v'
    if _rc {
        di as txt "Skip cointegration: `v' not found."
        local ok = 0
    }
}
if `ok' {
    tsset mdate, monthly
    quietly vecrank AEB MCSI, trend(constant) lags(1/6)
    estimates store RANKSEL

    quietly vec AEB MCSI, rank(1) lags(1/2)
    estimates store VECM1

    esttab RANKSEL VECM1 using "$TAB\T3_vecm_summary.txt", ///
        replace wide se label nonumber nogaps nodepvars
}

* ------------------ Wrap up ------------------
di as result "Step 3 analysis complete (Stata 19)."
di as txt     "Artifacts:"
di as txt     "  - $DOC\stationarity_decisions.md"
di as txt     "  - $PROC\monthly_panel_transformed.dta"
di as txt     "  - $TAB\T3_static_correlations.csv"
di as txt     "  - $TAB\T3_static_correlations_changes.csv"
di as txt     "  - $TAB\T3_vecm_summary.txt (if cointegration ran)"
di as txt     "  - $FIG\F3_rollcorr_AEB_EPU.png"
di as txt     "  - $FIG\F4_rollcorr_panel.png"
