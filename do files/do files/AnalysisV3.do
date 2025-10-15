*******************************************************
* 03_analysis_simple_v19.do — simple (Stata 19)
* Input: $PROC\aebcorrsv3.dta  (post-2014, +corn/soy, +Fourier)
*******************************************************

version 19.0
clear all
set more off

* 0) Load and set time
use "$PROC\aebcorrsv3.dta", clear
capture confirm variable mdate
if _rc {
    di as err "mdate missing"
    exit 459
}
format mdate %tm
tsset mdate, monthly

* Suffix for all outputs in this run
local SUF "_v3"

* Ensure AEB exists (canonical alias of AEB_aeb for convenience)
capture confirm variable AEB
if _rc {
    capture confirm variable AEB_aeb
    if !_rc {
        gen AEB = AEB_aeb
    }
}
capture confirm variable AEB
if _rc {
    di as err "AEB not found in the dataset."
    exit 111
}

* Identify numeric vars; exclude mdate and Fourier for analyses
ds, has(type numeric)
local allnum `r(varlist)'
local allnum : list allnum - mdate
ds F_s* F_c*, has(type numeric)
local four `r(varlist)'

* Keep a "no-Fourier" numeric set for transforms/analyses
local nums_nofour : list allnum - four

*******************************************************
* 1) Z-score all numeric vars (exclude mdate, Fourier, AEB alias)
*******************************************************
* Avoid duplicate z_ for AEB alias; we'll alias z_AEB later
local zbase `nums_nofour'
local zbase : list zbase - AEB
foreach v of local zbase {
    capture drop z_`v'
    quietly egen z_`v' = std(`v')
}
* Create z_AEB as alias of z_AEB_aeb if needed
capture confirm variable z_AEB
if _rc {
    capture confirm variable z_AEB_aeb
    if !_rc gen double z_AEB = z_AEB_aeb
}

**************************************************************
* 2) Pairwise correlations: z_AEB vs everyone else (no Fourier)
**************************************************************
* ---- Guardrails: build clean target list (no mdate/Fourier/AEB), require z_ partner
local others `nums_nofour'
local others : list others - mdate
local others : list others - AEB
local others : list others - AEB_aeb

* Keep only vars that actually have z_ counterparts
local safe_others
foreach v of local others {
    capture confirm variable z_`v'
    if !_rc {
        local safe_others `safe_others' `v'
    }
}

tempfile out
tempname ph
postfile `ph' str64 var double rho int N using "`out'", replace
foreach v of local safe_others {
    quietly corr z_AEB z_`v'
    post `ph' ("`v'") (r(rho)) (r(N))
}
postclose `ph'

use "`out'", clear
drop if missing(rho)
gen double abs_rho = abs(rho)
gsort -abs_rho
order var rho N
format rho %6.3f
format N   %9.0f
export delimited using "$TAB\T3_corr_AEB_pairwise`SUF'.csv", replace
drop abs_rho

*****************************************************************
* 3) (Optional) Rolling 24-month correlations (exclude Fourier)
*****************************************************************
cap which rangestat
if _rc ssc install rangestat, replace

preserve
    use "$PROC\aebcorrsv3.dta", clear
    tsset mdate, monthly

* ---- Guardrails: rebuild clean lists
ds, has(type numeric)
local allnum `r(varlist)'
local allnum : list allnum - mdate

ds F_s* F_c*, has(type numeric)
local four `r(varlist)'

local nums_nofour : list allnum - four

* Rebuild z_ for these only (skip AEB alias to avoid duplicates)
local zbase `nums_nofour'
local zbase : list zbase - AEB
foreach v of local zbase {
    capture drop z_`v'
    quietly egen z_`v' = std(`v')
}
capture confirm variable z_AEB
if _rc {
    capture confirm variable z_AEB_aeb
    if !_rc gen double z_AEB = z_AEB_aeb
}

* Rolling targets: no mdate/AEB; must have z_ partner
local others `nums_nofour'
local others : list others - mdate
local others : list others - AEB
local others : list others - AEB_aeb

local safe_others
foreach v of local others {
    capture confirm variable z_`v'
    if !_rc {
        local safe_others `safe_others' `v'
    }
}

* ---- Rolling correlations
foreach v of local safe_others {
    tempvar prod mx my mxy sx sy
    gen double `prod' = z_AEB * z_`v'
    rangestat (mean) `mx'=z_AEB `my'=z_`v' `mxy'=`prod' ///
              (sd)   `sx'=z_AEB `sy'=z_`v', interval(mdate -23 0)

    local base = strtoname("`v'")
    local maxbase = 32 - `=strlen("r_AEB_")'
    if strlen("`base'") > `maxbase' {
        local base = substr("`base'", 1, `maxbase')
    }
    local rname = "r_AEB_`base'"
    capture confirm variable `rname'
    local k = 1
    while !_rc {
        local rname = "r_AEB_`base'_`k'"
        local ++k
        capture confirm variable `rname'
    }

    capture drop `rname'
    gen double `rname' = (`mxy' - `mx'*`my') / (`sx'*`sy')
    replace    `rname' = . if (`sx'==0 | `sy'==0)
    drop `prod' `mx' `my' `mxy' `sx' `sy'
}

    order mdate, first
    compress
    save "$PROC\monthly_rollingcorrs_simple`SUF'.dta", replace
restore

**************************************************************
* 4) (Optional) Quick chart for one rolling series
**************************************************************
use "$PROC\monthly_rollingcorrs_simple`SUF'.dta", clear
local target r_AEB_News_Based_Policy_Uncert_Index
capture confirm variable `target'
if _rc {
    capture unab _cands : r_AEB_*
    if !_rc {
        local target : word 1 of `_cands'
    }
}
capture confirm variable `target'
if !_rc {
    twoway tsline `target', ///
        title("Rolling 24m corr: AEB vs selected") ytitle("corr") xtitle("")
    graph export "$FIG\F3_roll_AEB_selected`SUF'.png", replace width(1600)
}
* ---- Cleanup: drop any rolling-corr vars from the active dataset
capture unab rvars : r_*
if !_rc {
    drop `rvars'
}


**************************************************************
* 5) PCA + FA — auto-build varlists (include corn/sb vols)
**************************************************************
* Build RAW = all numeric except mdate, Fourier, z_* and the AEB alias
ds, has(type numeric)
local Xraw `r(varlist)'
local Xraw : list Xraw - mdate
local Xraw : list Xraw - four
local Xraw : list Xraw - AEB
capture unab zvars : z_*
if !_rc {
    local Xraw : list Xraw - zvars
}

* Build Z = all z_* constructed above
capture unab Z : z_*
if _rc local Z ""

local p_raw : word count `Xraw'
local p_z   : word count `Z'
if (`p_raw' < 2 & `p_z' < 2) {
    di as err "Not enough variables for PCA/FA."
    exit 111
}
local ncomp_raw = cond(`p_raw' >= 3, 3, `p_raw')
local ncomp_z   = cond(`p_z'   >= 3, 3, `p_z')

*******************************************************
* A) RAW (nominal) variables
*******************************************************
if `p_raw' >= 2 {
    quietly describe `Xraw'
    quietly summarize `Xraw'
    quietly corr `Xraw'

    * PCA (RAW)
    pca `Xraw'
    screeplot, name(G_scree_raw, replace)
    screeplot, yline(1) name(G_scree_raw_y1, replace)
    pca `Xraw', mineigen(1)
    pca `Xraw', comp(`ncomp_raw')
    pca `Xraw', comp(`ncomp_raw') blanks(.3)
    rotate, varimax

    estat loadings
    matrix L_pca_raw = e(L)
    preserve
        clear
        svmat double L_pca_raw, names(col)
        gen variable = ""
        local rn : rownames L_pca_raw
        local i = 1
        foreach r of local rn {
            replace variable = "`r'" in `i'
            local ++i
        }
        order variable
        export delimited using "$TAB\T_loadings_pca_raw_varimax`SUF'.csv", replace
    restore

    * PCA scores -> *_raw
    local pcs_raw
    forvalues i = 1/`ncomp_raw' {
        local pcs_raw `pcs_raw' pc`i'
    }
    capture drop `pcs_raw'
    predict `pcs_raw', score
    foreach v of local pcs_raw {
        rename `v' `v'_raw
    }

    * KMO
    estat kmo

    * FACTOR (RAW)
    factor `Xraw'
    screeplot, name(G_scree_raw_fa, replace)
    screeplot, yline(1) name(G_scree_raw_fa_y1, replace)
    factor `Xraw', mineigen(1)
    factor `Xraw', factor(`ncomp_raw')
    factor `Xraw', factor(`ncomp_raw') blanks(0.3)
    rotate, varimax

    estat common
    matrix L_fa_raw = e(L)
    preserve
        clear
        svmat double L_fa_raw, names(col)
        gen variable = ""
        local rn : rownames L_fa_raw
        local i = 1
        foreach r of local rn {
            replace variable = "`r'" in `i'
            local ++i
        }
        order variable
        export delimited using "$TAB\T_loadings_fa_raw_varimax`SUF'.csv", replace
    restore

    * FA scores -> *_raw
    local fs_raw
    forvalues i = 1/`ncomp_raw' {
        local fs_raw `fs_raw' f`i'
    }
    capture drop `fs_raw'
    predict `fs_raw'
    foreach v of local fs_raw {
        rename `v' `v'_raw
    }

    * Reliability + Bartlett (RAW)
    alpha `Xraw'
    cap which factortest
    if _rc ssc install factortest, replace
    factortest `Xraw'
}

*******************************************************
* B) STANDARDIZED (z_) variables
*******************************************************
if `p_z' >= 2 {
    quietly describe `Z'
    quietly summarize `Z'
    quietly corr `Z'

    * PCA (Z)
    pca `Z'
    screeplot, name(G_scree_z, replace)
    screeplot, yline(1) name(G_scree_z_y1, replace)
    pca `Z', mineigen(1)
    pca `Z', comp(`ncomp_z')
    pca `Z', comp(`ncomp_z') blanks(.3)
    rotate, varimax

    estat loadings
    matrix L_pca_z = e(L)
    preserve
        clear
        svmat double L_pca_z, names(col)
        gen variable = ""
        local rn : rownames L_pca_z
        local i = 1
        foreach r of local rn {
            replace variable = "`r'" in `i'
            local ++i
        }
        order variable
        export delimited using "$TAB\T_loadings_pca_z_varimax`SUF'.csv", replace
    restore

    * PCA scores -> *_z
    local pcs_z
    forvalues i = 1/`ncomp_z' {
        local pcs_z `pcs_z' pc`i'
    }
    capture drop `pcs_z'
    predict `pcs_z', score
    foreach v of local pcs_z {
        rename `v' `v'_z
    }

    * KMO
    estat kmo

    * FACTOR (Z)
    factor `Z'
    screeplot, name(G_scree_z_fa, replace)
    screeplot, yline(1) name(G_scree_z_fa_y1, replace)
    factor `Z', mineigen(1)
    factor `Z', factor(`ncomp_z')
    factor `Z', factor(`ncomp_z') blanks(0.3)
    rotate, varimax

    estat common
    matrix L_fa_z = e(L)
    preserve
        clear
        svmat double L_fa_z, names(col)
        gen variable = ""
        local rn : rownames L_fa_z
        local i = 1
        foreach r of local rn {
            replace variable = "`r'" in `i'
            local ++i
        }
        order variable
        export delimited using "$TAB\T_loadings_fa_z_varimax`SUF'.csv", replace
    restore

    * FA scores -> *_z
    local fs_z
    forvalues i = 1/`ncomp_z' {
        local fs_z `fs_z' f`i'
    }
    capture drop `fs_z'
    predict `fs_z'
    foreach v of local fs_z {
        rename `v' `v'_z
    }

    * Reliability + Bartlett (Z)
    alpha `Z'
    cap which factortest
    if _rc ssc install factortest, replace
    factortest `Z'
}

**************************************************************
* C) Save tidy score files — robust to missing scores
**************************************************************
* RAW scores
local have_raw ""
capture unab have_raw : pc*_raw f*_raw
if !_rc {
    preserve
        keep mdate `have_raw'
        compress
        save "$PROC\fa_pca_scores_raw`SUF'.dta", replace
    restore
}

* Z scores
local have_z ""
capture unab have_z : pc*_z f*_z
if !_rc {
    preserve
        keep mdate `have_z'
        compress
        save "$PROC\fa_pca_scores_z`SUF'.dta", replace
    restore
}

**************************************************************
* 6) Single–Factor "AEB-like" index (RAW and Z)
**************************************************************
* Use the discovered Xraw/Z lists above
local p_raw : word count `Xraw'
local p_z   : word count `Z'

* RAW single factor
if `p_raw' >= 2 {
    factor `Xraw', factor(1)
    matrix L_fa1_raw = e(L)
    preserve
        clear
        svmat double L_fa1_raw, names(col)
        gen variable = ""
        local rn : rownames L_fa1_raw
        local i = 1
        foreach r of local rn {
            replace variable = "`r'" in `i'
            local ++i
        }
        order variable
        export delimited using "$TAB\T_loadings_fa1_raw`SUF'.csv", replace
    restore

    capture drop F1_raw
    predict F1_raw
    capture confirm variable AEB_aeb
    if !_rc {
        quietly corr AEB_aeb F1_raw
        scalar s_raw = sign(r(rho))
        if s_raw < 0 {
            replace F1_raw = -F1_raw
        }
    }
    rename F1_raw AEB_like_raw
    preserve
        keep mdate AEB_like_raw
        compress
        save "$PROC\fa1_scores_raw`SUF'.dta", replace
    restore
}

* Z single factor
if `p_z' >= 2 {
    factor `Z', factor(1)
    matrix L_fa1_z = e(L)
    preserve
        clear
        svmat double L_fa1_z, names(col)
        gen variable = ""
        local rn : rownames L_fa1_z
        local i = 1
        foreach r of local rn {
            replace variable = "`r'" in `i'
            local ++i
        }
        order variable
        export delimited using "$TAB\T_loadings_fa1_z`SUF'.csv", replace
    restore

    capture drop F1_z
    predict F1_z
    capture confirm variable z_AEB_aeb
    if !_rc {
        quietly corr z_AEB_aeb F1_z
        scalar s_z = sign(r(rho))
        if s_z < 0 {
            replace F1_z = -F1_z
        }
    }
    rename F1_z AEB_like_z
    preserve
        keep mdate AEB_like_z
        compress
        save "$PROC\fa1_scores_z`SUF'.dta", replace
    restore
}

* Optional CSV exports of single-factor series
foreach which in raw z {
    capture confirm file "$PROC\fa1_scores_`which'`SUF'.dta"
    if !_rc {
        preserve
            use "$PROC\fa1_scores_`which'`SUF'.dta", clear
            export delimited using "$TAB\T_fa1_scores_`which'`SUF'.csv", replace
        restore
    }
}

*******************************************************
* 6) OLS: ΔAEB ~ all other first-difference vars (d_*)
*******************************************************
use "$PROC\aebcorrsv3_diff.dta", clear
format mdate %tm
tsset mdate, monthly

* Dependent variable (prefer d_AEB_aeb; fallback d_AEB)
local y ""
capture confirm variable d_AEB_aeb
if !_rc local y d_AEB_aeb
else {
    capture confirm variable d_AEB
    if !_rc local y d_AEB
}
if "`y'" == "" {
    di as err "Dependent variable d_AEB_aeb (or d_AEB) not found."
    exit 111
}

* Candidate predictors = all d_* vars
capture unab cand : d_*
if _rc {
    di as err "No d_* variables found."
    exit 111
}

* Exclusions: DV itself, d_mdate, any AEB variants, and differenced Fourier
local drop " `y' d_mdate d_AEB_aeb d_AEB "
local X ""
foreach v of local cand {
    * skip if in drop list or looks like differenced Fourier
    if strpos("`drop'"," `v' ")==0 & substr("`v'",1,4)!="d_F_" {
        local X `X' `v'
    }
}

* Safety: ensure predictors remain
local p : word count `X'
if `p'==0 {
    di as err "No predictors left after exclusions."
    exit 111
}

* Fit OLS with robust SEs
regress `y' `X', vce(robust)

* Multicollinearity check
estat vif

* ---- Top 5 most significant predictors (post-regress)
tempname b V
matrix `b' = e(b)
matrix `V' = e(V)
local cn : colnames `b'

tempfile _top5
postfile __pf str80 var double coef se t p using "`_top5'", replace
foreach v of local cn {
    if "`v'" != "_cons" {
        scalar __b  = `b'[1, colnumb(`b', "`v'")]
        scalar __se = sqrt(`V'[colnumb(`b', "`v'"), colnumb(`b', "`v'")])
        if __se < . & __se > 0 {
            scalar __t = __b / __se
            scalar __p = 2*ttail(e(df_r), abs(__t))
            post __pf ("`v'") (__b) (__se) (__t) (__p)
        }
    }
}
postclose __pf

use "`_top5'", clear
gsort +p
keep in 1/5
order var coef se t p
format coef se %9.3g
format t %8.2f
format p %6.4f
list, noobs

* Optional: export a CSV (uses your SUF if defined)
capture confirm local SUF
if _rc local SUF ""
export delimited using "$TAB\T_reg_top5`SUF'.csv", replace


*******************************************************
* 7) VAR & Granger causality — regular aebcorrsv3
*******************************************************
* You can run in levels or first differences (stationary recommended).
*  - USE_DIFF = 0 -> levels on aebcorrsv3.dta
*  - USE_DIFF = 1 -> uses differenced vars from aebcorrsv3_diff.dta
local USE_DIFF 0

* Optional file suffix from earlier; set empty if not defined
capture confirm local SUF
if _rc local SUF "_v3"

* ---------------------------
* A) Load data (levels/diff)
* ---------------------------
if `USE_DIFF'==0 {
    use "$PROC\aebcorrsv3.dta", clear
}
else {
    use "$PROC\aebcorrsv3_diff.dta", clear
}

format mdate %tm
tsset mdate, monthly

* ---------------------------
* B) Build variable set
* ---------------------------
* Choose AEB variable name
local AEBvar ""
capture confirm variable AEB
if !_rc local AEBvar AEB
else {
    capture confirm variable AEB_aeb
    if !_rc local AEBvar AEB_aeb
}
if "`AEBvar'"=="" {
    di as err "AEB (or AEB_aeb) not found."
    exit 111
}

* Core candidates (will include only those that exist)
local base_candidates ///
    News_Based_Policy_Uncert_Index ///
    CCI_monthly_cci ///
    UMCSENT_umc ///
    vix_vix ///
    corn_vol_month ///
    sb_vol_month

* If using differences, prefix d_ to names (AEB handled below)
local VARLIST ""
if `USE_DIFF'==0 {
    * Levels
    foreach c of local base_candidates {
        capture confirm variable `c'
        if !_rc local VARLIST `VARLIST' `c'
    }
    local DV "`AEBvar'"
}
else {
    * First differences
    local DV d_`AEBvar'
    capture confirm variable `DV'
    if _rc {
        di as err "Dependent variable `DV' not found in aebcorrsv3_diff.dta."
        exit 111
    }
    foreach c of local base_candidates {
        capture confirm variable d_`c'
        if !_rc local VARLIST `VARLIST' d_`c'
    }
}

* Final endogenous set = DV + available others
local ENDOG "`DV' `VARLIST'"

* Guard: drop Fourier if they slipped in
capture unab four : F_s* F_c* d_F_s* d_F_c*
if !_rc {
    local ENDOG : list ENDOG - four
}

* Ensure we have at least 2 endogenous variables
local k : word count `ENDOG'
if `k' < 2 {
    di as err "Not enough variables for VAR (need >=2). Endog: `ENDOG'"
    exit 111
}

* -----------------------------------
* C) Select lags via varsoc (max 6)
* -----------------------------------
* We'll parse SBIC to choose lags; fallback to 2 if parsing fails.
local LMAX = 6
capture noisily varsoc `ENDOG', maxlags(`LMAX')
local PSEL = 2
capture mata: st_local("PSEL", strofreal(panelsetup(st_matrix("r(stats)"))))  // harmless if r(stats) missing

* Programmatic pick (SBIC): scan info from r(stats); if not available, keep 2
capture matrix S = r(stats)
if !_rc {
    * r(stats): rows = lags, cols = ll aic hqic sbic
    local n = rowsof(S)
    local best = .
    local argbest = .
    forvalues p = 1/`n' {
        scalar sb = S[`p',4]
        if missing(`best') | sb < `best' {
            scalar best = sb
            local argbest = `p'
        }
    }
    if "`argbest'"!="" local PSEL = `argbest'
}

* -----------------------------------
* D) Estimate VAR and diagnostics
* -----------------------------------
quietly var `ENDOG', lags(1/`PSEL') dfk small

* Stability check (roots inside unit circle)
varstable, graph
graph export "$FIG\G_VAR_stability`SUF'.png", replace width(1600)

* Residual serial correlation (LM test)
varlmar, mlag(`PSEL')

* -----------------------------------
* E) Granger causality
* -----------------------------------
* Global table across all equations
vargranger

* Directed to AEB (AEB equation only): test each driver's lags jointly = 0
* Build list of driver vars (exclude the DV itself)
local DRIVERS : list ENDOG - DV

* Resolve equation name = the DV variable's equation
* In VAR, equation names equal variable names.
foreach x of local DRIVERS {
    * Test: do lags of x jointly Granger-cause DV in DV-equation
    test ([`DV']L(1/`PSEL').`x' = 0)
}

* -----------------------------------
* F) Optional: IRFs (levels only advisable if stationary/cointegrated)
* -----------------------------------
* Create IRFs for shocks to EPU/VIX/corn/soy into AEB (when present)
tempname has_epu has_vix has_corn has_sb
scalar `has_epu'  = 0
scalar `has_vix'  = 0
scalar `has_corn' = 0
scalar `has_sb'   = 0

if `USE_DIFF'==0 {
    capture confirm variable News_Based_Policy_Uncert_Index
    if !_rc scalar `has_epu' = 1
    capture confirm variable vix_vix
    if !_rc scalar `has_vix' = 1
    capture confirm variable corn_vol_month
    if !_rc scalar `has_corn' = 1
    capture confirm variable sb_vol_month
    if !_rc scalar `has_sb' = 1

    irf set "$PROC\IRF_aeb`SUF'", replace
    irf create var_lvl`SUF', step(24) replace

    local

