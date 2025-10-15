*******************************************************
* 03_analysis_simple_v19.do — keep it simple (Stata 19)
* Input: $PROC\firstcorrsv2_post2014.dta
*******************************************************

version 19.0
clear all
set more off

* 0) Load and set time
use "$PROC\firstcorrsv2_post2014.dta", clear
capture confirm variable mdate
if _rc {
    di as err "mdate missing"; exit 459
}
format mdate %tm
tsset mdate, monthly

* Ensure AEB is present
capture confirm variable AEB
if _rc {
    capture confirm variable AEB_aeb
    if !_rc gen AEB = AEB_aeb
}
capture confirm variable AEB
if _rc {
    di as err "AEB not found in the dataset."; exit 111
}

* 1) Z-score all numeric vars (except mdate)
ds, has(type numeric)
local nums `r(varlist)'
local nums : list nums - mdate

foreach v of local nums {
    capture drop z_`v'
    quietly egen z_`v' = std(`v')
}

* 2) Correlations: AEB vs everyone else (pairwise)
local others : list nums - AEB

tempfile out
tempname ph
postfile `ph' str64 var double rho int N using "`out'", replace
foreach v of local others {
    quietly corr z_AEB z_`v'   // pairwise by default
    post `ph' ("`v'") (r(rho)) (r(N))
}
postclose `ph'

* --- REOPEN results and sort by |rho| descending ---
use "`out'", clear
drop if missing(rho)
gen double abs_rho = abs(rho)
gsort -abs_rho

order var rho N
format rho %6.3f
format N   %9.0f
export delimited using "$TAB\T3_corr_AEB_pairwise.csv", replace
drop abs_rho
di as result "✅ Saved: $TAB\T3_corr_AEB_pairwise.csv"


* 3) (Optional) Rolling 24-month correlation of AEB vs everyone (safe names)
cap which rangestat
if _rc ssc install rangestat, replace

preserve
    use "$PROC\firstcorrsv2_post2014.dta", clear
    tsset mdate, monthly

    * Recreate z_ vars in this frame
    ds, has(type numeric)
    local nums `r(varlist)'
    local nums : list nums - mdate
    foreach v of local nums {
        capture drop z_`v'
        quietly egen z_`v' = std(`v')
    }
    local others : list nums - AEB

    foreach v of local others {
        * require z_ partner
        capture confirm variable z_`v'
        if _rc continue

        * windowed moments on z_ variables (consistent!)
        tempvar prod mx my mxy sx sy
        gen double `prod' = z_AEB * z_`v'
        rangestat (mean) `mx'=z_AEB `my'=z_`v' `mxy'=`prod' ///
                  (sd)   `sx'=z_AEB `sy'=z_`v', interval(mdate -23 0)

        * safe result name: r_AEB_<base> (<=32 chars, no spaces)
        local base = strtoname("`v'")
        local maxbase = 32 - `=strlen("r_AEB_")'
        if strlen("`base'") > `maxbase' local base = substr("`base'", 1, `maxbase')
        local rname = "r_AEB_`base'"

        * ensure uniqueness if truncation collides
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
    save "$PROC\monthly_rollingcorrs_simple.dta", replace
restore
di as result "✅ Saved: $PROC\monthly_rollingcorrs_simple.dta"


* 4) (Optional) Quick chart if EPU exists
use "$PROC\monthly_rollingcorrs_simple.dta", clear
capture confirm variable r_AEB_EPU
if !_rc {
    twoway tsline r_AEB_EPU, ///
        title("Rolling 24m corr: AEB vs EPU") ytitle("corr") xtitle("")
    graph export "$FIG\F3_roll_AEB_EPU.png", replace width(1600)
    di as result "✅ Saved: $FIG\F3_roll_AEB_EPU.png"
}

*---------------------------------------------------------------
* 5) PCA + FA (hard-coded varlists) on RAW and STANDARDIZED vars
*    Outputs:
*      $PROC\fa_pca_scores_raw.dta   (mdate + pc*_raw + f*_raw)
*      $PROC\fa_pca_scores_z.dta     (mdate + pc*_z   + f*_z)
*    NOTE: Assumes the dataset is already loaded and tsset mdate
*---------------------------------------------------------------

* 5a) Varlist setup for PCA/FA (Xraw vs Z)
* RAW (nominal) variables (exclude mdate)
local Xraw ///
    News_Based_Policy_Uncert_Index ///
    AEB_aeb ///
    SI_Dem_part ///
    SI_Ind_part ///
    SI_Rep_part ///
    CCI_monthly_cci ///
    UMCSENT_umc ///
    EOM_Close_de ///
    GEPU_current_gepu ///
    TPU_tpu ///
    SBOI_nfibo ///
    UncertaintyIndex_nfibu ///
    vix_vix

* STANDARDIZED (z_) variables (exclude z_mdate)
local Z ///
    z_News_Based_Policy_Uncert_Index ///
    z_AEB_aeb ///
    z_SI_Dem_part ///
    z_SI_Ind_part ///
    z_SI_Rep_part ///
    z_CCI_monthly_cci ///
    z_UMCSENT_umc ///
    z_EOM_Close_de ///
    z_GEPU_current_gepu ///
    z_TPU_tpu ///
    z_SBOI_nfibo ///
    z_UncertaintyIndex_nfibu ///
    z_vix_vix

* Guardrails
local p_raw : word count `Xraw'
local p_z   : word count `Z'
if (`p_raw' < 2 & `p_z' < 2) {
    di as err "Not enough variables for PCA/FA."
    exit 111
}

* Choose component/factor count (cap at 3 or #vars)
local ncomp_raw = cond(`p_raw' >= 3, 3, `p_raw')
local ncomp_z   = cond(`p_z'   >= 3, 3, `p_z')

*******************************************************
* A) RAW (nominal) variables
*******************************************************
if `p_raw' >= 2 {
    di as txt "---- RAW (nominal) variables ----"
    describe `Xraw'
    summarize `Xraw'
    corr `Xraw'

    * PCA (RAW)
    pca `Xraw'
    screeplot, name(G_scree_raw, replace)
    screeplot, yline(1) name(G_scree_raw_y1, replace)
    pca `Xraw', mineigen(1)
    pca `Xraw', comp(`ncomp_raw')
    pca `Xraw', comp(`ncomp_raw') blanks(.3)

    rotate, varimax
    rotate, varimax blanks(.3)
    rotate, clear
    rotate, promax
    rotate, promax blanks(.3)
    rotate, clear

    estat loadings
	
	* Save rotated PCA loadings (RAW) to CSV
matrix L_pca_raw = e(L)
preserve
    clear
    svmat double L_pca_raw, names(col)   // makes vars = component names
    gen variable = ""                    // add rownames (original vars)
    local rn : rownames L_pca_raw
    local i = 1
    foreach r of local rn {
        replace variable = "`r'" in `i'
        local ++i
    }
    order variable
    export delimited using "$TAB\T_loadings_pca_raw_varimax.csv", replace
restore


    * PCA scores (pc1..pc`ncomp_raw' -> *_raw)
    local pcs_raw
    forvalues i=1/`ncomp_raw' {
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
    factor `Xraw', factor(`ncomp_raw') pcf

    rotate, varimax
    rotate, varimax blanks(.3)
    rotate, clear
    rotate, promax
    rotate, promax blanks(.3)
    rotate, clear

    estat common
* Save rotated FA loadings (RAW) to CSV
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
    export delimited using "$TAB\T_loadings_fa_raw_varimax.csv", replace
restore

	
    * FA scores (f1..f`ncomp_raw' -> *_raw)
    local fs_raw
    forvalues i=1/`ncomp_raw' {
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
    di as txt "---- STANDARDIZED (z_) variables ----"
    describe `Z'
    summarize `Z'
    corr `Z'

    * PCA (Z)
    pca `Z'
    screeplot, name(G_scree_z, replace)
    screeplot, yline(1) name(G_scree_z_y1, replace)
    pca `Z', mineigen(1)
    pca `Z', comp(`ncomp_z')
    pca `Z', comp(`ncomp_z') blanks(.3)

    rotate, varimax
    rotate, varimax blanks(.3)
    rotate, clear
    rotate, promax
    rotate, promax blanks(.3)
    rotate, clear

    estat loadings
	* Save rotated PCA loadings (Z) to CSV
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
    export delimited using "$TAB\T_loadings_pca_z_varimax.csv", replace
restore


    * PCA scores (pc1..pc`ncomp_z' -> *_z)
    local pcs_z
    forvalues i=1/`ncomp_z' {
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
    factor `Z', factor(`ncomp_z') pcf

    rotate, varimax
    rotate, varimax blanks(.3)
    rotate, clear
    rotate, promax
    rotate, promax blanks(.3)
    rotate, clear

    estat common
	
	* Save rotated FA loadings (Z) to CSV
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
    export delimited using "$TAB\T_loadings_fa_z_varimax.csv", replace
restore


    * FA scores (f1..f`ncomp_z' -> *_z)
    local fs_z
    forvalues i=1/`ncomp_z' {
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

*******************************************************
* C) Save tidy score files — robust to missing scores
*******************************************************

* RAW scores
local have_raw ""
capture unab have_raw : pc*_raw f*_raw
if !_rc {
    preserve
        keep mdate `have_raw'
        compress
        save "$PROC\fa_pca_scores_raw.dta", replace
        di as result "✅ Saved: $PROC\fa_pca_scores_raw.dta"
    restore
} else {
    di as txt "ℹ️ No RAW scores found (pc*_raw / f*_raw) — skipping save."
}

* Z scores
local have_z ""
capture unab have_z : pc*_z f*_z
if !_rc {
    preserve
        keep mdate `have_z'
        compress
        save "$PROC\fa_pca_scores_z.dta", replace
        di as result "✅ Saved: $PROC\fa_pca_scores_z.dta"
    restore
} 
else {
    di as txt "ℹ️ No Z scores found (pc*_z / f*_z) — skipping save."
}

*---------------------------------------------------------------
* 6) Single–Factor "AEB-like" index (RAW and Z)
*    Idea: load all vars onto ONE factor; orient the score to move
*          in the same direction as AEB (or z_AEB_aeb for Z set).
*    Outputs:
*      $PROC\fa1_scores_raw.dta   (mdate + AEB_like_raw)
*      $PROC\fa1_scores_z.dta     (mdate + AEB_like_z)
*      $TAB\T_loadings_fa1_raw.csv , $TAB\T_loadings_fa1_z.csv
*---------------------------------------------------------------

* --- Ensure AEB variables exist for sign orientation ---
capture confirm variable AEB_aeb
if _rc {
    di as err "AEB_aeb not found for RAW orientation."; 
}
capture confirm variable z_AEB_aeb
if _rc {
    di as err "z_AEB_aeb not found for Z orientation.";
}

* --- If Section 5 locals are not in scope, (re)define hard-coded lists ---
if "`Xraw'" == "" {
    local Xraw ///
        News_Based_Policy_Uncert_Index ///
        AEB_aeb ///
        SI_Dem_part ///
        SI_Ind_part ///
        SI_Rep_part ///
        CCI_monthly_cci ///
        UMCSENT_umc ///
        EOM_Close_de ///
        GEPU_current_gepu ///
        TPU_tpu ///
        SBOI_nfibo ///
        UncertaintyIndex_nfibu ///
        vix_vix
}
if "`Z'" == "" {
    local Z ///
        z_News_Based_Policy_Uncert_Index ///
        z_AEB_aeb ///
        z_SI_Dem_part ///
        z_SI_Ind_part ///
        z_SI_Rep_part ///
        z_CCI_monthly_cci ///
        z_UMCSENT_umc ///
        z_EOM_Close_de ///
        z_GEPU_current_gepu ///
        z_TPU_tpu ///
        z_SBOI_nfibo ///
        z_UncertaintyIndex_nfibu ///
        z_vix_vix
}

*******************************************************
* A) RAW single-factor (factor(1)), orient to AEB_aeb
*******************************************************
quietly {
    local p_raw : word count `Xraw'
}
if `p_raw' >= 2 {
    di as txt "---- RAW single-factor (AEB-like) ----"

    * Extract ONE factor (no rotation needed for 1 factor)
    factor `Xraw', factor(1)

    * Save loadings to CSV
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
        export delimited using "$TAB\T_loadings_fa1_raw.csv", replace
    restore

    * Score and orient so that higher = more AEB-like
    capture drop F1_raw
    predict F1_raw
    capture confirm variable AEB_aeb
    if !_rc {
        quietly corr AEB_aeb F1_raw
        scalar s_raw = sign(r(rho))
        if s_raw < 0 replace F1_raw = -F1_raw
    }
    rename F1_raw AEB_like_raw
    order mdate AEB_like_raw, first

    * Save sidecar file
    preserve
        keep mdate AEB_like_raw
        compress
        save "$PROC\fa1_scores_raw.dta", replace
    restore
    di as result "✅ Saved: $PROC\fa1_scores_raw.dta"
} 
else di as txt "ℹ️ Skipping RAW single-factor (fewer than 2 vars)."

*******************************************************
* B) STANDARDIZED (Z) single-factor, orient to z_AEB_aeb
*******************************************************
quietly {
    local p_z : word count `Z'
}
if `p_z' >= 2 {
    di as txt "---- Z single-factor (AEB-like) ----"

    factor `Z', factor(1)

    * Save loadings to CSV
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
        export delimited using "$TAB\T_loadings_fa1_z.csv", replace
    restore

    * Score and orient to positive correlation with z_AEB_aeb
    capture drop F1_z
    predict F1_z
    capture confirm variable z_AEB_aeb
    if !_rc {
        quietly corr z_AEB_aeb F1_z
        scalar s_z = sign(r(rho))
        if s_z < 0 replace F1_z = -F1_z
    }
    rename F1_z AEB_like_z
    order mdate AEB_like_z, first

    * Save sidecar file
    preserve
        keep mdate AEB_like_z
        compress
        save "$PROC\fa1_scores_z.dta", replace
    restore
    di as result "✅ Saved: $PROC\fa1_scores_z.dta"
}
else di as txt "ℹ️ Skipping Z single-factor (fewer than 2 vars)."

*******************************************************
* C) (Optional) Export the single-factor series to CSVs
*******************************************************
foreach which in raw z {
    capture confirm file "$PROC\fa1_scores_`which'.dta"
    if !_rc {
        preserve
            use "$PROC\fa1_scores_`which'.dta", clear
            export delimited using "$TAB\T_fa1_scores_`which'.csv", replace
        restore
        di as result "✅ Saved: $TAB\T_fa1_scores_`which'.csv"
    }
}

clear
*******************************************************
* 6) OLS: pc_AEB_aeb ~ all other percent-change vars
*******************************************************

* Load pct-change dataset
use "$PROC\firstcorrsv2_post2014_pctchg.dta", clear
format mdate %tm
tsset mdate, monthly

* --- Hard-coded varlists ---
local y  pc_AEB_aeb

* Exclude mdate and pc_mdate; include all other pc_* vars
local X  ///
    pc_News_Based_Policy_Uncert_Inde ///
    pc_SI_Dem_part ///
    pc_SI_Ind_part ///
    pc_SI_Rep_part ///
    pc_CCI_monthly_cci ///
    pc_UMCSENT_umc ///
    pc_EOM_Close_de ///
    pc_GEPU_current_gepu ///
    pc_TPU_tpu ///
    pc_SBOI_nfibo ///
    pc_UncertaintyIndex_nfibu ///
    pc_vix_vix ///
	S1 ///
	C1 ///
	S2 ///
	C2 ///
	

* Quick check: drop AEB from RHS if it slipped in
local X : list X - `y'

* Run regression (HC1 robust SEs)
regress `y' `X', vce(robust)

* (Optional) multicollinearity check
estat vif


*---------------------------------------------------------------
* 7) Random Forest (rforest) on percent-change data (with seasonality & lags)
*     y = pc_AEB_aeb ; X = other pc_* vars + Fourier seasonality + target lags
*     - time-aware split: last 24 months = TEST
*     - handles missing y by restricting training/eval to observed y
*     - prints OOB RMSE, TEST MAE/RMSE, TEST R^2, TEST MAPE
*     - compares to naive-1M and seasonal-12M baselines
*     - saves predictions and variable importance to PROC/TAB
*---------------------------------------------------------------

* Load data and declare monthly time
use "$PROC\firstcorrsv2_post2014_pctchg.dta", clear
format mdate %tm
tsset mdate, monthly

* ---- Hard-coded varlists ----
local y  pc_AEB_aeb
local X  ///
    pc_News_Based_Policy_Uncert_Inde ///
    pc_SI_Dem_part pc_SI_Ind_part pc_SI_Rep_part ///
    pc_CCI_monthly_cci pc_UMCSENT_umc pc_EOM_Close_de ///
    pc_GEPU_current_gepu pc_TPU_tpu pc_SBOI_nfibo ///
    pc_UncertaintyIndex_nfibu pc_vix_vix

* ensure y is not on RHS
local X : list X - `y'

* ---- Fourier seasonality (K=1; robust to gaps by using mdate) ----
capture drop s1 c1
gen double _w = 2*_pi/12
gen double s1 = sin(_w*mdate)
gen double c1 = cos(_w*mdate)
drop _w
local X `X' s1 c1

* (Optional K=2) uncomment next two lines if needed:
* gen double s2 = sin(2*(2*_pi/12)*mdate)
* gen double c2 = cos(2*(2*_pi/12)*mdate)
* local X `X' s2 c2

* ---- Target lags (often boosts RF performance) ----
foreach k in 1 2 3 {
    capture drop pc_AEB_aeb_l`k'
    gen pc_AEB_aeb_l`k' = L`k'.pc_AEB_aeb
}
local X `X' pc_AEB_aeb_l1 pc_AEB_aeb_l2 pc_AEB_aeb_l3

* ---- Train/Test split (last 24 months = test) ----
quietly summarize mdate
local tmax = r(max)
capture drop is_test
gen byte is_test = (mdate > `tmax' - 23)
label define split 0 "train" 1 "test"
label values is_test split

* ---- Train only where y is observed (rforest requirement) ----
capture drop ok_train ok_test_eval
gen byte ok_train     = (is_test==0 & !missing(pc_AEB_aeb))
gen byte ok_test_eval = (is_test==1 & !missing(pc_AEB_aeb))

* ---- Install/verify rforest ----
cap which rforest
if _rc ssc install rforest, replace
set seed 12345

* Choose numvars ~ sqrt(#features)
local p : word count `X'
local mtry = ceil(sqrt(`p'))

* ---- Fit Random Forest (regression) ----
rforest pc_AEB_aeb `X' if ok_train, ///
    type(reg) iterations(800) numvars(`mtry') lsize(5) depth(0) seed(12345)

di as txt    "OOB RMSE (train build): " %9.4f e(OOB_Error)

* ---- Predict on TEST and print RF metrics ----
capture drop yhat_rf
predict yhat_rf if is_test==1
di as result "TEST MAE (RF):  " %9.4f e(MAE)  "   TEST RMSE (RF): " %9.4f e(RMSE)

* ---- TEST R^2 and MAPE ----
preserve
    keep if ok_test_eval
    gen resid = pc_AEB_aeb - yhat_rf
    quietly summarize pc_AEB_aeb
    local ybar = r(mean)
    gen sst = (pc_AEB_aeb - `ybar')^2
    gen sse = resid^2
    quietly summarize sse
    local SSE = r(sum)
    quietly summarize sst
    local SST = r(sum)
    di as result "TEST R^2 (RF): " %6.3f (1 - `SSE'/`SST')
    gen ape = 100*abs(resid/pc_AEB_aeb) if pc_AEB_aeb!=0
    quietly summarize ape
    di as result "TEST MAPE (RF): " %6.2f r(mean) " %"
restore

* ---- Baseline comparisons (Naive-1M & Seasonal-12M) ----
capture drop y_naive1 y_seas12
gen y_naive1 = L1.pc_AEB_aeb
gen y_seas12 = L12.pc_AEB_aeb

preserve
    keep if ok_test_eval
    * Naive-1M RMSE
    gen err1 = pc_AEB_aeb - y_naive1
    gen se1  = err1^2
    quietly summarize se1
    di as txt "Naive-1M RMSE:     " %9.4f sqrt(r(mean))
    * Seasonal-12M RMSE
    gen err12 = pc_AEB_aeb - y_seas12
    gen se12  = err12^2
    quietly summarize se12
    di as txt "Seasonal-12M RMSE: " %9.4f sqrt(r(mean))
restore

* ---- Save predictions (all rows) ----
preserve
    keep mdate pc_AEB_aeb yhat_rf is_test
    compress
    save "$PROC\rf_preds_pctchg_rforest.dta", replace
    export delimited using "$TAB\T_rf_preds_pctchg_rforest.csv", replace
restore
di as result "✅ Saved: $PROC\rf_preds_pctchg_rforest.dta and $TAB\T_rf_preds_pctchg_rforest.csv"

* ---- Variable importance -> CSV ----
capture matrix list e(importance)
if !_rc {
    matrix VI = e(importance)
    preserve
        clear
        svmat double VI, names(col)
        gen variable = ""
        local rn : rownames VI
        local i = 1
        foreach r of local rn {
            replace variable = "`r'" in `i'
            local ++i
        }
        order variable
        compress
        save "$PROC\rf_varimp_pctchg_rforest.dta", replace
        export delimited using "$TAB\T_rf_varimp_pctchg_rforest.csv", replace
    restore
    di as result "✅ Saved: $PROC\rf_varimp_pctchg_rforest.dta and $TAB\T_rf_varimp_pctchg_rforest.csv"
} 
else di as txt "ℹ️ e(importance) not found; skipping var-importance export."
