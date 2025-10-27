*******************************************************
* 03_analysis_simple_v19.do — simple (Stata 19)
* Input: $PROC\aebcorrsv3.dta  (post-2014, +corn/soy prices & vols, +Seasons)*******************************************************

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

* Identify numeric vars; exclude mdate and SEASON dummies (no Fourier terms exist)
ds, has(type numeric)
local allnum `r(varlist)'
local allnum : list allnum - mdate
ds se_*, has(type numeric)
local seasons `r(varlist)'

* Keep a "no-season" numeric set for transforms/analyses
local nums_noseason : list allnum - seasons

*******************************************************
* 1) Z-score all numeric vars (exclude mdate, seasons, AEB alias)
*******************************************************
local zbase `nums_noseason'
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
* 2) Correlations — include seasons + export full matrix
**************************************************************
* Build clean lists
ds se_*, has(type numeric)
local seasons `r(varlist)'

* Non-season numerics (used for z_ counterparts in pairwise)
local others `nums_noseason'
local others : list others - mdate
local others : list others - AEB
local others : list others - AEB_aeb

* Keep only vars that actually have z_ counterparts
local safe_others
foreach v of local others {
    capture confirm variable z_`v'
    if !_rc local safe_others `safe_others' `v'
}

* ------------ 2A) Pairwise: corr(z_AEB, z_* or se_*) ------------
preserve
    tempfile out
    tempname ph
    postfile `ph' str64 var double rho int N using "`out'", replace

    * Non-season variables (use z_ versions)
    foreach v of local safe_others {
        quietly corr z_AEB z_`v'
        post `ph' ("`v'") (r(rho)) (r(N))
    }

    * Season dummies (use raw se_* with z_AEB)
    if "`seasons'" != "" {
        foreach s of local seasons {
            quietly corr z_AEB `s'
            post `ph' ("`s'") (r(rho)) (r(N))
        }
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
restore

* ------------ 2B) Full correlation matrix (incl. seasons) ------------
* Build matrix varlist: AEB + (non-season numerics without AEB) + seasons
local tmp `nums_noseason'
local tmp : list tmp - AEB
local tmp : list tmp - AEB_aeb
local CMAT "AEB `tmp' `seasons'"

quietly corr `CMAT'
matrix C = r(C)

* --- Make column names legal, <=32 chars, and UNIQUE to satisfy svmat
local cols : colnames C
local safe ""
local used ""
local i = 0
foreach c of local cols {
    local ++i
    local base = strtoname("`c'")
    if strlen("`base'")>28 local base = substr("`base'",1,28)
    local nm "`base'"
    local k = 1
    while strpos(" `used' "," `nm' ") {
        local suffix _`k'
        local slen : length local suffix
        local blen = 32 - `slen'
        if `blen' < 1 local blen = 1
        local nm = substr("`base'",1,`blen')
        local nm "`nm'`suffix'"
        local ++k
    }
    local used `used' `nm'
    local safe `safe' `nm'
}
matrix colnames C = `safe'

preserve
    clear
    svmat double C, names(col)
    gen variable = ""
    local rn : rownames C
    local i = 1
    foreach r of local rn {
        replace variable = "`r'" in `i'
        local ++i
    }
    order variable
    export delimited using "$TAB\T_corr_full`SUF'.csv", replace
restore

*****************************************************************
* 3) (Optional) Rolling 24-month correlations (exclude seasons)
*****************************************************************
cap which rangestat
if _rc ssc install rangestat, replace

preserve
    use "$PROC\aebcorrsv3.dta", clear
    tsset mdate, monthly

    ds, has(type numeric)
    local allnum `r(varlist)'
    local allnum : list allnum - mdate
    ds se_*, has(type numeric)
    local seasons `r(varlist)'
    local nums_noseason : list allnum - seasons

    * Rebuild z_ for these only (skip AEB alias to avoid duplicates)
    local zbase `nums_noseason'
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
    local others `nums_noseason'
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
capture unab rvars : r_*
if !_rc {
    drop `rvars'
}

**************************************************************
* 5) PCA + FA — auto-build varlists (include corn/sb vols)
**************************************************************
* Build RAW = all numeric except mdate, seasons, z_* and the AEB alias
ds, has(type numeric)
local Xraw `r(varlist)'
local Xraw : list Xraw - mdate
local Xraw : list Xraw - seasons
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
* 6) OLS: ΔAEB ~ all other first-difference vars (d_*) + seasons
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

* Candidate differenced predictors = all d_* vars
capture unab cand : d_*
if _rc {
    di as err "No d_* variables found."
    exit 111
}

* Exclusions: DV itself, d_mdate, any AEB variants
local drop " `y' d_mdate d_AEB_aeb d_AEB "
local X ""
foreach v of local cand {
    if strpos("`drop'"," `v' ")==0 {
        local X `X' `v'
    }
}

* Add seasonal dummies as controls (omit se_winter to avoid dummy trap)
ds se_*, has(type numeric)
local SE `r(varlist)'
local SE_nobase `SE'
local SE_nobase : list SE_nobase - se_winter

* Safety: ensure predictors remain
local p : word count `X'
if `p'==0 & "`SE_nobase'"=="" {
    di as err "No predictors left after exclusions."
    exit 111
}

* Fit OLS with robust SEs
regress `y' `X' `SE_nobase', vce(robust)

* Multicollinearity check
estat vif

* ---- Quick print line + seasons joint test ----
di as res "Kitchen-sink OLS  |  N=" %9.0f e(N) "  R2=" %6.4f e(r2) ///
          "  adjR2=" %6.4f e(r2_a) "  RMSE=" %9.3g e(rmse)
if "`SE_nobase'" != "" {
    quietly testparm `SE_nobase'
    di as res "Seasons joint p = " %6.4f r(p)
}

* ---- Export tidy coefficient table WITH model stats: $TAB\kitchen_sink.csv ----
tempfile KS
tempname PF

* Capture model-level stats once
scalar R2    = e(r2)
scalar AR2   = e(r2_a)
scalar NN    = e(N)
scalar RMSE  = e(rmse)
scalar P_SES = .
capture confirm local SE_nobase
if !_rc & "`SE_nobase'" != "" {
    quietly testparm `SE_nobase'
    scalar P_SES = r(p)
}

postfile `PF' str64 var str32 short double b se t p str1 star ///
    double r2 adjr2 N rmse p_seasons using "`KS'", replace

matrix b = e(b)
matrix V = e(V)
local cn : colnames b
foreach v of local cn {
    scalar bb = b[1, colnumb(b, "`v'")]
    scalar ss = sqrt(V[colnumb(V, "`v'"), colnumb(V, "`v'")])
    scalar tt = bb/ss
    scalar pp = 2*ttail(e(df_r), abs(tt))
    local star ""
    if (pp < .05) local star "*"

    * Short name mapping (same scheme as horse-race)
    local short "`v'"
    local short : subinstr local short "News_Based_Policy_Uncert_Index" "EPU_news", all
    local short : subinstr local short "CCI_monthly_cci"               "CCI",       all
    local short : subinstr local short "UMCSENT_umc"                   "UMCSENT",   all
    local short : subinstr local short "EOM_Close_de"                  "EOM",       all
    local short : subinstr local short "GEPU_current_gepu"             "GEPU",      all
    local short : subinstr local short "TPU_tpu"                       "TPU",       all
    local short : subinstr local short "SBOI_nfibo"                    "SBOI",      all
    local short : subinstr local short "UncertaintyIndex_nfibu"        "NFIBU",     all
    local short : subinstr local short "SI_Dem_part"                   "SI_Dem",    all
    local short : subinstr local short "SI_Ind_part"                   "SI_Ind",    all
    local short : subinstr local short "SI_Rep_part"                   "SI_Rep",    all
    local short : subinstr local short "corn_vol_month"                "corn_mo",   all
    local short : subinstr local short "corn_vol_ann"                  "corn_ann",  all
	local short : subinstr local short "corn_close"               "corn_px",   all
    local short : subinstr local short "sb_vol_month"                  "soy_mo",    all
    local short : subinstr local short "sb_vol_ann"                    "soy_ann",   all
	  local short : subinstr local short "sb_close"                      "soy_px",    all
    local short : subinstr local short "vix_vix"                       "VIX",       all
    local short : subinstr local short "_cons"                         "Intercept", all
    local short : subinstr local short "se_spring"                     "Spring",    all
    local short : subinstr local short "se_summer"                     "Summer",    all
    local short : subinstr local short "se_fall"                       "Fall",      all
    local short : subinstr local short "se_winter"                     "Winter",    all

    post `PF' ("`v'") ("`short'") (bb) (ss) (tt) (pp) ("`star'") ///
        (R2) (AR2) (NN) (RMSE) (P_SES)
}
postclose `PF'

use "`KS'", clear␊
order var short b se t p star r2 adjr2 N rmse p_seasons
format b se %9.3g
format t %8.2f
format p p_seasons %6.4f
format r2 adjr2 %6.4f
format N %9.0f
format rmse %9.3g

export delimited using "$TAB\kitchen_sink.csv", replace



*******************************************************
* 6b) HORSE RACE — ΔAEB on each X one-at-a-time
*     (no seasons, N>=90, original order, star p<.05)
*******************************************************

* Load differenced data (consistent with later horse-race steps)
use "$PROC\aebcorrsv3_diff.dta", clear
format mdate %tm
tsset mdate, monthly

* DV (prefer d_AEB_aeb; fallback d_AEB)
local y ""
capture confirm variable d_AEB_aeb
if !_rc local y d_AEB_aeb
else {
    capture confirm variable d_AEB
    if !_rc local y d_AEB
}
if "`y'"=="" exit 111

* Build X list in DATASET ORDER: all d_* minus DV and d_mdate
ds d_*, has(type numeric)
local X `r(varlist)'
local X : list X - `y'
local X : list X - d_AEB_aeb
local X : list X - d_AEB
local X : list X - d_mdate
local K : word count `X'
if `K'==0 exit 111

* Collect results
tempfile HR
tempname PF
postfile `PF' str64 var str32 short double b se t p r2 N str1 star ///
    double const_b double const_se using "`HR'", replace


forvalues j = 1/`K' {
    local v : word `j' of `X'
    quietly regress `y' `v', vce(robust)
    if _rc | missing(e(N)) continue
    if (e(N) < 90) continue

    scalar b  = _b[`v']
    scalar se = _se[`v']
    scalar t  = b/se
    scalar p  = 2*ttail(e(df_r), abs(t))
    scalar r2 = e(r2)
    scalar NN = e(N)
	scalar b0  = _b[_cons]
	scalar se0 = _se[_cons]

    * Short names
    local short "`v'"
    local short : subinstr local short "News_Based_Policy_Uncert_Index" "EPU_news", all
    local short : subinstr local short "CCI_monthly_cci"               "CCI",       all
    local short : subinstr local short "UMCSENT_umc"                   "UMCSENT",   all
    local short : subinstr local short "EOM_Close_de"                  "EOM",       all
    local short : subinstr local short "GEPU_current_gepu"             "GEPU",      all
    local short : subinstr local short "TPU_tpu"                       "TPU",       all
    local short : subinstr local short "SBOI_nfibo"                    "SBOI",      all
    local short : subinstr local short "UncertaintyIndex_nfibu"        "NFIBU",     all
    local short : subinstr local short "SI_Dem_part"                   "SI_Dem",    all
    local short : subinstr local short "SI_Ind_part"                   "SI_Ind",    all
    local short : subinstr local short "SI_Rep_part"                   "SI_Rep",    all
    local short : subinstr local short "corn_vol_month"                "corn_mo",   all
    local short : subinstr local short "corn_vol_ann"                  "corn_ann",  all
	local short : subinstr local short "corn_close"                    "corn_px",   all
    local short : subinstr local short "sb_vol_month"                  "soy_mo",    all
    local short : subinstr local short "sb_vol_ann"                    "soy_ann",   all
	local short : subinstr local short "sb_close"                      "soy_px",    all
    local short : subinstr local short "vix_vix"                       "VIX",       all

    local star ""
    if (p < .05) local star "*"

post `PF' ("`v'") ("`short'") (b) (se) (t) (p) (r2) (NN) ("`star'") (b0) (se0)
}
postclose `PF'

use "`HR'", clear
order var short b se t p star r2 N const_b const_se 
format b se %9.3g
format const_b const_se %9.3g
format t %8.2f
format p %6.4f
format r2 %6.4f
format N  %9.0f

export delimited using "$TAB\horse race.csv", replace

*******************************************************
* 6c) HORSE RACE — non-season vars w/o season controls
*     + append season dummies as additional variables
*     (N>=90, original order, star p<.05, short names)
*******************************************************

* Load diffs and ensure seasons are available for the seasons-only model
use "$PROC\aebcorrsv3_diff.dta", clear
format mdate %tm
tsset mdate, monthly
capture ds se_*
if _rc | "`r(varlist)'"=="" {
    preserve
        tempfile _se
        use "$PROC\aebcorrsv3.dta", clear
        keep mdate se_*
        duplicates drop mdate, force
        save "`_se'"
    restore
    merge 1:1 mdate using "`_se'", nogen
}

* Dependent variable (prefer d_AEB_aeb; fallback d_AEB)
local y ""
capture confirm variable d_AEB_aeb
if !_rc local y d_AEB_aeb
else {
    capture confirm variable d_AEB
    if !_rc local y d_AEB
}
if "`y'"=="" {
    di as err "Horse-race: missing d_AEB_aeb/d_AEB."
    exit 111
}

* Candidate non-season predictors (dataset order): all d_* minus DV and d_mdate
ds d_*, has(type numeric)
local X `r(varlist)'
local X : list X - `y'
local X : list X - d_AEB_aeb
local X : list X - d_AEB
local X : list X - d_mdate
local K : word count `X'
if `K'==0 {
    di as err "Horse-race: no candidate d_* predictors."
    exit 111
}

* Seasons (for the seasons-only model); drop one base to avoid dummy trap
ds se_*, has(type numeric)
local SE `r(varlist)'
local SE_nobase `SE'
local SE_nobase : list SE_nobase - se_winter
local dropped "se_winter"
local nSE   : word count `SE'
local nSEnb : word count `SE_nobase'
if `nSE'==`nSEnb' & `nSE'>0 {
    local base : word 1 of `SE'
    local SE_nobase : list SE_nobase - `base'
    local dropped "`base'"
}

* Collector
tempfile HRs
tempname PFs
postfile `PFs' str64 var str32 short double b se t p r2 N str1 star ///
    double const_b double const_se using "`HRs'", replace

* ---- (A) One-at-a-time NON-SEASON regressions (NO season controls) ----
forvalues j = 1/`K' {
    local v : word `j' of `X'
    quietly regress `y' `v', vce(robust)
    if _rc | missing(e(N)) continue
    if (e(N) < 90) continue

    scalar b   = _b[`v']
    scalar se  = _se[`v']
    scalar t   = b/se
    scalar p   = 2*ttail(e(df_r), abs(t))
    scalar r2  = e(r2)
    scalar NN  = e(N)
    scalar b0  = _b[_cons]
    scalar se0 = _se[_cons]

    * Short names
    local short "`v'"
    local short : subinstr local short "News_Based_Policy_Uncert_Index" "EPU_news", all
    local short : subinstr local short "CCI_monthly_cci"               "CCI",       all
    local short : subinstr local short "UMCSENT_umc"                   "UMCSENT",   all
    local short : subinstr local short "EOM_Close_de"                  "EOM",       all
    local short : subinstr local short "GEPU_current_gepu"             "GEPU",      all
    local short : subinstr local short "TPU_tpu"                       "TPU",       all
    local short : subinstr local short "SBOI_nfibo"                    "SBOI",      all
    local short : subinstr local short "UncertaintyIndex_nfibu"        "NFIBU",     all
    local short : subinstr local short "SI_Dem_part"                   "SI_Dem",    all
    local short : subinstr local short "SI_Ind_part"                   "SI_Ind",    all
    local short : subinstr local short "SI_Rep_part"                   "SI_Rep",    all
    local short : subinstr local short "corn_vol_month"                "corn_mo",   all
    local short : subinstr local short "corn_vol_ann"                  "corn_ann",  all
	local short : subinstr local short "corn_close"                    "corn_px",   all
    local short : subinstr local short "sb_vol_month"                  "soy_mo",    all
    local short : subinstr local short "sb_vol_ann"                    "soy_ann",   all
	local short : subinstr local short "sb_close"                      "soy_px",    all
    local short : subinstr local short "vix_vix"                       "VIX",       all

    local star ""
    if (p < .05) local star "*"

    post `PFs' ("`v'") ("`short'") (b) (se) (t) (p) (r2) (NN) ("`star'") (b0) (se0)
}

* ---- (B) Append seasons as additional variables (from seasons-only model) ----
if "`SE_nobase'" != "" {
    quietly regress `y' `SE_nobase', vce(robust)
    if !_rc & e(N) >= 90 {
        scalar r2s  = e(r2)
        scalar NNs  = e(N)
        scalar b0s  = _b[_cons]
        scalar se0s = _se[_cons]

        foreach s of local SE_nobase {
            * human-friendly short label for season
            local sh "`s'"
            local sh : subinstr local sh "se_spring" "Spring", all
            local sh : subinstr local sh "se_summer" "Summer", all
            local sh : subinstr local sh "se_fall"   "Fall",   all
            local sh : subinstr local sh "se_winter" "Winter", all

            capture scalar bs  = _b[`s']
            capture scalar ses = _se[`s']
            if _rc continue
            scalar ts = bs/ses
            scalar ps = 2*ttail(e(df_r), abs(ts))

            local star ""
            if (ps < .05) local star "*"

            * NOTE: for season rows, use the season var/label and seasons-only stats
            post `PFs' ("`s'") ("`sh'") (bs) (ses) (ts) (ps) (r2s) (NNs) ("`star'") (b0s) (se0s)
        }
    }
}

postclose `PFs'

use "`HRs'", clear
order var short b se t p star r2 N const_b const_se
format b se %9.3g
format const_b const_se %9.3g
format t %8.2f
format p %6.4f
format r2 %6.4f
format N  %9.0f

export delimited using "$TAB\horse_race_seasons.csv", replace

*******************************************************
* 6d) KITCHEN SINK — NO seasonal dummies; drop DV diffs
*      + overfit diagnostics + tidy CSV export
*******************************************************

use "$PROC\aebcorrsv3_diff.dta", clear
format mdate %tm
tsset mdate, monthly

* ---- Dependent variable (prefer d_AEB_aeb; fallback d_AEB)
local y ""
capture confirm variable d_AEB_aeb
if !_rc local y d_AEB_aeb
else {
    capture confirm variable d_AEB
    if !_rc local y d_AEB
}
if "`y'"=="" {
    di as err "Kitchen sink (no seasons): missing d_AEB_aeb/d_AEB."
    exit 111
}

* ---- Candidate predictors: build from d_* while explicitly skipping DV aliases + d_mdate
ds d_*, has(type numeric)
local rawX `r(varlist)'
local X ""
foreach v of local rawX {
    if inlist("`v'", "`y'", "d_AEB_aeb", "d_AEB", "d_mdate") continue
    local X `X' `v'
}
local p0 : word count `X'
if `p0'==0 {
    di as err "Kitchen sink (no seasons): no candidate d_* predictors."
    exit 111
}

* ---- (Optional) remove zero-variance columns on the common sample
* Build common sample indicator across y and all X
tempvar T
gen byte `T' = !missing(`y')
foreach v of local X {
    replace `T' = `T' & !missing(`v')
}
count if `T'
local Ncommon = r(N)
* Drop predictors that are constant on the common sample
local Xclean
foreach v of local X {
    quietly summarize `v' if `T', meanonly
    if r(Var)>0 & r(N)>0 local Xclean `Xclean' `v'
}
local X "`Xclean'"

* Double-check: ensure the dependent variable (or its aliases) are not in X
local DV_present : list X & `y'
if "`DV_present'" != "" {
    di as err "Kitchen sink (no seasons): dropping DV `y' from predictors to avoid R2=1."
    local X : list X - `y'
}
* Also guard against alternate DV alias accidentally persisting
local DV_aliases "d_AEB_aeb d_AEB"
local drop_alias : list X & DV_aliases
if "`drop_alias'" != "" {
    local X : list X - `drop_alias'
}

* ---- Fit kitchen-sink (NO seasons), robust SEs
quietly regress `y' `X', vce(robust)

* ---- Quick diagnostics for overfitting
local p_used = e(rank) - 1   // parameters excluding intercept
local df_r   = e(df_r)
di as res "Kitchen-sink (NO seasons) | N=" %9.0f e(N) "  p_used≈" %3.0f `p_used' ///
          "  df_r=" %6.0f `df_r' "  R2=" %6.4f e(r2) "  adjR2=" %6.4f e(r2_a)
if (e(r2) >= .9999) di as err ">>> R2≈1: likely overfit or perfect linear combo present."
if (e(df_r) <= 5)   di as err ">>> Very low residual dof: p too large vs N (risk of overfit)."

* ---- Export tidy coefficient table WITH model stats: $TAB\kitchen_sink.csv
tempfile KS
tempname PF

scalar R2   = e(r2)
scalar AR2  = e(r2_a)
scalar NN   = e(N)
scalar RMSE = e(rmse)

postfile `PF' str64 var str32 short double b se t p str1 star ///
    double r2 adjr2 N rmse using "`KS'", replace

matrix b = e(b)
matrix V = e(V)
local cn : colnames b
foreach v of local cn {
    scalar bb = b[1, colnumb(b, "`v'")]
    scalar ss = sqrt(V[colnumb(V, "`v'"), colnumb(V, "`v'")])
    scalar tt = bb/ss
    scalar pp = 2*ttail(e(df_r), abs(tt))
    local star ""
    if (pp < .05) local star "*"

    * short names (same scheme as before)
    local short "`v'"
    local short : subinstr local short "News_Based_Policy_Uncert_Index" "EPU_news", all
    local short : subinstr local short "CCI_monthly_cci"               "CCI",       all
    local short : subinstr local short "UMCSENT_umc"                   "UMCSENT",   all
    local short : subinstr local short "EOM_Close_de"                  "EOM",       all
    local short : subinstr local short "GEPU_current_gepu"             "GEPU",      all
    local short : subinstr local short "TPU_tpu"                       "TPU",       all
    local short : subinstr local short "SBOI_nfibo"                    "SBOI",      all
    local short : subinstr local short "UncertaintyIndex_nfibu"        "NFIBU",     all
    local short : subinstr local short "SI_Dem_part"                   "SI_Dem",    all
    local short : subinstr local short "SI_Ind_part"                   "SI_Ind",    all
    local short : subinstr local short "SI_Rep_part"                   "SI_Rep",    all
    local short : subinstr local short "corn_vol_month"                "corn_mo",   all
    local short : subinstr local short "corn_vol_ann"                  "corn_ann",  all
	local short : subinstr local short "corn_close"                    "corn_px",   all
    local short : subinstr local short "sb_vol_month"                  "soy_mo",    all
    local short : subinstr local short "sb_vol_ann"                    "soy_ann",   all
	local short : subinstr local short "sb_close"                      "soy_px",    all
    local short : subinstr local short "vix_vix"                       "VIX",       all
    local short : subinstr local short "_cons"                         "Intercept", all

    post `PF' ("`v'") ("`short'") (bb) (ss) (tt) (pp) ("`star'") (R2) (AR2) (NN) (RMSE)
}
postclose `PF'

use "`KS'", clear
order var short b se t p star r2 adjr2 N rmse
format b se %9.3g
format t %8.2f
format p %6.4f
format r2 adjr2 %6.4f
format N %9.0f
format rmse %9.3g
export delimited using "$TAB\kitchen_sink_noseasons.csv", replace


*******************************************************
* 6e) HORSE RACE — with seasonal dummies as controls
*     (N>=90, original order, star p<.05, short names)
*******************************************************

* Always load the diff dataset for robustness
use "$PROC\aebcorrsv3_diff.dta", clear
format mdate %tm
tsset mdate, monthly

* If seasons aren't here, merge them in from levels file
capture ds se_*
if _rc | "`r(varlist)'"=="" {
    preserve
        tempfile _se
        use "$PROC\aebcorrsv3.dta", clear
        keep mdate se_*
        duplicates drop mdate, force
        save "`_se'"
    restore
    merge 1:1 mdate using "`_se'", nogen
}

* ---- Dependent variable (prefer d_AEB_aeb; fallback d_AEB)
local y ""
capture confirm variable d_AEB_aeb
if !_rc local y d_AEB_aeb
else {
    capture confirm variable d_AEB
    if !_rc local y d_AEB
}
if "`y'"=="" {
    di as err "Horse-race (seasons): could not find d_AEB_aeb or d_AEB. Build the diff file first."
    exit 111
}

* ---- Build X list in dataset order: all d_* minus DV and d_mdate
ds d_*, has(type numeric)
local rawX `r(varlist)'
local X ""
foreach v of local rawX {
    if inlist("`v'", "`y'", "d_AEB_aeb", "d_AEB", "d_mdate") continue
    local X `X' `v'
}
local K : word count `X'
if `K'==0 {
    di as err "Horse-race (seasons): no candidate d_* predictors found."
    exit 111
}

* ---- Season controls (drop one base to avoid dummy trap)
ds se_*, has(type numeric)
local SE `r(varlist)'
local SE_nobase `SE'
local SE_nobase : list SE_nobase - se_winter
local nSE   : word count `SE'
local nSEnb : word count `SE_nobase'
if `nSE'==`nSEnb' & `nSE'>0 {
    local base : word 1 of `SE'
    local SE_nobase : list SE_nobase - `base'
}

* ---- Collect results
tempfile HRs
tempname PFs
postfile `PFs' str64 var str32 short double b se t p r2 N str1 star ///
    double const_b double const_se double p_seasons using "`HRs'", replace

forvalues j = 1/`K' {
    local v : word `j' of `X'

    quietly regress `y' `v' `SE_nobase', vce(robust)
    if _rc | missing(e(N)) continue
     if (e(N) < 90) continue

    scalar pS = .
    if "`SE_nobase'" != "" {
        quietly testparm `SE_nobase'
        scalar pS = r(p)
    }
    scalar b   = _b[`v']
    scalar se  = _se[`v']
    scalar t   = b/se
    scalar p   = 2*ttail(e(df_r), abs(t))
    scalar r2  = e(r2)
    scalar NN  = e(N)
    scalar b0  = .
    capture scalar b0 = _b[_cons]
    scalar se0 = .
    capture scalar se0 = _se[_cons]

    * Short names
    local short "`v'"
    local short : subinstr local short "News_Based_Policy_Uncert_Index" "EPU_news", all
    local short : subinstr local short "CCI_monthly_cci"               "CCI",       all
    local short : subinstr local short "UMCSENT_umc"                   "UMCSENT",   all
    local short : subinstr local short "EOM_Close_de"                  "EOM",       all
    local short : subinstr local short "GEPU_current_gepu"             "GEPU",      all
    local short : subinstr local short "TPU_tpu"                       "TPU",       all
    local short : subinstr local short "SBOI_nfibo"                    "SBOI",      all
    local short : subinstr local short "UncertaintyIndex_nfibu"        "NFIBU",     all
    local short : subinstr local short "SI_Dem_part"                   "SI_Dem",    all
    local short : subinstr local short "SI_Ind_part"                   "SI_Ind",    all
    local short : subinstr local short "SI_Rep_part"                   "SI_Rep",    all
    local short : subinstr local short "corn_vol_month"                "corn_mo",   all
    local short : subinstr local short "corn_vol_ann"                  "corn_ann",  all
	    local short : subinstr local short "corn_close"                    "corn_px",   all
    local short : subinstr local short "sb_vol_month"                  "soy_mo",    all
    local short : subinstr local short "sb_vol_ann"                    "soy_ann",   all
	    local short : subinstr local short "sb_close"                      "soy_px",    all
    local short : subinstr local short "vix_vix"                       "VIX",       all

    local star ""
    if (p < .05) local star "*"

	  post `PFs' ("`v'") ("`short'") (b) (se) (t) (p) (r2) (NN) ("`star'") (b0) (se0) (pS)
}
postclose `PFs'

use "`HRs'", clear
order var short b se t p star r2 N const_b const_se p_seasons
format b se %9.3g
format const_b const_se %9.3g
format t %8.2f
format p %6.4f
format r2 %6.4f
format N  %9.0f

export delimited using "$TAB\horse_race_ctrlseasons.csv", replace

*******************************************************
* 6f) HORSE RACE — with seasonal dummies as controls
*     + append seasons as additional variables (rows)
*     (N>=90, original order, star p<.05, short names)
*******************************************************
clear all
* Load diffs and make sure seasons exist
use "$PROC\aebcorrsv3_diff.dta", clear
format mdate %tm
tsset mdate, monthly

capture ds se_*
if _rc | "`r(varlist)'"=="" {
    preserve
        tempfile _se
        use "$PROC\aebcorrsv3.dta", clear
        keep mdate se_*
        duplicates drop mdate, force
        save "`_se'"
    restore
    merge 1:1 mdate using "`_se'", nogen
}

* Dependent variable (prefer d_AEB_aeb; fallback d_AEB)
local y ""
capture confirm variable d_AEB_aeb
if !_rc local y d_AEB_aeb
else {
    capture confirm variable d_AEB
    if !_rc local y d_AEB
}
if "`y'"=="" {
    di as err "Horse-race (seasons): missing d_AEB_aeb/d_AEB."
    exit 111
}

* Candidate X in dataset order: all d_* minus DV and d_mdate
ds d_*, has(type numeric)
local rawX `r(varlist)'
local X ""
foreach v of local rawX {
    if inlist("`v'", "`y'", "d_AEB_aeb", "d_AEB", "d_mdate") continue
    local X `X' `v'
}
local K : word count `X'
if `K'==0 {
    di as err "Horse-race (seasons): no candidate d_* predictors."
    exit 111
}

* Seasons (drop one base to avoid dummy trap)
ds se_*, has(type numeric)
local SE `r(varlist)'
local SE_nobase `SE'
local SE_nobase : list SE_nobase - se_winter
local dropped "se_winter"
local nSE   : word count `SE'
local nSEnb : word count `SE_nobase'
if `nSE'==`nSEnb' & `nSE'>0 {
    local base : word 1 of `SE'
    local SE_nobase : list SE_nobase - `base'
    local dropped "`base'"
}

* Results collector
tempfile HRs
tempname PFs
postfile `PFs' str64 var str32 short double b se t p r2 N str1 star ///
    double const_b double const_se using "`HRs'", replace

* ---- (A) One-at-a-time predictors WITH season controls ----
forvalues j = 1/`K' {
    local v : word `j' of `X'
    quietly regress `y' `v' `SE_nobase', vce(robust)
    if _rc | missing(e(N)) continue
    if (e(N) < 90) continue

    scalar b   = _b[`v']
    scalar se  = _se[`v']
    scalar t   = b/se
    scalar p   = 2*ttail(e(df_r), abs(t))
    scalar r2  = e(r2)
    scalar NN  = e(N)
    scalar b0  = .
    capture scalar b0 = _b[_cons]
    scalar se0 = .
    capture scalar se0 = _se[_cons]

    * Short names
    local short "`v'"
    local short : subinstr local short "News_Based_Policy_Uncert_Index" "EPU_news", all
    local short : subinstr local short "CCI_monthly_cci"               "CCI",       all
    local short : subinstr local short "UMCSENT_umc"                   "UMCSENT",   all
    local short : subinstr local short "EOM_Close_de"                  "EOM",       all
    local short : subinstr local short "GEPU_current_gepu"             "GEPU",      all
    local short : subinstr local short "TPU_tpu"                       "TPU",       all
    local short : subinstr local short "SBOI_nfibo"                    "SBOI",      all
    local short : subinstr local short "UncertaintyIndex_nfibu"        "NFIBU",     all
    local short : subinstr local short "SI_Dem_part"                   "SI_Dem",    all
    local short : subinstr local short "SI_Ind_part"                   "SI_Ind",    all
    local short : subinstr local short "SI_Rep_part"                   "SI_Rep",    all
    local short : subinstr local short "corn_vol_month"                "corn_mo",   all
    local short : subinstr local short "corn_vol_ann"                  "corn_ann",  all
    local short : subinstr local short "corn_close"                    "corn_px",   all
    local short : subinstr local short "sb_vol_month"                  "soy_mo",    all
    local short : subinstr local short "sb_vol_ann"                    "soy_ann",   all
    local short : subinstr local short "sb_close"                      "soy_px",    all
    local short : subinstr local short "vix_vix"                       "VIX",       all

    local star ""
    if (p < .05) local star "*"

    post `PFs' ("`v'") ("`short'") (b) (se) (t) (p) (r2) (NN) ("`star'") (b0) (se0)
}

* ---- (B) Append the season dummies as additional variables (rows) ----
if "`SE_nobase'" != "" {
    quietly regress `y' `SE_nobase', vce(robust)
    if !_rc & e(N)>=90 {
        scalar r2s = e(r2)
        scalar NNs = e(N)
        scalar b0s = .
        capture scalar b0s = _b[_cons]
        scalar se0s = .
        capture scalar se0s = _se[_cons]
        foreach s of local SE_nobase {
            local sh "`s'"
            local sh : subinstr local sh "se_spring" "Spring", all
            local sh : subinstr local sh "se_summer" "Summer", all
            local sh : subinstr local sh "se_fall"   "Fall",   all
            local sh : subinstr local sh "se_winter" "Winter", all

            capture scalar bs  = _b[`s']
            capture scalar ses = _se[`s']
            if _rc continue
            scalar ts = bs/ses
            scalar ps = 2*ttail(e(df_r), abs(ts))
            local star ""
            if (ps < .05) local star "*"
            post `PFs' ("`s'") ("`sh'") (bs) (ses) (ts) (ps) (r2s) (NNs) ("`star'") (b0s) (se0s)
        }
    }
}

postclose `PFs'

use "`HRs'", clear
order var short b se t p star r2 N const_b const_se
format b se %9.3g
format const_b const_se %9.3g
format t %8.2f
format p %6.4f
format r2 %6.4f
format N  %9.0f

export delimited using "$TAB\horse_race_snlctrl+seasons.csv", replace
