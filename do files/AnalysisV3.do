*******************************************************
* 03_analysis_simple_v19.do — simple (Stata 19) - REFACTORED
* Input: $PROC\aebcorrsv3.dta  (post-2014, +corn/soy prices & vols, +Seasons)
*******************************************************

version 19.0
clear all
set more off

*******************************************************
* CONFIGURATION & CONSTANTS
*******************************************************
local SUF "_v3"
local MIN_N 90
local ALPHA 0.05
local ROLL_WINDOW 24
local BLANK_THRESH 0.3
local N_COMPONENTS 3

* Verify required global macros are defined
if "$PROC" == "" {
    di as err "Global macro PROC is not defined. Please set it before running this script."
    exit 111
}
if "$TAB" == "" {
    di as err "Global macro TAB is not defined. Please set it before running this script."
    exit 111
}
if "$FIG" == "" {
    di as err "Global macro FIG is not defined. Please set it before running this script."
    exit 111
}

*******************************************************
* HELPER PROGRAMS
*******************************************************

* Program: map_short_names
* Maps long variable names to short display names
capture program drop map_short_names
program define map_short_names
    syntax, var(string) [Return(string)]

    local short "`var'"
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
    local short : subinstr local short "se_spring"                     "Spring",    all
    local short : subinstr local short "se_summer"                     "Summer",    all
    local short : subinstr local short "se_fall"                       "Fall",      all
    local short : subinstr local short "se_winter"                     "Winter",    all

    if "`return'" != "" {
        c_local `return' "`short'"
    }
end

* Program: export_matrix
* Exports a Stata matrix to CSV with row names as first column
capture program drop export_matrix
program define export_matrix
    syntax, matrix(string) file(string)

    preserve
        clear
        svmat double `matrix', names(col)
        gen variable = ""
        local rn : rownames `matrix'
        local i = 1
        foreach r of local rn {
            replace variable = "`r'" in `i'
            local ++i
        }
        order variable
        export delimited using "`file'", replace
    restore
end

* Program: get_clean_varlist
* Gets numeric variables excluding specified types
capture program drop get_clean_varlist
program define get_clean_varlist
    syntax, Return(string) [EXclude(string) INClude(string) NOSeasons NOZ NOAEB]

    * Start with all numeric or specified inclusion pattern
    if "`include'" != "" {
        ds `include', has(type numeric)
    }
    else {
        ds, has(type numeric)
    }
    local varlist `r(varlist)'

    * Always exclude mdate
    local varlist : list varlist - mdate

    * Exclude seasons if requested
    if "`noseasons'" != "" {
        ds se_*, has(type numeric)
        if !_rc {
            local seasons `r(varlist)'
            local varlist : list varlist - seasons
        }
    }

    * Exclude z_* variables if requested
    if "`noz'" != "" {
        capture unab zvars : z_*
        if !_rc {
            local varlist : list varlist - zvars
        }
    }

    * Exclude AEB variants if requested
    if "`noaeb'" != "" {
        local aeb_vars "AEB_aeb AEB z_AEB_aeb z_AEB d_AEB_aeb d_AEB"
        local varlist : list varlist - aeb_vars
    }

    * Apply custom exclusions
    if "`exclude'" != "" {
        local varlist : list varlist - exclude
    }

    c_local `return' "`varlist'"
end

* Program: build_aebcorrsv3pca
* Helper: ensure an AEB-free snapshot exists for RAW PCA/FA routines
capture program drop build_aebcorrsv3pca
program define build_aebcorrsv3pca
    syntax

    preserve
        quietly use "$PROC\aebcorrsv3.dta", clear

        * Remove AEB-specific series and any standardized z_* remnants
        capture drop AEB_aeb AEB z_*

        compress
        save "$PROC\aebcorrsv3pca.dta", replace
    restore
end

* Program: run_pca_fa
* Parameterized PCA and Factor Analysis with rotation and export
capture program drop run_pca_fa
program define run_pca_fa
    syntax varlist, suffix(string) [Ncomp(integer 3) BLanks(real 0.3) ///
        EXportpca(string) EXportfa(string) NOScores]

    local p : word count `varlist'
    if `p' < 2 {
        di as txt "Insufficient variables for PCA/FA (need at least 2)"
        exit
    }

    * Determine number of components
    local ncomp = cond(`p' >= `ncomp', `ncomp', `p')

    quietly describe `varlist'
    quietly summarize `varlist'
    quietly corr `varlist'

    **** PCA ****
    pca `varlist'
    screeplot, name(G_scree`suffix', replace)
    screeplot, yline(1) name(G_scree`suffix'_y1, replace)
    pca `varlist', mineigen(1)
    pca `varlist', comp(`ncomp')
    pca `varlist', comp(`ncomp') blanks(`blanks')
    rotate, varimax
    rotate, varimax blanks(`blanks')
    rotate, clear
    rotate, promax
    rotate, promax blanks(`blanks')
    rotate, clear

    estat loadings
    matrix L_pca`suffix' = e(L)

    * Export PCA loadings if requested
    if "`exportpca'" != "" {
        export_matrix, matrix(L_pca`suffix') file("`exportpca'")
    }

    * Generate PCA scores if requested
    if "`noscores'" == "" {
        local pcs
        forvalues i = 1/`ncomp' {
            local pcs `pcs' pc`i'
        }
        capture drop `pcs'
        predict `pcs', score
        foreach v of local pcs {
            rename `v' `v'`suffix'
        }
    }

    estat kmo

    **** FACTOR ANALYSIS ****
    factor `varlist'
    screeplot, name(G_scree`suffix'_fa, replace)
    screeplot, yline(1) name(G_scree`suffix'_fa_y1, replace)
    factor `varlist', mineigen(1)
    factor `varlist', factor(`ncomp')
    factor `varlist', factor(`ncomp') blanks(`blanks')
    rotate, varimax
    rotate, varimax blanks(`blanks')
    rotate, clear
    rotate, promax
    rotate, promax blanks(`blanks')
    rotate, clear

    estat common
    matrix L_fa`suffix' = e(L)

    * Export FA loadings if requested
    if "`exportfa'" != "" {
        export_matrix, matrix(L_fa`suffix') file("`exportfa'")
    }

    * Generate FA scores if requested
    if "`noscores'" == "" {
        local fs
        forvalues i = 1/`ncomp' {
            local fs `fs' f`i'
        }
        capture drop `fs'
        predict `fs'
        foreach v of local fs {
            rename `v' `v'`suffix'
        }
    }

    * Reliability + Bartlett
    alpha `varlist'
    cap which factortest
    if _rc ssc install factortest, replace
    factortest `varlist'
end

* Program: horse_race_regression
* Run one-at-a-time regressions for all predictors
capture program drop horse_race_regression
program define horse_race_regression
    syntax varlist, Yvar(string) file(string) ///
        [SEasons(string) MINn(integer 90) ALpha(real 0.05) APpendseasons]

    * Count predictors
    local K : word count `varlist'
    if `K' == 0 {
        di as err "No predictors for horse race"
        exit 111
    }

    * Setup postfile structure based on options
    tempfile HR
    tempname PF

    if "`seasons'" != "" {
        postfile `PF' str64 var str32 short double b se t p r2 N str1 star ///
            double const_b double const_se double p_seasons using "`HR'", replace
    }
    else {
        postfile `PF' str64 var str32 short double b se t p r2 N str1 star ///
            double const_b double const_se using "`HR'", replace
    }

    * Run regressions for each predictor
    forvalues j = 1/`K' {
        local v : word `j' of `varlist'
        quietly regress `yvar' `v' `seasons', vce(robust)
        if _rc | missing(e(N)) continue
        if (e(N) < `minn') continue

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

        * Test seasons if present
        scalar pS = .
        if "`seasons'" != "" {
            quietly testparm `seasons'
            scalar pS = r(p)
        }

        * Get short name
        map_short_names, var("`v'") return(short)

        local star ""
        if (p < `alpha') local star "*"

        if "`seasons'" != "" {
            post `PF' ("`v'") ("`short'") (b) (se) (t) (p) (r2) (NN) ("`star'") (b0) (se0) (pS)
        }
        else {
            post `PF' ("`v'") ("`short'") (b) (se) (t) (p) (r2) (NN) ("`star'") (b0) (se0)
        }
    }

    * Append season dummies as separate rows if requested
    if "`appendseasons'" != "" & "`seasons'" != "" {
        quietly regress `yvar' `seasons', vce(robust)
        if !_rc & e(N) >= `minn' {
            scalar r2s  = e(r2)
            scalar NNs  = e(N)
            scalar b0s  = .
            capture scalar b0s = _b[_cons]
            scalar se0s = .
            capture scalar se0s = _se[_cons]

            foreach s of local seasons {
                map_short_names, var("`s'") return(sh)

                capture scalar bs  = _b[`s']
                capture scalar ses = _se[`s']
                if _rc continue
                scalar ts = bs/ses
                scalar ps = 2*ttail(e(df_r), abs(ts))

                local star ""
                if (ps < `alpha') local star "*"

                if "`seasons'" != "" {
                    post `PF' ("`s'") ("`sh'") (bs) (ses) (ts) (ps) (r2s) (NNs) ("`star'") (b0s) (se0s) (.)
                }
                else {
                    post `PF' ("`s'") ("`sh'") (bs) (ses) (ts) (ps) (r2s) (NNs) ("`star'") (b0s) (se0s)
                }
            }
        }
    }

    postclose `PF'

    * Format and export results
    use "`HR'", clear
    if "`seasons'" != "" {
        order var short b se t p star r2 N const_b const_se p_seasons
    }
    else {
        order var short b se t p star r2 N const_b const_se
    }
    format b se const_b const_se %9.3g
    format t %8.2f
    format p %6.4f
    format r2 %6.4f
    format N %9.0f
    if "`seasons'" != "" {
        format p_seasons %6.4f
    }

    export delimited using "`file'", replace
end

*******************************************************
* SECTION 0: Load and set time
*******************************************************
use "$PROC\aebcorrsv3.dta", clear
capture confirm variable mdate
if _rc {
    di as err "mdate missing"
    exit 459
}
format mdate %tm
tsset mdate, monthly

* Ensure AEB_aeb exists
capture confirm variable AEB_aeb
if _rc {
    di as err "AEB_aeb not found in the dataset."
    exit 111
}

* Pre-OLS safeguards: drop any lingering standardized AEB
capture drop z_AEB_aeb

* Build AEB-free snapshot for PCA/FA
quietly build_aebcorrsv3pca

*******************************************************
* SECTION 1: Z-score all numeric vars (exclude mdate, seasons, AEB_aeb)
*******************************************************
get_clean_varlist, return(zbase) noseasons noaeb
foreach v of local zbase {
    capture drop z_`v'
    quietly egen z_`v' = std(`v')
}

* Snapshot without AEB_aeb for PCA/FA work
tempfile aebcorrsv3pca
preserve
    drop AEB_aeb
    save `aebcorrsv3pca'
restore

*******************************************************
* SECTION 2: Correlations
*******************************************************

* Get variable lists
get_clean_varlist, return(nums_noseason) noseasons
ds se_*, has(type numeric)
local seasons `r(varlist)'

* Build safe list (only vars with z_ counterparts)
local others `nums_noseason'
local others : list others - AEB_aeb
local safe_others
foreach v of local others {
    capture confirm variable z_`v'
    if !_rc local safe_others `safe_others' `v'
}

* 2A) Pairwise correlations
preserve
    tempfile out
    tempname ph
    postfile `ph' str64 var double rho int N using "`out'", replace

    * Non-season variables
    foreach v of local safe_others {
        quietly corr AEB_aeb `v'
        post `ph' ("`v'") (r(rho)) (r(N))
    }

    * Season dummies
    if "`seasons'" != "" {
        foreach s of local seasons {
            quietly corr AEB_aeb `s'
            post `ph' ("`s'") (r(rho)) (r(N))
        }
    }

    postclose `ph'

    use "`out'", clear
    drop if missing(rho)

    quietly count
    if r(N) == 0 {
        di as txt "No valid pairwise correlations between AEB_aeb and the candidate variables."
    }
    else {
        gen double abs_rho = abs(rho)
        gsort -abs_rho
        order var rho N
        format rho %6.3f
        format N %9.0f
        export delimited using "$TAB\T_corr_pairwise`SUF'.csv", replace
        drop abs_rho
    }
restore

* 2B) Full correlation matrix (incl. seasons)
local tmp `nums_noseason'
local tmp : list tmp - AEB_aeb
local CMAT "`tmp' `seasons'"

if trim("`CMAT'") == "" {
    di as txt "No variables available for correlation matrix once AEB_aeb is excluded."
}
else {
    quietly corr `CMAT'
    matrix C = r(C)

    * Make column names legal, <=32 chars, and UNIQUE
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

    export_matrix, matrix(C) file("$TAB\T_corr_full`SUF'.csv")
}

*******************************************************
* SECTION 3: Rolling 24-month correlations (exclude seasons)
*******************************************************
cap which rangestat
if _rc ssc install rangestat, replace

preserve
    use "$PROC\aebcorrsv3.dta", clear
    tsset mdate, monthly

    get_clean_varlist, return(nums_noseason) noseasons

    * Rebuild z_ for these only
    foreach v of local nums_noseason {
        capture drop z_`v'
        quietly egen z_`v' = std(`v')
    }

    * Rolling targets: must have z_ partner, exclude mdate and AEB_aeb
    local others `nums_noseason'
    local others : list others - AEB_aeb

    local safe_others
    foreach v of local others {
        capture confirm variable z_`v'
        if !_rc {
            local safe_others `safe_others' `v'
        }
    }

    * Rolling correlations
    foreach v of local safe_others {
        tempvar prod mx my mxy sx sy
        gen double `prod' = z_AEB_aeb * z_`v'
        rangestat (mean) `mx'=z_AEB_aeb `my'=z_`v' `mxy'=`prod' ///
                  (sd)   `sx'=z_AEB_aeb `sy'=z_`v', interval(mdate -`ROLL_WINDOW'+1 0)

        local base = strtoname("`v'")
        local maxbase = 32 - `=strlen("r_AEB_aeb_")'
        if strlen("`base'") > `maxbase' {
            local base = substr("`base'", 1, `maxbase')
        }
        local rname = "r_AEB_aeb_`base'"
        capture confirm variable `rname'
        local k = 1
        while !_rc {
            local rname = "r_AEB_aeb_`base'_`k'"
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

*******************************************************
* SECTION 4: Quick chart for one rolling series
*******************************************************
use "$PROC\monthly_rollingcorrs_simple`SUF'.dta", clear
local target r_AEB_aeb_News_Based_Policy_Uncert_Index
capture confirm variable `target'
if _rc {
    capture unab _cands : r_AEB_aeb_*
    if !_rc {
        local target : word 1 of `_cands'
    }
}
capture confirm variable `target'
if !_rc {
    twoway tsline `target', ///
        title("Rolling `ROLL_WINDOW'm corr: AEB_aeb vs selected") ytitle("corr") xtitle("")
    graph export "$FIG\F3_roll_AEBaeb_selected`SUF'.png", replace width(1600)
}
capture unab rvars : r_*
if !_rc {
    drop `rvars'
}

*******************************************************
* SECTION 5: PCA + FA
*******************************************************

* Get variable lists
get_clean_varlist, return(Xraw) noseasons noz noaeb
capture unab Z : z_*
if _rc local Z ""
local Z : list Z - z_AEB_aeb

local p_z : word count `Z'
local ncomp_z = cond(`p_z' >= `N_COMPONENTS', `N_COMPONENTS', `p_z')

**** 5A) RAW (nominal) variables ****
preserve
    use "$PROC\aebcorrsv3pca.dta", clear

    * Build raw list from snapshot
    get_clean_varlist, return(Xsnap) noseasons noz noaeb

    local p_snap : word count `Xsnap'
    if `p_snap' < 2 {
        di as txt "No variables available for RAW PCA/FA snapshot once exclusions applied."
    }
    else {
        local ncomp_snap = cond(`p_snap' >= `N_COMPONENTS', `N_COMPONENTS', `p_snap')

        run_pca_fa `Xsnap', suffix(_raw) ncomp(`ncomp_snap') blanks(`BLANK_THRESH') ///
            exportpca("$TAB\T_loadings_pca_raw_varimax`SUF'.csv") ///
            exportfa("$TAB\T_loadings_fa_raw_varimax`SUF'.csv")

        * Save raw PCA/FA scores
        local have_raw ""
        capture unab have_raw : pc*_raw f*_raw
        if !_rc {
            preserve
                keep mdate `have_raw'
                compress
                save "$PROC\fa_pca_scores_raw`SUF'.dta", replace
            restore
            di as result "✅ Saved: $PROC\fa_pca_scores_raw`SUF'.dta"
        }
    }
restore

if (`p_snap' < 2 & `p_z' < 2) {
    di as err "Not enough variables for PCA/FA."
    exit 111
}

**** 5B) STANDARDIZED (z_) variables ****
if `p_z' >= 2 {
    run_pca_fa `Z', suffix(_z) ncomp(`ncomp_z') blanks(`BLANK_THRESH') ///
        exportpca("$TAB\T_loadings_pca_z_varimax`SUF'.csv") ///
        exportfa("$TAB\T_loadings_fa_z_varimax`SUF'.csv")

    * Save z scores
    local have_z ""
    capture unab have_z : pc*_z f*_z
    if !_rc {
        preserve
            keep mdate `have_z'
            compress
            save "$PROC\fa_pca_scores_z`SUF'.dta", replace
        restore
    }
}

**** 5D) PCA + FA WITHOUT corn_close and sb_close (nopxs) ****
preserve
    use "$PROC\aebcorrsv3pca.dta", clear
    drop corn_close sb_close mdate se_spring se_summer se_fall se_winter

    ds
    local varlist_nopxs `r(varlist)'

    di as result "=== PCA/FA without corn_close and sb_close ==="
    di as result "Variables included: `varlist_nopxs'"

    tempfile nopxs_data
    save `nopxs_data', replace

    * Run PCA
    pca `varlist_nopxs', mineigen(1)
    rotate, varimax blanks(`BLANK_THRESH')
    estat loadings
    matrix L_pca_nopxs = e(L)
    export_matrix, matrix(L_pca_nopxs) file("$TAB\T_loadings_pca_nopxs`SUF'.csv")
    di as result "Saved: $TAB\T_loadings_pca_nopxs`SUF'.csv"

    * Reload for FA
    use `nopxs_data', clear
    factor `varlist_nopxs', mineigen(1)
    rotate, varimax blanks(`BLANK_THRESH')
    estat common
    matrix L_fa_nopxs = e(L)
    export_matrix, matrix(L_fa_nopxs) file("$TAB\T_loadings_fa_nopxs`SUF'.csv")
    di as result "Saved: $TAB\T_loadings_fa_nopxs`SUF'.csv"
restore

**** 5E) PCA + FA for Z-SCORED variables WITHOUT corn/soy prices (z_nopxs) ****
preserve
    use "$PROC\aebcorrsv3.dta", clear
    keep mdate z_*
    drop z_AEB_aeb z_corn_close z_sb_close mdate

    ds
    local varlist_z_nopxs `r(varlist)'

    di as result "=== PCA/FA for z-scored variables without corn/soy prices ==="
    di as result "Variables included: `varlist_z_nopxs'"

    tempfile z_nopxs_data
    save `z_nopxs_data', replace

    * Run PCA
    pca `varlist_z_nopxs', mineigen(1)
    rotate, varimax blanks(`BLANK_THRESH')
    estat loadings
    matrix L_pca_z_nopxs = e(L)
    export_matrix, matrix(L_pca_z_nopxs) file("$TAB\T_loadings_pca_z_nopxs`SUF'.csv")
    di as result "Saved: $TAB\T_loadings_pca_z_nopxs`SUF'.csv"

    * Reload for FA
    use `z_nopxs_data', clear
    factor `varlist_z_nopxs', mineigen(1)
    rotate, varimax blanks(`BLANK_THRESH')
    estat common
    matrix L_fa_z_nopxs = e(L)
    export_matrix, matrix(L_fa_z_nopxs) file("$TAB\T_loadings_fa_z_nopxs`SUF'.csv")
    di as result "Saved: $TAB\T_loadings_fa_z_nopxs`SUF'.csv"
restore

*******************************************************
* SECTION 6: Single–Factor "AEB-like" index (RAW and Z)
*******************************************************
local p_raw : word count `Xraw'
local p_z   : word count `Z'

* RAW single factor
if `p_raw' >= 2 {
    factor `Xraw', factor(1)
    matrix L_fa1_raw = e(L)
    export_matrix, matrix(L_fa1_raw) file("$TAB\T_loadings_fa1_raw`SUF'.csv")

    capture drop F1_raw
    predict F1_raw
    di as txt "Skipping F1_raw sign alignment with AEB_aeb prior to OLS."
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
    export_matrix, matrix(L_fa1_z) file("$TAB\T_loadings_fa1_z`SUF'.csv")

    capture drop F1_z
    predict F1_z
    di as txt "Skipping F1_z sign alignment with z_AEB_aeb prior to OLS."
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
* SECTION 7: OLS - Kitchen Sink with Seasons
*******************************************************
use "$PROC\aebcorrsv3_diff.dta", clear
format mdate %tm
tsset mdate, monthly

* Dependent variable
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

* Add seasonal dummies as controls (omit se_winter)
ds se_*, has(type numeric)
local SE `r(varlist)'
local SE_nobase `SE'
local SE_nobase : list SE_nobase - se_winter

* Safety check
local p : word count `X'
if `p'==0 & "`SE_nobase'"=="" {
    di as err "No predictors left after exclusions."
    exit 111
}

* Fit OLS with robust SEs
regress `y' `X' `SE_nobase', vce(robust)

* Multicollinearity check
estat vif

* Quick summary
di as res "Kitchen-sink OLS  |  N=" %9.0f e(N) "  R2=" %6.4f e(r2) ///
          "  adjR2=" %6.4f e(r2_a) "  RMSE=" %9.3g e(rmse)
if "`SE_nobase'" != "" {
    quietly testparm `SE_nobase'
    di as res "Seasons joint p = " %6.4f r(p)
}

* Export coefficient table
tempfile KS
tempname PF

scalar R2    = e(r2)
scalar AR2   = e(r2_a)
scalar NN    = e(N)
scalar RMSE  = e(rmse)
scalar P_SES = .
if "`SE_nobase'" != "" {
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
    if (pp < `ALPHA') local star "*"

    map_short_names, var("`v'") return(short)

    post `PF' ("`v'") ("`short'") (bb) (ss) (tt) (pp) ("`star'") ///
        (R2) (AR2) (NN) (RMSE) (P_SES)
}
postclose `PF'

use "`KS'", clear
order var short b se t p star r2 adjr2 N rmse p_seasons
format b se %9.3g
format t %8.2f
format p p_seasons %6.4f
format r2 adjr2 %6.4f
format N %9.0f
format rmse %9.3g
export delimited using "$TAB\kitchen_sink.csv", replace

*******************************************************
* SECTION 8: Horse Races
*******************************************************

* Load differenced data
use "$PROC\aebcorrsv3_diff.dta", clear
format mdate %tm
tsset mdate, monthly

* Dependent variable
local y ""
capture confirm variable d_AEB_aeb
if !_rc local y d_AEB_aeb
else {
    capture confirm variable d_AEB
    if !_rc local y d_AEB
}
if "`y'"=="" exit 111

* Build predictor list
ds d_*, has(type numeric)
local X `r(varlist)'
local X : list X - `y'
local X : list X - d_AEB_aeb
local X : list X - d_AEB
local X : list X - d_mdate

* Get seasons (if available, merge if needed)
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

ds se_*, has(type numeric)
local SE `r(varlist)'
local SE_nobase `SE'
local SE_nobase : list SE_nobase - se_winter

**** 8A) Horse Race - No Seasons ****
horse_race_regression `X', yvar(`y') file("$TAB\horse race.csv") minn(`MIN_N') alpha(`ALPHA')

**** 8B) Horse Race - No Season Controls + Append Seasons ****
horse_race_regression `X', yvar(`y') file("$TAB\horse_race_seasons.csv") ///
    minn(`MIN_N') alpha(`ALPHA') appendseasons

**** 8C) Horse Race - With Season Controls ****
horse_race_regression `X', yvar(`y') file("$TAB\horse_race_ctrlseasons.csv") ///
    seasons(`SE_nobase') minn(`MIN_N') alpha(`ALPHA')

**** 8D) Horse Race - With Season Controls + Append Seasons ****
horse_race_regression `X', yvar(`y') file("$TAB\horse_race_snlctrl+seasons.csv") ///
    seasons(`SE_nobase') minn(`MIN_N') alpha(`ALPHA') appendseasons

*******************************************************
* SECTION 9: Kitchen Sink - NO Seasonal Dummies
*******************************************************
use "$PROC\aebcorrsv3_diff.dta", clear
format mdate %tm
tsset mdate, monthly

* Dependent variable
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

* Build predictors
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

* Remove zero-variance columns
tempvar T
gen byte `T' = !missing(`y')
foreach v of local X {
    replace `T' = `T' & !missing(`v')
}
count if `T'
local Ncommon = r(N)

local Xclean
foreach v of local X {
    quietly summarize `v' if `T', meanonly
    if r(Var)>0 & r(N)>0 local Xclean `Xclean' `v'
}
local X "`Xclean'"

* Fit kitchen-sink (NO seasons)
quietly regress `y' `X', vce(robust)

* Quick diagnostics
local p_used = e(rank) - 1
local df_r   = e(df_r)
di as res "Kitchen-sink (NO seasons) | N=" %9.0f e(N) "  p_used≈" %3.0f `p_used' ///
          "  df_r=" %6.0f `df_r' "  R2=" %6.4f e(r2) "  adjR2=" %6.4f e(r2_a)
if (e(r2) >= .9999) di as err ">>> R2≈1: likely overfit or perfect linear combo present."
if (e(df_r) <= 5)   di as err ">>> Very low residual dof: p too large vs N (risk of overfit)."

* Export coefficient table
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
    if (pp < `ALPHA') local star "*"

    map_short_names, var("`v'") return(short)

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

di as result "🎉 Analysis complete!"
