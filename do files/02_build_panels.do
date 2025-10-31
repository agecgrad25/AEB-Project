*******************************************************
* 02_Build_panels — FACTOR ANALYSIS ONLY (Stata 19)
* Outputs: CSV tables to $TABLES, figures to $FIGURES
*******************************************************

* ---------- Setup ----------
* Define output folders if not pre-defined; create if missing
capture confirm global TABLES
if _rc global TABLES "tables"
capture confirm global FIGURES
if _rc global FIGURES "figures"
cap mkdir "$TABLES"
cap mkdir "$FIGURES"

* Ensure esttab/estpost available for table exports
cap which esttab
if _rc di as error "NOTE: estout/esttab not installed. Run: ssc install estout"

* ---------- Load data ----------
use "$PROC/aeb_monthly.dta", clear
tsset mdate, monthly

* ---------- Guard: ensure items exist and keep complete cases ----------
capture noisily confirm variable Q1 Q2 Q3 Q4 Q5
if _rc {
    di as error "Q1–Q5 not found — skipping factor analysis."
    exit 111
}

preserve
    keep mdate Q1 Q2 Q3 Q4 Q5
    drop if missing(Q1, Q2, Q3, Q4, Q5)

    * Drop any constant items in the estimation sample
    local items Q1 Q2 Q3 Q4 Q5
    local keepitems
    foreach v of local items {
        quietly summarize `v'
        if (r(sd)>1e-8 & r(N)>0) local keepitems `keepitems' `v'
        else di as error "Dropping `v' (constant/no variance in FA sample)."
    }
    if "`keepitems'" == "" {
        di as error "No usable items remain — aborting FA."
        restore
        exit 498
    }
    * Retokenize for safety
    local keepitems : list retoken keepitems

    * ---------- Factor model (FORCE 2 factors), KMO ----------
    factor `keepitems', pcf factors(2) corr
    estat kmo

    * ---------- Manual Bartlett's test (no preserve) ----------
    tempvar touse
    mark `touse'
    markout `touse' `keepitems'      // drop rows with any missing item from the mark
    quietly corr `keepitems' if `touse'
    matrix Rb = r(C)
    scalar n   = r(N)
    local  p   = colsof(Rb)
    matrix symeigen evalRb VRb = Rb
    scalar minEig = evalRb[1,`p']
    if (minEig <= 0) {
        di as error "Bartlett chi2 undefined (R not positive definite)."
    }
    else {
        scalar logdet = 0
        forvalues i = 1/`p' {
            scalar logdet = logdet + ln(evalRb[1,`i'])
        }
        scalar df   = `p' * (`p' - 1) / 2
        scalar chi2 = - ( n - 1 - (2*`p' + 5)/6 ) * logdet
        scalar pval = chi2tail(df, chi2)
        di as txt "Bartlett test of sphericity: chi2(" df ") = " ///
            %9.3f chi2 ",  p = " %6.4f pval "   [n=" n ", p=" `p' "]"
    }

    * ---------- Rotate and capture results ----------
    rotate, promax

    * PATTERN loadings (regression-style) right after rotate
    matrix Lpat = e(L)
    mat rownames Lpat = `keepitems'

    * Factor names from k = number of columns in Lpat
    local k = colsof(Lpat)
    local fnames
    forvalues j=1/`k' {
        local fnames `fnames' F`j'
    }
    mat colnames Lpat = `fnames'

    * Phi (factor correlations). If absent, use I(k)
    capture matrix Phi = e(Phi)
    if (_rc) | (rowsof(Phi)!=`k') | (colsof(Phi)!=`k') {
        matrix Phi = I(`k')
    }
    mat rownames Phi = `fnames'
    mat colnames Phi = `fnames'

    * STRUCTURE loadings (= correlations) : Lstr = Lpat * Phi
    matrix Lstr = Lpat * Phi
    mat rownames Lstr = `keepitems'
    mat colnames Lstr = `fnames'

  * Communality = rowwise sum of squares of Lstr
tempname Cmn
matrix `Cmn' = Lstr
matrix `Cmn' = `Cmn''*`Cmn'    // this gives k x k, not what we want (so don't use)

* Better: loop rows to build a p×1 vector
local p = rowsof(Lstr)
local k = colsof(Lstr)
matrix Cmn = J(`p',1,.)
forvalues i=1/`p' {
    scalar rsq = 0
    forvalues j=1/`k' {
        scalar rsq = rsq + (Lstr[`i',`j']^2)
    }
    matrix Cmn[`i',1] = rsq
}
mat rownames Cmn = `keepitems'
mat colnames Cmn = Communality

matrix Unq = J(rowsof(Cmn), 1, 1)
matrix Unq = Unq - Cmn
mat rownames Unq = `keepitems'
mat colnames Unq = Uniqueness

cap which esttab
if !_rc {
    esttab matrix(Cmn) using "$TABLES/T_FA_communalities.csv",  replace nonumber nomtitles b(3)
    esttab matrix(Unq) using "$TABLES/T_FA_uniqueness.csv",    replace nonumber nomtitles b(3)
}



    * ---------- Export tables ----------
    cap which esttab
    if !_rc {
        esttab matrix(Lpat) using "$TABLES/T_FA_loadings_pattern.csv", ///
            replace nonumber nomtitles b(3)
        esttab matrix(Lstr) using "$TABLES/T_FA_loadings_structure.csv", ///
            replace nonumber nomtitles b(3)
        esttab matrix(Phi)  using "$TABLES/T_FA_factor_correlations.csv", ///
            replace nonumber nomtitles b(3)
        capture confirm matrix Unq
        if !_rc esttab matrix(Unq) using "$TABLES/T_FA_uniqueness.csv", replace b(3)
        capture confirm matrix Cmn
        if !_rc esttab matrix(Cmn) using "$TABLES/T_FA_communalities.csv", replace b(3)
    }
    else di as error "esttab not installed; tables not exported."
	
* ---------- Orthogonal rotation (VARIMAX) on same 2-factor model ----------
* Refit clean (unrotated) 2-factor model on the same keepitems
factor `keepitems', pcf factors(2) corr
estat kmo

* Varimax rotation (orthogonal)
rotate, varimax

* PATTERN loadings after varimax (for orthogonal rotations, pattern = structure)
matrix Lpat_v = e(L)
mat rownames Lpat_v = `keepitems'

* Factor names
local kv = colsof(Lpat_v)
local fnames_v
forvalues j=1/`kv' {
    local fnames_v `fnames_v' F`j'
}
mat colnames Lpat_v = `fnames_v'

* Phi (factor correlations) is identity under orthogonal rotation
matrix Phi_v = I(`kv')
mat rownames Phi_v = `fnames_v'
mat colnames Phi_v = `fnames_v'

* STRUCTURE loadings (= correlations); equals Lpat_v for orthogonal rotations
matrix Lstr_v = Lpat_v
mat rownames Lstr_v = `keepitems'
mat colnames Lstr_v = `fnames_v'

* Communalities and uniquenesses
local pv = rowsof(Lstr_v)
local kv2 = colsof(Lstr_v)
matrix Cmn_v = J(`pv',1,.)
forvalues i=1/`pv' {
    scalar rsq_v = 0
    forvalues j=1/`kv2' {
        scalar rsq_v = rsq_v + (Lstr_v[`i',`j']^2)
    }
    matrix Cmn_v[`i',1] = rsq_v
}
mat rownames Cmn_v = `keepitems'
mat colnames Cmn_v = Communality

matrix Unq_v = J(rowsof(Cmn_v), 1, 1)
matrix Unq_v = Unq_v - Cmn_v
mat rownames Unq_v = `keepitems'
mat colnames Unq_v = Uniqueness

* ---------- Export (VARIMAX) ----------
cap which esttab
if !_rc {
    esttab matrix(Lpat_v) using "$TABLES/T_FA_loadings_pattern_varimax.csv",   replace nonumber nomtitles b(3)
    esttab matrix(Lstr_v) using "$TABLES/T_FA_loadings_structure_varimax.csv", replace nonumber nomtitles b(3)
    esttab matrix(Phi_v)  using "$TABLES/T_FA_factor_correlations_varimax.csv", replace nonumber nomtitles b(3)
    esttab matrix(Cmn_v)  using "$TABLES/T_FA_communalities_varimax.csv",      replace nonumber nomtitles b(3)
    esttab matrix(Unq_v)  using "$TABLES/T_FA_uniqueness_varimax.csv",         replace nonumber nomtitles b(3)
}

* ---------- Tidy item→factor assignments (VARIMAX, 2-factor) ----------
clear
input str32 item
"Q1"
"Q2"
"Q3"
"Q4"
"Q5"
end

tempname Mv
matrix `Mv' = Lstr_v
svmat double `Mv', names(col)   // creates F1 F2 for varimax

gen double absF1 = abs(F1)
gen double absF2 = abs(F2)
gen str8  primary_factor  = cond(absF1 >= absF2, "F1", "F2")
gen double primary_loading = cond(absF1 >= absF2, F1,  F2)

local thresh = 0.30
local thresh_tag : display %02.0f 100*`thresh'
gen byte loads_ge_`thresh_tag' = (abs(primary_loading) >= `thresh')

gen byte cross_loading = (abs(absF1 - absF2) < 0.15)

order item F1 F2 primary_factor primary_loading loads_ge_`thresh_tag' cross_loading
format F1 F2 primary_loading %6.3f

export delimited using "$TABLES/T_FA_item_factor_assignments_varimax.csv", replace


   * ----- Tidy table (2-factor case) -----

    clear
    input str32 item
    "Q1"
    "Q2"
    "Q3"
    "Q4"
    "Q5"
    end

    * Bring STRUCTURE loadings into variables F1 and F2
    tempname M
    matrix `M' = Lstr
    svmat double `M', names(col)   // creates F1 F2

    * Primary factor = which absolute loading is larger
    gen double absF1 = abs(F1)
    gen double absF2 = abs(F2)
    gen str8  primary_factor  = cond(absF1 >= absF2, "F1", "F2")
    gen double primary_loading = cond(absF1 >= absF2, F1,  F2)

	* Threshold and a safe tag for the variable name
local thresh = 0.30
local thresh_tag : display %02.0f 100*`thresh'   // -> "30"

* Create the flag variable with a valid name
gen loads_ge_`thresh_tag' = (abs(primary_loading) >= `thresh')

	
   * ---- Flags: threshold + cross-loading (safe to re-run) ----
local thresh = 0.30
local thresh_tag : display %02.0f 100*`thresh'   // -> "30"

* Create-or-replace threshold flag
capture confirm variable loads_ge_`thresh_tag'
if !_rc {
    replace loads_ge_`thresh_tag' = (abs(primary_loading) >= `thresh')
}
else {
    gen byte loads_ge_`thresh_tag' = (abs(primary_loading) >= `thresh')
}

* Create-or-replace cross-loading flag (|F1|-|F2| within 0.15)
capture confirm variable cross_loading
if !_rc {
    replace cross_loading = (abs(absF1 - absF2) < 0.15)
}
else {
    gen byte cross_loading = (abs(absF1 - absF2) < 0.15)
}

* --- Self-healing tidy table footer (ensures all vars exist, then orders safely) ---

* 1) Ensure item list exists (from Lstr rownames)
capture confirm variable item
if _rc {
    local items : rowfullnames Lstr
    local n : word count `items'
    count
    if (r(N)==0) set obs `n'
    gen str32 item = "" in 1/`n'
    forvalues i=1/`n' {
        replace item = "``i''" in `i'
    }
}

* 2) Ensure F1/F2 exist (structure loadings into vars)
capture confirm variable F1
if _rc {
    tempname M
    matrix `M' = Lstr
    svmat double `M', names(col)   // creates F1 F2 (given 2 factors)
}

* 3) Ensure derived vars exist (primary_* and absF*)
capture confirm variable absF1
if _rc gen double absF1 = abs(F1)
capture confirm variable absF2
if _rc gen double absF2 = abs(F2)

capture confirm variable primary_factor
if _rc gen str8 primary_factor  = cond(absF1 >= absF2, "F1", "F2")

capture confirm variable primary_loading
if _rc gen double primary_loading = cond(absF1 >= absF2, F1,  F2)

* 4) Threshold flag with valid name (idempotent)
local thresh = 0.30
local thresh_tag : display %02.0f 100*`thresh'   // -> "30"

capture confirm variable loads_ge_`thresh_tag'
if !_rc {
    replace loads_ge_`thresh_tag' = (abs(primary_loading) >= `thresh')
}
else {
    gen byte loads_ge_`thresh_tag' = (abs(primary_loading) >= `thresh')
}

* 5) Cross-loading flag (idempotent)
capture confirm variable cross_loading
if !_rc {
    replace cross_loading = (abs(absF1 - absF2) < 0.15)
}
else {
    gen byte cross_loading = (abs(absF1 - absF2) < 0.15)
}

* 6) Order only the variables that actually exist (avoids errors)
local want item F1 F2 primary_factor primary_loading loads_ge_`thresh_tag' cross_loading
local have
foreach v of local want {
    capture confirm variable `v'
    if !_rc local have `have' `v'
}
capture noisily order `have'
format F1 F2 primary_loading %6.3f

* 7) Export
export delimited using "$TABLES/T_FA_item_factor_assignments.csv", 


* ---------- Scree / Elbow figures (no pca/screeplot dependency) ----------

* Use the items you actually factored. If Lstr exists, grab rownames from it.
capture confirm matrix Lstr
if _rc {
    * fallback: hard-code if needed
    local items Q1 Q2 Q3 Q4 Q5
}
else {
    local items : rowfullnames Lstr
}

* Eigenvalues from correlation matrix
corr `items'
matrix R = r(C)
matrix symeigen ev V = R

* Build a tiny dataset for plotting
clear
set obs `=colsof(ev)'
gen comp  = _n
gen eigen = .
forvalues i=1/`=colsof(ev)' {
    replace eigen = ev[1,`i'] in `i'
}

* Scree with Kaiser line
twoway (connected eigen comp, msym(circle)), ///
    title("Scree plot — Q1–Q5") ///
    xtitle("Component") ytitle("Eigenvalue") ///
    yline(1, lpattern(dash)) ///
    name(G_scree_items, replace)
graph export "$FIGURES/G_scree_items.png", replace width(1600)

* Elbow highlight (largest drop)
gen drop = eigen - eigen[_n+1]
quietly summarize drop in 1/`=_N-1', meanonly
local maxdrop = r(max)
gen byte is_elbow = (drop==`maxdrop' & _n<_N)

twoway ///
 (connected eigen comp, msym(circle)) ///
 (scatter   eigen comp if is_elbow, msym(diamond)), ///
 title("Eigenvalues (Elbow marked)") ///
 xtitle("Component") ytitle("Eigenvalue") yline(1, lpattern(dash)) ///
 name(G_scree_elbow, replace)
graph export "$FIGURES/G_scree_elbow.png", replace width(1600)




