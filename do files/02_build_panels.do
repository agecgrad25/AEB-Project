*02_Build_panels
*******************************************************
* AEB monthly workflow: correlations → FA (Q1–Q5) → event-study
* Expects macros: $PROC, $TABLES, $FIGURES
*******************************************************

*-----------------------------
* Load panel & time index
*-----------------------------
use "$PROC\aeb_monthly.dta", clear
tsset mdate, monthly

*-----------------------------
* 2.1 Descriptives & correlations (AEB, ICC, IFE)
*-----------------------------
summ AEB ICC IFE

pwcorr AEB ICC IFE, sig star(0.10)

* Export correlation table (no mat() option needed)
cap which esttab
if _rc di as error "NOTE: estout/esttab not installed. Run: ssc install estout"
estpost corr AEB ICC IFE, listwise
esttab using "$TABLES\T1_aeb_correlations.csv", ///
    replace unstack nonumber nomtitles b(3) label


*=====================================================
* Utility: Manual Bartlett's test of sphericity
* Usage:  bartlett_sphericity varlist
*=====================================================
capture program drop bartlett_sphericity
program define bartlett_sphericity
    version 19
    syntax varlist(min=2)

    preserve
        keep `varlist'
        drop if missing(`varlist')
        quietly count
        scalar n = r(N)
        if (n < 3) {
            di as error "Not enough complete cases for Bartlett's test."
            restore
            exit 2000
        }

        corr `varlist'
        matrix R = r(C)
        local p = colsof(R)

        * det(R) via eigenvalues (product of evals)
        matrix symeigen eval V = R
        scalar minEig = eval[1,`p']
        if (minEig <= 0) {
            di as error "Correlation matrix not positive definite; Bartlett chi2 undefined."
            restore
            exit 459
        }

        scalar logdet = 0
        forvalues i = 1/`p' {
            scalar logdet = logdet + ln(eval[1,`i'])
        }

        scalar df   = `p' * (`p' - 1) / 2
        scalar chi2 = - ( n - 1 - (2*`p' + 5)/6 ) * logdet
        scalar pval = chi2tail(df, chi2)

        di as txt "Bartlett test of sphericity: chi2(" df ") = " ///
            %9.3f chi2 ",  p = " %6.4f pval "   [n=" n ", p=" `p' "]"
    restore
end

*-----------------------------
* 2.2 Factor analysis on Q1–Q5 (force 2 factors if present & usable)
*-----------------------------
capture noisily confirm variable Q1 Q2 Q3 Q4 Q5
if _rc {
    di as result "Q1–Q5 not found — skipping factor analysis."
}
else {
    preserve
        keep mdate Q1 Q2 Q3 Q4 Q5
        drop if missing(Q1, Q2, Q3, Q4, Q5)

        * Drop any constant items in the estimation sample
        local items Q1 Q2 Q3 Q4 Q5
        local keepitems
        foreach v of local items {
            quietly summarize `v'
            if (r(sd)>1e-8 & r(N)>0) local keepitems `keepitems' `v'
            else di as error "Dropping `v' from FA (constant or no variance)."
        }
        if "`keepitems'" == "" {
            di as error "No non-constant items remain — skipping FA."
        }
        else {
            * Ensure clean spacing so varlist parses correctly
            local keepitems : list retoken keepitems

            * Factor model: FORCE two factors
            factor `keepitems', pcf factors(2) corr

            * KMO (informative, not required)
            estat kmo

            * ===== Manual Bartlett's test (prints χ2, df, p) =====
           * ===== Manual Bartlett's test (no preserve; always prints) =====
tempvar touse
mark `touse'
markout `touse' `keepitems'     // drop rows with any missing item from the mark

quietly corr `keepitems' if `touse'
matrix Rb = r(C)
scalar n   = r(N)
local p    = colsof(Rb)

matrix symeigen eval Vb = Rb
scalar minEig = eval[1,`p']
if (minEig <= 0) {
    di as error "Bartlett chi2 undefined (R not positive definite)."
}
else {
    scalar logdet = 0
    forvalues i = 1/`p' {
        scalar logdet = logdet + ln(eval[1,`i'])
    }
    scalar df   = `p' * (`p' - 1) / 2
    scalar chi2 = - ( n - 1 - (2*`p' + 5)/6 ) * logdet
    scalar pval = chi2tail(df, chi2)
    di as txt "Bartlett test of sphericity: chi2(" df ") = " ///
        %9.3f chi2 ",  p = " %6.4f pval "   [n=" n ", p=" `p' "]"
}
* ===== End manual Bartlett =====


           * Rotate obliquely (factors may correlate)
rotate, promax

* --- Export rotated PATTERN loadings (from e(L)) ---
matrix Lpat = e(L)
mat rownames Lpat = `keepitems'
mat colnames Lpat = F1 F2
cap which esttab
if !_rc {
    esttab matrix(Lpat) using "$TABLES\T2_FA_loadings_Q1_Q5_pattern.csv", ///
        replace nonumber nomtitles b(3)
}

            * Score factors (Bartlett) and label
            predict F1 F2, bartlett
            label var F1 "Factor 1 (Bartlett) — Q1–Q5"
            label var F2 "Factor 2 (Bartlett) — Q1–Q5"

            * Save scored subset to a tempfile
            tempfile fa_scores
            save `fa_scores', replace
        }
    restore

    * Merge factor scores back to the full panel (bring in all F* if present)
    capture confirm file `fa_scores'
    if !_rc {
        merge 1:1 mdate using `fa_scores', keepusing(F*) nogen
    }
    else {
        di as result "No factor score file found — skipping merge."
    }
}




*-----------------------------
* 2.3 Event-study plots (face validity)
*-----------------------------
twoway (tsline AEB), ///
    title("AEB over time") ///
    xtitle("") ytitle("Index") ///
    xline(`=ym(2016,11)' `=ym(2019,1)' `=ym(2020,3)', lpattern(dash)) ///
    note("Verticals: 2016 U.S. election | Trade tensions | COVID onset")
graph export "$FIGURES\F1_aeb_timeseries.png", replace width(1600)

twoway (tsline ICC IFE), ///
    title("ICC & IFE over time") ///
    xtitle("") ytitle("Index") ///
    xline(`=ym(2016,11)' `=ym(2019,1)' `=ym(2020,3)', lpattern(dash))
graph export "$FIGURES\F2_icc_ife_timeseries.png", replace width(1600)








