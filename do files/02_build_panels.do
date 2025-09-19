02_build_panels

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

* Export correlation table (fixes your prior mat() error)
cap which esttab
if _rc {
    di as error "NOTE: estout/esttab not installed. Run: ssc install estout"
}
estpost corr AEB ICC IFE, listwise
esttab using "$TABLES\T1_aeb_correlations.csv", ///
    replace unstack nonumber nomtitles b(3) label

*-----------------------------
* 2.2 Factor analysis on Q1–Q5 (only if they exist & usable)
*-----------------------------
capture noisily confirm variable Q1 Q2 Q3 Q4 Q5
if _rc {
    di as result "Q1–Q5 not found — skipping factor analysis."
}
else {
    * Keep only rows where all five items are present (clean EFA sample)
    preserve
    keep mdate Q1 Q2 Q3 Q4 Q5
    drop if missing(Q1, Q2, Q3, Q4, Q5)

    * Guardrails: drop any constant item (zero variance in the estimation sample)
    local items Q1 Q2 Q3 Q4 Q5
    local keepitems
    foreach v of local items {
        quietly summarize `v'
        if (r(sd)>1e-8 & r(N)>0) local keepitems `keepitems' `v'
        else di as error "Dropping `v' from FA (constant or no variance in sample)."
    }
    if "`keepitems'" == "" {
        di as error "No non-constant items remain — skipping FA."
    }
    else {
        * Check singularity of correlation matrix
        corr `keepitems'
        matrix R = r(C)
        matrix symeigen ev V = R
        * Count eigenvalues > 1 (Kaiser); fall back to 1 if none
        local nf = 0
        forvalues i = 1/`=colsof(ev)' {
            if (ev[1,`i']>1) local ++nf
        }
        if (`nf'==0) local nf = 1
        if (`nf'>= `=wordcount("`keepitems'")') local nf = `=wordcount("`keepitems'")' - 1

        * If smallest eigenvalue is ~0, warn (near-singular but proceed with FA)
        scalar min_eig = ev[1,`=colsof(ev)']
        if (min_eig < 1e-8) {
            di as error "Warning: correlation matrix near-singular; review items."
        }

        * KMO & Bartlett (postestimation-friendly order)
        factor `keepitems', pcf factors(`nf') corr
        estat kmo
        estat bartlett

        * Rotation: oblique (promax) if >1 factor; varimax if 1 factor
        if (`nf'>1) rotate, promax
        else rotate, varimax

        * Loadings table and export
        estat loadings, blanks(.30)
        matrix L = e(L)
        mat rownames L = `keepitems'
        * Name columns F1..Fk
        local cols : colfullnames L
        local k = colsof(L)
        local newcols
        forvalues j=1/`k' {
            local newcols `newcols' F`j'
        }
        mat colnames L = `newcols'

        if !c(rc) & ("`c(username)'" != "") {
            esttab matrix(L) using "$TABLES\T2_FA_loadings_Q1_Q5.csv", ///
                replace nonumber nomtitles b(3)
        }

        * Factor scores (Bartlett)
        predict F*, bartlett
        forvalues j=1/`nf' {
            label var F`j' "Factor `j' (Bartlett) — Q1–Q5"
        }

        * Reliability (Cronbach's alpha) for each factor bucket (optional):
        * After you inspect loadings, group items per factor and run alpha, e.g.:
        * alpha Q1 Q3, item detail
        * alpha Q2 Q4 Q5, item detail

        * Save scored data (optional)
        tempfile withFA
        save `withFA', replace
        restore
        merge 1:1 mdate using `withFA', nogen
    }
}

*-----------------------------
* 2.3 Event-study plots (face validity)
*-----------------------------
* AEB time series with key verticals (edit dates as needed)
twoway (tsline AEB), ///
    title("AEB over time") ///
    xtitle("") ytitle("Index") ///
    xline(ym(2016,11) 2019m1 2020m3, lpattern(dash)) ///
    note("Verticals: 2016 U.S. election | Trade tensions | COVID onset")
graph export "$FIGURES\F1_aeb_timeseries.png", replace width(1600)

* (Optional) ICC & IFE series
twoway (tsline ICC IFE), ///
    title("ICC & IFE over time") ///
    xtitle("") ytitle("Index") ///
    xline(ym(2016,11) 2019m1 2020m3, lpattern(dash))
graph export "$FIGURES\F2_icc_ife_timeseries.png", replace width(1600)







