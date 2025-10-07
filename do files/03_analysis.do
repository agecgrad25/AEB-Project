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
