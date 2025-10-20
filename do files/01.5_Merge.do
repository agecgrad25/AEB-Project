**************************************************************
* build_all_into_one_simple.do — Merge all series to one table
* Base timeline: USEPU_clean.dta (monthly mdate)
* Requires: Stata 19, $PROC set
**************************************************************
version 19.0
clear all
set more off

* -------------------------------
* 1) Load base (USEPU timeline)
* -------------------------------
use "$PROC\USEPU_clean.dta", clear

* Build/normalize mdate if needed
capture confirm variable mdate
if _rc {
    capture confirm variable Date
    if !_rc gen mdate = mofd(Date)
    else {
        capture confirm variable Year
        capture confirm variable Month
        if !_rc & !_rc gen mdate = ym(Year, Month)
        else {
            di as err "Cannot find or construct mdate in USEPU_clean.dta"
            exit 459
        }
    }
}
format mdate %tm
duplicates drop mdate, force
sort mdate
order mdate, first
label data "firstcorrsv2 — master timeline from USEPU (mdate)"

* -------------------------------------------------------
* 2) Helper: merge on mdate and suffix non-key variables
* -------------------------------------------------------
program define _merge_on_mdate, rclass
    version 19.0
    syntax , Path(string) Suffix(string)

    preserve
        quietly use "`path'", clear

        * Build/normalize mdate in using data
        capture confirm variable mdate
        if _rc {
            capture confirm variable Date
            if !_rc gen mdate = mofd(Date)
            else {
                capture confirm variable Year
                capture confirm variable Month
                if !_rc & !_rc gen mdate = ym(Year, Month)
                else {
                    di as err "Skipping `path' — cannot construct mdate."
                    restore
                    exit
                }
            }
        }
        format mdate %tm
        duplicates drop mdate, force
        sort mdate
        order mdate, first

        * Suffix all non-key vars to avoid name collisions
        ds mdate, not
        local vars `r(varlist)'
        if "`vars'" != "" {
            foreach v of local vars {
                capture confirm variable `v'
                if !_rc rename `v' `v'_`suffix'
            }
        }

        tempfile _u
        quietly save "`_u'"
    restore

    merge 1:1 mdate using "`_u'", keep(master match) nogen
end

* ------------------------------------
* 3) Merge all sources (one by one)
* ------------------------------------
* Assumes _merge_on_mdate is already defined above.

_merge_on_mdate, path("$PROC\aeb_monthly_AEBonly.dta")     suffix(aeb)
_merge_on_mdate, path("$PROC\partisan_monthly_SIonly.dta") suffix(part)
_merge_on_mdate, path("$PROC\cci_monthly.dta")             suffix(cci)
_merge_on_mdate, path("$PROC\umcsent_monthly.dta")         suffix(umc)
_merge_on_mdate, path("$PROC\DE_monthly.dta")              suffix(de)
_merge_on_mdate, path("$PROC\GEPU_clean.dta")              suffix(gepu)
_merge_on_mdate, path("$PROC\TPU_clean.dta")               suffix(tpu)
_merge_on_mdate, path("$PROC\NFIBO_clean.dta")             suffix(nfibo)
_merge_on_mdate, path("$PROC\NFIBU_clean.dta")             suffix(nfibu)
_merge_on_mdate, path("$PROC\vix_monthly.dta")             suffix(vix)

compress
save "$PROC\firstcorrsv2.dta", replace
export delimited using "$PROC\firstcorrsv2.csv", replace
di as result "✅ Saved: $PROC\firstcorrsv2.dta (+ CSV)."

* --------------------------
* 4) Post-2014 filtered copy
* --------------------------
preserve
    keep if mdate >= tm(2014m12)
    compress
    save "$PROC\firstcorrsv2_post2014.dta", replace
    export delimited using "$PROC\firstcorrsv2_post2014.csv", replace
restore
di as result "✅ Saved post-2014 copy."

* -----------------------------------------------
* 5) Add corn/soy vols to the post-2014 dataset
* -----------------------------------------------
use "$PROC\firstcorrsv2_post2014.dta", clear
format mdate %tm

merge 1:1 mdate using "$PROC\cornvol.dta", ///
    keep(master match) nogen keepusing(vol_cc_month vol_cc_ann)
rename vol_cc_month corn_vol_month
rename vol_cc_ann   corn_vol_ann

merge 1:1 mdate using "$PROC\sbvol.dta", ///
    keep(master match) nogen keepusing(vol_cc_month vol_cc_ann)
rename vol_cc_month sb_vol_month
rename vol_cc_ann   sb_vol_ann

order corn_vol_month corn_vol_ann sb_vol_month sb_vol_ann, after(mdate)
compress
save "$PROC\aebcorrsv3.dta", replace
di as result "✅ Final: $PROC\aebcorrsv3.dta (post-2014 + corn/soy vols)."

* -------------------------------------------------
* 6) Add simple seasonal dummies (Spring/Summer/Fall/Winter)
*     Meteorological seasons (first full month starts the season):
*     Spring = Mar–May; Summer = Jun–Aug; Fall = Sep–Nov; Winter = Dec–Feb
* -------------------------------------------------
use "$PROC\aebcorrsv3.dta", clear
format mdate %tm
sort mdate

gen byte se_spring = inlist(month(dofm(mdate)), 3, 4, 5)
label var se_spring "Season dummy: Spring (Mar–May)"

gen byte se_summer = inlist(month(dofm(mdate)), 6, 7, 8)
label var se_summer "Season dummy: Summer (Jun–Aug)"

gen byte se_fall   = inlist(month(dofm(mdate)), 9, 10, 11)
label var se_fall  "Season dummy: Fall (Sep–Nov)"

gen byte se_winter = inlist(month(dofm(mdate)), 12, 1, 2)
label var se_winter "Season dummy: Winter (Dec–Feb)"

order se_spring se_summer se_fall se_winter, after(mdate)
compress
save "$PROC\aebcorrsv3.dta", replace
export delimited using "$PROC\aebcorrsv3.csv", replace

* ------------------------------------
* 7) First-difference table (Δx_t)
*     Include seasonal dummies in output, but do NOT difference them
* ------------------------------------
use "$PROC\aebcorrsv3.dta", clear
format mdate %tm
tsset mdate, monthly

* All numeric vars except the key
ds, has(type numeric)
local nums `r(varlist)'
local nums : list nums - mdate

* Identify seasonal dummies
ds se_*, has(type numeric)
local seasons `r(varlist)'

* Vars to difference = all numeric minus seasonal dummies
local diffvars : list nums - seasons

foreach v of local diffvars {
    capture gen double d_`v' = D.`v'
    if !_rc label var d_`v' "First difference of `v'"
}

* Keep mdate, seasonal dummies, and all diffs
local keepvars "mdate d_*"
if "`seasons'" != "" local keepvars "mdate `seasons' d_*"
keep `keepvars'
order mdate `seasons' d_*, first

compress
save "$PROC\aebcorrsv3_diff.dta", replace
export delimited using "$PROC\aebcorrsv3_diff.csv", replace



