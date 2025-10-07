**************************************************************
* build_all_into_one.do  — Put ALL series onto one table
* Base timeline: USEPU_clean.dta  (monthly mdate)
* Requires: Stata 16+ (frames), $PROC (and $RAW if you use it) set
**************************************************************
clear all
set more off

*---------------------------------------
* 0) Quick sanity (optional)
*---------------------------------------
di as txt "PROC = " `"$PROC"'

*---------------------------------------
* 1) Base frame = USEPU timeline
*---------------------------------------
* Replace: frame drop _all
capture frame drop base
frame create base
frame change base
clear

use `"$PROC\USEPU_clean.dta"', clear

* Ensure monthly key exists
capture confirm variable mdate
if _rc {
    capture confirm variable Date
    if !_rc gen mdate = mofd(Date)
}
capture confirm variable mdate
if _rc {
    di as err "USEPU_clean.dta has no mdate and no Date to build it."
    exit 459
}
format mdate %tm
duplicates drop mdate, force
sort mdate
compress
label data "firstcorrsv2 — master timeline from USEPU (mdate)"

*------------------------------------------------------------
* Safe fetcher: link by mdate, frget vars, then rename in base
*------------------------------------------------------------
program drop _all
program define _bring_into_base3
    version 16
    syntax , Path(string) Suffix(string)

    * source frame
    frame create src_`suffix'
    frame src_`suffix': use "`path'", clear

    * ensure mdate
    frame src_`suffix': capture confirm variable mdate
    if _rc {
        frame src_`suffix': capture confirm variable Date
        if !_rc frame src_`suffix': gen mdate = mofd(Date)
        else {
            frame src_`suffix': capture confirm variable Year
            frame src_`suffix': capture confirm variable Month
            if !_rc & !_rc frame src_`suffix': gen mdate = ym(Year,Month)
        }
    }
    frame src_`suffix': capture confirm variable mdate
    if _rc {
        di as err "Skipping `path' — cannot construct mdate."
        frame drop src_`suffix'
        exit
    }

    * normalize & dedup
    frame src_`suffix': format mdate %tm
    frame src_`suffix': duplicates drop mdate, force
    frame src_`suffix': sort mdate

    * vars to bring (exclude key)
    frame src_`suffix': ds mdate, not
    local fetch `r(varlist)'
    if "`fetch'"=="" {
        frame drop src_`suffix'
        exit
    }

    * link and fetch into base
    frame change base
    frlink 1:1 mdate, frame(src_`suffix')
    frget `fetch', from(src_`suffix')
    capture drop src_`suffix'   // drop link var

    * rename fetched vars with suffix (avoid collisions later)
    foreach v of local fetch {
        capture confirm variable `v'
        if !_rc {
            local new = strtoname("`v'") + "_`suffix'"
            capture confirm variable `new'
            if _rc rename `v' `new'
            else {
                * if new already exists (unlikely now), append a numeric tail
                local k = 1
                while !`_rc' {
                    local new2 = strtoname("`v'") + "_`suffix'_" + string(`k')
                    capture confirm variable `new2'
                    if _rc {
                        rename `v' `new2'
                        continue, break
                    }
                    local ++k
                }
            }
        }
    }

    * cleanup source frame
    frame drop src_`suffix'
end



_bring_into_base3, path(`"$PROC\aeb_monthly_AEBonly.dta"')      suffix(aeb)
_bring_into_base3, path(`"$PROC\partisan_monthly_SIonly.dta"')  suffix(part)
_bring_into_base3, path(`"$PROC\cci_monthly.dta"')              suffix(cci)
_bring_into_base3, path(`"$PROC\umcsent_monthly.dta"')          suffix(umc)
_bring_into_base3, path(`"$PROC\DE_monthly.dta"')               suffix(de)
_bring_into_base3, path(`"$PROC\GEPU_clean.dta"')               suffix(gepu)
_bring_into_base3, path(`"$PROC\TPU_clean.dta"')                suffix(tpu)
_bring_into_base3, path(`"$PROC\NFIBO_clean.dta"')              suffix(nfibo)
_bring_into_base3, path(`"$PROC\NFIBU_clean.dta"')              suffix(nfibu)
_bring_into_base3, path(`"$PROC\vix_monthly.dta"')              suffix(vix)

order mdate, first
sort mdate
compress
save `"$PROC\firstcorrsv2.dta"', replace
export delimited using `"$PROC\firstcorrsv2.csv"', replace




*---------------------------------------
* 4) Finalize single wide table
*---------------------------------------
frame change base
order mdate, first
sort mdate
compress
label data "firstcorrsv2 — ALL monthly series on USEPU timeline"

* Save the "one new sheet"
save `"$PROC\firstcorrsv2.dta"', replace
export delimited using `"$PROC\firstcorrsv2.csv"', replace

di as result "✅ Built firstcorrsv2 with all series merged into one table."



* Make a post-2014 copy of the merged sheet
use `"$PROC\firstcorrsv2.dta"', clear

* Ensure monthly time index is set/consistent
capture confirm variable mdate
if _rc {
    di as err "mdate missing in firstcorrsv2.dta"
    exit 459
}
format mdate %tm

* Save a duplicate under a new name
save `"$PROC\firstcorrsv2_post2014.dta"', replace

* Reopen the duplicate and filter
use `"$PROC\firstcorrsv2_post2014.dta"', clear
drop if mdate < tm(2014m12)

* Housekeeping + exports
order mdate, first
sort mdate
compress
save `"$PROC\firstcorrsv2_post2014.dta"', replace
export delimited using `"$PROC\firstcorrsv2_post2014.csv"', replace

di as result "✅ Kept observations from 2014m12+ and saved to:"
di as txt     "  - $PROC\firstcorrsv2_post2014.dta"
di as txt     "  - $PROC\firstcorrsv2_post2014.csv"

