version 19
clear all
set more off

*==========================================================
* 01_clean_raw.do  —  Read raw files, clean, save to Processed
*   Uses globals from master.do:
*   $RAW  -> ...\Stata Files\Raw Datasets
*   $PROC -> ...\Stata Files\Processed
*==========================================================

*-----------------------------
* AEB + subindices (monthly) from Excel
*-----------------------------

import excel "$RAW\aebd.xlsx", firstrow clear
list in 1/5

version 19
* 1) Rename the first column (currently "A") to a clear name
rename A ddate                    // daily date
label var ddate "Daily date"

* 2) Create a monthly date from the daily date
gen mdate = mofd(ddate)
format mdate %tm
label var mdate "Month (YYYY–MM)"
order mdate, first

* 3) If you prefer only the monthly date, drop the daily
drop ddate



* 5) (Optional) Drop placeholder columns you don't need (I M Q R S T)
*   (Only run if these exist and you want them gone)
capture drop I M Q R S T

* 6) (Optional) Give friendly labels to key variables
label var AEB "Ag Economy Barometer"
label var ICC "Index of Current Conditions"
label var IFE "Index of Future Expectations"

* 7) Set time series structure on the monthly date
tsset mdate, monthly

save "$PROC\aeb_monthly.dta", replace

*-----------------------------
* (rest of your script continues…)
* Benchmarks (monthly): MCSI, CCI, partisan (R/D), EPU, TradeEPU
*-----------------------------
*MCSI 

	import delimited using "$RAW\UMCSENT.csv", ///
    varnames(1) case(preserve) clear
list in 1/4
*--- UMCSENT: clean to monthly time series

* 1) Make the date column explicit
capture confirm variable observa~e
if _rc==0 rename observa~e date_str
else rename observa* date_str   // fallback if the name truncates differently

* 2) Build a monthly date from the M/D/Y string
gen ddate = daily(date_str, "MDY")
format ddate %td
gen mdate = mofd(ddate)
format mdate %tm
label var mdate "Month (YYYY-MM)"
drop date_str

* 3) Ensure UMCSENT is numeric
capture confirm numeric variable UMCSENT
if _rc {
    destring UMCSENT, replace ignore(".")
}
label var UMCSENT "University of Michigan Consumer Sentiment (UMCSENT)"

* 4) Collapse to one row per month (safe even if already unique)
collapse (mean) UMCSENT, by(mdate)

* --- AFTER: collapse (mean) UMCSENT, by(mdate)
sort mdate

* flag rows that break a clean monthly tail:
* - missing UMCSENT
* - or a gap > 1 month vs previous row
gen bad = missing(UMCSENT) | ( (mdate - mdate[_n-1]) > 1 & _n>1 )

* compute reverse cumulative sum of 'bad' to see if any problems exist from this row to the end
gen __ord = _n
gen __rev = _N - _n + 1
sort __rev
gen __bad_tail = sum(bad)   // cumulate from the end backward
sort __ord
drop __rev __ord

* first row where the tail is clean (i.e., zero bads from here to the end)
quietly summarize mdate if __bad_tail==0, meanonly
local start = r(min)

* print the starting month
display as text "First uninterrupted, complete monthly segment starts at: " %tm `start'

*Removes time period when reporting was quarterly 
drop if __bad_tail > 0

*remove columns for filtering 
drop bad __bad_tail

* 5) Declare the time series structure
tsset mdate, monthly

* 6) Save to Processed (requires $PROC from master.do)
capture noisily save "$PROC\umcsent_monthly.dta", replace 


* CCI (minimal: keep only Period + Value up front)

import excel "$RAW\CCIMo.xlsx", firstrow clear

* --- Keep only the two needed columns right away ---
keep Period Value

* --- Make Value numeric (string -> double if needed) ---
if substr("`: type Value'",1,3)=="str" {
    replace Value = subinstr(Value, char(160), "", .)
    replace Value = strtrim(itrim(Value))
    replace Value = subinstr(Value, ",", "", .)
    replace Value = subinstr(Value, "$", "", .)
    replace Value = subinstr(Value, "%", "", .)
    replace Value = subinstr(Value, "(", "-", .)
    replace Value = subinstr(Value, ")", "", .)
    replace Value = "" if inlist(upper(Value), "NA","N/A","NULL",".","")
    destring Value, replace force
}
recast double Value

* --- Build monthly key from Period ---
capture confirm numeric variable Period
if !_rc {
    local f : format Period
    if strpos("`f'","%td") {
        gen mdate = mofd(Period)
    }
    else {
        gen mdate = mofd(Period + td(1899,12,30))
    }
}
else {
    gen double __dd = daily(Period,"DMY")
    replace __dd = daily(Period,"MDY") if missing(__dd)
    format __dd %td
    gen mdate = mofd(__dd)
    drop __dd
}
format mdate %tm

* --- Finalize monthly series ---
order mdate, first
duplicates drop mdate, force
sort mdate
drop if mdate < tm(1980m1)

* Keep just the monthly key and the value
keep mdate Value
tsset mdate, monthly
compress

* --- Rename + label (no spaces in Stata var names) ---
rename Value CCI_monthly
label var CCI_monthly "CCI monthly"

* Save to Processed
save "$PROC\cci_monthly.dta", replace





*Partisan MCSI
* 1) Import as strings so header rows (1–3) are preserved
import excel "$RAW\redbk05b.xls", allstring clear    // change path/sheet if needed
* e.g., add: sheet("Sheet1")
list in 1/4

* 2) Build new variable names from row 3, grouping party triplets as SI / CI / EI
capture unab vlist : _all
if _rc {
    di as err "No variables imported. Check sheet name/range."
    exit 198
}
local k : word count `vlist'
local partypos = 0   // counts DEM/IND/REP positions across row 3

forvalues j = 1/`k' {
    local v : word `j' of `vlist'
    local cell = `"`=strtrim(`v'[3])'"'    // text in row 3 of this column (DEM/IND/REP/Date...)

    if `j'==1 {
        rename `v' Period
        continue
    }

    * normalize party label
    local up = upper("`cell'")
    if inlist("`up'","DEM","IND","REP") {
        local ++partypos
        local grp = int((`partypos'-1)/3) + 1   // 1,2,3,...
        local party ""
        if "`up'"=="DEM" local party "Dem"
        else if "`up'"=="IND" local party "Ind"
        else if "`up'"=="REP" local party "Rep"

        local metric ""
        if `grp'==1 local metric "SI"
        else if `grp'==2 local metric "CI"
        else if `grp'==3 local metric "EI"
        else local metric "X`grp'"  // just in case there are extra sets

        local target "`metric'_`party'"
    }
    else {
        * if row 3 is blank/noise, give a generic name
        local target "col`j'"
    }

    * ensure uniqueness before renaming
    capture confirm variable `target'
    if !_rc local target "`target'_`j'"

    rename `v' `target'
}

* 3) Drop the first 3 header/spacer rows
drop in 1/3
list in 1/4
 
* --- Build monthly date from Period (month text) + col2 (year) ---
rename col2 Year

* make sure Month & Year are strings we can parse
capture confirm string variable Period
if _rc {  // if Period isn't string, coerce
    tostring Period, replace force
}
tostring Year, replace force

replace Period = strtrim(Period)
replace Year   = strtrim(Year)

* normalize month text (handles "aug", "AUGUST", etc.)
gen str20 _mon = strproper(Period)

* build a daily date for the first of the month, then convert to monthly
gen double _d = date("1 " + _mon + " " + Year, "DMY")
format _d %td

gen mdate = mofd(_d)
format mdate %tm
label var mdate "Month (YYYY-MM)"
order mdate, first

* (optional) clean up
drop _d
drop Period Year // uncomment if you don't need the originals

* (optional) declare monthly time scale like AEB
tsset mdate, monthly

* Convert all variables except mdate to numeric where possible
ds mdate, not
local tofix `r(varlist)'

foreach v of local tofix {
    capture confirm string variable `v'
    if !_rc {
        quietly replace `v' = subinstr(`v', char(160), "", .)
        quietly replace `v' = strtrim(`v')
        quietly replace `v' = subinstr(`v', ",", "", .)
        quietly replace `v' = subinstr(`v', "%", "", .)
        quietly replace `v' = subinstr(`v', "(", "-", .)
        quietly replace `v' = subinstr(`v', ")", "", .)
        capture destring `v', replace force
        * If you want proportions instead of percent points, also do:
        * capture confirm numeric variable `v'
        * if !_rc replace `v' = `v'/100
    }
}

drop col12 
drop col13 col14 col15 col16 col17 col18 col19 _mon 


save "$PROC\partisan_monthly.dta", replace

*USEPU

	version 19
*==============================================================
* US EPU cleaning — USEPU.xlsx
*==============================================================
clear
* If the workbook has multiple sheets, add: sheet("Sheet1")
import excel "$RAW\USEPU.xlsx", firstrow clear

* --- Sanity checks (optional) ---
describe
list in 1/10

* --- Ensure Year / Month exist (handle possible variants) ---
capture confirm variable Year
if _rc {
    capture rename year Year
}
capture confirm variable Month
if _rc {
    capture rename month Month
}

* --- Drop rows with missing in required fields (Year/Month + at least one EPU series) ---
* If you know the exact EPU column name, add it here; otherwise we drop on Year/Month
drop if missing(Year) | missing(Month)

* --- Coerce Year / Month to integers and validate ranges ---
capture confirm numeric variable Year
if _rc destring Year, replace force
capture confirm numeric variable Month
if _rc destring Month, replace force

replace Month = floor(Month)
replace Year  = floor(Year)

* Guard against out-of-range months
drop if Month < 1 | Month > 12

* --- Build YearMonth, Period (daily), and mdate (monthly) ---
tostring Year, replace force
tostring Month, replace force
replace Month = string(real(Month), "%02.0f")

gen str7 YearMonth = Year + "-" + Month
gen double Period = date(YearMonth + "-01", "YMD")
format Period %td

gen mdate = mofd(Period)
format mdate %tm
label var mdate "Month (YYYY-MM)"

* --- Drop helpers & reorder ---
drop Year Month YearMonth
order Period mdate, first

* --- Convert all other columns to numeric where possible ---
ds Period mdate, not
local tofix `r(varlist)'

foreach v of local tofix {
    capture confirm string variable `v'
    if !_rc {
        quietly replace `v' = subinstr(`v', char(160), "", .)
        quietly replace `v' = strtrim(`v')
        quietly replace `v' = subinstr(`v', ",", "", .)
        quietly replace `v' = subinstr(`v', "%", "", .)
        quietly replace `v' = subinstr(`v', "(", "-", .)
        quietly replace `v' = subinstr(`v', ")", "", .)
        capture destring `v', replace force
    }
}

* --- Optional: ensure one obs per month; if duplicates exist, keep last ---
bysort mdate (Period): keep if _n == _N

* --- Optional: quick validations ---
* assert !missing(Period) & !missing(mdate)
* duplicates report mdate
drop Period

* --- Export cleaned outputs ---
save "$PROC\USEPU_clean.dta", replace

display as text "US EPU cleaned: saved USEPU_clean.dta"
	

* GEPU

	version 19
*==============================================================
* Global EPU cleaning — GEPU.xlsx
*==============================================================
clear
* If needed, specify the sheet: sheet("Sheet1")
import excel "$RAW\GEPU.xlsx", firstrow clear

* --- Ensure Year / Month exist (handle lowercase variants) ---
capture confirm variable Year
if _rc capture rename year Year
capture confirm variable Month
if _rc capture rename month Month

* --- Drop rows missing required date fields ---
drop if missing(Year) | missing(Month)

* --- Coerce Year / Month to integers; validate month range ---
capture confirm numeric variable Year
if _rc destring Year, replace force
capture confirm numeric variable Month
if _rc destring Month, replace force
replace Month = floor(Month)
replace Year  = floor(Year)
drop if Month < 1 | Month > 12

* --- Build YearMonth, Period (daily), and mdate (monthly) ---
tostring Year, replace force
tostring Month, replace force
replace Month = string(real(Month), "%02.0f")

gen str7 YearMonth = Year + "-" + Month
gen double Period = date(YearMonth + "-01", "YMD")
format Period %td

gen mdate = mofd(Period)
format mdate %tm
label var mdate "Month (YYYY-MM)"

* --- Drop helpers & reorder ---
drop Year Month YearMonth
order Period mdate, first

* --- Convert all other columns to numeric where possible ---
ds Period mdate, not
local tofix `r(varlist)'
foreach v of local tofix {
    capture confirm string variable `v'
    if !_rc {
        quietly replace `v' = subinstr(`v', char(160), "", .)
        quietly replace `v' = strtrim(`v')
        quietly replace `v' = subinstr(`v', ",", "", .)
        quietly replace `v' = subinstr(`v', "%", "", .)
        quietly replace `v' = subinstr(`v', "(", "-", .)
        quietly replace `v' = subinstr(`v', ")", "", .)
        capture destring `v', replace force
    }
}

* --- Ensure one obs per month; if duplicates exist, keep last by source date ---
bysort mdate (Period): keep if _n == _N
drop Period GEPU_ppp
* --- Optional quick checks ---
* assert !missing(Period) & !missing(mdate)
* duplicates report mdate

* --- Export cleaned outputs ---
export delimited using "$PROC\GEPUClean.csv", replace
save "$PROC\GEPU_clean.dta", replace

display as text "Global EPU cleaned: saved to $PROC\GEPUClean.csv and GEPU_clean.dta"

* TPU
*==============================================================
* TPU cleaning — tpu_web_latest.xlsx (sheet 4)
*==============================================================
clear
import excel "$RAW\tpu_web_latest.xlsx", sheet(TPU_MONTHLY) firstrow clear

* --- Keep only DATE and TPU (handle case variants) ---
capture confirm variable DATE
if _rc {
    capture rename Date DATE
    capture rename date DATE
}
capture confirm variable TPU
if _rc {
    di as err "TPU column not found."
    exit 111
}
keep DATE TPU

* --- Parse DATE to daily %td robustly ---
capture confirm string variable DATE
if !_rc {
    gen double _d = daily(DATE, "YMD")
    replace _d = daily(DATE, "MDY") if missing(_d)
    replace _d = daily(DATE, "DMY") if missing(_d)
    format _d %td
    drop DATE
    rename _d Date
}
else {
    capture confirm numeric variable DATE
    if !_rc {
        local f : format DATE
        if !strpos("`f'", "%td") {
            gen double Date = DATE + td(1899,12,30)
            drop DATE
        }
        else {
            rename DATE Date
        }
        format Date %td
    }
    else {
        di as err "DATE could not be parsed."
        exit 198
    }
}

* --- Build monthly index and tidy order ---
gen mdate = mofd(Date)
format mdate %tm
label var mdate "Month (YYYY-MM)"
order Date mdate TPU

* --- Coerce TPU to numeric if needed ---
capture confirm numeric variable TPU
if _rc {
    quietly replace TPU = subinstr(TPU, char(160), "", .)
    quietly replace TPU = strtrim(TPU)
    quietly replace TPU = subinstr(TPU, ",", "", .)
    quietly replace TPU = subinstr(TPU, "%", "", .)
    quietly replace TPU = subinstr(TPU, "(", "-", .)
    quietly replace TPU = subinstr(TPU, ")", "", .)
    destring TPU, replace force
}

* --- Drop rows missing key fields ---
drop if missing(Date) | missing(TPU)

* --- Ensure one obs per month; keep last if duplicates ---
bysort mdate (Date): keep if _n == _N
drop Date
* --- Optional checks ---
* describe
* list in 1/10
* tsset mdate, monthly

* --- Export cleaned outputs ---
export delimited using "$PROC\TPUClean.csv", replace
save "$PROC\TPU_clean.dta", replace

display as text "TPU cleaned: saved to $PROC\TPUClean.csv and TPU_clean.dta"

*NFIB 
* ==============================================================
* NFIB Uncertainty — "NFIBUncertainty.xlsx"
* ==============================================================

import excel "$RAW\NFIBUncertainty.xlsx", firstrow clear

* Ensure we have a Date variable (common cases)
capture confirm variable Date
if _rc {
    capture rename DATE Date
    capture rename date Date
}

* Parse Date (handles "YYYY-MM-DD" or Excel serial)
capture confirm string variable Date
if !_rc {
    gen double _d = daily(Date, "YMD")
    format _d %td
    drop Date
    rename _d Date
}
else {
    capture confirm numeric variable Date
    if !_rc {
        local f : format Date
        if !strpos("`f'","%td") {
            replace Date = Date + td(1899,12,30)
            format Date %td
        }
    }
}

* Monthly index to match AEB scale
gen mdate = mofd(Date)
format mdate %tm
label var mdate "Month (YYYY-MM)"
order Date mdate, first
drop Date
* Export
export delimited using "$PROC\NFIBUClean.csv", replace
save "$PROC\NFIBU_clean.dta", replace


* ==============================================================
* NFIB Optimism Sentiment — "NFIB Optimism Sentiment.xlsx"
* ==============================================================

clear
import excel "$RAW\NFIB Optimism Sentiment.xlsx", firstrow clear

* Ensure Date exists and parse
capture confirm variable Date
if _rc {
    capture rename DATE Date
    capture rename date Date
}
capture confirm string variable Date
if !_rc {
    gen double _d = daily(Date, "YMD")
    format _d %td
    drop Date
    rename _d Date
}
else {
    capture confirm numeric variable Date
    if !_rc {
        local f2 : format Date
        if !strpos("`f2'","%td") {
            replace Date = Date + td(1899,12,30)
            format Date %td
        }
    }
}

* Monthly index
gen mdate = mofd(Date)
format mdate %tm
label var mdate "Month (YYYY-MM)"
order Date mdate, first

* Convert "Small Business Optimism Index" to numeric and rename SBOI
* (Covers both spaced and auto-renamed variants)
capture confirm variable "Small Business Optimism Index"
if !_rc {
    tempvar s
    gen strL `s' = "Small Business Optimism Index"
    replace `s' = subinstr(`s', char(160), "", .)
    replace `s' = strtrim(`s')
    replace `s' = subinstr(`s', ",", "", .)
    replace `s' = subinstr(`s', "%", "", .)
    replace `s' = subinstr(`s', "(", "-", .)
    replace `s' = subinstr(`s', ")", "", .)
    destring `s', gen(SBOI) force
    drop "Small Business Optimism Index"
}
else {
    * Try common auto-name from import excel
    capture rename Small_Business_Optimism_Index SBOI
    if _rc {
        * If a different name was used, try a wildcard
        capture unab sb : Small*Optimism*Index
        if !_rc {
            local cand : word 1 of `sb'
            capture confirm numeric variable `cand'
            if _rc {
                quietly replace `cand' = subinstr(`cand', char(160), "", .)
                quietly replace `cand' = strtrim(`cand')
                quietly replace `cand' = subinstr(`cand', ",", "", .)
                quietly replace `cand' = subinstr(`cand', "%", "", .)
                quietly replace `cand' = subinstr(`cand', "(", "-", .)
                quietly replace `cand' = subinstr(`cand', ")", "", .)
                destring `cand', replace force
            }
            rename `cand' SBOI
        }
    }
}
label var SBOI "Small Business Optimism Index"
drop Date

* Export
export delimited using "$PROC\NFIBOClean.csv", replace
save "$PROC\NFIBO_clean.dta", replace

display as text "NFIB cleaned: saved to $PROC\NFIBUClean.csv / NFIBU_clean.dta and $PROC\NFIBOClean.csv / NFIBO_clean.dta"

*-----------------------------
* VIX (daily) -> monthly mean
*-----------------------------
import excel "$RAW\VIX.xls", sheet("Chart Data") firstrow clear

* --- Parse PricingDate to daily %td, handling strings or Excel serials ---
capture confirm string variable PricingDate
if !_rc {
    gen double ddate = daily(PricingDate,"YMD")
    replace ddate = daily(PricingDate,"MDY") if missing(ddate)
    replace ddate = daily(PricingDate,"DMY") if missing(ddate)
    format ddate %td
}
else {
    capture confirm numeric variable PricingDate
    if !_rc {
        local f : format PricingDate
        if !strpos("`f'","%td") {
            gen double ddate = PricingDate + td(1899,12,30)   // Excel serial → %td
        }
        else {
            gen double ddate = PricingDate
        }
        format ddate %td
    }
    else {
        di as err "Could not parse PricingDate."
        exit 198
    }
}

* --- Identify the VIX value column and standardize its name to 'vix' ---
local vixvar ""
capture confirm variable vix
if !_rc local vixvar "vix"
if "`vixvar'"=="" {
    capture confirm variable VIX
    if !_rc local vixvar "VIX"
}
if "`vixvar'"=="" {
    capture unab _v : *VIX*
    if !_rc local vixvar : word 1 of `_v'
}
if "`vixvar'"=="" {
    * common alternates
    foreach cand in Close close Price price Value value {
        capture confirm variable `cand'
        if !_rc {
            local vixvar "`cand'"
            continue, break
        }
    }
}
if "`vixvar'"=="" {
    * fallback: first numeric column other than PricingDate
    ds PricingDate, not
    foreach c of local r(varlist) {
        capture confirm numeric variable `c'
        if !_rc {
            local vixvar "`c'"
            continue, break
        }
    }
}
if "`vixvar'"=="" {
    di as err "VIX value column not found."
    exit 111
}
rename `vixvar' vix

* Coerce vix to numeric if needed
capture confirm numeric variable vix
if _rc {
    quietly replace vix = subinstr(vix, char(160), "", .)
    quietly replace vix = strtrim(vix)
    quietly replace vix = subinstr(vix, ",", "", .)
    quietly replace vix = subinstr(vix, "%", "", .)
    quietly replace vix = subinstr(vix, "(", "-", .)
    quietly replace vix = subinstr(vix, ")", "", .)
    destring vix, replace force
}

* --- Daily -> monthly mean ---
gen mdate = mofd(ddate)
format mdate %tm
collapse (mean) vix, by(mdate)
tsset mdate, monthly

save "$PROC\vix_monthly.dta", replace
display as text "VIX monthly mean saved to $PROC\vix_monthly.dta"


*Deere Stock Price 
import excel "$RAW\DEStockPrice.xlsx"
list in 1/4

* Assume current vars are A and B, with row 1 holding header text

* 1) Rename columns using row-1 labels and drop the header row
rename A PricingDate
rename B SharePricing
drop in 1

* 2) Convert PricingDate ("02jan1968") to a proper daily date
capture confirm string variable PricingDate
if !_rc {
    gen double _d = date(PricingDate, "DMY")
    format _d %td
    drop PricingDate
    rename _d PricingDate
}
else {
    * if it's already numeric, just ensure date format
    format PricingDate %td
}

* 3) Ensure SharePricing is numeric double
capture confirm numeric variable SharePricing
if _rc destring SharePricing, replace force
recast double SharePricing
label var SharePricing "NYSE:DE - Share Pricing"
*Cleaning 
order PricingDate, first

* Assume: PricingDate = %td daily, SharePricing = double

* 1) Monthly index
gen mdate = mofd(PricingDate)
format mdate %tm
label var mdate "Month (YYYY-MM)"

* 2) Keep last trading day of each month (EOM)
bysort mdate (PricingDate): keep if _n == _N

* 3) (Optional) restrict to 2015m1+
drop if mdate < tm(2015m1)

* 4) Finalize
rename SharePricing EOM_Close
order mdate PricingDate EOM_Close
tsset mdate, monthly
drop PricingDate

 save "$PROC\DE_monthly.dta", replace
display as text "DE EOM Price saved to $PROC\DE_monthly.dta"

use `"$PROC\aeb_monthly.dta"', clear
keep mdate AEB
format mdate %tm
order mdate, first
duplicates drop mdate, force
sort mdate
compress
save `"$PROC\aeb_monthly_AEBonly.dta"', replace
export delimited using `"$PROC\aeb_monthly_AEBonly.csv"', replace

* Build a slim copy of partisan_monthly with mdate + SI_* only
use `"$PROC\partisan_monthly.dta"', clear
keep mdate SI_Dem SI_Ind SI_Rep   // note: if named differently, adjust here
format mdate %tm
order mdate, first
duplicates drop mdate, force
sort mdate
compress
save `"$PROC\partisan_monthly_SIonly.dta"', replace
export delimited using `"$PROC\partisan_monthly_SIonly.csv"', replace





