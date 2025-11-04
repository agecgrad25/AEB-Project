version 19
clear all
set more off
set rmsg on

*Run on Work PC
global PROJ "C:\Users\adolbowv\OneDrive - purdue.edu\Work Projects\Barometer correlations\Stata Files"

*Run on own pc
*global PROJ "C:\Users\acers\OneDrive - purdue.edu\Work Projects\Barometer correlations\Stata Files"

global CODE "$PROJ\do files"
global RAW  "$PROJ\Raw Datasets"
global PROC "$PROJ\Processed"
global FIG  "$PROJ\Figures"
global TAB  "$PROJ\Tables"
global DOC  "$PROJ\Docs"

cd "$PROJ"
cap mkdir "$PROC"
cap mkdir "$FIG"
cap mkdir "$TAB"
cap mkdir "$DOC"

log using "$DOC\runlog.smcl", replace

do "$CODE\00_setup_install.do"
do "$CODE\01_clean_raw.do"
do "$CODE\01.5_Merge.do"
do "$CODE\02_build_panels.do"
do "$CODE\03_analysis.do"
do "$CODE\04_export_outputs.do"

log close
display as text "Done!"

