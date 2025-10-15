*04_export_outputs

use "$PROC\aebcorrsv3.dta", clear
format mdate %tm
tsset mdate, monthly

* Fit AEB on Fourier terms only to isolate deterministic seasonality
capture drop AEB_season_fit
quietly regress AEB F_s* F_c*
predict double AEB_season_fit, xb

twoway (tsline AEB) (tsline AEB_season_fit), ///
    title("AEB vs Fourier seasonality fit") ///
    legend(order(1 "AEB" 2 "Fourier-only fit") pos(3) ring(0)) ///
    ytitle("") xtitle("")
graph export "$FIG\G_AEB_fourier_fit_v3.png", replace width(1600)
