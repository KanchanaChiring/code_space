*********************************************
** SAMPLE CODE: CREATING A NORMALISED DEISEL VARIABLE AND RUNNING
** PANEL FE REGRESSIONS OF CROP AREA AND EVAPOTRANSPIRATION ON
** FALLING GROWNDWATER TABLES
*********************************************

clear all
set more off
set rmsg on

* paths 
if "`c(username)'" == "kanchanachiring" {
	global dropbox = "/Users/kanchanachiring/Library/CloudStorage/Dropbox"
	global path = "$dropbox/water"
	/* global overleaf = "xx" */
} 
else if "`c(username)'" == "x" {
	global dropbox = "x"
	global path = "x"
	global overleaf = "x"
}
else {
	display as error "please specify root directory" ///
	"your username is: `c(username)'" ///
	"replace yourname with `c(username)'" 
	exit 
}

* global build path 
global binput "$path/build/input"
global bcode "$path/build/code"
global boutput "$path/build/output"
global btemp "$path/build/temp"
* global analysis path
global acode "$path/analysis/code"
global aoutput "$path/analysis/output"
global atemp "$path/analysis/temp"
global atables "$path/analysis/tables"
* other paths 
global paptab "$overleaf/tables"
global papfig "$overleaf/figures"

global estar3 star(* 0.10 ** 0.05 *** 0.01)
global estar4 star(+ 0.15 * 0.10 ** 0.05 *** 0.01)

*****************************
****

use "$binput/ic_bi_wb_1956_2011.dta", clear

rename dist ic_distcode
tab ic_distcode if ic_distcode == .
drop if ic_distcode == .
drop if year < 1966
drop if year > 2000

tab year if pop_ru !=.

isid ic_distcode year

keep ic_distcode year diesel_ps pop*

tab year if diesel_ps !=.

* Keep only 1991 values to use for imputation
preserve
keep if year == 1991
keep ic_distcode pop_ru
rename pop_ru pop_ru_1991
tempfile base
save `base'

* Restore full dataset and merge
restore
merge m:1 ic_distcode using `base', keep(match master) nogenerate

bysort ic_distcode: egen diesel_ps_max = max(diesel_ps)

gen dps_per_cap_b = diesel_ps_max/pop_ru_1991

// g pop_ru_91 = pop_ru if year == 1991

collapse (mean) dps_per_cap diesel_ps_max pop* , by(ic_distcode)

* bring the dcode61 which will be used to get geoid
merge m:1 ic_distcode ///
using "$btemp/icrisat_wbiac.dta", ///
keepusing(dcode_1961) gen(m_ic)

* flag: 10 districts that were not mapped to shpfiles
gdistinct ic_distcode if m_ic == 1

* drop the districts that were unsucc merged
drop if m_ic ~= 3 

isid dcode_1961

save "$boutput/dcode61_diesel.dta", replace

tempfile usingf_diesel
save `usingf_diesel'

* read outcome data 
use "$boutput/ready_wlcodexyear_for_analysis.dta", clear

xtset objectid year

* merge in source wise irr area 
merge m:1 dcode_1961 using `usingf_diesel', ///
gen(m_dies) 

* flag: 44 districts that were not mapped to wells
gdistinct dcode_1961 if m_dies != 3

drop if m_dies != 3

winsor2 dps_per_*, cuts (5 95)

* make median indicator for diesel pump sets (raw) and per capita
qui summ diesel_ps_max, det
gen diesel_x2 = r(p50)
summ diesel_x2

* diesel_above_p50 is a binary variable that indicates if the district's diesel pump sets is above the median
gen diesel_above_p50 = 1 if diesel_ps_max >= diesel_x2 & !missing(diesel_ps_max)
replace diesel_above_p50 = 0 if diesel_above_p50 != 1

qui summ dps_per_cap_b_w, det
gen dps_per_cap_b_x2 = r(p50)
summ dps_per_cap_b_x2

* dpscapb_above_p50 is a binary variable that indicates if the district's diesel pump sets per capita is above the median
gen dpscapb_above_p50 = 1 if dps_per_cap_b_w >= dps_per_cap_b_x2 & !missing(dps_per_cap_b_w)
replace dpscapb_above_p50 = 0 if dpscapb_above_p50 != 1

merge 1:1 wlcode year ///
using "$boutput/wlcodexyear_ca_cp_et_pet_rf.dta", gen(m_pet) keepusing(pet_win_)

rename pet_win_ pet_win

drop if m_pet !=3

xtset objectid year 

foreach i in ca cp {
	foreach j in 1 3 5 {
	g ln_`i'`j'_win_alt = ln(`i'`j'_win + 0.01)
}
}

foreach i in 1 3 5 {
	g ln_et`i'_win_alt = ln(et`i'_win + 1)
}

   g ln_gw_nov = ln(gw_nov)


***********************************************
** NO. DIESEL PUMPS (RAW NUMBERS)
***********************************************
** No. of diesel pump sets are dist x yr averages so running state x year FE regressions would be informative

*** Effect of GW on CA x Diesel pumps
* gw nov  
foreach j in 1 3 5 {
	* without rainfall control 
	eststo ca`j'_gnxd_r0_m3: reghdfe ca`j'_win gw_nov c.gw_nov##c.diesel_ps_max, absorb(wlcode i.year) vce(cluster tid)
	estadd ysumm
	estadd local rctrl "" 
	eststo ca`j'_gnxd_r0_m4: reghdfe ca`j'_win gw_nov c.gw_nov##c.diesel_ps_max, absorb(wlcode i.sid##i.year) vce(cluster tid) 
	estadd ysumm
	estadd local rctrl ""
   * with rainfall control 
	eststo ca`j'_gnxd_r1_m3: reghdfe ca`j'_win gw_nov rf_win c.gw_nov##c.diesel_ps_max, absorb(wlcode i.year) vce(cluster tid) 
	estadd ysumm
	estadd local rctrl "\checkmark"
	eststo ca`j'_gnxd_r1_m4: reghdfe ca`j'_win gw_nov rf_win c.gw_nov##c.diesel_ps_max, absorb(wlcode i.sid##i.year) vce(cluster tid) 
	estadd ysumm
	estadd local rctrl "\checkmark"
}

foreach j in 1 3 5 {
	esttab ca`j'_gnxd_r0_m* ca`j'_gnxd_r1_m*, keep(c.gw_nov#c.diesel_ps_max) p $estar3 compress ///
	mtitle("yr" "sxy" "yr + r" "sxy + r") ///
	title("winter crop area within `j' in km sq ~ gw + gw x diesel + fe") stats(mean N rctrl, fmt(%9.2gc %9.0gc) labels("Mean" "Obs" "Rainfall Ctrl"))
}


****************************************************************************
**** Below and Above median no. of diesel pumpsets per district across India
****************************************************************************

*** Effect of GW on CA x Diesel pumps above and below median
* gw nov  
foreach j in 1 3 5 {
	* without rainfall control 
	eststo ca`j'_gnxdm_r0_m3: reghdfe ca`j'_win c.gw_nov##i.diesel_above_p50, absorb(wlcode i.year) vce(cluster tid)
	estadd ysumm
	estadd local rctrl ""
	eststo ca`j'_gnxdm_r0_m4: reghdfe ca`j'_win c.gw_nov##i.diesel_above_p50, absorb(wlcode i.sid##i.year) vce(cluster tid) 
	estadd ysumm
	estadd local rctrl ""
}

foreach j in 1 3 5 {
	esttab ca`j'_gnxdm_r0_m*, keep(1.diesel_above_p50#c.gw_nov) p $estar3 compress ///
	mtitle("yr" "sxy") ///
	title("winter crop area within `j' in km sq ~ gw + no rf + gw x diesel above med + fe")
}

* gw nov  
foreach j in 1 3 5 {
	* with rainfall control 
	eststo ca`j'_gnxdm_r1_m3: reghdfe ca`j'_win c.gw_nov##i.diesel_above_p50 rf_win, absorb(wlcode i.year) vce(cluster tid) 
	estadd ysumm
	estadd local rctrl "\checkmark"
	eststo ca`j'_gnxdm_r1_m4: reghdfe ca`j'_win c.gw_nov##i.diesel_above_p50 rf_win, absorb(wlcode i.sid##i.year) vce(cluster tid) 
	estadd ysumm
	estadd local rctrl "\checkmark"
}

foreach j in 1 3 5 {
	esttab ca`j'_gnxdm_r1_m*, keep(1.diesel_above_p50#c.gw_nov) p $estar3 compress ///
	mtitle("yr" "sxy") ///
	title("winter crop area within `j' in km sq ~ gw + rf + gw x diesel above med + fe")
}

********************************************************************
** NO. DIESEL PUMPS PER  CAPITA BASELINE (to normalise across districts)
********************************************************************

** No. of diesel pump sets are dist x yr averages so running state x year FE regressions would be informative

*** Effect of GW x Diesel pumps per capita
* gw nov  
foreach j in 1 3 5 {
	foreach i in ca cp {
	*level-level
	eststo `i'`j'_dpslv_r1_m4: reghdfe `i'`j'_win rf_win c.gw_nov##c.dps_per_cap_b_w, absorb(wlcode i.sid##i.year) vce(cluster tid) 
	estadd ysumm
	estadd local whichm "level-level"
	eststo `i'`j'_dpslv_r1_m5: reghdfe `i'`j'_win rf_win c.gw_nov##c.dps_per_cap_b_w, absorb(wlcode i.did##i.year) vce(cluster tid) 
	estadd ysumm
	estadd local whichm "level-level"
	*log-level
	eststo `i'`j'_dpsllv_r1_m4: reghdfe ln_`i'`j'_win_alt rf_win c.gw_nov##c.dps_per_cap_b_w, absorb(wlcode i.sid##i.year) vce(cluster tid) 
	estadd ysumm
	estadd local whichm "log-level"
	eststo `i'`j'_dpsllv_r1_m5: reghdfe ln_`i'`j'_win_alt rf_win c.gw_nov##c.dps_per_cap_b_w, absorb(wlcode i.did##i.year) vce(cluster tid) 
	estadd ysumm
	estadd local whichm "log-level"
	*log-log
	eststo `i'`j'_dpslnln_r1_m4: reghdfe ln_`i'`j'_win_alt rf_win c.ln_gw_nov##c.dps_per_cap_b_w, absorb(wlcode i.sid##i.year) vce(cluster tid) 
	estadd ysumm
	estadd local whichm "log-log"
	eststo `i'`j'_dpslnln_r1_m5: reghdfe ln_`i'`j'_win_alt rf_win c.ln_gw_nov##c.dps_per_cap_b_w, absorb(wlcode i.did##i.year) vce(cluster tid) 
	estadd ysumm
	estadd local whichm "log-log"
	*level-log
	eststo `i'`j'_dpslvln_r1_m4: reghdfe `i'`j'_win rf_win c.ln_gw_nov##c.dps_per_cap_b_w, absorb(wlcode i.sid##i.year) vce(cluster tid) 
	estadd ysumm
	estadd local whichm "level-log"
	eststo `i'`j'_dpslvln_r1_m5: reghdfe `i'`j'_win rf_win c.ln_gw_nov##c.dps_per_cap_b_w, absorb(wlcode i.did##i.year) vce(cluster tid) 
	estadd ysumm
	estadd local whichm "level-log"
}
}


foreach j in 1 3 5 {
	foreach i in et {
	*level-level
	eststo `i'`j'_dpslv_r1_m4: reghdfe `i'`j'_win rf_win pet_win c.gw_nov##c.dps_per_cap_b_w, absorb(wlcode i.sid##i.year) vce(cluster tid) 
	estadd ysumm
	estadd local whichm "level-level"
	eststo `i'`j'_dpslv_r1_m5: reghdfe `i'`j'_win rf_win pet_win c.gw_nov##c.dps_per_cap_b_w, absorb(wlcode i.did##i.year) vce(cluster tid) 
	estadd ysumm
	estadd local whichm "level-level"
	*log-level
	eststo `i'`j'_dpsllv_r1_m4: reghdfe ln_`i'`j'_win_alt rf_win pet_win c.gw_nov##c.dps_per_cap_b_w, absorb(wlcode i.sid##i.year) vce(cluster tid) 
	estadd ysumm
	estadd local whichm "log-level"
	eststo `i'`j'_dpsllv_r1_m5: reghdfe ln_`i'`j'_win_alt rf_win pet_win c.gw_nov##c.dps_per_cap_b_w, absorb(wlcode i.did##i.year) vce(cluster tid) 
	estadd ysumm
	estadd local whichm "log-level"
	*log-log
	eststo `i'`j'_dpslnln_r1_m4: reghdfe ln_`i'`j'_win_alt rf_win pet_win c.ln_gw_nov##c.dps_per_cap_b_w, absorb(wlcode i.sid##i.year) vce(cluster tid) 
	estadd ysumm
	estadd local whichm "log-log"
	eststo `i'`j'_dpslnln_r1_m5: reghdfe ln_`i'`j'_win_alt rf_win pet_win c.ln_gw_nov##c.dps_per_cap_b_w, absorb(wlcode i.did##i.year) vce(cluster tid) 
	estadd ysumm
	estadd local whichm "log-log"
	*level-log
	eststo `i'`j'_dpslvln_r1_m4: reghdfe `i'`j'_win rf_win pet_win c.ln_gw_nov##c.dps_per_cap_b_w, absorb(wlcode i.sid##i.year) vce(cluster tid) 
	estadd ysumm
	estadd local whichm "level-log"
	eststo `i'`j'_dpslvln_r1_m5: reghdfe `i'`j'_win rf_win pet_win c.ln_gw_nov##c.dps_per_cap_b_w, absorb(wlcode i.did##i.year) vce(cluster tid) 
	estadd ysumm
	estadd local whichm "level-log"
}
}

** state-wise table
foreach i in et cp ca {
* 1 km 
esttab `i'1_dpslv_r1_m4 `i'1_dpsllv_r1_m4 `i'1_dpslnln_r1_m4 `i'1_dpslvln_r1_m4  ///
using "$atables/panel_regs_`i'_dies_st.tex", replace f ///
cells(b(fmt(4)) se(fmt(4) star par)) rename(gw_nov WT ln_gw_nov WT ///
c.gw_nov#c.dps_per_cap_b_w WTxDPS/Capita c.ln_gw_nov#c.dps_per_cap_b_w WTxDPS/Capita) $estar3 ///
keep(WT WTxDPS/Capita) noomitted ///
refcat(WT "\textbf{\emph{Panel A: Area 1 km}}", nolabel) /// 
booktabs nonotes collabels(none) nomtitles ///
stats(N ymean, fmt(%9.2gc %9.2gc) layout("\multicolumn{1}{c}{@}" "\multicolumn{1}{c}{@}") labels("N" "Mean")) 
* 3 km  
esttab `i'3_dpslv_r1_m4 `i'3_dpsllv_r1_m4 `i'3_dpslnln_r1_m4 `i'3_dpslvln_r1_m4 ///
using "$atables/panel_regs_`i'_dies_st.tex", append f ///
cells(b(fmt(4)) se(fmt(4) star par)) rename(gw_nov WT ln_gw_nov WT ///
c.gw_nov#c.dps_per_cap_b_w WTxDPS/Capita c.ln_gw_nov#c.dps_per_cap_b_w WTxDPS/Capita) $estar3 ///
keep(WT WTxDPS/Capita) noomitted ///
refcat(WT "\textbf{\emph{Panel B: Area 3 km}}", nolabel) /// 
booktabs nonotes collabels(none) nomtitles nonumber  ///
stats(N ymean, fmt(%9.2gc %9.2gc) layout("\multicolumn{1}{c}{@}" "\multicolumn{1}{c}{@}") labels("N" "Mean"))
* 5 km 
esttab `i'5_dpslv_r1_m4 `i'5_dpsllv_r1_m4 `i'5_dpslnln_r1_m4 `i'5_dpslvln_r1_m4 ///
using "$atables/panel_regs_`i'_dies_st.tex", append f ///
cells(b(fmt(4)) se(fmt(4) star par)) rename(gw_nov WT ln_gw_nov WT ///
c.gw_nov#c.dps_per_cap_b_w WTxDPS/Capita c.ln_gw_nov#c.dps_per_cap_b_w WTxDPS/Capita) $estar3 ///
keep(WT WTxDPS/Capita) noomitted ///
refcat(WT "\textbf{\emph{Panel B: Area 5 km}}", nolabel) /// 
booktabs nonotes collabels(none) nomtitles nonumber  ///
stats(N ymean whichm, fmt(%9.2gc %9.2gc) layout("\multicolumn{1}{c}{@}" "\multicolumn{1}{c}{@}" "\multicolumn{1}{c}{@}") labels("N" "Mean" "Model"))
}

** dist-wise table
foreach i in et cp ca {
* 1 km 
esttab `i'1_dpslv_r1_m5 `i'1_dpsllv_r1_m5 `i'1_dpslnln_r1_m5 `i'1_dpslvln_r1_m5  ///
using "$atables/panel_regs_`i'_dies_di.tex", replace f ///
cells(b(fmt(4)) se(fmt(4) star par)) rename(gw_nov WT ln_gw_nov WT ///
c.gw_nov#c.dps_per_cap_b_w WTxDPS/Capita c.ln_gw_nov#c.dps_per_cap_b_w WTxDPS/Capita) $estar3 ///
keep(WT WTxDPS/Capita) noomitted ///
refcat(WT "\textbf{\emph{Panel A: Area 1 km}}", nolabel) /// 
booktabs nonotes collabels(none) nomtitles ///
stats(N ymean, fmt(%9.2gc %9.2gc) layout("\multicolumn{1}{c}{@}" "\multicolumn{1}{c}{@}") labels("N" "Mean")) 
* 3 km  
esttab `i'3_dpslv_r1_m5 `i'3_dpsllv_r1_m5 `i'3_dpslnln_r1_m5 `i'3_dpslvln_r1_m5 ///
using "$atables/panel_regs_`i'_dies_di.tex", append f ///
cells(b(fmt(4)) se(fmt(4) star par)) rename(gw_nov WT ln_gw_nov WT ///
c.gw_nov#c.dps_per_cap_b_w WTxDPS/Capita c.ln_gw_nov#c.dps_per_cap_b_w WTxDPS/Capita) $estar3 ///
keep(WT WTxDPS/Capita) noomitted ///
refcat(WT "\textbf{\emph{Panel B: Area 3 km}}", nolabel) /// 
booktabs nonotes collabels(none) nomtitles nonumber  ///
stats(N ymean, fmt(%9.2gc %9.2gc) layout("\multicolumn{1}{c}{@}" "\multicolumn{1}{c}{@}") labels("N" "Mean"))
* 5 km 
esttab `i'5_dpslv_r1_m5 `i'5_dpsllv_r1_m5 `i'5_dpslnln_r1_m5 `i'5_dpslvln_r1_m5 ///
using "$atables/panel_regs_`i'_dies_di.tex", append f ///
cells(b(fmt(4)) se(fmt(4) star par)) rename(gw_nov WT ln_gw_nov WT ///
c.gw_nov#c.dps_per_cap_b_w WTxDPS/Capita c.ln_gw_nov#c.dps_per_cap_b_w WTxDPS/Capita) $estar3 ///
keep(WT WTxDPS/Capita) noomitted ///
refcat(WT "\textbf{\emph{Panel B: Area 5 km}}", nolabel) /// 
booktabs nonotes collabels(none) nomtitles nonumber  ///
stats(N ymean whichm, fmt(%9.2gc %9.2gc) layout("\multicolumn{1}{c}{@}" "\multicolumn{1}{c}{@}" "\multicolumn{1}{c}{@}") labels("N" "Mean" "Model"))
}

************************************************************************************
**** Below and Above median no. of diesel pumpsets per capita per district across India
************************************************************************************
*** Effect of GW on CA x Diesel pumps per capita above and below median
* gw nov  
foreach j in 1 3 5 {
	foreach i in ca cp {
	*level-level
	eststo `i'`j'_mdpslv_r1_m4: reghdfe `i'`j'_win rf_win c.gw_nov##c.dpscapb_above_p50, absorb(wlcode i.sid##i.year) vce(cluster tid) 
	estadd ysumm
	estadd local whichm "level-level"
	eststo `i'`j'_mdpslv_r1_m5: reghdfe `i'`j'_win rf_win c.gw_nov##c.dpscapb_above_p50, absorb(wlcode i.did##i.year) vce(cluster tid) 
	estadd ysumm
	estadd local whichm "level-level"
	*log-level
	eststo `i'`j'_mdpsllv_r1_m4: reghdfe ln_`i'`j'_win_alt rf_win c.gw_nov##c.dpscapb_above_p50, absorb(wlcode i.sid##i.year) vce(cluster tid) 
	estadd ysumm
	estadd local whichm "log-level"
	eststo `i'`j'_mdpsllv_r1_m5: reghdfe ln_`i'`j'_win_alt rf_win c.gw_nov##c.dpscapb_above_p50, absorb(wlcode i.did##i.year) vce(cluster tid) 
	estadd ysumm
	estadd local whichm "log-level"
	*log-log
	eststo `i'`j'_mdpslnln_r1_m4: reghdfe ln_`i'`j'_win_alt rf_win c.ln_gw_nov##c.dpscapb_above_p50, absorb(wlcode i.sid##i.year) vce(cluster tid) 
	estadd ysumm
	estadd local whichm "log-log"
	eststo `i'`j'_mdpslnln_r1_m5: reghdfe ln_`i'`j'_win_alt rf_win c.ln_gw_nov##c.dpscapb_above_p50, absorb(wlcode i.did##i.year) vce(cluster tid) 
	estadd ysumm
	estadd local whichm "log-log"
	*level-log
	eststo `i'`j'_mdpslvln_r1_m4: reghdfe `i'`j'_win rf_win c.ln_gw_nov##c.dpscapb_above_p50, absorb(wlcode i.sid##i.year) vce(cluster tid) 
	estadd ysumm
	estadd local whichm "level-log"
	eststo `i'`j'_mdpslvln_r1_m5: reghdfe `i'`j'_win rf_win c.ln_gw_nov##c.dpscapb_above_p50, absorb(wlcode i.did##i.year) vce(cluster tid) 
	estadd ysumm
	estadd local whichm "level-log"
}
}


foreach j in 1 3 5 {
	foreach i in et {
	*level-level
	eststo `i'`j'_mdpslv_r1_m4: reghdfe `i'`j'_win rf_win pet_win c.gw_nov##c.dpscapb_above_p50, absorb(wlcode i.sid##i.year) vce(cluster tid) 
	estadd ysumm
	estadd local whichm "level-level"
	eststo `i'`j'_mdpslv_r1_m5: reghdfe `i'`j'_win rf_win pet_win c.gw_nov##c.dpscapb_above_p50, absorb(wlcode i.did##i.year) vce(cluster tid) 
	estadd ysumm
	estadd local whichm "level-level"
	*log-level
	eststo `i'`j'_mdpsllv_r1_m4: reghdfe ln_`i'`j'_win_alt rf_win pet_win c.gw_nov##c.dpscapb_above_p50, absorb(wlcode i.sid##i.year) vce(cluster tid) 
	estadd ysumm
	estadd local whichm "log-level"
	eststo `i'`j'_mdpsllv_r1_m5: reghdfe ln_`i'`j'_win_alt rf_win pet_win c.gw_nov##c.dpscapb_above_p50, absorb(wlcode i.did##i.year) vce(cluster tid) 
	estadd ysumm
	estadd local whichm "log-level"
	*log-log
	eststo `i'`j'_mdpslnln_r1_m4: reghdfe ln_`i'`j'_win_alt rf_win pet_win c.ln_gw_nov##c.dpscapb_above_p50, absorb(wlcode i.sid##i.year) vce(cluster tid) 
	estadd ysumm
	estadd local whichm "log-log"
	eststo `i'`j'_mdpslnln_r1_m5: reghdfe ln_`i'`j'_win_alt rf_win pet_win c.ln_gw_nov##c.dpscapb_above_p50, absorb(wlcode i.did##i.year) vce(cluster tid) 
	estadd ysumm
	estadd local whichm "log-log"
	*level-log
	eststo `i'`j'_mdpslvln_r1_m4: reghdfe `i'`j'_win rf_win pet_win c.ln_gw_nov##c.dpscapb_above_p50, absorb(wlcode i.sid##i.year) vce(cluster tid) 
	estadd ysumm
	estadd local whichm "level-log"
	eststo `i'`j'_mdpslvln_r1_m5: reghdfe `i'`j'_win rf_win pet_win c.ln_gw_nov##c.dpscapb_above_p50, absorb(wlcode i.did##i.year) vce(cluster tid) 
	estadd ysumm
	estadd local whichm "level-log"
}
}

** state-wise table
foreach i in et cp ca {
* 1 km 
esttab `i'1_mdpslv_r1_m4 `i'1_mdpsllv_r1_m4 `i'1_mdpslnln_r1_m4 `i'1_mdpslvln_r1_m4  ///
using "$atables/panel_regs_`i'_diesm_st.tex", replace f ///
cells(b(fmt(4)) se(fmt(4) star par)) rename(gw_nov WT ln_gw_nov WT ///
c.gw_nov#c.dpscapb_above_p50 WTxMed_DPS/Capita c.ln_gw_nov#c.dpscapb_above_p50 WTxMed_DPS/Capita) $estar3 ///
keep(WT WTxMed_DPS/Capita) noomitted ///
refcat(WT "\textbf{\emph{Panel A: Area 1 km}}", nolabel) /// 
booktabs nonotes collabels(none) nomtitles ///
stats(N ymean, fmt(%9.2gc %9.2gc) layout("\multicolumn{1}{c}{@}" "\multicolumn{1}{c}{@}") labels("N" "Mean")) 
* 3 km  
esttab `i'3_mdpslv_r1_m4 `i'3_mdpsllv_r1_m4 `i'3_mdpslnln_r1_m4 `i'3_mdpslvln_r1_m4 ///
using "$atables/panel_regs_`i'_diesm_st.tex", append f ///
cells(b(fmt(4)) se(fmt(4) star par)) rename(gw_nov WT ln_gw_nov WT ///
c.gw_nov#c.dpscapb_above_p50 WTxMed_DPS/Capita c.ln_gw_nov#c.dpscapb_above_p50 WTxMed_DPS/Capita) $estar3 ///
keep(WT WTxMed_DPS/Capita) noomitted ///
refcat(WT "\textbf{\emph{Panel B: Area 3 km}}", nolabel) /// 
booktabs nonotes collabels(none) nomtitles nonumber  ///
stats(N ymean, fmt(%9.2gc %9.2gc) layout("\multicolumn{1}{c}{@}" "\multicolumn{1}{c}{@}") labels("N" "Mean"))
* 5 km 
esttab `i'5_mdpslv_r1_m4 `i'5_mdpsllv_r1_m4 `i'5_mdpslnln_r1_m4 `i'5_mdpslvln_r1_m4 ///
using "$atables/panel_regs_`i'_diesm_st.tex", append f ///
cells(b(fmt(4)) se(fmt(4) star par)) rename(gw_nov WT ln_gw_nov WT ///
c.gw_nov#c.dpscapb_above_p50 WTxMed_DPS/Capita c.ln_gw_nov#c.dpscapb_above_p50 WTxMed_DPS/Capita) $estar3 ///
keep(WT WTxMed_DPS/Capita) noomitted ///
refcat(WT "\textbf{\emph{Panel B: Area 5 km}}", nolabel) /// 
booktabs nonotes collabels(none) nomtitles nonumber  ///
stats(N ymean whichm, fmt(%9.2gc %9.2gc) layout("\multicolumn{1}{c}{@}" "\multicolumn{1}{c}{@}" "\multicolumn{1}{c}{@}") labels("N" "Mean" "Model"))
}

** dist-wise table
foreach i in et cp ca {
* 1 km 
esttab `i'1_mdpslv_r1_m5 `i'1_mdpsllv_r1_m5 `i'1_mdpslnln_r1_m5 `i'1_mdpslvln_r1_m5  ///
using "$atables/panel_regs_`i'_diesm_di.tex", replace f ///
cells(b(fmt(4)) se(fmt(4) star par)) rename(gw_nov WT ln_gw_nov WT ///
c.gw_nov#c.dpscapb_above_p50 WTxMed_DPS/Capita c.ln_gw_nov#c.dpscapb_above_p50 WTxMed_DPS/Capita) $estar3 ///
keep(WT WTxMed_DPS/Capita) noomitted ///
refcat(WT "\textbf{\emph{Panel A: Area 1 km}}", nolabel) /// 
booktabs nonotes collabels(none) nomtitles ///
stats(N ymean, fmt(%9.2gc %9.2gc) layout("\multicolumn{1}{c}{@}" "\multicolumn{1}{c}{@}") labels("N" "Mean")) 
* 3 km  
esttab `i'3_mdpslv_r1_m5 `i'3_mdpsllv_r1_m5 `i'3_mdpslnln_r1_m5 `i'3_mdpslvln_r1_m5 ///
using "$atables/panel_regs_`i'_diesm_di.tex", append f ///
cells(b(fmt(4)) se(fmt(4) star par)) rename(gw_nov WT ln_gw_nov WT ///
c.gw_nov#c.dpscapb_above_p50 WTxMed_DPS/Capita c.ln_gw_nov#c.dpscapb_above_p50 WTxMed_DPS/Capita) $estar3 ///
keep(WT WTxMed_DPS/Capita) noomitted ///
refcat(WT "\textbf{\emph{Panel B: Area 3 km}}", nolabel) /// 
booktabs nonotes collabels(none) nomtitles nonumber  ///
stats(N ymean, fmt(%9.2gc %9.2gc) layout("\multicolumn{1}{c}{@}" "\multicolumn{1}{c}{@}") labels("N" "Mean"))
* 5 km 
esttab `i'5_mdpslv_r1_m5 `i'5_mdpsllv_r1_m5 `i'5_mdpslnln_r1_m5 `i'5_mdpslvln_r1_m5 ///
using "$atables/panel_regs_`i'_diesm_di.tex", append f ///
cells(b(fmt(4)) se(fmt(4) star par)) rename(gw_nov WT ln_gw_nov WT ///
c.gw_nov#c.dpscapb_above_p50 WTxMed_DPS/Capita c.ln_gw_nov#c.dpscapb_above_p50 WTxMed_DPS/Capita) $estar3 ///
keep(WT WTxMed_DPS/Capita) noomitted ///
refcat(WT "\textbf{\emph{Panel B: Area 5 km}}", nolabel) /// 
booktabs nonotes collabels(none) nomtitles nonumber  ///
stats(N ymean whichm, fmt(%9.2gc %9.2gc) layout("\multicolumn{1}{c}{@}" "\multicolumn{1}{c}{@}" "\multicolumn{1}{c}{@}") labels("N" "Mean" "Model"))
}
