clear
set more off

global path "C:\Users\tsjac\Desktop\research\correlated_location\data/"

qui {
// Construct instrument for wages by constructing sector-occupation and city-sector exposure indices.
clear
import delimited "${path}city_sec_employment.csv"
keep if year == 1990
bysort comzone: egen E_c = total(employed)
gen pi_cs = employed / E_c 
keep comzone sector pi_cs 
duplicates drop pi_cs, force 
tempfile pics 
save `pics' 

// Generate HHI for shifts and report
clear
import delimited "${path}city_sec_employment.csv"
keep if year == 1990
egen e_t = total(employed)
bysort sector: egen e_st = total(employed)
duplicates drop sector, force 
keep sector e_t e_st 
gen s_st = e_st / e_t 
gen s2_st = s_st^2
egen hhi_shocks = total(s2_st)
local hhi = hhi_shocks[1]
}
di "HHI Index for shifts: `hhi'"
qui {
clear
import delimited "${path}sec_occ_wb.csv"
drop if indnaics == "Unemployed"
keep if year == 1990
duplicates drop indnaics occupation, force 
bysort occupation: egen wb_k = total(wage_bill)
gen delta_sk = wage_bill / wb_k 
keep indnaics occupation delta_sk 
sort occupation delta_sk 
rename indnaics sector 
tempfile deltask
save `deltask'

clear
import delimited "${path}master.csv"
drop if occupation == "other"
levelsof(occupation), local(levels)
local index = 1
foreach l of local levels {
    di "`l'"
	qui {
	use `deltask', clear
	keep if occupation == "`l'"
	merge 1:m sector using `pics'
	drop if _merge!=3
	drop _merge 
	tempfile temp 
	save `temp'

	foreach t in 1990 2000 2007 2019 {
		use "${path}tfp.dta", clear
		keep if year == `t' 
		merge 1:m sector using `temp' 
		drop if _merge!=3
		drop _merge 
	
		if `t' == 1990 {
			tempfile wiv
			save `wiv'
		}
		if `t' != 1990 {
			append using `wiv'
			sort comzone sector occupation year 
			tempfile wiv 
			save `wiv'
		}
	}
	use `wiv', clear
	
	bysort comzone occupation year: egen w_iv_ckt = total(pi_cs*delta_sk*tfp)
	duplicates drop comzone occupation year, force
	keep comzone year occupation w_iv_ckt 
	sort comzone year occupation 
	
	if `index' == 1 {
	    tempfile iv 
		save `iv'
	}
	if `index' > 1 {
	    append using `iv'
		sort comzone occupation year 
		tempfile iv 
		save `iv'
	}
	local index = `index'+1
	}
}
use `iv', clear 

tempfile iv
save `iv'

/////////////////////////////////////////////////////////////////////
clear
import delimited "${path}master.csv"

drop if occupation == "Other"
keep if year == 1990 | year == 2000 | year == 2007 | year == 2019

merge 1:1 year occupation comzone using `iv'
drop if _merge!=3
drop _merge 
gen temp = 1
bysort comzone: egen check = total(temp)
drop if check != 96 

// Calculate phi instrument as average wage shifter across all other cities within same occupation.
bysort occupation year: egen N = total(temp)
bysort occupation year: egen phi_temp = total(w_iv_ckt)
*gen phi_iv_ckt = w_iv_ckt - (1/(N-1))*(phi_temp-w_iv_ckt)
gen phi_iv_ckt = (phi_temp - w_iv_ckt) / (N-1)
drop temp N check phi_temp 

sort comzone occupation year
save "${path}est_master.dta", replace
}
////////////////////////////////////////////////////////////////
// Estimation 
use "${path}est_master.dta", clear

bysort year: egen e_t = total(employed)
gen pi_ckt = employed / e_t 
bysort year occupation: egen e_kt = total(employed)
gen phi_ckt = employed / e_kt 

gen lnpi = ln(pi_ckt)
gen lnw = ln(wage)
gen lnphi = ln(phi_ckt)

egen t = group(year)
egen c = group(comzone)
egen kt = group(occupation year)
egen ct = group(comzone year)
egen ck = group(comzone occupation)

label variable lnw "Log Wage"
label variable w_iv_ckt "Wage IV"
label variable lnphi "Log phi"
label variable phi_iv_ckt "Phi IV"
label variable lnpi "Log pi"


// First stag estimates
eststo a1: reghdfe lnw w_iv_ckt phi_iv_ckt, a(t) cluster(ck)

eststo a2: reghdfe lnw w_iv_ckt phi_iv_ckt, a(c t) cluster(ck)

eststo a3: reghdfe lnw w_iv_ckt phi_iv_ckt, a(ct) cluster(ck)

eststo a4: reghdfe lnphi w_iv_ckt phi_iv_ckt, a(t) cluster(ck)

eststo a5: reghdfe lnphi w_iv_ckt phi_iv_ckt, a(c t) cluster(ck)

eststo a6: reghdfe lnphi w_iv_ckt phi_iv_ckt, a(ct) cluster(ck)

estadd local tfe "\checkmark", replace: a1 a4
estadd local cplustfe "\checkmark", replace: a2 a5 
estadd local ctfe "\checkmark", replace: a3 a6 

esttab a1 a2 a3 a4 a5 a6 using "${path}wage_fs.tex", /*
*/	replace label b(3) se(3) nonotes keep(w_iv_ckt phi_iv_ckt) star(* 0.10 ** 0.05 *** 0.01) s(tfe cplustfe ctfe r2 N, fmt(%9.2f %9.2f %9.2f %9.2f %9.0f) label("Year FE" "City and Year FE" "City-Year FE" "R$^2$" "Observations")) lines gaps


// Headline estimates
eststo a1: reghdfe lnpi lnw lnphi, a(t) cluster(ck)

eststo a2: ivreghdfe lnpi (lnw lnphi = w_iv_ckt phi_iv_ckt), a(t) cluster(ck)

eststo a3: reghdfe lnpi lnw lnphi, a(c t) cluster(ck)

eststo a4: ivreghdfe lnpi (lnw lnphi = w_iv_ckt phi_iv_ckt), a(c t) cluster(ck)

eststo a5: reghdfe lnpi lnw lnphi, a(ct) cluster(ck)

eststo a6: ivreghdfe lnpi (lnw lnphi = w_iv_ckt phi_iv_ckt), a(ct) cluster(ck)

estadd local tfe "\checkmark", replace: a1 a2
estadd local cplustfe "\checkmark", replace: a3 a4 
estadd local ctfe "\checkmark", replace: a5 a6 
estadd local estimator "2SLS", replace: a2 a4 a6 
estadd local estimator "OLS", replace: a1 a3 a5 
estadd scalar fstat = e(widstat), replace: a2 a4 a6 

esttab a1 a2 a3 a4 a5 a6 using "${path}thetarho_estimates.tex", /*
*/	replace label b(3) se(3) nonotes keep(lnw lnphi) star(* 0.10 ** 0.05 *** 0.01) s(tfe cplustfe ctfe estimator fstat r2 N, fmt(%9.2f %9.2f %9.2f %9.2f %9.2f %9.2f %9.0f) label("Year FE" "City and Year FE" "City-Year FE" "Estimator" "K-P F-Statistic" "R$^2$" "Observations")) lines gaps



// Compare to estimates of theta/(1-rho) recovered using category-time fixed effects.
eststo a1: reghdfe lnpi lnw, a(kt) cluster(ck)

eststo a2: ivreghdfe lnpi (lnw = w_iv_ckt), a(kt) cluster(ck)

qui ivreghdfe lnpi (lnw lnphi = w_iv_ckt phi_iv_ckt), a(t) cluster(ck)
local est = _b[lnw] / (1-_b[lnphi])
estadd scalar baseline = `est', replace: a2 

eststo a3: reghdfe lnpi lnw, a(kt c) cluster(ck)

eststo a4: ivreghdfe lnpi (lnw = w_iv_ckt), a(kt c) cluster(ck)

qui ivreghdfe lnpi (lnw lnphi = w_iv_ckt phi_iv_ckt), a(c t) cluster(ck)
local est = _b[lnw] / (1-_b[lnphi])
estadd scalar baseline = `est', replace: a4 

eststo a5: reghdfe lnpi lnw, a(kt ct) cluster(ck)

eststo a6: ivreghdfe lnpi (lnw = w_iv_ckt), a(kt ct) cluster(ck)

qui ivreghdfe lnpi (lnw lnphi = w_iv_ckt phi_iv_ckt), a(ct) cluster(ck)
local est = _b[lnw] / (1-_b[lnphi])
estadd scalar baseline = `est', replace: a6

estadd local ktfe "\checkmark", replace: a1 a2 a3 a4 a5 a6
estadd local cfe "\checkmark", replace: a3 a4 
estadd local ctfe "\checkmark", replace: a5 a6 
estadd local estimator "2SLS", replace: a2 a4 a6 
estadd local estimator "OLS", replace: a1 a3 a5 
estadd scalar fstat = e(widstat), replace: a2 a4 a6 


esttab a1 a2 a3 a4 a5 a6 using "${path}thetarho_comparison.tex", /*
*/	replace label b(3) se(3) nonotes keep(lnw) star(* 0.10 ** 0.05 *** 0.01) s(ktfe cfe ctfe estimator fstat baseline r2 N, fmt(%9.2f %9.2f %9.2f %9.2f %9.2f %9.2f %9.2f %9.0f) label("Occupation-Year FE" "City FE" "City-Year FE" "Estimator" "K-P F-Statistic" "Baseline" "R$^2$" "Observations")) lines gaps



