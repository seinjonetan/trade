clear

// Set CW
local path = "C:\Users\tsjac\Downloads\"

//////////////////////////////////////////////////////////////
// This first section is just data cleaning. In order to use, you need to have employment, wages, and wage bill datafiles in the same folder. You might want to check the replication files for Galle, Rodrigues-Clare, and Yi (RESTUD). These are available through RESTUD, and may have all data ready to go.
local year = 1980

import delimited "`path'city_occ_e_`year'.csv"

drop other 
tempfile master 
save `master' 

local index = 1
foreach var of varlist occ2_cleric-occ3_shealth {
	di `index'
	qui {
	use `master'
	keep comzone `var' 
	rename `var' L_ck 
	gen k = `index'
	order k, before(L_ck) 
	gen k_desc = "`var'"
	
	if `index' == 1 {
		tempfile temp 
		save `temp'
	}
	else {
		append using `temp'
		sort comzone k
		tempfile temp 
		save `temp'
	}
	local index = `index'+1
	}

}

clear
import delimited "`path'city_occ_w_`year'.csv"

drop other 
tempfile master 
save `master' 

local index = 1
foreach var of varlist occ2_cleric-occ3_shealth {
	di `index'
	qui {
	use `master'
	keep comzone `var' 
	rename `var' w_ck 
	gen k = `index'
	order k, before(w_ck) 
	gen k_desc = "`var'"
	
	if `index' == 1 {
		tempfile temp_w
		save `temp_w'
	}
	else {
		append using `temp_w'
		sort comzone k
		tempfile temp_w
		save `temp_w'
	}
	local index = `index'+1
	}

}

clear
import delimited "`path'city_occ_wb_`year'.csv"

drop other 
tempfile master 
save `master' 

local index = 1
foreach var of varlist occ2_cleric-occ3_shealth {
	di `index'
	qui {
	use `master'
	keep comzone `var' 
	rename `var' y_ck 
	gen k = `index'
	order k, before(y_ck) 
	gen k_desc = "`var'"
	
	if `index' == 1 {
		tempfile temp_y
		save `temp_y'
	}
	else {
		append using `temp_y'
		sort comzone k
		tempfile temp_y
		save `temp_y'
	}
	local index = `index'+1
	}

}

use `temp', clear
merge 1:1 comzone k k_desc using `temp_w'
drop _merge 
merge 1:1 comzone k k_desc using `temp_y'
drop _merge 
sort comzone k 
order k_desc, after(y_ck)


egen c = group(comzone)
order comzone, before(k_desc)
order c, before(k)
gen t = `year'
order t, before(c)
foreach var of varlist L_ck-y_ck {
	replace `var' = 1 if `var' == .
}

save "`path'ipums_clean_1980.dta", replace
//////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////
// Set values for theta and rho. Note that theta>1 and 0<=rho<1.
local theta = 1.5
local rho = 0.5
//////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////
// Generate all necessary variables in order to calculate T_ck values for each city-occupation pair. 
egen L = total(L_ck)
gen pi_ck = L_ck / L 
bysort k: egen L_k = total(L_ck)
gen omega_k = L_k / L 

egen max_pi = max(pi_ck)
gen index = (pi_ck == max_pi)
egen max_omega = total(index*omega_k)

gen w_ck_model = ((pi_ck/max_pi)^((1-`rho')/`theta'))*((omega_k/max_omega)^(`rho'/`theta'))

bysort k: egen y_k = total(y_ck)
local alpha = 7
gen p_ck = (y_ck/y_k)^(1/(1-`alpha'))
gen T_ck = w_ck_model / p_ck 
sort T_ck 

reghdfe T_ck , a(T_c=i.c T_k=i.k) nocons 
gen t_ck = T_ck - T_k - T_c 

sort t_ck 
save "`path'ipums_clean_1980.dta", replace

bysort c: egen pi_c = total(pi_ck)
gen phi_ck = pi_ck / omega_k 
//////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////
// Generate technology elasticity and occupation elasticity.
gen sigma_ck = `theta'*omega_k*((phi_ck/pi_c)-1)
bysort c: egen temp = total(pi_ck*phi_ck)
gen sigma_cc = (`theta'/(1-`rho'))*(1-(pi_c*(1-`rho'))-((`rho'*temp)/pi_c))
drop temp 
//////////////////////////////////////////////////////////////

