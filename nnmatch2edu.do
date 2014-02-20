cap program drop nnmatch2edu
program def nnmatch2edu
	syntax varlist(min=3 numeric)
	
	tokenize `varlist'
	local y `1'
	local t `2'
	local x `3'
	
	// Generating ids (for future merging)
	tempvar id
	gen `id' = _n
	
	tempvar matchid
	qui gen `matchid' = .
	
	tempvar matchy
	qui gen `matchy' = .
	
	tempvar vari
	qui gen `vari' = .
	
	mata: st_store(., (st_varindex("`matchy'"),st_varindex("`matchid'"),st_varindex("`vari'")), bestmatch(st_data(., "`y'"), st_data(.,"`t'"),st_data(.,"`x'")))
	
	tempvar att
	gen `att' = (`y' - `matchy')*(`t' == 1) + (`matchy' - `y')*(`t' != 1)
	
	summ `att', meanonly
	di "{hline}" _n as text "ATE Estimator (not adjusted) = " as result r(mean) _n "{hline}"
	summ `vari', meanonly
	local se = r(mean)^0.5
	di "{hline}" _n as text "Standard error = " as result `se' _n "{hline}"

	
	/*
	// IF YOU WANT TO SAVE THE RESULTS
	cap drop att
	gen att = `att'
	
	cap drop matchid
	gen matchid = `matchid'
	
	cap drop matchy
	gen matchy = `matchy'*/
end

cap mata: mata drop bestmatch()
mata:
real matrix bestmatch(
	real colvector y,
	real colvector t,
	real matrix x
	)
	{
	
	// Objects declaration
	// vari contains the individual level variance
	// ki corresponds to the number of units each individual has been matched to
	// tau is the SATE
	// sigma2 is the individual variance under constant treatment and heterocedasticity
	real colvector ids, matchid, matchy, vari, ki
	real scalar nvars, nobs
	real scalar i, i2, tau, sigma2
	real scalar currentx, currentt, mindif
	
	// Defining constants and ids
	nvars = cols(x)
	nobs = rows(x)
	ids = 1::nobs
	matchid = J(nobs,1,.)
	matchy  = J(nobs,1,.)
	ki = J(nobs,1,0)
	
	// Adjusting for sd
	x = x :/diagonal(sqrt(variance(x)))'
	
	// First obsloop
	for(i=1;i<=nobs;i++) {
		
		// Starting values for i-th match
		currentx = x[i,1]
		currentt = t[i]
		mindif   = .
		
		// Second obsloop
		for(i2=1;i2<=nobs;i2++) {			
			// Matches only if ith != i2th
			if (i != i2 & currentt != t[i2]) {
				if (abs(currentx - x[i2,1]) < mindif) {
					mindif = abs(x[i,1] - x[i2,1])
					matchid[i] = ids[i2]
					matchy[i] = y[i2]
				}
			}
		}
	}
	// calculating ki (valid for only 1 match)
	for(i=1;i<=nobs;i++) {
		ki[i] = sum(matchid:==i)
	}
	//(t,ki)
	sum(ki)
	tau = 	1/nobs*sum(((2*t:-1):*(1:+ki)):*y)
	sigma2 = 1/(2*nobs)*(t:*(y-matchy:-tau)+(1:-t):*(matchy-y:-tau))'*(t:*(y-matchy:-tau)+(1:-t):*(matchy-y:-tau))
	vari = 1/nobs*((1:+ki):*(1:+ki))*sigma2
	//rows((matchy,matchid,vari))
	//cols((matchy,matchid,vari))
	//(matchy,matchid,vari)
	return((matchy,matchid,vari))	
	//return((matchy,matchid))	
}
end

sysuse auto, clear
expand 20

timer clear

gen idold = _n
qui timer on 1
nnmatch2edu price foreign weight
qui timer off 1
/*
qui timer on 2
nnmatch price foreign weight
qui timer off 2
*/
timer list
di 
di "{hline}" _n "N veces aumenta velocidad (con `=c(N)' obs) = x" `=r(t2)/r(t1)' _n "{hline}"
//sort price

