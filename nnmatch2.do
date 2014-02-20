clear all

cap cd i:/george/comandos_paquetes_librerias/stata/nnmatch2
*! vers 0.13.04 12april2013
*! auth: George Vega
*! cntb: Eduardo Fajnzylber
cap program drop nnmatch2
program def nnmatch2, eclass
	syntax varlist(min=3 numeric) [, METric(string) m(integer 1) Level(cilevel) tc(string)]
	
	vers 9
	
	// Reading variables
	gettoken yi left: varlist
	gettoken w x: left
	
	local nvars : word count `x'
	
	// Selecting metric
	if (length("`metric'") == 0) local metric = "inverse"
	else if (!regexm("`metric'","inverse|lineal|maha")) {
		di as error "Invalid metric `metric'"
		error 1
	}
	
	// Parsing estimator
	if (length("`tc'") == 0) local tc ate
	else if (!regexm("`tc'", "ate|att|atc")) {
		di as error "Invalid treatment effect `tc'"
		error 1
	}
	
	tempvar yl tau complete vartau
	quietly {
		qui gen `yl' = .
	}
	
	preserve
	
	// Keeping complete cases
	quietly {
		egen `complete' = rownonmiss(`yi' `w' `x')
		keep if `complete' == `=2 + `nvars''
	}
	
	marksample touse
	
	#delimit ;
	mata: st_store(., 
		(st_varindex("`yl'")), 
		nn_match(
			st_data(., "`yi'"), 
			st_data(.,"`w'"),
			st_data(.,"`x'"),
			"`metric'",
			`m'
			)
		);
	#delimit cr
	
	
	// Computing Tau
	qui gen `tau' = (`yi' - `yl')*(`w' == 1) + (`yl' - `yi')*(`w' != 1)
	if ("`tc'" == "ate") {         // Average Treatment Effect
		qui reg `tau'
		local tau_est = _b[_cons]
	}
	else if ("`tc'" == "att") {    // Average Treatment Effect for the Treated
		qui reg `tau' if `w' == 1
		local tau_est = _b[_cons]
	}
	else if ("`tc'" == "atc") {    // Average Treatment Effect for the Controls
		qui reg `tau' if `w' != 1
		local tau_est = _b[_cons]
	}
	
	// Computing Variance
	#delimit ;
	mata: st_local("variance", 
		strofreal(nn_variance(st_data(.,"`w'"), st_data(.,"`yi'"), st_data(.,"`yl'"), `tau_est', `m'))
		);
	#delimit cr

	// Printing and storing results
	nnmatch2_results, yi("`yi'") w("`w'") x(`=trim("`x'")') m(`m') tc("`tc'") nobs(`=c(N)') ///
		metric("`metric'") tau(`tau_est') touse(`touse') cmdline("`0'") level(`level') var(`variance')
	
end

// Stores results in e() and prints them as a table
cap program drop nnmatch2_results
program def nnmatch2_results, eclass

	syntax , yi(string) w(string) x(string) m(integer) tc(string) nobs(integer) ///
		metric(string) tau(real) touse(varname) cmdline(string) level(cilevel) var(real)
	
	// Picking estimator name
	if ("`tc'" == "ate") {
		local estname SATE
	}
	if ("`tc'" == "att") {
		local estname SATT
		local lastname "for the Treated"
	}
	if ("`tc'" == "atc") {
		local estname SATC
		local lastname "for the Controls"
	}
	
	tempname b V 
	
	// Returnling Estimators
	mat def `b' = `tau'
	mat colname `b' = `estname'
	mat rowname `b' = y1
	mat def `V' = `var'
	mat colname `V' = `estname'
	mat rowname `V' = `estname'

	ereturn post `b' `V', obs(`nobs') depname("`depn'") esample(`touse')
	
	// Returning scalar
	ereturn scalar m = `m'
	ereturn scalar N = `nobs'
	ereturn scalar se = `=sqrt(`var')'
	local depn "`yi'"
	
	// Returning locals
	ereturn local cmdline "nnmatch2 `cmdline'"
	ereturn local depvar "`depn'"
	ereturn local match_ind "`w'"
	ereturn local match_vars "`x'"
	ereturn local stat "`estname'"
	ereturn local metric "`metric'"
	ereturn local cmd "nnmatch2"

	// Picking weighting
	if `m' > 1 {
		if ("`metric'" == "inverse") {
				loc wmtxt "Weighting matrix: inverse variance"
		}
		else if ("`metric'" == "maha") {
				loc wmtxt "Weighting matrix: Mahalanobis"
		}
	}
	
	// Displaying the results
	di as text _n "Matching estimator: Average Treatment Effect `lastname'" _n
	di as text "`wmtxt' {col 45}Number of obs {col 68}=" as res %10.0f `nobs'
	di as text "`wvtxt' {col 45}Number of matches  (m) = " as res %9.0f `m' _n
	ereturn display, level(`level')
	di as text "{p 0 8 4}Matching variables: `x'{p_end}"
end

// Variance estimator (WORK IN PROGRESS!)
cap mata: mata drop nn_variance()
mata:
real scalar nn_variance(
	real colvector W,
	real colvector Yi,
	real colvector Yl,
	real scalar tau,
	real scalar M
	)
	{
	
	real scalar N
	real colvector output
	
	N = rows(W)
	
	output = sum((W:*(Yi - Yl :- tau) + (1 :- W):*(Yl - Yi :- tau)):^2)*(1/(2*N))
	
	return( ((2/N)*(output)) )
	
}
end

// General purpose matching in covariates.
// nn_match accepts single or multiple covariates as matching in one or more
// neighbours.
cap mata: mata drop nn_match()
mata:
real colvector nn_match(
	real colvector Yi,
	real colvector W,
	real matrix X,
	string scalar metric,
	real scalar M
	)
	{
	
	// Objects declaration
	real colvector Yl, Yl1, Yl0
	real scalar nvars, N, N1, N0, i
	real matrix S, V, Xl, Xl1, Xl0
	
	// Defining constants and ids
	nvars = cols(X)
	N = rows(X)
	Yl  = J(N,1,.)
	
	// Adjusting for sd
	S = quadvariance(X)
		
	// Choosing metric
	if (metric == "maha") { // Mahalanobis matrix
		V =  invsym(S)
	}
	else if (metric == "linear") { // Linear
		V = diag(J(nvars,1,1))
	}
	else { // diagonal matrix of inverse variances
		V = diag(1:/S)
	}
	
	N1 = sum(W)
	N0 = N - N1
	
	// Covariates
	Xl = X,Yi
	_collate(Xl, order(W, 1))
	Xl1 = Xl[|N0 + 1,1\N, nvars|]
	Xl0 = Xl[|1,1\N0, nvars|]
	
	Yl1 = Xl[|N0+1,nvars+1\N,nvars+1|]
	Yl0 = Xl[|1,nvars+1\N0,nvars+1|]

	for(i=1;i<=N;i++) {
		// Starting values for i-th match
		if (!W[|i|]) {
			// Mean match
			Yl[|i|] = mean(bestmatches(calc_distances(Xl1, V, Yl1, X[|i,.|], N1), N1, M))
		}
		else {
			// Mean match
			Yl[|i|] = mean(bestmatches(calc_distances(Xl0, V, Yl0, X[|i,.|], N0), N0, M))
		}
	}
	
	return(Yl)
}
end

// Calculates distances from vector Xi towards each row of matrix Xl using a
// V weighting square matrix (metric matrix such as Inverse matrix or mahalanobis).
// It assumes that Xi vector is control/treated individual and every row in Xl
// is treated/control individual respectively.
cap mata: mata drop calc_distances()
mata:
real matrix calc_distances(
	real matrix Xl,
	real matrix V,
	real colvector Yl,
	real rowvector Xi,
	real scalar Ni
	)
	{
	
	real scalar i
	real matrix xdif, candidates
	
	xdif = Xl - J(Ni,1,Xi)
	
	candidates = (J(Ni,1,.),Yl)
	
	for(i=1;i<=Ni;i++) {
			candidates[|i,1|] = xdif[|i,.|]*V*transposeonly(xdif[|i,.|])
	}
	
	return(candidates)
}
end

// Looks for at least the m best matches, in the case of two (or more) observations
// at the same distance from the current individual, bestmatch will keep both (or
// more) // so, just like nnmatch, nnmatch2's m() option could be more than what
// the user expects.
cap mata: mata drop bestmatches()
mata:
real matrix bestmatches(
	real matrix candidates,
	real scalar nobs,
	real scalar M
	)
	{
	
	// Sorting (descending) over distance
	_collate(candidates, order(candidates, 1))
	
	real scalar M2, i
	
	i  = 0
	M2 = 0
	while ((M2 < M) & (i++ < nobs)) {
		if (candidates[|i,1|] != candidates[|i+1,1|]) M2++
	}

	return(candidates[|1,2\i-1,2|])
}
end

timer clear

cap log close _all
log using results_nnmatch2.txt, text replace

// Lalonde
use ldw_exper, clear
//expand 5

/*
IMPORTANT:
After a couple of runs of nnmatch vs nnmatch2, it turns out that nnmatch gives
biased estimator depending in the number of observations in the dataset. A small
test running nnmatch(2) with original data and expanded (which should not imply
any change in the estimators) turns out in different estimators for nnmatch and
unchanged estimators for nnmatch2 (which I think should be right).
(SEE THE EXAMPLE BELOW)
*/

// Onematch
timer on 1
nnmatch2 re78 t age educ black hisp married re74 re75 u74 u75 , m(1)
timer off 1

timer on 2
nnmatch  re78 t age educ black hisp married re74 re75 u74 u75 , m(1)
timer off 2

qui timer list
di "{hline}" _n "N of speedups (with con `=c(N)' obs) " %5.2fc `=r(t2)' "/" %5.2fc `=r(t1)' " secs. = x" %5.2fc `=r(t2)/r(t1)' _n "{hline}"
timer clear

// Four matches
timer on 1
nnmatch2 re78 t age educ black hisp married re74 re75 u74 u75 , m(4)
timer off 1

timer on 2
nnmatch  re78 t age educ black hisp married re74 re75 u74 u75 , m(4)
timer off 2

qui timer list
di "{hline}" _n "N of speedups (with con `=c(N)' obs) " %5.2fc `=r(t2)' "/" %5.2fc `=r(t1)' " secs. = x" %5.2fc `=r(t2)/r(t1)' _n "{hline}"
timer clear

////////////////////////////////////////////////////////////////////////////////
// Expanding data (estimators should keep unchanged)
expand 10

// One matche
timer on 1
nnmatch2 re78 t age educ black hisp married re74 re75 u74 u75 , m(1)
timer off 1

timer on 2
nnmatch  re78 t age educ black hisp married re74 re75 u74 u75 , m(1)
timer off 2

qui timer list
di "{hline}" _n "N of speedups (with con `=c(N)' obs) " %5.2fc `=r(t2)' "/" %5.2fc `=r(t1)' " secs. = x" %5.2fc `=r(t2)/r(t1)' _n "{hline}"
timer clear

// Four matches
timer on 1
nnmatch2 re78 t age educ black hisp married re74 re75 u74 u75 , m(4)
timer off 1

timer on 2
nnmatch  re78 t age educ black hisp married re74 re75 u74 u75 , m(4)
timer off 2

qui timer list
di "{hline}" _n "N of speedups (with con `=c(N)' obs) " %5.2fc `=r(t2)' "/" %5.2fc `=r(t1)' " secs. = x" %5.2fc `=r(t2)/r(t1)' _n "{hline}"	

log close _all

// exit, clear

