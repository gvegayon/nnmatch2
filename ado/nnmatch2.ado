*! vers 0.13.04 18apr2013
*! auth: George Vega
*! auth: Eduardo Fajnzylber
*! A new version of nearest neighbour matching estimator in Stata
*! ON DEVELOPMENT VERSION
cap program drop nnmatch2
program def nnmatch2, eclass
	if !replay() { // If it is a replay, 
		nnmatch2main `0'
	}
	else if "`e(cmd)'" != "nnmatch2" { // If the last command runned wasn't nnmatch2
		di as error "Last estimation was not {bf:nnmatch2}, was {bf:`=e(cmd)'}" 
		exit 301
	}
	else { // if the last command was nnmatch2, it replays the results
		tempname b V
		mat def `b' = e(b)
		mat def `V' = e(V)
		
		tempvar touse
		gen `touse' = e(sample)
		
		// Printing and storing results
		nnmatch2_results, yi(`=e(depvar)') w(`=e(match_ind)') x(`=e(match_vars)') /// 
			m(`=e(m)') tc(`=lower(regexr(e(stat),"^[a-zA-Z]",""))') nobs(`=e(N)') ///
			metric(`=e(metric)') tau(`=`b'[1,1]') cmdline(`=e(cmdline)') level(`level') ///
			var(`=`V'[1,1]') sigma2(`=e(sigma2)') bias(`=e(bias)') touse(`touse')
	}
end

cap program drop nnmatch2main
program def nnmatch2main, eclass sortpreserve
	syntax varlist(min=3 numeric) [in] [if] [, METric(string) m(integer 1) Level(cilevel) tc(string) ///
		BIASadj(string) save(string)]
	
	vers 9
	
	// Reading variables
	gettoken yi left: varlist
	gettoken w x: left
	
	local nvars : word count `x'
	
	// Selecting metric
	if (length("`metric'") == 0) local metric = "inverse"
	else if (!regexm("`metric'","inverse|lineal|maha")) {
		cap confirm mat `metric'
		if (_rc != 0) {
			di as error "matrix `metric' not found"
		}
	}
	
	// Parsing estimator
	if (length("`tc'") == 0) local tc ate
	else if (!regexm("`tc'", "ate|att|atc")) {
		di as error "Invalid treatment effect `tc'"
		error 1
	}
	
	// Parsing bias-adjustment
	if (length("`biasadj'") > 0) {
		if ("`biasadj'" == "bias") local biasvars `x'
		else                       local biasvars `biasadj'
		
		foreach var in `biasvars' {
			if (!regexm("`x'", "(^|[\t ])`var'($|[\t ])")) {
				di as error "{bf:`var'} of biasadj varlist is not contained in {bf:x}"
				di as error "Bias-adj is not supported for other variables rather than those in {bf:x}"
				error 1
			}
		}
	}
	
	// Marking sample
	marksample touse
	markout `touse' `yi' `w' `x' // Marking complete cases
	
	
	// Counting the number of observations to be matched
	cap count if `touse'
	local ntouse = r(N)
	gsort -`touse'
	
	// Saving place for tempvars
	tempvar yl tau nmi yl2 kmi vari vari2 idi idl kl

	
	quietly {
		gen `idi' = _n if `touse'
		gen `yl' = . if `touse'  // Matched Y
		gen `nmi' = . if `touse' // Number of matched to i
		gen `yl2' = . if `touse' // Matched Y^2
		gen `kmi' = . if `touse' // Number of times i was matched to another obs
		
		// Average covariates of matched individuals
		foreach var of var `x' {
			tempvar `var'_l
			gen ``var'_l' = .
			local xls `xls' ``var'_l'
		}
	}
	
	// Matching algorithm	
	tempname matches	
	#delimit ;
	mata: st_store((1, `ntouse'), 
		(st_varindex(tokens("`yl' `nmi' `yl2' `kmi' `xls'"))), 
		nn_match(
			st_data((1,`ntouse'), "`yi'"), 
			st_data((1,`ntouse'),"`w'"),
			st_data((1,`ntouse'),"`x'"),
			"`metric'",
			`m',
			&(`matches' = J(1,5+`nvars'*2,.)) // Pointer to matches matrix (used for computing biasadj)
			)
		);
	#delimit cr 
	
	// If bias estimator
	if (length("`biasadj'") != 0) {
	
		tempname betas0 betas1			
		qui reg `yi' `biasvars' [aw=`kmi'] if `touse' & `w'
		mat def `betas1' = e(b)
		
		qui reg `yi' `biasvars' [aw=`kmi'] if `touse' & !`w'
		mat def `betas0' = e(b)
		
		// Computing mus predictions
		tempvar mu_0i mu_0l mu_1i mu_1l
		
		mat score `mu_1i' = `betas0'
		mat score `mu_0i' = `betas1'
		
		// Building list of Xls var names to be included in the bias-adjustment
		foreach var in `biasvars' {
			foreach var2 in `x' {
				if ("`var'"=="`var2'") {
					local xlsbias `xlsbias' ``var2'_l'
				}
			}
		}
		
		// Replacing betas names Xis with Xls
		mat colname `betas0' = `xlsbias' _cons
		mat colname `betas1' = `xlsbias' _cons
		
		mat score `mu_1l' = `betas0'
		mat score `mu_0l' = `betas1'
		
		// Getting back to the original names
		mat colname `betas0' = `biasvars' _cons
		mat colname `betas1' = `biasvars' _cons
		
		// Generating yl bias corrected of the form yl + u(xi) - u(xl)
		qui {
			replace `yl' = `yl' + `w'*(`mu_1i'-`mu_1l') + (1-`w')*(`mu_0i'-`mu_0l')
		}
	}
	if (length("`save'") > 0) save `save', replace
		
	
	// Computing Tau (Average Treatment Estimator)
	qui gen `tau' = (`yi' - `yl')*(`w' == 1) + (`yl' - `yi')*(`w' != 1) if `touse'
	tempvar touse2
	if      ("`tc'" == "ate") gen `touse2' = `touse' & `yl' != .
	else if ("`tc'" == "att") gen `touse2' = `w' == 1 & `touse' & `yl' != .
	else if ("`tc'" == "atc") gen `touse2' = `w' != 1 & `touse' & `yl' != .
	
	// Regression (MATA)
	mata: st_local("tau_est", strofreal(lm(st_data(.,"`touse2'"),st_data(.,"`tau'"),., 1)))
	
	
	// Calculates the adjusted average of squared yl
	if (length("`biasadj'") != 0) {
		
		#delimit ;
		mata: st_store(
			(1, `ntouse') , // Index
			st_varindex(tokens("`yl2'")),  // Varlist
			y2(st_matrix("`betas0'"),st_matrix("`betas1'"),&`matches', biasvars_index("`x'","`biasvars'"), `nvars') // Output
			);
		#delimit cr
	}

	// Dropping auxiliar mata matrix
	cap mata mata drop `matches'
	
	// Computing Variance
	qui gen `vari' = 0.5*(`yi'^2+`yl2'-2*`yi'*`yl'+`tau_est'^2-2*`tau_est'*(2*`w'-1)*(`yi'-`yl')) if `touse'
	summ `vari' if `touse', meanonly
	local sigma2 = r(mean)
	local n = r(N)
	qui gen `vari2' = (1+`kmi')^2*`sigma2'/`n' if `touse'
	summ `vari2' if `touse', meanonly
	local variance = r(mean)

	// Printing and storing results
	nnmatch2_results, yi("`yi'") w("`w'") x(`=trim("`x'")') m(`m') tc("`tc'") ///
		nobs(`ntouse') metric("`metric'") tau(`tau_est') touse(`touse') ///
		cmdline("`0'") level(`level') var(`variance') sigma2(`sigma2') bias(`biasvars')
	
end

// Stores results in e() and prints them as a table
cap program drop nnmatch2_results
program def nnmatch2_results, eclass

	syntax , yi(string) w(string) x(string) m(integer) tc(string) nobs(integer) ///
		metric(string) tau(real) cmdline(string) level(cilevel) touse(varname) ///
		var(real) sigma2(real) bias(string)
	
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
	
	if (length("`level'") == 0) local level = 95
	
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
	ereturn scalar N = `nobs'
	ereturn scalar m = `m'
	ereturn scalar se = `=sqrt(`var')'
	ereturn scalar sigma2 = `sigma2'
	local depn "`yi'"
	
	// Returning locals
	ereturn local cmdline "nnmatch2 `cmdline'"
	ereturn local depvar "`depn'"
	ereturn local match_ind "`w'"
	ereturn local match_vars "`x'"
	ereturn local stat "`estname'"
	ereturn local bias "`bias'"
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
	
	if (length("`bias'") >0) {
			di as text "{p 0 8 4}Bias-adj variables: "      ///
					"`bias'{p_end}"
	}

end

cap mata mata drop lm()
mata:
real colvector lm(
	real colvector touse,
	real colvector Y,
	real matrix X,
	| real scalar noconst,
	real colvector W
	) {
	
	real scalar nobs
	
	nobs = length(Y)
	
	if (X == .) X = J(nobs,1,1)
	
	// Checking noconst
	if (noconst == J(1,1,.)) noconst = 0
	if (!noconst) X = J(nobs,1,1),X
	
	// Checking weights
	if (W == J(0,1,.)) W = J(nobs,1,1)
	
	// Subseting
	if (touse != .) {
		Y = select(Y, touse)
		X = select(X, touse)
		W = select(W, touse)
	}

	// Returning LS
	return(invsym(quadcross(X,W,X))*quadcross(X,W,Y))
}

end


// General purpose matching in covariates.
// nn_match accepts single or multiple covariates as matching in one or more
// neighbours.
// cap mata: mata drop nn_match()
mata:
real matrix nn_match(
	real colvector Yi,
	real colvector W,
	real matrix X,
	string scalar metric,
	real scalar M, 
	pointer(real matrix) scalar xmatch
	)
	{
	
	// Objects declaration
	real colvector Yl, ids,ids_c,ids_t,idsm
	real scalar nvars, N, N1, N0, i, nmi, j
	real matrix S, V, Xl1, Xl0
	real scalar nmatchtot, outcols;
	
	// Defining constants and ids
	nvars = cols(X)
	N = rows(X)
	Yl  = J(N,4+nvars,0)
	
	// Adjusting for sd
	S = quadvariance(X)
		
	// Choosing metric
	if (regexm(metric,"maha|inverse|linear")) {
		if (metric == "maha") V =  invsym(S) // Mahalanobis matrix
		else if (metric == "linear") V = diag(J(nvars,1,1)) // Linear
		else V = diag(1:/S) // diagonal matrix of inverse variances
	}
	else V = st_matrix(metric)
	
	N1 = sum(W)
	N0 = N - N1
	
	// Colvector of indexes
	ids = 1::N
	_collate(ids, order(W, 1))
	
	// indexes of control units
	ids_c = ids[|1,1\N0,1|]
	_collate(ids_c, order(ids_c, 1))
	
	// indexes of treated units
	ids_t = ids[|(N0+1),1\N,1|]
	_collate(ids_t, order(ids_t, 1))
	
	// Covariates
	Xl1 = X[ids_t,.]
	Xl0 = X[ids_c,.]

	(*xmatch) = J(N*(max((N0,N1))-1),1,(*xmatch))
	nmatchtot = 0
	outcols = 5 + nvars*2 // cols((*xmatch))

	for(i=1;i<=N;i++) {
		// Starting values for i-th match
		if (!W[i]) {
		
			// ids of matched treated units
			idsm = bestmatches(calc_distances(Xl1, V, ids_t, X[i,.], N1), N1, M)
			
			// Mean matched Y
			Yl[i,1] = mean(Yi[idsm,.])
			
			// number of matches
			nmi = rows(idsm)
			Yl[i,2] = nmi
			
			// Mean matched Y^2
			Yl[i,3] = mean(Yi[idsm,.]:^2)
			
			// Number of units one is matched to
			Yl[idsm,4] = Yl[idsm,4] :+ 1/nmi
			
			// Average covariates among matched individuals
			Yl[|i,5\i,(4+nvars)|] = mean(X[idsm,.])
		
			// Adds a row to xmatch matrix throgh pointers
			// i, w(0), l, Yi, Yl, Xi1...Xin, <Xl1...Xln>
			(*xmatch)[|(nmatchtot + 1),1\(nmatchtot + nmi),outcols|] = (J(nmi,1,(i, 0)), idsm,J(nmi,1,Yi[i]),Yi[idsm,1], J(nmi, 1, X[i,.]), X[idsm,.]) 
			//(J(nmi,1,i), J(nmi, 1, 0), idsm,J(nmi,1,Yi[i]),Yi[idsm,1], J(nmi, 1, X[i,.]), X[idsm,.]) 
		}
		else {
			// ids of matched treated units
			idsm = bestmatches(calc_distances(Xl0, V, ids_c, X[i,.], N0), N0, M)
			
			// Mean matched Y
			Yl[i,1] = mean(Yi[idsm,.])
			
			// Number of matches
			nmi = rows(idsm)
			Yl[i,2] = nmi
			
			// Mean matched Y^2
			Yl[i,3] = mean(Yi[idsm,.]:^2)
			
			// Number of units one is matched to
			Yl[idsm,4] = Yl[idsm,4] :+ 1/nmi
			
			// Average covariates among matched individuals
			Yl[|i,5\i,(4+nvars)|] = mean(X[idsm,.])
		
			// Adds a row to xmatch matrix throgh pointers
			// i, w(0), l, Yl, Yl, kl, Xi1...Xin, <Xl1...Xln>
			(*xmatch)[|(nmatchtot + 1),1\(nmatchtot + nmi),outcols|] = (J(nmi,1,(i, 1)), idsm,J(nmi,1,Yi[i]),Yi[idsm,1], J(nmi, 1, X[i,.]), X[idsm,.])
		}
		
		// Increasing the number of matches made
		nmatchtot = nmatchtot + nmi
	}
	
	// Trimming the number of matches
	(*xmatch) = (*xmatch)[|1,1\nmatchtot,outcols|]
	
	return(Yl)
}
end

// Calculates distances from vector Xi towards each row of matrix Xl using a
// V weighting square matrix (metric matrix such as Inverse matrix or mahalanobis).
// It assumes that Xi vector is control/treated individual and every row in Xl
// is treated/control individual respectively.
// EF: returns a matrix with the distances in the first column and the Y values in the second
// cap mata: mata drop calc_distances()
mata:
real matrix calc_distances(
	real matrix Xl,
	real matrix V,
	real colvector ids,
	real rowvector Xi,
	real scalar Ni
	)
	{
	
	real scalar i
	real matrix xdif, xdif2, candidates
	
	xdif = Xl - J(Ni,1,Xi)
	
	candidates = J(Ni,1,.)
	// Computes distances
	for(i=1;i<=Ni;i++) {
			candidates[i,] = xdif[i,.]*V*transposeonly(xdif[i,.])
	}
	return((candidates, ids))
}
end

// Looks for at least the m best matches, in the case of two (or more) observations
// at the same distance from the current individual, bestmatch will keep both (or
// more) // so, just like nnmatch, nnmatch2's m() option could be more than what
// the user expects.
cap mata: mata drop bestmatches()
mata:
real colvector bestmatches(
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
		if ((candidates[i,1] != candidates[i+1,1]) | (M2 < M-1)) M2++
	}

	return(candidates[|1,2\(i-1),2|])
}
end

// Computes Yl and Yl^2 for variance estimator in case of Bias-adjusted estimator
cap mata: mata drop y2()
mata:
real colvector y2(
	real matrix betas0,
	real matrix betas1,
	pointer(real matrix) scalar xmatch,
	real rowvector biasvars,
	real scalar k
	)
{
	real scalar kbias, nobs, i, ny
	real matrix data
	real matrix xi, xl
	real colvector w, yl02, yl2, sumy2
	
	// Number of covariates
	kbias = length(biasvars)
	
	// Keeping only covariates for bias
	betas0 = betas0[1,1::kbias]
	betas1 = betas1[1,1::kbias]
	
	// Cooping matching database and deleting it from the memory
	data = *xmatch

	(*xmatch) = J(0,0,.)
	
	// Number of observations in the matching database
	nobs = rows(data)
	
	// Subseting covariates
	xi = data[|1,6\nobs,(5+k)|]
	xi = xi[.,biasvars]
	xl = data[|1,(6+k)\nobs,(5+2*k)|]
	
	xl = xl[.,biasvars]		

	// Subseting treatment
	w = data[|.,2|]
	
	// Adjusting for bias
	yl02 = (data[.,5]+w:*((xi-xl)*transposeonly(betas0))+(1:-w):*((xi-xl)*transposeonly(betas1))):^2
	
	// Starting empty colvectors of yl and yl^2
	yl2 = J(0,1,.)
	sumy2 = yl02[1,1]
	ny = 1
	for (i=2 ; i<=nobs ; i++) {
		if (data[i-1,1]==data[i,1]) {
			sumy2 = sumy2 + yl02[i,1]
			++ny
		}
		else {
			yl2 = yl2 \ (sumy2/ny)
			sumy2 = yl02[i,1]
			ny = 1
		}
	}
	return((yl2 \ (sumy2/ny)))
}
end

// Returns a permutation of X indexes which should be used in bias-adj
cap mata: mata drop biasvars_index()
mata:
real rowvector biasvars_index(
	string rowvector xvars,
	string rowvector biasvars
	)
	{
	
	real scalar i, j, k, kbias
	real rowvector indexes
	
	// Building string rowvectors (parsing)
	xvars = tokens(xvars)
	biasvars = tokens(biasvars)
	
	k = length(xvars)
	kbias = length(biasvars)
	indexes = J(1,kbias,.)
	
	// Building permutation
	for(i=1;i<=kbias;i++) {
		for(j=1;j<=k;j++) {
			if (xvars[j] == biasvars[i]) indexes[i] = j
		}
	
	}
	
	return(indexes)
}
	
end

/*

timer clear

cap log close _all
log using results_nnmatch2_v4.txt, text replace

// Lalonde
use ldw_exper, clear
//expand 5


// One match
timer on 1
nnmatch2 re78 t age educ black hisp married re74 re75 u74 u75 , m(1) bias(bias) /*save(prueba)*/
timer off 1

timer on 2
nnmatch  re78 t age educ black hisp married re74 re75 u74 u75 , m(1) bias(bias)
timer off 2

qui timer list
di "{hline}" _n "N of speedups (with con `=c(N)' obs) " %5.2fc `=r(t2)' "/" %5.2fc `=r(t1)' " secs. = x" %5.2fc `=r(t2)/r(t1)' _n "{hline}"
timer clear

// Four matches
timer on 1
nnmatch2 re78 t age educ black hisp married re74 re75 u74 u75 , m(2) bias(bias)
timer off 1

timer on 2
nnmatch  re78 t age educ black hisp married re74 re75 u74 u75 , m(2) bias(bias)
timer off 2

timer list
di "{hline}" _n "N of speedups (with con `=c(N)' obs) " %5.2fc `=r(t2)' "/" %5.2fc `=r(t1)' " secs. = x" %5.2fc `=r(t2)/r(t1)' _n "{hline}"
timer clear

////////////////////////////////////////////////////////////////////////////////
// Expanding data (estimators should keep unchanged)
expand 3

// One match
timer on 1
nnmatch2 re78 t age educ black hisp married re74 re75 u74 u75 if age > 20, m(1) bias(bias)
timer off 1

timer on 2
nnmatch  re78 t age educ black hisp married re74 re75 u74 u75 if age > 20, m(1) bias(bias)
timer off 2

qui timer list
di "{hline}" _n "N of speedups (with con `=c(N)' obs) " %5.2fc `=r(t2)' "/" %5.2fc `=r(t1)' " secs. = x" %5.2fc `=r(t2)/r(t1)' _n "{hline}"
timer clear

// Four matches
timer on 1
nnmatch2 re78 t age educ black hisp married re74 re75 u74 u75 , m(4) bias(bias)
timer off 1

timer on 2
nnmatch  re78 t age educ black hisp married re74 re75 u74 u75 , m(4) bias(bias)
timer off 2

timer list
di "{hline}" _n "N of speedups (with con `=c(N)' obs) " %5.2fc `=r(t2)' "/" %5.2fc `=r(t1)' " secs. = x" %5.2fc `=r(t2)/r(t1)' _n "{hline}"	

log close _all

// exit, clear

