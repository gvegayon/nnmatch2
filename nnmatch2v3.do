clear all
cap set mem 2g

cap cd i:/george/comandos_paquetes_librerias/stata/nnmatch2
cap cd "C:\Users\eduardo.fajnzylber\Dropbox\Trabajo\publicaciones\proyecto nnmatch"

*! vers 0.13.04 12april2013
*! auth: George Vega
*! auth: Eduardo Fajnzylber
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
		
		// Printing and storing results
		nnmatch2_results, yi(`=e(depvar)') w(`=e(match_ind)') x(`=e(match_vars)') /// 
			m(`=e(m)') tc(`=lower(regexr(e(stat),"^[a-zA-Z]",""))') nobs(`=e(N)') ///
			metric(`=e(metric)') tau(`=`b'[1,1]') cmdline(`=e(cmdline)') level(`level') ///
			var(`=`V'[1,1]') sigma2(`=e(sigma2)')
	}
end

cap program drop nnmatch2main
program def nnmatch2main, eclass sortpreserve
	syntax varlist(min=3 numeric) [in] [if] [, METric(string) m(integer 1) Level(cilevel) tc(string) ///
		biasadj(string) save(string)]
	
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
	
	// Marking sample
	marksample touse
	markout `touse' `yi' `w' `x' // Marking complete cases
	
	
	// Counting the number of observations to be matched
	cap count if `touse'
	local ntouse = r(N)
	gsort -`touse'
	
	// Saving place for tempvars
	tempvar yl tau nmi yl2 kmi vari vari2 idi idl
	
	quietly {
		gen `idi' = _n if `touse'
		gen `yl' = . if `touse'  // Matched Y
		gen `nmi' = . if `touse' // Number of matched to i
		gen `yl2' = . if `touse' // Matched Y^2
		gen `kmi' = . if `touse' // Number of times i was matched to another obs
	}
	
	// i, w(0), l, Yi, Yl, Xi1, ..., Xin, Xi1, ..., Xin
	tempname matches
	mata: `matches' = J(0,5+`nvars'*2,.)
	
	// Matching algorithm	
	#delimit ;
	mata: st_store((1, `ntouse'), 
		(st_varindex(tokens("`yl' `nmi' `yl2' `kmi'"))), 
		nn_match(
			st_data((1,`ntouse'), "`yi'"), 
			st_data((1,`ntouse'),"`w'"),
			st_data((1,`ntouse'),"`x'"),
			"`metric'",
			`m',
			&`matches' // Pointer to matches matrix (used for computing biasadj)
			)
		);
	#delimit cr

	// If bias estimator
	if (length("`biasadj'") != 0) {
		tempvar yl_tmp
		
		tempfile originalset
		save `originalset'
		
		drop _all
		
		mata: st_local("nobs", strofreal(rows(`matches')))		
		
		// Loading data from matches matrix
		quietly {
			gen `idi' = .
			gen `w' = .
			gen `idl' = .
			gen `yi' = .
			gen `yl' = .
			
			// Generation i control vars
			forval j = 1/`nvars' {
				gen xi`j' = .
				local xis `xis' xi`j'
			}
			
			// Generation l control vars
			forval j = 1/`nvars' {
				gen xl`j' = .
				local xls `xls' xl`j'
			}
			
			set obs `=`nobs''
		}
		
		mata: st_store(. , ., `matches')
		mata: mata drop `matches'
		
		// Computing mus coeficients
		quietly {
			// Computing weights
			tempvar kmi isfirst
			
			bysort `idi': gen `kmi' = 1/_N
			bysort `idl': replace `kmi' = sum(`kmi')
			by `idl': gen `isfirst' = _n == 1
			
			tempname betas0 betas1
			
			qui reg `yl' xl* [aw=`kmi'] if `w' & `isfirst'
			mat def `betas1' = e(b)
			
			qui reg `yl' xl* [aw=`kmi'] if !`w' & `isfirst'
			mat def `betas0' = e(b)
			
			// Computing mus predictions
			// tempvar mu_0i mu_0l mu_1i mu_1l
		
			mat score mu_1l = `betas0'
			mat score mu_0l = `betas1'
			
			forval j = 1/`nvars' {
				drop xl`j'
				ren xi`j' xl`j'
			}
			
			mat score mu_1i = `betas0'
			mat score mu_0i = `betas1'
		}
		
		if (length("`save'") > 0) save `save', replace
		
		// Generating yl bias corrected of the form yl + u(xi) - u(xl)
		qui {
			replace `yl' = `yl' + `w'*(mu_0i - mu_0l) + (1-`w')*(mu_1i - mu_1l)
			collapse (mean) `yi' `yl', by(`idi') fast			

			merge 1:1 `idi' using `originalset', nogen keep(3)
		}
	}

	// Computing Tau
	qui gen `tau' = (`yi' - `yl')*(`w' == 1) + (`yl' - `yi')*(`w' != 1) if `touse'
	if ("`tc'" == "ate") {         // Average Treatment Effect
		qui reg `tau' if `touse' & `yl' != .
		local tau_est = _b[_cons]
	}
	else if ("`tc'" == "att") {    // Average Treatment Effect for the Treated
		qui reg `tau' if `w' == 1 & `touse' & `yl' != .
		local tau_est = _b[_cons]
	}
	else if ("`tc'" == "atc") {    // Average Treatment Effect for the Controls
		qui reg `tau' if `w' != 1 & `touse' & `yl' != .
		local tau_est = _b[_cons]
	}
	
	
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
		cmdline("`0'") level(`level') var(`variance') sigma2(`sigma2')
	
end

// Stores results in e() and prints them as a table
cap program drop nnmatch2_results
program def nnmatch2_results, eclass

	syntax , yi(string) w(string) x(string) m(integer) tc(string) nobs(integer) ///
		metric(string) tau(real) cmdline(string) level(cilevel) touse(varname) ///
		var(real) sigma2(real)
	
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

// General purpose matching in covariates.
// nn_match accepts single or multiple covariates as matching in one or more
// neighbours.
cap mata: mata drop nn_match()
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
	
	// Defining constants and ids
	nvars = cols(X)
	N = rows(X)
	Yl  = J(N,4,0)
	
	// Adjusting for sd
	S = quadvariance(X)
		
	// Choosing metric
	if (regexm(metric,"maha|inverse|linear")) {
		if (metric == "maha") {        // Mahalanobis matrix
			V =  invsym(S)
		}
		else if (metric == "linear") { // Linear
			V = diag(J(nvars,1,1))
		}
		else {                         // diagonal matrix of inverse variances
			V = diag(1:/S)
		}
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
	ids_t = ids[|N0+1,1\N,1|]
	_collate(ids_t, order(ids_t, 1))
	
	// Covariates
	Xl1 = X[ids_t,.]
	Xl0 = X[ids_c,.]

	for(i=1;i<=N;i++) {
		// Starting values for i-th match
		if (!W[|i|]) {
		
			// ids of matched treated units
			idsm = bestmatches(calc_distances(Xl1, V, ids_t, X[|i,.|], N1), N1, M)
			
			// Mean matched Y
			Yl[|i,1|] = mean(Yi[idsm,.])
			
			// number of matches
			nmi = rows(idsm)
			Yl[|i,2|] = nmi
			
			// Mean matched Y^2
			Yl[|i,3|] = mean(Yi[idsm,.]:^2)
			
			// number of units one is matched to
			Yl[idsm,4] = Yl[idsm,4] :+ 1/nmi
			
			// Adds a row to xmatch matrix throgh pointers
			// i, w(0), l, Yl, Yl, Xi1...Xin
			(*xmatch) = (*xmatch)\(J(nmi,1,i), J(nmi, 1, 0), idsm,J(nmi,1,Yi[|i|]),Yi[idsm,1], J(nmi, 1, X[|i,.|]), X[idsm,.])
		}
		else {
			// ids of matched treated units
			idsm = bestmatches(calc_distances(Xl0, V, ids_c, X[|i,.|], N0), N0, M)
			
			// Mean matched Y
			Yl[|i,1|] = mean(Yi[idsm,.])
			
			// number of matches
			nmi = rows(idsm)
			Yl[|i,2|] = nmi
			
			// Mean matched Y^2
			Yl[|i,3|] = mean(Yi[idsm,.]:^2)
			
			// number of units one is matched to
			Yl[idsm,4] = Yl[idsm,4] :+ 1/nmi
			
			// Adds a row to xmatch matrix throgh pointers
			// i, w(0), l, Yl, Yl, Xi1...Xin
			(*xmatch) = (*xmatch)\(J(nmi,1,i), J(nmi, 1, 1), idsm,J(nmi,1,Yi[|i|]),Yi[idsm,1], J(nmi, 1, X[|i,.|]), X[idsm,.])
		}
	}
	
	return(Yl)
}
end

// Calculates distances from vector Xi towards each row of matrix Xl using a
// V weighting square matrix (metric matrix such as Inverse matrix or mahalanobis).
// It assumes that Xi vector is control/treated individual and every row in Xl
// is treated/control individual respectively.
// EF: returns a matrix with the distances in the first column and the Y values in the second
cap mata: mata drop calc_distances()
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
	real matrix xdif, candidates
	
	xdif = Xl - J(Ni,1,Xi)
	
	candidates = (J(Ni,1,.),ids)
	
	// Computes distances
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
		if ((candidates[|i,1|] != candidates[|i+1,1|]) | (M2 < M-1)) M2++
	}

	return(candidates[|1,2\i-1,2|])
}
end


timer clear

cap log close _all
log using results_nnmatch2_v3.txt, text replace

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

// One match
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
nnmatch2 re78 t age educ black hisp married re74 re75 u74 u75 , m(2)
timer off 1

timer on 2
nnmatch  re78 t age educ black hisp married re74 re75 u74 u75 , m(2)
timer off 2

qui timer list
di "{hline}" _n "N of speedups (with con `=c(N)' obs) " %5.2fc `=r(t2)' "/" %5.2fc `=r(t1)' " secs. = x" %5.2fc `=r(t2)/r(t1)' _n "{hline}"
timer clear

////////////////////////////////////////////////////////////////////////////////
// Expanding data (estimators should keep unchanged)
expand 3

// One matche
timer on 1
nnmatch2 re78 t age educ black hisp married re74 re75 u74 u75 if age > 20, m(1)
timer off 1

timer on 2
nnmatch  re78 t age educ black hisp married re74 re75 u74 u75 if age > 20, m(1)
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
