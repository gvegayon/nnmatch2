// this file does all the examples in the Stata Journal article
// "Implementing Matching Estimators for Average Treatment Effects in Stata"

clear
capture log close

// example 1
sjlog using match_ex1, replace
use artificial
nnmatch y w x
sjlog close, replace

// example 2
sjlog using match_ex2, replace
use ldw_exper
nnmatch re78 t age educ black hisp married re74 re75 u74 u75, m(4)
sjlog close, replace

// example 3
sjlog using match_ex2b, replace
nnmatch re78 t age educ black hisp married re74 re75 u74 u75, m(4) pop
sjlog close, replace

// example 4
sjlog using match_ex3, replace
nnmatch re78 t age educ black hisp married re74 re75 u74 u75, tc(att) m(4)
sjlog close, replace

// example 5
sjlog using match_ex4, replace
nnmatch re78 t age educ black hisp married re74 re75 u74 u75, tc(att) 
sjlog close, replace

// example 6
sjlog using match_ex5, replace
nnmatch re78 t age educ black hisp married re74 re75 u74 u75, tc(att) 	///
	m(4) bias(bias)
sjlog close, replace

// example 7
sjlog using match_ex6, replace
nnmatch re78 t age educ black hisp married re74 re75 u74 u75, tc(att) 	///
	m(4) robusth(4)
sjlog close, replace


