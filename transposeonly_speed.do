mata:
mata clear
real matrix transp(real matrix x) {

	real scalar ncols, nrows, i, j
	real matrix y
	
	ncols = cols(x)
	nrows = rows(x)
	y = J(ncols, nrows, .)
	
	for(i=1;i<=ncols;i++) {
		for(j=1;j<=nrows;j++) {
			y[|i,j|] = x[|j,i|]
		}
	}
	
	return(y)
}

X = (
1,2,3\
4,5,6)

timer_clear()
for(i=1;i<=1000000;i++) {
	timer_on(1)
	X=transposeonly(X)
	timer_off(1)
	timer_on(2)
	X=transp(X)
	timer_off(2)
}

timer()
end

