mata:
mata clear

void fun(pointer(real colvector) scalar x) {
	
	real scalar i, ncols

	ncols = rows((*x))

	for(i=1;i<=ncols;i++) {
		timer_on(1)	
		(*x)[i] = (*x)[i]*(-1)
		timer_off(1)
	}
	
	i = mean((*x))
}

real colvector fun2(real colvector x) {
	
	real scalar i, ncols

	ncols = rows(x)

	for(i=1;i<=ncols;i++) {
		timer_on(1)	
		x[i] = x[i]*(-1)
		timer_off(1)
	}
	
	i = mean(x)
	
	return(x)
}

x = 1::10000

timer_clear()
for(j=1;j<=1000;j++) {
	fun(&x)
}
timer()

timer_clear()
for(j=1;j<=1000;j++) {
	x = fun2(x)
}
timer()

end
