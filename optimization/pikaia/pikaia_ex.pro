function example_calc,n,params,_extra=extra
	; the fitness function will be a gaussian centered at x=extra.x and y = extra.y
	; params will be an array of size 2x?
	fitness = exp( -1.0*( (params[0,*]-extra.x)^2 + (params[1,*]-extra.y)^2 ) )

	; the result has to be transposed because the above expression
	; returns a column vector, but pikaia expects a row vector
	return,transpose(fitness)
end

function example_calc2,n,params
	; here's the highly optimized version
	return, exp( -1.0*total( (params-.5)^2, 1 ) )
end

function example_calc3,n,params
	; here's the highly non-optimized version
	npts = n_elements(params[0,*])
	fitness = fltarr(npts)

	for i=0,npts-1 do fitness[i] = exp( -1.0*( (params[0,i]-.5)^2 + (params[1,i]-.5)^2 ) )
	return,fitness
end

pro pikaia_ex

; the only thing we need to tell pikaia is the function name to use
; and the number of parameters to fit.  It will automatically load
; reasonable defaults for ctrl (population size: 100, generations: 200)
n = 2
t1 = systime(/seconds)
pikaia,n,ctrl,x,f,status,fname='example_calc',x=.4,y=.6
print,x

t2 = systime(/seconds)
pikaia,n,ctrl,x,f,status,fname='example_calc2'
print,x

t3 = systime(/seconds)
pikaia,n,ctrl,x,f,status,fname='example_calc3'
print,x

t4=systime(/seconds)
print,t2-t1
print,t3-t2
print,t4-t3
end