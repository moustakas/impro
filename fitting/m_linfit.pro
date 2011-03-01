;+
; NAME:
;	M_LINFIT
; VERSION:
;	3.0
; PURPOSE:
;	Linear fitting with an arbitrary number of parameters.
; CATEGORY:
;	Mathematical Function.
; CALLING SEQUENCE:
;	Result = M_LINFIT( X, Y [,W] [, keywords])
; INPUTS:
;    X
;	Numeric vector.  X coordinates of the data.
;    Y
;	Numeric vector.  Y coordinates of the data.  Same length as X.
; OPTIONAL INPUT PARAMETERS:
;    W
;	Numeric vector.  Weight parameters.  Same length as X.  Default is
;	constant weight (1).
; KEYWORD PARAMETERS:
;    ORDER
;	Specifies the order of the fit, i.e. one less than the number of free 
;	parameters (or base functions).  If not given, the order will be read 
;	from the number of entries in the BASE variable (see below).  If both
;	ORDER and BASE are specified, the higher one will determine the order.
;	Finally, if neither is specified, ORDER is set to 1 (meaning linear fit)
;    RESIDUAL
;	Optional output, see below.
;    BASE
;	Character array containing names of base functions.  Any blank entry 
;	will be replaced by a power of X.  For example, if the third entry 
;	(i = 2) is blank (or null string) the third base function will be X^2.
;    PARAMS
;	Array containing optional parameters (one per function) to be used in
;	the function calls.  
;    PARMASK
;	Parameter mask, specifies which of the parameters are to be used.  A 
;	nonzero entry means USE, zero means DON'T USE.  Default is USE for all
;	existing parameters.
; OUTPUTS:
;	Returns a vector containing the fit coefficients.
; OPTIONAL OUTPUT PARAMETERS:
;    RESIDUAL
;	The name of the variable to receive the residual Chi-Square value, 
;	normalized to the number of degrees of freedom.
; COMMON BLOCKS:
;	None.
; SIDE EFFECTS:
;	None.
; RESTRICTIONS:
;	None.
; PROCEDURE:
;	Standard linear optimization, using a Singular Value Decomposition
;	procedure to solve the linear system.  Uses DEFAULT, FPU_FIX, 
;	SOLVE_LINSYS and TYPE from MIDL.
; MODIFICATION HISTORY:
;	Created 1-JUN-93 by Mati Meron.
;	Renamed from LINFIT to M_LINFIT to avoid clashes with an IDL library
;	routine bearing the same name.
;	Modified 20-SEP-1998 by Mati Meron.  Underflow filtering added.
;-

    on_error, 1
    nv = n_elements(x)
    if n_elements(y) ne nv then message, 'X ,Y lengths must be equal!'
    wei = Default(w,replicate(1.,nv),/dtype)
    if n_elements(wei) ne nv then message, 'X ,W lengths must be equal!'

    nord = Default(nord,1)
    bas = Default(bas, strarr(nord + 1))
    if Type(bas) ne 7 then message, 'Function names must be strings!'
    nbas = n_elements(bas)
    if nord gt (nbas-1) then begin
	bas = [bas,strarr(nord-nbas+1)]
	nbas = nord + 1
    endif else nord = nbas - 1
    bmsk = strtrim(bas) ne ''

    npars = n_elements(pars) < nbas
    if npars gt 0 then begin
	pmsk = Default(pmsk,replicate(1,npars)) ne 0
	npmsk = n_elements(pmsk) < npars
	if npmsk lt nbas then pmsk = [pmsk,replicate(0,nbas-npmsk)]
    endif else pmsk = replicate(0,nbas)

    farr = make_array(nv, nord + 1, type = ((Type(x) > Type(y)) > 4))
    for i = 0l, nord do begin
	if bmsk(i) then begin
	    if pmsk(i) then farr(*,i) = call_function(bas(i),x,pars(i)) $
	    else farr(*,i) = call_function(bas(i),x)
	endif else farr(*,i) = x^i
	farr(*,i) = FPU_fix(farr(*,i)*wei)
    endfor
    yw = FPU_fix(y*wei)

    res = Solve_linsys(farr,yw,/svd)
    resid = FPU_fix(sqrt(total((yw - farr#res)^2)/(nv - nbas)))
    return, res
end
