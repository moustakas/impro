function hessian, x, func, h=h, func_args=func_args

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
; Compute the Hessian of a function numerically via
; finite-differencing.
;
; Author : Brandon C. Kelly
;
; INPUTS:
;   X - The independent variable, an NX-element vector or scalar.
;   FUNC - A string variable containing the name of the function to
;          have its derivatives evaluated.  This should be an IDL
;          function that takes as argument and NX-element vector or
;          scalar, X, and returns a scalar, Y.
;
; OPTIONAL INPUTS:
;   H - The perturbation to X used in the finite differencing, i.e., 
;       X -> X +/- h.  Can be an NX-element vector or a scalar. 
;
; OUTPUT:
;   The matrix of second partial derivatives, an [NX, NX] array. HESSIAN[i,j] contains
;   the partial derivative of Y w.r.t. X[i,j].
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

if n_params() lt 2 then begin
    print, 'Syntax- RESULT = HESSIAN(X, FUNC, H=H, FUNC_ARGS=FUNC_ARGS)'
    return, 0
endif

xinfo = size(x)
if xinfo[0] gt 1 then begin
    print, 'X must be a vector or a scalar.'
    return, 0
endif
nx = xinfo[1]
if n_elements(h) eq 0 then begin
    h = 1d-5
endif
if n_elements(h) eq 1 then h = replicate(h, nx)
if n_elements(h) ne nx then begin
    print, 'H must be scalar or a vector of size equal to X.'
    return, 0
endif

x = double(x)

if n_elements(func_args) gt 0 then pass = 1 else pass = 0

hessian = dblarr(nx, nx)
for i = 0, nx - 1 do begin
    
    for j = i, nx - 1 do begin

        dx = dblarr(nx)
        dx[i] = dx[i] + h[i]
        dx[j] = dx[j] + h[j]

        y1 = pass ? call_function(func, x + dx, _extra=func_args) : call_function(func, x + dx)
        
        dx[i] = dx[i] - 2 * h[i]
        y2 = pass ? call_function(func, x + dx, _extra=func_args) : call_function(func, x + dx)
        
        dx[i] = dx[i] + 2 * h[i]
        dx[j] = dx[j] - 2 * h[j]
        y3 = pass ? call_function(func, x + dx, _extra=func_args) : call_function(func, x + dx)
        
        dx[i] = dx[i] - 2 * h[i]
        y4 = pass ? call_function(func, x + dx, _extra=func_args) : call_function(func, x + dx)

        hessian[i,j] = (y1 - y2 - y3 + y4) / (4.0 * h[i] * h[j])

        if i ne j then hessian[j,i] = hessian[i,j]
        
    endfor
    
endfor

return, hessian
end
