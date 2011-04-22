FUNCTION BVLS_Solve_pxf, A, b, npoly
compile_opt idl2, hidden

; No need to enforce positivity constraints if fitting one single template:
; use faster linear least-squares solution instead of BVLS.
;
s = size(a)
if s[0] eq 1 then $ ; A is a vector, not an array
    soluz = total(A*b)/total(A^2) $
else if s[2] eq npoly+1 then $ ; Fitting a single template
    soluz = la_least_squares(transpose(A),b) $
else begin               ; Fitting multiple templates
    bnd = dblarr(2,s[2],/NOZERO)
    mx = (machar()).xmax
    if npoly gt 0 then bnd[0,0:npoly-1] = -mx ; No bounds on Legendre polynomials
    bnd[0,npoly:*] = 0d  ; Positivity constraints on the templates (and sky spectra)
    bnd[1,*] = mx
    BVLS, A, B, bnd, soluz, ITMAX=15*s[2], IERR=ierr
    if ierr ne 0 then message, 'BVLS Error n. ' + strtrim(ierr,2)
endelse

return, soluz
END
;----------------------------------------------------------------------------
