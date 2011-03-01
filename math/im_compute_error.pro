function im_compute_error, x1, x1err, x2, x2err, log=log, $
  quotient=quotient, product=product
; jm02jul22uofa
; takes two variables and computes the error added in quadrature
; correctly.  assumes product
; jm02nov8uofa - added log keyword; also can be passed just x1 and
;                x1err (no x2 information)

    if (n_elements(log) eq 0L) and (n_elements(quotient) eq 0L) and $
      (n_elements(product) eq 0L) then product = 1L

    if (n_elements(x2) eq 0L) or (n_elements(x2err) eq 0L) then begin
       x2 = 1.0
       x2err = 0.0
    endif
    
; assume 
; x = alog10(x1/x2) = alog10(x1) - alog10(x2) or 
; x = alog10(x1*x2) = alog10(x1) + alog10(x2) 

    if keyword_set(log) then begin 

       term1 = (x1err/x1/alog(10))^2.0
       term2 = (x2err/x2/alog(10))^2.0

    endif

; assume x = x1/x2
    
    if keyword_set(quotient) then begin 
       
       term1 = x1err^2.0/x2^2.0
       term2 = x2err^2.0*(x1/x2^2.0)^2.0
    
    endif

; assume x = x1*x2
    
    if keyword_set(product) then begin 

       term1 = (x2*x1err)^2.0
       term2 = (x1*x2err)^2.0

    endif
    
    xerr = sqrt(term1 + term2)

return, xerr
end
