;+
; NAME:
;   IM_COMPUTE_ERROR()
;
; PURPOSE:
;   Perform propagation of errors with one or two variables. 
;
; INPUTS: 
;   x1 - first variable
;   x1err - corresponding error
;
; OPTIONAL INPUTS: 
;   x2 - second variable
;   x2err - corresponding error
;   
; KEYWORD PARAMETERS: 
;   log - assume:
;     x = alog10(x1/x2) = alog10(x1) - alog10(x2) or 
;     x = alog10(x1*x2) = alog10(x1) + alog10(x2) 
;   quotient - assume x = x1/x2
;   product - assume x = x1*x2 (default)
;
; OUTPUTS: 
;   xerr - quadrature sum 
;
; COMMENTS:
;
; EXAMPLES:
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2002 Jul 22, U of A
;   jm02nov8uofa - added log keyword; also can be passed just x1 and
;     x1err (no x2 information) 
;
; Copyright (C) 2002, John Moustakas
; 
; This program is free software; you can redistribute it and/or modify 
; it under the terms of the GNU General Public License as published by 
; the Free Software Foundation; either version 2 of the License, or
; (at your option) any later version. 
; 
; This program is distributed in the hope that it will be useful, but 
; WITHOUT ANY WARRANTY; without even the implied warranty of
; MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
; General Public License for more details. 
;-

function im_compute_error, x1, x1err, x2, x2err, log=log, $
  quotient=quotient, product=product

    if (n_elements(log) eq 0) and (n_elements(quotient) eq 0) and $
      (n_elements(product) eq 0) then product = 1

    if (n_elements(x2) eq 0L) or (n_elements(x2err) eq 0L) then begin
       x2 = 1.0
       x2err = 0.0
    endif
    
; assume 
;   x = alog10(x1/x2) = alog10(x1) - alog10(x2) or 
;   x = alog10(x1*x2) = alog10(x1) + alog10(x2) 
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
