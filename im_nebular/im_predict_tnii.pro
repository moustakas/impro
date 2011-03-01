;+
; NAME:
;       IM_PREDICT_TNII()
;
; PURPOSE:
;       Predict T(N+) given T(O+).
;
; INPUTS:
;       toii - temperature of the T(O+) zone: [K, see IM_TEMDEN()]
;
; OUTPUTS:
;       tnii - T(N+) [K]
;
; OPTIONAL OUTPUTS:
;       tnii_err - error in TNII_OUT, if TOII_ERR is given 
;
; COMMENTS:
;       The (trivial!) relation T(N+) and T(O+) is given by Garnett
;       (1992), but see Nava et al. 2006 for an alternative argument. 
;
; EXAMPLES:
;
; MODIFICATION HISTORY:
;       J. Moustakas, 2007 Dec 07, NYU, written
;
; Copyright (C) 2007, John Moustakas
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

function im_predict_tnii, toii, toii_err=toii_err, tnii_err=tnii_err

    if (n_elements(toii) eq 0L) then begin
       doc_library, 'im_predict_tnii'
       return, -1L
    endif

    if (n_elements(toii_err) eq 0L) then toii_err = toii*0.0

    if (n_elements(toii) ne n_elements(toii_err)) then begin
       splog, 'Dimensions of TOII and TOII_ERR do not agree.'
       return, -1L
    endif
    
    tnii = toii*0.0-999.0
    tnii_err = toii*0.0-999.0
    good = where((toii gt 0.0),ngood)
    if (ngood ne 0L) then begin
       tnii[good] = toii[good]
       if arg_present(toii_err) then tnii_err[good] = toii_err[good]
    endif

return, tnii
end
