;+
; NAME:
;       IM_PREDICT_TOII()
;
; PURPOSE:
;       Predict T(O+) given T(O++).
;
; INPUTS:
;       toiii - temperature of the T(O++) zone: [K, see IM_TEMDEN()] 
;
; OUTPUTS:
;       toii - T(O+) [K]
;
; OPTIONAL OUTPUTS:
;       toii_err - error in TOII, if TOIII_ERR is given [K]
;
; COMMENTS:
;       The relation T(O+) and T(O++) used here was originally derived
;       by Campbell et al. (1986), based on a fit to the Stasinska
;       (1982) photoionization models, but was popularized by Garnett
;       (1992). 
;
;       Izotov et al. (1994) also give a relation between T(O+) and
;       T(O++), based on the Stasinska (1990) photoionization models,
;       given by T(O+) = 2E4/(1E4/T(O++)+0.8).  However, Bresolin et
;       al. 2004 note that if the Stasinska models are limited to
;       T_eff<40,000 K stars, then there is good agreement with the
;       Garnett (1992) relation.
;
;       Pagel et al. (1992) have also published a relation.
;
; EXAMPLES:
;
; MODIFICATION HISTORY:
;       J. Moustakas, 2004 July 5, U of A, written
;       jm07dec06nyu - removed STASINSKA relation and improved the
;                      documentation 
;
; Copyright (C) 2004, 2007, John Moustakas
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

function im_predict_toii, toiii, toiii_err=toiii_err, toii_err=toii_err

    if (n_elements(toiii) eq 0L) then begin
       doc_library, 'im_predict_toii'
       return, -1L
    endif

    if (n_elements(toiii_err) eq 0L) then toiii_err = toiii*0.0

    if (n_elements(toiii) ne n_elements(toiii_err)) then begin
       splog, 'Dimensions of TOIII and TOIII_ERR do not agree.'
       return, -1L
    endif
    
    toii = toiii*0.0-999.0
    toii_err = toiii*0.0-999.0
    good = where((toiii gt 0.0),ngood)
    if (ngood ne 0L) then begin
       toii[good] = 0.7*toiii[good] + 3000.0
       if arg_present(toii_err) then toii_err[good] = 0.7*toiii_err[good]
    endif
       
return, toii
end
