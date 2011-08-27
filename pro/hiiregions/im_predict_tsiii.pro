;+
; NAME:
;       IM_PREDICT_TSIII()
;
; PURPOSE:
;       Predict T(S++) given T(O++).
;
; INPUTS:
;       toiii - temperature of the T(O++) zone: [K, see IM_TEMDEN()]
;
; OUTPUTS:
;       siii - T(S++) [K]
;
; OPTIONAL OUTPUTS:
;       siii_err - error in SIII_OUT, if TOIII_ERR is given [K]
;
; COMMENTS:
;       The relation T(S++) and T(O++) is given by Garnett (1992).
;       Note, however, that Perez-Montero & Diaz (2003) argue that
;       because of updates to the [S III] collisional excitation rates
;       (Tayal & Gupta 1999) a better relation would be: T(O++) = 
;       9500*T(S++)+800.
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

function im_predict_tsiii, toiii, toiii_err=toiii_err, tsiii_err=tsiii_err

    if (n_elements(toiii) eq 0L) then begin
       doc_library, 'im_predict_tsiii'
       return, -1L
    endif

    if (n_elements(toiii_err) eq 0L) then toiii_err = toiii*0.0

    if (n_elements(toiii) ne n_elements(toiii_err)) then begin
       splog, 'Dimensions of TOIII and TOIII_ERR do not agree.'
       return, -1L
    endif
    
    tsiii = toiii*0.0-999.0
    tsiii_err = toiii*0.0-999.0
    good = where((toiii gt 0.0),ngood)
    if (ngood ne 0L) then begin
       tsiii[good] = 0.83*toiii[good] + 1700.0
       if arg_present(tsiii_err) then tsiii_err[good] = 0.83*toiii_err[good]
    endif

return, tsiii
end
