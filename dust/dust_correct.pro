;+
; NAME:
;   DUST_CORRECT()
;
; PURPOSE:
;   Correct a line flux measurement for dust extinction.  
;
; INPUTS:
;   lineflux  - emission line flux [NFLUX]
;   linewave  - emission line wavelength (scalar or [NFLUX] array)
;
; OPTIONAL INPUTS:
;   decrement     - Balmer decrement (H-alpha/H-beta or
;                   H-alpha/H-gamma) [NFLUX]
;   err_decrement - corresponding error in DECREMENT [NFLUX]
;   ebv           - color excess E(B-V) [NFLUX]
;   err_ebv       - corresponding error in EBV [NFLUX]
;
; KEYWORD PARAMETERS:
;   extra         - keywords for GET_EBV() and K_LAMBDA() 
;
; OUTPUTS:
;   newflux       - extinction-corrected LINEFLUX [NFLUX]
;
; OPTIONAL OUTPUTS:
;   err_newflux   - if either ERR_DECREMENT or ERR_EBV are passed 
;                   then compute the error in NEWFLUX [NFLUX]
;
; COMMENTS:
;   If both DECREMENT and EBV are passed then DECREMENT is used to
;   calculate and overwrite EBV.
;
; PROCEDURES USED:
;   GET_EBV(), K_LAMBDA(), SPLOG
;
; EXAMPLE:
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2002 May 28, U of A
;   jm04jul08uofa - bug fix in returning EBV_ERR from GET_EBV();
;                   variable name shuffle
;
; Copyright (C) 2002, 2004, John Moustakas
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

function dust_correct, lineflux, linewave, err_lineflux=err_lineflux, $
  decrement=decrement, err_decrement=err_decrement, ebv=ebv, $
  err_ebv=err_ebv, err_newflux=err_newflux, _extra=extra

    nflux = n_elements(lineflux)
    nlines = n_elements(linewave)

    if (nflux eq 0L) or (nlines eq 0L) then begin
       doc_library, 'dust_correct'
       return, -1L
    endif
    
    if (nlines ne 1L) and (nlines ne nflux) then begin
       splog, 'LINEWAVE must be a scalar quantity or an [NFLUX] array.'
       return, -1L
    endif

    nfluxerr = n_elements(err_lineflux)
    if nfluxerr eq 0L then err_lineflux = lineflux*0.0 else begin
       if nfluxerr ne nflux then begin
          splog, 'LINEFLUX and ERR_LINEFLUX have incompatible dimensions.'
          return, -1L
       endif 
    endelse
       
    ndec = n_elements(decrement)
    nebv = n_elements(ebv)

    if (ndec eq 0L) and (nebv eq 0L) then begin
       splog, 'Either DECREMENT or EBV must be specified.'
       return, -1L
    endif

    if (nebv ne 0L) then begin
       if nebv ne nflux then begin
          splog, 'EBV and LINEFLUX have incompatible dimensions.'
          return, -1L
       endif
    endif

    if (ndec ne 0L) then begin
       if ndec ne nflux then begin
          splog, 'DECREMENT and LINEFLUX have incompatible dimensions.'
          return, -1L
       endif
; calculate the color excess and error
       ebv = get_ebv(decrement,decrement_err=err_decrement,ebv_err=err_ebv,_extra=extra)
    endif

; boost the flux appropriately and compute the error

    kl = k_lambda(linewave,_extra=extra)
    newflux = lineflux*10D0^(0.4*ebv*kl)

; include the error in E(B-V)

    err_newflux = sqrt( (err_lineflux*10.0^(0.4*ebv*kl))^2.0 + $
      (lineflux*0.4*kl*alog(10)*10.0^(0.4*ebv*kl)*err_ebv)^2.0 )

; do not include the error in E(B-V)

;   newflux_err = err_lineflux*10.0^(0.4*ebv*kl)

return, newflux
end


