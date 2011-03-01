;+
; NAME:
;	COMPUTE_PARANGLE()
;
; PURPOSE:
;	Compute the parallactic angle.
;
; CALLING SEQUENCE:
;	pa = compute_parangle(latitude,dec,ha)
;
; INPUTS:
;	latitude - observatory latitude [decimal degrees]
;	dec      - declination [decimal degrees]
;	ha       - hour angle [decimal degrees]
;
; OUTPUTS:
;	parangle - parallactic angle [degrees]
;
; COMMENTS:
;	The parallactic angle is returned in the range [0-180]
;	degrees. 
;
; EXAMPLE:
;
; PROCEDURES USED:
;
; MODIFICATION HISTORY:
;	J. Moustakas, 2002 Feb 13, U of A
;-

function compute_parangle, latitude, dec, ha

    if n_params() ne 3L then begin
       print, 'Syntax - parangle = compute_parangle(latitude,dec,ha)'
       return, 0
    endif
    
    top = sin(ha*!dtor)
    bot = tan(latitude*!dtor)*cos(dec*!dtor)-sin(dec*!dtor)*cos(ha*!dtor)
    parangle = atan(top/bot)*!radeg
    
    paneg = where(parangle lt 0.0,nneg)
    if nneg ne 0L then parangle[paneg] = 360 + parangle[paneg] 
    papos = where(parangle gt 180.0,npos)
    if npos ne 0L then parangle[papos] = parangle[papos] - 180D

return, parangle
end
