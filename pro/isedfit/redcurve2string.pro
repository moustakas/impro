;+
; NAME:
;   REDCURVE2STRING()
; PURPOSE:
;   Map an iSEDfit reddening number to the appropriate reddening
;   curve. 
; COMMENTS:
;   ISEDFIT support routine.
; MODIFICATION HISTORY:
;   J. Moustakas, ???
;-

function redcurve2string, redcurve

    ncurve = n_elements(redcurve)
    if (ncurve eq 0) then begin
       doc_library, 'redcurve2string'
       return, -1
    endif

    if (ncurve gt 1) then begin
       str = strarr(ncurve)
       for ii = 0, ncurve-1 do str[ii] = redcurve2string(redcurve[ii])
       return, str
    endif
    
    case redcurve of
       -1: redcurvestring = 'none'
       0: redcurvestring = 'calzetti'
       1: redcurvestring = 'charlot'
       2: redcurvestring = 'odonnell'
       3: redcurvestring = 'smc'
       else: message, 'Unrecognized reddening curve'
    endcase
return, redcurvestring
end
