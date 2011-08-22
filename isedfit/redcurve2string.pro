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
    case redcurve of
       -1: redcurvestring = ''
       0: redcurvestring = 'calzetti'
       1: redcurvestring = 'charlot'
       2: redcurvestring = 'odonnell'
       3: redcurvestring = 'smc'
       else: message, 'Unrecognized reddening curve'
    endcase
return, redcurvestring
end
