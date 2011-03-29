function redcurve2string, redcurve
; reddening curves
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
