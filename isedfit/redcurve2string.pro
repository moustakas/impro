function redcurve2string, redcurve, params=params
; reddening curves

    case redcurve of
       0: redcurvestring = 'calzetti'
       1: redcurvestring = 'charlot'
       2: redcurvestring = 'odonnell'
       3: redcurvestring = 'smc'
       else: message, 'Unrecognized reddening curve'
    endcase

    if (n_elements(params) ne 0) then begin
       if (params.ebv[1] le 0) then redcurvestring = ''
    endif

return, redcurvestring
end
