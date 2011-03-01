function redcurve2string, redcurve, params=params
; reddening curves
    if (params.ebv[1] le 0) then redcurvestring = '' else begin
       case redcurve of
          0: redcurvestring = 'calzetti'
          1: redcurvestring = 'charlot'
          2: redcurvestring = 'odonnell'
          3: redcurvestring = 'smc'
          else: message, 'Unrecognized reddening curve'
       endcase
    endelse
return, redcurvestring
end
