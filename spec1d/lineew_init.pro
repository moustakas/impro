function lineew_init, linewave, linename
; function for IM_SPLOTEW()
    lineew = {$
      linename:        linename, $
      linewave:        linewave, $
      linecenter:      0.0,$
      linecenter_flux: 0.0,$
      linewidth:       0.0,$
      continuum:       [0.0,-1.0], $
      gaussflux:       [0.0,-1.0], $
      ew:              [0.0,-1.0]}
return, lineew
end
