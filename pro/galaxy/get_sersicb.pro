;+
; NAME:
;   GET_SERSICB()
;
; PURPOSE:
;   Compute the Sersic bn parameter.
;
; INPUTS: 
;   sersicn - Sersic "n" value
;
; COMMENTS:
;   See Graham & Driver (2005), equations 1 and 4.
;
; MODIFICATION HISTORY:
;   J. Moustakas ???
;-

function sersicb_func, bb
    common sersicb, nn
return, gamma(2.0*nn)-2D*igamma(2*nn,bb)*gamma(2*nn)
end

function get_sersicb, sersicn
    common sersicb, nn
    nn = (sersicn>0D)<20D
    if nn le 0D or nn ge 20D then stop
    bn = zbrent(0.0,20.0,func_name='sersicb_func',$
      max_iterations=50,tolerance=1E-2)
return, bn
end
