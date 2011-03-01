function sfr_09rieke, z, f24
; jm09nov23ucsd - compute the SFR from the observed 24-micron *flux*
;   density (this relation includes the K-correction, for various SED 
;   types, from 0<z<3)

    common rieke_table1, ztable, atable, btable
    
    nobj = n_elements(z)
    if (nobj ne n_elements(f24)) then begin
       splog, 'Z and F24 must be the same dimension'
       return, -1
    endif
    
;   if (nobj gt 1) then begin
;      sfr = fltarr(nobj)
;      for ii = 0, nobj-1 do sfr[ii] = $
;        sfr_09rieke(z[ii],f24[ii])
;      return, sfr
;   endif

; cache Table 1 from the paper    
    if (n_elements(ztable) eq 0) or (n_elements(atable) eq 0) or $
      (n_elements(btable) eq 0) then begin
       ztable = [0.0,0.2,0.4,0.6,0.8,1.0,1.2,1.4,$
         1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0]
       atable = [0.417,0.502,0.528,0.573,0.445,0.358,0.505,$
         0.623,0.391,0.072,0.013,0.029,0.053,0.162,0.281,0.371]
       btable = [1.032,1.169,1.272,1.270,1.381,1.565,1.745,$
         1.845,1.716,1.642,1.639,1.646,1.684,1.738,1.768,1.782]
    endif

; apply eq. (14); f24 should be in Jy
    dlum = dluminosity(z,/cm)
    az = interpol(atable,ztable,z)
    bz = interpol(btable,ztable,z)
    sfr = az + bz*(alog10(4.0*!dpi*dlum*dlum*f24)-53.0) ; log-SFR [M_sun/yr]

return, float(sfr)
end
    
