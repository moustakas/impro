pro test_ml_abundance, hii, _extra=extra
; jm05mar21uofa
; test the maximum likelihood abundance code

    if (n_elements(hii) eq 0L) then begin
       hii = read_hii_regions(/nolog)
;      keep = where(hii.zt_log12oh gt -900)
;      hii = hii[keep]
    endif
    nobj = n_elements(hii)
    
    rationames = ['OIII_5007_OII','OIII_5007_H_BETA','NII_6584_OII','R23']
;   rationames = ['OII_H_BETA','OIII_5007_H_BETA','NII_6584_H_BETA','SII_H_BETA']
    nratios = n_elements(rationames)

    fluxratios = fltarr(nratios,nobj)
    ferrratios = fltarr(nratios,nobj)
    for iratio = 0L, nratios-1L do begin
       fluxratios[iratio,*] = (struct_trimtags(hii,select=rationames)).(iratio)
       ferrratios[iratio,*] = (struct_trimtags(hii,select=rationames+'_ERR')).(iratio)
    endfor

    invvar = fluxratios*0.0+1D-32
    
    for i = 0L, nratios-1L do begin
       bad = where(fluxratios[i,*] lt -900.0,nbad,comp=good,ncomp=ngood)
       if (nbad ne 0L) then begin
          fluxratios[i,bad] = 0.0
          ferrratios[i,bad] = 1E16
       endif
       if (ngood ne 0L) then invvar[i,good] = 1.0/ferrratios[i,good]^2
    endfor

    result = ml_abundance(fluxratios,invvar,rationames,_extra=extra)

    good = where(hii.z_12oh_o3n2_pettini gt -900)
    plot, hii[good].z_12oh_o3n2_pettini, result[good].zml_12logoh, ps=4, $
      xr=[7.5,9.3], yr=[7.5,9.3], xsty=3, ysty=3
    oplot, !x.crange, !y.crange, line=0
    
stop    
    
return
end
