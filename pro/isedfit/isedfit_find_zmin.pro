function zmin_parabola, xx, yy
; fit a parabolla; support routine for ISEDFIT_FIND_ZMIN
    ss = sort(xx)
    xx = xx[ss]
    yy = yy[ss]

    a = xx[0]
    b = xx[1]
    c = xx[2]

    fa = yy[0]
    fb = yy[1]
    fc = yy[2]

    xmin = b - 0.5 * ( (b-a)^2*(fb-fc)-(b-c)^2*(fb-fa) ) / $
      ( (b-a)*(fb-fc) - (b-c) *(fb-fa))
return, xmin
end

function isedfit_find_zmin, z, ochi2, nmin=nmin
; jm14jun14siena - find the first N minima in a plot of redshift vs
; chi2; based on code written by Richard Cool for the PRIMUS survey

    if n_elements(nmin) eq 0 then nmin = 5
    
    minmask = intarr(n_elements(ochi2))
    npix = n_elements(z)
    dchi2 = deriv(findgen(npix),ochi2)
    
    delvarx, minval
    delvarx, minima
    
    while total(minmask) lt npix and n_elements(minval) lt nmin do begin
       chi2 = ochi2 + 1D64*minmask
       junk = min(chi2,jmin)
       lgrow = 0
       stoplgrow = 0
       while stoplgrow eq 0 do begin
          lgrow++
          if jmin-lgrow lt 0 then stoplgrow =1 else begin
             if dchi2[jmin-lgrow] ge 0 then stoplgrow = 1
          endelse
       endwhile
       minmask[jmin-(lgrow-1):jmin] = 1
       
       rgrow = 0
       stoprgrow = 0
       while stoprgrow eq 0 do begin
          rgrow++
          if jmin+rgrow gt npix-1 then stoprgrow = 1 else begin
             if dchi2[jmin+rgrow] lt 0 then stoprgrow = 1
          endelse
       endwhile
       minmask[jmin:jmin+(rgrow-1)] = 1
       
       range = [jmin-1,jmin+1]
       if range[0] lt 0 then range = [0,2]
       if range[1] gt npix-1 then range = [npix-3,npix-1]
       
       range1 = [jmin-2,jmin+2]
       if range1[0] lt 0 then range1 = [0,4]
       if range1[1] gt npix-1 then range1 = npix-[5,1]
       zmin = zmin_parabola(z[range[0]:range[1]], $
         ochi2[range[0]:range[1]])
       mval = interpol(ochi2[range[0]:range[1]], $
         z[range[0]:range[1]],zmin)
       err = interpol(abs(z[range1[0]:range1[1]]-zmin), $
         ochi2[range1[0]:range1[1]],mval+1.0)
       
       if jmin lt 2 or jmin gt n_elements(z)-2 then begin
          zmin = z[jmin]
          err = -1
          mval = ochi2[jmin]
       endif
       
       if n_elements(minima) eq 0 then begin
          minima = zmin
          minval = mval
          minerr = err
          rminima = z[jmin]
          rminval = ochi2[jmin]
       endif else begin
          minima = [minima,zmin]
          minval = [minval,mval]
          minerr = [minerr,err]
          rminima = [rminima,z[jmin]]
          rminval = [rminval,ochi2[jmin]]
       endelse
    endwhile

; pack the results into a structure    
    if n_elements(minval) gt 0 then begin
       ss = sort(minval)
       minval = minval[ss]
       minima = minima[ss]
       minerr = minerr[ss]
       
       output = {$
         chi2min:      0.0,$
         photoz:       0.0,$
         photoz_err:   0.0,$
         photoz_sigma: [0.0,0.0],$
         photoz_95:    [0.0,0.0],$
         edges:           0}
       output = replicate(output,n_elements(ss))
       output.chi2min = minval
       output.photoz = minima
       output.photoz_err = minerr

; temporary hack; drop redshift minima that are on or near the edge
; (the right way to do this is with an edge bit flag)
       keep = where(output.photoz_err gt 0.0)
       output = output[keep]

; get the asymmetric 1-sigma uncertainty on just the first minimum 
       onesig = errorf(1/sqrt(2))/2 ; =0.341
       pofz = exp(-0.5*ochi2)
       pofz = pofz/total(pofz)
       cumupofz = total(pofz,/cumu) ; assumed properly normalized
       for ii = 0, n_elements(output)-1 do begin
;         pofzmin = interpol(cumupofz,z,output[ii].photoz)
;         zup = interpolate(z,findex(cumupofz,pofzmin+onesig))
;         zlo = interpolate(z,findex(cumupofz,pofzmin-onesig))
; 1-sigma
          zup = interpolate(z,findex(cumupofz,0.5+onesig))
          zlo = interpolate(z,findex(cumupofz,0.5-onesig))
          output[ii].photoz_sigma = [zup-output[ii].photoz,output[ii].photoz-zlo]
; 95% CL
          zup = interpolate(z,findex(cumupofz,0.5+0.95/2))
          zlo = interpolate(z,findex(cumupofz,0.5-0.95/2))
          output[ii].photoz_95 = [zup-output[ii].photoz,output[ii].photoz-zlo]
       endfor
       return, output 
    endif else return, -1

end
