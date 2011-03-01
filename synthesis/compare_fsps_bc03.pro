pro compare_fsps_bc03
; jm10jan31ucsd - compare the FSPS and BC03 SSPs
    
    ff = im_read_fsps(/chab)
    bc = im_read_bc03()
    bcflux = interpolate(bc.flux,findex(bc.wave,ff.wave),$
      findex(bc.age,ff.age),/grid)

    psfile = getenv('IMPRO_DIR')+'/synthesis/compare_fsps_bc03.ps'
    im_plotconfig, 0, pos, psfile=psfile, height=5, width=6.5, xmargin=[1.8,0.2]
;   for ii = 30, n_elements(ff.age)-1 do begin
    for ii = 0, n_elements(ff.age)-1 do begin
;      bcflux[*,ii] = bcflux[*,ii]*ff.wave^2/im_light(/ang)
       djs_plot, ff.wave, ff.flux[*,ii], xrange=[85,7000], xsty=1, $
         ysty=1, /xlog, xtitle='Wavelength (\AA)', $ ; /ylog, $
         ytitle='Flux (erg s^{-1} \AA^{-1} M_{\odot}^{-1})', $
         position=pos, yrange=[min(ff.flux[*,ii])>min(bcflux[*,ii]),$
         max(ff.flux[*,ii])>max(bcflux[*,ii])]
       djs_oplot, ff.wave, bcflux[*,ii], color='red'
       im_legend, strtrim(string(ff.age[ii]/1D6,format='(F12.3)'),2)+' Myr', $
         /left, /top, box=0
;      cc = get_kbrd(1)
    endfor
    im_plotconfig, psfile=psfile, /psclose, /pdf
    
return
end
    
