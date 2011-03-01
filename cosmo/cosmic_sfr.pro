function cosmic_sfr, zmin=zmin, zmax=zmax, doplot=doplot
; jm05jun06uofa - integrate the cosmic SFR density as parameterized in
;                 Bell (2004)

    z = findgen(61)/10.0
    sfr = (0.006+0.072*z^1.35)/(1+(z/2.0)^2.4)
    
    if (n_elements(zmin) eq 0L) then zmin = 0.0
    if (n_elements(zmax) eq 0L) then zmax = 6.0
    
    int = qpint1d('(0.006+0.072*X^1.35)/(1+(X/2.0)^2.4)',zmin,zmax,/expression)

    if keyword_set(doplot) then begin
       im_window, 0, xratio=0.5, /square
       djs_plot, z, sfr, xr=[0,5], yr=[-2.5,-0.5], xsty=3, ysty=3, $
         thick=2.0, xthick=2.0, ythick=2.0, charsize=2.0, charthick=2.0, $
         xtitle='Redshift', ytitle='log \rho_{SFR}(z)'
       djs_oplot, z, alog10(sfr), line=0, thick=2.0
    endif
    
return, int
end
