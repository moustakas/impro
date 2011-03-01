pro refraction_plot, postscript=postscript
; jm02feb15uofa 
; generate a plot of atmospheric refraction with respect to 5000 
; Angstrom as a function of wavelength and zenith angle
; references: Filippenko 1982, PASP, 94, 715

; zenith angle vector
    
    seczmin = 1.0
    seczmax = 5.0
    secdza = 0.5
    secz = findgen((seczmax-seczmin)/secdza+1)*secdza+seczmin
    nza = n_elements(secz)
    
    refwave = 5000.0

    refract = compute_refraction(secz,refwave=refwave,lambda=lambda)

    if keyword_set(postscript) then begin
       dfpsplot, 'refraction_plot.ps', /landscape
    endif else window, 0, xs=550, ys=550
    
    yrange = minmax(refract)
    plot, lambda, refract[*,0], xsty=3, ysty=3, yrange=yrange, charsize=2.0, $
      charthick=2.0, xtitle='Wavelength (\AA)', ytitle='Refraction (arcsec)', $
      title='Atmospheric Refraction', xthick=2.0, ythick=2.0
    legend, ['H = 2 km','T = 7 C','P = 600 mm Hg','f = 8 mm Hg'], /right, /top, $
      box=0, charsize=2.0, charthick=2.0
    for j = 0L, nza-1L do oplot, lambda, refract[*,j]

    if keyword_set(postscript) then dfpsclose

return
end
