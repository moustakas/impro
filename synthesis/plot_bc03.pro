pro plot_bc03, bc, _extra=extra
; jm04feb26uofa

    window, 0, xs=550, ys=550
    
    nage = n_elements(bc.age)

    xrange = minmax(bc.wave)

    if (n_elements(extra) ne 0L) then begin
       if tag_exist(extra,'XRANGE') then xrange = extra.xrange
    endif
    
    for i = 0L, nage-1L do begin

;      if (i eq 0L) then $
       djs_plot, bc.wave, bc.flux[*,i], xstyle=3, $
         ystyle=3, xthick=2.0, ythick=2.0, charsize=1.5, charthick=2.0, $
         xrange=xrange, _extra=extra

       legend, string(bc.age[i]/1E6,format='(F8.2)')+' Myr', /left, /top, $
         box=0, charsize=1.5, charthick=2.0

       cc = get_kbrd(1)
       
    endfor
    
return
end
    
