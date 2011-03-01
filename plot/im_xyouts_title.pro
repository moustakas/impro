pro im_xyouts_title, xtitle=xtitle, ytitle=ytitle, charsize=charsize, $
  xspacing=xspacing, yspacing=yspacing, xlog=xlog, ylog=ylog, _extra=extra
; jm05apr21uofa
; jm07mar15nyu - added LOG keyword    

    if (n_elements(charsize) eq 0L) then charsize = 1.0
    if (n_elements(xspacing) eq 0L) then xspacing = 8.0
    if (n_elements(yspacing) eq 0L) then yspacing = 3.0
    
    if keyword_set(xlog) then xrange = 10^!x.crange else xrange = !x.crange
    if keyword_set(ylog) then yrange = 10^!y.crange else yrange = !y.crange

    if (n_elements(xtitle) ne 0L) then begin

       xypos = convert_coord(mean(xrange),yrange[1],/data,/to_normal)
       xpos = reform(xypos[0,*])
       ypos = reform(xypos[1,*]) + (!d.y_ch_size/float(!d.y_size))*(yspacing>charsize)

       xyouts, xpos, ypos, textoidl(xtitle), /normal, align=0.5, $
         charsize=charsize, _extra=extra

    endif

    if (n_elements(ytitle) ne 0L) then begin

       xypos = convert_coord(xrange[1],mean(yrange),/data,/to_normal)
       xpos = reform(xypos[0,*]) + (!d.x_ch_size/float(!d.x_size))*(xspacing>charsize)
       ypos = reform(xypos[1,*]) 

       xyouts, xpos, ypos, textoidl(ytitle), /normal, align=0.5, $
         orientation=90, charsize=charsize, _extra=extra

    endif

return
end
    
    
