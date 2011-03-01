pro im_oplot_box, xy, pa, xynew=xynew, _extra=extra

; box coordinates clock-wise from the lower left corner
; xy = [ [min(x),min(y)], [min(x),max(y)], [max(x),max(y)], [max(x),min(y)] ]
    
    pa = 45.0

    xy = fltarr(2,4) ; [x and y coordinate of each corner]
    xy = [ [-1,-1], [-1,1], [1,1], [1,-1] ]

    xynew = im_offset_and_rotate(xy,pa)
    
;   xynew = xy*0.0
;   for i = 0L, 3L do xynew[*,i] = im_offset_and_rotate(reform(xy[*,i]),pa)

    indx = [0,1,2,3,0]
    
    plot, [0], [0], /nodata, xrange=[-2,2], yrange=[-2,2], xsty=3, ysty=3
    for i = 0L, 3L do oplot, xy[*,indx[i]], xy[*,indx[i+1]], line=0
    for i = 0L, 3L do oplot, xynew[*,indx[i]], xynew[*,indx[i+1]], line=2
    
return
end    
