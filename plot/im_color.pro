;+
; NAME:
;   IM_COLOR()
; PURPOSE:
;   Trivial wrapper on CGCOLOR() to access all the colors in 'rgb.txt'.  
; INPUTS: 
;   colorname - color name (see ${IMPRO_DIR}/etc/rgb.txt)
; OPTIONAL INPUTS: 
;   colorindx - color table index number (default 100)
;   extra - additional inputs to CGCOLOR()
; MODIFICATION HISTORY:
;   J. Moustakas, 2011 Aug 16, UCSD
;-

function im_color, colorname, colorindx, _extra=extra
    ncolor = n_elements(colorname)
    if (ncolor eq 0) then begin
       doc_library, 'im_color'
       return, -1
    endif
    if (n_elements(colorindx) eq 0) then colorindx = 100
    if (ncolor gt 1) then begin
       col = lonarr(ncolor)
       for ii = 0, ncolor-1 do col[ii] = $
         im_color(colorname[ii],colorindx+ii)
       return, col
    endif
    if (strcompress(colorname,/rem) ne '') then cname = colorname else begin
       if (!d.name eq 'X') then cname = 'white' else cname = 'black'
    endelse
return, cgcolor(cname,file=getenv('IMPRO_DIR')+'/etc/rgb.txt',colorindx,_extra=extra)
end

    
