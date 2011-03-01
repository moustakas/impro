function dosums, linestruct, tags, linename, sum_err=sum_err, $
  index=index, ulimit=ulimit

    if strcompress(linename[0],/remove) eq '' then begin
       sum = 0.0
       sum_err = -1.0
       index = -1L
       return, sum
    endif
    
    nline = n_elements(linename)
    for i = 0L, nline-1L do begin

       matchline = match_string(linename[i],tags,/exact,index=w)
       
       if i eq 0L then begin
          sum = reform((linestruct.(w))[0,*])
          sum_err = reform((linestruct.(w))[1,*])
          index = where((sum_err gt 0.0) and finite(sum) and (finite(sum_err)))
          ulimit = where(sum_err eq -3.0,nulimit)
       endif else begin
          sum = sum + reform((linestruct.(w))[0,*])
          sum_err = sqrt( sum_err^2.0 + (reform((linestruct.(w))[1,*]))^2.0 )
          index = cmset_op(index,'AND',where((sum_err gt 0.0) and finite(sum) and (finite(sum_err))))
          ulimiti = where(sum_err eq -3.0,nulimiti)
          if (nulimit ne 0L) and (nulimiti ne 0l) then ulimit = cmset_op(ulimit,'OR',ulimiti)
       endelse

    endfor

return, sum
end

pro ilineratio, linefit, x, xerr, y, yerr, linex1=linex1, linex2=linex2, $
  liney1=liney1, liney2=liney2, xlog=xlog, ylog=ylog, cut=cut, newlinefit=newlinefit

    nspec = n_elements(linefit)
    if (nspec eq 0L) then begin
       print, 'Syntax - '
       return
    endif
    
    tags = tag_names(linefit)
    nx1 = n_elements(linex1)
    nx2 = n_elements(linex2)
    ny1 = n_elements(liney1)
    ny2 = n_elements(liney2)
    ncut = n_elements(cut)
    
    if (ncut ne 0L) then if (cut[0] eq -1L) then linecut = lindgen(nspec) else linecut = cut

; no X line ratio    
    
    if (nx1 eq 0L) and (nx2 eq 0L) then begin
       x = 0.0
       xerr = -1.0
    endif

; no Y line ratio    
    
    if (ny1 eq 0L) and (ny2 eq 0L) then begin
       y = 0.0
       yerr = -1.0
    endif

    newlinefit = linefit[linecut]
    for k = 0L, n_elements(newlinefit)-1L do begin

       if nx1 
       


    endfor
    
stop    
    
    x1 = dosums(linefit,tags,line1x,sum_err=x1err,index=indxx1,ulimit=ulimitx1)
    x2 = dosums(linefit,tags,line2x,sum_err=x2err,index=indxx2,ulimit=ulimitx2)
    y1 = dosums(linefit,tags,line1y,sum_err=y1err,index=indxy1,ulimit=ulimity1)
    y2 = dosums(linefit,tags,line2y,sum_err=y2err,index=indxy2,ulimit=ulimity2)

; ---------------------------------------------------------------------------
; x1, x2, y1, y2 defined
; ---------------------------------------------------------------------------
    
    if (indxx1[0] ne -1L) and (indxx2[0] ne -1L) and (indxy1[0] ne -1L) and (indxy2[0] ne -1L) then begin
       index = cmset_op(cmset_op(indxx1,'and',indxx2),'and',cmset_op(indxy1,'and',indxy2))
    endif

; ---------------------------------------------------------------------------
; x1, x2, y1 defined, no y2
; ---------------------------------------------------------------------------
    
    if (indxx1[0] ne -1L) and (indxx2[0] ne -1L) and (indxy1[0] ne -1L) and (indxy2[0] eq -1L) then begin
       index = cmset_op(cmset_op(indxx1,'and',indxx2),'and',indxy1)
    endif

; ---------------------------------------------------------------------------
; x1, x2, y2 defined, no y1
; ---------------------------------------------------------------------------

    if (indxx1[0] ne -1L) and (indxx2[0] ne -1L) and (indxy1[0] eq -1L) and (indxy2[0] ne -1L) then begin
       index = cmset_op(cmset_op(indxx1,'and',indxx2),'and',indxy2)
    endif

; ---------------------------------------------------------------------------
; x1, x2 defined, no y1, y2
; ---------------------------------------------------------------------------

    if (indxx1[0] ne -1L) and (indxx2[0] ne -1L) and (indxy1[0] eq -1L) and (indxy2[0] eq -1L) then begin
       index = cmset_op(indxx1,'and',indxx2)
    endif

; ---------------------------------------------------------------------------
; x1, y1, y2 defined, no x2
; ---------------------------------------------------------------------------

    if (indxx1[0] ne -1L) and (indxx2[0] eq -1L) and (indxy1[0] ne -1L) and (indxy2[0] ne -1L) then begin
       index = cmset_op(indxx1,'and',cmset_op(indxy1,'and',indxy2))
    endif

; ---------------------------------------------------------------------------
; x2, y1, y2 defined, no x1
; ---------------------------------------------------------------------------

    if (indxx1[0] eq -1L) and (indxx2[0] ne -1L) and (indxy1[0] ne -1L) and (indxy2[0] ne -1L) then begin
       index = cmset_op(indxx2,'and',cmset_op(indxy1,'and',indxy2))
    endif

; ---------------------------------------------------------------------------
; y1, y2 defined, no x1, x2
; ---------------------------------------------------------------------------

    if (indxx1[0] eq -1L) and (indxx2[0] eq -1L) and $
      (indxy1[0] ne -1L) and (indxy2[0] ne -1L) then begin
       
       index = cmset_op(indxy1,'AND',indxy2)

;       if (ulimity1[0] ne -1L) and (ulimity2[0] ne -1L) then begin
;
;          rem = cmset_op(ulimity1,'AND',ulimity2,/index,count=rcount) ; unable to form upper limit
;          if (rcount eq n_elements(ulimity1)) then ulimity1 = -1L else $
;            if (rem[0] ne -1L) then remove, rem, ulimity1
;
;       endif
       
    endif

; ---------------------------------------------------------------------------
; x1, y1 defined, no x2, y2
; ---------------------------------------------------------------------------

    if (indxx1[0] ne -1L) and (indxx2[0] eq -1L) and (indxy1[0] ne -1L) and (indxy2[0] eq -1L) then begin
       index = cmset_op(indxx1,'and',indxy1)
    endif

; ---------------------------------------------------------------------------
; x1 defined, no x2, y1, y2
; ---------------------------------------------------------------------------

    if (indxx1[0] ne -1L) and (indxx2[0] eq -1L) and (indxy1[0] eq -1L) and (indxy2[0] eq -1L) then begin
       index = indxx1
    endif

; ---------------------------------------------------------------------------
; y1 defined, no x1, x2, y2
; ---------------------------------------------------------------------------

    if (indxx1[0] eq -1L) and (indxx2[0] eq -1L) and (indxy1[0] ne -1L) and (indxy2[0] eq -1L) then begin
       index = indxy1
    endif

    nindex = n_elements(index)

; do not index the arrays to the well-measured lines; just return the
; indices 
    
    if not keyword_set(nosubscript) then begin
       
       if (indxx1[0] ne -1L) then begin
          x1 = x1[index] & x1err = x1err[index]
       endif
       if (indxx2[0] ne -1L) then begin
          x2 = x2[index] & x2err = x2err[index]
       endif
       if (indxy1[0] ne -1L) then begin
          y1 = y1[index] & y1err = y1err[index]
       endif
       if (indxy2[0] ne -1L) then begin
          y2 = y2[index] & y2err = y2err[index]
       endif

    endif

    if (indxx1[0] ne -1L) and (indxx2[0] ne -1L) then begin ; x1 defined, x2 defined

       if keyword_set(nolog) then begin
          x = x1/x2
          xerr = im_compute_error(x1,x1err,x2,x2err,/quotient)
       endif else begin 
          x = alog10(x1/x2)
          xerr = im_compute_error(x1,x1err,x2,x2err,/log)
       endelse

    endif

    if (indxx1[0] ne -1L) and (indxx2[0] eq -1L) then begin ; x1 defined, no x2

       if keyword_set(nolog) then begin
          x = x1
          xerr = x1err
       endif else begin 
          x = alog10(x1)
          xerr = im_compute_error(x1,x1err,/log)
       endelse
       
    endif

    if (indxx1[0] eq -1L) and (indxx2[0] ne -1L) then begin ; no x1, x2 defined

       if keyword_set(nolog) then begin
          x = x2 
          xerr = x2err
       endif else begin
          x = alog10(x2)
          xerr = im_compute_error(x2,x2err,/log)
       endelse

    endif

    if (indxy1[0] ne -1L) and (indxy2[0] ne -1L) then begin ; y1 defined, y2 defined

       if keyword_set(nolog) then begin
          y = y1/y2 
          yerr = im_compute_error(y1,y1err,y2,y2err,/quotient)
       endif else begin
          y = alog10(y1/y2)
          yerr = im_compute_error(y1,y1err,y2,y2err,/log)
       endelse
       
    endif

    if (indxy1[0] ne -1L) and (indxy2[0] eq -1L) then begin ; y1 defined, no y2

       if keyword_set(nolog) then begin
          y = y1 
          yerr = y1err
       endif else begin
          y = alog10(y1)
          yerr = im_compute_error(y1,y1err,/log)
       endelse

    endif

    if (indxy1[0] eq -1L) and (indxy2[0] ne -1L) then begin ; no y1, y2 defined
       
       if keyword_set(nolog) then begin
          y = y2
          yerr = y2err
       endif else begin
          y = alog10(y2)
          yerr = im_compute_error(y2,y2err,/log)
       endelse

    endif
    
return
end    
