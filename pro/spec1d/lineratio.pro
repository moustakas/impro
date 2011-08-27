;+
; NAME:
;       LINERATIO
;
; PURPOSE:
;       Compute emission-line ratios.
;
; CALLING SEQUENCE:
;
; INPUTS:
;
; OPTIONAL INPUTS:
;       xsyserr - systematic error in X [%]
;       ysyserr - systematic error in Y [%]
;
; KEYWORD PARAMETERS:
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;
; INTERNAL SUPPORT ROUTINES:
;       DOSUMS()
;
; PROCEDURES USED:
;       CMSET_OP(), IM_COMPUTE_ERROR(), SPLOG
;
; COMMENTS:
;       NB:  If NOSUBSCRIPT is set then the galaxies with poorly
;       measured line fluxes will return gibberish.  
;
; TODO:  
;       [1] Properly treat upper limits. 
;
; EXAMPLES:
;       The following call would return x = alog10(OIII_3727/H_ALPHA),
;       y = alog10(OIII_5007/H_BETA) and the associated errors:
;
;          IDL> lineratio, 'OII_3727', 'H_ALPHA', 'OIII_5007', $
;          IDL> 'H_BETA', x, xerr, y, yerr 
;
; MODIFICATION HISTORY:
;       J. Moustakas, 2002 Nov 07, U of A - written
;       jm04nov03uofa - added NOXLOG and NOYLOG keywords
;       jm05aug04uofa - added XSYSERR and YSYSERR optional inputs 
;
; Copyright (C) 2002, 2004-2005, John Moustakas
; 
; This program is free software; you can redistribute it and/or modify 
; it under the terms of the GNU General Public License as published by 
; the Free Software Foundation; either version 2 of the License, or
; (at your option) any later version. 
; 
; This program is distributed in the hope that it will be useful, but 
; WITHOUT ANY WARRANTY; without even the implied warranty of
; MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
; General Public License for more details. 
;-

function dosums, linestruct, tags, linename, sum_err=sum_err, $
  index=index, nindex=nindex, ulimit=ulimit, nogo=nogo, $
  snrcut=snrcut, syserr=syserr

    nspec = n_elements(linestruct)
    
    if (n_elements(syserr) eq 0L) then syserr = 0.0
    
    if strcompress(linename[0],/remove) eq '' then begin
       sum = 0D
       sum_err = -1D
       index = -1L
       nindex = 0L
       nogo = 1L       ; no line information
       return, sum
    endif else nogo = 0L ; line information
    
    nline = n_elements(linename)
    for i = 0L, nline-1L do begin

       matchline = match_string(linename[i],tags,/exact,index=w)
       
       if (i eq 0L) then begin

          newsum = reform((linestruct.(w))[0,*])
          newsum_err = reform((linestruct.(w))[1,*])

          sum = newsum
          sum_err = sqrt(newsum_err^2 + (newsum*syserr/100.0)^2)

          junk = where(sum_err le 0.0,njunk)
          if (njunk ne 0L) then stop
          
          index = where((newsum_err gt 0.0) and finite(newsum) and $
            (finite(newsum_err)) and (newsum/newsum_err gt snrcut),nindex)

          if (nindex eq 0L) then begin
;            message, 'Problem here!'
             index = replicate(-1L,nspec)
             return, sum
          endif

          ulimit = where(newsum_err eq -3.0,nulimit)

       endif else begin

          newsum = reform((linestruct.(w))[0,*])
          newsum_err = reform((linestruct.(w))[1,*])

          sum = sum + newsum
          sum_err = sqrt( sum_err^2.0 + newsum_err^2.0 )

          firstindex = where((newsum_err gt 0.0) and finite(newsum) and $
            (finite(newsum_err)) and (newsum/newsum_err gt snrcut),nfirstindex)

          if (nfirstindex eq 0L) then begin
             index = replicate(-1L,nspec)
             return, sum
          endif else begin
             index = cmset_op(index,'AND',firstindex)
;            if (index[0] eq -1L) then message, 'Problem here!'
          endelse

          ulimiti = where(newsum_err eq -3.0,nulimiti)
          if (nulimit ne 0L) and (nulimiti ne 0l) then ulimit = cmset_op(ulimit,'OR',ulimiti)

       endelse
       
    endfor

return, sum
end

pro lineratio, linestruct, line1x, line2x, line1y, line2y, $
  x, xerr, y, yerr, snrcut=snrcut, index=index, nindex=nindex, $
  nolog=nolog, noxlog=noxlog, noylog=noylog, nosubscript=nosubscript, $
  xsyserr=xsyserr, ysyserr=ysyserr
    
    if n_elements(linestruct) eq 0L then begin
       splog, 'LINESTRUCT not defined.'
       return
    endif

    if keyword_set(nolog) then begin
       noxlog = 1L
       noylog = 1L
    endif

    if n_elements(snrcut) eq 0L then snrcut = 1.0
    tags = tag_names(linestruct)

    x1 = dosums(linestruct,tags,line1x,snrcut=snrcut,sum_err=x1err,$
      index=indxx1,nindex=nindxx1,ulimit=ulimitx1,nogo=nogox1,$
      syserr=xsyserr)

    x2 = dosums(linestruct,tags,line2x,snrcut=snrcut,sum_err=x2err,$
      index=indxx2,nindex=nindxx2,ulimit=ulimitx2,nogo=nogox2,$
      syserr=xsyserr)

    y1 = dosums(linestruct,tags,line1y,snrcut=snrcut,sum_err=y1err,$
      index=indxy1,nindex=nindxy1,ulimit=ulimity1,nogo=nogoy1,$
      syserr=ysyserr)

    y2 = dosums(linestruct,tags,line2y,snrcut=snrcut,sum_err=y2err,$
      index=indxy2,nindex=nindxy2,ulimit=ulimity2,nogo=nogoy2,$
      syserr=ysyserr)

; initialize the INDEX variable; begin by clearing previously defined
; INDEX and NINDEX values    
    
    delvarx, index, nindex
    index = -1L

; ---------------------------------------------------------------------------
; x1, x2, y1, y2 defined
; ---------------------------------------------------------------------------
    
    if (nindxx1 gt 0L) and (nindxx2 gt 0L) and (nindxy1 gt 0L) and (nindxy2 gt 0L) then begin
       index = cmset_op(cmset_op(indxx1,'and',indxx2),'and',cmset_op(indxy1,'and',indxy2))
;      print, 'x1, x2, y1, y2 defined'
    endif

; ---------------------------------------------------------------------------
; x1, x2, y1 defined, no y2
; ---------------------------------------------------------------------------
    
    if (nindxx1 gt 0L) and (nindxx2 gt 0L) and (nindxy1 gt 0L) and (nindxy2 eq 0L) and nogoy2 then begin
       index = cmset_op(cmset_op(indxx1,'and',indxx2),'and',indxy1)
;      print, 'x1, x2, y1 defined, no y2'
    endif

; ---------------------------------------------------------------------------
; x1, x2, y2 defined, no y1
; ---------------------------------------------------------------------------

    if (nindxx1 gt 0L) and (nindxx2 gt 0L) and (nindxy1 eq 0L) and (nindxy2 gt 0L) and nogoy1 then begin
       index = cmset_op(cmset_op(indxx1,'and',indxx2),'and',indxy2)
;      print, 'x1, x2, y2 defined, no y1'
    endif

; ---------------------------------------------------------------------------
; x1, x2 defined, no y1, y2
; ---------------------------------------------------------------------------

    if (nindxx1 gt 0L) and (nindxx2 gt 0L) and (nindxy1 eq 0L) and (nindxy2 eq 0L) and nogoy1 and nogoy2 then begin
       index = cmset_op(indxx1,'and',indxx2)
;      print, 'x1, x2 defined, no y1, y2'
    endif

; ---------------------------------------------------------------------------
; x1, y1, y2 defined, no x2
; ---------------------------------------------------------------------------

    if (nindxx1 gt 0L) and (nindxx2 eq 0L) and (nindxy1 gt 0L) and (nindxy2 gt 0L) and nogox2 then begin
       index = cmset_op(indxx1,'and',cmset_op(indxy1,'and',indxy2))
;      print, 'x1, y1, y2 defined, no x2'
    endif

; ---------------------------------------------------------------------------
; x2, y1, y2 defined, no x1
; ---------------------------------------------------------------------------

    if (nindxx1 eq 0L) and (nindxx2 gt 0L) and (nindxy1 gt 0L) and (nindxy2 gt 0L) and nogox1 then begin
       index = cmset_op(indxx2,'and',cmset_op(indxy1,'and',indxy2))
;      print, 'x2, y1, y2 defined, no x1'
    endif

; ---------------------------------------------------------------------------
; y1, y2 defined, no x1, x2
; ---------------------------------------------------------------------------

    if (nindxx1 eq 0L) and (nindxx2 eq 0L) and (nindxy1 gt 0L) and (nindxy2 gt 0L) and nogox1 and nogox2 then begin
       index = cmset_op(indxy1,'AND',indxy2)
;      print, 'y1, y2 defined, no x1, x2'
    endif

; ---------------------------------------------------------------------------
; x1, y1 defined, no x2, y2
; ---------------------------------------------------------------------------

    if (nindxx1 gt 0L) and (nindxx2 eq 0L) and (nindxy1 gt 0L) and (nindxy2 eq 0L) and nogox2 and nogoy2 then begin
       index = cmset_op(indxx1,'and',indxy1)
;      print, 'x1, y1 defined, no x2, y2'
    endif

; ---------------------------------------------------------------------------
; x1 defined, no x2, y1, y2
; ---------------------------------------------------------------------------

    if (nindxx1 gt 0L) and (nindxx2 eq 0L) and (nindxy1 eq 0L) and (nindxy2 eq 0L) and nogox2 and nogoy1 and nogoy2 then begin
       index = indxx1
;      print, 'x1 defined, no x2, y1, y2'
    endif

; ---------------------------------------------------------------------------
; y1 defined, no x1, x2, y2
; ---------------------------------------------------------------------------

    if (nindxx1 eq 0L) and (nindxx2 eq 0L) and (nindxy1 gt 0L) and (nindxy2 eq 0L) and nogox1 and nogox2 and nogoy2 then begin
       index = indxy1
;      print, 'y1 defined, no x1, x2, y2'
    endif

    junk = where(index ne -1L,nindex)
    if (nindex gt 0L) then if (nindex ne n_elements(index)) then message, 'Houston, we have a problem.'

; if NOSUBSCRIPT=1 then do not index the arrays to the well-measured
; lines; just return the indices 
    
    if not keyword_set(nosubscript) then begin
       
       if (nindxx1 gt 0L) and (nindex ne 0L) then begin
          x1 = x1[index] & x1err = x1err[index]
       endif
       if (nindxx2 gt 0L) and (nindex ne 0L) then begin
          x2 = x2[index] & x2err = x2err[index]
       endif
       if (nindxy1 gt 0L) and (nindex ne 0L) then begin
          y1 = y1[index] & y1err = y1err[index]
       endif
       if (nindxy2 gt 0L) and (nindex ne 0L) then begin
          y2 = y2[index] & y2err = y2err[index]
       endif

    endif

    if (nindxx1 gt 0L) and (nindxx2 gt 0L) and (nindex ne 0L) then begin ; x1 defined, x2 defined

       if keyword_set(noxlog) then begin
          x = x1/x2
          xerr = im_compute_error(x1,x1err,x2,x2err,/quotient)
       endif else begin 
          x = alog10(x1/x2)
          xerr = im_compute_error(x1,x1err,x2,x2err,/log)
       endelse

    endif

    if (nindxx1 gt 0L) and (nindxx2 eq 0L) and (nindex ne 0L) then begin ; x1 defined, no x2

       if keyword_set(noxlog) then begin
          x = x1
          xerr = x1err
       endif else begin 
          x = alog10(x1)
          xerr = im_compute_error(x1,x1err,/log)
       endelse
       
    endif

    if (nindxx1 eq 0L) and (nindxx2 gt 0L) and (nindex ne 0L) then begin ; no x1, x2 defined

       if keyword_set(noxlog) then begin
          x = x2 
          xerr = x2err
       endif else begin
          x = alog10(x2)
          xerr = im_compute_error(x2,x2err,/log)
       endelse

    endif

    if (nindxy1 gt 0L) and (nindxy2 gt 0L) and (nindex ne 0L) then begin ; y1 defined, y2 defined

       if keyword_set(noylog) then begin
          y = y1/y2 
          yerr = im_compute_error(y1,y1err,y2,y2err,/quotient)
       endif else begin
          y = alog10(y1/y2)
          yerr = im_compute_error(y1,y1err,y2,y2err,/log)
       endelse
       
    endif

    if (nindxy1 gt 0L) and (nindxy2 eq 0L) and (nindex ne 0L) then begin ; y1 defined, no y2

       if keyword_set(noylog) then begin
          y = y1 
          yerr = y1err
       endif else begin
          y = alog10(y1)
          yerr = im_compute_error(y1,y1err,/log)
       endelse

    endif

    if (nindxy1 eq 0L) and (nindxy2 gt 0L) and (nindex ne 0L) then begin ; no y1, y2 defined
       
       if keyword_set(noylog) then begin
          y = y2
          yerr = y2err
       endif else begin
          y = alog10(y2)
          yerr = im_compute_error(y2,y2err,/log)
       endelse

    endif

return
end    
