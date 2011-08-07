;+
; NAME:
;	DKBOXSTATS()
;
; PURPOSE:
;	To take a one or two dimensional array and return each pixel
;	by either the minimum, maximum, mean, median, total (sum), or
;	product of the pixels in a rectangular box of size
;	[XWIDTH,YWIDTH] centered on that pixel.
;
; CALLING SEQUENCE:
;	outarray = dkboxstats(array,[xwidth=],[ywidth=],[boxstat=])
;
; INPUTS:
;	array   - one or two dimensional array
;
; OPTIONAL INPUTS:
;	xwidth  - 
;	ywidth  - 
;	boxstat - 
;       extra   - rejection parameter keywords for DJSIG() 
;	
; KEYWORD PARAMETERS:
;
; OUTPUTS:
;	outarray - output array
;
; COMMENTS:
;	Edges are not handled correctly (add a MISSING keyword?).
;
;       Assume that when computing the median statistics that the
;       array has been passed as sorted.
;
; EXAMPLE:
;
; PROCEDURES USED:
;	MIN(), MAX(), DJS_MEAN(), DJS_MEDIAN(), TOTAL(), PRODUCT(),
;	DJSIG()
;
; MODIFICATION HISTORY:
;	J. Moustakas, 2001 October 29, U of A
;       jm02nov15uofa - added TOTAL and PRODUCT statistics
;       jm04jan22uofa - added SIGMA statistics
;       jm04jul13uofa - when computing the median in one dimension
;                       also return the lower and upper one-sigma 
;
; Copyright (C) 2001-2002, 2004, John Moustakas
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

function dkboxstats, array, xwidth=xwidth, ywidth=ywidth, boxstat=boxstat, $
  lower=lower, upper=upper, _extra=extra

    if n_elements(array) eq 0L then begin
       print, 'Syntax - outarray = dkboxstats(array,xwidth=,ywidth=,$'
       print, '   boxstat=,lower=,upper=,_extra=extra)'
       return, -1
    endif
    
    dim = size(array,/dimension)
    ndim = n_elements(dim)

    if n_elements(boxstat) ne 0L then begin

       boxstat = strlowcase(strcompress(boxstat,/remove))
       
       case boxstat of
          'min':
          'max':
          'mean':
          'median':
          'total':
          'product':
          'sigma':
          else: begin
             print, 'Unsupported statistic:  select from MIN, MAX, MEAN, MEDIAN, TOTAL, PRODUCT, or SIGMA.'
             return, -1
          end
       endcase
    
    endif else boxstat = 'min'

    if n_elements(xwidth) eq 0L then xwidth = 3L ; default x-width
    wxhalf = fix(xwidth/2)                       ; half-width

    onesigma = errorf(1/sqrt(2.0))-0.5

    outarray = array*0.0 ; output array
    lower = array*0.0
    upper = array*0.0
       
    case ndim of

       1L: begin
          dim = dim[0]
          for i = 0L, dim-1L do begin
             subarr = array[(i-wxhalf)>0L:(i+wxhalf)<(dim-1L)]
             nsub = n_elements(subarr)
             case boxstat of
                'min': outarray[i] = min(subarr)
                'max': outarray[i] = max(subarr)
                'mean': outarray[i] = djs_mean(subarr)
                'total': outarray[i] = total(subarr)
                'product': outarray[i] = product(subarr)
                'sigma': outarray[i] = djsig(subarr,_extra=extra)
                'median': begin
                   srt = sort(subarr)
                   outarray[i] = interpol(subarr[srt],findgen(nsub)/nsub,0.5)
                   lower[i] = outarray[i] - interpol(subarr[srt],findgen(nsub)/nsub,0.5-onesigma)
                   upper[i] = interpol(subarr[srt],findgen(nsub)/nsub,0.5+onesigma) - outarray[i]
;                  outarray[i] = djs_median(subarr)
                end
             endcase
          endfor 
       end 

       2L: begin

          if n_elements(ywidth) eq 0L then ywidth = 3L ; default y-width
          wyhalf = fix(ywidth/2)

          dimx = dim[0]
          dimy = dim[1]
          for i = 0L, dimx-1L do for j = 0L, dimy-1L do begin
             subarr = array[(i-wxhalf)>0L:(i+wxhalf)<(dimx-1),(j-wyhalf)>0L:(j+wyhalf)<(dimy-1)]
             case boxstat of
                'min': outarray[i,j] = subarr[(where(subarr eq min(subarr)))[0]]
                'max': outarray[i,j] = subarr[(where(subarr eq max(subarr)))[0]]
                'mean': outarray[i,j] = djs_mean(subarr)
                'median': outarray[i,j] = djs_median(subarr)
                'total': outarray[i,j] = total(subarr)
                'product': outarray[i,j] = product(subarr)
                'sigma': outarray[i,j] = djsig(subarr,_extra=extra)
             endcase
          endfor 
       end 

       else: begin
          print, 'Unsupported number of dimensions!'
          return, -1
       end
    endcase 
       
return, outarray
end
