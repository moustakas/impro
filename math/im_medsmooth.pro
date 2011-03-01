;+
; NAME:
;    IM_MEDSMOOTH()
;
; PURPOSE:
;    Median smoothing of a vector, including points near its ends
;    with optional iterative sigma-clipping.
;
; INPUTS:
;    VECTOR  = The (1-d numeric) vector to be smoothed
;    WINDOW = Odd integer giving the full width of the window over which 
;         the median is determined for each point.     (If WINDOW is
;         specified as an even number, then the effect is the same as
;         using WINDOW+1)   
;    clipsig - reject pixels that deviate more than CLIPSIG from the
;      *median* (default 3.0)
;    maxiter - maximum number of clipping iterations (default 10)
;
; OUTPUT:
;    Function returns the smoothed vector
;
; PROCEDURE:
;    Each point is replaced by the median of the nearest WINDOW of points.
;    The width of the window shrinks towards the ends of the vector, so that
;    only the first and last points are not filtered. These points are 
;    replaced by forecasting from smoothed interior points.
;
; EXAMPLE:
;    Create a vector with isolated high points near its ends
;    IDL> a = randomn(seed,40) & a[1] = 10  & a[38] = 10
;    Now do median smoothing with a 7 point window 
;    IDL> b = medsmooth(a,7)
;    Note that, unlike MEDIAN(), that MEDSMOOTH will remove the isolated
;    high points near the ends.
;
; REVISION HISTORY:
;    Written, H. Freudenreich, STX, 12/89
;    H.Freudenreich, 8/90: took care of end-points by shrinking window.
;    Speed up using vector median when possible  W. Landsman
;    February 2002
;    J. Moustakas, 2007 Aug 30, NYU - optional weights
;    J. Moustakas, 2008 Aug 14, NYU - added optional sigma clipping 
;-

function im_medsmooth, array, window, weight=weight, $
  clip_outliers=clip_outliers

    lend = n_elements(array)-1L
    if (lend+1) lt window then begin
       message,/con, $
         'error - size of smoothing window must be smaller than array size'
       return,array
    endif

    if arg_present(weight) then begin
       if (n_elements(weight) ne n_elements(array)) then begin
          print, 'dimensions of array and weight must agree.'
          return, array
       endif
    endif else weight = array*0.0+1.0 ; equal weights

    if (n_elements(clipsig) eq 0L) then clipsig = 3.0
    if (n_elements(maxiter) eq 0L) then maxiter = 10
    
    smoothed = array*0.0
;   smoothed = median(array,window)
    if keyword_set(clip_outliers) then begin

stop       
       
       smoothed = im_quantile(array,weight,quant=0.5)
    endif else begin
       for jj = 0L, lend do smoothed[jj] = im_quantile($
         array[jj:jj+window-1L],weight[jj:jj+window-1L],quant=0.5)
    endelse


; fix the ends; don't bother clipping

    numloop = (window-1L)/2L - 1L
    if (numloop gt 0L) then begin
       for j = 1L, numloop do begin 
          len = 2L*j+1L
          smoothed[j] = im_quantile(array[0l:len-1l],weight[0l:len-1l],quant=0.5)
          smoothed[lend-j] = im_quantile(array[lend-len+1l:lend],weight[lend-len+1l:lend],quant=0.5)
       endfor
    endif

; now replace the very last and first points

    y0 = 3.*array[0]-2.*array[1] ; predicted value of point -1
    smoothed[0] = median([y0,array[0],array[1]])
    y0 = 3.*array[lend]-2.*array[lend-1] ; predicted value of point lend+1
    smoothed[lend] = median([y0,array[lend],array[lend-1]])
    
return, smoothed
end
