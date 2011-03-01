;+
; NAME: ARM_HISTOMAKER
;       
; CATEGORY: plotting
;
; PURPOSE: generate publication-quality histograms
;
; CALLING SEQUENCE:
;
; INPUTS: arr1 - array of values to create histogram of
;       
; OPTIONAL INPUTS:
;    arr[2-10] - additional arrays of values
;    color     - string array of line color names (DJS_OPLOT colors)
;    fillcolor - string array of fill color names (DJS_OPLOT colors)
;    linestyle - array of linestyles (default=[0,2,1,3,4...])
;    thick     - array of line thicknesses (default=!p.thick)
;    yscale    - array of values to scale histograms by
;    min       - lower limit of binning
;    max       - upper limit of binning
;    binsize   - width of bins
;    nbins     - number of bins (ignored if BINSIZE is defined)
;    ignore    - numerical value to ignore when binning data
;
; KEYWORDS:
;   cumulative - plots ascending cumulative distribution if > 0
;                plots descending cumulative distribution if < 0
;   noplot     - do not generate plot (meant to be used with HISTODATA)
;   oplot      - overplot histogram(s) onto existing plot
;   percent    - plot y-axis as fraction of distribution
;   grid       - overplots horizontal lines "above" histogram
;   compress   - remove bins with no entries
;   quiet      - suppress printed statistics
;   no_offset  - do not offset histograms from each other (primarily
;                useful when offsets are large due to many histograms)
; 
; OPTIONAL OUTPUTS: 
;   histodata  - structure containing x & y data for histograms
;
; EXAMPLE:
;
; PROCEDURES USED: ARM_ORDEROFMAG, ARM_MULTIPLE, ARM_STRTRIM,
;            DJS_PLOT, DJS_OPLOT 
;
; COMMENTS: accepts all standard parameters accepted by HISTOGRAM and
;           PLOT, FILLCOLOR does not work well when keyword COMPRESS
;           is set, x-axis labels overlap in some cases when keyword
;           COMPRESS is set, GRID lines not designed for when PERCENT
;           keyword is set
;
;           please report bugs to amarble@as.arizona.edu
;
; MODIFICATION HISTORY:
;    written by A.R.Marble, Steward Obs., November 2000
;    streamlined, ARM, July 2001
;    additional features added, ARM, August 2001
;    complete overhaul, ARM, December 2, 2003
;    IGNORE parameter added, ARM, December 4, 2003
;    increased precision of MIN/MAX calculation, ARM, December 10, 2003
;    changed default line thickness to !p.thick, ARM, January 3, 2004
;-

pro arm_histomaker, arr1, arr2, arr3, arr4, arr5, arr6, arr7, arr8, arr9, arr10, $
       binsize=binsize, nbins=nbins, color=color, yscale=yscale, linestyle=linestyle, $
       cumulative=cumulative, min=min, max=max, histodata=histodata, thick=thick, $
       COMPRESS=compress, PERCENT=percent, QUIET=quiet, GRID=grid, NOPLOT=noplot, $
       NO_OFFSET=no_offset, _EXTRA=extra, fillcolor=fillcolor, ignore=ignore, $
       OPLOT=oplot

    n = N_PARAMS()              ; number of datasets to bin
    
; error checking

    if n lt 1 then begin
       PRINT, 'SYNTAX: ARM_HISTOMAKER, arr1, [arr2, arr3, arr4, arr5, ... arr10], '
       PRINT, '            [binsize=, min=, max=, color=, linestyle=, fillcolor=, '
       PRINT, '             yscale=, thick=, histodata=, cumulative=], '
       PRINT, '             [/percent, /grid, /no_offset, /no_plot, /oplot, /quiet, /compress]'
       return
    endif
    
; define defaults
    
    nyscale = N_ELEMENTS(yscale)
    if nyscale ne n then begin
       if nyscale eq 0 then yscale = REPLICATE(1L, n) $
       else if nyscale eq 1 then yscale = REPLICATE(yscale, n) $
       else if nyscale gt n then yscale = yscale[0:n-1] $
       else yscale = [yscale, REPLICATE(1L, n-nyscale)]
    endif
    
    nthick = N_ELEMENTS(thick)
    if nthick ne n then begin
       if nthick eq 0 then thick = REPLICATE(!p.thick, n) $
       else if nthick eq 1 then thick = REPLICATE(thick, n) $
       else if nthick gt n then thick = thick[0:n-1] $
       else thick = [thick, REPLICATE(!p.thick, n-nthick)]
    endif
    
    nstyle = N_ELEMENTS(linestyle)
    styles = [0, 2, 1, 3, 4, 5, 6, 0, 0, 0]
    if nstyle ne n then begin
       if nstyle eq 0 then linestyle = styles[0:n-1] $
       else if nstyle eq 1 then linestyle = REPLICATE(linestyle, n) $
       else if nstyle gt n then linestyle = linestyle[0:n-1] $
       else linestyle = [linestyle, REPLICATE(0L, n-nstyle)]
    endif

    if N_ELEMENTS(cumulative) eq 0L then cumulative = 0L else percent = 1L
    if KEYWORD_SET(percent) then ytitle = 'Fraction' else ytitle = 'Number'
    
    ncolor = N_ELEMENTS(color)
    if ncolor ne n then begin
       if ncolor eq 0 then color = REPLICATE('', n) $
       else if ncolor eq 1 then color = REPLICATE(color, n) $
       else if ncolor gt n then color = color[0:n-1] $
       else color = [color, REPLICATE('', n-ncolor)]
    endif

; remove unwanted values

    if KEYWORD_SET(ignore) then begin
       for i=1,n do begin
          dummy = EXECUTE('arr = arr'+STRN(FIX(i)))
          good = WHERE(arr ne ignore, count)
          if count eq 0L then MESSAGE, 'Array '+STRN(FIX(i))+' has no valid entries!' else $
            dummy = EXECUTE('arr'+STRN(FIX(i))+' = arr'+STRN(FIX(i))+'[good]')
       endfor
    endif

; determine binning range
    
    ranges = MINMAX(arr1)
    for i=2,n do dummy = EXECUTE('ranges = [ranges, minmax(arr'+strn(fix(i))+')]')
    if MIN(ranges) ne MAX(ranges) then begin
       power = ARM_ORDEROFMAG([MIN(ranges), MAX(ranges)])
       if N_ELEMENTS(min) eq 0L then min = $
         (FLOOR(MIN(ranges)/(1d1^(power[0]-1))))*(1d1^(power[0]-1))
       if N_ELEMENTS(max) eq 0L then max = $
         (CEIL(MAX(ranges)/(1d1^(power[1]-1))))*(1d1^(power[1]-1))
    endif else begin
       binsize = 1d0
       min = ranges[0] - 100d0
       max = ranges[0] + 100d0
       xticks = 2
       xminor = -1
       xtickname = [' ', ARM_STRTRIM(ARM_STRTRIM(STRN(ranges[0]), '0'), '.'), ' ']
    endelse

    if N_ELEMENTS(binsize) eq 0L then begin
       if N_ELEMENTS(nbins) eq 0L then nbins = 10. else nbins = FLOAT(nbins)
       binsize = (max-min) / nbins
    endif

; create histogram(s)
    
    histo = HISTOGRAM(DOUBLE([arr1, arr1]), binsize=binsize, min=min, max=max)
    num = N_ELEMENTS(histo)
    totals = REPLICATE(DOUBLE(N_ELEMENTS(arr1)), num)
    offset = REPLICATE(0L, num)
    rr = DINDGEN(num) * binsize + min + binsize / 2d0
    if N_ELEMENTS(rr) eq 1 then rr = [rr, rr + binsize]
    
    for i=1,n do begin
       
       dummy = EXECUTE('arr = arr'+STRN(FIX(i)))
       if not KEYWORD_SET(quiet) then begin
          print
          print,'Array '+strn(fix(i))+': [Min, Max] = ['+STRN(MIN(arr))+', '+STRN(MAX(arr))+']'
          print,'          Mean      =  ' + STRN(TOTAL(arr) / N_ELEMENTS(arr))
          print,'          Sigma     =  ' + STRN(SQRT(VARIANCE(arr)))
       endif
       if i gt 1L then begin
          histo  = [[histo], [HISTOGRAM(DOUBLE([arr, arr]), binsize=binsize, min=min, max=max)]]
          totals = [[totals], [REPLICATE(DOUBLE(N_ELEMENTS(arr)), num)]]
          offset = [[offset], [REPLICATE(-1L, num)]]          
       endif
       
    endfor
    
    histo = histo / 2L
    
; convert to cumulative distributions if desired
    
    if cumulative ne 0 then begin
       cdf = DINDGEN(num, n)
       if cumulative gt 0 then begin 
          for k=0,n-1 do for j=0,num-1 do cdf[j,k] = $
            TOTAL(histo[0:j,k]) / TOTAL(histo[*,k])
       endif else begin
          for k=0,n-1 do for j=0,num-1 do cdf[j,k] = $
            TOTAL(histo[j:num-1,k]) / TOTAL(histo[*,k])
       endelse
       histo = cdf
       
; normalize area under histograms if desired
       
    endif else if KEYWORD_SET(percent) then histo = histo / totals
    
    for i=0,n-1 do histo[*,i] = histo[*,i] * yscale[i]
    
; compress histogram bins if desired
    
    if KEYWORD_SET(compress) then begin
       
       for i=0,n-1 do begin
          wh = where(histo[*,i] ne 0L)
          if i eq 0 then nonzero = wh else nonzero = [nonzero, wh]
       endfor
       nonzero = nonzero[SORT(nonzero)]
       nonzero = nonzero[UNIQ(nonzero)]
       rrtxt = STRING(rr[nonzero], f='(e9.2)')
       xtickname = STRARR(N_ELEMENTS(rrtxt)*2+3) + ' '
       if N_ELEMENTS(xtickname) gt 60 then begin
          MESSAGE, 'too many bins to compress histogram!', /continue
          xtickname = ' '
       endif else begin       
          xtickname[INDGEN(N_ELEMENTS(rrtxt))*2+2] = rrtxt
          xtickinterval = 0.5
          xminor = 1
          num = N_ELEMENTS(nonzero)
          rr = DINDGEN(num) + 0.5
          histo = histo[nonzero, *]
       endelse

    endif
    
; create structure containing histogram data
    
    text = "'x', rr"
    for i=0,n-1 do text = text+", 'y"+STRN(FIX(i+1))+"', histo[*,"+STRN(FIX(i))+"]"
    dummy = EXECUTE("histodata = CREATE_STRUCT("+text+")")
    
    if not KEYWORD_SET(noplot) then begin
       
;;; determine best y-axis tick values
;;       
;;       if KEYWORD_SET(percent) then ymax = 1.1 * MAX(histo) < 1 $
;;       else ymax = 1.1 * MAX(histo)
;;       pwr = ARM_ORDEROFMAG(ymax)    
;;       if ABS(pwr) gt 3 then begin
;;          histo = histo / (10.^pwr)
;;          ymax = CEIL(ymax / (10.^pwr))
;;          ytitle = ytitle + TEXTOIDL(' (x10^{'+STRN(pwr)+'})')
;;       endif
;;       if not KEYWORD_SET(percent) then ymax = CEIL(ymax)
;;       if ymax gt 40 then begin
;;          yints = [5d,10d,20d,25d]
;;          pwr = ARM_ORDEROFMAG(ymax)
;;          if pwr gt 1L then begin
;;             yints = [yints * 10^(pwr-2), yints * 10^(pwr-1), yints * 10^pwr]
;;          endif else yints = [1d, 2d, 4d, yints * 10^(pwr-1), yints * 10^pwr]
;;       endif else yints = [1d, 2d, 4d]
;;       whmin = WHERE(ABS(ymax/yints-10.) eq MIN(ABS(ymax/yints-10.)))
;;       yint=FIX(yints[whmin[0]])
;;       if not KEYWORD_SET(percent) then begin
;;          multcheck = ARM_MULTIPLE(yint, ymax)
;;          while multcheck eq 0 do begin
;;             ymax = ymax + 1
;;             multcheck = ARM_MULTIPLE(yint, ymax)
;;          endwhile
;;       endif
;;       yticks = ymax / yint
;;       if KEYWORD_SET(percent) then if yticks gt 1L then yticks = 1L
;;       if yint lt 5 then yminor = yint else yminor = 5

       ymax = MAX(histo) * 1.1

; create plot box
       
       if not KEYWORD_SET(oplot) then begin
          PLOT, [0, 0], /nodata, ytitle=ytitle, _EXTRA=extra, $
            xrange=[rr[0]-(rr[1]-rr[0]), rr[num-1]+rr[1]-rr[0]], xstyle=1, $
            xtickinterval=xtickinterval, xtickname=xtickname, xminor=xminor, $
            yrange=[0,ymax], $
            ystyle=1, yticks=yticks, yminor=yminor, xticks=xticks, $
            ytick_get=ytickvals
       endif

; get offset values for plotting multiple histograms
       
       normalcoords = CONVERT_COORD([[0.,0.],[1.,1.]], /normal, /to_data)
       xl = normalcoords[0,0]
       xr = normalcoords[0,1]
       
       devicecoords = CONVERT_COORD([[0,0], [1,1]], /device, /to_data)
       xpxl = devicecoords[0,1]-devicecoords[0,0]
       ypxl = devicecoords[1,1]-devicecoords[1,0]
       
       if not KEYWORD_SET(no_offset) then histo = histo + offset * ypxl
       
; overplot horizontal lines if grid keyword is set

       if KEYWORD_SET(grid) then begin
          if (SIZE(histo))[0] gt 1L then maxs = MAX(histo, dimension=2) else maxs = histo
          for i=0,N_ELEMENTS(ytickvals)-1 do begin 
             whless = where(maxs lt ytickvals[i], count)
             if count gt 0L then begin
                DJS_OPLOT,[xl, rr[0]-binsize/2.], ytickvals[[i,i]], thick=1.0, color='gray'
                DJS_OPLOT,[rr[num-1]+binsize/2., xr], ytickvals[[i,i]], thick=1.0, color='gray'
                for j=0,count-1 do DJS_OPLOT, thick=1.0, color='gray', $
                  [rr[whless[j]]-binsize/2.,rr[whless[j]]+binsize/2.], ytickvals[[i,i]]
             endif
          endfor
       endif
       
; plot histograms

       for j=0,n-1 do begin
          
          if not KEYWORD_SET(no_offset) then $
            r = rr + offset[0,j] * xpxl $ ;OFFSET HISTOGRAMS HORIZONTALLY BY ONE PIXEL
          else r = rr
          
          bnpx = binsize / xpxl / 2.0

          if N_ELEMENTS(fillcolor) gt j then for l=-bnpx,bnpx do for k=0,N_ELEMENTS(r)-1 do $
            DJS_OPLOT,r[k]*[1.,1.]+l*xpxl,[0,histo[k,j]],thick=1., color=fillcolor[j]
          
          DJS_OPLOT, r, histo[*,j], psym=10, color=color[j], linestyle=linestyle[j], $
            thick=thick[j]
          DJS_OPLOT, [r[0]-(r[1]-r[0])/2.,r[0]], [histo[0,j],histo[0,j]], $
            color=color[j], linestyle=linestyle[j], thick=thick[j]    
          DJS_OPLOT, [r[0]-(r[1]-r[0])/2.,r[0]-(r[1]-r[0])/2.], [0,histo[0,j]], $
            color=color[j], linestyle=linestyle[j], thick=thick[j]    
          DJS_OPLOT, [r[num-1],r[num-1]+(r[1]-r[0])/2.], $
            [histo[num-1,j],histo[num-1,j]], color=color[j], linestyle=linestyle[j], $
            thick=thick[j]    
          DJS_OPLOT, [r[num-1]+(r[1]-r[0])/2.,r[num-1]+(r[1]-r[0])/2.], $
            [0,histo[num-1,j]], color=color[j], linestyle=linestyle[j], thick=thick[j]    
          
       endfor

; re-draw plot box
       
       if not KEYWORD_SET(oplot) then begin
          PLOT, [0, 0], /noerase, /nodata, ytitle=ytitle, _EXTRA=extra, $
            xrange=[rr[0]-(rr[1]-rr[0]), rr[num-1]+rr[1]-rr[0]], xstyle=1, $
            xtickinterval=xtickinterval, xtickname=xtickname, xminor=xminor, $
            yrange=[0,ymax], $
            ystyle=1, yticks=yticks, yminor=yminor, xticks=xticks
       endif

    endif 
    
    return
    
 end

















