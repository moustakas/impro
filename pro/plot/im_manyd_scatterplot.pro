;+
; NAME:
;   im_manyd_scatterplot
; PURPOSE:
;   plot N-dimensional data sets
; INPUTS:
;   weight       [N] array of data-point weights
;   point        [d,N] array of data points - N vectors of dimension d
;   psfilename   name for PostScript file; if no filename is given, then the
;                  plots will simply be sent to the currently active device
; OPTIONAL INPUTS:
;   nsig         number of sigma for half-width of each plot; default 5
;   label        [d] array of axis labels; default 'x_i'
;   range        [2,d] array of plotting ranges
;   xdims,ydims  indices of data dimensions to use on each x and y axis
;   axis_char_scale size of characters on labels
;   default_font font command to send to set font for plotting
;   extrafun     name of procedure to call after each plot for overplotting
;                (procedure takes two inputs: d1 and d2 indicating
;                x and y dimensions of the current plot)
;   title        puts string title on top of page
; KEYWORDS:
;   nodata       don't plot anything at all, just axes!
;   meanweight      - plot the mean of the weight values that land in
;                     that pixel, rather than the total; don't use
;                     with /conditional!
;   [etc]        [options for hogg_scatterplot, see documentation]
; OUTPUTS:
; OPTIONAL OUTPUTS:
;   manyd        [nx,ny,NX,NY] data
; BUGS:
;   Can get infinite plot ranges.
; COMMENTS:
;   XDIMS and YDIMS can either be:
;     (a) [N] arrays, yielding an NxN set of plots of xdims[0] vs
;         ydims[0], xdims[0] vs ydims[1], etc up to xdims[n-1] vs
;         ydims[n-1]; OR
;     (b) [N,M] arrays, yielding an NxM set of plots of xdims[0,0] vs
;         ydims[0,0], xdims[0,1] vs ydims[0,1] etc.
; DEPENDENCIES:
;   hogg_plot_defaults
;   hogg_scatterplot
;   plus much, much more
; REVISION HISTORY:
;   2002-12-14  re-constructed from ex_max_plot -- Hogg
;   Moustakas - hacked version of hogg_manyd_scatterplot
;-

pro im_manyd_scatterplot, weight,point,psfilename,nsig=nsig, $
  label=label,range=range, xdims=xdims,ydims=ydims, axis_char_scale=axis_char_scale, $
  default_font=default_font, xnpix=xnpix,ynpix=ynpix, nodata=nodata,manyd=manyd,title=title, $
  extrafun=extrafun,meanweight=meanweight, _EXTRA=KeywordsForHoggScatterplot, $
  ccolor=ccolor, outcolor=outcolor, keynote=keynote, upper_triangle=upper_triangle, $
  in_nticks=in_nticks

    if(n_params() lt 2) then begin
       print,'Syntax - hogg_manyd_scatterplot, weight, point [, psfilename, '
       print,'           nsig=, label=, range=, xdims=, ydims=, axis_char_scale=,'
       print,'           default_font=, xnpix=, ynpix=, /nodata, manyd=, '
       print,'           title= ]'
       print,'(also takes keywords associated with hogg_scatterplot)'
       return
    endif

    if keyword_set(keynote) then keycolor = 'white' else keycolor = 'black'
    if n_elements(ccolor) eq 0 then ccolor = keycolor
    if n_elements(outcolor) eq 0 then outcolor = keycolor

; check dimensions
    ndata = n_elements(weight)                   ; N
    dimen = n_elements(point)/n_elements(weight) ; d
    splog, ndata,' data points,',dimen,' dimensions'

; set defaults
    if NOT keyword_set(label) then $
      label= 'x!d'+strcompress(string(lindgen(dimen)),/remove_all)
    if NOT keyword_set(nsig) then nsig = 5d
    if NOT keyword_set(axis_char_scale) then axis_char_scale= 1.75

; only plot the upper triangle
    if keyword_set(upper_triangle) then begin
       xdims = lonarr(dimen,dimen)-1
       ydims = lonarr(dimen,dimen)-1
       arr = lindgen(dimen)
       for iy = 0, dimen-1 do xdims[*,iy] = arr*(arr lt (iy+1))+(arr ge (iy+1))*(-1)
       for ix = 0, dimen-1 do ydims[ix,*] = arr*(arr ge ix)+(arr lt ix)*(-1)
    endif

;print, [[0,-1,-1],[0,1,-1],[0,1,2]]
;print, [[0,-1,-1],[1,1,-1],[2,2,2]]

; which dimensions should we look at?
    if(n_elements(xdims) eq 0 ) then xdims=lindgen(dimen)
    if(n_elements(ydims) eq 0 ) then ydims=lindgen(dimen)
    if(NOT keyword_set(default_font)) then default_font='!3'
    xdimen=n_elements(xdims) 
    ydimen=n_elements(ydims)

;; deal with explicit dimension case
    explicit=0L
    if((size(xdims))[0] eq 2) then begin
       explicit=1L
       if(((size(ydims))[0] ne 2) OR $
         (n_elements(xdims) ne n_elements(ydims))) then $
           message, 'Explicit 2d XDIMS and YDIMS input must be parallel!'
       xdimen= (size(xdims))[1]
       ydimen= (size(xdims))[2]
    endif

; cram inputs into correct format
    point= reform(double(point),dimen,ndata)
    weight= reform(double(weight),ndata)

; compute mean and variance of the whole sample for plot ranges
    if NOT keyword_set(range) then begin
       amp1= total(weight)
       mean1= total(weight##(dblarr(dimen)+1D)*point,2)/amp1
       var1= 0d
       for i=0L,ndata-1 do begin
          delta= point[*,i]-mean1
          var1= var1+weight[i]*product(finite(delta))*delta#delta
       endfor
       var1= var1/amp1
       range= dblarr(2,dimen)
       for d1= 0,dimen-1 do range[*,d1]= mean1[d1]+[-nsig,nsig]*sqrt(var1[d1,d1])
    endif

; save system plotting parameters for later restoration
    bangP= !P
    bangX= !X
    bangY= !Y

; setup postscript file
    xsize = 7.5 & ysize= 7.0
    if keyword_set(psfilename) then begin
       im_plotconfig, 0, psfile=psfilename, keynote=keynote
;      set_plot, 'PS'
;      device, file=psfilename, /portrait, xsize=xsize, ysize=ysize, $
;        xoffset=0.5, yoff=0.5+(10-ysize), /inch, encap=encap, /color, $
;        bits=8

;  dfpsplot, psfilename, /color, bits=8, xsize=xsize, ysize=ysize
;   set_plot, "PS"
;   device, file=psfilename,/inches,xsize=xsize,ysize=ysize, $
;     xoffset=0,yoffset=0.5+(10-ysize),/color, bits=8
;;    xoffset=(8.5-xsize)/2.0,yoffset=(11.0-ysize)/2.0,/color, bits=8
    endif
;   hogg_plot_defaults, axis_char_scale=axis_char_scale
    tiny = 1.d-4
    !X.CHARSIZE= tiny
    !Y.CHARSIZE= !X.CHARSIZE
    !P.MULTI= [xdimen*ydimen,xdimen,ydimen]
    if keyword_set(keynote) then !p.charthick=5.0 else !p.charthick=3.0
    xyouts, 0,0, default_font

    if n_elements(in_nticks) eq 0 then in_nticks = 4

; loop over all pairs of dimensions
    for id2=ydimen-1L,0L,-1 do begin
;  if keyword_set(upper_triangle) then xdimen_end = id2+1 else xdimen_end = xdimen
;  for id1=0L,xdimen_end-1 do begin
       for id1=0L,xdimen-1 do begin
          if(NOT keyword_set(explicit)) then begin
             d1=xdims[id1]
             d2=ydims[id2]
          endif else begin
             d1=xdims[id1,id2]
             d2=ydims[id1,id2]
          endelse
          if d1 lt 0 or d2 lt 0 then begin
             !P.MULTI[0]=!P.MULTI[0]-1L
          endif else begin 

; set axis label properties
             !X.CHARSIZE= tiny
             !Y.CHARSIZE= tiny
             !X.TITLE= ''
             !Y.TITLE= ''

; set plot range
             nticks=in_nticks/axis_char_scale
;           nticks=6/axis_char_scale
             !X.TICKINTERVAL= hogg_interval(range[*,d1],nticks=nticks)
             !Y.TICKINTERVAL= hogg_interval(range[*,d2],nticks=nticks)

; are we on one of the plot edges?
; NB: must run this check before plotting!
             xprevblank=0
             yprevblank=0
             xnextblank=0
             ynextblank=0
             if(keyword_set(explicit)) then begin
                if(id1 gt 0) then $
                  if(xdims[id1-1,id2] eq -1) then $
                    xprevblank=1
                if(id1 lt xdimen-1L) then $
                  if(xdims[id1+1,id2] eq -1) then $
                    xnextblank=1
                if keyword_set(upper_triangle) then begin
                   if(id2 lt ydimen-1L) then if(xdims[id1,id2+1] eq -1) then yprevblank=1
                   if(id2 gt 0) then if(xdims[id1,id2-1] eq -1) then ynextblank=1
                endif else begin
                   if(id2 gt 0) then if(xdims[id1,id2-1] eq -1) then yprevblank=1
                   if(id2 lt ydimen-1L) then if(xdims[id1,id2+1] eq -1) then ynextblank=1
                endelse
             endif else begin
                if (id1 gt 0) then $
                  if (xdims[id1-1] eq -1) then xprevblank=1
                if (id1 lt xdimen-1) then $
                  if (xdims[id1+1] eq -1) then xnextblank=1
                if (id2 gt 0) then $
                  if (ydims[id2-1] eq -1) then ynextblank=1
                if (id2 lt xdimen-1) then $
                  if (ydims[id2+1] eq -1) then yprevblank=1
             endelse
             leftside= 0B
             if (!P.MULTI[0] EQ 0) OR $
               (((!P.MULTI[0]-1) MOD xdimen) EQ (xdimen-1) OR $
               xprevblank eq 1) then leftside= 1B
             topside= 0B
             if (!P.MULTI[0] EQ 0) OR $
               (floor(float(!P.MULTI[0]-1)/xdimen) EQ (ydimen-1) OR $
               yprevblank eq 1) then topside= 1B
             rightside= 0B
             if (((!P.MULTI[0]-1) MOD xdimen) EQ 0 OR $
               xnextblank eq 1) then rightside= 1B
             bottomside= 0B
             if (floor(float(!P.MULTI[0]-1)/xdimen) EQ 0 OR $
               ynextblank eq 1) then bottomside= 1B
;if id1 eq 3 and id2 eq 4 then stop
; plot
             if d1 NE d2 then begin
                if keyword_set(nodata) then begin
                   plot, [0],[0], $
                     xrange=range[*,d1],yrange=range[*,d2], $
                     /nodata, color=im_color(keycolor)
                endif else begin
                   hogg_scatterplot, point[d1,*],point[d2,*], $
                     weight=weight, $
                     xrange=range[*,d1],yrange=range[*,d2], $
                     xnpix=xnpix,ynpix=ynpix, grid=grid, $
                     internal_weight=internal_weight,meanweight=meanweight, $
                     _EXTRA=KeywordsForHoggScatterplot, $
                     color=im_color(keycolor), ccolor=im_color(ccolor), $
                     outcolor=im_color(ccolor), cthick=!p.charthick
                   if(n_elements(manyd) ne xnpix*ynpix*xdimen*ydimen) then $
                     manyd=dblarr(xnpix,ynpix,xdimen,ydimen)
                   manyd[*,*,id1,id2]=grid
                endelse
             endif else begin
                if keyword_set(nodata) then begin
                   plot, [0],[0], $
                     xrange=range[*,d1],yrange=[0,1], $
                     /nodata,yticklen=1d-10
                endif else begin
                   im_plothist, point[d1,*],weight=weight, $
                     xrange=range[*,d1],npix=xnpix,yticklen=1d-10, $
                     meanweight=meanweight, thick=!p.charthick, $
                     _extra=KeywordsForHoggScatterplot
;                  hogg_plothist, point[d1,*],weight=weight, $
;                    xrange=range[*,d1],npix=xnpix,yticklen=1d-10, $
;                    meanweight=meanweight, thick=!p.charthick, $
;                    _extra=extra
                endelse
             endelse
;           print, id1, id2, bottomside, topside, leftside, rightside

; make axis labels afterwards
             if bottomside then begin
                axis,!X.CRANGE[0],!Y.CRANGE[0],xaxis=0, $
                  xtitle=label[d1],xcharsize=axis_char_scale*1.1, color=im_color(keycolor)
             endif
             if topside then begin
                axis,!X.CRANGE[0],!Y.CRANGE[1],xaxis=1, $
                  xtitle=label[d1],xcharsize=axis_char_scale*1.1, color=im_color(keycolor)
             endif
             if leftside AND (d1 NE d2) then begin
                axis,!X.CRANGE[0],!Y.CRANGE[0],yaxis=0, $
                  ytitle=label[d2],ycharsize=axis_char_scale*1.1, color=im_color(keycolor)
             endif
             if rightside AND (d1 NE d2) then begin
                axis,!X.CRANGE[1],!Y.CRANGE[0],yaxis=1, $
                  ytitle=label[d2],ycharsize=axis_char_scale*1.1, color=im_color(keycolor)
             endif

             if(keyword_set(extrafun)) then begin
                call_procedure, extrafun, d1, d2, point=point
             endif

; end loops and close file
          endelse 
       endfor 
    endfor

    if(keyword_set(title)) then begin
       !P.MULTI=[0,1,1]
       !X.RANGE=[0,1]
       !Y.RANGE=[0,1]
       !X.CRANGE=[0,1]
       !Y.CRANGE=[0,1]
       !X.S=[0,1]
       !Y.S=[0,1]
       xyouts,0.5,1.03,title,align=0.5
    endif

;if keyword_set(psfilename) then dfpsclose
;if keyword_set(psfilename) then device, /close
    if keyword_set(psfilename) then im_plotconfig, /psclose, keynote=keynote, psfile=psfilename

; restore system plotting parameters
    !P= bangP
    !X= bangX
    !Y= bangY
;set_plot,'x'

return
end
