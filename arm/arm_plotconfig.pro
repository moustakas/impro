;+
; NAME: ARM_PLOTCONFIG
;       
; CATEGORY: plotting
;
; PURPOSE: determine coordinates for desired plot arrangement
;
; CALLING SEQUENCE: ARM_PLOTCONFIG, [psfile, xpage, ypage, xmargin,
;                   ymargin, nx, ny, width, height, xspace, yspace, 
;                   coords=coords, /landscape]
;
; OPTIONAL INPUTS:
;    psfile  - name of postscript file to be opened
;    xpage   - dimension of plot area from left to right (inches,
;              default=8.5)
;    ypage   - dimension of plot area from top to bottom (inches,
;              default=11.0)
;    xmargin - left/right margin, array of 1 or 2 values (inches,
;              default=[1.25,0.75])
;    ymargin - top/bottom margin, array of 1 or 2 values (inches,
;              default=[0.75,0.75])
;    nx      - number of plots from left to right (default=1)
;    ny      - number of plots from top to bottom (default=2)                  
;    width   - width of plots, array of 1 or nx values (inches)
;    height  - height of plots, array of 1 or ny values (inches)
;    xspace  - horizontal spacing between plots, array of 1 or nx-1
;              values, (inches)
;    yspace  - vertical spacing between plots, array of 1 or ny-1
;              values (inches)
;
; KEYWORDS: 
;   landscape - rotates orientation of plots 90 degrees if set
;   show      - display schematic of plot arrangement
;   writeover - write over old postscript file without prompting user
;   bw        - black & white instead of color
;
; OPTIONAL OUTPUTS: 
;    coords - a N by (nx x ny) array of plot coordinates (normalized),
;             the N elements are [x0,y0,x2,y2] where the plot corners
;             are ordered 0-3 counter-clockwise from the lower-left
; EXAMPLE:
;    IDL> ARM_PLOTCONFIG, nx=3, ny=5, xsp=0.1, ysp=0.5, coords=crds
;    IDL> window, xsize=400, ysize=500
;    IDL> for i=0,14 do PLOT, [0,0], /nodata, position=crds[*,i]
;
; PROCEDURES USED:
;
; COMMENTS: specify as few or as many inputs as desired depending
;           on particular needs
;
; MODIFICATION HISTORY:
;    written by A.R.Marble, Steward Obs., January 2002
;    overhauled, January 2003, A.R.Marble
;    postscript device features added and name changed from PAGEMAKER, 
;       September 2003, A.R.Marble
;    graphical display of plot configuration added, ARM, October 31, 2003
;    XOFFSET/YOFFSET bug in LANDSCAPE mode fixed; ARM, 2005 Aug 02
;    BW keyword added; ARM, 2005 Oct 21
;    YOFFSET bug in LANDSCAPE/ENCAPSULATED mode fixed; ARM, 2005 Oct 21
;-

pro ARM_PLOTCONFIG, psfile=psfile, xpage=xpage, ypage=ypage, $
       xmargin=xmargin, ymargin=ymargin, nx=nx, ny=ny, width=width, $
       height=height, xspace=xspace, yspace=yspace, coords=coords, $
       landscape=LANDSCAPE, show=SHOW, writeover=WRITEOVER, bw=BW, _extra=extra

; specify default values

    dflt_xpage =  8.5           ; xpage default (inches)
    dflt_ypage = 11.0           ; ypage default (inches)
    dflt_xmargin = 1.0          ; xmargin default (inches)
    dflt_ymargin = 1.0          ; ymargin default (inches)
    dflt_nx      = 1            ; nx default
    dflt_ny      = 1            ; ny default
    dflt_xspace  = 0.75         ; xspace default (inches)
    dflt_yspace  = 0.75         ; yspace default (inches)

; error checking and default values

    if not KEYWORD_SET(xpage) then xpage = dflt_xpage $
    else if N_ELEMENTS(xpage) gt 1 then $
      MESSAGE, 'XPAGE must be a single value.'
    if not KEYWORD_SET(ypage) then ypage = dflt_ypage $
    else if N_ELEMENTS(ypage) gt 1 then $
      MESSAGE, 'YPAGE must be a single value.'
    if KEYWORD_SET(landscape) then begin
       xp = ypage & yp = xpage 
    endif else begin
       xp = xpage & yp = ypage
    endelse
       
    if N_ELEMENTS(xmargin) eq 0 then xm = [1,1] * dflt_xmargin $
    else if N_ELEMENTS(xmargin) eq 1 then xm = [1,1] * xmargin $
    else if N_ELEMENTS(xmargin) eq 2 then xm = xmargin $
    else if N_ELEMENTS(xmargin) gt 2 then $
      MESSAGE, 'XMARGIN must have at most 2 elements.'
    if N_ELEMENTS(ymargin) eq 0 then ym = [1,1] * dflt_ymargin $
    else if N_ELEMENTS(ymargin) eq 1 then ym = [1,1] * ymargin $
    else if N_ELEMENTS(ymargin) eq 2 then ym = ymargin $
    else if N_ELEMENTS(ymargin) gt 2 then $
      MESSAGE, 'YMARGIN must have at most 2 elements.'
 
    if not KEYWORD_SET(nx) then nx = dflt_nx $
    else if N_ELEMENTS(nx) gt 1 or nx le 0 then $
      MESSAGE, 'NX must be a single value which is > 0.'
    if not KEYWORD_SET(ny) then ny = dflt_ny $
    else if N_ELEMENTS(ny) gt 1 or ny le 0 then $
      MESSAGE, 'NY must be a single value which is > 0.'

    if nx gt 1 then begin 
       if N_ELEMENTS(xspace) eq 0 then xsp = $
         [0.0, REPLICATE(dflt_xspace, nx-1)] $
       else if N_ELEMENTS(xspace) eq 1 then xsp = $
         [0.0, REPLICATE(xspace, nx-1)] $
       else if N_ELEMENTS(xspace) eq nx-1 then xsp = [0.0, xspace] $
       else MESSAGE, 'XSPACE must have 1 or NX-1 elements.'
    endif else xsp = 0.0
    if ny gt 1 then begin 
       if N_ELEMENTS(yspace) eq 0 then ysp = $
         [0.0, REPLICATE(dflt_yspace, ny-1)] $
       else if N_ELEMENTS(yspace) eq 1 then ysp = $
         [0.0, REPLICATE(yspace, ny-1)] $
       else if N_ELEMENTS(yspace) eq ny-1 then ysp = [0.0, yspace] $
       else MESSAGE, 'YSPACE must have 1 or NY-1 elements.'
    endif else ysp = 0.0

    pmulti = !p.multi

    !p.multi = [nx*ny, nx, ny, 0, 0] ; specify when new pages needed

; calculate widths and heights of plot boxes, if not specified

    if not KEYWORD_SET(width) then w = $
      REPLICATE((xp-TOTAL(xm)-TOTAL(xsp))/FLOAT(nx), nx) $
    else if N_ELEMENTS(width) eq 1 then w = REPLICATE(width, nx) $
    else if N_ELEMENTS(width) ne nx then $
      MESSAGE, 'WIDTH must have 1 or NX elements.' else w = width
    if TOTAL(w)+TOTAL(xsp)+TOTAL(xm) gt xp then $
      MESSAGE, 'WIDTH and XSPACE values exceed usable area.'
    if not KEYWORD_SET(height) then h = $
      REPLICATE((yp-TOTAL(ym)-TOTAL(ysp))/FLOAT(ny), ny) $
    else if N_ELEMENTS(height) eq 1 then h = REPLICATE(height, ny) $
    else if N_ELEMENTS(height) ne ny then $
      MESSAGE, 'HEIGHT must have 1 or NY elements.' else h = height
    if TOTAL(h)+TOTAL(ysp)+TOTAL(ym) gt yp then $
      MESSAGE, 'HEIGHT and YSPACE values exceed usable area.'

; calculate corner coordinates for each plot box 

    coords = FLTARR(4, nx*ny)

    for i=0,ny-1 do begin 

       for j=0,nx-1 do begin

          x0 = xm[0] + TOTAL(w[0:j]) + TOTAL(xsp[0:j]) - w[j] 
          x2 = xm[0] + TOTAL(w[0:j]) + TOTAL(xsp[0:j])
          y0 = yp - ym[0] - TOTAL(h[0:i]) - TOTAL(ysp[0:i])
          y2 = yp - ym[0] - TOTAL(h[0:i]) - TOTAL(ysp[0:i]) + h[i] 
          
          coords[*, j+i*nx] = [x0, y0, x2, y2]

       endfor
       
    endfor   
          
; convert physical coordinates to normal coordinates

    coords[[0,2],*] = coords[[0,2],*] / xp
    coords[[1,3],*] = coords[[1,3],*] / yp

; display schematic of plot arrangement, if desired

    if KEYWORD_SET(show) then begin

;   open window

       wsize = 600. ; largest dimension of window
       if xp gt yp then xsize = wsize else $
         xsize = xp / yp * wsize
       WINDOW, 28, xsize=xsize, ysize=yp/xp*xsize, $
         title='ARM_PLOTCONFIG schematic'

;   plot outline of 8.5x11 page

       PLOT, [0,0], /nodata, position=[0,0,1,1], /normal, $
         xstyle=5, ystyle=5, xrange=[0,1], yrange=[0,1]
       OPLOT, [0,8.5/xp], [1,1]*(1-11.0/yp), line=1
       OPLOT, [1,1]*(8.5/xp), [1-11.0/yp,1], line=1
       OPLOT, [0,8.5/xp], [1,1], line=1
       OPLOT, [0,0], [1-11.0/yp,1], line=1
       XYOUTS, 0.5*8.5/xp, 1-11.0/yp, align=0.5, '8.5"'
       XYOUTS, 8.5/xp,1-0.5*11.0/yp, align=0.5, orientation=90, '11.0"'

;   overplot plot boxes

       for j=0,nx*ny-1 do begin

          if j eq 0 then noerase = 1 else noerase = 0
         
          PLOT, [0,0], /nodata, position=coords[*,j], /normal, $
            xticks=1, xtickname=REPLICATE(' ', 2), xminor=1, $
            yticks=1, ytickname=REPLICATE(' ', 2), yminor=1, noerase=noerase

       endfor

    endif

; open postscript file if desired

    if KEYWORD_SET(psfile) then begin
       
;   check to see if file already exists

       if not KEYWORD_SET(writeover) then begin
          
          check = FILE_SEARCH(psfile)
          if check[0] ne '' then begin
             
             PRINT
             PRINT, 'The file '+psfile+' already exists.'
             continue = ARM_YESORNO('Write over it?')          
             
             if not continue then $
               MESSAGE, 'You must specify an alternate file name.'
             
          endif
          
       endif

       if STRUPCASE(STRMID(psfile, 2, 3, /reverse)) eq 'EPS' $
         then encapsulated = 1 else encapsulated = 0

       SET_PLOT, 'ps'
       if not KEYWORD_SET(landscape) then begin
          xoffset = (8.5 - xp) / 2. > 0
          yoffset = (11.0 - yp) / 2. > 0
       endif else begin
          xoffset = (8.5 - yp) / 2. > 0
          if encapsulated then yoffset = xp $
          else yoffset = 11.0 - (11.0 - xp) / 2. > 0
       endelse
       
       if KEYWORD_SET(landscape) then landscape = 1 else landscape = 0

       if KEYWORD_SET(bw) then color = 0 else color = 1
       
       DEVICE, file=psfile, color=color, /inches, xsize=xp, ysize=yp, $
         xoffset=xoffset, yoffset=yoffset, landscape=landscape, $
         encapsulated=encapsulated, _extra=extra

    endif

    !p.multi = pmulti ; restore
    
    return
    
 end




