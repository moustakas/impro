;+
; NAME: ARM_OPLOTLINE
;       
; CATEGORY: plotting
;
; PURPOSE: overplots 1 or more lines on an existing plot
;
; CALLING SEQUENCE: ARM_OPLOTLINE, value, [linestyle, color, thick,
;                                          /horizontal, /identity]
;
; INPUTS: value - array of 1 or more line values (default=1.0)
;       
; OPTIONAL INPUTS: 
;    linestyle - line style number {0,1,2,3...} (PLOT keyword)
;    color     - string name of desired color (DJS_PLOT keyword)
;    thick     - line thickness (PLOT keyword)
;    xthick    - x-axis thickness of existing plot
;    ythick    - y-axis thickness of existing plot
;
; KEYWORDS: 
;    horizontal - plots a horizontal line rather than a vertical one
;    identity   - plots y = x
;
; EXAMPLE: 
;    IDL> plot, [0,0], /nodata, xr=[0,10], xst=1, yr=[-5,5], yst=1
;    IDL> ARM_OPLOTLINE, [2,4,7.5], color=['red','blue','red'], $
;         thick=[2,2,1]
;
; PROCEDURES USED: DJS_OPLOT
;
; COMMENTS: accepts all PLOT commands/keywords
;
; MODIFICATION HISTORY:
;    written by A.R.Marble, Steward Obs., 2002
;    defaults modified to fix bug, ARM, October 20, 2003
;
; TODO : don't plot over axis (done, but problem with log plots)
;        put lines behind plot?
;-

pro ARM_OPLOTLINE, value, linestyle=linestyle, color=color, $
       thick=thick, xthick=xthick, ythick=ythick, $
       horizontal=HORIZONTAL, identity=IDENTITY

    if N_ELEMENTS(xthick) ne 0 then xthick = round(xthick) else xthick = !x.thick
    if N_ELEMENTS(ythick) ne 0 then ythick = round(ythick) else ythick = !y.thick

    if N_ELEMENTS(value) eq 0 or KEYWORD_SET(identity) then value = 1.0
    n = N_ELEMENTS(value)       ; number of lines to plot

; default line styles

    nstyle = N_ELEMENTS(linestyle)
    if nstyle ne n then begin
       if nstyle eq 0 then linestyle = REPLICATE(0, n) $
       else if nstyle eq 1 then linestyle = REPLICATE(linestyle, n) $
       else if nstyle gt n then linestyle = linestyle[0:n-1] $
       else linestyle = [linestyle, REPLICATE(0, n-nstyle)] 
    endif
    
; default colors

    ncolor = N_ELEMENTS(color)
    if ncolor ne n then begin
       if ncolor eq 0 then color = REPLICATE('', n) $
       else if ncolor eq 1 then color = REPLICATE(color, n) $
       else if ncolor gt n then color = color[0:n-1] $
       else color = [color, REPLICATE('', n-ncolor)]
    endif

; default thicknesses

    nthick = N_ELEMENTS(thick)
    if nthick ne n then begin
       if nthick eq 0 then thick = REPLICATE(!p.thick, n) $
       else if nthick eq 1 then thick = REPLICATE(thick, n) $
       else if nthick gt n then thick = thick[0:n-1] $
       else thick = [thick, REPLICATE(!p.thick, n-nthick)] 
    endif

; get pixel size in data coordinates
    
    devicecoords = CONVERT_COORD([[0.,0.],[1.,1.]], /device, /to_data)
    xpxlsize = devicecoords[0,1] - devicecoords[0,0]
    ypxlsize = devicecoords[1,1] - devicecoords[1,0]

; get plot window range in data coordinates

    normalcoords = CONVERT_COORD([[0.,0.],[1.,1.]], /normal, /to_data)

    xl = normalcoords[0,0] + ROUND(xthick / 2.) * xpxlsize
    xr = normalcoords[0,1] - ROUND(xthick / 2.) * xpxlsize
    yb = normalcoords[1,0] + ROUND(ythick / 2.) * ypxlsize
    yt = normalcoords[1,1] - ROUND(ythick / 2.) * ypxlsize

;;    if KEYWORD_SET(xlog) then !x.crange = 10^!x.crange
;;    if KEYWORD_SET(ylog) then !y.crange = 10^!y.crange
;;    xl = !x.crange[0] + ROUND(xthick / 2.) * xpxlsize
;;    xr = !x.crange[1] - ROUND(xthick / 2.) * xpxlsize
;;    yb = !y.crange[0] + ROUND(ythick / 2.) * ypxlsize
;;    yt = !y.crange[1] - ROUND(ythick / 2.) * ypxlsize

    if KEYWORD_SET(identity) then begin

; output warning if keywords conflict

       if KEYWORD_SET(horizontal) then MESSAGE, /continue, $
         'IDENTITY keyword overrules HORIZONTAL keyword '

; y = x (identity line)

       y1 = MIN([xl,yb]) > yb
       y2 = MAX([xr,yt]) < yt

       x1 = (y1 / value) > xl
       x2 = (y2 / value) < xr

    endif else if KEYWORD_SET(horizontal) then begin

; y = value (horizontal line)

       x1 = REPLICATE(xl, n)
       x2 = REPLICATE(xr, n)
       
       y1 = (value > yb) < yt
       y2 = y1

    endif else begin

; x = value (vertical line)

       x1 = (value > xl) < xr
       x2 = x1

       y1 = REPLICATE(yb, n)
       y2 = REPLICATE(yt, n)
       
    endelse
    
; overplot desired lines

;; old = tvrd()

    for i=0,n-1 do DJS_OPLOT, [x1[i], x2[i]], [y1[i], y2[i]], $
      linestyle=linestyle[i], color=color[i], thick=thick[i]

;; new = tvrd()
;; copy=new
;; wh = where(new ne old and old ne 0, count)
;; if count gt 0 then new[wh] = old[wh]
;; 
;; stop

    return


 end
