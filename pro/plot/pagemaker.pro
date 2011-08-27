;+
;
;  NAME:  
;    PAGEMAKER
;
;  PURPOSE:
;    Returns plotting positions for desired page specifications,
;    pagemaker can be run with no inputs (all defaults) or as much
;    specification can be given as desired
;
;  OPTIONAL INPUT:  
;
;     XPAGE   - horizontal dimension of page in inches (default=8.5)
;     YPAGE   - vertical dimension of page in inches (default=11.0)
;     XMARGIN - left/right margin size in inches, can be an array (default=1.0)
;     YMARGIN - top/bottom margin size in inches, can be an array (default=1.0)
;     NX      - number of horizontal boxes across page (default=1)
;     NY      - number of vertical boxes across page (default=2)
;     WIDTH   - width of plot box(es) in inches, can be an array 
;               for left to right
;     HEIGHT  - height of plot box(es) in inches, can be an array
;               for top to bottom
;     XSPACE  - horizontal spacing between boxes in inches, can be
;               an array for left to right (default=0.1*WIDTH)
;     YSPACE  - vertical spacing between boxes in inches, can be
;               an array for top to bottom (default=0.1*HEIGHT)
;
;  OUTPUT: position - a 4 x (NX*NY) array, position coords for each
;                     plot box - [x0,y0,x1,y1]
;
;  KEYWORDS: 
;
;     landscape - switches from portrait mode to landscape mode
;     normal    - returns normalized coordinates
;
;   EXAMPLE: IDL> pagemaker, nx=3, ny=5, xspace=0., yspace=0., $
;                   position=pos, /normal
;            IDL> window, xsize=500, ysize=647
;            IDL> for i=0,14 do plot,findgen(10),position=pos[*,i],/noerase
;
;  HISTORY: written by Andy Marble (ARM), Steward Observatory,  Jan 2002
;           overhauled by ARM, Jan 2003 - currently still under revision
;-

pro pm_err, message
    print,'ERROR (pagemaker.pro): '+message
    retall
 end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro pm_swap, a, b
    tmp = a
    a   = b
    b   = a
 end
 
pro pagemaker, xpage=xpage, ypage=ypage, xmargin=xmargin, ymargin=ymargin, $
               nx=nx, ny=ny, width=width, height=height, $ 
               xspace=xspace, yspace=yspace, position=position, $
               landscape=LANDSCAPE, normal=NORMAL
;              postscript=postscript
;              ratio=ratio
;              show=SHOW
;              nodim=NODIM
    
;   DEFAULT VALUES AND ERROR CHECKING
    if not keyword_set(xpage) then xpage = 8.5 $
    else if n_elements(xpage) gt 1 then pm_err,'XPAGE must be a single value.'
    if not keyword_set(ypage) then ypage = 11.0 $
    else if n_elements(ypage) gt 1 then pm_err,'YPAGE must be a single value.'
    if keyword_set(landscape) then pm_swap, xpage, ypage

    if n_elements(xmargin) eq 0 then xmargin = [1.0, 1.0] $
    else if n_elements(xmargin) eq 1 then xmargin = [xmargin,xmargin] $
    else if n_elements(xmargin) gt 2 then pm_err,'XMARGIN must have at most 2 elements.'
    if n_elements(ymargin) eq 0 then ymargin = [1.0, 1.0] $
    else if n_elements(ymargin) eq 1 then ymargin = [ymargin,ymargin] $
    else if n_elements(ymargin) gt 2 then pm_err,'YMARGIN must have at most 2 elements.'

    if not keyword_set(nx) then nx = 1 $
    else if n_elements(nx) gt 1 or nx eq 0 then pm_err,'NX must be a single value which is > 0.'
    if not keyword_set(ny) then ny = 2 $
    else if n_elements(ny) gt 1 or ny eq 0 then pm_err,'NY must be a single value which is > 0.'

    if nx gt 1 then begin 
       if n_elements(xspace) eq 0 then $
         xspace = [0.0, replicate(0.1*(xpage-total(xmargin))/float(nx), nx-1)] $
       else if n_elements(xspace) eq 1 then xspace = [0.0, replicate(xspace, nx-1)] $
       else if n_elements(xspace) eq nx-1 then xspace = [0.0, xspace] $
       else pm_err,'XSPACE must have 1 or NX-1 elements.'
    endif else xspace = 0.0
    if ny gt 1 then begin 
       if n_elements(yspace) eq 0 then $
         yspace = [0.0, replicate(0.1*(ypage-total(ymargin))/float(ny), ny-1)] $
       else if n_elements(yspace) eq 1 then yspace = [0.0, replicate(yspace, ny-1)] $
       else if n_elements(yspace) eq ny-1 then yspace = [0.0, yspace] $
       else pm_err,'YSPACE must have 1 or NY-1 elements.'
    endif else yspace = 0.0

    if not keyword_set(width) then $
      width = replicate((xpage-total(xmargin)-total(xspace))/float(nx), nx) $
    else if n_elements(width) eq 1 then width = replicate(width,nx) $
    else if n_elements(width) ne nx then pm_err,'WIDTH must have 1 or NX elements.'
    if total(width)+total(xspace)+total(xmargin) gt xpage then $
      pm_err,'WIDTH and XSPACE values exceed usable area.'
    if not keyword_set(height) then $
      height = replicate((ypage-total(ymargin)-total(yspace))/float(ny), ny) $
    else if n_elements(height) eq 1 then height = replicate(height,nx) $
    else if n_elements(height) ne ny then pm_err,'HEIGHT must have 1 or NY elements.'
;   if total(height)+total(yspace)+total(ymargin) gt ypage then $
;     pm_err,'HEIGHT (or RATIO) and YSPACE values exceed usable area.'
    
;   CALCULATE COORDINATES OF LOWER-LEFT AND UPPER-RIGHT CORNERS FOR EACH PLOT BOX
    position = fltarr(4,nx*ny)
    for i=0,ny-1 do for j=0,nx-1 do position[*,j+i*nx] = $
      [xmargin[0]+total(width[0:j])-width[j]+total(xspace[0:j]), $
       ypage-ymargin[0]-total(height[0:i])-total(yspace[0:i]), $
       xmargin[0]+total(width[0:j])+total(xspace[0:j]), $
       ypage-ymargin[0]-total(height[0:i])+height[i]-total(yspace[0:i])]
          
;   IF DESIRED, NORMALIZE COORDINATES TO DIMENSIONS OF PAGE
    if keyword_set(normal) then begin
       position[0,*] = position[0,*] / xpage
       position[1,*] = position[1,*] / ypage
       position[2,*] = position[2,*] / xpage
       position[3,*] = position[3,*] / ypage
    endif

    return
    
 end




