;+
; NAME:
;       IM_OFFSET_AND_ROTATE()
;
; PURPOSE:
;       Apply an arbitrary offset and rotation to a position.  
;
; CALLING SEQUENCE:
;       xynew = im_offset_and_rotate(xy,pa,xoffset=,yoffset=)
;
; INPUTS:
;       xy - [XPOSITION,YPOSITION] to rotate and offset; can also be a
;            [2,*] element array
;       pa - rotation angle measured counter-clockwise from the y axis 
;            (identical to astronomical position angles) [degrees] 
;
; OPTIONAL INPUTS:
;       xoffset - shift each X position by XOFFSET (scalar, default 0)
;                 (perpendicular to the slit)
;       yoffset - shift each Y position by YOFFSET (scalar, default 0)
;                 (along the slit)
;
; OUTPUTS:
;       xynew   - offset and rotated XY positions
;
; COMMENTS & EXAMPLE:
;       The input XY is a position in the rotated (primed) coordinate
;       axis PA relative to an un-rotated (cardinal) coordinate axis.
;       For example, if your position is [-1,0] in a coordinate axis
;       that is rotated by PA=45 degrees, then your position in the
;       un-rotated frame is  
;
;       IDL> print, im_offset_and_rotate([-1.0,0],45)
;            -0.707107    -0.707107
;
;       Also see IM_OPLOT_BOX for another application of this
;       routine.  
;
; MODIFICATION HISTORY:
;       J. Moustakas, 2003 May 22, U of A
;       jm04mar16uofa - bug fix when using offsets
;-

function im_offset_and_rotate, xy, pa, xoffset=xoffset, yoffset=yoffset

    xysize = size(xy,/dimension)
    if (xysize[0] ne 2L) then begin
       print, 'XY must be at least a two-dimensional array.'
       return, -1L
    endif

    if n_elements(pa) eq 0L then begin
       print, 'No position angle PA provided . . . setting PA = 0.'
       pa = 0.0
    endif
    
    if n_elements(xoffset) eq 0L then xoffset = 0.0
    if n_elements(yoffset) eq 0L then yoffset = 0.0
    
    parad = pa*!dtor     ; [radians]

; counter-clockwise rotation matrix measured positive from the y axis     

    xyoffset = xy
    xyoffset[0,*] = xy[0,*] + xoffset
    xyoffset[1,*] = xy[1,*] + yoffset
    
    rot_ccw = [ [cos(parad),sin(parad)], [-sin(parad),cos(parad)] ]
;   xynew = float(xy) ## rot_ccw + ([xoffset,yoffset] # (fltarr(4)+1))
    xynew = float(xyoffset) ## rot_ccw
    
return, xynew
end    
