function get_ellipse_indices, x, y, major=major, minor=minor, $
  angle=angle, xcenter=xcenter, ycenter=ycenter, debug=debug
; jm10jul09ucsd - see SINGS_LOG12OH_HIIREGIONS for how to call this
; routine; in particular, the angle is defined as in TVELLIPSE;
; portions of the algorithm taken from DIST_ELLIPSE

    if (n_elements(xcenter) eq 0) then xcenter = 0.0
    if (n_elements(ycenter) eq 0) then ycenter = 0.0
    
; rotate pixels to match ellipse orientation
    ang = (90-angle)*!dtor
    rot_ccw = [ [cos(ang),sin(ang)], [-sin(ang),cos(ang)] ]
    xyprime = transpose([[(x-xcenter)],[(y-ycenter)]]) ## rot_ccw

    indx = where(sqrt((xyprime[0,*]*major/minor)^2 + xyprime[1,*]^2) lt major,nindx)

    if keyword_set(debug) then begin
       djs_plot, xyprime[0,*], xyprime[1,*], ps=6
       if (nindx ne 0) then djs_oplot, xyprime[0,indx], $
         xyprime[1,indx], ps=6, color='orange'
    endif

return, indx
end
