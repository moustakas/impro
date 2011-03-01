;+
; NAME:
; plot_earth
;
;
; PURPOSE:
; overplot the earth onto a plot measured in Re
;
;
; INPUTS:
; none
;
;
; OPTIONAL INPUTS:
; none
;
;
; KEYWORD PARAMETERS:
; FANCY - plot the earth as a globe
; BLANK - for non fancy output blank out the region, ingnored when
;         fancy is specified
; TIME - the time of the fancy earth to plot, defults to zero
;        This can be specifed in fractional hours (e.g. 12.345) or as
;        an EUVtime (e.g. 20012001234   yyyydddhhmmss)
;
;
; OUTPUTS:
; the earth onto the current plot
;
;
; OPTIONAL OUTPUTS:
; none
;
;
; COMMON BLOCKS:
; none
;
;
; SIDE EFFECTS:
; - puts the earth on the plot
; - messed up the colortable, need to reload ct after a run of this
;
; RESTRICTIONS:
; - Only works well for iso plots measured in Re from the earth.
;   Could be modified to work in meters or the like at a later date
; - Only works for directly over the North Pole as that I what I am
;   doing
; now
;
;
; EXAMPLE:
; loadct, 13
; contour, dist(10), findgen(10)-5, findgen(10)-5, /fill, nlevels=30, /iso
; plot_earth, /fancy, time=12.345
;
; MODIFICATION HISTORY:
;
;       Sun Jan 28 15:41:46 2007, Brian Larsen
;       <larsen@ssel.montana.edu>
;
;               written and tested with help from comp.lang.idl-pvwave
;             http://groups.google.com/group/comp.lang.idl-pvwave/
;             browse_thread/thread/4ba6a41393fc8f6f/?hl=en#
;
;-

pro plot_earth, BLANK=blank, FANCY=fancy, TIME=time

IF NOT KEYWORD_SET(fancy) THEN BEGIN
   IF KEYWORD_SET(blank)  THEN BEGIN
       circ_r = fltarr(200)
       circ_r[*] = 1
       circ = findgen(200)*2*!pi/200
       x_=circ_r * sin(circ)
       y_=circ_r * cos(circ)
       polyfill, x_, y_, color=0
   ENDIF
   rad_ = 2.*!pi*findgen(100)/100.
   earth_ = fltarr(100)
   earth_[*] = 1
   oplot,  /polar, earth_[*], 2.*rad_[*]
   plots, [0,1],[0,0], linestyle=2
ENDIF ELSE BEGIN
   ;; for time we are expecting either an euv time or factional hours
   ;; - an euvtime will be a string
   IF N_ELEMENTS(time) EQ 0 THEN rot = 0 ELSE BEGIN
       IF size(time, /type) EQ 7 THEN BEGIN
           rot = euv_date2arr( time)
           rot = (rot[2]+rot[3]/60.+rot[4]/3600.)/24.*360 - 270
       ENDIF ELSE rot = time/24. * 360 - 270
   ENDELSE

   ;; day side
   pos1 = convert_coord([0,1], [-1,1], /data, /to_normal)
   ;; night side
   pos2 = convert_coord([-1,0], [-1,1], /data, /to_normal)

   ;; day side
   loadct, 12
   MAP_SET,90,0,rot,/ORTHOGRAPHIC,/ISOTROPIC, /CONTINENTS,/HORIZON, $
     E_continents={FILL:1, color:23}, $
     position=[pos1[0,0], pos1[1,0], pos1[0,1], pos1[1,1]], $
     /noerase, /noborder, $
     E_HORIZON={FILL:1, COLOR:100}, $
     limit=[0,rot,90,rot+180]

   ;; night side
   loadct, 0
   MAP_SET,90,0,rot,/ORTHOGRAPHIC,/ISOTROPIC, /CONTINENTS,/HORIZON, $
     E_continents={FILL:1, color:119}, $
     position=[pos2[0,0], pos2[1,0], pos2[0,1], pos2[1,1]], $
     /noerase, /noborder, $
     E_HORIZON={FILL:1, COLOR:33}, $
     limit=[0,rot+180,90,rot]

ENDELSE

END
