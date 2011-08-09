;+
; NAME:
;	IM_WINDOW
;
; PURPOSE:
;	Open up an IDL window appropriate to the user's screen size.  
;
; INPUTS:
;	windx  - the index number of the window to open (default is to
;	         open the next available window)
;
; OPTIONAL INPUTS:
;	xratio - the ratio of the x-size of the window relative to the
;		 x-size of the monitor, which should be less than 1;
;		 default to 0.5
;	yratio - as xratio above but for the y-size (default is xratio)
;	xpos   - x-position of the window in device coordinates
;	         (default is the upper-right hand corner)
;	ypos   - y-position of the window in device coordinates
;
; KEYWORDS USED:
;       square - force YSIZE=XSIZE
;
; PROCEDURES USED:
;	WINDOW, GET_SCREEN_SIZE()
;
; MODIFICATION HISTORY:
;	John Moustakas, 2000 September 15, UofA
;       jm05jan21uofa - added SQUARE keyword; removed TITLE optional
;                       input 
;
; Copyright (C) 2000, 2005, John Moustakas
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

pro im_window, windx, xratio=xratio, yratio=yratio, xpos=xpos, $
  ypos=ypos, square=square

    if (n_elements(xratio) eq 0L) then xratio = 0.5
    if (n_elements(yratio) eq 0L) then yratio = xratio

    wsz = get_screen_size()
    xsize = wsz[0]*xratio
    ysize = wsz[1]*yratio

    if keyword_set(square) then ysize = xsize
    
    if (n_elements(windx) eq 0L) then begin
       windx = 0L
       free = 1L 
    endif

;   if (n_elements(xpos) eq 0L) then xpos = wsz[0]-xsize
;   if (n_elements(ypos) eq 0L) then ypos = wsz[1]-ysize

    window, windx, xsize=xsize, ysize=ysize, xpos=xpos, ypos=ypos, free=free

return
end
