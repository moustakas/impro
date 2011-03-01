;+
; NAME:
;       CIndex
;
; PURPOSE:
;       This is a program for viewing the current colors in the
;       colortable with their index numbers overlayed on each color.
;       On 24-bit systems you must click the cursor in the graphics window
;       to see the colors in the current color table.
;
; AUTHOR:
;       FANNING SOFTWARE CONSULTING
;       David Fanning, Ph.D.
;       1645 Sheely Drive
;       Fort Collins, CO 80526 USA
;       Phone: 970-221-0438
;       E-mail: davidf@dfanning.com
;       Coyote's Guide to IDL Programming: http://www.dfanning.com
;
; CATEGORY: Graphics
;
; CALLING SEQUENCE:  CIndex
;
; INPUTS:   None.
;
; Optional Inputs:   None
;
; OUTPUTS:  None
;
; OPTIONAL OUTPUTS:  None
;
; KEYWORD Parameters:   None
;
; COMMON BLOCKS:  None
;
; SIDE EFFECTS:   None
;
; RESTRICTIONS:   Reqires XCOLORS and TVIMAGE from the Coyote Library:
;
;                     http://www.dfanning.com/programs/xcolors.pro
;                     http://www.dfanning.com/programs/xtvimage.pro
;
; PROCEDURE:
;
;  Draws a 31x25 set of small rectangles in 256 different colors.
;  Writes the color index number on top of each rectangle.
;
; MODIFICATION HISTORY:  Written by David Fanning, May 1995
;
;  Widgetized and made it work in 24-bit color. Colors are
;     updated by clicking in window. 22 Oct 98. DWF
;  Replace POLYFILL with TV command to avoid underflow error in
;     Z-buffer. 8 March 99. DWF
;  Fixed a problem with 24-bit devices with color decomposition ON. 15 Feb 2000. DWF.
;-
;
;###########################################################################
;
; LICENSE
;
; This software is OSI Certified Open Source Software.
; OSI Certified is a certification mark of the Open Source Initiative.
;
; Copyright © 2000 Fanning Software Consulting.
;
; This software is provided "as-is", without any express or
; implied warranty. In no event will the authors be held liable
; for any damages arising from the use of this software.
;
; Permission is granted to anyone to use this software for any
; purpose, including commercial applications, and to alter it and
; redistribute it freely, subject to the following restrictions:
;
; 1. The origin of this software must not be misrepresented; you must
;    not claim you wrote the original software. If you use this software
;    in a product, an acknowledgment in the product documentation
;    would be appreciated, but is not required.
;
; 2. Altered source versions must be plainly marked as such, and must
;    not be misrepresented as being the original software.
;
; 3. This notice may not be removed or altered from any source distribution.
;
; For more information on Open Source Software, visit the Open Source
; web site: http://www.opensource.org.
;
;###########################################################################


PRO CIndex_Colors, event

Widget_Control, event.top, Get_UValue=info, /No_Copy

   ; What kind of event is this?

thisEvent = Tag_Names(event, /Structure_Name)
CASE thisEvent OF

   'WIDGET_BUTTON': XColors, Group_Leader=event.top, $
      NotifyID=[event.id, event.top]

   'XCOLORS_LOAD': BEGIN
      Device, Get_Visual_Depth=thisDepth
      IF thisDepth GT 8 THEN BEGIN
         thisID = !D.Window
         WSet, info.wid
         TV, info.snap
         WSet, thisID
      ENDIF

      ENDCASE
ENDCASE
Widget_Control, event.top, Set_UValue=info, /No_Copy
END
;---------------------------------------------------------------------



PRO CIndex_Event, event
IF event.type NE 1 THEN RETURN
Widget_Control, event.top, Get_UValue=info, /No_Copy
thisID = !D.Window
WSet, info.wid

   ; Use TVIMAGE if you can. If not, use TV.

Catch, theError
IF theError NE 0 THEN BEGIN
   Catch, /Cancel
   Device, Decomposed=0
   TV, info.snap
   GOTO, skip_tvimage
ENDIF

TVImage, info.snap

skip_tvimage:

WSet, thisID
Widget_Control, event.top, Set_UValue=info, /No_Copy
END

PRO CIndex
oldWindowID = !D.Window
thisDevice = !D.Name
Set_Plot, 'Z'
Device, Set_Resolution=[496,400]

   ; Set the starting index for the polygons.

xindex = 0
yindex = 0

   ; Start drawing. There are 16 rows and 16 columns of colors.

FOR i=0,15 DO BEGIN

    y = [yindex, yindex+25, yindex+25, yindex, yindex]
    yindex = yindex+25
    xindex = 0

    FOR j=0,15 DO BEGIN

        x = [xindex, xindex, xindex+31, xindex+31, xindex]
        color = j+(i*16)

           ; Draw the polygon in a specfic color.

        ;Polyfill, x, y, /Device, Color=color
        TV, Replicate(color,31,25), j*31, i*25 ; To avoid bug in Z-buffer.
        output = StrTrim(j+(i*16), 2)

           ; Draw the index number in the "opposite" color.

        XYOutS, xindex+8, yindex-15, output, Color=Byte(255-color), $
           /Device, Charsize=0.75

           ; Reset the xindex number.

        xindex = xindex+31

    ENDFOR

ENDFOR

   ; Take a snapshot of the Z-Buffer.

snap = TVRD()

Set_Plot, thisDevice

tlb = Widget_Base(Title='Color Index Numbers', $
   TLB_Frame_Attr=1, MBar=menuID)
controlID = Widget_Button(menuID, Value='Change Colors')
xcolorsID = Widget_Button(controlID, Value='XColors', $
   Event_Pro='CIndex_Colors')
drawID = Widget_Draw(tlb, XSize=496, YSize=400, /Button_Events)

   ; Center the top-level base.

Device, Get_Screen_Size=screenSize
xCenter = screenSize(0) / 2
yCenter = screenSize(1) / 2

geom = Widget_Info(tlb, /Geometry)
xHalfSize = geom.Scr_XSize / 2
yHalfSize = geom.Scr_YSize / 2

Widget_Control, tlb, XOffset = xCenter-xHalfSize, $
   YOffset = yCenter-yHalfSize

   ; Realize the top-level base.

Widget_Control, tlb, /Realize
Widget_Control, drawID, Get_Value=wid
WSet, wid

   ; Use TVIMAGE if you can. If not, use TV.

Catch, theError
IF theError NE 0 THEN BEGIN
   Device, Decomposed=0
   TV, snap
   GOTO, skip_tvimage
ENDIF

TVImage, snap

Catch, /Cancel
Skip_TVImage:

info = {snap:snap, wid:wid}
Widget_Control, tlb, Set_UValue=info, /No_Copy
WSet, oldWindowID
XManager, 'cindex', tlb, /No_Block
END
