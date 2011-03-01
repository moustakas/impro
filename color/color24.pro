;+
; NAME:
;       COLOR24
;
; PURPOSE:
;
;       The purpose of this function is to convert a RGB color triple
;       into the equivalent 24-bit long integer. The 24-bit integer
;       can be decomposed into the appropriate color by interpreting
;       the lowest 8 bits as red, the middle 8 bits as green, and the
;       highest 8 bits as blue.
;
; AUTHOR:
;
;       FANNING SOFTWARE CONSULTING
;       David Fanning, Ph.D.
;       1645 Sheely Drive
;       Fort Collins, CO 80526 USA
;       Phone: 970-221-0438
;       E-mail: davidf@dfanning.com
;       Coyote's Guide to IDL Programming: http://www.dfanning.com/
;
; CATEGORY:
;
;       Graphics, Color Specification.
;
; CALLING SEQUENCE:
;
;       color = COLOR24(rgb_triple)
;
; INPUTS:
;
;       RGB_TRIPLE: A three-element column or row array representing
;       a color triple. Or an N-by-three element array of color triples.
;       The values of the elements must be between 0 and 255.
;
; KEYWORD PARAMETERS:
;       None.
;
; COMMON BLOCKS:
;       None.
;
; SIDE EFFECTS:
;       None.
;
; RESTRICTIONS:
;       None.
;
; EXAMPLE:
;
;       To convert the color triple for the color YELLOW,
;       (255, 255, 0), to the hexadecimal value '00FFFF'x
;       or the decimal number 65535, type:
;
;       color = COLOR24([255, 255, 0])
;
;       This routine was written to be used with device-independent
;       color programs like GETCOLOR.
;
; MODIFICATION HISTORY:
;
;       Written by:  David Fanning, 3 February 96.
;       Completely revised the algorithm to accept color arrays. 19 October 2000. DWF.
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



FUNCTION COLOR24, color

   ; This FUNCTION accepts a [red, green, blue] triple that
   ; describes a particular color and returns a 24-bit long
   ; integer that is equivalent to (can be decomposed into)
   ; that color. The triple can be either a row or column
   ; vector of 3 elements or it can be an N-by-3 array of
   ; color triples.

ON_ERROR, 2

s = Size(color)

IF s[0] EQ 1 THEN BEGIN
   IF s[1] NE 3 THEN Message, 'Input color parameter must be a 3-element vector.'
   RETURN, color[0] + (color[1] * 2L^8) + (color[2] * 2L^16)
ENDIF ELSE BEGIN
   IF s[2] GT 3 THEN Message, 'Input color parameter must be an N-by-3 array.'
   RETURN, color[*,0] + (color[*,1] * 2L^8) + (color[*,2] * 2L^16)
ENDELSE

END ;--------------------------------------------------------------------------------------------
