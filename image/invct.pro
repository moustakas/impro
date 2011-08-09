;+
; NAME:
;	INVCT
;
; PURPOSE:
;	Invert the current and original color tables.
;
; CATEGORY:
;	Image display.
;
; CALLING SEQUENCE:
;	INVCT
;
; KEYWORD PARAMETERS:
;	/QUIET:	Set this switch to suppress printing of a confirming message.
;
; COMMON BLOCKS:
;	Colors:	Common block that contains current and original R, G, and B
;		colors.  Used by LOADCT, HSV, and HLS routines.
;
; SIDE EFFECTS:
;	Image display color tables are altered.
;	Original color table in common block is inverted.
;
; RESTRICTIONS:
;	The current color table is read from the device and not from the
;	common block.
;
; PROCEDURE:
;	Straight forward.  :-)
;
; EXAMPLE:
;	LOADCT, 0		;load a black to white color table
;	INVCT			;make it a white to back color table (better for
;				;  PostScript
;	image = DIST(200,200)	;make an image
;	TVSCL, image		;display the image
;
; MODIFICATION HISTORY:
; 	Written by:	Fen Tamanaha, July 10, 1993.  Revision 1.1
;-

Pro	InvCT, $
		QUIET=quiet

SccsId = '@(#)invct.pro 1.1 7/12/93 Fen Tamanaha'

    Common colors, r,g,b,cur_red,cur_green,cur_blue

    On_Error,2

    nc = !D.Table_Size
    If ( nc LE 0 ) Then Begin
	Message, 'Error: This device does not support a color table.'
    EndIf

    If ( N_Elements(r) LE 0) Then Begin
	r = Fix(FIndGen(nc)*256.0/nc)		;load black to white
	g = r
	b = r
    EndIf

    TVLCT, cur_red, cur_green, cur_blue, /Get	;get current color table

    cur_red = Reverse(cur_red)			;invert current color table
    cur_green = Reverse(cur_green)
    cur_blue = Reverse(cur_blue)

    TVLCT, cur_red, cur_green, cur_blue		;load inverted color table

    r = Reverse(r)				;invert original color table
    g = Reverse(g)
    b = Reverse(b)

    If ( Not Keyword_Set(quiet) ) Then Begin
	Message, 'Color table inverted.', /Info
    EndIf

    Return
End
