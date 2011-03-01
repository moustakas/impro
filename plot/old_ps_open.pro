Pro   PS_Open, filename, $
		PORTRAIT=portrait, $
		COLOR=color, $
		THICK=thick, $
		PS_FONTS=ps_fonts, $
		ENCAPSULATED=encapsulated, $
		ADOBESTANDARD=adobe, $
		HELP=help

SccsId = '@(#)ps_open.pro 25 Jan 1995 2.4 Fen Tamanaha'
;+
; NAME:
;	PS_OPEN
;
; PURPOSE:
;	This procedure switches the current output device to PostScript.
;
; CATEGORY:
;	Device management.
;
; CALLING SEQUENCE:
;	PS_OPEN [,FName]
;
; OPTIONAL INPUTS:
;	FName:	File name of output PostScript file.  Either '.ps' or '.eps'
;		will automatically be appended.  If not present 'idl.ps' or
;		'idl.eps' is used.
;
; KEYWORD PARAMETERS:
;	/HELP:		Use this switch so see the USAGE line with all the
;			possible parameters and keywords.
;			
;	/PORTRAIT:	Forces portrait orientation.  Landscape is the
;			default.
;
;	/COLOR:		Set this switch if you want to use ANY other color
;			table besides the default black to white (low to
;			high values).  The current color table will be
;			copied to the PostScript device.  The color table
;			will be interpolated to fill the color table of
;			the PostScript device (typically 256).  The system
;			variables !P.Color and !P.Background are set to the
;			darkest and lightest colors in the table,
;			respectively.  Reset these, if desired.  Other color
;			tables can be loaded after calling PS_OPEN.  If you
;			load another color table, be sure to reset the
;			!P.Color system variable to color table entry to
;			to use in PLOT annotations.
;
;	THICK=		Set this keyword to the value by which all line
;			thicknesses system variables will be multiplied.
;			For full page camera ready figures a value of
;			THICK=4 or 5 is reasonable.
;			THICK increase will not override any explicit line
;			thickness definitions like "PLOT, [0,1], THICK=2".
;			To use the THICK keyword effectively code,
;			"PLOT, [0,1], THICK=2*!P.THICK", instead.
;
;	/PS_FONTS:	Set this switch if you want to use the PostScript
;			hardware fonts in plot annotations.  This is identical
;			to entering "!P.FONTS=0".  It's just a little easier.
;			ISOLatin1 versions of the PostScript fonts are used.
;
;	/ENCAPSULATED:	Set this switch is you want to generate an Encapsulated
;			PostScript file.  The file extension will be '.eps'
;			rather than the default '.ps'.  Since this plot is
;			assumed to become part of another document, no rotation
;			of the figure will occur even in the default landscape
;			mode.  Portrait and landscape aspect rations will be
;			maintained.
;
;	/ADOBESTANDARD	Set this switch to use the standard Adobe fonts.
;			By default, PS_OPEN uses ISOLatin1 Adobe fonts which
;			are similar to the ISOLatin1 Hershey font !3 (Simplex
;			Roman).  For example, the ISOLatin1 encoding for
;			the Angstrom symbol is STRING(197B), and plus-or-
;			minus is STRING(177B).  See any ISOLatin1 chart for
;			other symbols.
;
; COMMON BLOCKS:
;	fen_ps_common:	The common block hold the PostScript file name and
;			variables for the original device.  It is shared with
;			PS_CLOSE so that the environment of the original
;			device and its variables can be restored and the
;			output file closed.
;
; SIDE EFFECTS:
;	The plotting device is changed to PostScript until PS_CLOSE is
;	called.
;
;	Various system variables are changed in the PostScript environment,
;	but all the ones from the previous plotting environment are saved and
;	will be restored by PS_CLOSE.
;
;	A 1 inch border is set in both landscape and portrait orientations.
;
;	The BITS_PER_PIXEL is set to 8 to give the maximum number of colors
;	(and thus the largest possible file size).  This is only of importance
;	if you're using a color table.
;
;	Since PostScript plots black ink on white paper, the system variables
;	!P.BACKGROUND and !P.FOREGROUND are set to the lightest and darkest
;	entries in the current color table if /COLOR is set.  Otherwise the
;	lightest and darkest entries are obtained from the default color
;	table.  Changing !P.BACKGROUND does NOT change the actual PostScript
;	background.  This is a bug/feature of IDL.  But, the variable is set
;	so that any routines that wish to plot in the background (white)
;	color can do so by accessing this system variable.
;
; RESTRICTIONS:
;	PS_OPEN and PS_CLOSE are complimentary functions.  Always use PS_CLOSE
;	to end a PS_OPEN session regardless of whether you want or don't want
;	the output.  PS_CLOSE does not automatically send output to the printer,
;	so you don't have to worry about that happening.
;
; PROCEDURE:
;	Straight forward.  :-)
;
; EXAMPLE:
;	PS_OPEN, 'empty', /PS_Fonts, /Port, Thick=5
;	PLOT, [0,1], /NoData, XTitle='!15Day Number', YTitle='Hour Angle'
;	PLOT, day, test1	 		;normal thickness (5)
;	PLOT, day, test2, Thick=0.5*!P.Thick	;half as thick (2.5)
;	PS_CLOSE
;
;	LOADCT, 5
;	PS_OPEN, 'im22', /Port, /Color
;	image(Where(image EQ 0)) = !P.Background    ;map low values to white
;	TV, image
;	PS_CLOSE
;
; MODIFICATION HISTORY:
; 	Written by:	Fen Tamanaha, 07/09/93 Release 2.1.
;	940804 Fen - ISOLatin1 encoding made the default, and the
;			/AdobeStandard keyword added.
;	940805 Fen - Cleaned up handling of color tables.
;	950125 Fen - Fixed /Portrait keyword.
;-

    Common 	fen_ps_common, orig_device, ps_file, $
				orig_p, orig_x, orig_y, orig_z

    On_Error, 2

;
; Check parameters.
;
    If ( Keyword_Set(help) ) Then Begin
	Message, 'Usage: PS_Open [,FName] [,/PORTRAIT] [,/COLOR] [,THICK=]', /Info
	Message, '               [,/PS_FONTS] [,/ENCAPSULATED] [,/ADOBESTANDARD]', /Info
	Return
    EndIf

;
; Initialize the common block if UNDEFINED.
;
    dev_sz = Size(orig_device)
    If ( dev_sz(dev_sz(0)+1) NE 7 ) Then orig_device = !D.Name

    file_sz = Size(ps_file)
    If ( file_sz(file_sz(0)+1) NE 7 ) Then ps_file = ''

;
; Set suffix to '.ps' unless /ENCAPSULATED is set, then set suffix to '.eps'.
;
    If ( Keyword_Set(encapsulated) ) Then Begin
	suffix = '.eps'
    EndIf Else Begin
	suffix = '.ps'
    EndElse

;
; Check if a plot file is already open.
;
    If ( N_Elements(filename) GT 0 ) Then Begin
	fn = filename + suffix
    EndIf Else Begin
	fn = 'idl' + suffix
    EndElse
    If ( ps_file NE '' ) Then Begin
	msg = 'Plot file ' + ps_file + ' is already open.'
	Message, msg, /Cont
	Message, 'Use PS_Close to close the file.', /Info
	Return
    EndIf

;
; Load common block.
;
    ps_file = fn
    orig_device = !D.Name
    orig_p = !P
    orig_x = !X
    orig_y = !Y
    orig_z = !Z

;
; Change graphics device to a PostScript file.
;
; If the /Color keyword is not specified, no color map is placed in the
;	PostScript file and the hardware color table (black to white)
;	is used.  To use your own color map (even grayscale one) the
;	/Color keyword must be set.  The current color table is
;	interpolated into the PostScript color table by default.
;
    If ( Keyword_Set(color) ) Then Begin
;Set_Plot, 'PS';, /Interpolate ; jm01dec18uofa
	Set_Plot, 'PS', /Interpolate
    EndIf Else Begin
	Set_Plot, 'PS'
    EndElse
    Device, File=ps_file
    Device, Color=Keyword_Set(color)

;
; The bits per pixel in
;	the color map are maximized at 8.  This avoids the image size
;	restrictions discussed in the IDL Reference manual.
; The default encoding is set to ISOLatin1 unless the keyword
;	/ADOBESTANDARD is set.
;
    Device, Bits_Per_Pixel=8
    Device, ISOLatin1=(Not Keyword_Set(adobe))

;
; Use encapsulated PostScript format if the /ENCAPSULATED keyword has
;	been set.
;
    Device, Encapsulated=Keyword_Set(encapsulated)

;
; Increase thickness if /Thick keyword set.  Default values of the
;	thickness <= 1 are assumed to be 1.
;
    If ( Keyword_Set(thick) ) Then Begin
	!P.Thick = (!P.Thick > 1) * thick
	!P.CharThick = (!P.CharThick > 1) * thick
	!X.Thick = (!X.Thick > 1) * thick
	!Y.Thick = (!Y.Thick > 1) * thick
	!Z.Thick = (!Z.Thick > 1) * thick
    EndIf

;
; If the /PS_FONTS keyword is set then use hardware/PostScript fonts.
;
    If ( Keyword_Set(ps_fonts) ) Then !P.Font = 0

;
; SET_PLOT automatically sets !P.Color to 0 for a PostScript Device.
;	!P.Background is ignored in PostScript printing.
; If the color map is to be copied to the PostScript output, set the
;	background and foreground colors to the lightest and the darkest.
;	The background is set the lightest to be most like the white
;	paper and the foreground is set to the darkest so it shows up
;	well on the white paper.  These are only good for the current
;	color table.  If you load your own color table after running
;	PS_OPEN, you must reset these values yourself.
;
    If ( Keyword_Set(color) ) Then Begin
	TVLCT, r, g, b, /Get
	Color_Convert, r, g, b, h, l, s, /RGB_HLS	;hue/lightness/sat.
	ct_min = Min(l, min_index)
	ct_max = Max(l, max_index)
	!P.Background = max_index			;closest to white
	!P.Color = min_index				;closest to black
    EndIf

;
; In portrait orientation the page margins are set to 1.0 inch.  In 
;	landscape orientation the page margins are set to 1.0 inch.
;	For landscape, first the box defined by the origin and size
;	is moved and scaled, then it's rotated counterclockwise by 270
;	degrees about the origin.  This is why the YOffSet must be
;	so large.
;
    If ( Keyword_Set(portrait) ) Then Begin
	xoff = 1.0
	yoff = 1.0
	xsize = 6.5
	ysize = 9.0
	Device, /Portrait
    EndIf Else Begin		;/Landscape is the default
	xoff = 1.0
	yoff = 10.0
	xsize = 9.0
	ysize = 6.5
	Device, /Landscape
    EndElse

;
; If the output is encapsulated PostScript, then the offsets are
;	unimportant.  The inclusion macros used to insert the
;	figure into a '.dvi' file will handle the scaling.  The
;	assumption is made that even a figure with a landscape
;	aspect ratio will be maintained but scaled down.
; For rotated output, either print the image out normally or use
;	the macro flags and switches to perform the translation
;	and rotation.  This is painful, but less commonly needed.
;
    If ( Keyword_Set(encapsulated) ) Then Begin
	xoff = 0.0
	yoff = 0.0
	Device, /Portrait
    EndIf

    Device, /Inches, XOffSet=xoff, YOffSet=yoff, XSize=xsize, YSize=ysize

    Return
End
