Pro	PS_Close, $
		PRINT=prt, $
		HELP=help

SccsId = '@(#)ps_close.pro 28 Apr 1994 2.2 Fen Tamanaha'
;+
; NAME:
;	PS_CLOSE
;
; PURPOSE:
;	This procedure the PostScript output file and restore the output
;	device to its value before invocation of PS_OPEN.
;
; CATEGORY:
;	Device management.
;
; CALLING SEQUENCE:
;	PS_CLOSE
;
; KEYWORD PARAMETERS:
;	/HELP:		Use this switch so see the USAGE line with all the
;			possible parameters and keywords.
;			
;	/PRINT:		Set this switch if you want the PostScript file to
;			be closed and sent the printer named 'PostScript'.
;
; COMMON BLOCKS:
;	fen_ps_common:	The common block hold the PostScript file name and
;			variables for the original device.  It is shared with
;			PS_CLOSE so that the environment of the original
;			device and its variables can be restored and the
;			output file closed.
;
; SIDE EFFECTS:
;	The plotting device is changed from PostScirpt back to its pre-
;	PS_OPEN destination.
;
;	System variables changed in the PostScript environment, are
;	restored by PS_CLOSE.
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
;	PLOT, [0,1], /NoData, Title='!15', $
;			XTitle='!15Day Number', YTitle='Hour Angle'
;	PLOT, day, test1	 		;normal thickness (5)
;	PLOT, day, test2, Thick=0.5*!P.Thick	;half as thick (2.5)
;	PS_CLOSE
;
;	LOADCT, 5
;	PS_OPEN, 'im22', /Port, Color
;	image(Where(image EQ 0)) = !P.Background    ;map low values to white
;	TV, image
;	PS_CLOSE, /Print
;
; MODIFICATION HISTORY:
; 	Written by:	Fen Tamanaha, July 10, 1993  Release 2.1.
;	April 28, 1994	Fen: (2.2) Use "lp" rather than "lpr" to print.
;-

    Common 	fen_ps_common, orig_device, ps_file, $
				orig_p, orig_x, orig_y, orig_z

    On_Error, 2

;
; Print usage if /HELP set.
;
    If ( Keyword_Set(help) ) Then Begin
	Message, 'Usage: PS_Close [,/PRINT]', /Info
	Return
    EndIf

;
; Close the output file not matter what.  This should be called after
;	checking for an open file.  It is safer to try closing the file
;	first since errors in the file checking code might cause an
;	error and the file might remain open.
;
    If ( StrUpCase(!D.Name) EQ 'PS' ) Then Begin
	Device, /Close
    EndIf Else Begin
	Message, 'The current device is not "PS" (PostScript).', /Cont
	Return
    EndElse

;
; Check if a file is open.
;
    If ( ps_file EQ '' ) Then Begin
	Message, 'No PostScript file is currently opened.', /Cont
	Return
    EndIf

    If ( Keyword_Set(prt) ) Then Begin
	cmd = 'lp ' + ps_file
	Spawn, cmd, /Sh
	Print, ps_file + ' closed and sent to the default printer.'
    EndIf Else Begin
	Print, ps_file + ' closed.'
    EndElse

;
; Reset common block and plot settings.
;
    ps_file = ''
    Set_Plot, orig_device		;restore original plotting device
    !P = orig_p
    !X = orig_x
    !Y = orig_y
    !Z = orig_z
   
    Return
End
