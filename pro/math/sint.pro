;+
; NAME:
;    sint
; PURPOSE:
;    Sinc interpolation of a 1-D vector of data.
; DESCRIPTION:
;
; CATEGORY:
;    Numerical
; CALLING SEQUENCE:
;    result = sint( x, f )
; INPUTS:
;    x  : Independent variable values for which f is to be interpolated.
;         Note: The implied independent variable values for f are the indicies
;         of the vector f.
;    f  : Vector of function values (dependent variable).
;
; OPTIONAL INPUT PARAMETERS:
;
; KEYWORD PARAMETERS:
;
; OUTPUTS:
;    Interpolated value(s).
; COMMON BLOCKS:
;
; SIDE EFFECTS:
;
; RESTRICTIONS:
;
; PROCEDURE:
;
; MODIFICATION HISTORY:
;    Written by Doug Loucks, Lowell Observatory, September, 1993.
;    Adapted from the IDL function sshift.pro written by Marc Buie.
;    01/14/94, DWL, Documentation update.
;-
FUNCTION sint, x, f

dampfac = 3.25
ksize   = 21

nx = N_ELEMENTS( x )
nf = N_ELEMENTS( f )

ix = FIX( x )
fx = x - ix
z = WHERE( fx EQ 0, countz )
i = WHERE( fx NE 0, counti )

r = x * 0

IF countz NE 0 THEN BEGIN
   ;There are integral values of the independent variable. Select the function
   ;values directly for these points.
   r[ z ] = f[ ix[z] ]
ENDIF

IF counti NE 0 THEN BEGIN
   ;Use sinc interpolation for the points having fractional values.
   FOR point=0, counti-1 DO BEGIN
      xkernel = ( FINDGEN( ksize ) - 10 ) - fx[ i[ point ] ]
      u1 = xkernel / dampfac
      u2 = !pi * xkernel
      sinc = EXP( -( u1*u1 ) ) * SIN( u2 ) / u2
      lobe = ( INDGEN( ksize ) - 10 ) + ix[ i[point] ]
      vals = FLTARR( ksize )
      w = WHERE( lobe LT 0, count )
      IF count NE 0 THEN vals[w] = f[0]
      w = WHERE( lobe GE 0 AND lobe LT nf, count )
      IF count NE 0 THEN vals[w] = f[ lobe[w] ]
      w = WHERE( lobe GE nf, count )
      IF count NE 0 THEN vals[w] = f[ nf-1 ]
      r[ i[ point ] ] = TOTAL( sinc * vals )
   ENDFOR
ENDIF

RETURN, r

END
