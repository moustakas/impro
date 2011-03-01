FUNCTION PLANEFIT,XIN,YIN,ZIN,W,  YFIT
;+
; NAME:
;	PLANEFIT
; PURPOSE:
;	Least-squares fit of a plane to a set of (X,Y,Z) points: 
;		Z=R(0)+R(1)*X+R(2)*Y
; CALLING SEQUENCE:
;	coef = PLANEFIT( x,y,z,0., yfit )  ; for equal weighting
;	coef = PLANEFIT( x,y,z, weights_vector, yfit ) ;otherwise
; INPUT:
;	XIN, YIN, ZIN = the points to be fit
;	Weights_vector = their weights. (For equal weights, let W=0.)
; RETURNS:
;	coef = the fit coefficients
; OPTIONAL OUTPUT:
;	YFIT = the Z values at (X,Y), calculated from the fit
;
; NOTE:
;	If a fit is not possible, R will be returned as a scalar, the average 
;	of Z.
; REVISION HISTORY:
;	Written by H.T. Freudenreich, HSTX, 5/18/93
;-

R=FLTARR(3)
EPS=1.0E-24
NUMFAIL=0
WENT_DOUBLE=0

N = N_ELEMENTS(ZIN)
IF N LT 3 THEN BEGIN
   PRINT,'PLANEFIT: too few points!'
   R=AVG(Z)
   RETURN,R
ENDIF

; Shift X, Y, Z to the center of gravity frame:
X0 = TOTAL(Xin)/N & Y0 = TOTAL(Yin)/N & Z0 = TOTAL(Zin)/N
X  = XIN-X0       & Y  = YIN-Y0       & Z  = ZIN-Z0 

TRY_AGAIN:

; Do we have weights?
IF N_ELEMENTS(W) LT 3 THEN BEGIN
   S = N*1.
   SX =TOTAL(X)    &  SY =TOTAL(Y)    &  SZ =TOTAL(Z)
   SXX=TOTAL(X*X)  &  SYY=TOTAL(Y*Y)  
   SXY=TOTAL(X*Y)  &  SXZ=TOTAL(X*Z)  &  SYZ=TOTAL(Y*Z)
ENDIF ELSE BEGIN
   S = TOTAL(W)
   SX =TOTAL(W*X)    &  SY =TOTAL(W*Y)    &  SZ =TOTAL(W*Z)
   SXX=TOTAL(W*X*X)  &  SYY=TOTAL(W*Y*Y)  
   SXY=TOTAL(W*X*Y)  &  SXZ=TOTAL(W*X*Z)  &  SYZ=TOTAL(W*Y*Z)
ENDELSE

D = S*(SXX*SYY-SXY*SXY)+SX*(SXY*SY-SX*SYY)+SY*(SX*SXY-SXX*SY)

IF ABS(D) LT EPS THEN BEGIN
;  If this is the 1st failure, convert the vectors to double precision if they
;  aren't already and try again.
   NUMFAIL = NUMFAIL+1
   IF NUMFAIL EQ 1 THEN BEGIN
      SYZ=SIZE(ZIN)
      IF SYZ(2) NE 5 THEN BEGIN
         X=DOUBLE(X) & Y=DOUBLE(Y) & Z=DOUBLE(Z)
         WENT_DOUBLE = 1
         GOTO,TRY_AGAIN
      ENDIF
   ENDIF
   PRINT,'PLANEFIT: Unable to invert matrix! Taking mean instead.'
   R=TOTAL(Z)/S
   RETURN,R
ENDIF

R(0)=SZ*(SXX*SYY-SXY*SXY) + SX*(SXY*SYZ-SXZ*SYY) + SY*(SXZ*SXY-SXX*SYZ)
R(1)=S *(SXZ*SYY-SXY*SYZ) + SZ*(SXY*SY-SX*SYY)   + SY*(SX*SYZ-SXZ*SY)
R(2)=S *(SXX*SYZ-SXZ*SXY) + SX*(SXZ*SY-SX*SYZ)   + SZ*(SX*SXY-SXX*SY)
R=R/D

; Now shift (X,Y,Z) back:
R(0) = R(0) + Z0 - R(1)*X0 - R(2)*Y0

IF N_PARAMS(0) GT 4 THEN YFIT=R(0)+R(1)*Xin+R(2)*Yin

IF WENT_DOUBLE THEN BEGIN
;  Go back to single precision.
   R=FLOAT(R)
   IF N_PARAMS(0) GT 4 THEN YFIT=FLOAT(YFIT)
ENDIF

RETURN,R
END
