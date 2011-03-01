PRO BIN2D,NX,NY, X0,X1, Y0,Y1, X,Y,DATA, MAP,NUM, minbinnum,$
                   weights, tot=tot, mask=mask, notused = notused
;+
; NAME:
;	BIN2D (modified version of ROBUST_BIN2D) MHS 4/2002, CAT 10/2003
;
; PURPOSE:
;	Bins data in a 2-dimensional array suitable for surface or contour 
;	plots.
;
; CALLING SEQUENCE:
;	BIN2D,NX,NY,X0,X1,Y0,Y1,X,Y,DATA,Z,NUM,minbinnum,weights,/tot,
;             notused = notused
;
; INPUTS:
;	NX = Number of bins in the X direction
;	NY = Number of bins in the Y direction
;	X0 = Lower bound on X
;	X1 = Higher bound on X
;	Y0 = Lower bound on Y
;	Y1 = Higher bound on Y
;	X = Independent variable vector, floating-point
;	Y = Dependent variable vector
;	DATA = Vector of quantity to be binned. X, Y, DATA must have same 
;		length.
;       minbinnum = minimum number a data points needed to qualify for a bin 
;       weights = for weigthed average
;       tot = causes data to be totaled instead of averaged
;
; OUTPUTS:
;	MAP   = 2-dimensional array containing the average Z per bin
;	NUM  = 2-dimensional array containing the number of entries per bin.
;       MASK = NaN where number in bin is zero. Equals zero otherwise
;
; OPTIONAL OUTPUTS:
;       NOTUSED = index values of data points not included in the bins
;                 (added by CAT)
;
; SUBROUTINES CALLED:
;
;
; REVISION HISTORY:
;	Written by H.T. Freudenreich, ?/1990
;	Modified: H.F., 3/94 to return SIGGMA
;       modified MHS 4/2002
;- 

ON_ERROR,2

BINMAX = 2500
A    = FLTARR(BINMAX)
Z    = FLTARR(BINMAX,NX,NY)
zindex = FLTARR(BINMAX,NX,NY)
mask = FLTARR(NX,NY)
NUM  = INTARR(NX,NY)
MAP  = FLTARR(NX,NY)

IF N_PARAMS() GT 12 THEN BEGIN
   doweight=1
   W  = FLTARR(BINMAX,NX,NY)
ENDIF ELSE begin 
   doweight=0
   WEIGHTS=DATA *0. + 1.
ENDELSE

; Bin the data:
DELX = float(X1-X0)/NX   
DELY = float(Y1-Y0)/NY   

notused = 0

FOR I = 0L, N_ELEMENTS(DATA)-1 DO BEGIN
  IX =  (X(I)-X0)/DELX
  IY =  (Y(I)-Y0)/DELY 
  IF( (IX GE 0) AND (IY GE 0) AND (IX LT NX) AND (IY LT NY) )THEN BEGIN
     IF( NUM(IX,IY) LT BINMAX )THEN BEGIN
       Z(NUM(IX,IY),IX,IY) = DATA(I)
       zindex(NUM(IX,IY),IX,IY) = I
       IF doweight THEN W(NUM(IX,IY),IX,IY) = WEIGHTS(I)
       NUM(IX,IY)          = NUM(IX,IY)+1
     ENDIF ELSE PRINT,'BIN2D: Exceeded maximum entries per bin'
  ENDIF
ENDFOR

; Average or total each non-empty bin:
FOR I = 0, NX-1 DO BEGIN
  FOR J = 0, NY-1 DO BEGIN
    ;MAP(I,J) = Z(I,J,0)
    MAP(I,J) = Z(0,I,J) ;mhs
    
    IF( NUM(I,J) GE (minbinnum > 1) )THEN BEGIN

      MAP(I,J) = avg(Z[0: NUM(I,J)-1,I,J]) ;mhs

      IF doweight EQ 1 THEN BEGIN 
        MAP(I,J) = total(Z[0: NUM(I,J)-1,I,J] * W[0: NUM(I,J)-1,I,J] ) / $
                   total(W[0: NUM(I,J)-1,I,J])
      ENDIF

      IF keyword_set(tot) THEN MAP(I,J) = total(Z[0: NUM(I,J)-1,I,J])
    ENDIF ELSE BEGIN
      if num(I,J) ge 1 then notused = [notused, zindex[0:NUM(I,J)-1,I,J]]
      MAP(I,J) = 0
      num(I,J) = 0
      mask(I,J) = !VALUES.F_NAN
    ENDELSE
 
  ENDFOR
ENDFOR

if n_elements(notused) gt 1 then notused = notused[1:*]

RETURN
END
