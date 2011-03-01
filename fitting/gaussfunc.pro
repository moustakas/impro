PRO GAUSSFUNC,X,A,YFIT,PDER

; Part of FITAGAUSS
; Evaluates a Gausian and its partial derivative for FITAGAUSS

NUMPTS = N_ELEMENTS(X)
NTERMS = N_ELEMENTS(A)

; Make sure the width is > 0
IF A(2) LT 0. THEN A(2) = -A(2)
IF ABS(A(2)) LT .000001 THEN A(2)=.000001

; Evaluate the function:
DEL = X-A(1)
YFIT=A(0)*EXP( -.5*(DEL/A(2))^2 ) 
IF NTERMS EQ 4 THEN YFIT = YFIT + A(3)

; If requested, provide the partial derivatives:
IF N_PARAMS() EQ 4 THEN BEGIN
   PDER   = FLTARR(NUMPTS,NTERMS)
   IF NTERMS EQ 3 THEN BEGIN
      PDER(*,0) = YFIT/A(0)
      PDER(*,1) = YFIT*DEL/A(2)^2
      PDER(*,2) = PDER(*,1)*DEL/A(2)
   ENDIF ELSE BEGIN
      UFIT = YFIT-A(3)
      PDER(*,0) = UFIT/A(0)
      PDER(*,1) = UFIT*DEL/A(2)^2
      PDER(*,2) = PDER(*,1)*DEL/A(2)
      PDER(*,3) = 1.
   ENDELSE
ENDIF

RETURN
END
