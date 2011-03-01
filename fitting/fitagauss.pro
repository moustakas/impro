FUNCTION FITAGAUSS,X,Y,A,SIGMAA
;+
; NAME:
;	FITAGAUSS
; PURPOSE:
;	A version of CURVEFIT that fits a Gaussian to the data, using the
;	mean absolute deviation rather than the chi^2. DESIGNED ONLY FOR
;	USE BY HISTOGAUSS AND HALFAGAUSS. EMPLOY ELSEWHERE AT USER'S RISK.
; CALLING SEQUNENCE: 
;	value = FITAGAUSS( X, Y, A, SIGMAA )
; INPUT:
;	X, Y = the points to fit
;	A = the vector of parameters: [height,center,width,Y offset]. The
;	offset is optional. If the A input has 3 element, no offset is returned.
;       If A has four elements, the last parameter is a constant offset:
;       y = G(x) + A(3)
; OUTPUT:
;	Value -	The value of the Gaussian fit at the input points.
;	SIGMAA = uncertainties of parameters
; REVISION HISTORY:
;	Written,   ; H.T. Freudenreich, ?/?
;-
	ON_ERROR,2		;RETURN TO CALLER IF ERROR
        W=Y
        W(*)=1.
        EPS=1.0E-8
	A = 1.*A		;MAKE PARAMS FLOATING
	NTERMS = N_ELEMENTS(A)	;# OF PARAMS.
        NPTS  = N_ELEMENTS(Y)
        M0 = (NPTS-1)/2
	NFREE = (N_ELEMENTS(Y)<N_ELEMENTS(X))-NTERMS ;Degs of freedom
	IF NFREE LE 0 THEN STOP,'Curvefit - not enough data points.'
	FLAMBDA = 0.001		;Initial lambda
	DIAG = INDGEN(NTERMS)*(NTERMS+1) ;SUBSCRIPTS OF DIAGONAL ELEMENTS
;
	FOR ITER = 1,30 DO BEGIN	;Iteration loop
;
;		EVALUATE ALPHA AND BETA MATRICES.
;
        GAUSSFUNC,X,A,YFIT,PDER 

	BETA = (Y-YFIT)*W # PDER
	ALPHA = TRANSPOSE(PDER) # (W # (FLTARR(NTERMS)+1)*PDER)

        slant = alpha(diag)
        q = where(abs(slant) lt 1.0e-20)
        if q(0) ge 0 then begin
           slant(q) = 1.0e-20
           alpha(diag) = slant
        endif

        INR=WHERE( ABS(X-A(1)) LE (1.5*A(2)), COUNT )
        DEV=ABS(Y(INR)-YFIT(INR))
        CHISQ1=TOTAL(DEV)/COUNT
;
;	INVERT MODIFIED CURVATURE MATRIX TO FIND NEW PARAMETERS.
;
	REPEAT BEGIN
		C = SQRT(ALPHA(DIAG) # ALPHA(DIAG))
		ARRAY = ALPHA/C
		ARRAY(DIAG) = 1.+FLAMBDA
		ARRAY = INVERT(ARRAY)
		B = A+ ARRAY/C # TRANSPOSE(BETA) ;NEW PARAMS
                GAUSSFUNC,X,B,YFIT 
                INR=WHERE( ABS(X-B(1)) LE (1.5*B(2)), COUNT )
                DEV=ABS(Y(INR)-YFIT(INR))
                CHISQR=TOTAL(DEV)/COUNT
		FLAMBDA = FLAMBDA*10.	;ASSUME FIT GOT WORSE
                CHEESE_DIP = CHISQR-CHISQ1
		ENDREP UNTIL CHEESE_DIP LE EPS
    	        FLAMBDA = FLAMBDA/10.	;DECREASE FLAMBDA BY FACTOR OF 10
	        A=B			;SAVE NEW PARAMETER ESTIMATE.
        	IF ((CHISQ1-CHISQR)/CHISQ1) LE .0010 THEN GOTO,DONE ;Finished?
	ENDFOR			;ITERATION LOOP
;
	PRINT,'FITAGAUSS - Failed to converge'
;
DONE:	
        SIGMAA = SQRT(ARRAY(DIAG)/ALPHA(DIAG)) ;RETURN SIGMA'S
	RETURN,YFIT		;RETURN RESULT
END
