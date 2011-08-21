;+
; NAME:
;       IM_INTEGRAL()
;
; PURPOSE:
;       Routine to perform trapezoidal integration in X,Y between
;       limits xmin to xmax. 
;
; CALLING SEQUENCE:
;	result = integral(x,y,xmin,xmax)
;
; INPUTS:
;	x,y - vectors to be integrated
;	xmin,xmax - vectors with lower and upper integral limits
; OUTPUTS:
;	the integrations between xmin and xmax are returned
;	as the function result
; RESTRICTIONS:
;	The values in XMIN must be greater than or equal to the minimum
;	of X. The values in XMAX must be less than or equal to the 
;	maximum of X. X must be in ascending order.
;
; HISTORY:
;	Version 1,  D. Lindler  (a long time ago)
;	Version 2,  JKF/ACC	28-jan-1992	- moved to IDL V2.
; 	Version 3,  DJL 	17-jun-1992	- FIXED BUG AT ENDPOINTS
;	version 4,  DJL		27-Jul-1992	- fixed previous change to
;						  work with vector
;						  inputs
;
;       J. Moustakas, 2003 November 25, U of A, made double precision,
;          added square brackets, and error checking
;       jm08apr14nyu - additional error checking; default XMIN,XMAX to
;                      min and max of X array
;
;-

FUNCTION IM_INTEGRAL, X, Y, XMIN1, XMAX1

    nx = n_elements(x)
    ny = n_elements(y)
    if (nx eq 0L) or (ny eq 0L) or (nx ne ny) then begin
       doc_library, 'im_integral'
       return, -1L
    endif

    if (n_elements(xmin1) eq 0) then xmin = min(x) else xmin = xmin1>min(x)
    if (n_elements(xmax1) eq 0) then xmax = max(x) else xmax = xmax1<max(x)
    
    TABINV,X,XMIN,RMIN
    TABINV,X,XMAX,RMAX
    n = n_elements(x)
;
; CHECK FOR VALID LIMITS
;
    A=MAX(XMIN)>MAX(XMAX)
    B=MIN(XMIN)<MIN(XMAX)
    D=MIN(XMAX-XMIN)
    IF (A GT MAX(X)) OR (B LT MIN(X)) OR (D LT 0.0) THEN $
      message,'INVALID INTEGRAL LIMITS SUPPLIED TO INTEGRAL FUNCTION'
;
; COMPUTE DIFFERENCES IN X AND Y
;
    DX=SHIFT(X,-1)-X
    DY=SHIFT(Y,-1)-Y
;
; COMPUTE INTEGRALS FOR EACH FULL INTERVAL IN X
;
    DINT=(SHIFT(Y,-1)+Y)/2.0*DX
;
; COMPUTE FULL INTERVALS TO INTEGRATE BETWEEN
;
    IMIN=FIX(RMIN)
    IMAX=FIX(RMAX)
;
; COMPUTE FUNCTION VALUES AT XMIN AND XMAX
;                                                 
    DXMIN=XMIN-X[IMIN]
    YMIN=Y[IMIN]+DXMIN*(Y[IMIN+1]-Y[IMIN])/DX[IMIN]
    DXMAX=XMAX-X[IMAX]
    YMAX=Y[IMAX]+DXMAX*(Y[(IMAX+1L)<(n-1L)] - Y[IMAX])/DX[IMAX]
;
; COMPUTE INTEGRAL FROM IMIN TO IMAX
;
    NOUT=N_ELEMENTS(XMIN)
    INT = xmin*0.0
    FOR I=0L, NOUT-1L DO $
       IF IMAX[I] NE IMIN[I] THEN INT[I]=TOTAL(DINT[IMIN[I]:IMAX[I]-1L])
;
; SUBTRACT INTEGRAL FROM IMIN TO RMIN
;
    INT=INT - (Y[IMIN]+YMIN)/2.*DXMIN
;
; ADD INTEGRAL FROM IMAX TO RMAX
;
    INT=INT + (Y[IMAX]+YMAX)/2.0*DXMAX

RETURN,INT
END
