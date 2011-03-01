 FUNCTION MED,A       
;+
; NAME
;	 MED
; PURPOSE
;	Compute the median of an array, which may be of even length. 
;
; CALLING SEQUENCE:  
; 	MID_VALUE = MED(A)
;
; OUTPUTS 
; 	The median of array A
; REVISION HISTORY:
;	H.T. Freudenreich, ?/89
;-
ON_ERROR,2
NUM = N_ELEMENTS(A)
IF NUM MOD 2 EQ 0 THEN BEGIN  ; even # points. Can't call MEDIAN.
   B = A   &  B = B( SORT(B) )  &  I0  = (NUM-1)/2  & MED =.5*(B(I0)+B(I0+1)) 
ENDIF ELSE MED=MEDIAN(A)
RETURN, MED
END
