; Copyright (c) 1998, Forschungszentrum Juelich GmbH ICG-3
; All rights reserved.
; Unauthorized reproduction prohibited.
;
;+
; NAME:
;   is_string
;
; PURPOSE:
;   Test if variable is of type string
;
; CATEGORY:
;   PROG_TOOLS
;
; CALLING SEQUENCE:
;   Result = is_string(var)
;
; INPUTS:
;   var: The variable to be tested.
;
; OUTPUTS:
;   This function returns 1 if var is of type string else it returns 0
;
;
; KEYWORD PARAMETERS:
;   SCALAR: if this keword is set the function returns 1 only, if the
;           input var is a scalar string
;
; EXAMPLE:
;   A = 'Hallo'
;   PRINT, is_string(A)
;   -> 1
;
; MODIFICATION HISTORY:
;   Written by: Frank Holland, 19.01.98
;   27.01.98: debug statement auskommentiert
;   07.09.98: Return statement vereinfacht
;   15.11.00: T.B. scalar keyword
;-

FUNCTION is_string, var, SCALAR=sc
    a = SIZE(var)
    n = N_ELEMENTS(a)
    IF KEYWORD_SET(sc) THEN $
        RETURN, (a[n - 2] EQ 7) AND (a[n-1] EQ 1)
    RETURN, a[n - 2] EQ 7
END