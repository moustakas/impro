;+
; NAME: ARM_SWAP
;       
; CATEGORY: miscellaneous
;
; PURPOSE: swap the values of two variables
;
; CALLING SEQUENCE: ARM_SWAP, variable1, variable2
;
; INPUTS: 
;    variable1 - variable to be swapped with variable2
;    variable2 - variable to be swapped with variable1         
;
; EXAMPLE:
;    IDL> x = 3
;    IDL> y = ['a', 'b', 'c']
;    IDL> ARM_SWAP, x, y
;
;    x -> a b c, y -> 3
;
; PROCEDURES: ARM_ERROR
;
; COMMENTS: very simple!
;
; MODIFICATION HISTORY:
;    written by A.R.Marble, Steward Obs., September 2003
;-

pro ARM_SWAP, variable1, variable2

    if variable1 eq variable2 then $
      ARM_ERROR, 'trying to swap identical variables!'

    temp = variable1
    variable1 = variable2
    variable2 = temp
    
    return
    
 end
