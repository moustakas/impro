;+
; NAME: ARM_REARRANGE
;       
; CATEGORY: miscellaneous
;
; PURPOSE: apply new ordering to multiple arrays
;
; CALLING SEQUENCE: ARM_REARRANGE, a1, [a2, ... a100, order=order]
;
; INPUTS: 
;   a1    - array of values to rearrange
;   order - desired ordering [eg, output of SORT()]
;       
; OPTIONAL INPUTS:
;   a2-a100 - additional arrays to rearrange
; 
; KEYWORDS:
;
; OUTPUTS: a1-a100 - original arrays replaced by with new ordering
; 
; OPTIONAL OUTPUTS:
;
; EXAMPLE:
;
; PROCEDURES USED:
;
; COMMENTS: 
;
; MODIFICATION HISTORY:
;    written by A.R.Marble, Steward Obs., 2004 Feb 10
;-

pro arm_rearrange, order=order, $
       a1,  a2,  a3,  a4,  a5,  a6,  a7,  a8,  a9,  a10, $
       a11, a12, a13, a14, a15, a16, a17, a18, a19, a20, $
       a21, a22, a23, a24, a25, a26, a27, a28, a29, a30, $
       a31, a32, a33, a34, a35, a36, a37, a38, a39, a40, $
       a41, a42, a43, a44, a45, a46, a47, a48, a49, a50, $
       a51, a52, a53, a54, a55, a56, a57, a58, a59, a60, $
       a61, a62, a63, a64, a65, a66, a67, a68, a69, a70, $
       a71, a72, a73, a74, a75, a76, a77, a78, a79, a80, $
       a81, a82, a83, a84, a85, a86, a87, a88, a89, a90, $
       a91, a92, a93, a94, a95, a96, a97, a98, a99, a100

    n = N_PARAMS()

; error checking

    size = N_ELEMENTS(a1)

    if N_ELEMENTS(order) ne size then MESSAGE, $
      'A1 and ORDER and incompatible dimensions.'

    for i = 1, n do begin

       ndx = STRCOMPRESS(STRING(i), /remove_all)
       dummy = EXECUTE('a = a'+ndx)

       if N_ELEMENTS(a) ne size then $
         MESSAGE, 'A'+ndx+' has incompatible dimensions.'

    endfor

; rearrange arrays

    for i = 1, n do begin
       
       ndx = STRCOMPRESS(STRING(i), /remove_all)
       dummy = EXECUTE('a'+ndx+' = a'+ndx+'[order]')
    
    endfor

 end
