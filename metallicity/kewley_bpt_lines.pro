;+
; NAME:
;       KEWLEY_BTP_LINES()
;
; PURPOSE:
;       Return the Kewley et al 2001 theoretical starburst
;       classification lines corresponding to the BPT diagrams.  
;
; INPUTS:
;       None.
;
; OPTIONAL INPUTS:
;
; KEYWORD PARAMETERS:
;       kauffmann - use the Kauffmann et al. (2003) dividing line for
;                   AGN
;
; OUTPUTS:
;       models - output data structure
;
; COMMENTS:
;       See p. 184 of L. Kewley's thesis or Eqs. (5-7) in Kewley et al
;       2001, ApJS, 132, 37.
;
; PROCEDURES USED:
;
; EXAMPLE:
;
; MODIFICATION HISTORY:
;       J. Moustakas, 2002 April 19, U of A, written
;       jm03jun14uofa - added upper and lower bounds
;       jm04jul26uofa - to prevent the second hyperbola from being
;                       drawn require NII/Ha < 0.4
;       jm04jul29uofa - added KAUFFMANN keyword
;
; Copyright (C) 2002-2004, John Moustakas
; 
; This program is free software; you can redistribute it and/or modify 
; it under the terms of the GNU General Public License as published by 
; the Free Software Foundation; either version 2 of the License, or
; (at your option) any later version. 
; 
; This program is distributed in the hope that it will be useful, but 
; WITHOUT ANY WARRANTY; without even the implied warranty of
; MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
; General Public License for more details. 
;-

function kewley_bpt_lines, ratio_nii=ratio_nii, ratio_sii=ratio_sii, $
  ratio_oi=ratio_oi, nii=nii, sii=sii, oi=oi, kauffmann=kauffmann

; ---------------------------------------------------------------------------
; NII/H-alpha theoretical classification line    
; ---------------------------------------------------------------------------

    if keyword_set(nii) and (n_elements(ratio_nii) ne 0L) then begin

       if keyword_set(kauffmann) then $
         x_nii = ratio_nii < 0.4 else $ ; <-- NOTE
         x_nii = ratio_nii < 0.4

    endif else begin

       if keyword_set(kauffmann) then $
         x_nii = findgen((0.0-(-2.7))/0.01+1)*0.01+(-2.7) else $ ; <-- NOTE
         x_nii = findgen((0.36-(-2.7))/0.01+1)*0.01+(-2.7)
       
    endelse
       
    if keyword_set(kauffmann) then begin

       y_nii = 0.61 / (x_nii - 0.05) + 1.3
       y_nii_upper = y_nii*0.0
       y_nii_lower = y_nii*0.0

    endif else begin

       y_nii = 0.61/(x_nii-0.47)+1.19
       y_nii_upper = 0.61/(x_nii-0.37)+1.09
       y_nii_lower = 0.61/(x_nii-0.57)+1.29

    endelse

; ---------------------------------------------------------------------------
; SII/H-alpha theoretical classification line    
; ---------------------------------------------------------------------------

    if keyword_set(sii) and (n_elements(ratio_sii) ne 0L) then $
      x_sii = ratio_sii else $
      x_sii = findgen((0.2-(-2.7))/0.01+1)*0.01+(-2.7)
    
    y_sii = 0.72/(x_sii-0.32)+1.30
    y_sii_upper = 0.72/(x_sii-0.22)+1.20
    y_sii_lower = 0.72/(x_sii-0.42)+1.40

; ---------------------------------------------------------------------------
; OI/H-alpha theoretical classification line    
; ---------------------------------------------------------------------------
    
    if keyword_set(oi) and (n_elements(ratio_oi) ne 0L) then $
      x_oi = ratio_oi else $
      x_oi = findgen(((-0.7)-(-3.5))/0.01+1)*0.01+(-3.5)
    
    y_oi = 0.73/(x_oi+0.59)+1.33
    y_oi_upper = 0.73/(x_oi+0.69)+1.23
    y_oi_lower = 0.73/(x_oi+0.49)+1.43

; ---------------------------------------------------------------------------
; fill the data structure and return
; ---------------------------------------------------------------------------
    
    models = {$
      x_nii:       x_nii, $
      y_nii:       y_nii, $
      y_nii_upper: y_nii_upper, $
      y_nii_lower: y_nii_lower, $
      x_sii:       x_sii, $
      y_sii:       y_sii, $
      y_sii_upper: y_sii_upper, $
      y_sii_lower: y_sii_lower, $
      x_oi:        x_oi, $
      y_oi:        y_oi, $
      y_oi_upper:  y_oi_upper, $
      y_oi_lower:  y_oi_lower}
    
return, models
end
