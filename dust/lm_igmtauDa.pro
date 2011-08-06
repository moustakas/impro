;+
; NAME:
;   LM_IGMTAUDA()
;
; PURPOSE:
;   Compute IGM depression in the Lyman alpha forest at 1050-1170. 
;
; INPUTS:
;   z  - redshift
;
; KEYWORD PARAMETERS:
;
; OUTPUTS:
;   DA = lm_igmtauda(z)
;
; MODIFICATION HISTORY:
;   Leonidas A. Moustakas, 2003 March
;-

FUNCTION igmtauDaint, ww
   return,exp(-3.6d-3 * (ww/1216.d0)^3.46)
END

FUNCTION lm_igmtauDa, z

   zl1 = 1050.*(1.+z)
   zl2 = 1170.*(1.+z)

return,1.-1./(120.*(1.+z))*qromb('igmtauDaint',zl1,zl2,/double)
END 
