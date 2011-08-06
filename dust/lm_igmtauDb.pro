;+
; NAME:
;   LM_IGMTAUDB()
;
; PURPOSE:
;   Compute IGM depression in the Lyman beta forest at 920-1015. 
;
; INPUTS:
;   z  - redshift
;
; KEYWORD PARAMETERS:
;
; OUTPUTS:
;   DB = lm_igmtauda(z)
;
; MODIFICATION HISTORY:
;   Leonidas A. Moustakas, 2003 March
;-

FUNCTION igmtauDbint, ww
   lyser = [1216.,1026.,973.,950.]*1.d0
   Alyser = [3.6e-3,1.7e-3,1.2e-3,9.3e-4]*1.d0
   tmptau = 0.d0
   FOR ii=1,n_elements(lyser)-1 DO tmptau = tmptau+Alyser[ii]*(ww/lyser[ii])^3.46
   return,exp(-tmptau)
END

FUNCTION lm_igmtauDb, z
   zl1 = 920.*(1.+z)
   zl2 = 1015.*(1.+z)
return,1.-1./(95.*(1.+z))*qromb('igmtauDbint',zl1,zl2,/double)
END 
