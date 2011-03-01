FUNCTION igmtauDaint, ww
   return,exp(-3.6d-3 * (ww/1216.d0)^3.46)
END

FUNCTION lm_igmtauDa, z
;+
; igmtauDa, z
; 
; returns the IGM depression in the Lyman alpha forest at 1050-1170
; for redshift z
;-
   zl1 = 1050.*(1.+z)
   zl2 = 1170.*(1.+z)
   return,1.-1./(120.*(1.+z))*qromb('igmtauDaint',zl1,zl2,/double)
END 
