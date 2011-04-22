
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; 
;; Routines for computing the Alcock-Paczynski distortion
;;
;; 1998/10  Jonathan Baker <jbaker@astro.berkeley.edu>
;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


; H(z) in units where H(0) = 1

FUNCTION hubble, z

   COMMON cosmol, OmegaM, OmegaL, OmegaQ, wQ

   x = 1. + z
   OmegaK = 1. - OmegaM - OmegaL - OmegaQ
   h2 = OmegaM * x^3 + OmegaK * x^2 + OmegaL + OmegaQ * x^(3.+3.*wQ)
   
   RETURN, sqrt(h2)

END


; Compute d(chi)/dz where d(chi)^2 = dr^2 / (1-kr^2)

FUNCTION dchi_dz, z

   RETURN, 1./hubble(z)

END


; Compute coordinate distance r(z) for FRW metric

FUNCTION rz, z, $
            OmegaMatter = OmegaM, $
            OmegaLambda = OmegaL, $
            OmegaQ = OmegaQ, $
            wQ = wQ, $
            eps = eps

   COMMON cosmol, cOmegaM, cOmegaL, cOmegaQ, cwQ
   
   IF NOT keyword_set(OmegaM) THEN OmegaM = 1.d0
   IF NOT keyword_set(OmegaL) THEN OmegaL = 0.d0
   IF NOT keyword_set(OmegaQ) THEN OmegaQ = 0.d0
   IF NOT keyword_set(wQ) THEN wQ = 0.d0
   
   cOmegaM = OmegaM
   cOmegaL = OmegaL
   cOmegaQ = OmegaQ
   cwQ = wQ
   
   kurv = OmegaM + OmegaL + OmegaQ - 1.d0
   
   nz = n_elements(z)
   dchi = dblarr(nz)
   chi = dblarr(nz)
   r = dblarr(nz)
   
   z2 = z[0]
   IF z2 EQ 0. THEN $
     dchi[0] = 0.d0 $
   ELSE $
     dchi[0] = qromb('dchi_dz', 0., z2, /double, eps=eps)
   
   FOR i = 1, nz-1 do begin
      z1 = z[i-1]
      z2 = z[i]
      dchi[i] = qromb('dchi_dz', z1, z2, /double, eps=eps)
   ENDFOR
   
   chi[0] = dchi[0]
   FOR i = 1, nz-1 do $
     chi[i] = chi[i-1] + dchi[i]
   
   IF abs(kurv) LT 1.e-4 THEN $ ; flat
     r = chi $
   ELSE IF kurv GT 0.d0 THEN $  ; closed
     r = sin(chi*sqrt(kurv))/sqrt(kurv) $
   ELSE $                       ; open
     r = sinh(chi*sqrt(-kurv))/sqrt(-kurv)
   
   RETURN, r
   
END

