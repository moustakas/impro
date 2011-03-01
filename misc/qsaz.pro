function ngam, z
     common qsed, q
     q=qsaz()
     q.e = q.e 
     if n_elements(z) eq 0 then z = 5.8
     opz = 1.d0 + z
     Emax = 1.8 * sqrt(opz/15.d0) * 1.d3 ; eV
     sol = 2.9979247d10         ; cm/s
     return,4.*!dpi / sol * qpint1d('FEoverE', 13.6/opz, emax/opz)
 end 
function FEoverE, x
     common qsed, q
     return, interpol(q.e,q.w,x/1d3)/(x/1d3)
 end 

function sfe 
     common qsed, q
     q=qsaz()
     return,qpint1d('myfunc',0.5d,2.0d) / qpint1d('myfunc',1d-5,1d6)
 end
; 
function myfunc, x 
     common qsed, q 
     return, interpol(q.e,q.w,x)
 end 

;-

function qsaz, plot=plot
     common E, energy
     common qsed, q 
     
; Sazonov, ostriker, sunyaev 2004, MNRAS, 347, 144
; equation 8

     alpha = 0.24
     beta  = 1.60
     gamma = 1.06
     E1    = 83.0               ; keV
     k     = 4.1d-3

     E0    = (beta - alpha) * E1
     
     A = 2.^alpha * exp(2./E1)
     B = E0^(beta - alpha) * exp(-(beta-alpha)) / (1+k*E0^(beta-gamma))

;     np = 20000l
     np = 5000l
     emin = alog10(1d-5)        ; 0.01eV = 1e6 A
     emax = alog10(1d+6)        ; 1e6 keV
     Earr = 10^((findgen(np)/(np-1)*(1-emin/emax)+emin/emax)*emax)

     lx = where(Earr ge 2. and Earr lt E0)
     hx = where(Earr gt E0)

     FE = dblarr(n_elements(Earr))
     FE[lx] = A * Earr[lx]^(-alpha) * exp(-Earr[lx]/E1)
     FE[hx] = A * B * (1.+k*Earr[hx]^(beta-gamma)) * Earr[hx]^(-beta)

; equation 14

     smx = where(Earr lt 1d-3)
     slx = where(Earr ge 1d-3 and Earr lt 1d-2)
     shx = where(Earr ge 1d-2 and Earr lt 2.)

     for w=0l,n_elements(smx)-1 do begin
         energy = Earr[w]
;         FE[w] = 1.5d42 * qromb('irwavs',100.d0,1400.d0,/double)
         FE[w] = 1.5d42 * qpint1d('irwavs',100.d0,1400.d0)
     endfor 

     FE[slx] = 159.d0 * Earr[slx]^(-0.6)            ; * 1.2
     FE[shx] = Earr[shx]^(-1.7) * exp(Earr[shx]/2.) ; * 1.2

; normalize to Eddington luminosity for the entire SED
     G      = 6.67259d-8
     mp     = 1.67262d-24
     sigmaT = 6.65246d-25
     Msun   = 1.9891d33
     sol    = 2.99792d10
     Ledd   = 4. * !dpi * G * mp * sol / sigmaT * Msun
     q = {w: earr, E: FE}

     q.e = q.e / qpint1d('myfunc',1d-5,1d6) * Ledd

;     q.e = q.e / int_tabulated(q.w,q.e,/double) ; * Ledd

     if n_elements(plot) ne 0 then $
       plot,q.w,q.e*q.w,/xlog,/ylog,$
            xr=[1d-6,1e7],/xst,$
            xtit='Energy [keV]',ytit='E * F!DE!N'

     return, q
 end 

;     k = 1.380658d-16 
; k = 1.3807d-16 erg/K = 1.3807d-23 J/K 
; 1eV = 1.1604d4 K

; want kT to be in keV, where T is Kelvin.  so k needs to be in
; keV/K. so need to convert erg to keV

; 1.6021779e-09 erg/keV

function irwavs,temp
     common E, energy
     k = 1.380658d-16           ; erg/K
     k = k / 1.6022d-9          ; [erg/K] / [erg/keV] = [keV/K]
     return,energy^5.d0 * temp^(-7.d0) / (exp(energy/(k*temp))-1.d0)
 end 
