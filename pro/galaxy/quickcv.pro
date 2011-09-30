;+
; NAME:
;
;  QUICKCV
;
; PURPOSE:
;
;   calculate cosmic variance for some geometry
;
; CATEGORY:
;
;   cosmic variance
;
;
; CALLING SEQUENCE:
;
;   fracerror=quickcv(side1,side2 [zarr,deltaz,sig8=sig8,npoints=npoints]) 
;
; INPUTS:
;
;   side1: length of 1 side of field on sky, in degrees (rect. geom)
;   side2: length of other side of field on sky, in degrees
;
; OPTIONAL INPUTS:
;   zarr: array of z's for which cosmic variance will be calculated; default:1
;   deltaz: total length of volume in z direction; default: 0.1
;
; KEYWORD PARAMETERS:
;   sig8: desired value of sigma_8; default: 0.77
;   npoints: desired number of integration points for first part of 
;     integration; default: 200
;
;
; OUTPUTS:
;
;  fracerror: array containing the fractional error in a count,
;  sigma/N, due to cosmic variance for objects with bias b=1 
;  in a volume side1 x side2 degrees x delta_z centered at redshifts
;  zarr.  NOTE THIS IS SIGMA, NOT VARIANCE!!!!
;
; RESTRICTIONS:
;
;  uses power spectrum in pofk.pro (gamma=0.1872, omega0=0.26), distance
;  relation in rz.pro 
;
;
; MODIFICATION HISTORY:
; modified by rss july 2006 to include Dlin scaling with z
; modified by bpm sep  2009 to wmap3 cosmology
;-

function pofk,karr,gamma=gamma, n=n

if n_elements(n) eq 0 then n = 0.95
if n_elements(gamma) eq 0 then gamma=0.1872

        omega0=0.26
        sigma8 = 0.58*omega0^(-0.47+0.16*omega0)
        h=1.0

; determine sigma8 from model
; following Peacock, p. 528

        keff=(0.172+0.011*(alog(gamma/0.34))^2)*1.
        q=keff/(h*gamma)
        tk=alog(1+2.34*q)/(2.34*q)*(1+3.89*q+(16.1*q)^2+(5.46*q)^3+(6.71*q)^4)^(-0.25)
        sigma8squn=keff^(3.0+n)*tk^2

        q=karr/(h*gamma)   ; k/h Gamma
        tk=alog(1+2.34*q)/(2.34*q)*(1+3.89*q+(16.1*q)^2+(5.46*q)^3+(6.71*q)^4)^(-0.25)
        delsq=karr^n*tk^2       ; took out factor of k^3
        delsq=delsq*(sigma8)^2/sigma8squn        ;normalize to Borgani et al. sigma 8, omegam=0.26,omegal=0.74

        pofk=2*!Pi^2*delsq

return,pofk
end



function intpk4,nsteps=nsteps,rangemul=rangemul,nmc=nmc,x1=x1,x2=x2,x3=x3

        if n_elements(rangemul) eq 0 then rangemul=1.0
        

        if n_elements(nsteps) eq 0 then nsteps=200
        if n_elements(nmc) eq 0 then nmc=floor(nsteps/3.)


        nsteps=long(nsteps)
        
        if n_elements(x1) eq 0 then x1=(20.20)/2.
        if n_elements(x2) eq 0 then x2=(80.78)/2.
        if n_elements(x3) eq 0 then x3=(1301.4)/2.

        ;print,'in intpk4: ', x1,x2,x3

        k1max=(0.4+16/x1)*1.25*rangemul
        k2max=(0.4+16/x2)*1.25*rangemul
        k3max=(0.4+16/x3)*1.25*rangemul

; sufficient to encompass entire region where sin^2=1/2 not adequate:
        k1maxb=4.*!pi/x1
        k2maxb=4.*!pi/x2
        k3maxb=4.*!pi/x3

        nstep2=nmc
;floor(nsteps/5)

        dk1=(k1max-k1maxb)/(nstep2-1)
        dk2=(k2max-k2maxb)/(nstep2-1)
        dk3=(k3max-k3maxb)/(nstep2-1)

        dk1b=k1maxb/(nsteps-1)
        dk2b=k2maxb/(nsteps-1)
        dk3b=k3maxb/(nsteps-1)

        k1arr=findgen(nsteps+1)
        k2arr=k1arr*dk2b
        k3arr=k1arr*dk3b
        k1arr=k1arr*dk1b

        w1 = k1arr*0+1
        gt0 = where(k1arr gt 0) ; jm11sep23ucsd
        w1[gt0]=sin(k1arr[gt0]*x1)/k1arr[gt0]/x1

        w2 = k2arr*0+1
        gt0 = where(k2arr gt 0) ; jm11sep23ucsd
        w2[gt0]=sin(k2arr[gt0]*x2)/k2arr[gt0]/x2

        w3 = k3arr*0+1
        gt0 = where(k3arr gt 0) ; jm11sep23ucsd
        w3[gt0]=sin(k3arr[gt0]*x3)/k3arr[gt0]/x3
        
;       w1=sin(k1arr*x1)/k1arr/x1
;       w1(where(finite(w1) eq 0))=1.
;
;       w2=sin(k2arr*x2)/k2arr/x2
;       w2(where(finite(w2) eq 0))=1.
;
;       w3=sin(k3arr*x3)/k3arr/x3
;       w3(where(finite(w3) eq 0))=1.

        w1=w1^2
        w2=w2^2
        w3=w3^2

        k1b=exp(alog(k1max/k1maxb)/nstep2*(findgen(nstep2)+1.))*k1maxb
        k1arr=[k1arr,k1b]
        w1=[w1,(0.5)/x1^2/k1b^2]
        k2b=exp(alog(k2max/k2maxb)/nstep2*(findgen(nstep2)+1.))*k2maxb
        k2arr=[k2arr,k2b]
        w2=[w2,(0.5)/x2^2/k2b^2]
        k3b=exp(alog(k3max/k3maxb)/nstep2*(findgen(nstep2)+1.))*k3maxb
        k3arr=[k3arr,k3b]
        w3=[w1,(0.5)/x3^2/k3b^2]

        integral=0.

        nsteps=nsteps+nstep2
;plot,k3arr

        for i=0,nsteps  do begin
                k1=k1arr(i)
                windowf=w1(i)
                dk1=(k1arr((i+1)<nsteps)-k1arr((i-1)>0))/2.
;               if i eq 0 or i eq nsteps then dk1=dk1*2.
                ktot=k1^2

                ind2=((lindgen(nsteps+1,nsteps+1)) mod (nsteps+1))
                dk2=(k2arr((ind2+1)<nsteps)-k2arr((ind2-1)>0))/2.
;               wh=where(ind2 eq 0 or ind2 eq nsteps)
;               dk2(wh)=dk2(wh)*2.
                ktot=ktot+k2arr(ind2)^2
                windowf=windowf*w2(ind2)
                ind2=1.
                
                ind3=(lindgen(nsteps+1,nsteps+1) / (nsteps+1))
                dk3=(k3arr((ind3+1)<nsteps)-k3arr((ind3-1)>0))/2.
;               wh=where(ind3 eq 0 or ind3 eq nsteps)
;               dk3(wh)=dk3(wh)*2.
                ktot=ktot+k3arr(ind3)^2
                windowf=windowf*w3(ind3)
                ind3=1.
                ktot=sqrt(ktot)

                pk=pofk(ktot)
                if i eq 0 then pk(0,0)=0.
                integral=integral+total(double(pk*windowf*dk1*dk2*dk3))
        endfor



        return,integral*1/!pi^3
end


function Dlin, Omega_0, z
;assumes a flat cosmology
  Omega = Omega_0*(1+z)^3/(1-Omega_0 + (1+z)^3*Omega_0)
  gz = 2.5*Omega/(1./70. + 209.0*Omega/140. + Omega^(4./7.))
  g0 = 2.5*Omega_0/(1./70. + 209.0*Omega_0/140. + Omega_0^(4./7.))
  D = gz/(g0*(1+z))
  return, D
end


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


; Make Fig. 13.5 from Peebles (1993) book: 
; angular size / physical size vs. z

PRO alcock_test1, $
                  zMin = zMin, $
                  zMax = zMax, $
                  nz = nz
   
   IF NOT keyword_set(zMin) THEN zMin = 0.001d0
   IF NOT keyword_set(zMax) THEN zMax = 10.d0
   IF NOT keyword_set(nz) THEN nz = 200
   
   dz = (zMax-zMin) / (nz-1.d0)
   z = zMin + dz*findgen(nz)
   
   lincolr
   
; models with no curvature
   
   window, 0
   r = rz(z) & plot, z, (1+z)/r, yr=[0,10], /ys, color=1
   r = rz(z, OmegaM=2, OmegaL=-1) & oplot, z, (1+z)/r, color=2
   r = rz(z, OmegaM=0.5, OmegaL=0.5) & oplot, z, (1+z)/r, color=3
   r = rz(z, OmegaM=0.2, OmegaL=0.8) & oplot, z, (1+z)/r, color=4
   r = rz(z, OmegaM=0.1, OmegaL=0.9) & oplot, z, (1+z)/r, color=5
   r = rz(z, OmegaM=0.05, OmegaL=0.95) & oplot, z, (1+z)/r, color=6
   
; models with no lambda
   
   window, 1
   r = rz(z) & plot, z, (1+z)/r, yr=[0,10], /ys, color=1
   r = rz(z, OmegaM=2) & oplot, z, (1+z)/r, color=2
   r = rz(z, OmegaM=0.5) & oplot, z, (1+z)/r, color=3
   r = rz(z, OmegaM=0.2) & oplot, z, (1+z)/r, color=4
   r = rz(z, OmegaM=0.1) & oplot, z, (1+z)/r, color=5
   r = rz(z, OmegaM=0.05) & oplot, z, (1+z)/r, color=6
   
END


; Plot comoving size / angular size vs. z

PRO alcock_plot1, $
                  zMin = zMin, $
                  zMax = zMax, $
                  nz = nz, $
                  _extra = extra
   
   IF NOT keyword_set(zMin) THEN zMin = 0.d0
   IF NOT keyword_set(zMax) THEN zMax = 4.d0
   IF NOT keyword_set(nz) THEN nz = 200
   
   dz = (zMax-zMin) / (nz-1.d0)
   z = zMin + dz*findgen(nz)
   
   lincolr
   !p.multi = [0,1,2]
   
   r = rz(z) 
   mplot, z, r, _extra=extra, $
     mxtitle="z", mytitle=textoidl("H_0 r(z) / c"), $
     mtitle="Comoving size / Angular size", $
     charsize=1.5
   
   r = rz(z, OmegaM=0.2, OmegaL=0.8) & oplot, z, r, linesty=4, color=3
   r = rz(z, OmegaM=0.4, OmegaL=0.6) & oplot, z, r, linesty=5, color=3
   r = rz(z, OmegaM=0.2) & oplot, z, r, linesty=2, color=2
   r = rz(z, OmegaM=0.4) & oplot, z, r, linesty=3, color=2
   oplot, [0,10], [0,10], linesty=1
   
   xr = [0.55, 0.65]
   yr = [0.4, 0.32, 0.24, 0.16]
   reg2dev, xr, yr[0]*[1,1], x, y
   plots, x, y, linesty=4, color=3, /device
   reg2dev, xr, yr[1]*[1,1], x, y
   plots, x, y, linesty=5, color=3, /device
   reg2dev, xr, yr[2]*[1,1], x, y
   plots, x, y, linesty=2, color=2, /device
   reg2dev, xr, yr[3]*[1,1], x, y
   plots, x, y, linesty=3, color=2, /device

   xr = [0.67]
   yr = yr - 0.02
   reg2dev, xr, [yr[0]], x, y
   xyouts, x, y, textoidl("\Omega_m=0.2, \Omega_\Lambda=0.8"), /device, $
     color=3, charsize=1.4
   reg2dev, xr, [yr[1]], x, y
   xyouts, x, y, textoidl("\Omega_m=0.4, \Omega_\Lambda=0.6"), /device, $
     color=3, charsize=1.4
   reg2dev, xr, [yr[2]], x, y
   xyouts, x, y, textoidl("\Omega_m=0.2, \Omega_\Lambda=0"), /device, $
     color=2, charsize=1.4
   reg2dev, xr, [yr[3]], x, y
   xyouts, x, y, textoidl("\Omega_m=0.4, \Omega_\Lambda=0"), /device, $
     color=2, charsize=1.4
   xyouts, 2.3, 0.95, "E-dS", orientation=15, chars=1.2
   xyouts, 0.8, 0.9, "Euclidean", orientation=55, chars=1.2

   r = rz(z)
   mplot, z, r, _extra=extra, $
     charsize=1.5

   r = rz(z, OmegaM=0.4, OmegaQ=0.6, wQ=-1.) 
   oplot, z, r, color=3
   r = rz(z, OmegaM=0.4, OmegaQ=0.6, wQ=-0.6667) 
   oplot, z, r, color=2, lines=2
   r = rz(z, OmegaM=0.4, OmegaQ=0.6, wQ=-0.5) 
   oplot, z, r, color=4, lines=3
   r = rz(z, OmegaM=0.4, OmegaQ=0.6, wQ=-0.3333) 
   oplot, z, r, color=7, lines=4
   r = rz(z, OmegaM=0.4, OmegaQ=0.6, wQ=-0.1667) 
   oplot, z, r, color=6, lines=5
   oplot, [0,10], [0,10], linesty=1

   xr = [0.55, 0.65]
   yr = [0.4, 0.32, 0.24, 0.16]
   reg2dev, xr, yr[0]*[1,1], x, y
   plots, x, y, linesty=2, color=2, /device
   reg2dev, xr, yr[1]*[1,1], x, y
   plots, x, y, linesty=3, color=4, /device
   reg2dev, xr, yr[2]*[1,1], x, y
   plots, x, y, linesty=4, color=7, /device
   reg2dev, xr, yr[3]*[1,1], x, y
   plots, x, y, linesty=5, color=6, /device
   
   xr = [0.67]
   yr = yr - 0.02
   reg2dev, xr, [yr[0]], x, y
   xyouts, x, y, textoidl("w=-2/3"), /device, color=2, charsize=1.4
   reg2dev, xr, [yr[1]], x, y
   xyouts, x, y, textoidl("w=-1/2"), /device, color=4, charsize=1.4
   reg2dev, xr, [yr[2]], x, y
   xyouts, x, y, textoidl("w=-1/3"), /device, color=7, charsize=1.4
   reg2dev, xr, [yr[3]], x, y
   xyouts, x, y, textoidl("w=-1/6"), /device, color=6, charsize=1.4
   xyouts, 2.3, 0.8, "E-dS (w=0)", orientation=15, chars=1.2
   xyouts, 0.8, 0.9, "Euclidean", orientation=55, chars=1.2
   xyouts, 0.3, 1.3, textoidl("\Omega_m=0.4"), chars=1.5
   xyouts, 0.3, 1.2, textoidl("\Omega_Q=0.6"), chars=1.5
   xyouts, 1.7, 1.08, textoidl("\Lambda (w=-1)"), orient=25, chars=1.2, color=3
   
END


FUNCTION squish_factor, z, $
                        _extra=extra
   
   f = z / (rz(z, _extra=extra) * hubble(z))

   RETURN, f
END





function quickcv,side1,side2,za,deltaz,sig8=sig8,npoints=npoints

if n_elements(npoints) eq 0 then npoints=200.

; desired sigma8 
if n_elements(sig8) eq 0 then sig8=0.77

if n_elements(deltaz) eq 0 then deltaz=0.10


if n_elements(za) eq 0 then za=1.
nz=n_elements(za)

nfields=1


; assume LCDM
; uses jonathan baker's old code

;rza=rz(za,omegam=0.3,omegal=0.7) ; jm11sep23ucsd
;rzmin=rz(za-deltaz/2,omegam=0.3,omegal=0.7)
;rzmax=rz(za+deltaz/2,omegam=0.3,omegal=0.7)
rza=rz(za,omegam=0.26,omegal=0.74)
rzmin=rz(za-deltaz/2,omegam=0.26,omegal=0.74)
rzmax=rz(za+deltaz/2,omegam=0.26,omegal=0.74)
cvi=fltarr(nz)


x1a=3000./2.*(rza * side1)/!radeg
x2a=3000./2.*(rza * side2)/!radeg
x3a=3000.*(rzmax-rzmin)/2.

for i=0,nz-1 do  $
  cvi(i)=intpk4(x1=x1a(i),x2=x2a(i),x3=x3a(i),nsteps=npoints)

;cvi is the fractional variance in a count in a redshift bin
;cvar is the fractional error (sigma)

cvi=cvi*sig8^2/1.033^2
  d = Dlin(0.26, za)
  ;print, 'dlin: ', d
  cvar=sqrt(cvi)*d

;print,'Fractional sigma: ',cvar

  return,cvar

end
