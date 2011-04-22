function quickcv,side1,side2,za,deltaz,sig8=sig8,npoints=npoints
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
;   sig8: desired value of sigma_8; default: 0.85.  
;
;     IMPORTANT NOTE: This should be sigma_8 at the redshift of
;     interest, not sigma_8 today.
;
;   A simple solution is to use a measured correlation function to
;   estimate the net sigma_8 of the sample (hence, the cosmic variance
;   estimate will be for that sample, not for an unbiased tracer of
;   dark matter).  
;   Following Adelberger's paper on clustering from counts-in-cells :
;    sig8=sqrt(72.*(r0/8)^gamma/(3.-gamma)/(4.-gamma)/(6.-gamma)/2.^gamma)
;   where xi(r)=(r/r0)^(-gamma) and r0 is in comoving units.
;
;   npoints: desired number of integration points for first part of 
;     integration; default: 100 (for ~1 part in 1000 accuracy; 200 is
;     better if high accuracy is required).
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
;  uses power spectrum in pofk.pro (gamma=0.21, omega0=0.3), distance
;  relation in rz.pro , integration in intpk4.pro
;
; EXAMPLE: 
;   print,quickcv(0.5,0.7,1.,0.1)
;    will print out the fractional error in a count of an unbiased
;    tracer of mass over a DEEP2 pointing (0.5x0.7 deg) and
;    delta_z=0.1, with central redshift 1
;
; MODIFICATION HISTORY:
;   released JAN 6/16/04
;-

if n_elements(npoints) eq 0 then npoints=100.

; desired sigma8 
if n_elements(sig8) eq 0 then sig8=0.85

if n_elements(deltaz) eq 0 then deltaz=0.10


if n_elements(za) eq 0 then za=1.
nz=n_elements(za)

nfields=1


; assume LCDM
; uses jonathan baker's old code

rza=rz(za,omegam=0.3,omegal=0.7)
rzmin=rz(za-deltaz/2,omegam=0.3,omegal=0.7)
rzmax=rz(za+deltaz/2,omegam=0.3,omegal=0.7)
cvi=fltarr(nz)


x1a=3000./2.*(rza * side1)/!radeg

x2a = 3000./2.*(rza * side2)/!radeg
x3a=3000.*(rzmax-rzmin)/2.

for i=0,nz-1 do $
  cvi(i)=intpk4(x1=x1a(i),x2=x2a(i),x3=x3a(i),nsteps=npoints)

;cvi is the fractional variance in a count in a redshift bin
;cvar is the fractional error (sigma)

cvi=cvi*sig8^2/.964^2
cvar=sqrt(cvi)

;print,'Fractional sigma: ',cvar

return,cvar

end

