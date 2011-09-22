;+
; NAME:
;   K_LAMBDA()
;
; PURPOSE:
;   Return a variety of extinction curves.
;
; INPUTS:
;   wave  - wavelength vector [Angstroms]
;
; OPTIONAL INPUTS:
;   r_v   - total to selective extinction ratio (default 3.1, but
;           poorly known for the LMC and SMC, and inappropriate
;           for the Calzetti attenuation curve); also hard-wired
;           to be 3.1 for Li & Draine (2001) since it was only
;           calibrated for the diffuse MW ISM
;
; KEYWORD PARAMETERS:
;   calzetti - Calzetti (2000) continuum attenuation curve for IUE
;              starburst galaxies
;   charlot  - Charlot & Fall (2000) attenuation curve
;   ccm      - Cardelli, Clayton, & Mathis (1989) [Milky Way]
;   odonnell - CCM with O'Donnell (1994) coefficients [Milky Way]
;   fm       - Fitzpatrick (1999) [Milky Way]
;   avglmc   - average Large Magellanic Cloud from Gordon et
;              al. 2003, ApJ, 594, 279; R_V = 3.41+/-0.06
;   lmc2     - Large Magellanic Cloud Supershell (LMC2) field,
;              including 30 Dor from Gordon et al. 2003, ApJ, 594,
;              279; R_V = 2.76+/-0.09
;   smc      - Small Magellanic Cloud Bar region from Gordon et
;              al. 2003, ApJ, 594, 279; R_V = 2.74+/-0.13
;   li       - Li & Draine (2001) Milky Way extinction curve
;   silent   - suppress warning messages
; OUTPUTS:
;   k_lambda - defined as A(lambda)/E(B-V)
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   Much of the code for this routine has been excised from the
;   Goddard routines CCM_UNRED(), FM_UNRED(), and CALZ_UNRED() for
;   my own nefarious purposes.  
;
; PROCEDURES USED:
;   READ_SMC()
;
; EXAMPLE:
;   Plot the Calzetti and the CCM extinction curves in the
;   optical and UV.
;
;   IDL> wave = findgen(6000)+1000.0 ; Angstrom
;   IDL> plot, wave, k_lambda(wave,/calzetti), xsty=3, ysty=3
;   IDL> oplot, wave, k_lambda(wave,/ccm)
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2002 September 19, U of A
;   jm03mar2uofa - added the Charlot & Fall (2000) curve
;   jm03sep2uofa - general updates
;   jm03sep28uofa - added Seaton (1979) extinction curve
;   jm04jan26uofa - updated SMC, LMC2, and AVGLMC curves from
;                   Gordon et al. (2003)
;   jm04apr26uofa - the Charlot & Fall R_V=5.9, not 3.1 
;   jm06feb24uofa - added SILENT keyword
;   jm07feb19nyu  - added Li & Draine (2001) extinction curve 
;   jm09aug19ucsd - ensure that the Calzetti, O'Donnell and
;     SMC extinction curves behave properly at very long and short 
;     wavelengths  
;   jm11aug05ucsd - made READ_GORDON_2003() an internal support
;     routine 
;
; Copyright (C) 2002-2004, 2006-2007, 2009, 2011, John Moustakas
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

function read_gordon_2003, wave, smc=smc, lmc2=lmc2, avglmc=avglmc, r_v=r_v
; internal support routine for K_LAMBDA()
; jm04jan26uofa - written
; jm09aug20ucsd - force k(lambda) to be positive at long wavelengths 

    path = filepath('',root_dir=getenv('IMPRO_DIR'),subdirectory='etc')
    file = '2003_gordon.dat'

    if (file_test(path+file) eq 0L) then begin
       splog, 'Gordon et al. (2003) data file not found.'
       return, -1
    endif

    readcol, path+file, l, x, Al_AV_smc, err, Al_AV_lmc2, err, $
      Al_Av_avglmc, err, comment='#', /silent

    R_V_smc    = 2.74
    R_V_lmc2   = 2.76
    R_V_avglmc = 3.41

    if keyword_set(smc) then begin
       if n_elements(r_v) eq 0 then r_v = R_V_smc
       good = where(Al_AV_smc ne -9.999,ngood)
       Al_AV = Al_AV_smc[good]
       l = l[good]
    endif
    if keyword_set(lmc2) then begin
       if n_elements(r_v) eq 0 then r_v = R_V_lmc2
       good = where(Al_AV_lmc2 ne -9.999,ngood)
       Al_AV = Al_AV_lmc2[good]
       l = l[good]
    endif
    if keyword_set(avglmc) then begin
       if n_elements(r_v) eq 0 then r_v = R_V_avglmc
       good = where(Al_AV_avglmc ne -9.999,ngood)
       Al_AV = Al_AV_avglmc[good]
       l = l[good]
    endif
    
    lambda = 1D4*l
    k_lambda = R_V*Al_AV

    if (n_elements(wave) ne 0L) then k_lambda = $
      interpol(k_lambda,lambda,wave)>0.0 else wave = lambda
    
return, k_lambda
end    

function k_lambda, wave, r_v=r_v, calzetti=calzetti, charlot=charlot, $
  ccm=ccm, odonnell=odonnell, seaton=seaton, fm=fm, avglmc=avglmc, $
  lmc2=lmc2, smc=smc, li=li, silent=silent

    if (n_elements(wave) eq 0L) then begin
       doc_library, 'k_lambda'
       return, -1
    endif

    k_lambda = -1.0

    zero = where(wave le 0.0,nzero)
    if (nzero ne 0L) then begin
       splog, 'WAVE contains zero or negative values!'
       return, k_lambda
    endif
    
; select the default extinction curve (CCM+O'Donnell coefficients)
    if $
      (n_elements(calzetti) eq 0) and $
      (n_elements(charlot) eq 0)  and $
      (n_elements(ccm) eq 0)      and $
      (n_elements(odonnell) eq 0) and $
      (n_elements(seaton) eq 0)   and $
      (n_elements(fm) eq 0)       and $
      (n_elements(avglmc) eq 0)   and $
      (n_elements(lmc2) eq 0)     and $
      (n_elements(li) eq 0)       and $
      (n_elements(smc) eq 0)      then odonnell = 1 ; ccm = 1L

    if n_elements(odonnell) ne 0 then ccm = 1
;   if (n_elements(avglmc) ne 0L) or (n_elements(lmc2) ne 0L) then fm = 1L

; ---------------------------------------------------------------------------    
; Calzetti et al. (2000)
    if keyword_set(calzetti) then begin 
       if n_elements(r_v) eq 0L then r_v = 4.05
       w1 = where((wave ge 6300) and (wave le 22000),c1)
       w2 = where((wave ge  912) and (wave lt  6300),c2)
       x  = 10000.0/wave        ; wavelength in inverse microns

       k_lambda = 0.0*wave

       if c1 gt 0 then $
         k_lambda[w1] = 2.659*(-1.857 + 1.040*x[w1])+r_v
       if c2 gt 0 then $
         k_lambda[w2] = 2.659*(poly(x[w2], [-2.156, 1.509d0, -0.198d0, 0.011d0])) + r_v

; if necessary, extrapolate to longer and shorter wavelengths
       inrange = where((wave ge 912.0) and (wave le 22000.0),ninrange)
       lowave = where(wave lt 912.0,nlowave)
       hiwave = where(wave gt 22000.0,nhiwave)
       if (nlowave ne 0) or (nhiwave ne 0) then $
         if (keyword_set(silent) eq 0) then splog, 'Extrapolating '+$
         'the Calzetti k(lambda)'
       if (nlowave ne 0) then k_lambda[lowave] = k_lambda[inrange[0]]
       if (nhiwave ne 0) then k_lambda[hiwave] = interpol(k_lambda[inrange],$
         wave[inrange],wave[hiwave])>0.0
    endif

; ---------------------------------------------------------------------------    
; Charlot & Fall (2000)
    if keyword_set(charlot) then begin 
       if n_elements(r_v) eq 0 then r_v = 5.9 ; S. Charlot, private communication
       k_lambda = r_v*(wave/5500.0)^(-0.7)
    endif

; ---------------------------------------------------------------------------    
; Cardelli, Clayton, & Mathis (1989)
    if keyword_set(ccm) then begin 
       if n_elements(r_v) eq 0 then r_v = 3.1

       x = 10000./ wave         ; convert to inverse microns 
       npts = n_elements(x)
       a = fltarr(npts)  
       b = fltarr(npts)

       good = where( (x gt 0.3) and (x lt 1.1),ngood ) ; infrared
       if ngood gt 0 then begin
          a[good] =  0.574 * x[good]^(1.61)
          b[good] = -0.527 * x[good]^(1.61)
       endif

       good = where( (x ge 1.1) and (x lt 3.3),ngood) ; optical/NIR
       if ngood GT 0 then begin

          y = x[good] - 1.82

          if keyword_set(odonnell) then begin ; new coefficients from O'Donnell (1994)

             c1 = [ 1. , 0.104,   -0.609,    0.701,  1.137,    $
               -1.718,   -0.827,    1.647, -0.505 ]
             c2 = [ 0.,  1.952,    2.908,   -3.989, -7.985,    $
               11.102,    5.491,  -10.805,  3.347 ]

          endif else begin ; original coefficients

             c1 = [ 1. , 0.17699, -0.50447, -0.02427,  0.72085,$
               0.01979, -0.77530,  0.32999 ]
             c2 = [ 0.,  1.41338,  2.28305,  1.07233, -5.38434,$
               -0.62251,  5.30260, -2.09002 ]

          endelse

          a[good] = poly(y,c1)
          b[good] = poly(y,c2)

       endif

       good = where( (x ge 3.3) and (x lt 8),ngood) ; mid-UV
       if ngood gt 0 then begin

          y = x[good]
          F_a = fltarr(Ngood)    & F_b = fltarr(Ngood)
          good1 = where( (y GT 5.9), Ngood1 )
          if Ngood1 GT 0 then begin
             y1 = y[good1] - 5.9
             F_a[ good1] = -0.04473 * y1^2 - 0.009779 * y1^3
             F_b[ good1] =   0.2130 * y1^2  +  0.1207 * y1^3
          endif
          
          a[good] =  1.752 - 0.316*y - (0.104 / ( (y-4.67)^2 + 0.341 )) + F_a
          b[good] = -3.090 + 1.825*y + (1.206 / ( (y-4.62)^2 + 0.263 )) + F_b

       endif

       good = where( (x GE 8) and (x LE 11),ngood ) ; far-UV
       if ngood gt 0 then begin

          y = x[good] - 8.
          c1 = [ -1.073, -0.628,  0.137, -0.070 ]
          c2 = [ 13.670,  4.257, -0.420,  0.374 ]
          a[good] = poly(y, c1)
          b[good] = poly(y, c2)

       endif

       k_lambda = r_v * (a + b/r_v)

; if necessary, make K_LAMBDA a smooth function by extrapolating
       inrange = where((wave ge 912.0) and (wave le 10000.0/0.3),ninrange)
       lowave = where(wave lt 912.0,nlowave)
       hiwave = where(wave gt 10000.0/0.3,nhiwave)
       if (nlowave ne 0) or (nhiwave ne 0) then $
         if (keyword_set(silent) eq 0) then splog, "Extrapolating "+$
         "the O'Donnell k(lambda)"
       if (nlowave ne 0) then k_lambda[lowave] = k_lambda[inrange[0]]
       if (nhiwave ne 0) then k_lambda[hiwave] = interpol(k_lambda[inrange],$
         wave[inrange],wave[hiwave])>0.0

    endif

; ---------------------------------------------------------------------------    
; Seaton (1979) - old Milky Way
    if keyword_set(seaton) then begin

       if n_elements(r_v) eq 0L then r_v = 3.1

       x = 10000.0/wave         ; convert to inverse microns 
       k_lambda = x*0.0

; Table 2 (adopted from Nandy et al. (1975)       
       
       xprime = findgen(18)*0.1+1.0
       kprime = [1.36,1.44,1.84,2.04,2.24,2.44,2.66,2.88,3.14,$
         3.36,3.56,3.77,3.96,4.15,4.26,4.40,4.52,4.64]
              
       w1 = where((x lt 2.70),nw1)
       w2 = where((x ge 2.70) and (x le 3.65),nw2)
       w3 = where((x gt 3.65) and (x le 7.14),nw3)
       w4 = where((x gt 7.14) and (x le 10.0),nw4)

       if nw1 ne 0L then k_lambda[w1] = interpol(kprime,xprime,x[w1])
       if nw2 ne 0L then k_lambda[w2] = 1.56+1.048*x[w2]+1.01/((x[w2]-4.60)^2+0.280)
       if nw3 ne 0L then k_lambda[w3] = 2.29+0.848*x[w3]+1.01/((x[w3]-4.60)^2+0.280)
       if nw4 ne 0L then k_lambda[w4] = 16.17-3.20*x[w4]+0.2975*x[w4]^2

       k_lambda = r_v*k_lambda/3.2
       
    endif

; ---------------------------------------------------------------------------    
; Fitzpatrick (1999) - new Milky Way, AVGLMC, and LMC2
    if keyword_set(fm) then begin

       if n_elements(r_v) eq 0 then r_v = 3.1

       x = 10000./ wave         ; convert to inverse microns 
       k_lambda = x*0.

       if keyword_set(lmc2) then  begin
          if N_elements(x0) EQ 0 then x0    =  4.626
          if N_elements(gamma) EQ 0 then gamma =  1.05	
          if N_elements(c4) EQ 0 then c4   =  0.42   
          if N_elements(c3) EQ 0 then c3    =  1.92	
          if N_elements(c2) EQ 0 then c2    = 1.31
          if N_elements(c1) EQ 0 then c1    =  -2.16
       endif else if keyword_set(avglmc) then begin
          if N_elements(x0) EQ 0 then x0 = 4.596  
          if N_elements(gamma) EQ 0 then gamma = 0.91
          if N_elements(c4) EQ 0 then c4   =  0.64  
          if N_elements(c3) EQ 0 then c3    =  2.73	
          if N_elements(c2) EQ 0 then c2    = 1.11
          if N_elements(c1) EQ 0 then c1    =  -1.28
       endif else begin
          if N_elements(x0) EQ 0 then x0    =  4.596  
          if N_elements(gamma) EQ 0 then gamma =  0.99	
          if N_elements(c3) EQ 0 then c3    =  3.23	
          if N_elements(c4) EQ 0 then c4   =  0.41    
          if N_elements(c2) EQ 0 then c2    = -0.824 + 4.717/R_V
          if N_elements(c1) EQ 0 then c1    =  2.030 - 3.007*c2
       endelse

; Compute UV portion of A(lambda)/E(B-V) curve using FM fitting function and 
; R-dependent coefficients
       
       xcutuv = 10000.0/2700.0
       xspluv = 10000.0/[2700.0,2600.0]
       iuv = where(x ge xcutuv, N_UV)
       IF (N_UV GT 0) THEN xuv = [xspluv,x[iuv]] ELSE  xuv = xspluv

       yuv = c1  + c2*xuv
       yuv = yuv + c3*xuv^2/((xuv^2-x0^2)^2 +(xuv*gamma)^2)
       yuv = yuv + c4*(0.5392*((xuv>5.9)-5.9)^2+0.05644*((xuv>5.9)-5.9)^3)
       yuv = yuv + R_V
       yspluv  = yuv[0:1]       ; save spline points

       IF (N_UV GT 0) THEN k_lambda[iuv] = yuv[2:*] ; remove spline points
       
; Compute optical portion of A(lambda)/E(B-V) curve
; using cubic spline anchored in UV, optical, and IR

       xsplopir = [0,10000.0/[26500.0,12200.0,6000.0,5470.0,4670.0,4110.0]]
       ysplir   = [0.0,0.26469,0.82925]*R_V/3.1 
       ysplop   = [poly(R_V, [-4.22809e-01, 1.00270, 2.13572e-04] ), $
         poly(R_V, [-5.13540e-02, 1.00216, -7.35778e-05] ), $
         poly(R_V, [ 7.00127e-01, 1.00184, -3.32598e-05] ), $
         poly(R_V, [ 1.19456, 1.01707, -5.46959e-03, 7.97809e-04, $ 
         -4.45636e-05] ) ]
       
       ysplopir = [ysplir,ysplop]

       iopir = where(x lt xcutuv, Nopir)
       if (Nopir GT 0) then $
         k_lambda[iopir] = CSPLINE([xsplopir,xspluv],[ysplopir,yspluv],x[iopir])

    endif
       
; ---------------------------------------------------------------------------    
; Small Magellanic Cloud
    if keyword_set(smc) then k_lambda = read_gordon_2003(wave,/smc,r_v=r_v)

; ---------------------------------------------------------------------------    
; Large Magellanic Cloud
    if keyword_set(lmc2) then k_lambda = read_gordon_2003(wave,/lmc2,r_v=r_v)
    
; ---------------------------------------------------------------------------    
; Large Magellanic Cloud Average
    if keyword_set(avglmc) then k_lambda = read_gordon_2003(wave,/avglmc,r_v=r_v)

; ---------------------------------------------------------------------------    
; Li & Draine (2001)
    if keyword_set(li) then begin 
       path = filepath('',root_dir=getenv('IMPRO_DIR'),subdirectory='etc')
       data = rsex(path+'li_draine01.dat')
       if n_elements(r_v) eq 0 then r_v = 3.1
       k_lambda = r_v*2.146D21*interpol(data.sigma_ext,data.wave*1D4,wave,/spline)
    endif
    
    if n_elements(k_lambda) eq 1 then k_lambda = k_lambda[0]
    
return, k_lambda
end
