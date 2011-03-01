;+
; NAME:
;       IABUNDANCE()
;
; PURPOSE:
;       Calculate strong-line abundances using abundance diagnostics
;       in Kewley & Dopita (2002, ApJS, 142, 35).
;
; CALLING SEQUENCE:
;
;
; INPUTS:
;       line - reddening-corrected emission-line flux structure
;
; OPTIONAL INPUTS:
;
;
; KEYWORD PARAMETERS:
;
;
; OUTPUTS:
;
;
; OPTIONAL OUTPUTS:
;
;
; COMMENTS:
;       Use IUNRED_LINEFIT() to correct the data for dust extinction. 
;
; PROCEDURES USED:
;
;
; EXAMPLE:
;
;
; MODIFICATION HISTORY:
;       Written by Lisa Kewley (Harvard), unknown date
;       J. Moustakas, 2003 June 13, U of A - adopted original program
;       for use with ISPEC.  major re-write of I/O
;
;-

function parse_line, line

; set line fluxes that were not measured or that are upper limits to
; zero

    tags = tag_names(line)
    linename = strcompress(line[0].linename,/remove)
    nline = n_elements(linename)
    nspec = n_elements(line)

    for j = 0L, nline-1L do if (j eq 0L) then $
      out = create_struct(linename[j],0.0) else $
      out = create_struct(out,linename[j],0.0)
    out = replicate(out,nspec)

    for i = 0L, nline-1L do begin

       match = where(linename[i] eq tags)
       array = line.(match)
       flux = reform(array[0,*])
       ferr = reform(array[1,*])
       
       neg = where(ferr lt 0.0,nneg)
       if nneg ne 0L then flux[neg] = 0.0
       out.(i) = flux

    endfor
    
return, out
end

function iabundance, line

    num = n_elements(line)
    if num eq 0L then begin
       print, 'Syntax - abundance = iabundance(line)
       return, -1L
    endif

; ---------------------------------------------------------------------------    
; L. Kewley's original comments
; ---------------------------------------------------------------------------    
    
; ## This program does extinction correction, abundance determination, 
; ## and outputs results in "filename_abund.txt", "filename_abund2.txt", and
; ## "filename__corrfluxes.txt" where filename is given in the filename
; ## parameter below.  Extinction correction uses the Whitford reddening
; ## curve as parameterized by Miller & Mathews (1972).

; ## Outputs:
; ## --------
; ## "filename_abund.txt" contains a list of the abundance estimates for each ID
; ## using the diagnostics given in Kewley & Dopita (2002).  These include our
; ## recommended combined diagnostic, [NII]/[OII], R23, R23 combined, 
; ## [NII]/[SII], McGaugh (1991; M91), Zaritsky (1994; Z94), 
; ## Charlot (2001; C01), and the average of M91, Z94 and C01 (referred to as 
; ## the comparison average or comp ave).

; ## "filename_abund2.txt" is a shorter version of "filename_abund.txt" and only
; ## contains the abundance estimated from out combined diagnostic, R23 
; ## combined diagnostic, and the comparison average.

; ## "filename_corrfluxes.txt" contains a list of the extinction corrected
; ## emission-line flux ratios used in the Baldwin, Phillips & Terlevich
; ## AGN-starburst diagnostic diagrams.  To determine which galaxies contain
; ## AGN, we recommend the use of the theoretical diagnostic line derived 
; ## using the same models as for our abundance diagnostics.  This theoretical 
; ## line is defined in Kewley et al. (2001a, ApJ, 556, 121) and used in 
; ## Kewley et al. (2001b ApJS, 132, 37).  An IDL script which reads 
; ## filename_corrfluxes.txt, makes the classification according to 
; ## Kewley et al. 2001b and plots the data on the diagnostic diagrams is
; ## available from lkewley@cfa.harvard.edu.  Abundance estimates for galaxies 
; ## containing AGN will not be reliable.

; ## Inputs:
; ## -------
; ## Accepts any non-extinction corrected or extinction corrected data set 
; ## input in 12 columns (ID number, fluxes) corresponding to:
; ## ID number, [OII]3727, Hb, [OIII]4959, [OIII]5007, [OI]6300, Ha, [NII]6584, 
; ## [SII]6717, [SII] 6731,[SIII]9069, [SIII]9532
; ## If any of these line fluxes are unavailable, please put 0.0 in that
; ## column.  ID number should be an integer.

; ---------------------------------------------------------------------------    

; parse the input data structure

    parsed = parse_line(line)    

; rename the data structure fields to correspond to Lisa's convention 

    name           = lindgen(num)

    if tag_exist(parsed,'OII_3727')  then OII3727_raw    = parsed.OII_3727  else OII3727_raw = name*0.0
    if tag_exist(parsed,'H_BETA')    then Hb_raw         = parsed.H_BETA    else Hb_raw = name*0.0
    if tag_exist(parsed,'OIII_4959') then OIII4959_raw   = parsed.OIII_4959 else OIII4959_raw = name*0.0
    if tag_exist(parsed,'OIII_5007') then OIII5007_raw   = parsed.OIII_5007 else OIII5007_raw = name*0.0
    if tag_exist(parsed,'OI_6300')   then OI6300_raw     = parsed.OI_6300   else OI6300_raw = name*0.0
    if tag_exist(parsed,'H_ALPHA')   then Ha_raw         = parsed.H_ALPHA   else Ha_raw = name*0.0
    if tag_exist(parsed,'NII_6584')  then NII6584_raw    = parsed.NII_6584  else NII6584_raw = name*0.0
    if tag_exist(parsed,'SII_6716')  then SII6717_raw    = parsed.SII_6716  else SII6717_raw = name*0.0
    if tag_exist(parsed,'SII_6731')  then SII6731_raw    = parsed.SII_6731  else SII6731_raw = name*0.0
    if tag_exist(parsed,'SIII_9069') then SIII9069_raw   = parsed.SIII_9069 else SIII9069_raw = name*0.0
    if tag_exist(parsed,'SIII_9532') then SIII9532_raw   = parsed.SIII_9532 else SIII9532_raw = name*0.0

    OIII49595007_raw=fltarr(num)
    SII67176731_raw=fltarr(num)
    SIII90699532_raw=fltarr(num)

; initialize the output data structure

    abundance = {$
      KD_id:                  0L, $
      KD_Z_combined:         0.0, $
      KD_NII_OII:            0.0, $
      KD_R23_combined:       0.0, $
      KD_R23:                0.0, $
      KD_R23_flag:            0L, $ ; corresponds to R23_CHECK.TXT
      KD_NII_SII:            0.0, $
      KD_S23:                0.0, $
      KD_C01:                0.0, $
      KD_M91:                0.0, $
      KD_Z94:                0.0, $
      KD_Comparison_average: 0.0, $
      KD_q:                  0.0}
    abundance = replicate(abundance,num)

    red_corr = 'no'
    
;;; ## enter the 3 parameters below:

;;   num=10                       ; number of entries in data
;;    filename='Fluxes.txt'       ; input data filename 
;;    red_corr='yes'              ; if data need to be reddening corrected, this 
                                  ;  should be 'yes', otherwise 'no'

;;; ## the following should not need changing: ##
;;
;;    data=fltarr(12,num)
;;    openr,1,filename
;;    readf,1,data
;;    close,1
;;
;;    name=intarr(num)
;;    name=floor(data(0,*))
;;    OII3727_raw=data(1,*)
;;    Hb_raw=data(2,*)
;;    OIII4959_raw=data(3,*)
;;    OIII5007_raw=data(4,*)
;;    OI6300_raw=data(5,*)
;;    Ha_raw=data(6,*)
;;    NII6584_raw=data(7,*)
;;    SII6717_raw=data(8,*)
;;    SII6731_raw=data(9,*)
;;    SIII9069_raw=data(10,*)
;;    SIII9532_raw=data(11,*)
;;    OIII49595007_raw=fltarr(num)
;;    SII67176731_raw=fltarr(num)
;;    SIII90699532_raw=fltarr(num)

    for i=0,num-1 do begin
       if (OII3727_raw(i) lt 0.0) then begin
          OII3727_raw(i)=0.0
       endif
       if (OIII4959_raw(i) lt 0.0) then begin
          OIII4959_raw(i)=0.0
       endif
       if (Hb_raw(i) lt 0.0) then begin
          Hb_raw(i)=0.0
       endif
       if (OIII5007_raw(i) lt 0.0) then begin
          OIII5007_raw(i)=0.0
       endif
       if (NII6584_raw(i) lt 0.0) then begin
          NII6584_raw(i)=0.0
       endif
       if (Ha_raw(i) lt 0.0) then begin
          Ha_raw(i)=0.0
       endif

       if (OIII4959_raw(i) eq 0.0) and (OIII5007_raw(i) ne 0.0) then begin
          OIII4959_raw(i)=0.347*OIII5007_raw(i)
       endif
       if (OIII4959_raw(i) ne 0.0) and (OIII5007_raw(i) ne 0.0) then begin
          OIII49595007_raw(i)=OIII4959_raw(i)+OIII5007_raw(i)
       endif 
       if (SII6717_raw(i) ne 0.0) and (SII6731_raw(i) ne 0.0) then begin
          SII67176731_raw(i)=SII6717_raw(i)+SII6731_raw(i)
       endif
       if (SIII9069_raw(i) ne 0.0) and (SIII9532_raw(i) ne 0.0) then begin
          SIII90699532_raw(i)=SIII9069_raw(i)+SIII9532_raw(i)
       endif
    end

; extinction correction:

; correction constants from Miller & Mathews 1972

    fOII=1.29945
    fHb=1.0
    fOIII=0.961067              ; for OIII5007, but is only 0.01 different from OIII4959 (0.973617)
    fSII=0.63068
    fNII=0.650591
    fSIII=0.380485
    fOI=0.695
    fHa=0.655 

; OIII=4959+5007, else OIII5007 is specified
    EB_V=fltarr(num)
    c_EB_V=fltarr(num)
    logNIIOII=fltarr(num)
    logNIIHa=fltarr(num)
    logSIIHa=fltarr(num)
    logOIHa=fltarr(num)
    logOIIIHb=fltarr(num)
    logNIISII=fltarr(num)
    logR23=fltarr(num)
    logOIIIOII=fltarr(num)
    logS23=fltarr(num)
    logSIIISII=fltarr(num)
    logOIII5007OII=fltarr(num)
    logOIIOIII5007=fltarr(num)
    logOIIOIII5007Hb=fltarr(num)
    R23=fltarr(num)
    OIII5007OII=fltarr(num)
    OIIOIII5007=fltarr(num)
    OIIIOII=fltarr(num)
    logOIIIOII=fltarr(num)
    R23_5007=fltarr(num)
    OIIOIII5007Hb=fltarr(num)
    logOIIOIII5007Hb=fltarr(num)
    
    if red_corr eq 'yes' then begin 
       for i=0,num-1 do begin
          if (Hb_raw(i) ne 0.0) and (Ha_raw(i) ne 0.0) then begin
             EB_V(i)=2.1575*alog10(Ha_raw(i)/(Hb_raw(i)*2.85)) ; E(B-V)
             c_EB_V(i)=EB_V(i)/0.77 ; constant, c
             c=c_EB_V(i)
             if (NII6584_raw(i) ne 0.0) and (OII3727_raw(i) ne 0.0) then begin
                logNIIOII(i)=alog10(NII6584_raw(i)/OII3727_raw(i))+EB_V(i)/0.77*(fNII-fOII) 
             endif 
             if (OIII49595007_raw(i) ne 0.0) and (OII3727_raw(i) ne 0.0) then begin    
                logR23(i)=alog10((OII3727_raw(i)/Hb_raw(i))*10^(c*(fOII-fHb)) +  $
                  (OIII49595007_raw(i)/Hb_raw(i))*10^(c*(fOIII-fHb)))
                R23(i)=((OII3727_raw(i)/Hb_raw(i))*10^(c*(fOII-fHb)) +  $
                  (OIII49595007_raw(i)/Hb_raw(i))*10^(c*(fOIII-fHb)))
                logOIIIOII(i)=alog10(OIII49595007_raw(i)/OII3727_raw(i))+EB_V(i)/0.77*(fOIII-fOII)        
             endif    
             if (SII67176731_raw(i) ne 0.0) and (SIII90699532_raw(i) ne 0.0) then begin    
                logS23(i)=alog10((SII67176731_raw(i)/Hb_raw(i))*10^(c*(fSII-fHb)) + $
                  (SIII90699532_raw(i)/Hb_raw(i))*10^(c*(fSIII-fHb)))
                logSIIISII(i)=alog10(SIII90699532_raw(i)/SII67176731_raw(i))+EB_V(i)/0.77*(fSIII-fSII)
             endif    
             if (OIII5007_raw(i) ne 0.0) and (OII3727_raw(i) ne 0.0) then begin
                logOIII5007OII(i)=alog10(OIII5007_raw(i)/OII3727_raw(i))+EB_V(i)/0.77*(fOIII-fOII)
                logOIIOIII5007(i)=alog10(OII3727_raw(i)/OIII5007_raw(i))+EB_V(i)/0.77*(fOII-fOIII)
                logOIIOIII5007Hb(i)=alog10((OII3727_raw(i) + OIII5007_raw(i))/Hb_raw(i))
                                ; ratios for other diagnostics - slightly different ratios needed
                OIII5007OII(i)=10^(logOIII5007OII(i))
                OIIOIII5007(i)=10^(logOIIOIII5007(i))
                OIIIOII(i)=1.347*OIII5007OII(i)
                logOIIIOII(i)=alog10(OIIIOII(i))
                                ; R23 but without [OIII]4959
                R23_5007(i)=(1./OIII5007OII(i) + 1.)/(1./OIII5007OII(i) + 1.347)*R23(i)  
                OIIOIII5007Hb(i)=R23_5007(i)
                logOIIOIII5007Hb(i)=alog10(OIIOIII5007Hb(i))
             endif
             if (NII6584_raw(i) ne 0.0) and (SII67176731_raw(i) ne 0.0) then begin  
                logNIISII(i)=alog10(NII6584_raw(i)/SII67176731_raw(i))+EB_V(i)/0.77*(fNII-fSII) 
             endif
             if (NII6584_raw(i) ne 0.0) and (Ha_raw(i) ne 0.0) then begin  
                logNIIHa(i)=alog10(NII6584_raw(i)/Ha_raw(i))+EB_V(i)/0.77*(fNII-fHa) 
             endif
             if (SII67176731_raw(i) ne 0.0) and (Ha_raw(i) ne 0.0) then begin  
                logSIIHa(i)=alog10(SII67176731_raw(i)/Ha_raw(i))+EB_V(i)/0.77*(fSII-fHa) 
             endif
             if (OI6300_raw(i) ne 0.0) and (Ha_raw(i) ne 0.0) then begin  
                logOIHa(i)=alog10(OI6300_raw(i)/Ha_raw(i))+EB_V(i)/0.77*(fOI-fHa) 
             endif
             if (OIII5007_raw(i) ne 0.0) and (Hb_raw(i) ne 0.0) then begin
                logOIIIHb(i)=alog10(OIII5007_raw(i)/Hb_raw(i))+EB_V(i)/0.77*(fOIII-fHb)
             endif
          endif
       end
    endif else begin
       for i=0,num-1 do begin 
          if (NII6584_raw(i) ne 0.0) and (OII3727_raw(i) ne 0.0) then begin
             logNIIOII(i)=alog10(NII6584_raw(i)/OII3727_raw(i))
          endif
          if (OIII49595007_raw(i) ne 0.0) and (OII3727_raw(i) ne 0.0) then begin    
             logR23(i)=alog10((OII3727_raw(i) + OIII49595007_raw(i))/Hb_raw(i))
             R23(i)=((OII3727_raw(i) + OIII49595007_raw(i))/Hb_raw(i))
             logOIIIOII(i)=alog10(OIII49595007_raw(i)/OII3727_raw(i))     
          endif    
          if (SII67176731_raw(i) ne 0.0) and (SIII90699532_raw(i) ne 0.0) then begin    
             logS23(i)=alog10((SII67176731_raw(i) + SIII90699532_raw(i))/Hb_raw(i))
             logSIIISII(i)=alog10(SIII90699532_raw(i)/SII67176731_raw(i))
          endif    
          if (OIII5007_raw(i) ne 0.0) and (OII3727_raw(i) ne 0.0) then begin
             logOIII5007OII(i)=alog10(OIII5007_raw(i)/OII3727_raw(i))
             logOIIOIII5007(i)=alog10(OII3727_raw(i)/OIII5007_raw(i))
             logOIIOIII5007Hb(i)=alog10((OII3727_raw(i) + OIII5007_raw(i))/Hb_raw(i))
             OIII5007OII(i)=10^(logOIII5007OII(i))
             OIIOIII5007(i)=10^(logOIIOIII5007(i))
                                ; ratios for other diagnostics - slightly different ratios needed
             OIIIOII(i)=1.347*OIII5007OII(i)
             logOIIIOII(i)=alog10(OIIIOII(i))
                                ; R23 but without [OIII]4959
             R23_5007(i)=(1./OIII5007OII(i) + 1.)/(1./OIII5007OII(i) + 1.347)*R23(i) 
             OIIOIII5007Hb(i)=R23_5007(i)
             logOIIOIII5007Hb(i)=alog10(OIIOIII5007Hb(i))
          endif
          if (NII6584_raw(i) ne 0.0) and (SII67176731_raw(i) ne 0.0) then begin  
             logNIISII(i)=alog10(NII6584_raw(i)/SII67176731_raw(i))
          endif
          if (NII6584_raw(i) ne 0.0) and (Ha_raw(i) ne 0.0) then begin  
             logNIIHa(i)=alog10(NII6584_raw(i)/Ha_raw(i))
          endif
          if (SII67176731_raw(i) ne 0.0) and (Ha_raw(i) ne 0.0) then begin  
             logSIIHa(i)=alog10(SII67176731_raw(i)/Ha_raw(i))
          endif
          if (OI6300_raw(i) ne 0.0) and (Ha_raw(i) ne 0.0) then begin  
             logOIHa(i)=alog10(OI6300_raw(i)/Ha_raw(i))
          endif
          if (OIII5007_raw(i) ne 0.0) and (Hb_raw(i) ne 0.0) then begin
             logOIIIHb(i)=alog10(OIII5007_raw(i)/Hb_raw(i))
          endif
       end
    endelse

; ## calculating z from Kobulnicky parametrization of Zaritzky et al. (1994) ##

    Z94_Z=fltarr(num)

    x=logR23
    for i=0,num-1 do begin
       if x(i) ne 0.0 then begin
          Z94_Z(i)=9.265-0.33*x(i)-0.202*x(i)^2-0.207*x(i)^3-0.333*x(i)^4 
       endif
    end

; ## calculating Charlot 01 calibration: (case A) ##

    NIISII=10^logNIISII
    x2=OIIOIII5007/1.5
    x6=NIISII/0.85
    C01_Z=fltarr(num)

    for i=0,num-1 do begin 
       if (x2(i) ne 0.0) and (x6(i) ne 0.0) then begin
          C01_Z(i)=alog10(5.09e-4*(x2(i)^0.17)*(x6(i)^1.17))+12
       endif
    end

; ## calculating M91 calibration using Z94 as initial estimate of abundance ##:
; this initial estimate can be changed by replacing OH_init by another guess, 
; eg C01_Z 

    Z_init=Z94_Z

;Z_init=KD02_NIIOII_Z

    x=logR23
    y=logOIIIOII
    M91_Z=fltarr(num)

    for i=0,num-1 do begin
       if (x(i) ne 0.0) and (y(i) ne 0.0) then begin
          if (Z_init(i) le 8.4) then begin
             M91_Z(i)=12.0-4.944+0.767*x(i)+0.602*x(i)^2-y(i)*(0.29+0.332*x(i)-0.331*x(i)^2)
          endif else begin
             M91_Z(i)=12.0-2.939-0.2*x(i)-0.237*x(i)^2-0.305*x(i)^3-0.0283*x(i)^4-y(i)*$
               (0.0047-0.0221*x(i)-0.102*x(i)^2-0.0817*x(i)^3-0.00717*x(i)^4)
          endelse
       endif
    end

; ## calculating Kewley & Dopita (2002) (KD02) estimates of abundance ##

; KD02 [NII]/[OII] estimate (can be used for whole log(O/H)+12 range, but rms 
; scatter increases to 0.11 rms for log(O/H)+12 < 8.6;  rms = 0.04 for
; log(O/H)+12 > 8.6
; uses equation (4) from paper:

    NIIOII_roots=complexarr(4,num)
    KD02_NIIOII_Z=fltarr(num)

    for i=0,num-1 do begin
       if logNIIOII(i) ne 0.0 then begin
          NIIOII_coef=[1106.87,-532.154,96.3733,-7.81061,0.239282] ; q=2e7 line (approx average)
          NIIOII_coef(0)=NIIOII_coef(0)-logNIIOII(i)
          NIIOII_roots(*,i)=fz_roots(NIIOII_coef,/double) ; finding roots for eq (4)
          for k=0,3 do begin    ; root must be real and
             if imaginary(NIIOII_roots(k,i)) eq 0.0 then begin ; between 7.5 and 9.4
                if (abs(NIIOII_roots(k,i)) ge 7.5) and (abs(NIIOII_roots(k,i)) le 9.4) then begin
                   KD02_NIIOII_Z(i)=abs(NIIOII_roots(k,i))
                endif
             endif
          end
       endif
    end  

; ### KD02 [NII]/[OII] estimate ###
; (can be used for for log(O/H)+12 > 8.6 only)
; uses equation (5) from paper, this should be identical 
; to the estimate above for abundances log(O/H)+12 > 8.6

    KD02_NIIOII2_Z=fltarr(num)

    for i=0,num-1 do begin
       if logNIIOII(i) ne 0.0 then begin
          KD02_NIIOII2_Z(i)=alog10(8.511e-4*(1.54020+1.26602*logNIIOII(i)+$
            0.167977*logNIIOII(i)^2))+12.
       endif
    end

; if [NII]/[OII] after extinction correction is less than -1.5, then check the data.
; if it is only slightly less than 1.5, then this can be a result of either noisy
; data, inaccurate fluxes or extinction correction, or a higher ionization parameter
; than modelled.  For these cases, the average of the M91,Z94 and C01 should be used.


; HERE!!!!!!!!!!!!!!!!!!!


; KD02 R23 estimate (not reliable for  8.4 < log(O/H)+12 < 8.8)
; uses [NII]/[OII] estimate as initial guess - this can be changed below

    KD02_R23_Z=fltarr(num)

; initializing:

    Zi=[0.05,0.1,0.2,0.5,1.0,1.5,2.0,3.0] ; model grid abundances in solar units
    ZiOH=alog10(Zi*8.511e-4)+12 ; log(O/H)+12 units
    Zstep=[0.025,0.075,0.15,0.35,0.75,1.25,1.75,2.5,3.5] ; middle of model grid abundances
    ZstepOH=alog10(Zstep*8.511e-4)+12
    qstep=alog10([3.5e6,7.5e6,1.5e7,3e7,6e7,1.16e8,2.25e8,3.25e8]) ; model grid ionization parameters
    n_iter=3                    ; number of iterations for abundance determination
    tol=1.0e-2                  ; tolerance for convergance 
    R23_roots=complexarr(4,num) ; all possible roots of R23 diagnostic
    q_roots=complexarr(3,num)   ; possible roots of q diagnostic
    q_R23=fltarr(num,n_iter+1) 
    q=fltarr(num,n_iter+1)      ; actual q value
    OIIIOII_coef=fltarr(4,8)    ; coefficients from model grid fits
    R23_coef=fltarr(5,7)        ; coefficients from model grid fits
    R23_Z=fltarr(num,n_iter+1)  ; Z value for each iteration

    R23_Z(*,0)=KD02_NIIOII_Z    ; use [NII]/[OII] abundance as initial estimate
                                ; may be changed

; occasionally, for noisy data or badly fluxed [OII],[OIII] or Hb, 
; or for high ionization parameter galaxies, R23 is slightly higher
; than the curves in our models - this will result in all complex roots of
; the R23 curve unless a slightly lower R23 is used.  These should
; be checked individually to make sure that it is just noise etc in the
; data that is causing the problem, rather than wrong fluxes input.
; the R23 ratio should be close to 0.95, not much more than 1.0 for the
; data to be physical.

; ---------------------------------------------------------------------------    
; jm03jun13uofa - added
; ---------------------------------------------------------------------------    
    check = where(logR23 gt 0.95,ncheck)
    if ncheck ne 0L then abundance[check].KD_R23_flag = 1L
; ---------------------------------------------------------------------------    

;    openw,1,'R23_check.txt'
;    for i=0,num-1 do begin
;       if logR23(i) gt 0.95 then begin
;          logR23(i)=0.95
;          printf,1,i,name(i),logR23(i)
;       endif
;    end
;    close,1

    for i=0,num-1 do begin
       if (logOIII5007OII(i) ne 0.0) and (logR23(i) ne 0.0) then begin
                                ; coefficients from KD02 paper:

          R23_coef(*,0)=[-3267.93,1611.04,-298.187,24.5508,-0.758310] ; q=5e6
          R23_coef(*,1)=[-3727.42,1827.45,-336.340,27.5367,-0.845876] ; q=1e7
          R23_coef(*,2)=[-4282.30,2090.55,-383.039,31.2159,-0.954473] ; q=2e7
          R23_coef(*,3)=[-4745.18,2309.42,-421.778,34.2598,-1.04411] ; q=4e7
          R23_coef(*,4)=[-4516.46,2199.09,-401.868,32.6686,-0.996645] ; q=8e7
          R23_coef(*,5)=[-3509.63,1718.64,-316.057,25.8717,-0.795242] ; q=1.5e8
          R23_coef(*,6)=[-1550.53,784.262,-149.245,12.6618,-0.403774] ; q=3e8

          OIIIOII_coef(*,0)=[-36.9772,10.2838 ,-0.957421,0.0328614] ;z=0.05 
          OIIIOII_coef(*,1)=[-74.2814,24.6206,-2.79194,0.110773] ; z=0.1
          OIIIOII_coef(*,2)=[-36.7948,10.0581,-0.914212,0.0300472] ; z=0.2
          OIIIOII_coef(*,3)=[-81.1880,27.5082,-3.19126,0.128252] ; z=0.5
          OIIIOII_coef(*,4)=[-52.6367,16.0880,-1.67443,0.0608004] ; z=1.0
          OIIIOII_coef(*,5)=[-86.8674,28.0455,-3.01747,0.108311] ; z=1.5
          OIIIOII_coef(*,6)=[-24.4044,2.51913,0.452486,-0.0491711] ; z=2.0
          OIIIOII_coef(*,7)=[49.4728,-27.4711,4.50304,-0.232228] ; z=3.0

          OIIIOII_coef(0,*)=OIIIOII_coef(0,*)-logOIII5007OII(i)
          R23_coef(0,*)=R23_coef(0,*)-logR23(i)

          for iter=1,n_iter do begin ; iterate if tolerance level not met
             if abs(R23_Z(i,iter)-R23_Z(i,iter-1)) ge tol then begin

                                ;   calculate ionization parameter using [OIII]/[OII] with
                                ;   [NII]/[OII] abundance for the first iteration, and the R23
                                ;   abundance in consecutive iterations

                for j=0,7 do begin   
                   if R23_Z(i,iter-1) gt ZstepOH(j) then begin
                      if R23_Z(i,iter-1) le ZstepOH(j+1) then begin
                         q_roots(*,i)=fz_roots(OIIIOII_coef(*,j),/double)
                         q_roots(*,i)=fz_roots(OIIIOII_coef(*,j),/double)        
                      endif
                   endif
                end

                                ;   q must be between 3.5e6 and 3.5e8 cm/s because that is the range it
                                ;   is defined over by the model grids, and it must be real.

                for k=0,2 do begin
                   if imaginary(q_roots(k,i)) eq 0.0 then begin
                      if float(q_roots(k,i)) ge 6.54407 then begin ;log units (q ge 1e6 cm/s) 
                         if float(q_roots(k,i)) le 8.30103 then begin ;log units (q le 2e8 cm/s)
                            q(i,iter)=float(q_roots(k,i))
                            q_R23(i,iter)=float(q_roots(k,i))
                         endif
                      endif
                   endif
                end
                                ;   calculate abundance using ionization parameter:

                for j=0,6 do begin   
                   if q(i,iter) gt qstep(j) then begin
                      if q(i,iter) le qstep(j+1) then begin
                         R23_roots(*,i)=fz_roots(R23_coef(*,j),/double)
                         R23_qstepno=j
                      endif
                   endif
                end

                
                                ;   There will be four roots, two complex ones, and two real ones.
                                ;   use previous R23 value (or [NII]/[OII] if first iteration) 
                                ;   and q to find which real root to use (ie. which side of R23 curve 
                                ;   to use).  

                nn=0
;   Rmax=[1.04967,1.06497,1.06684,1.06329,1.03844,0.991261,0.91655]
                Smax=[8.69020,8.65819,8.61317,8.58916,8.49012,8.44109,8.35907]

                for k=0,3 do begin
                   if imaginary(R23_roots(k,i)) eq 0.0 then begin
                      if (R23_Z(i,iter-1) ge Smax(R23_qstepno)) and (float(R23_roots(k,i)) ge Smax(R23_qstepno)) then begin
                         R23_Z(i,iter)=float(R23_roots(k,i))
                      endif 
                      if (R23_Z(i,iter-1) le Smax(R23_qstepno)) and (float(R23_roots(k,i)) le Smax(R23_qstepno)) then begin
                         R23_Z(i,iter)=float(R23_roots(k,i))
                      endif
                   endif
                end

; around maximum of R23 sometimes the R23 ratio will be slightly higher than
; that available for the closest q value.  This will depend on noise addded to
; data.  If this happens, step up in ionization parameter to find abundance
; using that one instead. Around local maximum, the actual ionization parameter
; used is not significant compared to the errors associated with the lack of
; abundance sensitivity of the R23 ratio in this region.

                while (R23_Z(i,iter) eq 0.0) and (R23_qstepno le 5) do begin
                   R23_roots(*,i)=fz_roots(R23_coef(*,R23_qstepno+1),/double)
                   for k=0,3 do begin
                      if imaginary(R23_roots(k,i)) eq 0.0 then begin
                         if (R23_Z(i,iter-1) ge Smax(R23_qstepno)) and (float(R23_roots(k,i)) ge Smax(R23_qstepno)) then begin
                            R23_Z(i,iter)=float(R23_roots(k,i))
                         endif 
                         if (R23_Z(i,iter-1) le Smax(R23_qstepno)) and (float(R23_roots(k,i)) le Smax(R23_qstepno)) then begin
                            R23_Z(i,iter)=float(R23_roots(k,i))
                            
                         endif
                      endif
                   end
                   R23_qstepno=R23_qstepno+1
                endwhile

             endif else begin
                R23_Z(i,iter)=R23_Z(i,iter-1)  
                q(i,iter)=q(i,iter-1)
             endelse
          end
       endif
    end

    KD02_R23_Z=R23_Z(*,n_iter)
    q_OIIIOII=q(*,n_iter)

;  ### Combined R23 method outlined in KD02 paper Section 7. ###
;  ie for objects with only [OII], [OIII], Hb available

    KD02_R23comb_Z=fltarr(num)

    for i=0,num-1 do begin
       if Z94_Z(i) ge 9.0 then begin 
          KD02_R23comb_Z(i)=(KD02_R23_Z(i)+M91_Z(i)+Z94_Z(i))/3 ; my technique averaged with M91 and Z94
       endif 
       if (KD02_R23comb_Z(i) le 9.0) and (KD02_R23comb_Z(i) ge 8.5) then begin
          KD02_R23comb_Z(i)=(M91_Z(i)+Z94_Z(i))/2 ; average of M91 and Z94
       endif
       if (Z94_Z(i) le 9.0) and (Z94_Z(i) ge 8.5) then begin
          KD02_R23comb_Z(i)=(M91_Z(i)+Z94_Z(i))/2 ; average of M91 and Z94
       endif 
       if KD02_R23comb_Z(i) le 8.5 then begin
          KD02_R23comb_Z(i)=KD02_R23_Z(i) ; my technique
       endif
       if Z94_Z(i) le 8.5 then begin
          KD02_R23comb_Z(i)=KD02_R23_Z(i) ; my technique
       endif
    end

; HERE!!!!!!!!!!!!!!!!


; ### [NII]/[SII] method outlined in KD02 paper ###
; this method produces a systematic shift of 0.2 dex in log(O/H)+12
; compared with the average of M91, Z94, and C01.  We believe this
; is a result of inaccurate abundances or depletion factors, which are known 
; problems in sulfur modelling.  Initial guess of [NII]/[OII] used
; can be changed.  Initial guess is not critical except for high
; ionization parameters.  ionization parameter diagnostic is [OIII]/[OII]

    KD02_NIISII_Z=fltarr(num)


    n_iter=3                    ; number of iterations for abundance determination
    tol=1.0e-2                  ; tolerance for convergance 
    NIISII_roots=complexarr(4,num) ; all possible roots of NIISII diagnostic
    q_roots=complexarr(3,num)   ; possible roots of q diagnostic
    q=fltarr(num,n_iter+1)      ; actual q value
    q_NIISII=fltarr(num,n_iter+1) ; actual q value
    OIIIOII_coef=fltarr(4,8)    ; coefficients from model grid fits
    NIISII_coef=fltarr(5,7)     ; coefficients from model grid fits
    NIISII_Z=fltarr(num,n_iter+1) ; Z value for each iteration

; initializing:

    NIISII_Z(*,0)=KD02_NIIOII_Z ; use [NII]/[OII] abundance as initial estimate
                                ; may be changed

    for i=0,num-1 do begin
       if (logOIII5007OII(i) ne 0.0) and (logNIISII(i) ne 0.0) then begin
                                ; coefficients from KD02 paper:

          NIISII_coef(*,0)=[-1042.47,521.076,-97.1578,8.00058,-0.245356]
          NIISII_coef(*,1)=[-1879.46,918.362,-167.764,13.5700,-0.409872]
          NIISII_coef(*,2)=[-2027.82,988.218,-180.097,14.5377,-0.438345]
          NIISII_coef(*,3)=[-2080.31,1012.26,-184.215,14.8502,-0.447182]
          NIISII_coef(*,4)=[-2162.93,1048.97,-190.260,15.2859,-0.458717]
          NIISII_coef(*,5)=[-2368.56,1141.97,-205.908,16.4451,-0.490553]
          NIISII_coef(*,6)=[-2910.63,1392.18,-249.012,19.7280,-0.583763]

          OIIIOII_coef(*,0)=[-36.9772,10.2838 ,-0.957421,0.0328614] ;z=0.05 
          OIIIOII_coef(*,1)=[-74.2814,24.6206,-2.79194,0.110773] ; z=0.1
          OIIIOII_coef(*,2)=[-36.7948,10.0581,-0.914212,0.0300472] ; z=0.2
          OIIIOII_coef(*,3)=[-81.1880,27.5082,-3.19126,0.128252] ; z=0.5
          OIIIOII_coef(*,4)=[-52.6367,16.0880,-1.67443,0.0608004] ; z=1.0
          OIIIOII_coef(*,5)=[-86.8674,28.0455,-3.01747,0.108311] ; z=1.5
          OIIIOII_coef(*,6)=[-24.4044,2.51913,0.452486,-0.0491711] ; z=2.0
          OIIIOII_coef(*,7)=[49.4728,-27.4711,4.50304,-0.232228] ; z=3.0

          OIIIOII_coef(0,*)=OIIIOII_coef(0,*)-logOIII5007OII(i)
          NIISII_coef(0,*)=NIISII_coef(0,*)-logNIISII(i)

          for iter=1,n_iter do begin ; iterate if tolerance level not met
             if abs(NIISII_Z(i,iter)-NIISII_Z(i,iter-1)) ge tol then begin

                                ;   calculate ionization parameter using [OIII]/[OII] with
                                ;   [NII]/[OII] abundance for the first iteration, and the [NII]/[SII]
                                ;   abundance in consecutive iterations

                for j=0,7 do begin   
                   if NIISII_Z(i,iter-1) gt ZstepOH(j) then begin
                      if NIISII_Z(i,iter-1) le ZstepOH(j+1) then begin
                         q_roots(*,i)=fz_roots(OIIIOII_coef(*,j),/double)
                      endif
                   endif
                end

                                ;   q must be between 3.5e6 and 3.5e8 cm/s because that is the range it
                                ;   is defined over by the model grids, and it must be real.

                for k=0,2 do begin
                   if imaginary(q_roots(k,i)) eq 0.0 then begin
                      if float(q_roots(k,i)) ge 6.54407 then begin ;log units (q ge 1e6 cm/s) 
                         if float(q_roots(k,i)) le 8.30103 then begin ;log units (q le 2e8 cm/s)
                            q(i,iter)=float(q_roots(k,i))
                            q_NIISII(i,iter)=float(q_roots(k,i))
                         endif
                      endif
                   endif
                end

                                ;   calculate abundance using ionization parameter:

                for j=0,6 do begin   
                   if q(i,iter) gt qstep(j) then begin
                      if q(i,iter) le qstep(j+1) then begin
                         NIISII_roots(*,i)=fz_roots(NIISII_coef(*,j),/double)
                         q_NIISII_used=j
                      endif
                   endif
                end
                
                                ;   There will be four roots, two complex ones, and two real ones.
                                ;   use previous NIISII value (or [NII]/[OII] if first iteration) 
                                ;   and q to find which real root to use (ie. which side of R23 curve 
                                ;   to use).  

                nn=0

                for k=0,3 do begin
                   if imaginary(NIISII_roots(k,i)) eq 0.0 then begin
                      if (float(NIISII_roots(k,i)) ge 8.0) and (float(NIISII_roots(k,i)) le 9.35) then begin
                         NIISII_Z(i,iter)=float(NIISII_roots(k,i))
                      endif 
                   endif
                end

                if NIISII_Z(i,iter) eq 0.0 then begin
                   NIISII_roots(*,i)=fz_roots(NIISII_coef(*,q_NIISII_used+1),/double)
                   for k=0,3 do begin
                      if imaginary(NIISII_roots(k,i)) eq 0.0 then begin
                         if (float(NIISII_roots(k,i)) ge 8.0) and (float(NIISII_roots(k,i)) le 9.35) then begin
                            NIISII_Z(i,iter)=float(NIISII_roots(k,i))
                         endif 
                      endif
                   end
                endif  


             endif else begin
                NIISII_Z(i,iter)=NIISII_Z(i,iter-1)  
                q(i,iter)=q(i,iter-1)
             endelse
          end
       endif
    end

    KD02_NIISII_Z=NIISII_Z(*,n_iter)

; ### KD02 S23 estimate (not reliable at all)  ###
; We believe the inaccuracies in this method are a result of inaccurate 
; abundances or depletion factors, and systematically overestimating the
; [SIII]/[SII] ratio, all of which are known problems in sulfur modelling. 
; uses [NII]/[OII] estimate as initial guess - this can be changed below
; used [SIII]/[SII] as ionization parameter diagnostic

    KD02_S23_Z=fltarr(num)

    n_iter=3                    ; number of iterations for abundance determination
    tol=1.0e-2                  ; tolerance for convergance 
    S23_roots=complexarr(4,num) ; all possible roots of S23 diagnostic
    q_roots=complexarr(3,num)   ; possible roots of q diagnostic
    q=fltarr(num,n_iter+1)      ; actual q value
    q_S23=fltarr(num,n_iter+1)  ; actual q value
    SIIISII_coef=fltarr(4,8)    ; coefficients from model grid fits
    S23_coef=fltarr(5,7)        ; coefficients from model grid fits
    S23_Z=fltarr(num,n_iter+1)  ; Z value for each iteration

; initializing:

    S23_Z(*,0)=KD02_NIIOII_Z    ; use [NII]/[OII] abundance as initial estimate
                                ; may be changed

    for i=0,num-1 do begin
       if (logSIIISII(i) ne 0.0) and (logS23(i) ne 0.0) then begin
                                ; coefficients from KD02 paper:

          S23_coef(*,0)=[-1543.68, 761.018, -141.061, 11.6389, -0.360280] ; q=5e6
          S23_coef(*,1)=[-1542.15, 758.664, -140.351, 11.5597, -0.357261] ; q=1e7
          S23_coef(*,2)=[-1749.48, 855.280, -157.198, 12.8626, -0.394978] ; q=2e7
          S23_coef(*,3)=[-1880.06, 914.362, -167.192, 13.6119, -0.415997] ; q=4e7
          S23_coef(*,4)=[-1627.10, 790.891, -144.699, 11.7988, -0.361423] ; q=8e7
          S23_coef(*,5)=[-1011.65, 497.017, -92.2429, 7.64915, -0.238660] ; q=1.5e8
          S23_coef(*,6)=[-81.6519, 55.3453, -13.7783, 1.46716, -0.0563760] ; q=3e8

          SIIISII_coef(*,0)=[30.0116, -14.8970, 2.30577, -0.112314] ;z=0.05 
          SIIISII_coef(*,1)=[16.8569,-9.62876,1.59938,-0.0804552] ; z=0.1
          SIIISII_coef(*,2)=[32.2358,-15.2438,2.27251,-0.106913] ; z=0.2
          SIIISII_coef(*,3)=[-3.06247,-0.864092,0.328467,-0.0196089] ; z=0.5
          SIIISII_coef(*,4)=[-2.94394,-0.546041,0.239226,-0.0136716 ] ; z=1.0
          SIIISII_coef(*,5)=[-38.1338,13.0914,-1.51014,0.0605926] ; z=1.5
          SIIISII_coef(*,6)=[-21.0240,6.55748,-0.683584,0.0258690] ; z=2.0
          SIIISII_coef(*,7)=[-6.61131,1.36836,-0.0717560,0.00225792] ; z=3.0

          
          SIIISII_coef(0,*)=SIIISII_coef(0,*)-logSIIISII(i)
          S23_coef(0,*)=S23_coef(0,*)-logS23(i)

          for iter=1,n_iter do begin ; iterate if tolerance level not met
             if abs(S23_Z(i,iter)-S23_Z(i,iter-1)) ge tol then begin

                                ;   calculate ionization parameter using [SIII]/[SII] with
                                ;   [NII]/[OII] abundance for the first iteration, and the S23
                                ;   abundance in consecutive iterations

                for j=0,7 do begin   
                   if S23_Z(i,iter-1) gt ZstepOH(j) then begin
                      if S23_Z(i,iter-1) le ZstepOH(j+1) then begin
                         q_roots(*,i)=fz_roots(SIIISII_coef(*,j),/double)
                      endif
                   endif
                end

                                ;   q must be between 3.5e6 and 3.5e8 cm/s because that is the range it
                                ;   is defined over by the model grids, and it must be real.

                for k=0,2 do begin
                   if imaginary(q_roots(k,i)) eq 0.0 then begin
                      if float(q_roots(k,i)) ge 6.54407 then begin ;log units (q ge 1e6 cm/s) 
                         if float(q_roots(k,i)) le 8.54407 then begin ;log units (q le 3.5e8 cm/s)
                            q(i,iter)=float(q_roots(k,i))
                            q_S23(i,iter)=float(q_roots(k,i))
                         endif
                      endif
                   endif
                end

                                ;   calculate abundance using ionization parameter:

                for j=0,6 do begin   
                   if q(i,iter) gt qstep(j) then begin
                      if q(i,iter) le qstep(j+1) then begin
                         S23_roots(*,i)=fz_roots(S23_coef(*,j),/double)
                         S23_qstepno=j
                      endif
                   endif
                end

                
                                ;   There will be four roots, two complex ones, and two real ones.
                                ;   use previous S23 value (or [NII]/[OII] if first iteration) 
                                ;   and q to find which real root to use (ie. which side of S23 curve 
                                ;   to use).  

                nn=0

                for k=0,3 do begin
                   if imaginary(S23_roots(k,i)) eq 0.0 then begin
                      if (S23_Z(i,iter-1) ge 8.9) and (float(S23_roots(k,i)) ge 8.9) then begin
                         S23_Z(i,iter)=float(S23_roots(k,i))
                      endif 
                      if (S23_Z(i,iter-1) le 8.9) and (float(S23_roots(k,i)) le 8.9) then begin
                         S23_Z(i,iter)=float(S23_roots(k,i))
                      endif
                   endif
                end
                while S23_Z(i,iter) eq 0.0 do begin
                   S23_qstepno=S23_qstepno-1
                   S23_roots(*,i)=fz_roots(S23_coef(*,S23_qstepno),/double)
                   for k=0,3 do begin
                      if imaginary(S23_roots(k,i)) eq 0.0 then begin
                         if (S23_Z(i,iter-1) ge 8.9) and (float(S23_roots(k,i)) ge 8.9) then begin
                            S23_Z(i,iter)=float(S23_roots(k,i))
                         endif 
                         if (S23_Z(i,iter-1) le 8.9) and (float(S23_roots(k,i)) le 8.9) then begin
                            S23_Z(i,iter)=float(S23_roots(k,i))
                         endif
                      endif
                   end
                endwhile


             endif else begin
                S23_Z(i,iter)=S23_Z(i,iter-1)  
                q(i,iter)=q(i,iter-1)
             endelse
          end
       endif
    end

    KD02_S23_Z=S23_Z(*,n_iter)
    q_SIIISII=q(*,n_iter)

; KD01 combined method (uses [NII], [OII], [OIII], [SII]):
; uses KD02 [NII]/[OII] method if [NII]/[OII] gives 8.6 < log(O/H)+12
; uses average of M91 and Z94 if 8.5 < log(O/H)+12 < 8.6
; uses average of C01 and KD02 if  log(O/H)+12 < 8.5
; Also calculates comparison average.  The combined method is
; the preferred abundance estimator.  See Kewley et al. 2002 for
; more details.

    KD02_comb_Z=fltarr(num)
    KD02C01_ave=fltarr(num)
    M91Z94C01_ave=fltarr(num)
    M91Z94_ave=fltarr(num)

    for i=0,num-1 do begin
       if (M91_Z(i) ne 0.0) and (Z94_Z(i) ne 0.0) then begin
          M91Z94_ave(i)=(M91_Z(i)+Z94_Z(i))/2 
          if (C01_Z(i) ne 0.0) then begin
             M91Z94C01_ave(i)=(M91_Z(i)+Z94_Z(i)+C01_Z(i))/3 ; comparison average
          endif
       endif
       if (KD02_R23_Z(i) ne 0.0) and (C01_Z(i) ne 0.0) then begin
          KD02C01_ave(i)=(KD02_R23_Z(i)+C01_Z(i))/2
       endif
       if KD02_NIIOII_Z(i) le 8.6 then begin
          if M91Z94_ave(i) ge 8.5 then begin
             KD02_comb_Z(i)=M91Z94_ave(i) ; average of M91 and Z94
          endif else begin
             KD02_comb_Z(i)=KD02C01_ave(i)
          endelse
       endif else begin
          KD02_comb_Z(i)=KD02_NIIOII_Z(i)
       endelse
    end

; ---------------------------------------------------------------------------    
; store the outputs in a structure and return
; ---------------------------------------------------------------------------    

    abundance.KD_id = name
    abundance.KD_Z_combined = KD02_comb_Z
    abundance.KD_NII_OII = KD02_NIIOII_Z
    abundance.KD_R23_combined = KD02_R23comb_Z
    abundance.KD_R23 = KD02_R23_Z
    abundance.KD_NII_SII = KD02_NIISII_Z
    abundance.KD_S23 = KD02_S23_Z
    abundance.KD_C01 = C01_Z
    abundance.KD_M91 = M91_Z
    abundance.KD_Z94 = Z94_Z
    abundance.KD_Comparison_average = M91Z94C01_ave
    abundance.KD_q = q_OIIIOII

; ---------------------------------------------------------------------------    
    
;;;; print out abundances derived:
;;;
;;;    outfilename=filename+'_abund.txt'
;;;    openw,1,outfilename
;;;    printf,1,'-----------------------------------------------------------------------'
;;;    printf,1,'### Abundances found for data set ',filename
;;;    printf,1,'KD02 = Kewley, L. J., & Dopita, M. A., 2002, ApJS, 142,35 '
;;;    printf,1,'C01 = Charlot, S., & Longhetti, M., 2001, MNRAS, 323, 887'
;;;    printf,1,'M91 = McGaugh, S.S., 1991, ApJ, 380, 140'
;;;    printf,1,'Z94 = Zaritsky, D., Kennicutt, R. C., & Huchra, J. P. 1994, ApJ, 420, 87'
;;;    printf,1,'------------------------------------------------------------------------'
;;;    printf,1,'Method:#       ID          Combined   [NII]/[OII]    R23 Comb      R23  '
;;;    printf,1,'  [NII]/[SII]          S23         C01           M91           Z94    Comp Ave'
;;;
;;;    for i=0,num-1 do begin
;;;       printf,1,i,name(i),KD02_comb_Z(i),KD02_NIIOII_Z(i),KD02_R23comb_Z(i),KD02_R23_Z(i),KD02_NIISII_Z(i),KD02_S23_Z(i),C01_Z(i),M91_Z(i),Z94_Z(i),M91Z94C01_ave(i)
;;;    end
;;;    close,1
;;;
;;;; print out abundances derived but only if KD02_comb, KD02_R23comb, or Comp Ave are not zero:
;;;
;;;    outfilename=filename+'_abund2.txt'
;;;    openw,1,outfilename
;;;    printf,1,'-----------------------------------------------------------------------'
;;;    printf,1,'### Abundances found for data set ',filename
;;;    printf,1,'KD02 = Kewley, L. J., & Dopita, M. A., 2002, ApJS, 142,35 '
;;;    printf,1,'C01 = Charlot, S., & Longhetti, M., 2001, MNRAS, 323, 887'
;;;    printf,1,'M91 = McGaugh, S.S., 1991, ApJ, 380, 140'
;;;    printf,1,'Z94 = Zaritsky, D., Kennicutt, R. C., & Huchra, J. P. 1994, ApJ, 420, 87'
;;;    printf,1,'------------------------------------------------------------------------'
;;;    printf,1,'Method:#       ID          Combined     R23 Comb    Comp Ave'
;;;
;;;    for i=0,num-1 do begin
;;;       zerotest=KD02_comb_Z(i)+KD02_R23comb_Z(i)+M91Z94C01_ave(i)
;;;
;;;       if (zerotest ne 0.0) then begin
;;;          printf,1,i,name(i),KD02_comb_Z(i),KD02_R23comb_Z(i),M91Z94C01_ave(i)
;;;       endif
;;;    end
;;;    close,1
;;;
;;;; print out extinction and extinction-corrected line ratios for classification:
;;;
;;;    if red_corr eq 'yes' then begin 
;;;       outfilename=filename+'_corrfluxes.txt'
;;;       openw,1,outfilename
;;;       printf,1,'#    ID   E(B-V)   log([OIII]/Hb)  log([OI]/Ha)   log([NII]/Ha)  log([SII]/Ha)'
;;;       for i=0,num-1 do begin
;;;          zerotest=logOIIIHb(i)+logOIHa(i)+logNIIHa(i)+logSIIHa(i)
;;;;      if zerotest ne 0.0 then begin
;;;          printf,1,i,name(i),EB_V(i),logOIIIHb(i),logOIHa(i),logNIIHa(i),logSIIHa(i)
;;;;      endif
;;;       end
;;;       close,1
;;;    endif
;;;
;;;    if red_corr eq 'yes' then begin 
;;;       outfilename=filename+'_corrfluxes.txt.short'
;;;       openw,1,outfilename
;;;       for i=0,num-1 do begin
;;;          zerotest=logOIIIHb(i)+logOIHa(i)+logNIIHa(i)+logSIIHa(i)
;;;;      if zerotest ne 0.0 then begin
;;;          printf,1,i,name(i),EB_V(i),logOIIIHb(i),logOIHa(i),logNIIHa(i),logSIIHa(i)
;;;;      endif
;;;       end
;;;       close,1
;;;    endif

return, abundance
end
