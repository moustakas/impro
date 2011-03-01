;+
; NAME:
;   BUILD_PEGASE_BASEMODELS
;
; PURPOSE:
;   Generate a grid of basic PEGASE models that will be used by 
;   BUILD_ISEDFIT_SFHGRID to build a more complex set of SFHs.
;
; INPUTS:
;
; OPTIONAL INPUTS:
;
; KEYWORD PARAMETERS:
;   chabrier - use the CHABRIER IMF (default is the SALPETER IMF)
;
; OUTPUTS:
;   An information structure is written to
;   getenv('ISEDFIT_SFHGRID_DIR')+'/basemodels/' and the models
;   themselves are written to the 'pegase' subdirectory.
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   Rather than using the PEGASE FORTRAN code, these base models could
;   be built using my own SFH convolution code.
;
; EXAMPLES:
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2010 Mar 19, UCSD - based on BUILD_BC03_BASEMODELS
;
; Copyright (C) 2009, John Moustakas
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

pro build_pegase_basemodels, remake_ssps=remake_ssps, test=test

    splog, 'Building the PEGASE basemodels'

    pegpath = getenv('PEGASE_HR_DIR')+'/data/user_defined/'
    rootpath = getenv('ISEDFIT_SFHGRID_DIR')+'/basemodels/'
    outpath = rootpath+'pegase/'
    if (file_test(outpath,/dir) eq 0) then spawn, 'mkdir -p '+outpath, /sh
    
    iofilespath = '/tmp/' ; could optionally keep these files
    isedoutpath = '/tmp/'
    
; read the ages.dat file to retrieve the default Pegase model ages
    agesfile = pegpath+'ages.dat'
    if (file_test(agesfile) eq 0) then begin
       splog, 'File ages.dat not found!'
       return
    endif
    readcol, agesfile, ages, format='A', comment='!', /silent
    nage = n_elements(ages)

;   ages = round(range(1,20000,150,/log))
;   ages = ages[uniq(ages,sort(ages))]
;   junk = replicate({ages: 0L},n_elements(ages))
;   junk.ages = ages
;   struct_print, junk, file='ages.dat', /no_hea
    
; read the IMFs file
    imffile = pegpath+'list_IMFs.dat'
    splog, 'Reading '+imffile
    readcol, imffile, imf, format='A', comment='!', /silent
    nimf = n_elements(imf)
    splog, 'Building BASEMODELS for the following IMFs:'
    niceprint, '    '+imf
    splog, "If this is OK, type 'Y' otherwise edit the list_IMFs.dat file:"
;   cc = 'y'
    cc = get_kbrd(1)
    if (strupcase(cc) ne 'Y') then return

    timf = systime(1)
    for ii = 0, nimf-1 do begin

       imfroot = repstr(repstr(imf[ii],'.dat',''),'IMF_','')
       case imfroot of 
          'Salpeter': imfstr = 'salp'
          'Kroupa02': imfstr = 'kroupa'
          else: if strmatch(imfroot,'*cosmicimf*') then begin ; for my AGES/IMF project
             imfstr = imfroot
          endif else message, 'Fix me'
       endcase

; ---------------------------------------------------------------------------
; check whether the SSPs files exist; if not, make them; also make
; them if /REMAKE_SSPs (which you would only need to do if the IMF
; changes)

       sspsfile = outpath+imfroot+'_SSPs.dat'
       if (file_test(sspsfile) eq 0) or keyword_set(remake_ssps) then begin
          splog, 'Building SSPs for IMF: '+imf[ii]
; clean up old files
          allfiles = [file_search(outpath+imfroot+'*_SSPs.dat'),$
            file_search(outpath+imfroot+'*_tracks*.dat')]
          nall = n_elements(allfiles)
          if (nall ne 0) then rmfile, allfiles
          infile = repstr(sspsfile,'.dat','.input')
          logfile = repstr(sspsfile,'.dat','.log')

; we have to read the IMF itself to get the lowest and higest stellar
; mass for this IMF
          junk = djs_readlines(pegpath+imf[ii])
          junk = junk[where(strcompress(junk,/remove) ne '')]
          lomass = strtrim(string((strsplit(junk[1],' ',/extract))[0],format='(F12.2)'),2)
          himass = strtrim(string(junk[n_elements(junk)-1],format='(F12.2)'),2)

          openw, lun, infile, /get_lun            ; SSPs_HR.f input file
          printf, lun, string(ii+1,format='(I0)') ; IMF number
          printf, lun, lomass   ; lowest mass M_sun
          printf, lun, himass   ; highest mass M_sun
          printf, lun, 'B'      ; SNII ejecta model of Woosley & Weaver
          printf, lun, 'y'      ; stellar winds
          printf, lun, '1'      ; Basel stellar library
          printf, lun, imfroot  ; prefix
          printf, lun, 'end'
          free_lun, lun

          pushd, outpath
          spawn, 'SSPs_HR < '+infile+' > '+logfile, /sh ; run SSPs_HR.f
          popd
       endif

; ---------------------------------------------------------------------------
; metallicity grid
       if keyword_set(test) then $
         Z = [0.02] else $
         Z = [0.0004,0.004,0.008,0.02,0.05]
;        Z = [0.0001,0.0004,0.004,0.008,0.02,0.05,0.1]
       Zstr = strtrim(string(Z,format='(F12.4)'),2)
       nZ = n_elements(Z)

; ---------------------------------------------------------------------------
; define the tau-values between 0 Gyr (SSP) and 13.5 Gyr, the maximum
; age of the universe
       if keyword_set(test) then begin
          tau = [3.0D,7.0D]
       endif else begin
          tau1 = dindgen(((13.5D)-0.0D)/(0.5D)+1D0)*(0.5D)+0.0D ; 0.5 Gyr spacing [0-13.5]
          tau2 = dindgen(((20.0D)-14.0D)/(1.0D)+1D0)*(1.0D)+14.0D ; 1 Gyr spacing [14-20]
          tau = [tau1,tau2,100.0D]                                ; add 100 Gyr to simulate continuous SF (<0.1% error)
       endelse
       ntau = n_elements(tau)

       taufile = strarr(ntau,nZ)
       for iZ = 0, nZ-1 do taufile[*,iZ] = imfstr+'_'+Z2string(Z[iZ])+$
         '_tau_'+tau2string(tau)+'.fits'
       if (nZ eq 1) then taufile = reform(taufile,ntau,nZ)

; ---------------------------------------------------------------------------
; define the constant star formation duration parameters; select the
; durations uniformly between 0.1 and 1 Gyr; append a model with
; continuous star formation for 13.5 Gyr to simulate a constant star
; formation model of the entire age of the universe
       if keyword_set(test) then begin
          const = [0.1D,13.5D]  ; the 13.5 needs to be here
       endif else begin
; 0.01-0.1 Gyr
          constmin1 = 0.01D & constmax1 = 0.09D & dconst1 = 0.01D ; [Gyr]
          const1 = dindgen((constmax1-constmin1)/dconst1+1)*$
            dconst1+constmin1
; 0.1-1 Gyr
          constmin2 = 0.1D & constmax2 = 0.9D & dconst2 = 0.1D ; [Gyr]
          const2 = dindgen((constmax2-constmin2)/dconst2+1)*$
            dconst2+constmin2
; 1-10 Gyr
          constmin3 = 1.0D & constmax3 = 10.0D & dconst3 = 1.0D ; [Gyr]
          const3 = dindgen((constmax3-constmin3)/dconst3+1)*$
            dconst3+constmin3
          const = [const1,const2,const3,13.5D] ; add 13.5 Gyr of continuous star formation
       endelse
       nconst = n_elements(const)

       constfile = strarr(nconst,nZ)
       for iZ = 0, nZ-1 do constfile[*,iZ] = imfstr+'_'+Z2string(Z[iZ])+$
         '_const_'+const2string(const)+'.fits'
       if (nZ eq 1) then constfile = reform(constfile,nconst,nZ)

; define the characteristic time for star formation for each model
       allsfhfiles = [taufile,constfile]
       time = [tau,const] 
       time_str = strtrim(string(1E3*time,format='(F12.1)'),2) ; [Myr]
       istau = fix([tau*0+1,const*0])
       nsfh = n_elements(time)

; we want the constant SFR models for 1 M_sun/yr but that causes
; numerical problems; so use 1E-4 M_sun/Myr = 1E-10 M_sun/yr when
; building the grid, and then scale the output quantities by 1E4
       sfr_const = 1D-4
       sfr_const_str = strtrim(string(sfr_const,format='(E12.1)'),2)
    
; ---------------------------------------------------------------------------
; write out an information structure
       info = {$
         imf:               imfstr,$
         Z:                    0.0,$
         tau:           float(tau),$
         const:       float(const),$
         taufile:     strarr(ntau),$
         constfile: strarr(nconst),$
         ages:      float(ages)*1E6} ; [yr]
       info = replicate(info,nZ)
       info.Z = Z
       info.taufile = taufile+'.gz'
       info.constfile = constfile+'.gz'

; write out the info structure
       infofile = rootpath+'info_pegase_'+imfstr+'.fits'
       im_mwrfits, info, infofile, /clobber
       
; ---------------------------------------------------------------------------
; loop on each metallicity
       tz = systime(1)
       for iZ = 0, nZ-1 do begin ; loop on metallicity

; initialize the ISED and input/output file names       
          splog, 'IMF='+imfstr+', Z='+string(Z[iZ],format='(G0.0)')

; --------------------------------------------------
; build the scenarios file
          t0 = systime(1)

          scenarios_infile = outpath+'scenarios_'+imfroot+'.input'        ; scenarios.f input
          scenarios_outfile = repstr(scenarios_infile,'.input','.output') ; scenarios.f output
          scenarios_logfile = repstr(scenarios_infile,'.input','.log')
          spectra_logfile = repstr(scenarios_logfile,'scenarios','spectra')
    
          openw, lun, scenarios_infile, /get_lun ; scenarios.f input file
          printf, lun, scenarios_outfile         ; scenarios.f output file
          printf, lun, sspsfile                  ; SSPs pre-computed by SSPs.f
          printf, lun, '0.05'                    ; default binary fraction
          printf, lun, '1'                       ; Basel stellar library
    
          for jj = 0, nsfh-1 do begin
             sfhfile = outpath+allsfhfiles[jj,iZ]
             if file_test(sfhfile) then spawn, '/bin/rm '+sfhfile, /sh
             printf, lun, sfhfile ; spectra.f output file name
; treat the tau- and continuous star formation scenarios differently 
             if istau[jj] then begin
                if (time[jj] eq 0.0) then begin ; tau=0
                   printf, lun, Zstr[iZ] ; metallicity of the ISM at t=0 (assume "solar")
                   printf, lun, 'n'      ; infall - no
                   printf, lun, '0'      ; star formation scenario: instantaneous burst
                endif else begin
                   printf, lun, Zstr[iZ]     ; metallicity of the ISM at t=0 (for no consistent evolution)
                   printf, lun, 'n'          ; infall - no
                   printf, lun, '2'          ; star formation scenario: e.g., exponentially declining
                   printf, lun, time_str[jj] ; characteristic time for star formation [tau or p1, Myr]
                   printf, lun, '1.0'        ; star formation law pre-factor [p2]
                endelse
             endif else begin ; continuous star formation
                printf, lun, Zstr[iZ]      ; metallicity of the ISM at t=0
                printf, lun, 'n'           ; infall - no
                printf, lun, '1'           ; star formation scenario: continuous star formation
                printf, lun, sfr_const_str ; star formation rate [p1, M_sun/Myr]
                printf, lun, time_str[jj]  ; star formation duration [p2, Myr]
             endelse
             printf, lun, 'n'      ; consistent evolution of the stellar metallicity - no
             printf, lun, Zstr[iZ] ; (fixed!) stellar metallicity
             printf, lun, '0.0'    ; mass fraction of substellar objects formed
             printf, lun, 'n'      ; no galactic winds
             printf, lun, 'y'      ; include nebular emission
             printf, lun, '0'      ; no extinction
;            printf, lun, '2'      ; inclination-averaged extinction for a disk geometry
          endfor 
          printf, lun, 'end'    ; exit
          free_lun, lun

; clean up files    
          if file_test(scenarios_outfile) then spawn, '/bin/rm '+scenarios_outfile, /sh

; generate the models
          t0 = systime(1)
          splog, 'Building the SFH grid for IMF '+imfroot
          spawn, 'scenarios_HR < '+scenarios_infile+' > '+scenarios_logfile, /sh     ; run scenarios.f
          spawn, 'echo "'+scenarios_outfile+'" | '+'spectra_HR > '+spectra_logfile, /sh ; run spectra.f
          splog, format='("Time = ",G0," minutes.")', (systime(1)-t0)/60.0

; now go back through all the models, convert the spectra to an
; isedfit-compatible format, and then overwrite
          for itau = 0, ntau-1 do begin
             peg = im_read_peg(outpath+taufile[itau,iZ])
             ised = peg2isedfit(peg)
             ised = struct_addtags({imf: imfstr, Z: float(Z[iZ]), $
               tau: float(tau[itau])},ised)
             im_mwrfits, ised, outpath+taufile[itau,iZ], /clobber
          endfor
          for iconst = 0, nconst-1 do begin
             peg = im_read_peg(outpath+constfile[iconst,iZ])
             ised = peg2isedfit(peg,sfr_const=sfr_const)
             ised = struct_addtags({imf: imfstr, Z: float(Z[iZ]), $
               const: float(const[iconst])},ised)
             im_mwrfits, ised, outpath+constfile[iconst,iZ], /clobber
          endfor
       endfor ; close metallicity loop
       splog, format='("Time for all metallicities = '+$
         '",G0," minutes.")', (systime(1)-tZ)/60.0
    endfor ; close IMF loop
    splog, format='("Time for all IMFs = '+$
      '",G0," minutes.")', (systime(1)-timf)/60.0
       
return
end
