;+
; NAME:
;   CREATE_PEGASE_SSPS
;
; PURPOSE:
;   Pegase does not come with any SSPs, so build a standard set for a
;   variety of IMFs.  Use IM_READ_PEGASE() to read these files. 
;
; INPUTS: 
;
; KEYWORD PARAMETERS: 
;   cosmic_imf - generate SSPs for my PRIMUS cosmic-IMF project 
;
; OUTPUTS: 
;
; COMMENTS:
;   To change the ages.dat file (must be integers!)
;      pushd, '$PEGASE_HR_DIR/data/user_defined/'
;      ages = round([0D,range(1D,15D3,149,/log)])
;      ages = ages[uniq(ages)]
;      openw, lun, 'ages.dat', /get_lun
;      for ii = 0, n_elements(ages)-1 do printf, lun, ages[ii]
;      free_lun, lun
; 
;   For now do not include nebular emission, but that is on my ToDo
;   list. 
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2011 Mar 29, UCSD
;
; Copyright (C) 2011, John Moustakas
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

pro create_pegase_ssps, cosmic_imf=cosmic_imf, dossps=dossps, $
  doscenarios=doscenarios

    hrpath = getenv('PEGASE_HR_DIR')+'/data/user_defined/'
    ssppath = getenv('PEGASE_HR_DIR')+'/SSPs/'
    if (file_test(ssppath,/dir) eq 0) then $
      spawn, 'mkdir -p '+ssppath, /sh

; --------------------------------------------------    
; [1] specify (and create if necessary) the IMFs 

; special case for my cosmic IMF project; consider a grid of high-mass
; slopes with a fixed low-mass slope 
    if keyword_set(cosmic_imf) then begin
       alpha2 = primus_cosmicimf_slope() ; high-mass slope
       mass = [0.1,0.5,120.0] ; stellar mass ranges

       imf = 'cosmicimf_'+string(alpha2,format='(F4.2)')
       imfstr = imf
       nimf = n_elements(imf)
       minmass = replicate(mass[0],nimf)
       maxmass = replicate(mass[2],nimf)

       imffile = hrpath+'IMF_'+imf+'.dat' 
       for ii = 0, nimf-1 do begin 
          slope = 1.0-[1.3,alpha2[ii]] ; fixed low-mass slope
          write_pegase_imf, mass, slope, outfile=imffile[ii]
       endfor
    endif else begin 
       imf = ['Salpeter_100','Kroupa01_100'] ; fiducial set
       imfstr = ['salp','kroupa01']
       minmass = [0.1,0.1]
       maxmass = [100,100]
       imffile = hrpath+'IMF_'+imf+'.dat' 
       nimf = n_elements(imf)

       write_pegase_imf, [0.1,100.0], 1.0-2.35, outfile=imffile[0]          ; Salpeter
       write_pegase_imf, [0.1,0.5,100.0], 1.0-[1.3,2.3], outfile=imffile[1] ; Kroupa01
    endelse

; update the list_IMFs.dat file; be sure a copy of the *original* IMF
; file exists because we are going to be overwriting it
    if (file_test(hrpath+'list_IMFs.original.dat') eq 0) then begin
       splog, 'Please manually copy/backup the *original* list_IMF.dat file to: '
       splog, '   '+hrpath+'list_IMFs.original.dat'
       return
    endif
     
    outfile = hrpath+'list_IMFs.dat'
    splog, 'Writing '+outfile
    openw, lun, outfile, /get_lun
    printf, lun, '! File written by CREATE_PEGASE_SSPs on '+im_today()
    for ii = 0, nimf-1 do printf, lun, file_basename(imffile[ii])
    free_lun, lun

; --------------------------------------------------    
; [2] convolve the stellar models with each IMF with a default set of
; optional parameters, if necessary
    if keyword_set(dossps) then begin
       pushd, ssppath
       for ii = 0, nimf-1 do begin
          splog, 'Working on IMF '+imf[ii]
; clean up old files
          allfiles = file_search(ssppath+imf[ii]+'_*.dat',count=nall)
          if (nall ne 0) then rmfile, allfiles
          sspsfile = imf[ii]+'_SSPs.dat'
          infile = repstr(sspsfile,'.dat','.input')
          logfile = repstr(sspsfile,'.dat','.log')
          
          openw, lun, infile, /get_lun         ; SSPs_HR.f input file
          printf, lun, string(ii+1,format='(I0)') ; IMF number
          printf, lun, strtrim(minmass[ii],2)     ; lowest mass M_sun
          printf, lun, strtrim(maxmass[ii],2)     ; highest mass M_sun
          printf, lun, 'B'                        ; SNII ejecta model of Woosley & Weaver
          printf, lun, 'y'                        ; stellar winds
          printf, lun, '1'                        ; BaSeL stellar library
          printf, lun, imf[ii]                    ; prefix
          printf, lun, 'end'
          free_lun, lun
          
          t0 = systime(1)
          splog, 'Building SSPs for IMF: '+imf[ii]
;         pushd, ssppath
          spawn, 'SSPs_HR < '+infile+' > '+logfile, /sh ; run SSPs_HR.f
;         popd
          splog, format='("Time = ",G0," minutes")', (systime(1)-t0)/60.0
       endfor
       popd
    endif

; --------------------------------------------------    
; [3] call scenarios with tau=0 to build the actual SSPs
    if keyword_set(doscenarios) then begin
;      Zgrid = ['0.02']
       Zgrid = ['0.0001','0.0004','0.004','0.008','0.02','0.05','0.1']
       nZ = n_elements(Zgrid)

       pushd, ssppath
       for ii = 0, nimf-1 do begin
          sspsfile = imf[ii]+'_SSPs.dat'
          if (file_test(sspsfile) eq 0) then begin
             splog, sspsfile+' not found!'
             splog, 'Rerun with /DOSSPS!'
             continue
          endif

          scenarios_infile = 'scenarios_'+imf[ii]+'.input'        ; scenarios.f input
          scenarios_outfile = repstr(scenarios_infile,'.input','.output') ; scenarios.f output
          scenarios_logfile = repstr(scenarios_infile,'.input','.log')
          spectra_logfile = repstr(scenarios_logfile,'scenarios','spectra')
    
          openw, lun, scenarios_infile, /get_lun ; scenarios.f input file
          printf, lun, scenarios_outfile         ; scenarios.f output file
          printf, lun, sspsfile                  ; SSPs pre-computed by SSPs.f
          printf, lun, '0.05'                    ; default binary fraction
          printf, lun, '1'                       ; BaSeL stellar library
    
          for iZ = 0, nZ-1 do begin ; loop on each metallicity
             sfhfile = 'SSP_'+imfstr[ii]+'_Z'+Zgrid[iZ]+'.fits'
             if file_test(sfhfile) then rmfile, sfhfile
             printf, lun, sfhfile ; spectra.f output file name
             printf, lun, Zgrid[iZ]  ; metallicity of the ISM at t=0
             printf, lun, 'n'       ; infall - no
             printf, lun, '0'       ; star formation scenario: instantaneous burst
             printf, lun, 'n'       ; consistent evolution of the stellar metallicity - no
             printf, lun, Zgrid[iZ] ; (fixed!) stellar metallicity
             printf, lun, '0.0'     ; mass fraction of substellar objects formed
             printf, lun, 'n'       ; no galactic winds
             printf, lun, 'n'       ; include nebular emission - not for now
             printf, lun, '0'       ; no extinction
;            printf, lun, '2'       ; inclination-averaged extinction for a disk geometry
          endfor 
          printf, lun, 'end'    ; exit
          free_lun, lun

; clean up files and then generate the models
          if file_test(scenarios_outfile) then spawn, '/bin/rm '+scenarios_outfile, /sh

          t0 = systime(1)
          splog, 'Building the SSPs for IMF '+imf[ii]
          spawn, 'scenarios_HR < '+scenarios_infile+' > '+scenarios_logfile, /sh     ; run scenarios.f
          spawn, 'echo "'+scenarios_outfile+'" | '+'spectra_HR > '+spectra_logfile, /sh ; run spectra.f
          splog, format='("Time = ",G0," minutes.")', (systime(1)-t0)/60.0
       endfor 
    endif 
    
return
end
    
