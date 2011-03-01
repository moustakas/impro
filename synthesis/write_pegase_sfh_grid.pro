; ###########################################################################
; THIS VERSION WORKS WITH PEGASE-HR, WHICH SUPERSEDES PEGASE.2 (see
; WRITE_PEGASE_SFH_GRID_LOWRES) 
; ###########################################################################
;+
; NAME:
;       WRITE_PEGASE_SFH_GRID
;
; PURPOSE:
;       Generate a star formation history grid of Pegase models.
;
; INPUTS:
;
; OPTIONAL INPUTS:
;       gridpath - output path [default getenv('PEGASE_HR_SFHGRID_DIR')] 
;
; KEYWORD PARAMETERS:
;       dotau    - generate the tau models
;       doinfall - generate the infall models
;       measure  - measure quantities of interest
;
; OUTPUTS:
;       Models and an information structure are written to GRIDPATH. 
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;       Note that MEASURE_PEGASE_SFH_GRID is run at the end. 
;
; MODIFICATION HISTORY:
;       J. Moustakas, 2006 Apr 19, U of A - based on
;          WRITE_SFH_BASE_MODELS 
;       jm07feb26nyu - updated to use PEGASE-HR
;
; Copyright (C) 2006-2007, John Moustakas
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

pro write_pegase_sfh_grid, elodie=elodie, dotau=dotau, doinfall=doinfall, $
  measure=measure, compute_ssps=compute_ssps

    gridpath = getenv('PEGASE_HR_SFHGRID_DIR')+'/' ; SFH grid path

    if keyword_set(elodie) then begin
       library = '2'            ; 2=ELODIE
       library_str = 'elodie'
    endif else begin
       library = '1'            ; 1=Basel
       library_str = 'basel'
    endelse

    imf = 'salp'
    sspsfile = gridpath+imf+'_'+library_str+'_SSPs.dat'

    zsun = 0.0122 ; the tau=0 model starts out with solar metallicity
    
; ---------------------------------------------------------------------------
; run SSPs_HR - compute for both libraries
; ---------------------------------------------------------------------------

    if keyword_set(compute_ssps) then begin

       lib = ['1','2']
       lib_str = ['basel','elodie']
       
       for ilib = 0L, 1L do begin

          temp_sspsfile = gridpath+imf+'_'+lib_str[ilib]+'_SSPs.dat'
          
          if file_test(temp_sspsfile,/regular) then begin
             splog, 'Deleting '+temp_sspsfile+'.'
             spawn, ['/bin/rm '+temp_sspsfile], /sh
          endif

          ssps_infile = repstr(temp_sspsfile,'.dat','.input')
          ssps_logfile = repstr(temp_sspsfile,'.dat','.logfile')
          
          openw, lun, ssps_infile, /get_lun  ; SSPs_HR.f input file
          printf, lun, '4'                   ; Salpeter IMF
          printf, lun, '0.1'                 ; 0.1 M_sun
          printf, lun, '100.0'               ; 100.0 M_sun
          printf, lun, 'B'                   ; SNII ejecta model of Woosley & Weaver
          printf, lun, 'y'                   ; stellar winds
          printf, lun, lib[ilib]             ; stellar library
          printf, lun, imf+'_'+lib_str[ilib] ; prefix
          printf, lun, 'end'
          free_lun, lun

          t0 = systime(1)
          splog, 'Generating the SSPs for stellar library '+strupcase(lib_str[ilib])+'.'
          spawn, ['SSPs_HR < '+ssps_infile+' > '+ssps_logfile], /sh ; run SSPs_HR.f
          splog, format='("Time to generate the SSPs = ",G0," minutes.")', (systime(1)-t0)/60.0

       endfor
       return
       
    endif

    if (file_test(sspsfile,/regular) eq 0L) then begin
       splog, sspsfile+' not found!'
       return
    endif
    
; ---------------------------------------------------------------------------
; define the tau values
; ---------------------------------------------------------------------------

;   taumin = 0.0D & taumax = 20.0D & dtau = 0.5D ; [Gyr]
;   ntau = fix((taumax-taumin)/dtau)+1
;   tau = dindgen(ntau)*dtau+taumin

    tau = [0.0,0.5,1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0,15.0,999.0] ; [Gyr]
    ntau = n_elements(tau)

    taustring = strarr(ntau)
    flag = where((tau lt 10.0),nflag)
    if (nflag ne 0L) then taustring[flag] = '00'+string(tau[flag],format='(F3.1)')+'Gyr'
    flag = where((tau ge 10.0) and (tau lt 100.0),nflag)
    if (nflag ne 0L) then taustring[flag] = '0'+string(tau[flag],format='(F4.1)')+'Gyr'
    flag = where((tau ge 100.0),nflag)
    if (nflag ne 0L) then taustring[flag] = string(tau[flag],format='(F5.1)')+'Gyr'

    if (max(tau) ge 1000.0) then message, 'TAU value not supported!'

; ---------------------------------------------------------------------------
; define the infall values
; ---------------------------------------------------------------------------

;   tinfall = [1.0,5.0] ; [Gyr]
;   tinfall = [0.5,1.0,2.0,3.0,5.0,10.0] ; [Gyr]
    tinfall = [0.5,1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0,15.0,100.0] ; [Gyr]
    ntinfall = n_elements(tinfall)

    tinfallstring = strarr(ntinfall)
    flag = where((tinfall lt 10.0),nflag)
    if (nflag ne 0L) then tinfallstring[flag] = '00'+string(tinfall[flag],format='(F3.1)')+'Gyr'
    flag = where((tinfall ge 10.0) and (tinfall lt 100.0),nflag)
    if (nflag ne 0L) then tinfallstring[flag] = '0'+string(tinfall[flag],format='(F4.1)')+'Gyr'
    flag = where((tinfall ge 100.0),nflag)
    if (nflag ne 0L) then tinfallstring[flag] = string(tinfall[flag],format='(F5.1)')+'Gyr'

    if (max(tinfall) ge 1000.0) then message, 'TINFALL value not supported!'

; ---------------------------------------------------------------------------
; define the Kennicutt law efficiencies
; ---------------------------------------------------------------------------

; eta is analogous to p2 [Myr/M_sun] on p. 10 of the Pegase
; documentation; e.g., 1000 [Myr/M_sun] = 0.1 (10%) per 10^8 yr 

    eta = [0.01D,0.1D,1.0D] ; efficiency per 10^8 yr
    p2 = 1.0/(eta*1D-2)     ; inverse efficiency per Myr
    neta = n_elements(eta)

    etastring = ['0.01','0.10','1.00'] ; per 10^8 yr

; --------------------------------------------------
; generate the tau models with no infall
; --------------------------------------------------

    if keyword_set(dotau) then begin
    
       pegfile = gridpath+imf+'_tau_'+taustring+'.fits'
    
       scenarios_infile = gridpath+'scenarios_'+imf+'_tau.input'        ; scenarios.f input
       scenarios_outfile = repstr(scenarios_infile,'.input','.output')  ; scenarios.f output
       scenarios_logfile = repstr(scenarios_infile,'.input','.logfile')
       spectra_logfile = repstr(scenarios_logfile,'scenarios','spectra')
    
       openw, lun, scenarios_infile, /get_lun  ; scenarios.f input file
       printf, lun, scenarios_outfile          ; scenarios.f output file
       printf, lun, sspsfile                   ; SSPs pre-computed by SSPs.f
       printf, lun, '0.05'                     ; default binary fraction
       printf, lun, library                    ; stellar library
    
       for itau = 0L, ntau-1L do begin
    
;         splog, pegfile[itau]
          if file_test(pegfile[itau],/regular) then begin
             splog, 'Deleting '+pegfile[itau]+'.'
             spawn, ['/bin/rm '+pegfile[itau]], /sh
          endif
          
          tau_str = strtrim(string(1E3*tau[itau],format='(F12.1)'),2) ; [Myr]
    
          printf, lun, pegfile[itau] ; spectra.f output file name
          case tau[itau] of
             0.0:   begin
                printf, lun, string(zsun,format='(G0.0)') ; metallicity of the ISM at t=0 (solar)
                printf, lun, 'n'   ; infall - no
                printf, lun, '0'   ; star formation scenario, instantaneous burst
             end
             999.0: begin
; these values will yield a ~1 M_sun galaxy after 20 Gyr
                printf, lun, '0.0'    ; metallicity of the ISM at t=0
                printf, lun, 'n'      ; infall - no
                printf, lun, '1'      ; star formation scenario, constant star formation
                printf, lun, '5.0E-5' ; star formation rate [p1, M_sun/Myr]
                printf, lun, '20000.0'; star formation duration [p2, Myr]
             end
             else: begin
                printf, lun, '0.0'    ; metallicity of the ISM at t=0
                printf, lun, 'n'      ; infall - no
                printf, lun, '2'      ; star formation scenario (e.g., exponentially declining)
                printf, lun, tau_str  ; characteristic time for star formation [tau or p1, Myr]
                printf, lun, '1.0'    ; star formation law pre-factor [p2]
             end
          endcase
          printf, lun, 'y'           ; consistent evolution of the stellar metallicity
          printf, lun, '0.0'         ; mass fraction of substellar objects formed
          printf, lun, 'n'           ; no galactic winds
          printf, lun, 'y'           ; include nebular emission
          printf, lun, '0'           ; no extinction
;         printf, lun, '2'           ; inclination-averaged extinction for a disk geometry
       
       endfor
       printf, lun, 'end'          ; exit
       free_lun, lun

; clean up files    
    
       if file_test(scenarios_outfile,/regular) then begin
          splog, 'Deleting '+scenarios_outfile+'.'
          spawn, ['/bin/rm '+scenarios_outfile], /sh
       endif

; generate the models

       t0 = systime(1)
       splog, 'Generating the tau models with no infall.'
       spawn, ['scenarios_HR < '+scenarios_infile+' > '+scenarios_logfile], /sh     ; run scenarios.f
       spawn, ['echo "'+scenarios_outfile+'" | '+'spectra_HR > '+spectra_logfile], /sh ; run spectra.f
       splog, format='("Time to generate tau models with no infall = ",G0," minutes.")', (systime(1)-t0)/60.0

    endif
       
; --------------------------------------------------
; generate the infall models with a Kennicutt SF law
; --------------------------------------------------

    if keyword_set(doinfall) then begin
       
       scenarios_infile = gridpath+'scenarios_'+imf+'_infall.input' ; scenarios.f input
       scenarios_outfile = repstr(scenarios_infile,'.input','.output') ; scenarios.f output
       scenarios_logfile = repstr(scenarios_infile,'.input','.logfile')
       spectra_logfile = repstr(scenarios_logfile,'scenarios','spectra')

       openw, lun, scenarios_infile, /get_lun ; scenarios.f input file
       printf, lun, scenarios_outfile         ; scenarios.f output file
       printf, lun, sspsfile                  ; SSPs pre-computed by SSPs.f
       printf, lun, '0.05'                    ; default binary fraction
       printf, lun, library                   ; stellar library

       for ieta = 0L, neta-1L do begin          

          p2_str = strtrim(string(p2[ieta],format='(F12.1)'),2)

          for itinfall = 0L, ntinfall-1L do begin

             tinfall_str = strtrim(string(1E3*tinfall[itinfall],format='(F12.1)'),2) ; [Myr]
             pegfile = gridpath+imf+'_kennlaw_'+etastring[ieta]+'_infall_'+tinfallstring[itinfall]+'.fits'

;            splog, pegfile
             if file_test(pegfile,/regular) then begin
                splog, 'Deleting '+pegfile+'.'
                spawn, ['/bin/rm '+pegfile], /sh
             endif

             printf, lun, pegfile     ; spectra.f output file name
             printf, lun, '0.0'       ; metallicity of the ISM at t=0
             printf, lun, 'y'         ; infall - yes
             printf, lun, tinfall_str ; infall time scale [Myr]
             printf, lun, '0.0'       ; metallicity of the infalling gas
             printf, lun, '3'         ; star formation scenario, Kennicutt Law
             printf, lun, '1.0'       ; Kennicutt law exponent
             printf, lun, p2_str      ; Kennicutt law efficiency factor
             printf, lun, 'y'         ; consistent evolution of the stellar metallicity
             printf, lun, '0.0'       ; mass fraction of substellar objects formed
             printf, lun, 'n'         ; no galactic winds
             printf, lun, 'y'         ; include nebular emission
             printf, lun, '0'         ; no extinction
;            printf, lun, '2'         ; inclination-averaged extinction for a disk geometry

          endfor
          
       endfor
       
       printf, lun, 'end'
       free_lun, lun

; clean up files    
       
       if file_test(scenarios_outfile,/regular) then begin
          splog, 'Deleting '+scenarios_outfile+'.'
          spawn, ['/bin/rm '+scenarios_outfile], /sh
       endif

; generate the models

       t0 = systime(1)
       splog, 'Generating the Kennicutt-Schmidt models with infall.'
       spawn, ['scenarios_HR < '+scenarios_infile+' > '+scenarios_logfile], /sh     ; run scenarios.f
       spawn, ['echo "'+scenarios_outfile+'" | '+'spectra_HR > '+spectra_logfile], /sh ; run spectra.f
       splog, format='("Time to generate tau models with infall = ",G0," minutes.")', (systime(1)-t0)/60.0

    endif
       
; --------------------------------------------------
; now measure everything
; --------------------------------------------------

    if keyword_set(measure) then measure_pegase_sfh_grid, /write
    
return
end
