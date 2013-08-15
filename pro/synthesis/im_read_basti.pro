;+
; NAME:
;   IM_READ_BASTI()
;
; PURPOSE:
;   Read a BaSTI SSP into a structure.
;
; INPUTS:
;   None required.
;
; OPTIONAL INPUTS:
;   metallicity - stellar metallicity (default 'zsun')
;     z104 -> Z = 0.0001
;     z304 -> Z = 0.0003
;     z604 -> Z = 0.0006
;     z103 -> Z = 0.001
;     z203 -> Z = 0.002
;     z403 -> Z = 0.004
;     z803 -> Z = 0.008
;     z102 -> Z = 0.01
;     zsun -> Z = 0.0198 (default)
;     z302 -> Z = 0.03
;     z402 -> Z = 0.04
;
;   eta - mass loss rate, either '04' (eta=0.4; normal), or '02'
;     (eta=0.2; low)
;
; KEYWORD PARAMETERS:
;   hires - read the theoretical high-resolution models (default is to
;     read the low-resolution models, which have greater wavelength
;     coverage)  
;   enhanced - read the alpha-enhanced models (default is to use the
;     solar-scaled models)
;
; OUTPUTS:
;   basti - output data structure
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   See http://albione.oa-teramo.inaf.it for additional relevant
;   details. 
; 
;   Only the Kroupa+01 IMF is available. 
;
; File name key:
;   solar scaled models:
;     eta=0.4
;       sso - scaled solar, eta=0.4, with overshooting (for t<11 Gyr)
;       sss - scaled solar, eta=0.4, without overshooting (for t>11 Gyr)
;     eta=0.2
;       so2 - scaled solar, eta=0.2, with overshooting (for t<11 Gyr)
;       ss2 - scaled solar, eta=0.2, without overshooting (for t>11 Gyr)
;
;   alpha-enhanced models:
;     eta=0.4
;       aeo - alpha enhanced, eta=0.4, with overshooting (for t<14 Gyr)
;       aes - elpha enhanced, eta=0.4, without overshooting (for t>14 Gyr)
;     eta=0.2
;       ao2 - alpha enhanced, eta=0.2, with overshooting (for t<14 Gyr)
;       ae2 - elpha enhanced, eta=0.2, without overshooting (for t>14 Gyr)
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2011 Jan 23, UCSD
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

function im_read_basti, metallicity=metallicity, eta=eta, hires=hires, $
  enhanced=enhanced, silent=silent

    ssppath = getenv('IM_RESEARCH_DIR')+'/synthesis/basti/'

    if (n_elements(metallicity) eq 0) then metallicity = 'zsun'
    case metallicity of
       'z104': begin
          zz = 0.0001
          yy = 'y245'
       end
       'z304': begin
          zz = 0.0003
          yy = 'y245'
       end
       'z604': begin
          zz = 0.0006
          yy = 'y246'
       end
       'z103': begin
          zz = 0.001
          yy = 'y246'
       end
       'z203': begin
          zz = 0.002
          yy = 'y248'
       end
       'z403': begin
          zz = 0.004
          yy = 'y251'
       end
       'z803': begin
          zz = 0.008
          yy = 'y256'
       end
       'z102': begin
          zz = 0.01
          yy = 'y259'
       end
       'zsun': begin
          zz = 0.0198
          yy = 'ysun' ; y=0.2734
       end
       'z302': begin
          zz = 0.03
          yy = 'y288'
       end
       'z402': begin
          zz = 0.04
          yy = 'y303'
       end
       else: begin
          splog, 'Metallicity '+metallicity+' not recognized!'
          return, -1
       end
    endcase
    
    if (n_elements(eta) eq 0) then eta = '04'
    if (eta ne '04') and (eta ne '02') then begin
       splog, 'ETA parameter not recognized'
       return, -1
    endif

    imf = 'Kroupa01'
    if keyword_set(hires) then begin
       prefix = 'HRSPEC' 
       wavefactor = 1.0
       fluxfactor = 1.0
    endif else begin
       prefix = 'SPEC'
       wavefactor = 10.0 ; convert micron-->A
       fluxfactor = 1D/4.3607D-33/1D10 ; convert to [erg/s/A]
    endelse

    if keyword_set(enhanced) then begin
       suffix = 'ae'
       if eta eq '04' then models = ['aeo','aes'] else $
         models = ['ao2','as2']
    endif else begin
       suffix = 'ss'
       if eta eq '04' then models = ['sso','sss'] else $
         models = ['so2','ss2']
    endelse

    massfile = ssppath+'MLR'+metallicity+yy+suffix+'_e'+eta+'.dat'
    sspfile = file_search(ssppath+prefix+metallicity+$
      yy+models+'_agb.t6?????',count=nage)
    if (nage eq 0) then message, 'Problem finding files'

    age = strmid(file_basename(sspfile),4,/rev)*1D-3 ; [Gyr]

; read the mass file
    if (file_test(massfile) eq 0) then begin
       splog, 'Mass file '+massfile+' not found!'
       return, -1
    endif
    readfast, massfile, allmass, skip=1
    
; read all the files and pack into a data structure; note the native
; flux units for the hires models are arbitrary, while for the lores
; models the units are 4.3607D-33 erg/s/m!
    for ii = 0L, nage-1 do begin
       if keyword_set(silent) eq 0 then splog, 'Reading '+sspfile[ii]
       readfast, sspfile[ii], data, skip=0
       npix = n_elements(data[0,*])
       if (ii eq 0) then begin
          basti = {$
            imf:                imf,$
            Z:                   zz,$
            age:            age*1D9,$ ; [yr]
            mstar:     fltarr(nage),$ ; [Msun]
            wave: wavefactor*reform(data[0,*]),$ ; [A]
            flux:      fltarr(npix,nage)} ; [erg/s/A]
       endif
       basti.flux[*,ii] = fluxfactor*reform(data[1,*])
;      if ii eq 20 then stop
; get the total stellar mass (living stars plus remnants)  
       basti.mstar = interpol(total(allmass[10:13,*],1),allmass[0,*]*1D6,basti.age)
    endfor

return, basti
end
