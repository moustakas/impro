;+
; NAME:
;       NED_WEBGET_PHOTOMETRY
;
; PURPOSE:
;       Collects photometry data from NED.
;
; INPUTS:
;       galaxy - NED-compatible list of galaxies [NGALAXY]
;
; OPTIONAL INPUTS:
;       outfile - output file name if WRITE=1
;
; KEYWORD PARAMETERS:
;       write - write PHOTOMETRY as a binary FITS table
;
; OUTPUTS:
;       photo - output data structure
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;       This routine is a more convenient drop-in replacement routine
;       for NED's batch-retrieval system.  See also
;       PARSE_NED_PHOTOMETRY. 
;
; EXAMPLES:
;       IDL> ned_webget_photometry, 'ngc5194', p
;       IDL> struct_print, p
;
; MODIFICATION HISTORY:
;       J. Moustakas, 2007 Dec 01, NYU - written based on
;          NED_WEBGET_DIAMETERS 
;
; Copyright (C) 2007, John Moustakas
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

pro ned_webget_photometry, galaxy1, photo, outfile=outfile, write=write
    
    ngalaxy = n_elements(galaxy1)
    if (ngalaxy eq 0L) then begin
       doc_library, 'ned_webget_photometry'
       return
    endif

    if (n_elements(outfile) eq 0L) then outfile = 'ned_webget_photometry.fits'
    
; initializez the output structure    

    photo = {$
      galaxy:               '...',  $
      ned_galaxy:           '...',  $
      U_RC3:                -999.0, $
      U_RC3_ERR:            -999.0, $
      U_RC3_FLAG:           '...',  $
      B_RC3:                -999.0, $
      B_RC3_ERR:            -999.0, $
      B_RC3_FLAG:           '...',  $
      V_RC3:                -999.0, $
      V_RC3_ERR:            -999.0, $
      V_RC3_FLAG:           '...',  $

      J_2MASS:              -999.0, $
      J_2MASS_ERR:          -999.0, $
      J_2MASS_REF:          '...',  $
      H_2MASS:              -999.0, $
      H_2MASS_ERR:          -999.0, $
      H_2MASS_REF:          '...',  $
      K_2MASS:              -999.0, $
      K_2MASS_ERR:          -999.0, $
      K_2MASS_REF:          '...',  $

      J20_2MASS:            -999.0, $
      J20_2MASS_ERR:        -999.0, $
      J20_2MASS_REF:        '...',  $
      H20_2MASS:            -999.0, $
      H20_2MASS_ERR:        -999.0, $
      H20_2MASS_REF:        '...',  $
      K20_2MASS:            -999.0, $
      K20_2MASS_ERR:        -999.0, $
      K20_2MASS_REF:        '...',  $

      IRAS_12:              -999.0, $
      IRAS_12_ERR:          -999.0, $
      IRAS_12_REF:          '...',  $
      IRAS_25:              -999.0, $
      IRAS_25_ERR:          -999.0, $
      IRAS_25_REF:          '...',  $
      IRAS_60:              -999.0, $
      IRAS_60_ERR:          -999.0, $
      IRAS_60_REF:          '...',  $
      IRAS_100:             -999.0, $
      IRAS_100_ERR:         -999.0, $
      IRAS_100_REF:         '...',  $

      galex_nuv:            -999.0, $
      galex_nuv_err:        -999.0, $
      galex_nuv_ref:         '...', $
      galex_fuv:            -999.0, $
      galex_fuv_err:        -999.0, $
      galex_fuv_ref:         '...', $

      irac_ch1:             -999.0, $
      irac_ch1_err:         -999.0, $
      irac_ch1_ref:          '...', $
      irac_ch2:             -999.0, $
      irac_ch2_err:         -999.0, $
      irac_ch2_ref:          '...', $
      irac_ch3:             -999.0, $
      irac_ch3_err:         -999.0, $
      irac_ch3_ref:          '...', $
      irac_ch4:             -999.0, $
      irac_ch4_err:         -999.0, $
      irac_ch4_ref:          '...', $

      mips_24:              -999.0, $
      mips_24_err:          -999.0, $
      mips_24_ref:           '...', $
      mips_70:              -999.0, $
      mips_70_err:          -999.0, $
      mips_70_ref:           '...', $
      mips_160:             -999.0, $
      mips_160_err:         -999.0, $
      mips_160_ref:          '...'  $
      }
    photo = replicate(photo,ngalaxy)

    phototags = tag_names(photo[0])
    
    galaxy = strcompress(strupcase(galaxy1),/remove)
    photo.galaxy = galaxy

    for igal = 0L, ngalaxy-1L do begin

       print, format='("Querying NED photometry for galaxy ",I0,"/",I0,".",A10,$)', igal+1L, ngalaxy, string(13b)
       
       nedgalaxy = repstr(galaxy[igal],'+','%2B')
       
       webinfo = webget('http://nedwww.ipac.caltech.edu/cgi-bin/nph-datasearch?objname='+nedgalaxy+$
         '&meas_type=bot&ebars_spec=ebars&label_spec=no&x_spec=freq&y_spec=Fnu_jy&xr=-1&of=ascii_bar&search_type=Photometry')
       text = webinfo.text

       notblank = where(strcompress(text,/remove) ne '',nnotblank)
       if (nnotblank eq 0L) then begin
          splog, 'Galaxy '+galaxy[igal]+' not recognized by NED.'
          continue
       endif
       
       text = text[notblank]
              
       noned = where(strmatch(text,'*not currently recognized*',/fold),nnoned)
       if (nnoned ne 0L) then begin
          splog, 'Galaxy '+galaxy[igal]+' not recognized by NED.'
          continue
       endif

       noned = where(strmatch(text,'*there is no object with this name*',/fold),nnoned)
       if (nnoned ne 0L) then begin
          splog, 'Galaxy '+galaxy[igal]+' not recognized by NED.'
          continue
       endif

       noned = where(strmatch(text,'*no photometric data points*',/fold),nnoned)
       if (nnoned ne 0L) then begin
          splog, 'Galaxy '+galaxy[igal]+' contains no photometry in NED.'
          continue
       endif

; NED name

       line = (where(strmatch(text,'*photometric data for*',/fold) eq 1B,nline))[0]
       if (nline eq 0L) then begin
          splog, 'Galaxy '+galaxy[igal]+' contains no photometry in NED.'
          continue
       endif

       i1 = strpos(text[line],'Photometric data for ')+21
       photo[igal].ned_galaxy = strcompress(strmid(text[line],i1,30),/remove)

; parse the photometry by generating a structure of everything

       ascii_photo = text[line+1L:n_elements(text)-1L]
       ascii_photo = repstr(ascii_photo,'||','| |')
       nphoto = n_elements(ascii_photo)
       
       tags = strsplit(ascii_photo[0],'|',/extract)
       ntags = n_elements(tags)
       data = create_struct(idl_validname(tags[0],/convert_all), '')
       for it = 1L, ntags-1L do data = create_struct(data,$
         idl_validname(tags[it],/convert_all), '')
       data = replicate(data,nphoto)

       for ip = 1L, nphoto-1L do begin
          str = strsplit(ascii_photo[ip],'|',/extract)
          nstr = n_elements(str)
          if (nstr ne ntags) then str = [str,replicate('...',ntags-nstr)]
          for it = 0L, ntags-1L do data[ip-1L].(it) = str[it]
       endfor

       data = im_struct_trimtags(data,select=['OBSERVED_PASSBAND','PHOTOMETRY_MEASUREMENT',$
         'UNCERTAINTY','REFCODE'],newtags=['band','flux','ferr','refcode'])
       
; now parse the bandpasses of interest (see PARSE_NED_PHOTOMETRY)

       allbands = strupcase(strtrim(data.band,2))
       
       thesebands = ['U (U_T)','B (B_T)','V (V_T)',$
         'B (m_B)',$
         'J*2MASS','H*2MASS','K*2MASS',$
         'IRAS 12 microns','IRAS 25 microns','IRAS 60 microns','IRAS 100 microns',$
         '3.6 microns (IRAC)','4.5 microns (IRAC)','5.8 microns (IRAC)','8.0 microns (IRAC)',$
         '24 microns (MIPS)','70 microns (MIPS)','160 microns (MIPS)',$
         'FUV (GALEX)','NUV (GALEX)']
       thesebands = strupcase(thesebands)

       for ib = 0L, n_elements(thesebands)-1L do begin
          match = where(strmatch(allbands,+'*'+thesebands[ib]+'*'),nmatch)
          if (nmatch ge 1L) then if (n_elements(mydata) eq 0L) then $
            mydata = data[match] else mydata = [mydata,data[match]]
       endfor

       if (n_elements(mydata) eq 0L) then begin
          splog, 'Galaxy '+galaxy[igal]+' contains no photometry of interest.'
          continue              ; no photometry of interest
       endif
       
; clean up the "uncertainty" measurements

       for ib = 0L, n_elements(mydata)-1L do begin
          flux = strtrim(mydata[ib].flux,2)
          ferr = repstr(strtrim(mydata[ib].ferr,2),'+/-','')
          if (flux eq '') then flux = '-999.0'
          if (ferr eq '') then ferr = '-999.0'
          if strmatch(ferr,'*%*') then ferr = repstr(ferr,'%','')/100.0
          if strmatch(ferr,'*<*') then begin
             flux = repstr(ferr,'<','-')
             ferr = '-999.0'
          endif
          mydata[ib].flux = flux
          mydata[ib].ferr = ferr
       endfor

;      struct_print, mydata

; now create a map between the bandpass, the reference code, and the
; appropriate PHOTO structure tag

; RC3       
       
       mybands = ['U (U_T)','B (B_T)','V (V_T)'] ; NED-compatible bandpass name
       mytags = ['u','b','v']+'_rc3'
       refcode = '1991RC3.9.C...0000d'

       for ib = 0L, n_elements(mybands)-1L do begin
          match = where(strmatch(mydata.band,'*'+mybands[ib]+'*',/fold) and $
            strmatch(mydata.refcode,'*'+refcode+'*',/fold),nmatch)
          if (nmatch gt 1L) then message, 'Should not happen.'
          if (nmatch eq 1L) then begin
             photo[igal].(where(strupcase(mytags[ib]) eq phototags)) = mydata[match].flux
             photo[igal].(where(strupcase(mytags[ib])+'_ERR' eq phototags)) = mydata[match].ferr
             photo[igal].(where(strupcase(mytags[ib])+'_FLAG' eq phototags)) = mydata[match].band
          endif
       endfor

; add additional B-band measurements, where necessary

       mybands = ['B (m_B)'] ; NED-compatible bandpass name
       mytags = ['b']+'_rc3'
       refcode = '1991RC3.9.C...0000d'

       if (photo[igal].b_rc3 lt -900.0) then begin
          for ib = 0L, n_elements(mybands)-1L do begin
             match = where(strmatch(mydata.band,'*'+mybands[ib]+'*',/fold) and $
               strmatch(mydata.refcode,'*'+refcode+'*',/fold),nmatch)
             if (nmatch gt 1L) then message, 'Should not happen.'
             if (nmatch eq 1L) then begin
                photo[igal].(where(strupcase(mytags[ib]) eq phototags)) = mydata[match].flux
                photo[igal].(where(strupcase(mytags[ib])+'_ERR' eq phototags)) = mydata[match].ferr
                photo[igal].(where(strupcase(mytags[ib])+'_FLAG' eq phototags)) = mydata[match].band
             endif
          endfor
       endif
       
;      help, struct_trimtags(photo[igal],sel=['*galaxy*','*rc3*']), /st 

; IRAS; sometimes the Moshir et al. photometry appears twice, one in
; exponential format and one in floating point; choose the exponential
; one, which contains more significant figures

       mybands = 'IRAS '+['12','25','60','100']+' microns' ; NED-compatible bandpass name
       allrefcode = [$ ; prioritized list
         '2003AJ....126.1607S',$ ; Sanders et al.
         '1988ApJS...68...91R',$ ; Rice et al. 
         '1989AJ.....98..766S',$ ; Soifer et al.
         '1990IRASF.C...0000M']  ; Moshir et al.
       allmyref = [$
         'sanders03',$
         'rice88',$
         'soifer89',$
         'moshir90']
       mytags = 'IRAS_'+['12','25','60','100']

       for ib = 0L, n_elements(mybands)-1L do begin
          for ir = n_elements(allrefcode)-1L, 0L, -1L do begin ; loop in reverse
             match = where(strmatch(mydata.band,'*'+mybands[ib]+'*',/fold) and $
               strmatch(mydata.refcode,'*'+allrefcode[ir]+'*',/fold),nmatch)
             if (nmatch gt 1L) then begin
                keep = where(strmatch(mydata[match].flux,'*E*'))
                match = match[keep] & nmatch = 1L
             endif
             if (nmatch eq 1L) then begin
                photo[igal].(where(strupcase(mytags[ib]) eq phototags)) = mydata[match].flux
                photo[igal].(where(strupcase(mytags[ib])+'_ERR' eq phototags)) = mydata[match].ferr
                photo[igal].(where(strupcase(mytags[ib])+'_REF' eq phototags)) = allmyref[ir]
             endif
          endfor
       endfor

;      help, struct_trimtags(photo[igal],sel='*iras*'), /st 

; MIPS; in the future (jm07dec02nyu), look for Engelbracht et al. 2008
; photometry 

       mybands = ['24','70','160']+' microns (MIPS)' ; NED-compatible bandpass name
       allrefcode = [$ ; prioritized list
         '2007ApJ...655..863D',$ ; Dale et al.  (SINGS)
         '2007AJ....133..791S']  ; Smith et al. (Spirals, Bridges, & Tails)
       allmyref = [$
         'dale07',$
         'smith07']
       mytags = 'mips_'+['24','70','160']

       for ib = 0L, n_elements(mybands)-1L do begin
          for ir = n_elements(allrefcode)-1L, 0L, -1L do begin ; loop in reverse
             match = where(strmatch(mydata.band,'*'+mybands[ib]+'*',/fold) and $
               strmatch(mydata.refcode,'*'+allrefcode[ir]+'*',/fold),nmatch)
             if (nmatch gt 1L) then message, 'Should not happen.'
             if (nmatch eq 1L) then begin
                photo[igal].(where(strupcase(mytags[ib]) eq phototags)) = mydata[match].flux
                photo[igal].(where(strupcase(mytags[ib])+'_ERR' eq phototags)) = mydata[match].ferr
                photo[igal].(where(strupcase(mytags[ib])+'_REF' eq phototags)) = allmyref[ir]
             endif
          endfor
       endfor

;      help, struct_trimtags(photo[igal],sel='*mips*'), /st 
       
; IRAC; in the future (jm07dec02nyu), look for Engelbracht et al. 2008
; photometry 

       mybands = ['3.6','4.5','5.8','8.0']+' microns (IRAC)' ; NED-compatible bandpass name
       allrefcode = [$ ; prioritized list
         '2007ApJ...655..863D',$ ; Dale et al.  (SINGS)
         '2007AJ....133..791S']  ; Smith et al. (Spirals, Bridges, & Tails)
       allmyref = [$
         'dale07',$
         'smith07']
       mytags = 'irac_'+['ch1','ch2','ch3','ch4']

       for ib = 0L, n_elements(mybands)-1L do begin
          for ir = n_elements(allrefcode)-1L, 0L, -1L do begin ; loop in reverse
             match = where(strmatch(mydata.band,'*'+mybands[ib]+'*',/fold) and $
               strmatch(mydata.refcode,'*'+allrefcode[ir]+'*',/fold),nmatch)
             if (nmatch gt 1L) then message, 'Should not happen.'
             if (nmatch eq 1L) then begin
                photo[igal].(where(strupcase(mytags[ib]) eq phototags)) = mydata[match].flux
                photo[igal].(where(strupcase(mytags[ib])+'_ERR' eq phototags)) = mydata[match].ferr
                photo[igal].(where(strupcase(mytags[ib])+'_REF' eq phototags)) = allmyref[ir]
             endif
          endfor
       endfor

;      help, struct_trimtags(photo[igal],sel='*irac*'), /st 
       
; GALEX; in the future (jm07dec02nyu), look for the GALEX UV Atlas of
; Nearby Galaxies by Gil de Paz et al. 2007, ApJS and also the data
; from 11HUGS

       mybands = ['FUV','NUV']+' (GALEX)' ; NED-compatible bandpass name
       allrefcode = [$ ; prioritized list
         '2007ApJ...655..863D'] ; Dale et al.  (SINGS)
       allmyref = [$
         'dale07']
       mytags = 'galex_'+['fuv','nuv']

       for ib = 0L, n_elements(mybands)-1L do begin
          for ir = n_elements(allrefcode)-1L, 0L, -1L do begin ; loop in reverse
             match = where(strmatch(mydata.band,'*'+mybands[ib]+'*',/fold) and $
               strmatch(mydata.refcode,'*'+allrefcode[ir]+'*',/fold),nmatch)
             if (nmatch gt 1L) then begin
;               splog, 'Multiple 2MASS '+strupcase(mytags[ib])+' magnitudes for '+galaxy[igal]
                brighter = min(mydata[match].flux,keep)
                match = match[keep] & nmatch = 1L
             endif
             if (nmatch eq 1L) then begin
                photo[igal].(where(strupcase(mytags[ib]) eq phototags)) = mydata[match].flux
                photo[igal].(where(strupcase(mytags[ib])+'_ERR' eq phototags)) = mydata[match].ferr
                photo[igal].(where(strupcase(mytags[ib])+'_REF' eq phototags)) = allmyref[ir]
             endif
          endfor
       endfor

;      help, struct_trimtags(photo[igal],sel='*galex*'), /st 

; 2MASS "total"

       mybands = ['J','H','K']+'_*tot* (2MASS*)' ; NED-compatible bandpass name
       allrefcode = [$ ; prioritized list
         '2003AJ....125..525J',$ ; Jarrett et al. 2003
         '20032MASX.C.......:']  ; XSC
       allmyref = [$
         'jarrett03',$
         'jarrett01']
       mytags = ['j','h','k']+'_2mass'

       for ib = 0L, n_elements(mybands)-1L do begin
          for ir = n_elements(allrefcode)-1L, 0L, -1L do begin ; loop in reverse
             match = where(strmatch(mydata.band,'*'+mybands[ib]+'*',/fold) and $
               strmatch(mydata.refcode,'*'+allrefcode[ir]+'*',/fold),nmatch)
             if (nmatch gt 1L) then begin
;               splog, 'Multiple 2MASS '+strupcase(mytags[ib])+' magnitudes for '+galaxy[igal]
                brighter = min(mydata[match].flux,keep)
                match = match[keep] & nmatch = 1L
             endif
             if (nmatch eq 1L) then begin
                photo[igal].(where(strupcase(mytags[ib]) eq phototags)) = mydata[match].flux
                photo[igal].(where(strupcase(mytags[ib])+'_ERR' eq phototags)) = mydata[match].ferr
                photo[igal].(where(strupcase(mytags[ib])+'_REF' eq phototags)) = allmyref[ir]
             endif
          endfor
       endfor

; 2MASS: 20th isophotal limit; sometimes multiple XSC photometric
; measurements are returned by NED; see PARSE_NED_PHOTOMETRY for
; details 

       mybands = ['J','H','K']+'_*20* (2MASS*)' ; NED-compatible bandpass name
       allrefcode = [$ ; prioritized list
         '2003AJ....125..525J',$ ; Jarrett et al. 2003
         '20032MASX.C.......:']  ; XSC
       allmyref = [$
         'jarrett03',$
         'jarrett01']
       mytags = ['j20','h20','k20']+'_2mass'

       for ib = 0L, n_elements(mybands)-1L do begin
          for ir = n_elements(allrefcode)-1L, 0L, -1L do begin ; loop in reverse
             match = where(strmatch(mydata.band,'*'+mybands[ib]+'*',/fold) and $
               strmatch(mydata.refcode,'*'+allrefcode[ir]+'*',/fold),nmatch)
             if (nmatch gt 1L) then begin
;               splog, 'Multiple 2MASS '+strupcase(mytags[ib])+' magnitudes for '+galaxy[igal]
                brighter = min(mydata[match].flux,keep)
                match = match[keep] & nmatch = 1L
             endif
             if (nmatch eq 1L) then begin
                photo[igal].(where(strupcase(mytags[ib]) eq phototags)) = mydata[match].flux
                photo[igal].(where(strupcase(mytags[ib])+'_ERR' eq phototags)) = mydata[match].ferr
                photo[igal].(where(strupcase(mytags[ib])+'_REF' eq phototags)) = allmyref[ir]
             endif
          endfor
       endfor

;      help, struct_trimtags(photo[igal],sel='*2mass*'), /st 
       
       delvarx, mydata
;      help, photo[igal], /str & stop
       
    endfor
    print

; write out    
    
    if keyword_set(write) then begin

       splog, 'Writing '+outfile
       mwrfits, photo, outfile, /create
       spawn, 'gzip -f '+outfile, /sh

    endif
    
return
end
    
;;; IRAS
;;
;;       mybands = 'IRAS '+['12','25','60','100']+' MICRONS' ; NED-compatible bandpass name
;;       allname = ['SANDERS','RICE','SOIFER','MOSHIR'] ; the order must match ALLREFCODE
;;       allrefcode = ['2003AJ....126.1607S','1988ApJS...68...91R',$
;;         '1989AJ.....98..766S','1990IRASF.C...0000M']
;;       allmytags = 'IRAS_'+['12','25','60','100']
;;
;;       for in = 0L, n_elements(allname)-1L do begin
;;
;;          name = allname[in]
;;          mytags = name+'_'+allmytags
;;          refcode = allrefcode[in]
;;          
;;          for ib = 0L, n_elements(mybands)-1L do begin
;;             match = where(strmatch(mydata.band,'*'+mybands[ib]+'*',/fold) and $
;;               strmatch(mydata.refcode,'*'+refcode+'*',/fold),nmatch)
;;             if (nmatch gt 1L) then message, 'Should not happen.'
;;             if (nmatch eq 1L) then begin
;;                photo[igal].(where(strupcase(mytags[ib]) eq phototags)) = mydata[match].flux
;;                photo[igal].(where(strupcase(mytags[ib])+'_ERR' eq phototags)) = mydata[match].ferr
;;             endif
;;          endfor
;;
;;       endfor
;;       
