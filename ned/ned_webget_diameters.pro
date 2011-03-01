;+
; NAME:
;       NED_WEBGET_DIAMETERS
;
; PURPOSE:
;       Collect size and position angle information from NED.
;
; INPUTS:
;       galaxy - NED-compatible list of galaxies [NGALAXY]
;
; OPTIONAL INPUTS:
;       outfile - output file name if WRITE=1
;
; KEYWORD PARAMETERS:
;       write - write DIAMETERS as a binary FITS table
;
; OUTPUTS:
;       diameters - output data structure
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;       This routine is a work-around to the fact that NED does not
;       return diameters in batch mode.
;
; EXAMPLES:
;
; MODIFICATION HISTORY:
;       J. Moustakas, 2005 Jul 22, U of A
;       jm06jan19uofa - if multiple outer MCG diameters, select the
;                       first one 
;       jm06feb09uofa - also return the 2MASS K20 isophotal diameters
;
; Copyright (C) 2005-2006, John Moustakas
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

pro ned_webget_diameters, galaxy, diameters, outfile=outfile, write=write
    
    ngalaxy = n_elements(galaxy)
    if (ngalaxy eq 0L) then begin
       doc_library, 'ned_webget_diameters'
       return
    endif

    if (n_elements(outfile) eq 0L) then outfile = 'ned_webget_diameters.fits'
    
; initializez the output structure    
    
    diameters = {$
      galaxy:                         ' ', $
      ned_galaxy:                     ' ', $
      rc3_reference:                  ' ', $
      rc3_major_axis:              -999.0, $
      rc3_minor_axis:              -999.0, $
      rc3_posangle:                -999.0, $
      rc3_posangle_equinox:           ' ', $
      rc3_diameter_reflevel:          ' ', $
      rc3_diameter_frequency:         ' ', $
      rc3_measured_quantity:          ' ', $
      rc3_ned_comment:                ' ', $
                                   
      eso_reference:                  ' ', $
      eso_major_axis:              -999.0, $
      eso_minor_axis:              -999.0, $
      eso_posangle:                -999.0, $
      eso_posangle_equinox:           ' ', $
      eso_diameter_reflevel:          ' ', $
      eso_diameter_frequency:         ' ', $
      eso_measured_quantity:          ' ', $
      eso_ned_comment:                ' ', $
                                   
      mcg_reference:                  ' ', $
      mcg_major_axis:              -999.0, $
      mcg_minor_axis:              -999.0, $
      mcg_posangle:                -999.0, $
      mcg_posangle_equinox:           ' ', $
      mcg_diameter_reflevel:          ' ', $
      mcg_diameter_frequency:         ' ', $
      mcg_measured_quantity:          ' ', $
      mcg_ned_comment:                ' ', $
                                   
      ugc_reference:                  ' ', $
      ugc_major_axis:              -999.0, $
      ugc_minor_axis:              -999.0, $
      ugc_posangle:                -999.0, $
      ugc_posangle_equinox:           ' ', $
      ugc_diameter_reflevel:          ' ', $
      ugc_diameter_frequency:         ' ', $
      ugc_measured_quantity:          ' ', $
      ugc_ned_comment:                ' ', $
                                   
      twomass_reference:              ' ', $
      twomass_major_axis:          -999.0, $
      twomass_minor_axis:          -999.0, $
      twomass_k20_major_axis:      -999.0, $
      twomass_k20_minor_axis:      -999.0, $
      twomass_posangle:            -999.0, $
      twomass_posangle_equinox:       ' ', $
      twomass_diameter_reflevel:      ' ', $
      twomass_k20_diameter_reflevel:  ' ', $
      twomass_diameter_frequency:     ' ', $
      twomass_measured_quantity:      ' ', $
      twomass_ned_comment:            ' '}
    diameters = replicate(diameters,ngalaxy)
    diameters.galaxy = galaxy

    refinfo_template = {$
      reference:             ' ', $
      frequency_targeted:    ' ', $
      major_axis:         -999.0, $
      major_axis_units:      ' ', $
      minor_axis:         -999.0, $
      minor_axis_units:      ' ', $
      axis_ratio:         -999.0, $
      posangle:           -999.0, $
      posangle_equinox:      ' ', $
      diameter_reflevel:     ' ', $
      diameter_frequency:    ' ', $
      detector_type:         ' ', $
      fitting_technique:     ' ', $
      measured_quantity:     ' ', $
      measurement_qualifiers:' ', $
      data_qualifiers:       ' ', $
      ned_comment:           ' '}

    for igalaxy = 0L, ngalaxy-1L do begin

       nedgalaxy = repstr(galaxy[igalaxy],'+','%2B')
      
       webinfo = webget('http://nedwww.ipac.caltech.edu/cgi-bin/nph-datasearch?objname='+$
         nedgalaxy+'&search_type=Diameters&search_type=Diameters')
       text = webinfo.text

       noned = where(strmatch(text,'*not currently recognized*',/fold),nnoned)
       if (nnoned ne 0L) then begin
          splog, 'Object '+galaxy[igalaxy]+' not recognized by NED.'
          continue
       endif

       istart = where(strmatch(text,'*reference code*:*',/fold) eq 1B,nstart)
       iend = where(strmatch(text,'*NED Comment*',/fold) eq 1B,nend)
       if (nstart ne nend) then message, 'Problem, Dave.'

       if (nstart eq 0L) then begin
          splog, galaxy[igalaxy]+': no diameter data.'
          continue
       endif
       
; remove particular references, based on the NED comment
       
       retain = where($
;        (strmatch(text[iend],'*Diameters flagged uncertain*',/fold) eq 0B) and $
         (strmatch(text[iend],'*From magnitude-aperture curves*',/fold) eq 0B) and $
         (strmatch(text[iend],'*Extinction and inclination corrected*',/fold) eq 0B),nretain)

       if (nretain eq 0L) then begin
          splog, 'Need to handle this case.'
          stop
       endif else begin
          istart = istart[retain]
          iend = iend[retain]
       endelse
       
; parse data on each individual reference and store the results        

       refinfo = replicate(refinfo_template,nretain)
       for i = 0L, nretain-1L do begin

          subtext = text[istart[i]:iend[i]]

          line = where(strmatch(subtext,'*reference code*',/fold) eq 1B)
          i1 = strpos(subtext[line],'refcode=')+8L
          refinfo[i].reference = strtrim(strmid(subtext[line[0]],i1,19),2)

          line = where(strmatch(subtext,'*frequency targeted*',/fold) eq 1B)
          refinfo[i].frequency_targeted = strtrim(strmid(subtext[line[0]],strpos(subtext[line],':')+1L),2)

          line = where(strmatch(subtext,'*major axis*',/fold) eq 1B)
          subline = strmid(subtext[line[0]],68)
          if (strcompress(subline,/remove) ne '') and (strmatch(subline,'*...*') eq 0B) then begin
             i1 = strpos(subline,'2a =')+4L & i2 = strpos(subline,'arcsec')
             refinfo[i].major_axis       = strmid(subline,i1,i2-i1)
             refinfo[i].major_axis_units = strtrim(strmid(subline,i2),2)
          endif
             
          line = where(strmatch(subtext,'*minor axis*',/fold) eq 1B,nline)
          if (nline ne 0L) then begin
             subline = strmid(subtext[line[0]],68)
             if (strcompress(subline,/remove) ne '') and (strmatch(subline,'*...*') eq 0B) then begin
                i1 = strpos(subline,'2b =')+4L & i2 = strpos(subline,'arcsec')
                refinfo[i].minor_axis       = strmid(subline,i1,i2-i1)
                refinfo[i].minor_axis_units = strtrim(strmid(subline,i2),2)
             endif
          endif
             
          line = where(strmatch(subtext,'*axis ratio*',/fold) eq 1B,nline)
          if (nline ne 0L) then begin
             subline = strmid(subtext[line[0]],68)
             if (strcompress(subline,/remove) ne '') and (strmatch(subline,'*...*') eq 0B) then begin
                i1 = strpos(subline,'b/a =')+5L
                refinfo[i].axis_ratio = strmid(subline,i1)
             endif
          endif
             
          if (refinfo[i].major_axis gt -900.0) and (refinfo[i].minor_axis lt -900.0) and $
            (refinfo[i].axis_ratio gt -900.0) then begin
             refinfo[i].minor_axis       = refinfo[i].major_axis*refinfo[i].axis_ratio
             refinfo[i].minor_axis_units = refinfo[i].major_axis_units
          endif
             
          line = where(strmatch(subtext,'*position angle*',/fold) eq 1B)
          subline = strmid(subtext[line[0]],68)
          if (strcompress(subline,/remove) ne '') and (strmatch(subline,'*...*') eq 0B) then begin
             i1 = 68L & i2 = strpos(subline,'deg')
             refinfo[i].posangle = strmid(subline,0L,i2)
             if (strmatch(subline,'*J2000*') or strmatch(subline,'*B1950*')) then $
               refinfo[i].posangle_equinox = strmid(subline,strpos(subline,'[')+1,5)
          endif
          
          line = where(strmatch(subtext,'*Diameter reference level*',/fold) eq 1B)
          refinfo[i].diameter_reflevel = strtrim(strmid(subtext[line[0]],strpos(subtext[line],':')+1L),2)

          line = where(strmatch(subtext,'*Frequency  *',/fold) eq 1B)
          i1 = strpos(subtext[line],':')+1L & i2 = 67L ; <-- NOTE!
          refinfo[i].diameter_frequency = strtrim(strmid(subtext[line[0]],i1,i2-i1),2)

          line = where(strmatch(subtext,'*detector type*',/fold) eq 1B)
          refinfo[i].detector_type = strtrim(strmid(subtext[line[0]],strpos(subtext[line],':')+1L),2)

          line = where(strmatch(subtext,'*fitting technique*',/fold) eq 1B)
          refinfo[i].fitting_technique = strtrim(strmid(subtext[line[0]],strpos(subtext[line],':')+1L),2)

          line = where(strmatch(subtext,'*measured quantity*',/fold) eq 1B)
          refinfo[i].measured_quantity = strtrim(strmid(subtext[line[0]],strpos(subtext[line],':')+1L),2)

          line = where(strmatch(subtext,'*qualifiers*',/fold) eq 1B,nline)
          i1 = strpos(subtext[line[0]],':')+1L & i2 = strlen(subtext[line[0]])

          if (nline eq 1L) then refinfo[i].data_qualifiers = strtrim(strmid(subtext[line[0]],i1,i2-i1),2)

          if (nline eq 2L) then begin
             refinfo[i].measurement_qualifiers = strtrim(strmid(subtext[line[0]],i1,i2-i1),2)
             i1 = strpos(subtext[line[1]],':')+1L & i2 = strlen(subtext[line[1]])
             refinfo[i].data_qualifiers = strtrim(strmid(subtext[line[1]],i1,i2-i1),2)
          endif

          line = where(strmatch(subtext,'*ned comment*',/fold) eq 1B)
          i1 = strpos(subtext[line],':')+1L & i2 = strlen(subtext[line])
          refinfo[i].ned_comment = strtrim(strmid(subtext[line[0]],i1,i2-i1),2)

       endfor

; now retrieve only the references of interest
       
       rc3 = where($
         (strmatch(refinfo.reference,'*1991RC3.9.C...0000d*') eq 1B) and $
         (strmatch(refinfo.measured_quantity,'*isophotal diameter*',/fold) eq 1B),nrc3)

       if (nrc3 ne 0L) then begin

          if (nrc3 eq 1L) then begin

             diameters[igalaxy].rc3_reference          = refinfo[rc3].reference
             diameters[igalaxy].rc3_major_axis         = refinfo[rc3].major_axis
             diameters[igalaxy].rc3_minor_axis         = refinfo[rc3].minor_axis
             diameters[igalaxy].rc3_posangle           = refinfo[rc3].posangle
             diameters[igalaxy].rc3_posangle_equinox   = refinfo[rc3].posangle_equinox
             diameters[igalaxy].rc3_diameter_reflevel  = refinfo[rc3].diameter_reflevel
             diameters[igalaxy].rc3_diameter_frequency = refinfo[rc3].diameter_frequency
             diameters[igalaxy].rc3_measured_quantity  = refinfo[rc3].measured_quantity
             diameters[igalaxy].rc3_ned_comment        = refinfo[rc3].ned_comment

          endif else begin

             splog, 'Multiple RC3 diameters for '+galaxy[igalaxy]+'.'

          endelse
             
       endif

; ESO-Uppsala "Quick Blue"
       
       eso = where($
         (strmatch(refinfo.reference,'*1982ESOU..C...0000L*') eq 1B) or $
         (strmatch(refinfo.reference,'*1989ESOLV.C...0000L*') eq 1B),neso)

       if (neso ne 0L) then begin

          eso89 = where((strmatch(refinfo[eso].reference,'*1989ESOLV.C...0000L*') eq 1B) and $
            (strmatch(refinfo[eso].diameter_reflevel,'*25.0*') eq 1B),neso89)
          
          if (neso89 ne 0L) then begin
             if (neso89 eq 1L) then begin
                diameters[igalaxy].eso_reference          = refinfo[eso[eso89]].reference
                diameters[igalaxy].eso_major_axis         = refinfo[eso[eso89]].major_axis
                diameters[igalaxy].eso_minor_axis         = refinfo[eso[eso89]].minor_axis
                diameters[igalaxy].eso_posangle           = refinfo[eso[eso89]].posangle
                diameters[igalaxy].eso_posangle_equinox   = refinfo[eso[eso89]].posangle_equinox
                diameters[igalaxy].eso_diameter_reflevel  = refinfo[eso[eso89]].diameter_reflevel
                diameters[igalaxy].eso_diameter_frequency = refinfo[eso[eso89]].diameter_frequency
                diameters[igalaxy].eso_measured_quantity  = refinfo[eso[eso89]].measured_quantity
                diameters[igalaxy].eso_ned_comment        = refinfo[eso[eso89]].ned_comment
             endif else begin
                splog, 'Multiple ESO/1989 25.0 SB diameters for '+galaxy[igalaxy]+'.'
             endelse
          endif

          eso82 = where((strmatch(refinfo.reference,'*1982ESOU..C...0000L*') eq 1B),neso82)

          if (neso82 ne 0L) and (neso89 eq 0L) then begin
             if (neso82 eq 1L) then begin
                diameters[igalaxy].eso_reference          = refinfo[eso[eso82]].reference
                diameters[igalaxy].eso_major_axis         = refinfo[eso[eso82]].major_axis
                diameters[igalaxy].eso_minor_axis         = refinfo[eso[eso82]].minor_axis
                diameters[igalaxy].eso_posangle           = refinfo[eso[eso82]].posangle
                diameters[igalaxy].eso_posangle_equinox   = refinfo[eso[eso82]].posangle_equinox
                diameters[igalaxy].eso_diameter_reflevel  = refinfo[eso[eso82]].diameter_reflevel
                diameters[igalaxy].eso_diameter_frequency = refinfo[eso[eso82]].diameter_frequency
                diameters[igalaxy].eso_measured_quantity  = refinfo[eso[eso82]].measured_quantity
                diameters[igalaxy].eso_ned_comment        = refinfo[eso[eso82]].ned_comment
             endif else begin
                splog, 'Multiple ESO/1982 diameters for '+galaxy[igalaxy]+'.'
             endelse
          endif 
             
       endif 

; MCG

       mcgall = where($
         (strmatch(refinfo.reference,'*1962MCG1..C...0000V*') eq 1B) or $
         (strmatch(refinfo.reference,'*1964MCG2..C...0000V*') eq 1B) or $
         (strmatch(refinfo.reference,'*1963MCG3..C...0000V*') eq 1B) or $
         (strmatch(refinfo.reference,'*1968MCG4..C...0000V*') eq 1B),nmcgall)
       
       if (nmcgall ne 0L) then begin

          mcg1 = where(strmatch(refinfo[mcgall].reference,'*1962MCG1..C...0000V*') eq 1B,nmcg1)
          mcg2 = where(strmatch(refinfo[mcgall].reference,'*1964MCG2..C...0000V*') eq 1B,nmcg2)
          mcg3 = where(strmatch(refinfo[mcgall].reference,'*1963MCG3..C...0000V*') eq 1B,nmcg3)
          mcg4 = where(strmatch(refinfo[mcgall].reference,'*1968MCG4..C...0000V*') eq 1B,nmcg4)

          if (nmcg1 ne 0L) and (nmcg2 eq 0L) and (nmcg3 eq 0L) and (nmcg4 eq 0L) then mcgindx = mcgall[mcg1]
          if (nmcg1 eq 0L) and (nmcg2 ne 0L) and (nmcg3 eq 0L) and (nmcg4 eq 0L) then mcgindx = mcgall[mcg2]
          if (nmcg1 eq 0L) and (nmcg2 eq 0L) and (nmcg3 ne 0L) and (nmcg4 eq 0L) then mcgindx = mcgall[mcg3]
          if (nmcg1 eq 0L) and (nmcg2 eq 0L) and (nmcg3 eq 0L) and (nmcg4 ne 0L) then mcgindx = mcgall[mcg4]

          nmcg = n_elements(mcgindx)
          
          if (nmcg gt 1L) then begin
             outer = where(strmatch(refinfo[mcgindx].measured_quantity,'*"outer" diameter*',/fold) eq 1B,nouter)
             if (nouter gt 1L) then begin
                splog, 'Multiple outer MCG diameters for '+galaxy[igalaxy]+'.'
                outer = outer[0] & nouter = 1L
             endif
             if (nouter ne 0L) then mcgindx = mcgindx[outer]
             inner = where(strmatch(refinfo[mcgindx].measured_quantity,'*"inner" diameter*',/fold) eq 1B,ninner)
             if (ninner ne 0L) then mcgindx = mcgindx[inner]
          endif

          diameters[igalaxy].mcg_reference          = refinfo[mcgindx].reference
          diameters[igalaxy].mcg_major_axis         = refinfo[mcgindx].major_axis
          diameters[igalaxy].mcg_minor_axis         = refinfo[mcgindx].minor_axis
          diameters[igalaxy].mcg_posangle           = refinfo[mcgindx].posangle
          diameters[igalaxy].mcg_posangle_equinox   = refinfo[mcgindx].posangle_equinox
          diameters[igalaxy].mcg_diameter_reflevel  = refinfo[mcgindx].diameter_reflevel
          diameters[igalaxy].mcg_diameter_frequency = refinfo[mcgindx].diameter_frequency
          diameters[igalaxy].mcg_measured_quantity  = refinfo[mcgindx].measured_quantity
          diameters[igalaxy].mcg_ned_comment        = refinfo[mcgindx].ned_comment
             
       endif else nmcg = 0L

; UGC

       ugc = where((strmatch(refinfo.reference,'*1973UGC...C...0000N*') eq 1B),nugc)
         
       if (nugc ne 0L) then begin

          if (nugc gt 1L) then begin
             outer = where(strmatch(refinfo[ugc].measured_quantity,'*"outer" diameter*',/fold) eq 1B,nouter)
             if (nouter ne 0L) then ugc = ugc[outer]
             inner = where(strmatch(refinfo[ugc].measured_quantity,'*"inner" diameter*',/fold) eq 1B,ninner)
             if (ninner ne 0L) then ugc = ugc[inner]
          endif

          diameters[igalaxy].ugc_reference          = refinfo[ugc].reference
          diameters[igalaxy].ugc_major_axis         = refinfo[ugc].major_axis
          diameters[igalaxy].ugc_minor_axis         = refinfo[ugc].minor_axis
          diameters[igalaxy].ugc_posangle           = refinfo[ugc].posangle
          diameters[igalaxy].ugc_posangle_equinox   = refinfo[ugc].posangle_equinox
          diameters[igalaxy].ugc_diameter_reflevel  = refinfo[ugc].diameter_reflevel
          diameters[igalaxy].ugc_diameter_frequency = refinfo[ugc].diameter_frequency
          diameters[igalaxy].ugc_measured_quantity  = refinfo[ugc].measured_quantity
          diameters[igalaxy].ugc_ned_comment        = refinfo[ugc].ned_comment
             
       endif

; 2MASS       
       
       twomass = where($
         ((strmatch(refinfo.reference,'*2003AJ....125..525J*') eq 1B) or $
         (strmatch(refinfo.reference,'*20032MASX.C.......*') eq 1B)) and $
         (strmatch(refinfo.measured_quantity,'*isophotal diameter*',/fold) eq 1B),ntwomass)
       
       if (ntwomass ne 0L) then begin

          this_tot = where((strmatch(refinfo[twomass].frequency_targeted,'*total*',/fold) eq 1B),nthis_tot)
          this_iso = where((strmatch(refinfo[twomass].frequency_targeted,'*isophotal*',/fold) eq 1B),nthis_iso)
          if (nthis_tot ne nthis_iso) then message, 'Need to handle this case.'

          if (nthis_tot ne 0L) then begin
             if (nthis_tot eq 1L) then begin
                diameters[igalaxy].twomass_reference             = refinfo[twomass[this_tot]].reference
                diameters[igalaxy].twomass_major_axis            = refinfo[twomass[this_tot]].major_axis
                diameters[igalaxy].twomass_minor_axis            = refinfo[twomass[this_tot]].minor_axis
                diameters[igalaxy].twomass_k20_major_axis        = refinfo[twomass[this_iso]].major_axis
                diameters[igalaxy].twomass_k20_minor_axis        = refinfo[twomass[this_iso]].minor_axis
                diameters[igalaxy].twomass_posangle              = refinfo[twomass[this_tot]].posangle
                diameters[igalaxy].twomass_posangle_equinox      = refinfo[twomass[this_tot]].posangle_equinox
                diameters[igalaxy].twomass_diameter_reflevel     = refinfo[twomass[this_tot]].diameter_reflevel
                diameters[igalaxy].twomass_k20_diameter_reflevel = refinfo[twomass[this_iso]].diameter_reflevel
                diameters[igalaxy].twomass_diameter_frequency    = refinfo[twomass[this_tot]].diameter_frequency
                diameters[igalaxy].twomass_measured_quantity     = refinfo[twomass[this_tot]].measured_quantity
                diameters[igalaxy].twomass_ned_comment           = refinfo[twomass[this_tot]].ned_comment
             endif else begin
                splog, 'Multiple 2MASS diameters for '+galaxy[igalaxy]+'.'
             endelse
          endif

       endif 

       if (nrc3 eq 0L) and (neso eq 0L) and (nmcg eq 0L) and (nugc eq 0L) and $
         (ntwomass eq 0L) then splog, 'No preferred diameter data for '+galaxy[igalaxy]+'.'
       
    endfor

; write out    
    
    if keyword_set(write) then begin

       splog, 'Writing '+outfile
       mwrfits, diameters, outfile, /create
       spawn, ['gzip -f '+outfile], /sh

    endif
    
return
end
    
