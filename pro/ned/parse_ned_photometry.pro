;+
; NAME:
;       PARSE_NED_PHOTOMETRY
;
; PURPOSE:
;       Parse galaxy photometry from a NED batch request.
;
; INPUTS:
;       nedfile      - photometry text file from NED
;
; OPTIONAL INPUTS:
;       inputnedfile - input text file sent to NED with a list of
;                      objects, one per line (spaces and funny
;                      characters are okay), for which to get
;                      photometry
;       nedpath      - path name to NEDFILE and INPUTNEDFILE
;       outfile      - output file name (either OUTFILE.FITS or
;                      OUTFILE.TXT, depending on whether TEXTFILE=1)
;       outpath      - output path name for OUTFILE
;
; KEYWORD PARAMETERS:
;       textfile     - write a columnular text file rather than a
;                      binary FITS file
;
; OUTPUTS:
;       data         - output photometry data structure
;          GALAXY        - input galaxy name from INPUTNEDFILE
;          NEDGALAXY     - NED galaxy name from NEDFILE
;          U_RC3         - U-band RC3 magnitude
;          U_RC3_ERR     - U-band RC3 magnitude error
;          U_FLAG        - U-band magnitude flag (U_T or U_T^0)
;          B_RC3         - B-band RC3 magnitude
;          B_RC3_ERR     - B-band RC3 magnitude error
;          B_FLAG        - B-band magnitude flag (B_T, m_B, or m_B^0)
;          V_RC3         - V-band RC3 magnitude
;          V_RC3_ERR     - V-band RC3 magnitude error
;          V_FLAG        - V-band magnitude flag (V_T or V_T^0)
;          J_2MASS       - J-band 2MASS magnitude
;          J_2MASS_ERR   - J-band 2MASS magnitude error
;          J_2MASS_FLAG  - J-band 2MASS flag [not used]
;          H_2MASS       - H-band 2MASS magnitude
;          H_2MASS_ERR   - H-band 2MASS magnitude error
;          H_2MASS_FLAG  - H-band 2MASS flag [not used]
;          K_2MASS       - K-band 2MASS magnitude
;          K_2MASS_ERR   - K-band 2MASS magnitude error
;          K_2MASS_FLAG  - K-band 2MASS flag [not used]
;          IRAS_12       - IRAS 12 micron flux density [Jy]
;          IRAS_12_ERR   - IRAS 12 micron flux density error [Jy]
;          IRAS_12_FLAG  - IRAS 12 micron reference
;          IRAS_25       - IRAS 25 micron flux density [Jy]
;          IRAS_25_ERR   - IRAS 25 micron flux density error [Jy]
;          IRAS_25_FLAG  - IRAS 25 micron reference
;          IRAS_60       - IRAS 60 micron flux density [Jy]
;          IRAS_60_ERR   - IRAS 60 micron flux density error [Jy]
;          IRAS_60_FLAG  - IRAS 60 micron reference
;          IRAS_100      - IRAS 100 micron flux density [Jy]
;          IRAS_100_ERR  - IRAS 100 micron flux density error [Jy]
;          IRAS_100_FLAG - IRAS 100 micron reference
;          UV_1650       - UV(1650) flux [erg/s/cm2/A]
;          UV_1650_ERR   - UV(1650) flux error [erg/s/cm2/A]
;          UV_1650_REF   - UV(1650) reference
;          UV_2500       - UV(2500) flux [erg/s/cm2/A]
;          UV_2500_ERR   - UV(2500) flux error [erg/s/cm2/A]
;          UV_2500_REF   - UV(2500) reference
;          UV_3150       - UV(3150) flux [erg/s/cm2/A]
;          UV_3150_ERR   - UV(3150) flux error [erg/s/cm2/A]
;          UV_3150_REF   - UV(3150) reference
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;       Upper limits on IRAS fluxes are negative.  The UV data are
;       from Rifatto, Longo, & Capaciol 1995, A&AS, 114, 527 (see also
;       Rifatto et al. 1995, A&AS, 109, 341).  We exclude upper limits
;       on the UV fluxes.
; 
;       Sometimes the 2MASS data appears twice.  According to NED, the
;       brighter of the two should be used:
;       
;          Date: Mon, 29 Nov 2004 15:38:11 -0800 (PST)
;          From: Harold Corwin <hgcjr@ipac.caltech.edu>
;          To: jmoustakas@as.arizona.edu
;          Cc: jarrett@ipac.caltech.edu, zb4ms@ipac.caltech.edu
;          Subject: 2MASS photometry
;          
;          Hi, John,
;          
;           > I am using the NED batch photometry retrieval system to compile 2MASS
;           > photometry for my galaxy sample.  I have discovered that several objects
;           > in my sample have *two* measurements corresponding to the "total"  
;           > magnitudes in JHKs from the Extended Source Catalog.  These magnitudes are
;           > similar in some objects but very different in others.  Below is the list
;           > of objects in my sample with the two measurements.  If you can help me
;           > resolve which measurement I should be using I would appreciate the help.
;          
;          In all of these cases, there are two 2MASS XSC names attached to 
;          the NED object; these "interloper" XSC names pulled in the 
;          additional magnitudes your NED searches returned.
;          
;          I've checked the images for each of these and have found that the 
;          brighter of the two sets of magnitudes is the appropriate one to 
;          use in all of these cases (some of these, however, will need to 
;          be revised by Tom Jarrett using his Large Galaxy Atlas software; 
;          NGC 4138 in particular is affected).  I've also appended a note 
;          briefly explaining the source of each problem.
;
;       NED provides 2MASS photometry for some objects from the Large
;       Galaxy Atlas (LGA).  This catalog includes total photometry
;       that is "better" than the XSC.  It also provides isophotal
;       photometry (at 20 mag/arcsec2) in each of the three bands
;       (JHKs), whereas the XSC provides isophotal photometry at the
;       K20 radius.  For uniformity we only take the XSC isophotal
;       photometry.
;
;       IRAS fluxes exist from three sources: Rice et al. 1988, ApJS,
;       68, 91 for large optical galaxies; Soifer et al. 1989, AJ, 98,
;       766 for the Bright Galaxy Sample; and Moshir et al. 1990 for
;       the IRAS Faint Source Catalog.  All three of these fluxes are
;       stored.  Problems are sometimes encountered whereby the fluxes
;       are very different, or one source has a flux, but no flux
;       error, etc..  
;
; EXAMPLE:
;       IDL> parse_ned_photometry, 'ned.phot.dat', outfile='ned.phot.fits'
;       IDL> parse_ned_photometry, 'ned.phot.dat', /textfile, outfile='ned.phot.txt'
;
; MODIFICATION HISTORY:
;       J. Moustakas, 2002 May 14, U of A
;       jm03feb10uofa - documented and streamlined
;       jm03mar15uofa - grabbed IRAS errors and distinguished upper
;                       limits; more error checking
;       jm03apr30uofa - parse 2MASS photometry newly added to NED 
;       jm03dec29uofa - added error checking for multiple 2MASS
;                       measurements 
;       jm04nov29uofa - smarter reading of the 2MASS magnitudes;
;                       account for total magnitudes from the Large
;                       Galaxy Atlas; also read both the 20th
;                       isophotal radius magnitudes and the "total"
;                       magnitudes 
;       jm04dec07uofa - generalized the IRAS photometry to return all
;                       three of the RICE, SOIFER, and MOSHIR
;                       photometry 
;       jm07jun25nyu  - parse GALEX and MIPS photometry
;
; Copyright (C) 2002-2004, 2007, John Moustakas
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

pro parse_ned_photometry, nedfile, data, inputnedfile=inputnedfile, $
  nedpath=nedpath, textfile=textfile, outfile=outfile, outpath=outpath
    
    if (n_elements(nedfile) eq 0L) then begin
       doc_library, 'parse_ned_photometry'
       return
    endif

; read the NED photometry file

    if n_elements(nedpath) eq 0L then nedpath = cwd()
    
    if file_test(nedpath+nedfile,/regular) eq 0L then begin
       splog, 'File '+nedpath+nedfile+' not found.'
       return
    endif
    
    neddata = djs_readlines(nedpath+nedfile)
    lstart = (where(strmatch(neddata,'*SEARCH RESULTS*') eq 1B))[0]
    neddata = neddata[lstart:n_elements(neddata)-1L]

    objstart = where(strmatch(neddata,'*PHOTO *') eq 1B,nobj)
    objstart = [objstart,where(strmatch(neddata,'*ABOUT THE DATA*') eq 1B)]
    splog, 'There are '+strn(nobj)+' objects in the NED photometry file.'

; read the input text file to NED, if it was given

    if n_elements(inputnedfile) ne 0L then begin
       
       if file_test(nedpath+inputnedfile,/regular) eq 0L then begin
          splog, 'Input ned file '+nedpath+inputnedfile+' not found.'
          return
       endif

       inputfile = strcompress(djs_readlines(nedpath+inputnedfile),/remove)
       galaxy = inputfile[where(inputfile ne '')]
       ngalaxy = n_elements(galaxy)

       splog, 'There are '+strn(ngalaxy)+' objects in '+inputnedfile+'.'

       if (nobj ne ngalaxy) then begin
          splog, 'NED did not return photometry for all the objects . . . returning.'
          return
       endif

    endif
    
; initialize the output data structure

    data = {$
      galaxy:               '...',  $
      nedgalaxy:            '...',  $
      U_RC3:                -999.0, $
      U_RC3_ERR:            -999.0, $
      U_FLAG:               '...',  $
      B_RC3:                -999.0, $
      B_RC3_ERR:            -999.0, $
      B_FLAG:               '...',  $
      V_RC3:                -999.0, $
      V_RC3_ERR:            -999.0, $
      V_FLAG:               '...',  $
      J_2MASS:              -999.0, $
      J_2MASS_ERR:          -999.0, $
      J20_2MASS:            -999.0, $
      J20_2MASS_ERR:        -999.0, $
      J_2MASS_FLAG:         '...',  $
      H_2MASS:              -999.0, $
      H_2MASS_ERR:          -999.0, $
      H20_2MASS:            -999.0, $
      H20_2MASS_ERR:        -999.0, $
      H_2MASS_FLAG:         '...',  $
      K_2MASS:              -999.0, $
      K_2MASS_ERR:          -999.0, $
      K20_2MASS:            -999.0, $
      K20_2MASS_ERR:        -999.0, $
      K_2MASS_FLAG:         '...',  $
      SANDERS_IRAS_12:      -999.0, $
      SANDERS_IRAS_12_ERR:  -999.0, $
      RICE_IRAS_12:         -999.0, $
      RICE_IRAS_12_ERR:     -999.0, $
      SOIFER_IRAS_12:       -999.0, $
      SOIFER_IRAS_12_ERR:   -999.0, $
      MOSHIR_IRAS_12:       -999.0, $
      MOSHIR_IRAS_12_ERR:   -999.0, $
      SANDERS_IRAS_25:      -999.0, $
      SANDERS_IRAS_25_ERR:  -999.0, $
      RICE_IRAS_25:         -999.0, $
      RICE_IRAS_25_ERR:     -999.0, $
      SOIFER_IRAS_25:       -999.0, $
      SOIFER_IRAS_25_ERR:   -999.0, $
      MOSHIR_IRAS_25:       -999.0, $
      MOSHIR_IRAS_25_ERR:   -999.0, $
      SANDERS_IRAS_60:      -999.0, $
      SANDERS_IRAS_60_ERR:  -999.0, $
      RICE_IRAS_60:         -999.0, $
      RICE_IRAS_60_ERR:     -999.0, $
      SOIFER_IRAS_60:       -999.0, $
      SOIFER_IRAS_60_ERR:   -999.0, $
      MOSHIR_IRAS_60:       -999.0, $
      MOSHIR_IRAS_60_ERR:   -999.0, $
      SANDERS_IRAS_100:     -999.0, $
      SANDERS_IRAS_100_ERR: -999.0, $
      RICE_IRAS_100:        -999.0, $
      RICE_IRAS_100_ERR:    -999.0, $
      SOIFER_IRAS_100:      -999.0, $
      SOIFER_IRAS_100_ERR:  -999.0, $
      MOSHIR_IRAS_100:      -999.0, $
      MOSHIR_IRAS_100_ERR:  -999.0, $
      UV_1650:              -999.0, $
      UV_1650_ERR:          -999.0, $
      UV_1650_REF:           '...', $
      UV_2500:              -999.0, $
      UV_2500_ERR:          -999.0, $
      UV_2500_REF:           '...', $
      UV_3150:              -999.0, $
      UV_3150_ERR:          -999.0, $
      UV_3150_REF:          '...' , $
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
    data = replicate(data,nobj)
    
; go through each object and pull out the photometry

    for i = 0L, nobj-1L do begin

       subtext = neddata[objstart[i]:objstart[i+1L]-1L]

       if n_elements(ngalaxy) ne 0L then data[i].galaxy = strcompress(galaxy[i],/remove)
       data[i].nedgalaxy = strcompress(strmid(subtext[0],6,strlen(subtext[0])),/remove)

; MIPS references to think about: Dale et al. 2007, UV-IR atlas
; (2007ApJ...655..863D) (IRAC+MIPS), and Smith et al. 2007, Spiral,
; Bridges, and Tails Spitzer Survey (2007AJ....133..791S) (only
; IRAC+MIPS/24)
       
; MIPS 24/70/160

       mips24 = where((strmatch(subtext,'*24*(MIPS)*') eq 1B) and $
         (strmatch(subtext,'*2007ApJ...655..863D*') eq 0B),nmips24)
       if nmips24 ne 0L then begin
;         data[i].mips_24 = 
       endif

;; --------------------------------------------------
;; 24 microns
;; --------------------------------------------------
;
;       mips24 = where((strmatch(subtext,'*24*(MIPS)*') eq 1B) and $
;         (strmatch(subtext,'*1991RC3.9.C...0000d*') eq 1B),nmips24)
;       if nmips24 ne 0L then data[i].U_flag = 'U_T'
;
;       if nmips24 eq 0L then begin
;          mips24 = where((strmatch(subtext,'*U (U_T^0)*') eq 1B) and $
;            (strmatch(subtext,'*1991RC3.9.C...0000d*') eq 1B),nmips24)
;          if nmips24 ne 0L then data[i].U_flag = 'U_T^0'
;       endif
;
;       if nmips24 ne 0L then begin
;          data[i].U_RC3 = strmid(subtext[mips24],27,11)
;          err = strmid(subtext[mips24],41,6)
;          if strcompress(err,/remove) ne '' then data[i].U_RC3_ERR = err
;       endif
       
       
; RC3 catalog [mag] - the B magnitudes are prioritized as below
       
; --------------------------------------------------
; U-band
; --------------------------------------------------
       
       urc3 = where((strmatch(subtext,'*U (U_T)*') eq 1B) and $
         (strmatch(subtext,'*1991RC3.9.C...0000d*') eq 1B),nurc3)
       if nurc3 ne 0L then data[i].U_flag = 'U_T'

       if nurc3 eq 0L then begin
          urc3 = where((strmatch(subtext,'*U (U_T^0)*') eq 1B) and $
            (strmatch(subtext,'*1991RC3.9.C...0000d*') eq 1B),nurc3)
          if nurc3 ne 0L then data[i].U_flag = 'U_T^0'
       endif

       if nurc3 ne 0L then begin
          data[i].U_RC3 = strmid(subtext[urc3],27,11)
          err = strmid(subtext[urc3],41,6)
          if strcompress(err,/remove) ne '' then data[i].U_RC3_ERR = err
       endif
       
; --------------------------------------------------
; B-band
; --------------------------------------------------
       
       brc3 = where((strmatch(subtext,'*B (B_T)*') eq 1B) and $
         (strmatch(subtext,'*1991RC3.9.C...0000d*') eq 1B),nbrc3)
       if nbrc3 ne 0L then data[i].B_flag = 'B_T'
       
       if nbrc3 eq 0L then begin
          brc3 = where((strmatch(subtext,'*B (m_B)*') eq 1B) and $
            (strmatch(subtext,'*1991RC3.9.C...0000d*') eq 1B),nbrc3)
          if nbrc3 ne 0L then data[i].B_flag = 'm_B'
       endif

       if nbrc3 eq 0L then begin
          brc3 = where((strmatch(subtext,'*B (B_T^0)*') eq 1B) and $
            (strmatch(subtext,'*1991RC3.9.C...0000d*') eq 1B),nbrc3)
          if nbrc3 ne 0L then data[i].B_flag = 'B_T^0'
       endif

       if nbrc3 eq 0L then begin
          brc3 = where((strmatch(subtext,'*B (m_B^0)*') eq 1B) and $
            (strmatch(subtext,'*1991RC3.9.C...0000d*') eq 1B),nbrc3)
          if nbrc3 ne 0L then data[i].B_flag = 'm_B^0'
       endif

       if nbrc3 ne 0L then begin
          data[i].B_RC3 = strmid(subtext[brc3],27,11)
          err = strmid(subtext[brc3],41,6)
          if strcompress(err,/remove) ne '' then data[i].B_RC3_ERR = err
       endif

; --------------------------------------------------
; V-band
; --------------------------------------------------
       
       vrc3 = where((strmatch(subtext,'*V (V_T)*') eq 1B) and $
         (strmatch(subtext,'*1991RC3.9.C...0000d*') eq 1B),nvrc3)
       if nvrc3 ne 0L then data[i].V_flag = 'V_T'

       if nvrc3 eq 0L then begin
          vrc3 = where((strmatch(subtext,'*V (V_T^0)*') eq 1B) and $
            (strmatch(subtext,'*1991RC3.9.C...0000d*') eq 1B),nvrc3)
          if nvrc3 ne 0L then data[i].V_flag = 'V_T^0'
       endif
       
       if nvrc3 ne 0L then begin
          data[i].V_RC3 = strmid(subtext[vrc3],27,11)
          err = strmid(subtext[vrc3],41,6)
          if strcompress(err,/remove) ne '' then data[i].V_RC3_ERR = err
       endif

; --------------------------------------------------
; 2MASS J-band
; --------------------------------------------------

; total magnitude; first check if this object is in the LGA
       
       J2mass = where($
         ((strmatch(subtext,'*J_tot (2MASS LGA)*') eq 1B) and $
         (strmatch(subtext,'*2003AJ....125..525J*') eq 1B)),nJ2mass)

       if (nJ2mass ne 0L) then begin

          data[i].J_2MASS_FLAG = 'LGA'

       endif else begin

          J2mass = where((strmatch(subtext,'*J_total (2MASS)*') eq 1B) and $
            (strmatch(subtext,'*20032MASX.C.......:*') eq 1B),nJ2mass)
          if (nJ2mass ne 0L) then data[i].J_2MASS_FLAG = 'XSC'

       endelse

       if (nJ2mass ne 0L) then begin

          mags = strmid(subtext[J2mass],27,11)
          mag = min(mags,minindx)

          if (nJ2mass gt 1L) then begin
             splog, 'Multiple 2MASS J-band photometry: '+data[i].galaxy+': '+strjoin(strtrim(mags,2),', ')
          endif
          
          data[i].J_2MASS = mag
          err = strmid(subtext[J2mass[minindx]],41,6)
          if strcompress(err,/remove) ne '' then data[i].J_2MASS_ERR = err

       endif
       
; isophotal magnitude; first check if this object is in the LGA

       J2mass = where($
         ((strmatch(subtext,'*J_20 (2MASS LGA)*') eq 1B) and $
         (strmatch(subtext,'*2003AJ....125..525J*') eq 1B)),nJ2mass)

       if (nJ2mass ne 0L) then begin

       endif else begin

          J2mass = where((strmatch(subtext,'*J_K20 (2MASS)*') eq 1B) and $
            (strmatch(subtext,'*20032MASX.C.......:*') eq 1B),nJ2mass)

       endelse

       if (nJ2mass ne 0L) then begin

          mags = strmid(subtext[J2mass],27,11)
          mag = min(mags,minindx)

          if (nJ2mass gt 1L) then begin
             splog, 'Multiple 2MASS J-band photometry: '+data[i].galaxy+': '+strjoin(strtrim(mags,2),', ')
          endif
          
          data[i].J20_2MASS = mag
          err = strmid(subtext[J2mass[minindx]],41,6)
          if strcompress(err,/remove) ne '' then data[i].J20_2MASS_ERR = err

       endif
       
; --------------------------------------------------
; 2MASS H-band
; --------------------------------------------------

; total magnitude; first check if this object is in the LGA
       
       H2mass = where($
         ((strmatch(subtext,'*H_tot (2MASS LGA)*') eq 1B) and $
         (strmatch(subtext,'*2003AJ....125..525J*') eq 1B)),nH2mass)

       if (nH2mass ne 0L) then begin

          data[i].H_2MASS_FLAG = 'LGA'

       endif else begin

          H2mass = where((strmatch(subtext,'*H_total (2MASS)*') eq 1B) and $
            (strmatch(subtext,'*20032MASX.C.......:*') eq 1B),nH2mass)
          if (nH2mass ne 0L) then data[i].H_2MASS_FLAG = 'XSC'

       endelse

       if (nH2mass ne 0L) then begin

          mags = strmid(subtext[H2mass],27,11)
          mag = min(mags,minindx)

          if (nH2mass gt 1L) then begin
             splog, 'Multiple 2MASS H-band photometry: '+data[i].galaxy+': '+strjoin(strtrim(mags,2),', ')
          endif
          
          data[i].H_2MASS = mag
          err = strmid(subtext[H2mass[minindx]],41,6)
          if strcompress(err,/remove) ne '' then data[i].H_2MASS_ERR = err

       endif

; isophotal magnitude; first check if this object is in the LGA
       
       H2mass = where($
         ((strmatch(subtext,'*H_20 (2MASS LGA)*') eq 1B) and $
         (strmatch(subtext,'*2003AJ....125..525J*') eq 1B)),nH2mass)

       if (nH2mass ne 0L) then begin

       endif else begin

          H2mass = where((strmatch(subtext,'*H_K20 (2MASS)*') eq 1B) and $
            (strmatch(subtext,'*20032MASX.C.......:*') eq 1B),nH2mass)

       endelse

       if (nH2mass ne 0L) then begin

          mags = strmid(subtext[H2mass],27,11)
          mag = min(mags,minindx)

          if (nH2mass gt 1L) then begin
             splog, 'Multiple 2MASS H-band photometry: '+data[i].galaxy+': '+strjoin(strtrim(mags,2),', ')

          endif
          
          data[i].H20_2MASS = mag
          err = strmid(subtext[H2mass[minindx]],41,6)
          if strcompress(err,/remove) ne '' then data[i].H20_2MASS_ERR = err

       endif

; --------------------------------------------------
; 2MASS K-band
; --------------------------------------------------

; total magnitude; first check if this object is in the LGA
       
       K2mass = where($
         ((strmatch(subtext,'*K_tot (2MASS LGA)*') eq 1B) and $
         (strmatch(subtext,'*2003AJ....125..525J*') eq 1B)),nK2mass)

       if (nK2mass ne 0L) then begin

          data[i].K_2MASS_FLAG = 'LGA'
          
       endif else begin

          K2mass = where((strmatch(subtext,'*K_s_total (2MASS)*') eq 1B) and $
            (strmatch(subtext,'*20032MASX.C.......:*') eq 1B),nK2mass)
          if (nK2mass ne 0L) then data[i].K_2MASS_FLAG = 'XSC'

       endelse

       if (nK2mass ne 0L) then begin

          mags = strmid(subtext[K2mass],27,11)
          mag = min(mags,minindx)

          if (nK2mass gt 1L) then begin
             splog, 'Multiple 2MASS K-band photometry: '+data[i].galaxy+': '+strjoin(strtrim(mags,2),', ')
          endif
          
          data[i].K_2MASS = mag
          err = strmid(subtext[K2mass[minindx]],41,6)
          if strcompress(err,/remove) ne '' then data[i].K_2MASS_ERR = err

       endif

; isophotal magnitude; first check if this object is in the LGA
       
       K2mass = where($
         ((strmatch(subtext,'*K_20 (2MASS LGA)*') eq 1B) and $
         (strmatch(subtext,'*2003AJ....125..525J*') eq 1B)),nK2mass)

       if (nK2mass ne 0L) then begin

       endif else begin

          K2mass = where((strmatch(subtext,'*K_s_20 (2MASS)*') eq 1B) and $
            (strmatch(subtext,'*20032MASX.C.......:*') eq 1B),nK2mass)

       endelse

       if (nK2mass ne 0L) then begin

          mags = strmid(subtext[K2mass],27,11)
          mag = min(mags,minindx)

          if (nK2mass gt 1L) then begin
             splog, 'Multiple 2MASS K-band photometry: '+data[i].galaxy+': '+strjoin(strtrim(mags,2),', ')
          endif
          
          data[i].K20_2MASS = mag
          err = strmid(subtext[K2mass[minindx]],41,6)
          if strcompress(err,/remove) ne '' then data[i].K20_2MASS_ERR = err

       endif

; --------------------------------------------------
; IRAS 12 micron       
; --------------------------------------------------

; ##################################################       
; SANDERS
; ##################################################       
       
       iras12 = where((strmatch(subtext,'*IRAS 12 microns*') eq 1B) and $
         (strmatch(subtext,'*2003AJ....126.1607S*') eq 1B),niras12)

       if (niras12 ne 0L) then begin

          flux = strmid(subtext[iras12],26,11)

          if (niras12 gt 1L) then begin
             splog, 'WARNING: '+data[i].galaxy+' has multiple IRAS 12-micron [SANDERS] measurements:', flux
             iras12 = iras12[0]
             flux = flux[0]
          endif
          
          if strcompress(flux,/remove) ne '' then begin 
             if strmatch(flux,'*<*') eq 1B then begin ; upper limit
                flux = -float(strmid(flux,strpos(flux,'<')+1,strlen(flux)))
             endif
             data[i].sanders_iras_12 = float(flux)
          endif

          err = strmid(subtext[iras12],41,7)
          if strcompress(err,/remove) ne '' then begin
             if strmatch(err,'*%*') eq 1B then err = (strmid(err,0,strpos(err,'%'))/100.0)*abs(data[i].sanders_iras_12)
             data[i].sanders_iras_12_err = float(err)
          endif

          units = strcompress(strmid(subtext[iras12],48,11),/remove)
          if strmatch(units,'Jy',/fold) ne 1B then splog, 'IRAS 12 UNITS WARNING: ', i, data[i].galaxy
          
       endif

; ##################################################       
; RICE
; ##################################################       
       
       iras12 = where((strmatch(subtext,'*IRAS 12 microns*') eq 1B) and $
         (strmatch(subtext,'*1988ApJS...68...91R*') eq 1B),niras12)

       if (niras12 ne 0L) then begin

          flux = strmid(subtext[iras12],26,11)

          if (niras12 gt 1L) then begin
             splog, 'WARNING: '+data[i].galaxy+' has multiple IRAS 12-micron [RICE] measurements:', flux
             iras12 = iras12[0]
             flux = flux[0]
          endif
          
          if strcompress(flux,/remove) ne '' then begin 
             if strmatch(flux,'*<*') eq 1B then begin ; upper limit
                flux = -float(strmid(flux,strpos(flux,'<')+1,strlen(flux)))
             endif
             data[i].rice_iras_12 = float(flux)
          endif

          err = strmid(subtext[iras12],41,7)
          if strcompress(err,/remove) ne '' then begin
             if strmatch(err,'*%*') eq 1B then err = (strmid(err,0,strpos(err,'%'))/100.0)*abs(data[i].rice_iras_12)
             data[i].rice_iras_12_err = float(err)
          endif

          units = strcompress(strmid(subtext[iras12],48,11),/remove)
          if strmatch(units,'Jy',/fold) ne 1B then splog, 'IRAS 12 UNITS WARNING: ', i, data[i].galaxy
          
       endif

; ##################################################       
; SOIFER
; ##################################################       
       
       iras12 = where((strmatch(subtext,'*IRAS 12 microns*') eq 1B) and $
         (strmatch(subtext,'*1989AJ.....98..766S*') eq 1B),niras12)

       if (niras12 ne 0L) then begin

          flux = strmid(subtext[iras12],26,11)

          if (niras12 gt 1L) then begin
             splog, 'WARNING: '+data[i].galaxy+' has multiple IRAS 12-micron [SOIFER] measurements:', flux
             iras12 = iras12[0]
             flux = flux[0]
          endif
          
          if strcompress(flux,/remove) ne '' then begin 
             if strmatch(flux,'*<*') eq 1B then begin ; upper limit
                flux = -float(strmid(flux,strpos(flux,'<')+1,strlen(flux)))
             endif
             data[i].soifer_iras_12 = float(flux)
          endif

          err = strmid(subtext[iras12],41,7)
          if strcompress(err,/remove) ne '' then begin
             if strmatch(err,'*%*') eq 1B then err = (strmid(err,0,strpos(err,'%'))/100.0)*abs(data[i].soifer_iras_12)
             data[i].soifer_iras_12_err = float(err)
          endif

          units = strcompress(strmid(subtext[iras12],48,11),/remove)
          if strmatch(units,'Jy',/fold) ne 1B then splog, 'IRAS 12 UNITS WARNING: ', i, data[i].galaxy
          
       endif

; ##################################################       
; MOSHIR
; ##################################################       
       
       iras12 = where((strmatch(subtext,'*IRAS 12 microns*') eq 1B) and $
         (strmatch(subtext,'*1990IRASF.C...0000M*') eq 1B),niras12)

       if (niras12 ne 0L) then begin

          flux = strmid(subtext[iras12],26,11)

          if (niras12 gt 1L) then begin
             splog, 'WARNING: '+data[i].galaxy+' has multiple IRAS 12-micron [MOSHIR] measurements:', flux
             iras12 = iras12[0]
             flux = flux[0]
          endif
          
          if strcompress(flux,/remove) ne '' then begin 
             if strmatch(flux,'*<*') eq 1B then begin ; upper limit
                flux = -float(strmid(flux,strpos(flux,'<')+1,strlen(flux)))
             endif
             data[i].moshir_iras_12 = float(flux)
          endif

          err = strmid(subtext[iras12],41,7)
          if strcompress(err,/remove) ne '' then begin
             if strmatch(err,'*%*') eq 1B then err = (strmid(err,0,strpos(err,'%'))/100.0)*abs(data[i].moshir_iras_12)
             data[i].moshir_iras_12_err = float(err)
          endif

          units = strcompress(strmid(subtext[iras12],48,11),/remove)
          if strmatch(units,'Jy',/fold) ne 1B then splog, 'IRAS 12 UNITS WARNING: ', i, data[i].galaxy
          
       endif

; --------------------------------------------------
; IRAS 25 micron       
; --------------------------------------------------

; ##################################################       
; SANDERS
; ##################################################       
       
       iras25 = where((strmatch(subtext,'*IRAS 25 microns*') eq 1B) and $
         (strmatch(subtext,'*2003AJ....126.1607S*') eq 1B),niras25)

       if (niras25 ne 0L) then begin

          flux = strmid(subtext[iras25],26,11)

          if (niras25 gt 1L) then begin
             splog, 'WARNING: '+data[i].galaxy+' has multiple IRAS 25-micron [SANDERS] measurements:', flux
             iras25 = iras25[0]
             flux = flux[0]
          endif
          
          if strcompress(flux,/remove) ne '' then begin 
             if strmatch(flux,'*<*') eq 1B then begin ; upper limit
                flux = -float(strmid(flux,strpos(flux,'<')+1,strlen(flux)))
             endif
             data[i].sanders_iras_25 = float(flux)
          endif

          err = strmid(subtext[iras25],41,7)
          if strcompress(err,/remove) ne '' then begin
             if strmatch(err,'*%*') eq 1B then err = (strmid(err,0,strpos(err,'%'))/100.0)*abs(data[i].sanders_iras_25)
             data[i].sanders_iras_25_err = float(err)
          endif

          units = strcompress(strmid(subtext[iras25],48,11),/remove)
          if strmatch(units,'Jy',/fold) ne 1B then splog, 'IRAS 25 UNITS WARNING: ', i, data[i].galaxy
          
       endif

; ##################################################       
; RICE
; ##################################################       
       
       iras25 = where((strmatch(subtext,'*IRAS 25 microns*') eq 1B) and $
         (strmatch(subtext,'*1988ApJS...68...91R*') eq 1B),niras25)

       if (niras25 ne 0L) then begin

          flux = strmid(subtext[iras25],26,11)

          if (niras25 gt 1L) then begin
             splog, 'WARNING: '+data[i].galaxy+' has multiple IRAS 25-micron [RICE] measurements:', flux
             iras25 = iras25[0]
             flux = flux[0]
          endif
          
          if strcompress(flux,/remove) ne '' then begin 
             if strmatch(flux,'*<*') eq 1B then begin ; upper limit
                flux = -float(strmid(flux,strpos(flux,'<')+1,strlen(flux)))
             endif
             data[i].rice_iras_25 = float(flux)
          endif

          err = strmid(subtext[iras25],41,7)
          if strcompress(err,/remove) ne '' then begin
             if strmatch(err,'*%*') eq 1B then err = (strmid(err,0,strpos(err,'%'))/100.0)*abs(data[i].rice_iras_25)
             data[i].rice_iras_25_err = float(err)
          endif

          units = strcompress(strmid(subtext[iras25],48,11),/remove)
          if strmatch(units,'Jy',/fold) ne 1B then splog, 'IRAS 25 UNITS WARNING: ', i, data[i].galaxy
          
       endif

; ##################################################       
; SOIFER
; ##################################################       
       
       iras25 = where((strmatch(subtext,'*IRAS 25 microns*') eq 1B) and $
         (strmatch(subtext,'*1989AJ.....98..766S*') eq 1B),niras25)

       if (niras25 ne 0L) then begin

          flux = strmid(subtext[iras25],26,11)

          if (niras25 gt 1L) then begin
             splog, 'WARNING: '+data[i].galaxy+' has multiple IRAS 25-micron [SOIFER] measurements:', flux
             iras25 = iras25[0]
             flux = flux[0]
          endif
          
          if strcompress(flux,/remove) ne '' then begin 
             if strmatch(flux,'*<*') eq 1B then begin ; upper limit
                flux = -float(strmid(flux,strpos(flux,'<')+1,strlen(flux)))
             endif
             data[i].soifer_iras_25 = float(flux)
          endif

          err = strmid(subtext[iras25],41,7)
          if strcompress(err,/remove) ne '' then begin
             if strmatch(err,'*%*') eq 1B then err = (strmid(err,0,strpos(err,'%'))/100.0)*abs(data[i].soifer_iras_25)
             data[i].soifer_iras_25_err = float(err)
          endif

          units = strcompress(strmid(subtext[iras25],48,11),/remove)
          if strmatch(units,'Jy',/fold) ne 1B then splog, 'IRAS 25 UNITS WARNING: ', i, data[i].galaxy
          
       endif

; ##################################################       
; MOSHIR
; ##################################################       
       
       iras25 = where((strmatch(subtext,'*IRAS 25 microns*') eq 1B) and $
         (strmatch(subtext,'*1990IRASF.C...0000M*') eq 1B),niras25)

       if (niras25 ne 0L) then begin

          flux = strmid(subtext[iras25],26,11)

          if (niras25 gt 1L) then begin
             splog, 'WARNING: '+data[i].galaxy+' has multiple IRAS 25-micron [MOSHIR] measurements:', flux
             iras25 = iras25[0]
             flux = flux[0]
          endif
          
          if strcompress(flux,/remove) ne '' then begin 
             if strmatch(flux,'*<*') eq 1B then begin ; upper limit
                flux = -float(strmid(flux,strpos(flux,'<')+1,strlen(flux)))
             endif
             data[i].moshir_iras_25 = float(flux)
          endif

          err = strmid(subtext[iras25],41,7)
          if strcompress(err,/remove) ne '' then begin
             if strmatch(err,'*%*') eq 1B then err = (strmid(err,0,strpos(err,'%'))/100.0)*abs(data[i].moshir_iras_25)
             data[i].moshir_iras_25_err = float(err)
          endif

          units = strcompress(strmid(subtext[iras25],48,11),/remove)
          if strmatch(units,'Jy',/fold) ne 1B then splog, 'IRAS 25 UNITS WARNING: ', i, data[i].galaxy
          
       endif

; --------------------------------------------------
; IRAS 60 micron       
; --------------------------------------------------

; ##################################################       
; SANDERS
; ##################################################       
       
       iras60 = where((strmatch(subtext,'*IRAS 60 microns*') eq 1B) and $
         (strmatch(subtext,'*2003AJ....126.1607S*') eq 1B),niras60)

       if (niras60 ne 0L) then begin

          flux = strmid(subtext[iras60],26,11)

          if (niras60 gt 1L) then begin
             splog, 'WARNING: '+data[i].galaxy+' has multiple IRAS 60-micron [SANDERS] measurements:', flux
             iras60 = iras60[0]
             flux = flux[0]
          endif
          
          if strcompress(flux,/remove) ne '' then begin 
             if strmatch(flux,'*<*') eq 1B then begin ; upper limit
                flux = -float(strmid(flux,strpos(flux,'<')+1,strlen(flux)))
             endif
             data[i].sanders_iras_60 = float(flux)
          endif

          err = strmid(subtext[iras60],41,7)
          if strcompress(err,/remove) ne '' then begin
             if strmatch(err,'*%*') eq 1B then err = (strmid(err,0,strpos(err,'%'))/100.0)*abs(data[i].sanders_iras_60)
             data[i].sanders_iras_60_err = float(err)
          endif

          units = strcompress(strmid(subtext[iras60],48,11),/remove)
          if strmatch(units,'Jy',/fold) ne 1B then splog, 'IRAS 60 UNITS WARNING: ', i, data[i].galaxy
          
       endif

; ##################################################       
; RICE
; ##################################################       
       
       iras60 = where((strmatch(subtext,'*IRAS 60 microns*') eq 1B) and $
         (strmatch(subtext,'*1988ApJS...68...91R*') eq 1B),niras60)

       if (niras60 ne 0L) then begin

          flux = strmid(subtext[iras60],26,11)

          if (niras60 gt 1L) then begin
             splog, 'WARNING: '+data[i].galaxy+' has multiple IRAS 60-micron [RICE] measurements:', flux
             iras60 = iras60[0]
             flux = flux[0]
          endif
          
          if strcompress(flux,/remove) ne '' then begin 
             if strmatch(flux,'*<*') eq 1B then begin ; upper limit
                flux = -float(strmid(flux,strpos(flux,'<')+1,strlen(flux)))
             endif
             data[i].rice_iras_60 = float(flux)
          endif

          err = strmid(subtext[iras60],41,7)
          if strcompress(err,/remove) ne '' then begin
             if strmatch(err,'*%*') eq 1B then err = (strmid(err,0,strpos(err,'%'))/100.0)*abs(data[i].rice_iras_60)
             data[i].rice_iras_60_err = float(err)
          endif

          units = strcompress(strmid(subtext[iras60],48,11),/remove)
          if strmatch(units,'Jy',/fold) ne 1B then splog, 'IRAS 60 UNITS WARNING: ', i, data[i].galaxy
          
       endif

; ##################################################       
; SOIFER
; ##################################################       
       
       iras60 = where((strmatch(subtext,'*IRAS 60 microns*') eq 1B) and $
         (strmatch(subtext,'*1989AJ.....98..766S*') eq 1B),niras60)

       if (niras60 ne 0L) then begin

          flux = strmid(subtext[iras60],26,11)

          if (niras60 gt 1L) then begin
             splog, 'WARNING: '+data[i].galaxy+' has multiple IRAS 60-micron [SOIFER] measurements:', flux
             iras60 = iras60[0]
             flux = flux[0]
          endif
          
          if strcompress(flux,/remove) ne '' then begin 
             if strmatch(flux,'*<*') eq 1B then begin ; upper limit
                flux = -float(strmid(flux,strpos(flux,'<')+1,strlen(flux)))
             endif
             data[i].soifer_iras_60 = float(flux)
          endif

          err = strmid(subtext[iras60],41,7)
          if strcompress(err,/remove) ne '' then begin
             if strmatch(err,'*%*') eq 1B then err = (strmid(err,0,strpos(err,'%'))/100.0)*abs(data[i].soifer_iras_60)
             data[i].soifer_iras_60_err = float(err)
          endif

          units = strcompress(strmid(subtext[iras60],48,11),/remove)
          if strmatch(units,'Jy',/fold) ne 1B then splog, 'IRAS 60 UNITS WARNING: ', i, data[i].galaxy
          
       endif

; ##################################################       
; MOSHIR
; ##################################################       
       
       iras60 = where((strmatch(subtext,'*IRAS 60 microns*') eq 1B) and $
         (strmatch(subtext,'*1990IRASF.C...0000M*') eq 1B),niras60)

       if (niras60 ne 0L) then begin

          flux = strmid(subtext[iras60],26,11)

          if (niras60 gt 1L) then begin
             splog, 'WARNING: '+data[i].galaxy+' has multiple IRAS 60-micron [MOSHIR] measurements:', flux
             iras60 = iras60[0]
             flux = flux[0]
          endif
          
          if strcompress(flux,/remove) ne '' then begin 
             if strmatch(flux,'*<*') eq 1B then begin ; upper limit
                flux = -float(strmid(flux,strpos(flux,'<')+1,strlen(flux)))
             endif
             data[i].moshir_iras_60 = float(flux)
          endif

          err = strmid(subtext[iras60],41,7)
          if strcompress(err,/remove) ne '' then begin
             if strmatch(err,'*%*') eq 1B then err = (strmid(err,0,strpos(err,'%'))/100.0)*abs(data[i].moshir_iras_60)
             data[i].moshir_iras_60_err = float(err)
          endif

          units = strcompress(strmid(subtext[iras60],48,11),/remove)
          if strmatch(units,'Jy',/fold) ne 1B then splog, 'IRAS 60 UNITS WARNING: ', i, data[i].galaxy
          
       endif

; --------------------------------------------------
; IRAS 100 micron       
; --------------------------------------------------

; ##################################################       
; SANDERS
; ##################################################       
       
       iras100 = where((strmatch(subtext,'*IRAS 100 microns*') eq 1B) and $
         (strmatch(subtext,'*2003AJ....126.1607S*') eq 1B),niras100)

       if (niras100 ne 0L) then begin

          flux = strmid(subtext[iras100],26,11)

          if (niras100 gt 1L) then begin
             splog, 'WARNING: '+data[i].galaxy+' has multiple IRAS 100-micron [SANDERS] measurements:', flux
             iras100 = iras100[0]
             flux = flux[0]
          endif
          
          if strcompress(flux,/remove) ne '' then begin 
             if strmatch(flux,'*<*') eq 1B then begin ; upper limit
                flux = -float(strmid(flux,strpos(flux,'<')+1,strlen(flux)))
             endif
             data[i].sanders_iras_100 = float(flux)
          endif

          err = strmid(subtext[iras100],41,7)
          if strcompress(err,/remove) ne '' then begin
             if strmatch(err,'*%*') eq 1B then err = (strmid(err,0,strpos(err,'%'))/100.0)*abs(data[i].sanders_iras_100)
             data[i].sanders_iras_100_err = float(err)
          endif

          units = strcompress(strmid(subtext[iras100],48,11),/remove)
          if strmatch(units,'Jy',/fold) ne 1B then splog, 'IRAS 100 UNITS WARNING: ', i, data[i].galaxy
          
       endif

; ##################################################       
; RICE
; ##################################################       
       
       iras100 = where((strmatch(subtext,'*IRAS 100 microns*') eq 1B) and $
         (strmatch(subtext,'*1988ApJS...68...91R*') eq 1B),niras100)

       if (niras100 ne 0L) then begin

          flux = strmid(subtext[iras100],26,11)

          if (niras100 gt 1L) then begin
             splog, 'WARNING: '+data[i].galaxy+' has multiple IRAS 100-micron [RICE] measurements:', flux
             iras100 = iras100[0]
             flux = flux[0]
          endif
          
          if strcompress(flux,/remove) ne '' then begin 
             if strmatch(flux,'*<*') eq 1B then begin ; upper limit
                flux = -float(strmid(flux,strpos(flux,'<')+1,strlen(flux)))
             endif
             data[i].rice_iras_100 = float(flux)
          endif

          err = strmid(subtext[iras100],41,7)
          if strcompress(err,/remove) ne '' then begin
             if strmatch(err,'*%*') eq 1B then err = (strmid(err,0,strpos(err,'%'))/100.0)*abs(data[i].rice_iras_100)
             data[i].rice_iras_100_err = float(err)
          endif

          units = strcompress(strmid(subtext[iras100],48,11),/remove)
          if strmatch(units,'Jy',/fold) ne 1B then splog, 'IRAS 100 UNITS WARNING: ', i, data[i].galaxy
          
       endif

; ##################################################       
; SOIFER
; ##################################################       
       
       iras100 = where((strmatch(subtext,'*IRAS 100 microns*') eq 1B) and $
         (strmatch(subtext,'*1989AJ.....98..766S*') eq 1B),niras100)

       if (niras100 ne 0L) then begin

          flux = strmid(subtext[iras100],26,11)

          if (niras100 gt 1L) then begin
             splog, 'WARNING: '+data[i].galaxy+' has multiple IRAS 100-micron [SOIFER] measurements:', flux
             iras100 = iras100[0]
             flux = flux[0]
          endif
          
          if strcompress(flux,/remove) ne '' then begin 
             if strmatch(flux,'*<*') eq 1B then begin ; upper limit
                flux = -float(strmid(flux,strpos(flux,'<')+1,strlen(flux)))
             endif
             data[i].soifer_iras_100 = float(flux)
          endif

          err = strmid(subtext[iras100],41,7)
          if strcompress(err,/remove) ne '' then begin
             if strmatch(err,'*%*') eq 1B then err = (strmid(err,0,strpos(err,'%'))/100.0)*abs(data[i].soifer_iras_100)
             data[i].soifer_iras_100_err = float(err)
          endif

          units = strcompress(strmid(subtext[iras100],48,11),/remove)
          if strmatch(units,'Jy',/fold) ne 1B then splog, 'IRAS 100 UNITS WARNING: ', i, data[i].galaxy
          
       endif

; ##################################################       
; MOSHIR
; ##################################################       
       
       iras100 = where((strmatch(subtext,'*IRAS 100 microns*') eq 1B) and $
         (strmatch(subtext,'*1990IRASF.C...0000M*') eq 1B),niras100)

       if (niras100 ne 0L) then begin

          flux = strmid(subtext[iras100],26,11)

          if (niras100 gt 1L) then begin
             splog, 'WARNING: '+data[i].galaxy+' has multiple IRAS 100-micron [MOSHIR] measurements:', flux
             iras100 = iras100[0]
             flux = flux[0]
          endif
          
          if strcompress(flux,/remove) ne '' then begin 
             if strmatch(flux,'*<*') eq 1B then begin ; upper limit
                flux = -float(strmid(flux,strpos(flux,'<')+1,strlen(flux)))
             endif
             data[i].moshir_iras_100 = float(flux)
          endif

          err = strmid(subtext[iras100],41,7)
          if strcompress(err,/remove) ne '' then begin
             if strmatch(err,'*%*') eq 1B then err = (strmid(err,0,strpos(err,'%'))/100.0)*abs(data[i].moshir_iras_100)
             data[i].moshir_iras_100_err = float(err)
          endif

          units = strcompress(strmid(subtext[iras100],48,11),/remove)
          if strmatch(units,'Jy',/fold) ne 1B then splog, 'IRAS 100 UNITS WARNING: ', i, data[i].galaxy
          
       endif

; --------------------------------------------------       
; UV 1650
; --------------------------------------------------       

       uv1650 = where((strmatch(subtext,'*UV_1650 (m_T)*') eq 1B) and $
         (strmatch(subtext,'*1995A&AS..114..527R*') eq 1B),nuv1650)
       if (nuv1650 ne 0L) then begin

          mag = float(strmid(subtext[uv1650],27,11))
          if strmatch(mag,'*>*') ne 1B then begin ; exclude upper limits

             data[i].UV_1650_REF = 'Rifatto et al. 1995'

             flux = 10.0^(-0.4*(mag+21.1))
             data[i].UV_1650 = flux ; erg/s/cm2/A
             
             magerr = strmid(subtext[uv1650],41,6)
             if (strcompress(magerr,/remove) ne '') then begin
                fluxerr = flux*magerr*alog(10.0)/2.5 ; erg/s/cm2/A
                data[i].UV_1650_ERR = fluxerr
             endif

          endif
          
       endif 
       
; --------------------------------------------------       
; UV 2500
; --------------------------------------------------       

       uv2500 = where((strmatch(subtext,'*UV_2500 (m_T)*') eq 1B) and $
         (strmatch(subtext,'*1995A&AS..114..527R*') eq 1B),nuv2500)
       if (nuv2500 ne 0L) then begin

          mag = float(strmid(subtext[uv2500],27,11))
          if strmatch(mag,'*>*') ne 1B then begin ; exclude upper limits

             data[i].UV_2500_REF = 'Rifatto et al. 1995'

             flux = 10.0^(-0.4*(mag+21.1))
             data[i].UV_2500 = flux ; erg/s/cm2/A
             
             magerr = strmid(subtext[uv2500],41,6)
             if (strcompress(magerr,/remove) ne '') then begin
                fluxerr = flux*magerr*alog(10.0)/2.5 ; erg/s/cm2/A
                data[i].UV_2500_ERR = fluxerr
             endif

          endif
          
       endif
       
; --------------------------------------------------       
; UV 3150
; --------------------------------------------------       

       uv3150 = where((strmatch(subtext,'*UV_3150 (m_T)*') eq 1B) and $
         (strmatch(subtext,'*1995A&AS..114..527R*') eq 1B),nuv3150)
       if (nuv3150 ne 0L) then begin

          mag = float(strmid(subtext[uv3150],27,11))
          if strmatch(mag,'*>*') ne 1B then begin ; exclude upper limits

             data[i].UV_3150_REF = 'Rifatto et al. 1995'

             flux = 10.0^(-0.4*(mag+21.1))
             data[i].UV_3150 = flux ; erg/s/cm2/A
             
             magerr = strmid(subtext[uv3150],41,6)
             if (strcompress(magerr,/remove) ne '') then begin
                fluxerr = flux*magerr*alog(10.0)/2.5 ; erg/s/cm2/A
                data[i].UV_3150_ERR = fluxerr
             endif

          endif 

       endif
       
    endfor

; write a binary fits table
    
    if not keyword_set(outpath) then outpath = cwd()
    if keyword_set(textfile) then begin
       if not keyword_set(outfile) then outfile = 'outfile.txt'
       splog, 'Writing '+outpath+outfile+'.'
       print_struct, data, file=outpath+outfile
    endif else begin
       if not keyword_set(outfile) then outfile = 'outfile.fits'
       splog, 'Writing '+outpath+outfile+'.gz.'
       mwrfits, data, outpath+outfile, /create
       spawn, ['gzip -f '+outpath+outfile], /sh
    endelse
    
return
end
