;+
; NAME:
;   SPECTRAL_INDICES()
;
; PURPOSE:
;   Measure a variety of spectral indices from a composite stellar
;   population spectrum.
;
; INPUTS:
;   wave     - wavelength vector in Angstrom [NPIX]
;   flux     - flux density in erg/s/cm2/Angstrom [NPIX]
;
; OPTIONAL INPUTS:
;   ivar      - corresponding inverse variance spectrum [NPIX] (set
;     IVAR=0 for crummy pixels)
;   indexfile - full path name of the file listing the indices
;               (default '${ISPEC_DIR}/etc/indexlist.dat') 
;
; KEYWORD PARAMETERS:
;   nobreaks - do not measure any of the spectral breaks like D(4000) 
;   nolick   - do not measure any of the Lick indices
;   debug    - generate a plot of the index measurements and wait
;              for a keystroke
;   silent   - suppress messages to STDOUT
;
; OUTPUTS:
;   indices   - data structure with all the results and errors
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   We measure D4000 (Bruzual 1983, ApJ, 273, 105), D4000_narrow
;   (Balogh et al. 1999, ApJ, 527, 54), the 41-50 color index
;   (Kennicutt 1992, ApJS, 79, 255), the Balmer break (Leonardi &
;   Rose 2003, AJ, 126, 1811), all the continuum indices defined
;   in INDEXFILE, and several composite indices (see below for the
;   references). 
;
; EXAMPLE:
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2002 October 31, U of A
;   jm03oct30uofa - major overhaul and updates
;   jm03nov25uofa - incorporated IABSLINEEW() and cleaned up the
;                   code and documentation
;   jm04mar10uofa - bug fix in call to IABSLINEEW() (was not
;                   passing the line limits); added SILENT keyword 
;   jm05jan28uofa - bug fix when computing c41-50 to avoid NaN
;   jm07mar13nyu  - additional error checking when computing
;                   extended Lick indices
;   jm08aug14nyu  - when computing MgFe, Mg1Fe, and Mg2Fe check to
;                   make sure everything is positive!
;   jm08aug22nyu  - removed INDEXPATH variable; now INDEXLIST must
;                   include the full path name
;   jm08nov24nyu  - added NOBREAKS keyword
;   jm09aug20ucsd - added NOLICK keyword
;   jm09dec08ucsd - changed FERR optional input to IVAR (inverse
;     variance), so that bad pixels could be tracked (IVAR=0)
;
; Copyright (C) 2002-2009, John Moustakas
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

function spectral_indices, wave, flux, ivar=ivar1, indexfile=indexfile, $
  nobreaks=nobreaks, nolick=nolick, debug=debug, silent=silent, $
  empty_structure=empty_structure, _extra=extra

; index definitions    
    if (n_elements(indexfile) eq 0L) then indexfile = getenv('IMPRO_DIR')+$
      '/etc/indexlist.dat'
    if (file_test(indexfile,/regular) eq 0L) then begin
       splog, 'Lick index file '+indexfile+' not found'
       return, -1L
    endif

; initialize the output data structure; note that NAMES needs to
; include *all* indices measured by this routine, including all the
; extended and additional indices, below
    readcol, indexfile, licknames, w1, w2, wb1, wb2, $
      wr1, wr2, units, format='A,D,D,D,D,D,D,A', /silent, comment='#'
    nlick = n_elements(licknames)
    names = ['DBalmer','D4000','D4000_Narrow',licknames]
;   names = ['DBalmer','D4000','D4000_Narrow','CaII_Strength','c41_50']

; add composite index names
    if (total(strmatch(licknames,'Lick_Fe5270',/fold)+$
      strmatch(licknames,'Lick_Fe5335',/fold)) eq 2) then $
        names = [names,'Lick_Fe']
    
    if (total(strmatch(licknames,'Lick_Fe5270',/fold) and $
      strmatch(licknames,'Lick_Fe5335',/fold) and $
      strmatch(licknames,'Lick_Mgb',/fold)) eq 3) then $
        names = [names,'Lick_MgFe']
    if (total(strmatch(licknames,'Lick_Fe4531',/fold) and $
      strmatch(licknames,'Lick_Fe5015',/fold) and $
      strmatch(licknames,'Lick_Mg1',/fold) and $
      strmatch(licknames,'Lick_Mg2',/fold)) eq 4) then $
       names = [names,'Lick_Mg1Fe','Lick_Mg2Fe']
    indices = create_struct('indices', names)
    for j = 0L, n_elements(names)-1L do indices = create_struct($
      temporary(indices),names[j],[0.0,-2.0])
    if keyword_set(empty_structure) then return, indices ; return

; now error-check the input spectrum    
    npix = n_elements(wave)
    nflux = n_elements(flux)

    if (npix eq 0L) or (nflux eq 0L) then begin
       doc_library, 'spectral_indices'
       return, -1L
    endif

    if (n_elements(ivar1) eq 0) then begin
       sig = djsig(flux)
       if (sig le 0) then sig = 1.0
       ivar = flux*0.0+1.0/sig^2
    endif else ivar = ivar1

    neg = where(ivar lt 0,nneg)
    if (nneg ne 0) then message, 'IVAR array contains negative values'
    nivar = n_elements(ivar)
    
    if (npix ne nflux) or (npix ne nivar) then begin
       print, 'Dimensions of WAVE, FLUX, and IVAR do not agree'
       return, -1L
    endif

; convenient variables to convert to the flux vector (erg/s/cm2/A) to
; F_nu (erg/s/cm2/Hz); where IVAR is unity then we need to take that
; into account when converting to FNU_IVAR
    light = 2.99792458D18
    factor = wave^2/light
    fnu = flux*factor
    fnu_ivar = ivar/factor^2
    unity = where(ivar eq 1.0,nunity)
    if (nunity ne 0) then fnu_ivar[unity] = 1.0
    nu = light/wave

; measure the Balmer Break (Leonardi & Rose 2003, AJ, 126, 1811); the
; following two definitions of the 4000-Angstrom break: Bruzual 1983,
; ApJ, 273, 105 and Balogh et al. 1999, ApJ, 527, 54; the ratio of the
; strength of the Ca II absorption lines (Stoughton et al. 2002, AJ,
; 123, 485; and the 41-50 continuum color index (Kennicutt 1992, ApJS,
; 79, 255)
    breaks = {$
      DBalmer:       [0.0,-2.0],$
      D4000:         [0.0,-2.0],$
      D4000_Narrow:  [0.0,-2.0]}
;     CaII_Strength: [0.0,-2.0],$
;     c41_50:        [0.0,-2.0]}

    breaknames = tag_names(breaks)
    nbreaks = n_elements(breaknames)

    wbvec = [$             ; blue wavelength interval
      [3525.0,3600.0],$
      [3750.0,3950.0],$
      [3850.0,3950.0]]
;     [3921.0,3946.0],$
;     [4050.0,4250.0]]
    llimit = total(wbvec,1)/2.0
    lwidth = wbvec[1,*]-wbvec[0,*]
    
    wrvec = [ $             ; red wavelength interval
      [3700.0,3825.0],$
      [4050.0,4250.0],$
      [4000.0,4100.0]]
;     [3956.0,3981.0],$
;     [4900.0,5100.0]]
    ulimit = total(wrvec,1)/2.0
    uwidth = wrvec[1,*]-wrvec[0,*]

    midwave = total([[llimit],[ulimit]],2)/2.0 ; central wavelength
    blabel = repstr(breaknames,'_',' ')

    if (keyword_set(nobreaks) eq 0) then begin
;      if (keyword_set(silent) eq 0) then splog, 'Measuring continuum breaks and color indices'
       cbreak = iabslineew(wave,fnu,midwave,ivar=fnu_ivar,llimit=llimit,$
         lwidth=lwidth,ulimit=ulimit,uwidth=uwidth,label=blabel,/noline,$
         absplot=cbreakplot,debug=debug,/fnu,silent=silent,_extra=extra)

; fill the structure and copy the results to INDICES
       for ibreak = 0L, nbreaks-1L do begin
          if strmatch(breaknames[ibreak],'*41*50*') eq 1B then begin
             if (cbreak[ibreak].cratio gt 0.0) then begin
                bk = 2.5*alog10(cbreak[ibreak].cratio)
                bkerr = (2.5/alog(10.0))*cbreak[ibreak].cratio_err/cbreak[ibreak].cratio
             endif else begin
                bk = 0.0
                bkerr = -2.0
             endelse
          endif else begin
             bk = cbreak[ibreak].cratio
             bkerr = cbreak[ibreak].cratio_err
          endelse
          breaks.(ibreak) = [bk,bkerr]
       endfor
       struct_assign, breaks, indices, /nozero
    endif

; ---------------------
; extended Lick indices
    if (keyword_set(nolick) eq 0) then begin

; define a temporary lick structure
       lick = create_struct(licknames[0],[0.0,-2.0])
       for j = 1L, nlick-1L do lick = create_struct(lick,licknames[j],[0.0,-2.0])
       
; measure the indices; redefine some variables for use in IABSLINEEW() 
       llimit = total([[wb1],[wb2]],2)/2.0
       ulimit = total([[wr1],[wr2]],2)/2.0
       linewave = total([[w1],[w2]],2)/2.0
       lwidth = wb2-wb1 & uwidth = wr2-wr1
       lline = w1 & uline = w2
       llabel = repstr(licknames,'_',' ')
; custom labels
       llabel = repstr(llabel,'Lick Hd A','Lick H\delta_{A}')
       llabel = repstr(llabel,'Lick Hg A','Lick H\gamma_{A}')

       mags = where(units eq 'mag',nmags,comp=EWs,ncomp=nEWs)
       if nmags ne 0L then begin
;         if (keyword_set(silent) eq 0) then splog, 'Measuring spectral index equivalent widths [magnitude].'
          lick_mags = iabslineew(wave,flux,linewave[mags],ivar=ivar,llimit=llimit[mags],$
            lwidth=lwidth[mags],ulimit=ulimit[mags],uwidth=uwidth[mags],$
            lline=lline[mags],uline=uline[mags],label=llabel[mags],absplot=lickplot,$
            debug=debug,/magnitude,silent=silent,_extra=extra)
          for ilick = 0L, nmags-1L do lick.(mags[ilick]) = $
            [lick_mags[ilick].lineew,lick_mags[ilick].lineew_err]
       endif
       
       if nEWs ne 0L then begin
;         if (keyword_set(silent) eq 0) then splog, 'Measuring spectral index equivalent widths [Angstrom].'
          lick_ews = iabslineew(wave,flux,linewave[EWs],ivar=ivar,llimit=llimit[EWs],$
            lwidth=lwidth[EWs],ulimit=ulimit[EWs],uwidth=uwidth[EWs],lline=lline[EWs],$
            uline=uline[EWs],label=llabel[EWs],absplot=lickplot,debug=debug,$
            silent=silent,_extra=extra)
          for ilick = 0L, nEWs-1L do lick.(EWs[ilick]) = $
            [lick_ews[ilick].lineew,lick_ews[ilick].lineew_err]
       endif

; now copy over the results
       struct_assign, lick, indices, /nozero
    endif

; -----------------------------------------
; finally compute several additional useful indices; see Thomas,
; Maraston, & Bender 2003, MNRAS, 339, 897 and Bruzual & Charlot 2003,
; MNRAS, 344, 1000
    if tag_exist(indices,'LICK_FE5270') and tag_exist(indices,'LICK_FE5335') and $
      tag_exist(indices,'LICK_Fe') then begin
       if (indices.lick_fe5270[1] gt 0.0) and (indices.lick_fe5335[1] gt 0.0) then begin
          indices.Lick_Fe[0] = 0.5*(indices.lick_fe5270[0] + indices.lick_fe5335[0])
          indices.Lick_Fe[1] = 0.5*sqrt(indices.lick_fe5270[1]^2 + indices.lick_fe5335[1]^2)
       endif
    endif

    if tag_exist(indices,'LICK_MGB') and tag_exist(indices,'LICK_FE5270') and $
      tag_exist(indices,'LICK_FE5335') and tag_exist(indices,'LICK_MgFe') then begin

       u1 = indices.lick_mgb[0]
       u1_err = indices.lick_mgb[1]
       u2 = 0.72*indices.lick_fe5270[0] + 0.28*indices.lick_fe5335[0]
       u2_err = sqrt((0.72*indices.lick_fe5270[1])^2 + (0.28*indices.lick_fe5335[1])^2)

       if (u1 gt 0.0) and (u2 gt 0.0) and (u1_err gt 0.0) and (u2_err gt 0.0) then begin
          indices.Lick_MgFe[0] = sqrt(u1*u2)
          indices.Lick_MgFe[1] = sqrt((u1*u2_err)^2 + (u2*u1_err)^2) / (2*indices.Lick_MgFe[0])
       endif
    endif

    if tag_exist(indices,'LICK_MG1') and tag_exist(indices,'LICK_MG2') and $
      tag_exist(indices,'LICK_FE4531') and tag_exist(indices,'LICK_FE5015') and $
      tag_exist(indices,'LICK_Mg1Fe') and tag_exist(indices,'LICK_Mg2Fe') then begin
       if (indices.lick_fe4531[1] gt 0.0) and (indices.lick_fe5015[1] gt 0.0) and $
         (indices.lick_mg1[1] gt 0.0) and (indices.lick_mg2[1] gt 0.0) then begin
          
          u = indices.lick_fe4531[0] + indices.lick_fe5015[0]
          u_err = sqrt(indices.lick_fe4531[1]^2 + indices.lick_fe5015[1]^2)

          if (u gt 0.0) then begin
             indices.Lick_Mg1Fe[0] = 0.6*indices.lick_mg1[0] + 0.4*alog10(u)
             indices.Lick_Mg1Fe[1] = sqrt((0.6*indices.lick_mg1[1])^2 + ((0.4/alog(10.0))*u_err/u)^2)

             indices.Lick_Mg2Fe[0] = 0.6*indices.lick_mg2[0] + 0.4*alog10(u)
             indices.Lick_Mg2Fe[1] = sqrt((0.6*indices.lick_mg2[1])^2 + ((0.4/alog(10.0))*u_err/u)^2)
          endif
       endif
    endif

; quick error check    
    if (n_tags(struct_trimtags(indices,except=['indices'])) ne $
      n_elements(indices.indices)) then splog, 'WARNING: Incomplete INDICES vector!'
    
return, indices
end
