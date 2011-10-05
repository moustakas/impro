;+
; NAME:
;   IM_READ_BC03()
;
; PURPOSE:
;   Read the Bruzual & Charlot (2003) binary format population
;   synthesis models into a convenient data structure.
;
; CALLING SEQUENCE:
;   bc03 = im_read_bc03(isedfile=,isedpath,metallicity=,$
;      age=,minwave=,maxwave=,endian=,bc03_extras=,$
;      /salpeter,/lr,/readzero,/silent,/array_extras,/float)
;
; INPUTS:
;   None required.  By default this routine reads the high
;   resolution (hr), solar metallicity (m62), Chabrier SSP
;   models. 
;
; OPTIONAL INPUTS:
;   isedfile - read this binary SED rather than the default SSP
;              models 
;   isedpath - full data path (with trailing slash) to ISEDFILE
;   metallicity - SSP metallicity
;      0 - Z=0.0001 (m22)
;      1 - Z=0.0004 (m32)
;      2 - Z=0.004  (m42)
;      3 - Z=0.008  (m52)
;      4 - Z=0.02   (m62) [Solar, default]
;      5 - Z=0.05   (m72)
;   age - return the SSP(s) corresponding to this scalar or vector
;         age(s) [Gyr]
;   minwave - crop the SSP spectra to this minimum wavelength
;             [Angstrom] 
;   maxwave - crop the SSP spectra to this maximum wavelength
;             [Angstrom]
;   endian  - byte ordering for the binary file; see documentation
;             for READ_BINARY(); usually, if the binary files have
;             been created on a Linux box then endian='little',
;             which is the default; on a SUN, endian='big' is
;             appropriate (for example, if you use 'bin_ised'
;             provided with BC03 to create binary-format SSP's) 
;
; KEYWORD PARAMETERS:
;   salpeter - read the Salpeter IMF models (default is to read
;              the Chabrier models)
;   lr       - read the low resolution models (default is to read
;              the high resolution models)
;   readzero - read the zeroth age SSP (has not been tested with
;              convolved star-formation histories!)
;   silent   - do not print any messages to STDOUT
;   array_extras - return the BC03_EXTRAS structure in a slightly
;                  different format
;   float    - read the AGE, WAVE, and FLUX arrays in single (not 
;              double) precision
;   abmag - convert the output spectra to AB mag at 10 pc
;
; OUTPUTS:
;   bc03 - data structure with the following fields:
;      age  - vector of SSP ages [NAGE] [yr]
;      wave - wavelength array [NPIX] [Angstrom]
;      flux - SSP models [NPIX,NAGE] [L_sun/M_sun/A]
;
; OPTIONAL OUTPUTS:
;   bc03_extras - data structure containing the extra parameters
;                 associated with each SSP (see the BC03
;                 documentation) 
;
; COMMENTS:
;   N.B.  The environment variable ${bc03_dir} must be defined in
;   your '.cshrc' indicating the *root* directory of the BC03
;   models.  For example, in my '.cshrc' I have
;
;          setenv bc03_dir ${HOME}/synthesis/bc03
;
;   Note that this routine only works with the Padova (1994)
;   binary isochrones and not the Pickles isochrones.
;   Alternatively, you can use the ISEDFILE and ISEDPATH optional
;   inputs to read in an arbitrary BC03 SED.
;
; BUGS:
;
; INTERNAL SUPPORT ROUTINES:
;   READ_EXTRAS(), GET_ELEMENT
;
; PROCEDURES USED:
;   READFAST, STRUCT_TRIMTAGS(), STRUCT_ADDTAGS(), MATCH,
;   MRD_STRUCT() 
;
; EXAMPLES:
;   [1] Read the high-resolution, Salpeter IMF, Solar metallicity
;       (Z=0.02) SSP models:
;
;      IDL> bc03 = im_read_bc03()
;      IDL> help, bc03, /str
;
;   [2] Read the low-resolution, Salpeter IMF, LMC metallicity
;       (Z=0.004) SSP models:
;
;      IDL> bc03 = im_read_bc03(metallicity=2,/lr,/salpeter)
;
;   [3] Read the high-resolution, Chabrier IMF, twice-solar
;       metallicity (Z=0.05) SSP models and plot the 10 Gyr model: 
;
;      IDL> bc03 = im_read_bc03(metallicity=5)
;      IDL> indx = where(bc03.age eq 1E10)
;      IDL> plot, bc03.wave, bc03.flux[*,indx], /xlog, /ylog
;
;   [4] Read the high-resolution, Chabrier IMF, solar metallicity
;       (Z=0.02) SSP models and extras and plot D4000 versus
;       H-delta_A:  
;
;      IDL> bc03 = im_read_bc03(bc03_extras=bce)
;      IDL> plot, bce.d_4000_, bce.h_delta_a, xsty=3, ysty=3
;
;   [5] Read a non-SSP SED in my temp subdirectory:
;
;      IDL> bc03 = im_read_bc03(isedfile='mysed.ised',isedpath='temp/')
;
;   [6] Retrieve the 12 Gyr Salpeter IMF, high-resolution SSP
;       between 3000 and 10000 Angstroms.
;
;      IDL> bc03 = im_read_bc03(/salpeter,age=12.0,minwave=3000,maxwave=1E4)
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2003 October 29, U of A, based in small part on
;      IDL code originally written by C. Papovich
;   jm03nov12uofa - added ISEDFILE and ISEDPATH optional inputs
;                   and various bug fixes
;   jm04mar01uofa - added AGE, MINWAVE, and MAXWAVE optional
;                   inputs; changed wave, flux, and age to double
;                   precision 
;   jm04sep05uofa - generalized to read non-SSP models; NAGE is no
;                   longer fixed to be 221 age bins; improved
;                   error checking when reading the "extras"
;                   files; do not read AGE=0
;   jm04oct31uofa - added ENDIAN optional keyword (thanks to
;                   M. Blanton) 
;   jm05jul26uofa - added READZERO optional input
;   jm06feb19uofa - added ARRAY_EXTRAS and FLOAT keywords 
;
; Copyright (C) 2003-2005, John Moustakas
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

function read_extras, extrafile, extrapath=extrapath, array_extras=array_extras, silent=silent

    if file_test(extrapath+extrafile) eq 0L then begin
       print, 'Extras file '+extrapath+extrafile+' not found.'
       return, -1L
    endif
       
    if not keyword_set(silent) then print, 'Reading extras file '+extrapath+extrafile

; read the first 50 lines of the file to figure out the column names
; and also to determine on which row the data begin 

    temp = strarr(50L)
    tmpstr = ''
    openr, lun1, extrapath+extrafile, /get_lun
    for i = 0L, 49L do begin
       readf, lun1, tmpstr
       temp[i] = tmpstr
    endfor
    free_lun, lun1

    problem = 0L
    
    headlines = where(strmatch(temp,'*#*') eq 1B,nheadlines)
    if (nheadlines eq 0L) then problem = 1L

    agerow = where(strmatch(temp,'*log-age*',/fold) eq 1B,nagerow)
    if (nagerow ne 1L) then problem = 1L

    if problem then begin
    
       splog, 'WARNING: There was a problem reading with '+extrapath+extrafile+'.'
       return, -1L

    endif

; read the column names and convert them to valid structure tag names
    
    cols = strupcase(strcompress(strsplit(strmid(temp[agerow],1),' ',/extract),/remove))
    cols[0] = 'LOGAGE'
    ncols = n_elements(cols)
    for j = 0L, ncols-1L do cols[j] = idl_validname(cols[j],/convert_all)

; read the data
    
    readfast, extrapath+extrafile, data, header, skipline=nheadlines, $
      nlines=nage, ncols=ncheck

    if (ncheck ne ncols) then begin
       splog, 'There was a problem reading '+extrapath+extrafile+'.'
       return, -1L
    endif

; initialize the output data structure and fill it

    if keyword_set(array_extras) then begin

       extras = mrd_struct(cols,replicate('fltarr('+string(nage,format='(I0)')+')',ncols),1.0)
       for icol = 0L, ncols-1L do extras.(icol) = data[icol,*]
       
    endif else begin
       
       extras = create_struct(cols[0],0.0)
       for k = 1L, ncols-1L do extras = create_struct(extras,cols[k],0.0)
       extras = replicate(extras,nage)

       for iage = 0L, ncols-1L do extras.(iage) = reform(data[iage,*])

    endelse
       
return, extras
end

pro get_element, x, value, position
; jm01jan28uofa
; a generalization of GETELEMENT_VECTOR, this routine will also accept
; an array of positions, and return an array of indices

    position = long(value-value)
    for i = 0L, n_elements(value)-1L do begin
       array_value = min((abs(x-value[i])),temp)
       position[i] = temp
    endfor
      
return
end

function im_read_bc03, isedfile=isedfile, isedpath=isedpath, metallicity=metallicity, $
  age=age, minwave=minwave, maxwave=maxwave, endian=endian, salpeter=salpeter, lr=lr, $
  readzero=readzero, silent=silent, array_extras=array_extras, bc03_extras=bc03_extras, $
  float=float, abmag=abmag

    if n_elements(metallicity) eq 0L then metallicity = 4L
    if n_elements(endian) eq 0L then endian = 'little'

    case long(metallicity) of
       0L: begin
          Zstr = 'm22'
          Zinfo = 'Z=0.0001'
       end
       1L: begin
          Zstr = 'm32'
          Zinfo = 'Z=0.0004'
       end
       2L: begin
          Zstr = 'm42'
          Zinfo = 'Z=0.004'
       end
       3L: begin
          Zstr = 'm52'
          Zinfo = 'Z=0.008'
       end
       4L: begin
          Zstr = 'm62'
          Zinfo = 'Z=0.02'
       end
       5L: begin
          Zstr = 'm72'
          Zinfo = 'Z=0.05'
       end
       else: begin
          Zstr = 'm62'
          Zinfo = 'Z=0.02'
       end
    endcase

    if n_elements(isedpath) eq 0L then begin
       isedpath = getenv('bc03_dir')+'/models/Padova1994/'
       if keyword_set(salpeter) then $
         isedpath = isedpath+'salpeter/' else $
         isedpath = isedpath+'chabrier/'
    endif

    if keyword_set(salpeter) then begin
       imfstr = 'salp'
       imfinfo = 'Salpeter IMF' 
    endif else begin
       imfstr = 'chab'
       imfinfo = 'Chabrier IMF'
    endelse
    
    if keyword_set(lr) then begin ; low resolution

       npix = 1221L        
       resstr = 'lr'
       resinfo = 'Low Resolution'

       if keyword_set(salpeter) then begin
          ifs = 306L   ; starting index of the wavelength vector
          ifs2 = 2807L ; starting index of the first spectrum
       endif else begin
          ifs = 300L
          ifs2 = 2801L
       endelse

    endif else begin              ; high resolution

       npix = 6900L 
       resstr = 'hr'
       resinfo = 'High Resolution'

       if keyword_set(salpeter) then begin
          ifs = 306L   
          ifs2 = 7209L 
       endif else begin
          ifs = 300L
          ifs2 = 7203L
       endelse
       
    endelse
      
; read the binary file

    if n_elements(isedfile) eq 0L then $
      isedfile = 'bc2003_'+resstr+'_'+Zstr+'_'+imfstr+'_ssp.ised' else $
      Zinfo = 'Z=Unknown'

    if file_test(isedpath+isedfile) then begin
       if not keyword_set(silent) then begin
          print, 'Reading SSP file '+isedpath+isedfile
          print, imfinfo+', '+resinfo+', '+Zinfo
       endif
       tempbin = read_binary(isedpath+isedfile,data_type=4,endian=endian)
    endif else begin
       print, 'SSP file '+isedpath+isedfile+' not found.'
       return, -1L
    endelse
        
; initialize some indices
    
    offa = 2L   ; offset from beginning of file to first age index
    offs = 56L  ; space between spectra
;   nage = 221L ; number of age bins (and models)

; retrieve the number of age records in the file (requires testing!)
    
    ntimesteps = read_binary(isedpath+isedfile,data_type=2,data_dims=1,$
      data_start=4,endian=endian)
    ntimesteps = ntimesteps[0]

; this additional offset accounts for non-SSP models, which have more
; than 221 age bins

    ifs = ifs + ntimesteps - 221L
    ifs2 = ifs2 + ntimesteps - 221L

; initialize the output data structure; do not read age zero unless
; READZERO=1

    if keyword_set(readzero) then begin
       nage = ntimesteps
       istart = 0L
    endif else begin
       nage = ntimesteps - 1L
       istart = 1L              ; offset from age zero
    endelse
    
    if not keyword_set(silent) then splog, 'Found '+string(nage,format='(I0)')+' age bins.'

    if keyword_set(float) then double = 0L else double = 1L
    
    bc03 = {$
      age:         make_array(nage,double=double), $
      wave:        make_array(npix,double=double), $
      flux:        make_array(npix,nage,double=double)}
;     bolflux:     dblarr(nage),      $
;     stellarmass: dblarr(nage)}

    bc03.age = (tempbin[offa:offa+ntimesteps-1L])[istart:ntimesteps-1L] ; age vector
    bc03.wave = tempbin[ifs:ifs+npix-1L]                                ; wavelength vector

    for i = istart, ntimesteps-1L do begin

       i1 = ifs2 + i*(npix+offs)
       i2 = ifs2 + i*(npix+offs) + npix-1L
;      print, i1, i2, (i2-i1)+1

       bc03.flux[*,i-istart] = tempbin[i1:i2]

    endfor

;; read the bolometric flux and stellar mass as a function of time 
;
;    j1 = i2 + offs + 1L + 1L       ; the second "1" is to offset from age zero 
;    j2 = i2 + offs + ntimesteps
;    bc03.bolflux = tempbin[j1:j2]
;
;    k1 = j2 + 4L + 1L              ; offset from zero age
;    k2 = j2 + 4L + ntimesteps - 1L
;    bc03.stellarmass = tempbin[k1:k2]

    ninputage = n_elements(age)
    if (ninputage ne 0L) then begin

       get_element, bc03.age/1E9, age, ageindx
       if not keyword_set(silent) then begin
          splog, 'Retrieving the following SSP(s):'
          for j = 0L, ninputage-1L do print, '   '+string(bc03.age[ageindx[j]]/1E9,format='(F12.5)')+' Gyr'
       endif

       bc03 = {age: bc03.age[ageindx], wave: bc03.wave, flux: bc03.flux[*,ageindx]}
;      bolflux: bc03.bolflux[ageindx], stellarmass: bc03.stellarmass[ageindx]}
    
    endif

    nminwave = n_elements(minwave)
    nmaxwave = n_elements(maxwave)
    
    if (nminwave ne 0L) or (nmaxwave ne 0L) then begin

       if (nminwave eq 0L) then minwave = min(bc03.wave)
       if (nmaxwave eq 0L) then maxwave = max(bc03.wave)

       get_element, bc03.wave, [minwave,maxwave], ww
       if not keyword_set(silent) then begin
          splog, 'Reducing wavelength range to ['+string(minwave,format='(I0)')+','+$
            string(maxwave,format='(I0)')+'] Angstrom.'
       endif

       bc03 = {age: bc03.age, wave: bc03.wave[(ww[0]-1L)>0L:(ww[1]+1L)<(npix-1L)], $
         flux: bc03.flux[(ww[0]-1L)>0L:(ww[1]+1L)<(npix-1L),*]}
    
    endif

; convert to AB magnitudes, if desired    
    if keyword_set(abmag) then begin
       pc10 = 10.0*3.085678D18  ; =10 pc
       light = 2.99792458D18    ; [Angstrom/s]
       for ii = 0, nage-1 do begin
          bc03.flux[*,ii] = 3.826D33*bc03.flux[*,ii]/(4.0*!dpi*pc10^2)*bc03.wave^2/light
          gd = where(bc03.flux[*,ii] gt 0.0,ngd)
          if (ngd ne 0) then bc03.flux[gd,ii] = -2.5*alog10(bc03.flux[gd,ii])-48.6
       endfor
    endif
    
; read the extra parameters for this SSP
    if arg_present(bc03_extras) then begin

       if (n_elements(isedfile) eq 0L) then $       
         base = 'bc2003_'+resstr+'_'+Zstr+'_'+imfstr+'_ssp.' else $
         base = repstr(isedfile,'ised','')
       
       colorfiles = base+string(lindgen(5)+1,format='(I0)')+'color'
       ABmagfile = base+'1ABmag'
          
       if keyword_set(lr) then begin

          extrafiles = [colorfiles,ABmagfile]

       endif else begin

          indxfiles6 = base+'6lsindx_'+'sed' ; ['ffn','sed','sed_lick_system']
          indxfiles7 = base+'7lsindx_'+'sed' ; ['ffn','sed','sed_lick_system']

          extrafiles = [colorfiles,indxfiles6,indxfiles7,ABmagfile]

       endelse
          
       nextra = n_elements(extrafiles)

       for k = 0L, nextra-1L do begin

          extras1 = read_extras(extrafiles[k],extrapath=isedpath,silent=silent,array_extras=array_extras)
          if (size(extras1,/type) eq 8L) then begin
             
             if k eq 0L then bc03_extras = extras1 else begin

; remove repeated structure tag names before appending the data             
                
                oldtags = tag_names(bc03_extras)
                newtags = tag_names(extras1)

                match, oldtags, newtags, oldindx, newindx, count=count                
                if count ne 0L then extras1 = struct_trimtags(extras1,except=newtags[newindx])

                bc03_extras = struct_addtags(bc03_extras,extras1)

             endelse

          endif
             
       endfor

; crop the extras structure to the output ages requested
    
       if (ninputage ne 0L) then bc03_extras = bc03_extras[ageindx]

    endif 
    
return, bc03
end
