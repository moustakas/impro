;+
; NAME:
;   BUILD_HIIREGIONS_DATABASE
;
; PURPOSE:
;   Parse the individual HII-region datafiles and write out a
;   structure containing the fluxes in a couple different formats. 
;
; INPUTS:
;
; KEYWORD PARAMETERS:
;
; OUTPUTS:
;   ${CATALOGS_DIR}/hiiregions/hiiregions_database_VERSION.fits.gz
;   ${CATALOGS_DIR}/hiiregions/hiiregions_database_nolog_VERSION.fits.gz
;   ${CATALOGS_DIR}/hiiregions/hiiregions_linefit_VERSION.fits.gz
;
; COMMENTS:
;
; TODO:
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2010 Feb 22, UCSD - modularized version of
;     WRITE_HII_REGIONS (which had a long, long history)
;
; Copyright (C) 2010, John Moustakas
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

pro build_hiiregions_database, data, data_nolog, linefit, $
  nmonte=nmonte, clobber=clobber

    hiipath = hiiregions_path()
    version = hiiregions_version()

    outfile = hiipath+'hiiregions_database_'+version+'.fits'
    outfile_nolog = hiipath+'hiiregions_database_nolog_'+version+'.fits'
    linefit_outfile = hiipath+'hiiregions_linefit_'+version+'.fits'
    if file_test(outfile+'.gz',/regular) and $
      (keyword_set(clobber) eq 0) then begin
       splog, 'Output file '+outfile+' exists; use /CLOBBER'
       return
    endif

; see BUILD_HIIREGIONS_ANCILLARY
    nedfile = hiipath+'hii_region_ned_basic_'+version+'.fits.gz' 
    diamfile = hiipath+'hii_region_ned_diameters_'+version+'.fits.gz'
    if (file_test(nedfile,/regular) eq 0) then begin
       splog, 'NED basic file '+nedfile+' not found'
       return
    endif
    if (file_test(diamfile,/regular) eq 0) then begin
       splog, 'NED diameters file '+diamfile+' not found'
       return
    endif

    if (n_elements(nmonte) eq 0) then nmonte = 250
    
; ---------------------------------------------------------------------------    
; read and parse all the data files; treat the individual HII regions
; separately from the HII galaxies
    splog, 'Reading '+hiipath+'datafiles_hii_regions.txt'
    readcol, hiipath+'datafiles_hii_regions.txt', hii_datafile, hii_reference, $
      hii_texref, format='A,A,A', comment='#', delimiter='|', /silent
    hii_datafile = strtrim(hii_datafile,2)
    hii_reference = strtrim(hii_reference,2)
    hii_texref = strtrim(hii_texref,2)
    
    splog, 'Reading '+hiipath+'datafiles_hii_galaxies.txt'
    readcol, hiipath+'datafiles_hii_galaxies.txt', galaxy_datafile, galaxy_reference, $
      galaxy_texref, format='A,A,A', comment='#', delimiter='|', /silent
    galaxy_datafile = strtrim(galaxy_datafile,2)
    galaxy_reference = strtrim(galaxy_reference,2)
    galaxy_texref = strtrim(galaxy_texref,2)

; parse everything       
    splog, 'Parsing the HII-region data files'
    hii_data = parse_hiiregions_datafile(hiipath+hii_datafile,hii_reference,$
      hii_texref,linefit=hii_linefit,data_nolog=hii_data_nolog)
    newtags = replicate({hiiregion: 1, hiigalaxy: 0},n_elements(hii_data))
    hii_data = struct_addtags(hii_data,newtags)
    hii_data_nolog = struct_addtags(hii_data_nolog,newtags)

    splog, 'Parsing the HII-galaxy data files'
    galaxy_data = parse_hiiregions_datafile(hiipath+galaxy_datafile,galaxy_reference,$
      galaxy_texref,linefit=galaxy_linefit,data_nolog=galaxy_data_nolog)
    newtags = replicate({hiiregion: 0, hiigalaxy: 1},n_elements(galaxy_data))
    galaxy_data = struct_addtags(galaxy_data,newtags)
    galaxy_data_nolog = struct_addtags(galaxy_data_nolog,newtags)

; combine...
    data = struct_append(hii_data,galaxy_data)
    data_nolog = struct_append(hii_data_nolog,galaxy_data_nolog)
    linefit = struct_append(hii_linefit,galaxy_linefit)

; ---------------------------------------------------------------------------    
; merge the ancillary data from NED

    splog, 'Storing NED properties'

    splog, 'Reading '+nedfile
    splog, 'Reading '+diamfile
    ned = mrdfits(nedfile,1)
    diam = mrdfits(diamfile,1)

    doit = match_string(data.hii_galaxy,ned.galaxy,index=index,/exact)
    good = where(index ne -1L,ngood)
    niceprint, ned[index[good]].galaxy, doit[good], data[good].hii_galaxy

    data[good].ned_galaxy          = strupcase(ned[index[good]].ned_galaxy)
    data[good].galaxy_ra           = ned[index[good]].ra
    data[good].galaxy_dec          = ned[index[good]].dec
    data[good].galaxy_rc3_pa       = diam[index[good]].rc3_posangle
    data[good].galaxy_twomass_pa   = diam[index[good]].twomass_posangle
    data[good].galaxy_rc3_incl     = compute_inclination(diam[index[good]])
    data[good].galaxy_twomass_incl = compute_inclination(diam[index[good]],/twomass)

    indx = where(diam[index[good]].rc3_major_axis gt -900.0,nindx)
    if (nindx ne 0L) then data[good[indx]].galaxy_rc3_r25 = $
      diam[index[good[indx]]].rc3_major_axis/60.0/2.0

    indx = where(diam[index[good]].twomass_k20_major_axis gt -900.0,nindx)
    if (nindx ne 0L) then data[good[indx]].galaxy_twomass_k20 = $
      diam[index[good[indx]]].twomass_k20_major_axis/60.0/2.0

; set the inclination angle for M51 to 20 (Tully 1974)
    m51 = where(strmatch(data.hii_galaxy,'*5194*'),nm51)
;   niceprint, data[m51].hii_galaxy, data[m51].galaxy_rc3_pa, data[m51].galaxy_twomass_pa, $
;     data[m51].galaxy_rc3_incl, data[m51].galaxy_twomass_incl
    if (nm51 ne 0L) then begin
       splog, 'Setting M51 inclination angle equal to 20!'
       data[m51].galaxy_rc3_incl = 20.0
       data[m51].galaxy_twomass_incl = 20.0
    endif
    
; set the inclination and position angleS for NGC6822 (Lee & Skillman 2006)
    n6822 = where(strmatch(data.hii_galaxy,'*6822*'),nn6822)
;   niceprint, data[n6822].hii_galaxy, data[n6822].galaxy_rc3_pa, data[n6822].galaxy_twomass_pa, $
;     data[n6822].galaxy_rc3_incl, data[n6822].galaxy_twomass_incl
    if (nn6822 ne 0L) then begin
       splog, 'Setting NGC 6822 inclination and position angles!'
       data[n6822].galaxy_rc3_incl = 50.1
       data[n6822].galaxy_twomass_incl = 50.1
       data[n6822].galaxy_rc3_pa = 122.0
       data[n6822].galaxy_twomass_pa = 122.0
    endif
    
; needed quantitites
    need = where((data.galaxy_rc3_pa lt -900.0) and $
      (data.galaxy_twomass_pa gt -900.0),nneed)
    if (nneed ne 0L) then begin
;      splog, 'Replacing RC3 with 2MASS position angles for the following objects'
;      niceprint, data[need].hii_galaxy, data[need].ned_galaxy, $
;        data[need].galaxy_rc3_pa, data[need].galaxy_twomass_pa
       data[need].galaxy_rc3_pa = data[need].galaxy_twomass_pa
    endif

    need = where((data.galaxy_rc3_incl lt -900.0) and $
      (data.galaxy_twomass_incl gt -900.0),nneed)
    if (nneed ne 0L) then begin
;      splog, 'Replacing RC3 with 2MASS inclination angles for the following objects'
;      niceprint, data[need].hii_galaxy, data[need].ned_galaxy, $
;        data[need].galaxy_rc3_incl, data[need].galaxy_twomass_incl
       data[need].galaxy_rc3_incl = data[need].galaxy_twomass_incl
    endif

; de-project using the RC3 and 2MASS angles       
    splog, 'Computing deprojected galactocentric radii'

    hii_region = strtrim(data.hii_galaxy,2)+'/'+strtrim(data.hii_region,2)

    data.hii_rc3_radius = im_hiiregion_deproject(data.galaxy_rc3_incl,$
      data.galaxy_rc3_pa,data.hii_raoffset,data.hii_deoffset,hii_phi=hii_phi,$
      hii_region=hii_region)
    data.hii_rc3_phi = hii_phi
    good = where((data.galaxy_rc3_r25 gt -900.0) and (data.hii_rc3_radius gt -900.0),ngood)
    if (ngood ne 0L) then data[good].hii_rc3_rr25 = $
      data[good].hii_rc3_radius/(data[good].galaxy_rc3_r25*60.0)

    data.hii_twomass_radius = im_hiiregion_deproject(data.galaxy_twomass_incl,$
      data.galaxy_twomass_pa,data.hii_raoffset,data.hii_deoffset,hii_phi=hii_phi,$
      hii_region=hii_region)
    data.hii_twomass_phi = hii_phi
    good = where((data.galaxy_twomass_k20 gt -900.0) and (data.hii_twomass_radius gt -900.0),ngood)
    if (ngood ne 0L) then data[good].hii_twomass_rr25 = $
      data[good].hii_twomass_radius/(data[good].galaxy_twomass_k20*60.0)

; fill the "nolog" structure       
    data_nolog.ned_galaxy          = data.ned_galaxy
    data_nolog.galaxy_ra           = data.galaxy_ra
    data_nolog.galaxy_dec          = data.galaxy_dec
    data_nolog.galaxy_rc3_pa       = data.galaxy_rc3_pa
    data_nolog.galaxy_twomass_pa   = data.galaxy_twomass_pa
    data_nolog.galaxy_rc3_incl     = data.galaxy_rc3_incl
    data_nolog.galaxy_twomass_incl = data.galaxy_twomass_incl
    data_nolog.galaxy_rc3_r25      = data.galaxy_rc3_r25
    data_nolog.galaxy_twomass_k20  = data.galaxy_twomass_k20

    data_nolog.hii_rc3_radius      = data.hii_rc3_radius
    data_nolog.hii_rc3_phi         = data.hii_rc3_phi
    data_nolog.hii_rc3_rr25        = data.hii_rc3_rr25
    data_nolog.hii_twomass_radius  = data.hii_twomass_radius
    data_nolog.hii_twomass_phi     = data.hii_twomass_phi
    data_nolog.hii_twomass_rr25    = data.hii_twomass_rr25

; ---------------------------------------------------------------------------    
; compute electron temperatures; possibly apply a nominal S/N>1 cut 
    splog, 'Computing electron temperatures'
    t0 = systime(1)
    te = im_compute_te(linefit,snrcut=0.0,nmonte=nmonte)
    splog, format='("Total time = ",G0," minutes.")', $
      (systime(1)-t0)/60.0

    data = struct_addtags(temporary(data),te)
    data_nolog = struct_addtags(temporary(data_nolog),te)

; ---------------------------------------------------------------------------
; compute electron-temperature abundances; ToDo: some objects (from,
; e.g., the Kniazev et al. 2004; Izotov et al. 2006 catalogs of
; metal-poor galaxies in the SDSS) do not have a measure of [OII]
; 3727, but do have a measure of [OII] 7325; compute the total ionic
; abundances of those guys here:
   splog, 'Computing Te abundances'

   t0 = systime(1)
   oh12te = im_direct_abundance(linefit,oii_ion='3727',oiii_ion='5007',$
     oii_temp=data.zt_toii,oiii_temp=data.zt_toiii,err_oii_temp=data.zt_toii_err,$
     err_oiii_temp=data.zt_toiii_err,dens=dens,err_dens=err_dens,nmonte=nmonte)
   splog, format='("Total time = ",G0," minutes.")', (systime(1)-t0)/60.0

   data = struct_addtags(temporary(data),oh12Te)
   data_nolog = struct_addtags(temporary(data_nolog),oh12Te)

; ---------------------------------------------------------------------------    
; compute strong-line abundances; any S/N cut other than 0.0 screws up
; NGC2541 

   splog, 'Computing strong-line abundances'
   strong = im_abundance(linefit,snrcut=0.0,$ 
     nmonte=nmonte,/justflux) 
   strong = struct_trimtags(strong,except='*ZSTRONG*EW*')

   data = struct_addtags(temporary(data),strong)
   data_nolog = struct_addtags(temporary(data_nolog),strong)

; ---------------------------------------------------------------------------
; write the final data structure in the zeroth extension and the
; LINEFIT structure in the first extension
   
   splog, 'Writing '+outfile
   mwrfits, data, outfile, /create
   mwrfits, linefit, outfile
   spawn, 'gzip -f '+outfile, /sh
   
   splog, 'Writing '+outfile_nolog
   mwrfits, data_nolog, outfile_nolog, /create
   mwrfits, linefit, outfile_nolog
   spawn, 'gzip -f '+outfile_nolog, /sh

return
end    
