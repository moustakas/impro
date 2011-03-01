;+
; NAME:
;   BUILD_GALEX_CATALOG
;
; PURPOSE:
;   Given an input catalog, build the corresponding line-matched GALEX
;   catalog, including upper limits.
;
; INPUTS: 
;   ra, dec - coordinates of the galaxies in the survey of interest
;     (decimal degrees)
;
; OPTIONAL INPUTS: 
;   tilefile - name of the tile information structure written out by
;     BUILD_GALEX_TILELIST
;   outfile - output file name
;   galex_snrcut - minimum S/N for a GALEX source to be matched
;     (default 5)
;
; KEYWORD PARAMETERS: 
;   clobber - overwrite OUTFILE, if it exists
;   debug - render some simple debugging plots (and wait for a
;     keystroke) 
;
; OUTPUTS: 
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   Currently hard-coded to use the gr4/gr5 data releases.  The upper
;   limits bit of the code currently does not work - it's too slow! 
;
; EXAMPLES:
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2010 Feb 23, UCSD - largely based on code written by
;     R. Cool
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

function resolve_galex_duplicates, ingalex
; identify duplicates given an input GALEX catalog; among repeat
; observations, choose the object with the high S/N in the NUV 

    nobj = n_elements(ingalex)
    outgalex = struct_addtags(ingalex,$
      replicate({isbest: 1, ndup: 0},nobj))

; run spheregroup and then loop through every unique group to check
; for repeats; note that given the 2" matching distance most of
; these repeats will be repeat observations (rather than overlaps from
; the same tile)
    group = spheregroup(ingalex.alpha_j2000,ingalex.delta_j2000,$
      2D/3600D,firstg=firstg,multg=multg,nextg=nextg)
    
    ugroup = group[uniq(group[sort(group)])]
    ngroup = n_elements(ugroup)
    
    for igroup = 0L, ngroup-1L do begin
       print, "Group ", igroup, " of ", ngroup, string(13b), $
         format='(A0,I0,A0,I0,A1,$)'
       these = where(group eq ugroup[igroup],nthese)
       if (nthese ge 2) then begin
          bestsnr = max(ingalex[these].nuv_s2n,ibest)
          outgalex[these].ndup = nthese ; number of duplicates
          outgalex[these].isbest = 0
          outgalex[these[ibest]].isbest = 1
       endif
     endfor
    
return, outgalex
end

function galex_intile, ra, dec, tilestr, objobserved=objobserved
; simple wrapper to identify the subset of tiles that overlap the
; input set of objects
    
    nobj = n_elements(ra)
    ntile = n_elements(tilestr)
    tileused = lonarr(ntile)
    objobserved = lonarr(nobj)

    for itile = 0L, ntile-1 do begin
       dist = djs_diff_angle(ra,dec,tilestr[itile].tile_ra,$
         tilestr[itile].tile_dec)
       kobs = where(dist le 0.6, nmatch)
       if nmatch gt 0 then objobserved[kobs] = 1
       if nmatch gt 0 then tileused[itile] = 1
    endfor

return, tileused
end

function select_galextags
; return the list of all pipeline GALEX tags we want to keep for our
; output catalogs

    galextags = [$
      'E_BV', $
      'NUV_S2N', $
      'FUV_S2N', $ 
      'NUV_EXPTIME', $
      'NUV_A_IMAGE', $
      'NUV_B_IMAGE', $
      'NUV_KRON_RADIUS', $
      'NUV_FLUX_AUTO', $
      'NUV_FLUXERR_AUTO', $ 
      'NUV_SKYBG', $
      'NUV_WEIGHT', $ 
      'NUV_MAG_AUTO', $ 
      'NUV_MAGERR_AUTO', $
      'NUV_ZPMAG', $
      'NUV_FLUX_APER_1', $
      'NUV_FLUX_APER_2', $
      'NUV_FLUX_APER_3', $
      'NUV_FLUX_APER_4', $
      'NUV_FLUX_APER_5', $
      'NUV_FLUX_APER_6', $
      'NUV_FLUX_APER_7', $
      'NUV_FLUXERR_APER_1', $ 
      'NUV_FLUXERR_APER_2', $ 
      'NUV_FLUXERR_APER_3', $ 
      'NUV_FLUXERR_APER_4', $ 
      'NUV_FLUXERR_APER_5', $ 
      'NUV_FLUXERR_APER_6', $ 
      'NUV_FLUXERR_APER_7', $ 
      'FUV_A_IMAGE', $
      'FUV_B_IMAGE', $
      'FUV_KRON_RADIUS', $
      'FUV_FLUX_AUTO', $
      'FUV_FLUXERR_AUTO', $ 
      'FUV_SKYBG', $
      'FUV_WEIGHT', $ 
      'FUV_MAG_AUTO', $ 
      'FUV_MAGERR_AUTO', $
      'FUV_ZPMAG', $
      'FUV_FLUX_APER_1', $
      'FUV_FLUX_APER_2', $
      'FUV_FLUX_APER_3', $
      'FUV_FLUX_APER_4', $
      'FUV_FLUX_APER_5', $
      'FUV_FLUX_APER_6', $
      'FUV_FLUX_APER_7', $
      'FUV_FLUXERR_APER_1', $ 
      'FUV_FLUXERR_APER_2', $ 
      'FUV_FLUXERR_APER_3', $ 
      'FUV_FLUXERR_APER_4', $ 
      'FUV_FLUXERR_APER_5', $ 
      'FUV_FLUXERR_APER_6', $ 
      'FUV_FLUXERR_APER_7', $ 
      'alpha_j2000', $
      'delta_j2000']
    
return, galextags
end

pro build_galex_catalog, ra, dec, tilefile=tilefile, outfile=outfile, $
  galex_snrcut=galex_snrcut, out_galex=out_galex, clobber=clobber, $
  debug=debug, nowrite=nowrite
    
    nobj = n_elements(ra)
    if (nobj eq 0) or (n_elements(dec) eq 0) then begin
       doc_library, 'build_galex_catalog'
       return
    endif
    if (nobj ne n_elements(dec)) then begin
       splog, 'Dimensions of RA and DEC do not match'
       return
    endif

; output file name    
    if (n_elements(outfile) eq 0) then outfile = 'galex_matched.fits'
    if file_test(outfile+'.gz') and (keyword_set(clobber) eq 0) then begin
       splog, 'Output file '+outfile+' exists; use /CLOBBER'
       return
    endif

    if (n_elements(galex_snrcut) eq 0) then galex_snrcut = 5.0

; check for the tilefile written out by BUILD_GALEX_TILELIST
    if (n_elements(tilefile) eq 0) then tilefile = $
      getenv('CATALOGS_DIR')+'/galex/galex_tiles_gr4_gr5.fits.gz'
    if (file_test(tilefile) eq 0) then begin
       splog, 'Tilefile structure '+tilefile+' not found'
       return
    endif
    splog, 'Reading '+tilefile
    tiles = mrdfits(tilefile,1)    

; only consider the GALEX tiles that overlap the input survey
    splog, 'Finding matching tiles'
    t0 = systime(1)
    tileused = galex_intile(ra,dec,tiles,objobserved=objobserved)
    splog, 'Total time = ', (systime(1)-t0)/60.0

    if keyword_set(debug) then begin
       these = where(tileused)
       djs_plot, ra, dec, ps=3, xsty=3, ysty=3
       tvcircle, 0.6, tiles.tile_ra, tiles.tile_dec, $
         color=djs_icolor('red'), /data
       tvcircle, 0.6, tiles[these].tile_ra, tiles[these].tile_dec, $
         color=djs_icolor('green'), thick=3, /data
       cc = get_kbrd(1)
    endif

    these_tiles = where(tileused,nthese_tiles)
    if (nthese_tiles eq 0L) then begin
       splog, 'No GALEX photometry available'
       return
    endif
    splog, 'Found '+strtrim(nthese_tiles,2)+' matching tiles'
    struct_print, tiles[these_tiles]
    
; build a merged catalog of all the GALEX observations we care about    
    splog, 'Building merged GALEX catalog'
    galextags = select_galextags() ; tags we care about
    delvarx, galex
    t0 = systime(1)
    for itile = 0L, nthese_tiles-1L do begin
       print, "Tile ", itile, " of ", nthese_tiles, string(13b), $
         format='(A0,I0,A0,I0,A1,$)'
       tilename = galex_tilename(tiles[these_tiles[itile]])
       galex1 = mrdfits(tilename,1,/silent,columns=galextags)
; apply S/N cut
       highsnr = where(galex1.nuv_s2n gt galex_snrcut,ngalex1)
       if (ngalex1 ne 0L) then begin
          galex1 = galex1[highsnr]
          galex1 = struct_addtags(replicate(tiles[itile],n_elements(galex1)),galex1)
          if (n_elements(galex) eq 0L) then galex = galex1 else $
            galex = [temporary(galex),temporary(galex1)]
       endif
    endfor
    splog, 'Total time = ', (systime(1)-t0)/60.0

; find and resolve duplicates
    splog, 'Resolving duplicates'
    t0 = systime(1)
    ugalex = resolve_galex_duplicates(galex)
    best = where(ugalex.isbest,ncomp=ndup)
    splog, 'Removing '+string(ndup,format='(I0)')+' inferior observations'
    ugalex = galex[best]
    splog, 'Total time = ', (systime(1)-t0)/60.0

;; remove low-S/N observations
;    highsnr = where(ugalex.nuv_s2n gt galex_snrcut,ngalex,ncomp=ncrap)
;    splog, 'Removing '+string(ncrap,format='(I0)')+' observations with '+$
;      'S/N<'+strtrim(string(galex_snrcut,format='(F12.1)'),2)
;    splog, 'Final GALEX catalog contains '+$
;      string(ngalex,format='(I0)')+' objects'
;    galex = ugalex[highsnr]

; finally, spherematch!  use a 2" search radius
    splog, 'Sphere-matching'
    t0 = systime(1)
    spherematch, ra, dec, galex.alpha_j2000, $
      galex.delta_j2000, 2D/3600D, m1, m2, dist
    nmatch = long(total(m1 ne -1))
    splog, 'Found '+string(nmatch,format='(I0)')+' matching sources'
    out_galex = im_empty_structure(galex,$
      empty_value=-999.0,ncopies=nobj)
    if (m1[0] ne -1) then out_galex[m1] = galex[m2]

; add tags of convenience and relating to the upper limit
; calculations, below
    moretags = replicate({ra: -999.0D, dec: -999.0D, $
      galex_detect: 0, galex_limit: 0, galex_insurvey: 0, $
      nuv_flux_aper3_limit: 0.0, fuv_flux_aper3_limit: 0.0, $
      nuv_flux_aper4_limit: 0.0, fuv_flux_aper4_limit: 0.0, $
      nuv_flux_aper5_limit: 0.0, fuv_flux_aper5_limit: 0.0},nobj)
    out_galex = struct_addtags(moretags,temporary(out_galex))
    out_galex.ra = ra
    out_galex.dec = dec
    if (m1[0] ne -1) then out_galex[m1].galex_detect = 1
    splog, 'Total time = ', (systime(1)-t0)/60.0

; flag objects that are not in the GALEX footprint (and therefore
; would never had a good match)    
    splog, 'Constructing the GALEX window function'
    t0 = systime(1)
    radius = 0.6                ; GALEX FOV [degrees]
    circles = construct_polygon(ncaps=1,nelem=nthese_tiles)
    for ii = 0L, nthese_tiles-1L do begin 
       circles[ii].caps = ptr_new(circle_cap(ra=tiles[these_tiles[ii]].tile_ra,$
         dec=tiles[these_tiles[ii]].tile_dec,radius))
       circles[ii].ncaps = 1 
       circles[ii].use_caps = 1L
       circles[ii].str = garea(circles[ii]) 
       circles[ii].weight = 1.0
    endfor
    
    nodetect = where((out_galex.galex_detect eq 0),nlimit,$
      comp=detect,ncomp=ndetect)
    if (ndetect ne 0L) then out_galex[detect].galex_insurvey = 1
    if (nlimit ne 0L) then begin
       inwindow = is_in_window(ra=out_galex[nodetect].ra, $
         dec=out_galex[nodetect].dec,circles)
       insurvey = where(inwindow,ninsurvey)
       if (ninsurvey ne 0L) then $
         out_galex[nodetect[insurvey]].galex_insurvey = 1
    endif
    splog, 'Total time = ', (systime(1)-t0)/60.0

; the following code chunk is for computing upper limits, but it is
; ridiculously slow; needs to be rewritten!!    
    
;;;  Now find the objects that were not detect but that are in the
;;;  survey
;;    nodetect = where(out_galex.galex_detect eq 0 and $
;;      out_galex.galex_insurvey eq 1, nlimit)
;;          
;;;  Sadly, where we need to loop here
;;    limitradius = 120.0 ;; arcseconds
;;    nuv_limit = fltarr(nlimit)
;;    fuv_limit = fltarr(nlimit)
;;    npartner = lindgen(nlimit)*0
;;    
;;    print, "Performing duplicate spherematch"
;;    time0 = systime(/sec)
;;    spherematch, out_galex[nodetect].parent_ra, out_galex[nodetect].parent_dec, $
;;      galexcat_full.alpha_j2000, galexcat_full.delta_j2000, limitradius/3600.0, $
;;      m1, m2, dist, maxmatch=0
;;    print, "Time for spherematch is ", systime(/sec)-time0
;;    time0 = systime(/sec)
;;    print, "Calculating upper limits"
;;    for i = 0l, nlimit -1l do begin
;;       print, "Fluxing object ", i, " of ", nlimit, string(13b), $
;;         format='(a20, i6, a5, i6, a1, $)'
;;       keep = where(m1 eq i, nkeep)
;;       if nkeep ge 3 then begin
;;          subset = m2[keep]
;;          out_galex[nodetect[i]].nuv_flux_aper3_limit = median(galexcat_full[subset].nuv_fluxerr_aper_3)
;;          out_galex[nodetect[i]].fuv_flux_aper3_limit = median(galexcat_full[subset].fuv_fluxerr_aper_3)
;;          out_galex[nodetect[i]].nuv_flux_aper4_limit = median(galexcat_full[subset].nuv_fluxerr_aper_4)
;;          out_galex[nodetect[i]].fuv_flux_aper4_limit = median(galexcat_full[subset].fuv_fluxerr_aper_4)
;;          out_galex[nodetect[i]].nuv_flux_aper5_limit = median(galexcat_full[subset].nuv_fluxerr_aper_5)
;;          out_galex[nodetect[i]].fuv_flux_aper5_limit = median(galexcat_full[subset].fuv_fluxerr_aper_5)
;;       endif
;;    endfor
;;    out_galex[nodetect].galex_limit = 1
;;    print, "Time to calculate upper limits ", systime(/sec) - time0

; finally write out!    
    if (keyword_set(nowrite) eq 0) then $
      im_mwrfits, out_galex, outfile, clobber=clobber

return
end
