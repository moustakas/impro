;+
; NAME:
;   BUILD_GALEX_TILELIST
;
; PURPOSE:
;   Produce summary data structures for all the available gr4 and gr5
;   tiles. 
;
; INPUTS: 
;
; OPTIONAL INPUTS: 
;   tilefile - output file name, which will be needed by
;     BUILD_GALEX_CATALOG; default:
;     $CATALOGS_DIR/galex/galex_tiles_gr4_gr5.fits
;
; KEYWORD PARAMETERS: 
;
; OUTPUTS: 
;   Data structure with some basic information (ra, dec, etc.) for
;   each tile in GR4 and GR5.
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   The code assumes $GALEX_DIR is set to the location of the galex
;   data products.  Also, note that the code is not particularly
;   fast. 
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2010 Apr 21, UCSD - based on original code by
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

pro build_galex_tilelist, tilefile=tilefile

    galex_dir = getenv('GALEX_DIR')
    if (strtrim(galex_dir,2) eq '') then message, $
      '$GALEX_DIR environment variable required'

; default output filename    
    if (n_elements(tilefile) eq 0) then tilefile = $
      getenv('CATALOGS_DIR')+'/galex/galex_tiles_gr4_gr5.fits'

; build the list of MCATS files for each GR release and both versions;
; 02-vsn contains data from the AIS, while the 01-vsn directory seems
; to include some of the deeper imaging
    mcats1 = file_search(galex_dir+'/gr4/pipe/01-vsn/*/d/01-main/*/*/*mcat*fits*')
    mcats2 = file_search(galex_dir+'/gr4/pipe/02-vsn/*/d/01-main/*/*/*mcat*fits*')
    mcats3 = file_search(galex_dir+'/gr5/pipe/01-vsn/*/d/01-main/*/*/*mcat*fits*')
    mcats4 = file_search(galex_dir+'/gr5/pipe/02-vsn/*/d/01-main/*/*/*mcat*fits*')
    help, mcats1, mcats2, mcats3, mcats4
    mcats = [mcats1,mcats2,mcats3]
;   mcats = [mcats1,mcats2,mcats3,mcats4]
    nfiles = n_elements(mcats)

    out = {gr: '', vsn: 0, tile_num: 0L, tile_name: '', $
      tile_ra: 0.0D, tile_dec: 0.0D}
    out = replicate(out,nfiles)

;   for ifile = 5000L, 5100L do begin
    for ifile = 0L, nfiles-1L do begin
;   for ifile = 0L, 10 do begin
       print, "Reading file ", ifile, " of ", nfiles, string(13B), $
         format='(A15, I6, A6, I6, A1, $)'
       mcats0 = mcats[ifile]
       data = mrdfits(mcats0,1,/silent)

       out[ifile].gr = strmid(repstr(mcats0,galex_dir,''),1,3)
       out[ifile].vsn = data[0].vsn
       out[ifile].tile_name = strtrim((strsplit(file_basename(mcats0),'-',/extract))[0],2)
       out[ifile].tile_num = data[0].tilenum
       out[ifile].tile_ra = djs_median(data.alpha_j2000)
       out[ifile].tile_dec = djs_median(data.delta_j2000)
    endfor

    im_mwrfits, out, tilefile, /clobber
    
return
end
