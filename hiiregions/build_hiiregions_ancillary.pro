;+
; NAME:
;   BUILD_HIIREGIONS_ANCILLARY
;
; PURPOSE:
;   Retrieve ancillary data (basic NED data and diameters) for the
;   HII regions database. 
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2007 Dec 17, NYU
;   jm10feb22ucsd - modified to use the galaxy list output from
;     BUILD_HIIREGIONS_GALAXY_LIST 
;-

pro build_hiiregions_ancillary, basic=basic, diam=diam

    hiipath = hiiregions_path()
    version = hiiregions_version()

    listfile = hiipath+'hiiregions_galaxy_list_'+version+'.txt'
    readcol, listfile, galaxy, format='A', comment='#', /silent
    galaxy = strtrim(galaxy,2)
    galaxy = galaxy[uniq(galaxy,sort(galaxy))]

    galaxy = [$
      'KDG061',$
      'NGC7552',$
      'UGC04305',$
      'UGC04459',$
      'UGC05139',$
      'UGC05336',$
      'UGC05666',$
      'UGC05692',$
      'UGC05918',$
      'UGC08201']
    
; basic data    
    if keyword_set(basic) then begin
       ned_webget_basic, galaxy, basic
       fitsfile = 'hii_region_ned_basic_'+version+'.fits'
       fitsfile = 'junk_basic.fits'
       im_mwrfits, basic, hiipath+fitsfile
    endif

; diameters
    if keyword_set(diam) then begin
       ned_webget_diameters, galaxy, diam
       fitsfile = 'hii_region_ned_diameters_'+version+'.fits'
       fitsfile = 'junk_diameters.fits'
       im_mwrfits, diam, hiipath+fitsfile
    endif

return
end
