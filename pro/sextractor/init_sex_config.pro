;+
; NAME:
;   INIT_SEX_CONFIG()
; PURPOSE:
;   Initialize a sex configuration structure.
; OPTIONAL INPUTS: 
;   ncopies - number of times to copy the data structure
; OUTPUTS: 
;   config - sex configuration structure with all the data needed to
;     run sex using IM_SEX
; MODIFICATION HISTORY:
;   J. Moustakas, 2008 Aug 06, NYU
;   jm09mar05nyu - added default path to PARAMETERS, FILTER, and NNW
;     files  
;   jm11mar02ucsd - updated to SExtractor 2.8.6
;-

function init_sex_config, ncopies
; jm08aug06nyu - support routine for IM_SEX

    if (n_elements(ncopies) eq 0) then ncopies = 1

    etcpath = getenv('IMPRO_DIR')+'/etc/'
    config = {$
; catalog
      CATALOG_NAME:     'test.cat'        ,$
      CATALOG_TYPE:     'ASCII_HEAD'      ,$
      PARAMETERS_NAME:  etcpath+'default.sex.param'   ,$
; extraction
      DETECT_TYPE:      'CCD'             ,$
      DETECT_MINAREA:   5                 ,$
      THRESH_TYPE:      'RELATIVE'        ,$
      DETECT_THRESH:    1.5               ,$
      ANALYSIS_THRESH:  1.5               ,$
      FILTER:           'Y'               ,$
      FILTER_NAME:      etcpath+'default.conv'    ,$
      FILTER_THRESH:    0                 ,$
      DEBLEND_NTHRESH:  32                ,$
      DEBLEND_MINCONT:  0.005             ,$
      CLEAN:            'Y'               ,$
      CLEAN_PARAM:      1.0               ,$
      MASK_TYPE:        'CORRECT'         ,$
; weighting
      WEIGHT_TYPE:      'NONE'            ,$
      WEIGHT_IMAGE:     'weight.fits'     ,$
      WEIGHT_GAIN:      'Y'               ,$
      WEIGHT_THRESH:    ''                ,$
; flaging
      FLAG_IMAGE:       'flag.fits'       ,$
      FLAG_TYPE:        'OR'              ,$
; photometry
      PHOT_APERTURES:   '5'               ,$ ; should be string
      PHOT_AUTOPARAMS:  '2.5,3.5'         ,$
      PHOT_PETROPARAMS: '2.0,3.5'         ,$

      PHOT_AUTOAPERS:   '0.0,0.0'         ,$

      PHOT_FLUXFRAC:    '0.5'             ,$ ; string
      SATUR_LEVEL:      50000.0           ,$
;     SATUR_KEY:        'SATURATE'        ,$
      MAG_ZEROPOINT:    0.0               ,$
      MAG_GAMMA:        4.0               ,$
      GAIN:             '0.0'             ,$ ; should be string, to handle MEF
;     GAIN_KEY:         'GAIN'            ,$
      PIXEL_SCALE:      0.0               ,$
; star/galaxy separation
      SEEING_FWHM:      1.2               ,$
      STARNNW_NAME:     etcpath+'default.nnw'     ,$
; background
      BACK_TYPE:        'AUTO'            ,$
      BACK_VALUE:       0.0               ,$
      BACK_SIZE:        64                ,$
      BACK_FILTERSIZE:  3                 ,$
      BACKPHOTO_TYPE:   'GLOBAL'          ,$
      BACKPHOTO_THICK:  24                ,$
      BACK_FILTTHRESH:  0.0               ,$
; check image
      CHECKIMAGE_TYPE:  'NONE'            ,$
      CHECKIMAGE_NAME:  'check.fits'      ,$
; memory
      MEMORY_OBJSTACK:  3000              ,$
      MEMORY_PIXSTACK:  300000            ,$
      MEMORY_BUFSIZE:   1024              ,$
; association
      ASSOC_NAME:       'sky.list'        ,$
      ASSOC_DATA:       '2,3,4'           ,$
      ASSOC_PARAMS:     '2,3,4'           ,$
      ASSOC_RADIUS:     2.0               ,$
      ASSOC_TYPE:       'NEAREST'         ,$
      ASSOCSELEC_TYPE:  'MATCHED'         ,$
; miscellaneous
      VERBOSE_TYPE:     'NORMAL'          ,$
      WRITE_XML:        'N'               ,$
      XML_NAME:         'sex.xml'         ,$
      XSL_URL:          'file://sex.xsl'  , $
      NTHREADS:         0                 ,$
      FITS_UNSIGNED:    'N'               ,$
      INTERP_MAXXLAG:   16                ,$
      INTERP_MAXYLAG:   16                ,$
      INTERP_TYPE:      'ALL'             ,$
; experimental stuff
      PSF_NAME:         'default.psf'     ,$
      PSF_NMAX:         9                 ,$
      PSFDISPLAY_TYPE:  'SPLIT'           ,$
;     PATTERN_TYPE:     'RINGS-HARMONIC'  ,$
      SOM_NAME:         'default.som'     }
    config = replicate(config,ncopies)

return, config
end
