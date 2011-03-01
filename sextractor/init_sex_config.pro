function init_sex_config, ncopies
; jm08aug06nyu - support routine for IM_SEX
; jm09mar05nyu - added default path to PARAMETERS, FILTER, and NNW
; files 

    if (n_elements(ncopies) eq 0L) then ncopies = 1L

    path = getenv('IMPRO_DIR')+'/sextractor/'
    config = {$
      CATALOG_NAME:     'test.cat'        ,$
      CATALOG_TYPE:     'ASCII_HEAD'      ,$
      PARAMETERS_NAME:  path+'default.sex.param'   ,$
      DETECT_TYPE:      'CCD'             ,$
      DETECT_MINAREA:   5                 ,$
      THRESH_TYPE:      'RELATIVE'        ,$
      DETECT_THRESH:    1.5               ,$
      ANALYSIS_THRESH:  1.5               ,$
      FILTER:           'Y'               ,$
      FILTER_NAME:      path+'default.conv'    ,$
      FILTER_THRESH:    0                 ,$
      DEBLEND_NTHRESH:  32                ,$
      DEBLEND_MINCONT:  0.005             ,$
      CLEAN:            'Y'               ,$
      CLEAN_PARAM:      1.0               ,$
      MASK_TYPE:        'CORRECT'         ,$
      WEIGHT_TYPE:      'NONE'            ,$
      WEIGHT_IMAGE:     'weight.fits'     ,$
      WEIGHT_GAIN:      'Y'               ,$
      WEIGHT_THRESH:    ''                ,$
      FLAG_IMAGE:       'flag.fits'       ,$
      FLAG_TYPE:        'OR'              ,$
      PHOT_APERTURES:   '5'               ,$ ; should be string
      PHOT_AUTOPARAMS:  '2.5,3.5'         ,$
      PHOT_PETROPARAMS: '2.0,3.5'         ,$
      PHOT_AUTOAPERS:   '0.0,0.0'         ,$
      PHOT_FLUXFRAC:    '0.5'             ,$ ; string
      SATUR_LEVEL:      50000.0           ,$
      MAG_ZEROPOINT:    0.0               ,$
      MAG_GAMMA:        4.0               ,$
      GAIN:             0.0               ,$
      PIXEL_SCALE:      1.0               ,$
      SEEING_FWHM:      1.2               ,$
      STARNNW_NAME:     path+'default.nnw'     ,$
      BACK_TYPE:        'AUTO'            ,$
      BACK_VALUE:       0.0               ,$
      BACK_SIZE:        64                ,$
      BACK_FILTERSIZE:  3                 ,$
      BACKPHOTO_TYPE:   'GLOBAL'          ,$
      BACKPHOTO_THICK:  24                ,$
      BACK_FILTTHRESH:  0.0               ,$
      CHECKIMAGE_TYPE:  'NONE'            ,$
      CHECKIMAGE_NAME:  'check.fits'      ,$
      MEMORY_OBJSTACK:  3000              ,$
      MEMORY_PIXSTACK:  300000            ,$
      MEMORY_BUFSIZE:   1024              ,$
      ASSOC_NAME:       'sky.list'        ,$
      ASSOC_DATA:       '2,3,4'           ,$
      ASSOC_PARAMS:     '2,3,4'           ,$
      ASSOC_RADIUS:     2.0               ,$
      ASSOC_TYPE:       'MAG_SUM'         ,$
      ASSOCSELEC_TYPE:  'MATCHED'         ,$
      VERBOSE_TYPE:     'NORMAL'          ,$
      WRITE_XML:        'N'               ,$
      XML_NAME:         'sex.xml'         ,$
      XSL_URL:          'file://sex.xsl'  , $
      NTHREADS:         1                 ,$
      FITS_UNSIGNED:    'N'               ,$
      INTERP_MAXXLAG:   16                ,$
      INTERP_MAXYLAG:   16                ,$
      INTERP_TYPE:      'ALL'             ,$
      PSF_NAME:         'default.psf'     ,$
      PSF_NMAX:         9                 ,$
      PSFDISPLAY_TYPE:  'SPLIT'           ,$
      SOM_NAME:         'default.som'     }
    config = replicate(config,ncopies)

return, config
end
