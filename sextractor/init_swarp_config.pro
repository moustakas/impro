function init_swarp_config, ncopies
; jm08aug06nyu - support routine for IM_SWARP

    if (n_elements(ncopies) eq 0L) then ncopies = 1L

    config = {$
      IMAGEOUT_NAME:          'coadd.fits'      , $
      WEIGHTOUT_NAME:       'coadd.weight.fits' , $
      HEADER_ONLY:            'N'               , $
      HEADER_SUFFIX:          '.head'           , $
      WEIGHT_TYPE:            'NONE'            , $
      WEIGHT_SUFFIX:          '.weight.fits'    , $
      WEIGHT_IMAGE:           ''                , $
      WEIGHT_THRESH:          ''                , $
      COMBINE:                'Y'               , $
      COMBINE_TYPE:           'MEDIAN'          , $
      BLANK_BADPIXELS:        'N'               , $
      CELESTIAL_TYPE:         'NATIVE'          , $
      PROJECTION_TYPE:        'TAN'             , $
      PROJECTION_ERR:         0.001             , $
      CENTER_TYPE:            'ALL'             , $
      CENTER:         '00:00:00.0, +00:00:00.0' , $
      PIXELSCALE_TYPE:        'MEDIAN'          , $
      PIXEL_SCALE:            0.0               , $
      IMAGE_SIZE:             '0'               , $ ; needs to be a string
      RESAMPLE:               'Y'               , $
      RESAMPLE_DIR:           '.'               , $
      RESAMPLE_SUFFIX:        '.resamp.fits'    , $
      RESAMPLING_TYPE:        'LANCZOS3'        , $
      OVERSAMPLING:           0                 , $
      INTERPOLATE:            'N'               , $
      FSCALASTRO_TYPE:        'FIXED'           , $
      FSCALE_KEYWORD:         'FLXSCALE'        , $
      FSCALE_DEFAULT:         1.0               , $
      GAIN_KEYWORD:           'GAIN'            , $
      GAIN_DEFAULT:           0.0               , $
      SUBTRACT_BACK:          'Y'               , $
      BACK_TYPE:              'AUTO'            , $
      BACK_DEFAULT:           0.0               , $
      BACK_SIZE:              128               , $
      BACK_FILTERSIZE:        3                 , $
      BACK_FILTTHRESH:        0.0               , $
      VMEM_DIR:               '.'               , $
      VMEM_MAX:               2047              , $
      MEM_MAX:                128               , $
      COMBINE_BUFSIZE:        64                , $
      DELETE_TMPFILES:        'Y'               , $
      COPY_KEYWORDS:          'OBJECT'          , $
      WRITE_FILEINFO:         'N'               , $
      WRITE_XML:              'N'               , $
      XML_NAME:               'swarp.xml'       , $
      XSL_URL:                'file://swarp.xsl', $
      VERBOSE_TYPE:           'NORMAL'          , $
      NNODES:                 1                 , $
      NODE_INDEX:             0                 , $
      NTHREADS:               0                 }
    config = replicate(config,ncopies)

return, config
end
