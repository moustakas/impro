function init_scamp_config, ncopies
; jm08aug06nyu - support routine for IM_SCAMP

    if (n_elements(ncopies) eq 0L) then ncopies = 1L

    config = {$
      FGROUP_RADIUS:          1.0             , $
      REF_SERVER:       'cocat1.u-strasbg.fr' , $
      REF_PORT:               1660            , $
      CDSCLIENT_EXEC:       'aclient'         , $
      ASTREF_CATALOG:       'USNO-B1'         , $
      ASTREF_BAND:          'DEFAULT'         , $
      ASTREFCAT_NAME:       'astrefcat.cat'   , $
      ASTREFCENT_KEYS:      'X_WORLD,Y_WORLD' , $
      ASTREFERR_KEYS:       'ERRA_WORLD, ERRB_WORLD, ERRTHETA_WORLD', $
      ASTREFMAG_KEY:        'MAG'             , $
      SAVE_REFCATALOG:      'N'               , $
      REFOUT_CATPATH:       '.'               , $
      MERGEDOUTCAT_TYPE:    'NONE'            , $
      MATCH:                  'Y'             , $
      MATCH_NMAX:             0               , $
      PIXSCALE_MAXERR:        1.2             , $
      POSANGLE_MAXERR:        5.0             , $
      POSITION_MAXERR:        1.0             , $
      MATCH_RESOL:            0               , $
      MATCH_FLIPPED:          'N'               , $
      MOSAIC_TYPE:            'UNCHANGED'       , $
      FIXFOCALPLANE_NMIN:     1               , $
      CROSSID_RADIUS:         2.0             , $
      SOLVE_ASTROM:           'Y'               , $
      ASTRINSTRU_KEY:         'FILTER,QRUNID'   , $
      STABILITY_TYPE:         'INSTRUMENT'      , $
      CENTROID_KEYS:          'XWIN_IMAGE,YWIN_IMAGE' , $
      CENTROIDERR_KEYS:       'ERRAWIN_IMAGE,ERRBWIN_IMAGE,ERRTHETAWIN_IMAGE', $
      DISTORT_KEYS:           'XWIN_IMAGE,YWIN_IMAGE' , $
      DISTORT_GROUPS:         '1,1'           , $
      DISTORT_DEGREES:        '3'             , $ ; should be string to allow comma-separated values
      ASTREF_WEIGHT:          1.0             , $
      ASTRCLIP_NSIGMA:        3.0             , $
      CORRECT_COLOURSHIFTS:   'N'             , $
      SOLVE_PHOTOM:           'Y'             , $
      MAGZERO_OUT:            0.0             , $
      MAGZERO_INTERR:         0.01            , $
      MAGZERO_REFERR:         0.03            , $
      PHOTINSTRU_KEY:         'FILTER'        , $
      MAGZERO_KEY:            'PHOT_C'        , $
      EXPOTIME_KEY:           'EXPTIME'       , $
      AIRMASS_KEY:            'AIRMASS'       , $
      EXTINCT_KEY:            'PHOT_K'        , $
      PHOTOMFLAG_KEY:         'PHOTFLAG'      , $
      PHOTFLUX_KEY:           'FLUX_AUTO'     , $
      PHOTFLUXERR_KEY:        'FLUXERR_AUTO'  , $
      PHOTCLIP_NSIGMA:        3.0             , $
      CHECKPLOT_CKEY:         'SCAMPCOL'      , $
      CHECKPLOT_DEV:          'PNG'           , $
      CHECKPLOT_RES:          '1000,1000'     , $ ; originally 0 (i.e., 800,600)
      CHECKPLOT_ANTIALIAS:    'Y'             , $
      CHECKPLOT_TYPE:         'FGROUPS,DISTORTION,ASTR_INTERROR2D,ASTR_INTERROR1D,ASTR_REFERROR2D,ASTR_REFERROR1D,ASTR_CHI2,PHOT_ERROR', $
      CHECKPLOT_NAME:         'fgroups,distort,astr_interror2d,astr_interror1d,astr_referror2d,astr_referror1d,astr_chi2,psphot_error' , $
      CHECKIMAGE_TYPE:        'NONE'          , $
      CHECKIMAGE_NAME:        'check.fits'    , $
      SN_THRESHOLDS:          '10.0,100.0'    , $
      FWHM_THRESHOLDS:        '0.0,100.0'     , $
      FLAGS_MASK:             '0x00f0'        , $
      WEIGHTFLAGS_MASK:       '0x00ff'        , $
      IMAFLAGS_MASK:          '0x0'           , $
      AHEADER_GLOBAL:         'scamp.ahead'   , $
      AHEADER_SUFFIX:         '.ahead'        , $
      HEADER_SUFFIX:          '.head'         , $
      HEADER_TYPE:            'NORMAL'        , $
      VERBOSE_TYPE:           'NORMAL'        , $
      WRITE_XML:              'Y'             , $
      XML_NAME:               'scamp.xml'     , $
      XSL_URL:                'file://scamp.xsl' , $
      NTHREADS:               0               }
    config = replicate(config,ncopies)

return, config
end
