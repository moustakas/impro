;+
; NAME:
;   ISEDFIT_QUICK
;
; PURPOSE:
;   Quick start script that goes through the process.
;   Overwrites previous run results
;   This is just to test that you can run isedfit quickly.
;-

PRO isedfit_quick, extras=extras
  
  ;; output directory -- depending on parameters may be large [few Gb]
  prefix =  'example'
  
  isedfit_dir = '/tmp/isedfit_example/'
  montegrids_dir =  isedfit_dir+'montegrids/'
  file_mkdir, montegrids_dir

  isedfit_paramfile =  isedfit_dir+prefix+'_paramfile.par'

  ; --------------------------------------------------
  ;; catalog to be processed
  example = im_mrdfits('cosmos_example.fits.gz')
  filterlist = strtrim(example[0].filterlist, 2)

  
  ; --------------------------------------------------
  ; choose your priors: write the iSEDfit parameter file 
  write_isedfit_paramfile,  params=params,  isedfit_dir=isedfit_dir,  $
         prefix=prefix,  filterlist=filterlist,  spsmodels='fsps_v2.4_miles',  $
         imf='chab',  redcurve='charlot',  /igm,  zminmax=[0.1, 1.3],  zbin=0.05,  $
         nmodel=500L,  age=[0.1, 13.0],  tau=[0.1, 5.0],  Zmetal=[0.004, 0.03],  $
         AV=[0.35, 2.0],  mu=[0.1, 4.0],  pburst=0.2,  interval_pburst=2.0,  $
         tburst=[0.1, 13.0],  /delayed,  galchunksize=250L,  /clobber

  ; --------------------------------------------------
  ; build the Monte Carlo grids [time to go read astro-ph]
  isedfit_montegrids,  isedfit_paramfile,  isedfit_dir=isedfit_dir,  $
         montegrids_dir=montegrids_dir,  /clobber
  
  ; --------------------------------------------------
  ; calculate the model photometry 
  isedfit_models,  isedfit_paramfile,  isedfit_dir=isedfit_dir,  $
         montegrids_dir=montegrids_dir,  /clobber

  ; --------------------------------------------------
  ; fit the catalog
  isedfit,  isedfit_paramfile,  example.maggies,  example.ivarmaggies,  $
         example.z,  ra=example.ra,  dec=example.dec,  $
         isedfit_dir=isedfit_dir, /clobber, $
         isedfit_results=mass
  
  ;; simple plot
  plot, mass.mstar_50, mass.sfr100_50, psym=6, $
        xtitle='Mass [M_sun]', ytitle='SFR [M_sun / yr]'


  print, 'Huzza! iSEDfit finished.'
  print, ' Check the plot, and make sure that these sfr/masses make sense'
  

  ;; Here be dragons...

  IF keyword_set(extras) THEN BEGIN
    ;; Get some fancy plots [beta] that show the properties of the
    ;; models and the sort.
    
    isedfit_qaplot_models,  isedfit_paramfile,  example.maggies,  $
         example.ivarmaggies,  example.z,  isedfit_dir=isedfit_dir,  $
         thesefilters=filterlist, /clobber
    
    isedfit_kcorrect,  isedfit_paramfile,  isedfit_dir=isedfit_dir,  $
         montegrids_dir=montegrids_dir,  /clobber
    
    isedfit_qaplot_sed,  isedfit_paramfile,  nrandom=50,  $
         isedfit_dir=isedfit_dir,  montegrids_dir=montegrids_dir,  $
         /xlog, /clobber
    
  ENDIF



END
