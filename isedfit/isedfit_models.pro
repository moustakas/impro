;+
; THE DOCUMENTATION IS OUT OF DATE!
;
; NAME:
;   ISEDFIT_MODELS
;
; PURPOSE:
;   Compute model photometry on a grid of redshift, star-formation
;   history, age, and stellar metallicity.
;
; INPUTS:
;   filterlist  - compute model photometry for this list of
;                 filters; the filters must be loaded properly into
;                 the KCORRECT database, otherwise a disaster
;                 beyond your imagination will occur [NFILTER]
;
; OPTIONAL INPUTS:
;   datapath     - I/O path [default to CWD()]
;   modelsprefix - prefix to prepend to the output data
;                  structures (default 'isedfit_models') 
;   minredshift  - compute model photometry beginning at redshift
;                  MINREDSHIFT (default 0.01)
;   maxredshift  - compute model photometry ending at MAXREDSHIFT
;                  (default 1.0) 
;   dredshift    - redshift spacing between MINREDSHIFT and
;                  MAXREDSHIFT (default 0.01)
;   minage       - minimum model age (default 0.1 Gyr)
;   maxage       - maximum model age; in practice, should be set to
;                  the age of the universe at the lowest redshift
;                  of the sample (default 13.5 Gyr) 
;   dlogage      - logarithmic interval between MINAGE and MAXAGE
;                  (default 0.1) 
;   Zstellar     - grid of stellar metallicities; this input must
;                  be an integer array ranging between 0 and
;                  5, inclusive; the default is [2,3,4,5],
;                  corresponding to stellar metallicities
;                  Z=[0.004,0.008,0.02,0.05] 
;   tau          - input tau values (stronger than MINTAU, MAXTAU,
;                  and DTAU; Gyr)
;   mintau       - minimum tau value (default 0.0 Gyr)
;   maxtau       - maximum tau value (default 20.0 Gyr)
;   dtau         - tau interval (default 1.0 Gyr)
;   minebv       - minimum reddening (default 0.0 mag)
;   maxebv       - maximum reddening (default 1.0 mag)
;   debv         - reddening interval (default 0.05 mag)
;
; KEYWORD PARAMETERS:
;   allages - use the ages inherent to BC03, but still bounded by
;             MINAGE and MAXAGE (recommended)
;   alltau  - use all the available TAU values between MINTAU and
;             MAXTAU (see WRITE_SFH_BASE_MODELS)
;   debug   - generate a simple QA plot showing the each model
;             and the broadband photometry overlaid; not very
;             useful, but fun to watch
;   write   - write out the results to DATAPATH, using the
;             specified MODELSPREFIX
;
; OUTPUTS:
;   isedfit      - output data structure containing the model
;                  photometry for every model [NFILTER, NTAU,
;                  NAGE, NZ, NREDSHIFT]
;   isedfit_info - information structure specifying the parameters 
;                  of the models
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   This routine draws on a subset of models pre-computed by
;   WRITE_SFH_BASE_MODELS. 
;
;   The initial mass function has been hard-wired to Salpeter
;   evaluated between 0.1 and 100 solar masses.  
;
;   The ages of the models will not be at the exact age grid
;   specified by MINAGE, MAXAGE, and DLOGAGE since the BC03
;   SSP's have been computed only at pre-determined ages.  The
;   output ages are replaced by the actual BC03 ages.
; 
;   The tau values specified by TAU, or MINTAU and MAXTAU, must
;   have been pre-computed in WRITE_SFH_BASE_MODELS, otherwise the
;   routine exits with an error.  I suggest using tau-values that
;   are spaced quasi-logarithmically, since the SFH is 
;   characterized by the ratio AGE/TAU.
;
;   The cosmology has been hard-wired to h=0.7, Omega_0=0.3, and
;   Omega_Lambda=0.7, although in principal these can be
;   optional. 
;
;   IGM absorption is properly accounted for.
;
;   Model photometry is also synthesized for a solar metallicity,
;   SSP with an age equal to the age of the universe at each value
;   of REDSHIFT.  This model is used in the two-component fits
;   implemented in ISEDFIT_CHI2.
;
; EXAMPLES:
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2005 Feb 10-12, U of A - written
;   jm05mar22uofa - cleaned up documentation; added IGM
;                   attenuation; various extinction curves (Milky
;                   Way, SMC, Calzetti, or Charlot & Fall) now
;                   possible; IGM attenuation incorporated 
;   jm06feb21uofa - updates incorporated to be compatible with 
;                   WRITE_SFH_BASE_MODELS; variables containing
;                   "TAGE" now changed to "AGE"; obsolete MINWAVE
;                   and MAXWAVE parameters removed
;   jm06mar01uofa - excised from ISEDFIT_MODELS; dust is now
;                   treated as a free parameter in ISEDFIT; also,
;                   the burst model photometry is computed along
;                   with the tau-model photometry; NTAU input
;                   changed to DTAU; TAU optional input added 
;   jm06mar16uofa - synthesize photometry for a maximally old
;                   second component
;   jm06jul14uofa - added ALLAGES, ALLTAU keywords; use
;                   IM_FILTERMAG() instead of K_PROJECT_FILTERS() 
;   jm07jun20nyu  - speeded up by a factor of 10 or more, but now
;                   using Blanton's filter code; there's still a
;                   discrepancy between IM_FILTERMAG() and
;                   K_PROJECT_FILTERS()!
;
; Copyright (C) 2005-2007, John Moustakas
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

pro isedfit_models, paramfile, params=params, iopath=iopath, $
  sfhgrid_basedir=sfhgrid_basedir, clobber=clobber

    if (n_elements(paramfile) eq 0) and $
      (n_elements(params) eq 0) then begin
       doc_library, 'isedfit_models'
       return
    endif

; read the parameter file; parse to get the relevant path and
; filenames
    if (n_elements(iopath) eq 0) then iopath = './'
    if (n_elements(params) eq 0) then params = $
      read_isedfit_paramfile(paramfile)

; SFHGRID and REDCURVE can be vectors; however, if SFHGRID=3 (no
; reddening) then ignore REDCURVE    
    nsfhgrid = n_elements(params.sfhgrid)
    nredcurve = n_elements(params.redcurve)
    if (nsfhgrid gt 1) or (nredcurve gt 1) then begin
       for ii = 0, nsfhgrid-1 do begin
          newparams1 = struct_trimtags(params,except='sfhgrid')
          newparams1 = struct_addtags(newparams1,{sfhgrid: params.sfhgrid[ii]})
          if (newparams1.sfhgrid eq 3) then nredcurve = 1
          for jj = 0, nredcurve-1 do begin
             newparams2 = struct_trimtags(newparams1,except='redcurve')
             newparams2 = struct_addtags(newparams2,{redcurve: params.redcurve[jj]})
             isedfit_models, params=newparams2, iopath=iopath, $
               sfhgrid_basedir=sfhgrid_basedir, clobber=clobber
          endfor
       endfor 
       return 
    endif 

    fp = isedfit_filepaths(params,iopath=iopath,sfhgrid_basedir=sfhgrid_basedir)
    chunkfile = fp.modelspath+fp.isedfit_models_chunkfiles[0] ; check the first file
    if file_test(chunkfile+'.gz',/regular) and $
      (keyword_set(clobber) eq 0) then begin
       splog, 'First ChunkFile '+chunkfile+' exists; use /CLOBBER'
       return
    endif
    
    splog, 'SYNTHMODELS='+params.synthmodels+', '+$
      'REDCURVE='+strtrim(params.redcurve,2)+', IMF='+$
      params.imf+', '+'SFHGRID='+$
      string(params.sfhgrid,format='(I2.2)')

; filters and redshift grid
    filterlist = strtrim(params.filterlist,2)
    filtinfo = im_filterspecs(filterlist=filterlist)
    nfilt = n_elements(filterlist)

    splog, 'Synthesizing photometry in '+$
      string(nfilt,format='(I0)')+' bandpasses:'
    niceprint, replicate('  ',nfilt), filterlist

    redshift = params.redshift
    nredshift = n_elements(redshift)
    splog, 'Redshift grid: '
    splog, '  Nz   = '+string(nredshift,format='(I0)')
    splog, '  zmin = '+string(min(redshift),format='(G0)')
    splog, '  zmax = '+string(max(redshift),format='(G0)')

    if (min(redshift) le 0.0) then begin
       splog, 'REDSHIFT should be positive and non-zero'
       return
    endif
    
    dist = 10.0*3.085678D18 ; fiducial distance [10 pc in cm]
    dlum = dluminosity(redshift,/cm) ; luminosity distance [cm]

; IGM attenuation    
    if params.igm then begin
       splog, 'Reading IGM attenuation lookup table'
       igmgrid = mrdfits(getenv('IMPRO_DIR')+'/dust/igmtau_grid.fits.gz',1)
    endif else begin
       splog, 'Neglecting IGM absorption'
    endelse 

; now loop on each "chunk" of models and build modelmaggies
    nchunk = n_elements(fp.sfhgrid_chunkfiles)
    splog, 'NCHUNK = '+string(nchunk,format='(I0)')
    t1 = systime(1)
    for ichunk = 0L, nchunk-1L do begin
       t0 = systime(1)
       splog, 'Reading '+fp.sfhgrid_chunkfiles[ichunk]
       chunk = mrdfits(fp.sfhgrid_chunkfiles[ichunk],1,/silent)
       ndim = size(chunk.flux,/n_dim)
       sz = size(chunk.flux,/dim)
       npix = sz[0] & nage = sz[1]
       if (ndim eq 2L) then nmodel = 1L else nmodel = sz[2]
       distfactor = rebin(reform((dist/dlum)^2.0,nredshift,1),nredshift,nmodel)
; initialize the output structure
       isedfit_models = struct_trimtags(chunk,except=['WAVE','FLUX'])
       isedfit_models = struct_addtags(temporary(isedfit_models),$
         replicate({modelmaggies: fltarr(nfilt,nage,nredshift)},nmodel))
; build the IGM absorption vector
       if (params.igm eq 1) then begin
          igm = fltarr(npix,nredshift)
          for iz = 0L, nredshift-1L do begin
             zwave = chunk[0].wave*(1.0+redshift[iz])
             windx = findex(igmgrid.wave,zwave)
             zindx = findex(igmgrid.zgrid,redshift[iz])
             igm[*,iz] = interpolate(igmgrid.igm,windx,zindx,/grid,missing=1.0)
;            plot, zwave, igm[*,iz], xr=[0,9000], xsty=3, ysty=3
;            get_element, igmgrid.zgrid, redshift[iz], jj                                   
;            plot, igmgrid.wave, igmgrid.igm[*,jj], ysty=3
          endfor
       endif
       nwave = n_elements(chunk[0].wave)
       for iage = 0L, nage-1L do begin
          if ((iage mod 5) eq 0) then print, format='("Chunk ",I0,"/",I0,", '+$
            'Age ",I0,"/",I0,A10,$)', ichunk+1L, nchunk, $
            iage+1L, nage, string(13b)
; jm09may16nyu - added IGM attenuation
; project the filters onto each SED, including IGM attenuation, and
; scale by the distance ratio
          flux = chunk.flux[*,iage]
; fast code, but doesn't include IGM attenuation
          if (params.igm eq 0) then begin
             k_projection_table, rmatrix, flux, k_lambda_to_edges(chunk[0].wave), $
               redshift, filterlist, /silent
             for ff = 0L, nfilt-1L do rmatrix[*,*,ff] = rmatrix[*,*,ff]*distfactor
             isedfit_models.modelmaggies[*,iage,*] = $
               reform(transpose(rmatrix,[2,0,1]),nfilt,1,nredshift,nmodel)
          endif else begin
; reasonably fast code that includes IGM absorption
;            for iz = 207, 207 do begin
             for iz = 0L, nredshift-1 do begin
                bigigm = rebin(igm[*,iz],nwave,nmodel)
                k_projection_table, rmatrix1, flux*bigigm, $
                  k_lambda_to_edges(chunk[0].wave), $
                  redshift[iz], filterlist, /silent
                isedfit_models.modelmaggies[*,iage,iz] = $
                  reform(transpose(rmatrix1*(dist/dlum[iz])^2.0,[2,0,1]))
                endfor
          endelse
;; very slow!
;         for iz = 0L, nredshift-1 do begin
;            for it = 0L, ngood-1L do begin
;               zwave = chunk[0].wave*(1.0+redshift[iz])
;               igm = interpolate(igmgrid.igm,findex(igmgrid.wave,zwave),$
;                 findex(igmgrid.zgrid,redshift[iz]),/grid)
;               isedfit_models[good[it]].modelmaggies[*,iage,iz] = $
;                 k_project_filters(k_lambda_to_edges(zwave),$
;                 igm*flux[*,it]/(1.0+redshift[iz]),filterlist=filterlist,$
;                 /silent)*(dist/dlum[iz])^2.0
;            endfor
;         endfor
; test that I'm doing the distance ratio right
;         zindx = 10 & mindx = 5
;         ff = flux[*,mindx]/(1.0+redshift[zindx])*distfactor[zindx]
;         ww = k_lambda_to_edges(chunk[0].wave)*(1.0+redshift[zindx])
;         mm = k_project_filters(ww,ff,filterlist=filterlist)
;         niceprint, mm, rmatrix[zindx,mindx,*]
       endfor 
       if (ichunk eq 0L) then splog, format='("Total time for CHUNK-001 '+$
         '= ",G0," minutes")', (systime(1)-t0)/60.0

       outfile = fp.modelspath+fp.isedfit_models_chunkfiles[ichunk]
       splog, 'Writing '+outfile
       mwrfits, isedfit_models, outfile, /create
       spawn, 'gzip -f '+outfile, /sh
    endfor
    splog, format='("Total time for all models '+$
         '= ",G0," minutes")', (systime(1)-t1)/60.0    

return
end
