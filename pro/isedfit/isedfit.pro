;+
; NAME:
;   ISEDFIT
;
; PURPOSE:
;   Infer the physical properties of galaxies (probabalistically) by
;   modeling their observed spectral energy distributions.
;
; INPUTS:
;   isedfit_paramfile - iSEDfit parameter file
;   maggies - input galaxy photometry [NFILT,NGAL]
;   ivarmaggies - inverse variance array for MAGGIES [NFILT,NGAL]  
;   z - input galaxy redshifts [NGAL] 
;
; OPTIONAL INPUTS:
;   params - data structure with the same information contained in
;     ISEDFIT_PARAMFILE (over-rides ISEDFIT_PARAMFILE)
;   thissfhgrid - if ISEDFIT_PARAMFILE contains multiple grids then
;     build this SFHgrid (may be a vector)
;   isedfit_dir - full directory path where the results should be
;     written (default PWD=present working directory)  
;
;   outprefix - optionally write out files with a different prefix
;     from that specified in ISEDFIT_PARAMFILE (or PARAMS) (very
;     useful for fitting various subsets of the input sample with
;     different assumptions)
;   index - use this optional input to fit a zero-indexed subset of
;     the full sample (default is to fit everything)
;
;   ra,dec - galaxy right ascension and declination which will be
;     copied to the output structure for the user's convenience
;     [NGAL] (decimal degrees)
;
; KEYWORD PARAMETERS:
;   allages - allow solutions with ages that are older than the age of
;     the universe at the redshift of the object 
;   silent - suppress messages to STDOUT
;   nowrite - do not write out any of the output files (generally not
;     recommended but can be useful in certain situations) 
;   clobber - overwrite existing files of the same name (the default
;     is to check for existing files and if they exist to exit
;     gracefully)  
;
; OUTPUTS:
;   Binary FITS tables containing the fitting results are written
;   to ISEDFIT_DIR. 
;
; OPTIONAL OUTPUTS:
;   isedfit_results - output data structure containing all the
;     results; see the iSEDfit documentation for a detailed breakdown
;     and explanation of all the outputs
;   isedfit_post - output data structure containing the random draws
;     from the posterior distribution function, which can be used to
;     rebuild the posterior distributions of any of the output
;     parameters (using ISEDFIT_RECONSTRUCT_POSTERIOR) 
;   isedfit_outfile - output file name of the iSEDfit results
;
; COMMENTS:
;   Better documentation of the output data structures would be
;   helpful here.
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2011 Sep 01, UCSD - I began writing iSEDfit in 2005
;     while at the U of A, adding updates off-and-on through 2007;
;     however, the code has evolved so much that the old modification
;     history became obsolete!  Future changes to the officially
;     released code will be documented here.
;   jm13jan13siena - documentation rewritten and updated to reflect
;     many major changes
;   jm13aug09siena - updated to conform to a new and much simpler data
;     model; documentation updated
;
; Copyright (C) 2011, 2013, John Moustakas
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

function init_isedfit, ngal, nfilt, params=params, ra=ra, dec=dec, $
  isedfit_post=isedfit_post, photoz=photoz
; ISEDFIT support routine - initialize the output structure 

    ndraw = params.ndraw ; number of random draws
    if (params.nmaxburst eq 0) then burstarray1 = -1.0 else $
      burstarray1 = fltarr(params.nmaxburst)-1.0

    if keyword_set(photoz) then begin
       isedfit1 = {$
         isedfit_id:            -1L,$ ; unique ID number
         ra:                    -1D,$ ; RA [decimal degrees]
         dec:                   -1D,$ ; Dec [decimal degrees]
         maggies:     fltarr(nfilt),$ ; observed maggies
         ivarmaggies: fltarr(nfilt),$ ; corresponding inverse variance
         bestmaggies: fltarr(nfilt)}  ; best-fitting model photometry
    endif else begin
       isedfit1 = {$
         isedfit_id:            -1L,$ ; unique ID number
         ra:                    -1D,$ ; RA [decimal degrees]
         dec:                   -1D,$ ; Dec [decimal degrees]
         z:                    -1.0,$ ; redshift
         maggies:     fltarr(nfilt),$ ; observed maggies
         ivarmaggies: fltarr(nfilt),$ ; corresponding inverse variance
         bestmaggies: fltarr(nfilt)}  ; best-fitting model photometry
    endelse

; best-fit values (at the chi2 minimum); see ISEDFIT_MONTEGRIDS
    best = {$
      chunkindx:                 -1,$
      modelindx:                 -1,$
      delayed:       params.delayed,$ ; delayed tau model?
      bursttype:   params.bursttype,$ ; burst type
      chi2:                     1E6,$ ; chi2 minimum
      totalmass:               -1.0,$ 
      totalmass_err:           -1.0,$ 

      mstar:                   -1.0,$ 
      age:                     -1.0,$
      sfrage:                  -1.0,$
      tau:                     -1.0,$
      Zmetal:                  -1.0,$
      AV:                       0.0,$ ; initialize with zero to accommodate dust-free models!
      mu:                       1.0,$ ; always default to 1.0!
      oiiihb:                  -1.0,$ 
      nlyc:                    -1.0,$ 
      sfr:                     -1.0,$ ; instantaneous
      sfr100:                  -1.0,$ ; 100 Myr timescale
      b100:                    -1.0,$ ; birthrate parameter [100 Myr]
      b1000:                   -1.0,$ ; birthrate parameter [1 Gyr]
      ewoii:                   -1.0,$ 
      ewoiiihb:                -1.0,$ 
      ewniiha:                 -1.0,$ 
      nburst:                     0,$
      trunctau:                -1.0,$ ; burst truncation time scale
      tburst:           burstarray1,$
      dtburst:          burstarray1,$ ; 
      fburst:           burstarray1}  ; burst mass fraction

; median quantities and PDF quantiles
    qmed = {$
      mstar_50:    -1.0,$
      age_50:      -1.0,$
      sfrage_50:   -1.0,$
      tau_50:      -1.0,$
      Zmetal_50:   -1.0,$
      AV_50:       -1.0,$
      mu_50:       -1.0,$
      oiiihb_50:   -1.0,$
      sfr_50:      -1.0,$ ; instantaneous
      sfr100_50:   -1.0,$ ; 100 Myr
      b100_50:     -1.0,$
      b1000_50:    -1.0,$
      ewoii_50:    -1.0,$
      ewoiiihb_50: -1.0,$
      ewniiha_50:  -1.0,$

      mstar_avg:   -1.0,$
      age_avg:     -1.0,$
      sfrage_avg:  -1.0,$
      tau_avg:     -1.0,$
      Zmetal_avg:  -1.0,$
      AV_avg:      -1.0,$
      mu_avg:      -1.0,$
      oiiihb_avg:  -1.0,$
      sfr_avg:     -1.0,$ ; instantaneous
      sfr100_avg:  -1.0,$ ; 100 Myr
      b100_avg:    -1.0,$
      b1000_avg:   -1.0,$
      ewoii_avg:   -1.0,$
      ewoiiihb_avg:-1.0,$
      ewniiha_avg: -1.0,$

      mstar_err:    -1.0,$
      age_err:      -1.0,$
      sfrage_err:   -1.0,$
      tau_err:      -1.0,$
      Zmetal_err:   -1.0,$
      AV_err:       -1.0,$
      mu_err:       -1.0,$
      oiiihb_err:   -1.0,$
      sfr_err:      -1.0,$
      sfr100_err:   -1.0,$
      b100_err:     -1.0,$
      b1000_err:    -1.0,$
      ewoii_err:    -1.0,$
      ewoiiihb_err: -1.0,$
      ewniiha_err:  -1.0}

; add the photometric redshift info    
    if keyword_set(photoz) then begin
       pofz = {$
         photoz:      fltarr(5),$ ; three best redshifts
         photoz_err:  fltarr(5),$ ; three best redshifts
         photoz_chi2: fltarr(5),$ ; corresponding chi2 values
         z:              -1.0,$ ; maximum likelihood redshift
         z_err:          -1.0}
       isedfit = struct_addtags(struct_addtags(temporary(isedfit1),$
         struct_addtags(best,qmed)),pofz)
    endif else begin
       isedfit = struct_addtags(temporary(isedfit1),$
         struct_addtags(best,qmed))
    endelse
    isedfit = replicate(temporary(isedfit),ngal)

    if n_elements(ra) ne 0L and n_elements(dec) ne 0L then begin
       isedfit.ra = ra
       isedfit.dec = dec
    endif

; initialize the posterior distribution structure
    isedfit_post = {$
      draws:         lonarr(ndraw)-1,$
      chi2:          fltarr(ndraw)-1,$
      totalmass:     fltarr(ndraw)-1,$
      totalmass_err: fltarr(ndraw)-1}
    if keyword_set(photoz) then isedfit_post = $
      struct_addtags(isedfit_post, {pofz:  fltarr(params.nzz)})
    isedfit_post = replicate(temporary(isedfit_post),ngal)
    
return, isedfit
end

function zmin_parabola, xx, yy
; fit a parabolla; support routine for ISEDFIT_FIND_ZMIN
    ss = sort(xx)
    xx = xx[ss]
    yy = yy[ss]

    a = xx[0]
    b = xx[1]
    c = xx[2]

    fa = yy[0]
    fb = yy[1]
    fc = yy[2]

    xmin = b - 0.5 * ( (b-a)^2*(fb-fc)-(b-c)^2*(fb-fa) ) / $
      ( (b-a)*(fb-fc) - (b-c) *(fb-fa))
return, xmin
end

function isedfit_find_zmin, z, ochi2, nmin=nmin
; find the first N minima in a plot of redshift vs chi2; based on code
; written by Richard Cool for the PRIMUS survey

    if n_elements(nmin) eq 0 then nmin = 5
    
    minmask = intarr(n_elements(ochi2))
    npix = n_elements(z)
    dchi2 = deriv(findgen(npix),ochi2)
    
    delvarx, minval
    delvarx, minima
    
    while total(minmask) lt npix and n_elements(minval) lt nmin do begin
       chi2 = ochi2 + 1D64*minmask
       junk = min(chi2,jmin)
       lgrow = 0
       stoplgrow = 0
       while stoplgrow eq 0 do begin
          lgrow++
          if jmin-lgrow lt 0 then stoplgrow =1 else begin
             if dchi2[jmin-lgrow] ge 0 then stoplgrow = 1
          endelse
       endwhile
       minmask[jmin-(lgrow-1):jmin] = 1
       
       rgrow = 0
       stoprgrow = 0
       while stoprgrow eq 0 do begin
          rgrow++
          if jmin+rgrow gt npix-1 then stoprgrow = 1 else begin
             if dchi2[jmin+rgrow] lt 0 then stoprgrow = 1
          endelse
       endwhile
       minmask[jmin:jmin+(rgrow-1)] = 1
       
       range = [jmin-1,jmin+1]
       if range[0] lt 0 then range = [0,2]
       if range[1] gt npix-1 then range = [npix-3,npix-1]
       
       range1 = [jmin-2,jmin+2]
       if range1[0] lt 0 then range1 = [0,4]
       if range1[1] gt npix-1 then range1 = npix-[5,1]
       zmin = zmin_parabola(z[range[0]:range[1]], $
         ochi2[range[0]:range[1]])
       mval = interpol(ochi2[range[0]:range[1]], $
         z[range[0]:range[1]],zmin)
       err = interpol(abs(z[range1[0]:range1[1]]-zmin), $
         ochi2[range1[0]:range1[1]],mval+1.0)
       
       if jmin lt 2 or jmin gt n_elements(z)-2 then begin
          zmin = z[jmin]
          err = -1
          mval = ochi2[jmin]
       endif
       
       if n_elements(minima) eq 0 then begin
          minima = zmin
          minval = mval
          minerr = err
          rminima = z[jmin]
          rminval = ochi2[jmin]
       endif else begin
          minima = [minima,zmin]
          minval = [minval,mval]
          minerr = [minerr,err]
          rminima = [rminima,z[jmin]]
          rminval = [rminval,ochi2[jmin]]
       endelse
    endwhile
        
    if n_elements(minval) gt 0 then begin
       ss = sort(minval)
       minval = minval[ss]
       minima = minima[ss]
       minerr = minerr[ss]
       
       output = {minval:0.0, $
         minima : 0.0, $
         minerr : 0.0, $
         edges : 0}
       output = replicate(output, n_elements(ss))
       output.minval = minval
       output.minima = minima
       output.minerr = minerr

; temporary hack; drop redshift minima that are on or near the edge
; (the right way to do this is with an edge bit flag)
       keep = where(output.minerr gt 0.0)
       output = output[keep]
       return, output 
    endif else return, -1
end

pro isedfit, isedfit_paramfile, maggies, ivarmaggies, z, params=params, $
  thissfhgrid=thissfhgrid, isedfit_dir=isedfit_dir, outprefix=outprefix, $
  index=index, ra=ra, dec=dec, isedfit_results=isedfit_results, $
  isedfit_post=isedfit_post, isedfit_outfile=isedfit_outfile, allages=allages, $
  maxold=maxold, silent=silent, nowrite=nowrite, photoz=photoz, clobber=clobber

    if n_elements(isedfit_paramfile) eq 0 and n_elements(params) eq 0 then begin
       doc_library, 'isedfit'
       return
    endif

; read the parameter file; parse to get the relevant path and
; filenames
    if (n_elements(params) eq 0) then params = $
      read_isedfit_paramfile(isedfit_paramfile,thissfhgrid=thissfhgrid)
    if (n_elements(isedfit_dir) eq 0) then isedfit_dir = get_pwd()

; error checking on the input photometry    
    ndim = size(maggies,/n_dim)
    dims = size(maggies,/dim)
    if (ndim eq 1) then ngal = 1 else ngal = dims[1]  ; number of galaxies
    nfilt = dims[0] ; number of filters

    nmaggies = n_elements(maggies)
    nivarmaggies = n_elements(ivarmaggies)

; allow for photoz mode    
    if keyword_set(photoz) then begin
       nz = ngal
       z = replicate(0.0,nz)
    endif else nz = n_elements(z)
    
    if (nmaggies eq 0L) or (nivarmaggies eq 0L) or $
      (nz eq 0L) then begin
       doc_library, 'isedfit'
       return
    endif

    ndim = size(maggies,/n_dimension)
    if (ndim ne 2L) then begin ; reform into a 2D array
       maggies = reform(maggies,n_elements(maggies),1)
       ivarmaggies = reform(ivarmaggies,n_elements(maggies),1)
    endif

    if (n_elements(maggies) ne n_elements(ivarmaggies)) then begin
       splog, 'Dimensions of MAGGIES and IVARMAGGIES do not match'
       return
    endif
    if (nz ne ngal) then begin
       splog, 'Dimensions of MAGGIES and Z do not match'
       return
    endif
    if (total(finite(maggies) eq 0B) ne 0.0) or $
      (total(finite(ivarmaggies) eq 0B) ne 0.0) then begin
       splog, 'MAGGIES and IVARMAGGIES cannot have infinite values!'
       return
    endif
    if keyword_set(photoz) eq 0 then begin
       if (total(z le 0.0) ne 0.0) then begin
          splog, 'Z should all be positive'
          return
       endif
    endif

; check for RA,DEC    
    if n_elements(ra) ne 0L then begin
       if size(ra,/type) ne 5 then splog, 'Warning: RA should be type double!'
       if n_elements(ra) ne ngal then begin
          splog, 'Dimensions of RA must match the input number of galaxies!'
          return
       endif
    endif else ra = dblarr(ngal)-1D

    if n_elements(dec) ne 0L then begin
       if size(dec,/type) ne 5 then splog, 'Warning: DEC should be type double!'
       if n_elements(dec) ne ngal then begin
          splog, 'Dimensions of DEC must match the input number of galaxies!'
          return
       endif
    endif else dec = dblarr(ngal)-1D
    
; treat each SFHgrid separately
    ngrid = n_elements(params)
    if ngrid gt 1 then begin
       for ii = 0, ngrid-1 do begin
          isedfit, isedfit_paramfile1, maggies, ivarmaggies, z, params=params[ii], $
            isedfit_dir=isedfit_dir, outprefix=outprefix, index=index, $
            ra=ra, dec=dec, isedfit_results=isedfit_results, isedfit_post=isedfit_post, $
            isedfit_outfile=isedfit_outfile, allages=allages, maxold=maxold, $
            silent=silent, nowrite=nowrite, clobber=clobber, photoz=photoz
       endfor 
       return
    endif

    fp = isedfit_filepaths(params,isedfit_dir=isedfit_dir,outprefix=outprefix)

    isedfit_outfile = fp.isedfit_dir+fp.isedfit_outfile
    if file_test(isedfit_outfile+'.gz') and (keyword_set(clobber) eq 0) and $
      (keyword_set(nowrite) eq 0) then begin
       splog, 'Output file '+isedfit_outfile+' exists; use /CLOBBER'
       return
    endif

; fit the requested subset of objects and return
    if (n_elements(index) ne 0L) then begin
       isedfit, isedfit_paramfile1, maggies[*,index], ivarmaggies[*,index], $
         z[index], params=params, isedfit_dir=isedfit_dir, outprefix=outprefix, $
         ra=ra[index], dec=dec[index], isedfit_results=isedfit_results1, $
         isedfit_post=isedfit_post1, isedfit_outfile=isedfit_outfile, allages=allages, $
         maxold=maxold, silent=silent, /nowrite, clobber=clobber, photoz=photoz
       isedfit_results = init_isedfit(ngal,nfilt,params=params,$
         ra=ra,dec=dec,isedfit_post=isedfit_post,photoz=photoz)
       isedfit_results[index] = isedfit_results1
       isedfit_post[index] = isedfit_post1
       if (keyword_set(nowrite) eq 0) then begin
          im_mwrfits, isedfit_results, isedfit_outfile, /clobber, silent=silent
          im_mwrfits, isedfit_post, fp.isedfit_dir+fp.post_outfile, /clobber, silent=silent
       endif
       return
    endif

    if (keyword_set(silent) eq 0) then begin
       splog, 'SPSMODELS='+strtrim(params.spsmodels,2)+', '+$
         'REDCURVE='+strtrim(params.redcurve,2)+', IMF='+$
         strtrim(params.imf,2)+', '+'SFHGRID='+$
         string(params.sfhgrid,format='(I2.2)')
    endif

; filters and redshift grid
    filterlist = strtrim(params.filterlist,2)
    nfilt = n_elements(filterlist)

    redshift = params.redshift
    nredshift = params.nzz
    if keyword_set(photoz) eq 0 then begin
       if (min(z)-min(redshift) lt -1E-3) or $
         (max(z)-max(redshift) gt 1E-3) then begin
          splog, 'Need to rebuild model grids using a wider redshift grid!'
          return
       endif
    endif

;   redshift1 = reform(rebin(reform(redshift,1,nredshift),$
;     params.nmodel,nredshift),params.nmodel*nredshift)
;   redshift1 = rebin(reform(redshift,1,nredshift),params.nmodel,nredshift)

; if REDSHIFT is not monotonic then FINDEX(), below, can't be
; used to interpolate the model grids properly; this only really
; matters if USE_REDSHIFT is passed
    if params.nzz gt 1 then begin
       if monotonic(redshift) eq 0 then begin
          splog, 'REDSHIFT should be a monotonically increasing or decreasing array!'
          return
       endif
    endif

; initialize the output structure(s)
    isedfit_results = init_isedfit(ngal,nfilt,params=params,$
      ra=ra,dec=dec,isedfit_post=isedfit_post,photoz=photoz)
    isedfit_results.isedfit_id = lindgen(ngal)
    isedfit_results.maggies = maggies
    isedfit_results.ivarmaggies = ivarmaggies
    if keyword_set(photoz) eq 0 then isedfit_results.z = z

; read all the model photometry and modelgrid parameters into memory;
; even for a significant number of models (~100,000) this should be
; manageable on a modern machine
    ngalchunk = ceil(ngal/float(params.galchunksize))    
    nchunk = params.nmodelchunk

    splog, 'Reading all the model parameters into memory.'
    for ichunk = 0, nchunk-1 do begin
       print, format='("ISEDFIT: Chunk ",I0,"/",I0, A10,$)', ichunk+1, nchunk, string(13b)
       chunkfile = fp.models_chunkfiles[ichunk]
;      if (keyword_set(silent) eq 0) then splog, 'Reading '+chunkfile
       modelgrid1 = gz_mrdfits(chunkfile,1,/silent)
       modelmaggies1 = modelgrid1.modelmaggies
       modelgrid1 = struct_trimtags(temporary(modelgrid1),except='modelmaggies')
       if ichunk eq 0 then begin
          modelmaggies = temporary(modelmaggies1)
          modelgrid = temporary(modelgrid1)
       endif else begin
          modelmaggies = [[[temporary(modelmaggies)]],[[temporary(modelmaggies1)]]]
          modelgrid = [temporary(modelgrid),temporary(modelgrid1)]
       endelse
    endfor

; in photoz mode we need to loop on each galaxy because the
; bookkeeping is a little more manageable; in the standard mode we
; loop on each galaxy chunk
    t1 = systime(1)
    mem1 = memory(/current)
    if keyword_set(photoz) then begin
       zindx = lindgen(nredshift)
       maxage = lf_z2t(redshift,omega0=params.omega0,$ ; [Gyr]
         omegal0=params.omegal)/params.h100

; we need an NMODELxNREDSHIFT copy of all the model parameters       
       modelgrid = reform(cmreplicate(temporary(modelgrid),nredshift),$
         params.nmodel*nredshift)

       for igal = 0L, ngal-1 do begin
          if ((igal mod 5) eq 0) then print, format='("ISEDFIT: '+$
            'Galaxy ",I0,"/",I0, A10,$)', igal, ngal, string(13b)
          maggies1 = rebin(maggies[*,igal],nfilt,nredshift) ; copy the photometry 
          ivarmaggies1 = rebin(ivarmaggies[*,igal],nfilt,nredshift)

; compute chi2 over all models and redshifts
          galaxygrid = isedfit_chi2(maggies1,ivarmaggies1,modelmaggies,$
            modelage=modelgrid[0:params.nmodel-1].age,maxage=maxage,$
            zindx=zindx,nzz=nredshift,nminphot=params.nminphot,/silent)

; marginalize over all the models and then compute the first three
; redshift minima; the factor in the numerator is to account for
; models/templates that violated the maxage prior
          pofz = total(exp(-0.5*(galaxygrid.chi2-min(galaxygrid.chi2)))/rebin(reform($
            total(galaxygrid.chi2 lt 1E6,1),1,nredshift),params.nmodel,nredshift),1)
          if total(pofz,/double) eq 0.0 then message, 'Deal with zero total probability.'
          pofz = pofz/total(pofz,/double)
          isedfit_post[igal].pofz = pofz

          zmin = isedfit_find_zmin(redshift,-2.0*alog(pofz>1D-20),nmin=nmin)
          isedfit_results[igal].photoz = zmin.minima
          isedfit_results[igal].photoz_err = zmin.minerr
          isedfit_results[igal].photoz_chi2 = zmin.minval
          isedfit_results[igal].z = isedfit_results[igal].photoz[0]
          isedfit_results[igal].z_err = isedfit_results[igal].photoz_err[0]

;         djs_plot, redshift, pofz, xsty=3, ysty=3, psym=-8
;         for bb = 0, n_elements(zmin)-1 do djs_oplot, zmin[bb].minima*[1,1], !y.crange, color='red'
;         struct_print, zmin
;         cc = get_kbrd(1)
          
; get the posterior probabilities on all the parameters
          galaxygrid = reform(temporary(galaxygrid),params.nmodel*nredshift)

          temp_isedfit_post = isedfit_post[igal]
          isedfit_results[igal] = isedfit_posterior(isedfit_results[igal],$
            modelgrid=modelgrid,galaxygrid=galaxygrid,params=params,$
            isedfit_post=temp_isedfit_post,/silent)
          isedfit_post[igal] = temporary(temp_isedfit_post) ; pass-by-value

; now recompute and then minimize chi2 to get the model parameters at
; the best redshift 
          galaxygrid = isedfit_chi2(maggies[*,igal],ivarmaggies[*,igal],$
            modelmaggies,modelage=modelgrid[0:params.nmodel-1].age,$
            maxage=lf_z2t(isedfit_results[igal].z,omega0=params.omega0,$ ; [Gyr]
            omegal0=params.omegal)/params.h100,nminphot=params.nminphot,$
            nzz=nredshift,zindx=findex(redshift,isedfit_results[igal].z),$
            /silent)
          isedfit_results[igal] = isedfit_posterior(isedfit_results[igal],$
            modelgrid=modelgrid[0:params.nmodel-1],galaxygrid=galaxygrid,$
            params=params,/silent,/bestfitonly)
       endfor                   ; close galaxy loop
    endif else begin
; if NZZ=1 we don't use ZINDX for the interpolation
       if params.nzz gt 1 then zindx = findex(redshift,z)
; do not allow the galaxy to be older than the age of the universe at
; redshift Z
       maxage = lf_z2t(z,omega0=params.omega0,$ ; [Gyr]
         omegal0=params.omegal)/params.h100
       for gchunk = 0L, ngalchunk-1 do begin
          g1 = gchunk*params.galchunksize
          g2 = ((gchunk*params.galchunksize+params.galchunksize)<ngal)-1
          gnthese = g2-g1+1
          gthese = lindgen(gnthese)+g1
; loop on each "chunk" of output from ISEDFIT_MODELS
          t0 = systime(1)
          mem0 = memory(/current)
; compute chi2
          galaxygrid = isedfit_chi2(maggies[*,gthese],ivarmaggies[*,gthese],$
            modelmaggies,modelage=modelgrid.age,maxage=maxage[gthese],$
            zindx=zindx[gthese],gchunk=gchunk,ngalchunk=ngalchunk,ichunk=1,$
            nchunk=1,nminphot=params.nminphot,nzz=params.nzz,allages=allages,$
            silent=silent,maxold=maxold)
; minimize chi2
;         if (keyword_set(silent) eq 0) then splog, 'Building the posterior distributions...'
          temp_isedfit_post = isedfit_post[gthese]
          isedfit_results[gthese] = isedfit_posterior(isedfit_results[gthese],$
            modelgrid=modelgrid,galaxygrid=galaxygrid,params=params,$
            isedfit_post=temp_isedfit_post)
          isedfit_post[gthese] = temporary(temp_isedfit_post) ; pass-by-value
          if keyword_set(silent) eq 0 and gchunk eq 0 then begin
             dt0 = systime(1)-t0
             if dt0 lt 60.0 then begin
                tfactor = 1.0
                suff = 'sec'
             endif else begin
                tfactor = 60.0
                suff = 'min'
             endelse
             splog, 'First GalaxyChunk = '+string(dt0/tfactor,format='(G0)')+$
               ' '+suff+', '+strtrim(string((memory(/high)-mem0)/$
               1.07374D9,format='(F12.3)'),2)+' GB'
          endif 
       endfor                   ; close GalaxyChunk
       if keyword_set(silent) eq 0 then begin
          dt1 = systime(1)-t1
          if dt0 lt 60.0 then begin
             tfactor = 1.0
             suff = 'sec'
          endif else begin
             tfactor = 60.0
             suff = 'min'
          endelse
          splog, 'All GalaxyChunks = '+string(dt1/tfactor,format='(G0)')+$
            ' '+suff+', '+strtrim(string((memory(/high)-mem1)/$
            1.07374D9,format='(F12.3)'),2)+' GB'
       endif
    endelse

;; some legacy code
;; loop on each "chunk" of output from ISEDFIT_MODELS
;          t0 = systime(1)
;          mem0 = memory(/current)
;          for ichunk = 0, nchunk-1 do begin
;             print, format='("ISEDFIT: Chunk ",I0,"/",I0, A10,$)', ichunk+1, nchunk, string(13b)
;             chunkfile = fp.models_chunkfiles[ichunk]
;;            if (keyword_set(silent) eq 0) then splog, 'Reading '+chunkfile
;             modelchunk = gz_mrdfits(chunkfile,1,/silent)
;; compute chi2
;             galaxychunk = isedfit_chi2(maggies[*,gthese],ivarmaggies[*,gthese],$
;               modelchunk,maxage,zindx,gchunk=gchunk,ngalchunk=ngalchunk,ichunk=ichunk,$
;               nchunk=nchunk,nminphot=params.nminphot,nzz=params.nzz,allages=allages,$
;               silent=silent,maxold=maxold)
;             if (ichunk eq 0) then begin
;                galaxygrid = temporary(galaxychunk)
;                modelgrid = temporary(modelchunk)
;             endif else begin
;                galaxygrid = [temporary(galaxygrid),temporary(galaxychunk)]
;                modelgrid = [temporary(modelgrid),temporary(modelchunk)]
;             endelse
;          endfor                ; close ModelChunk
    
; write out the final structure and the full posterior distributions
    if keyword_set(nowrite) eq 0 then begin
       im_mwrfits, isedfit_results, isedfit_outfile, silent=silent, /clobber
       im_mwrfits, isedfit_post, fp.isedfit_dir+fp.post_outfile, silent=silent, /clobber
    endif

return
end
