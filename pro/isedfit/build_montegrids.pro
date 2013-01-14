;+
; NAME:
;   BUILD_MONTEGRIDS
;
; PURPOSE:
;   Build the Monte Carlo based star formation history grids for use
;   with iSEDfit.  
;
; INPUTS: 
;   sfhgrid_paramfile - parameter file describing the star formation
;     history (SFH) priors
;
; OPTIONAL INPUTS: 
;   supergrid_paramfile - file describing the supergrid parameters 
;     (parameters can be overridden using SFHGRID, SYNTHMODELS, IMF,
;     and REDCURVE optional inputs) 
;
;   thissupergrid - if SUPERGRID_PARAMFILE contains multiple
;     supergrids then build this SUPERGRID (may be a vector)
;
;   sfhgrid - SFH grid number to build
;   synthmodels - population synthesis models to use (see the
;     corresponding BUILD_*_SSP routine for details)
;       bc03 - 
;       bc03_lowres - low-resolution BC03 models
;       basti - solar-scaled models
;       basti_ae - alpha-enhanced models
;       pegase - 
;       maraston05 - 
;
;   imf - initial mass function; the available IMFs depend on which
;     SYNTHMODELS are adopted, but as of 2011 Aug 24 they are:
;       chab = Chabrier (2003): bc03, bc03_lowres, fsps
;       salp = Salpeter (1955): bc03, bc03_lowres, fsps, maraston05,
;         pegase 
;       kroupa01 = Kroupa (2001): basti, basti_ae, fsps, maraston05,
;         pegase 
;
;   redcurve - reddening curve; current options are: 
;    -1 = none
;     0 = Calzetti 2000 
;     1 = Charlot & Fall 2000
;     2 = O'Donnell 1994 (i.e., standard Milky Way)
;     3 = SMC
; 
;   isedfit_dir - base pathname for iSEDfit files; the Monte Carlo
;     grids themselves are written to ISEDFIT_DIR/MONTEGRIDS; note
;     that a directory will be created if none exists
;   montegrids_dir - override ISEDFIT_DIR+'/'+MONTEGRIDS
;
; KEYWORD PARAMETERS: 
;   clobber - delete old files from a previous call to this routine,
;     including all the Monte Carlo distributions of parameter values
;   debug - make some debugging plots and wait for a keystroke 
;
; OUTPUTS: 
;   The grids get written out in a data model that ISEDFIT
;   understands, and which should be transparent to the user. 
;
; COMMENTS:
;   Using SUPERGRID_PARAMFILE is the recommended way to call this
;   routine. 
;
; TODO:
;   Build diagnostic plots showing the range of parameters spanned by
;   the models.
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2009 Feb 20, NYU - written
;   jm09may22nyu - generalized how the reddening curves are treated 
;   jm10jul11ucsd - include SFR measures on various timescales and the
;     birthrate parameter at each age
;   jm10nov12ucsd - streamlined and parameter file adopted
;   jm11aug24ucsd - documentation updated
;   jm13jan13siena - semi-major changes in preparation for public
;     release 
;
; Copyright (C) 2009-2011, 2013, John Moustakas
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

function init_montegrid, nmodel, nage, imf=imf, nmaxburst=nmaxburst
; initialize the structure containing the parameters of the Monte
; Carlo grid

    burstarray1 = -1.0
    if (nmaxburst gt 1) then burstarray1 = fltarr(nmaxburst)-1.0
    
    montegrid = {$
      chunkindx:                 -1,$
      modelindx:                 -1,$
      delayed:                    0,$ ; delayed tau model?
      bursttype:                  0,$ ; burst type: 0=step function (default); 1=gaussian; 2=step function with exponential wings
      tau:                     -1.0,$
      Z:                       -1.0,$
      av:                       0.0,$ ; always initialize with zero to accommodate the dust-free models!
      mu:                       1.0,$ ; always default to 1.0!
      nburst:                     0,$
      trunctau:                -1.0,$ ; burst truncation time scale
      minage:                  -1.0,$
      maxage:                  -1.0,$
      mintburst:               -1.0,$
      maxtburst:               -1.0,$
      tburst:           burstarray1,$
      dtburst:          burstarray1,$ ; 
      fburst:           burstarray1,$ ; burst mass fraction
;     aburst:           burstarray1,$ ; burst amplitude
      age:             fltarr(nage)}
    montegrid = replicate(montegrid,nmodel)

return, montegrid
end

function get_bracket_files, info, sspinfo=sspinfo
; get the two SSPs that bracket the desired metallicity, to allow for
; interpolation
    ninfo = n_elements(info)
    bracket_files = strarr(2,ninfo)

    for ii = 0L, ninfo-1 do begin
       indx = findex(sspinfo.Z,info[ii].Z)
       below = fix(indx)
       above = ceil(indx)
       bracket_files[0,ii] = sspinfo.sspfile[below]
       bracket_files[1,ii] = sspinfo.sspfile[above]
    endfor
    
return, bracket_files
end

function read_and_interpolate, files, info=info, $
  ssppath=ssppath, fits_grid=fits_grid
; read and interpolate the base model file onto the desired tau grid

    files = strtrim(files,2)
    junk1 = file_search(ssppath+files[0],count=c1)
    junk2 = file_search(ssppath+files[1],count=c2)
    if (c1 eq 0) or (c2 eq 0) then message, 'Problem finding files'
    
    nobj = n_elements(info)
    if (strtrim(files[0],2) eq strtrim(files[1],2)) then begin
       fits = mrdfits(ssppath+files[0],1,/silent)
       fits = replicate(temporary(fits),nobj) ; identical for all objects
    endif else begin
       fits_grid = [mrdfits(ssppath+files[0],1,/silent),$
         mrdfits(ssppath+files[1],1,/silent)]
       fits = replicate(fits_grid[0],nobj)
       indx = findex(fits_grid.Z,info.Z)
       fits.flux = interpolate(fits_grid.flux,indx)
       fits.mstar = interpolate(fits_grid.mstar,indx)
       fits.Z = interpolate(fits_grid.Z,indx)
    endelse       

return, fits
end

pro build_grid, montegrid, chunkinfo, ssppath=ssppath, $
  sspinfo=sspinfo, redcurve=redcurve, params=params, debug=debug
; this is the main driver routine for building the grid

; reddening curves
    if (n_elements(redcurve) gt 0) then begin
       case redcurve of
         -1: 
          0: calzetti = 1
          1: charlot = 1
          2: odonnell = 1
          3: smc = 1
          else: message, 'Unrecognized reddening curve'
       endcase
    endif    

; build the fiducial output age grid
    nmodel = chunkinfo.nmodel
    nchunk = chunkinfo.nchunk
    chunksize = chunkinfo.chunksize

    for ichunk = 0, nchunk-1 do begin
       i1 = ichunk*chunksize
       i2 = ((ichunk*chunksize+chunksize)<nmodel)-1L
       nthese = i2-i1+1L & these = lindgen(nthese)+i1
       outinfo = montegrid[these]
       outinfo.chunkindx = ichunk

; get the two SSPs that bracket the desired metallicity, to allow for
; interpolation; do this here so that we're not reading in the same
; set of models multiple times
       bracket_sspfile = get_bracket_files(outinfo,sspinfo=sspinfo)
       uindx = uniq(bracket_sspfile[0,*],sort(bracket_sspfile[0,*]))
       ubracket_sspfile = bracket_sspfile[*,uindx]
       nssp = n_elements(uindx)

       if (n_elements(tbc) eq 0) then tbc = 0.01D9 ; dust dispersal timescale [10 Myr]

       for issp = 0, nssp-1 do begin ; loop on SSP metallicity
          indx1 = where(ubracket_sspfile[0,issp] eq bracket_sspfile[0,*],nindx1)
          outinfo[indx1].modelindx = indx1
          sspfits = read_and_interpolate(ubracket_sspfile[*,issp],$
            info=outinfo[indx1],ssppath=ssppath)
          npix = n_elements(sspfits[0].wave)

; initialize the output SED data structure
          if (issp eq 0) then outinfo = struct_addtags(temporary(outinfo),$
            replicate({mstar: fltarr(params.nage), $
            wave: float(sspfits[0].wave), $
            flux: fltarr(npix,params.nage)},nthese))

          delvarx, rv
          klam = k_lambda(outinfo[indx1[0]].wave,r_v=rv,charlot=charlot,$
            calzetti=calzetti,odonnell=odonnell,smc=smc,/silent)
          klam = rebin(reform(klam,npix,1),npix,n_elements(sspfits[0].age)) ; [npix,nage]
;         if keyword_set(charlot) then nsmth = total((sspfits[0].age gt 0.9*tbc) and (sspfits[0].age lt 1.1*tbc))

; convolve each model with the specified SFH
          for jj = 0L, nindx1-1 do begin
             print, format='("Chunk=",I4.4,"/",I4.4,", SSP=",I4.4,"/",I4.4,", '+$
               'Model=",I4.4,"/",I4.4,"   ",A5,$)', ichunk+1, nchunk, issp+1, $
               nssp, jj+1, nindx1, string(13b)
; attenuate
             alam = klam*(outinfo[indx1[jj]].av/rv)
             if keyword_set(charlot) then begin
                old = where(sspfits[jj].age gt tbc,nold)
                if (nold ne 0) then alam[*,old] = outinfo[indx1[jj]].mu*alam[*,old]
             endif 
             sspfits[jj].flux = sspfits[jj].flux*10.0^(-0.4*alam)

; do the convolution; tau=0 and no bursts is a special case
             outage = outinfo[indx1[jj]].age
             if (outinfo[indx1[jj]].tau eq 0.0) and (outinfo[indx1[jj]].nburst eq 0) then begin ; special case
                ageindx = findex(sspfits[jj].age,outage*1D9)
                outflux = interpolate(sspfits[jj].flux,lindgen(npix),ageindx,/grid)
                outmstar = interpolate(sspfits[jj].mstar,ageindx)
                outmgal = outmstar*0+1
                outsfr = isedfit_reconstruct_sfh(outinfo[indx1[jj]],$
                  outage=im_double(outage),debug=debug)
             endif else begin
                outflux = isedfit_convolve_sfh(sspfits[jj],info=outinfo[indx1[jj]],$
                  time=im_double(outage),mstar=sspfits[jj].mstar,cspmstar=outmstar,$
                  sfh=outsfr,nsamp=1.0,debug=debug);,/ylog,xr=[4,7],yrange=[1D-14,1D-6])
;               test = isedfit_reconstruct_sfh(outinfo[indx1[jj]],outage=outage,/debug,xr=[3.5,8])
             endelse

             inf = where(finite(outflux) eq 0)
             if (inf[0] ne -1) then message, 'Bad bad bad'

             outinfo[indx1[jj]].age = outage
             outinfo[indx1[jj]].mstar = outmstar
             outinfo[indx1[jj]].flux = outflux
          endfor 
       endfor    ; close SSP loop 
       im_mwrfits, outinfo, chunkinfo.chunkfiles[ichunk], /clobber
    endfor 
    
return
end    

pro build_montegrids, sfhgrid_paramfile, supergrid_paramfile=supergrid_paramfile, $
  thissupergrid=thissupergrid, sfhgrid=sfhgrid, synthmodels=synthmodels, imf=imf, $
  redcurve=redcurve, isedfit_dir=isedfit_dir, montegrids_dir=montegrids_dir, $
  clobber=clobber, debug=debug

; read the SFHGRID parameter file
    if n_elements(sfhgrid_paramfile) eq 0 then begin
       doc_library, 'build_montegrids'
       return
    endif

; read the SUPERGRID parameter file, if given    
    if n_elements(supergrid_paramfile) ne 0 then begin
       super = read_supergrid_paramfile(supergrid_paramfile,supergrid=thissupergrid)
       if n_elements(sfhgrid) eq 0 then sfhgrid = super.sfhgrid
       if n_elements(synthmodels) eq 0 then synthmodels = super.synthmodels
       if n_elements(imf) eq 0 then imf = super.imf
       if n_elements(redcurve) eq 0 then redcurve = super.redcurve
    endif else begin
       if (n_elements(sfhgrid) eq 0) or (n_elements(synthmodels) eq 0) or $
         (n_elements(imf) eq 0) or (n_elements(redcurve) eq 0) then begin
          splog, 'You must either provide SUPERGRID_PARAMFILE or *all* of '+$
            'SFHGRID, SYNTHMODELS, IMF, and REDCURVE'
          return
       endif
    endelse
       
; call this routine iteratively
    nsfh = n_elements(sfhgrid)
    if nsfh gt 1 then begin
       if nsfh ne n_elements(synthmodels) or nsfh ne n_elements(imf) or $
         nsfh ne n_elements(redcurve) then message, 'SFHGRID, SYNTHMODELS, '+$
         'IMF, and REDCURVE must have the same number of elements!'
       for ii = 0, n_elements(sfhgrid)-1 do begin
          build_montegrids, sfhgrid_paramfile, sfhgrid=sfhgrid[ii], $
            synthmodels=synthmodels[ii], imf=imf[ii], redcurve=redcurve[ii], $
            isedfit_dir=isedfit_dir, montegrids_dir=montegrids_dir, $
            clobber=clobber, debug=debug
       endfor
       return
    endif

; read the SFHGRID parameter file and get the reddening curve
    params = read_sfhgrid_paramfile(sfhgrid_paramfile,thissfhgrid=sfhgrid)
    redcurvestring = redcurve2string(redcurve)

; read the SSP information structure    
    ssppath = getenv('ISEDFIT_SSP_DIR')
    if (file_test(ssppath,/dir) eq 0) then begin
       splog, 'Verify that ${ISEDFIT_SSP_DIR} environment variable is defined!'
       return
    endif
    ssppath=ssppath+'/'
    sspinfofile = ssppath+'info_'+synthmodels+'_'+imf+'.fits.gz'
    if (file_test(sspinfofile) eq 0) then begin
       splog, 'SSP info file '+sspinfofile+' not found!'
       splog, 'Run the appropriate BUILD_*_SSP code!'
       return
    endif
    sspinfo = mrdfits(sspinfofile,1,/silent)
    
; make directories and delete old files
    if (n_elements(isedfit_dir) eq 0) then isedfit_dir = './'
    if (n_elements(montegrids_dir) eq 0) then montegrids_dir = isedfit_dir+'montegrids/'
;   if (n_elements(montegrids_dir) eq 0) then montegrids_dir = $
;     '${MONTEGRIDS_DIR}/'
    sfhgridstring = 'sfhgrid'+string(sfhgrid,format='(I2.2)')
    sfhgridpath = montegrids_dir+sfhgridstring+$
      '/'+synthmodels+'/'+redcurvestring+'/'

    if (file_test(sfhgridpath,/dir) eq 0) then begin
       splog, 'Making directory '+sfhgridpath
       spawn, 'mkdir -p '+sfhgridpath
    endif

; ---------------------------------------------------------------------------
; first major step: build the Monte Carlo grid, if it doesn't
; already exist 
    montefile = sfhgridpath+imf+'_montegrid.fits'
    if (file_test(montefile+'.gz') eq 0) or keyword_set(clobber) then begin
       cc = 'N'
       if (keyword_set(clobber) eq 0) then begin
          splog, 'Delete all *'+imf+'* files from '+sfhgridpath+' [Y/N]?'
          cc = get_kbrd(1)
       endif
       if keyword_set(clobber) or (strupcase(cc) eq 'Y') then $
         spawn, '/bin/rm -f '+sfhgridpath+'*'+imf+'*.fits.gz', /sh
       
       splog, 'Building SFHGRID='+sfhgridstring+' REDCURVE='+redcurvestring+$
         ' NMODEL='+string(params.nmonte*params.nage,format='(I0)')

; compute the maximum number of bursts, if any
       if (params.pburst le 0.0) then nmaxburst = 0 else $
         nmaxburst = ceil(long(100D*(params.maxtburst-params.mintburst)/params.interval_pburst)/100D)
;      if (params.pburst le 0.0) then nmaxburst = 0 else $
;        nmaxburst = ceil(long(100D*(params.maxage-params.minage)/params.interval_pburst)/100D)

       montegrid = init_montegrid(params.nmonte,params.nage,imf=imf,nmaxburst=nmaxburst)

       montegrid.delayed = params.delayed ; delayed SFH?
       montegrid.minage = params.minage
       montegrid.maxage = params.maxage

; draw uniformly from linear TAU, or 1/TAU?
       tau = randomu(seed,params.nmonte)*(params.tau[1]-params.tau[0])+params.tau[0]
       if params.oneovertau eq 1 and params.delayed eq 1 then $
         message, 'DELAYED and ONEOVERTAU may not work well together.'
       if params.oneovertau eq 1 then tau = 1D/tau
       montegrid.tau = tau

; metallicity; check to make sure that the prior boundaries do not
; exceed the metallicity range available from the chosen SYNTHMODELS
       if (params.Z[0] lt min(sspinfo.Z)) then begin
          splog, 'Adjusting minimum prior metallicity!'
          params.Z[0] = min(sspinfo.Z)
       endif
       if (params.Z[1] gt max(sspinfo.Z)) then begin
          splog, 'Adjusting maximum prior metallicity!'
          params.Z[1] = max(sspinfo.Z)
       endif
       montegrid.Z = randomu(seed,params.nmonte)*(params.Z[1]-params.Z[0])+params.Z[0]

; age; unfortunately I think we have to loop to sort
;      splog, 'TESTING WITH A UNIFORM AGE GRID!'
;      for ii = 0L, params.nmonte-1 do montegrid[ii].age = range(params.minage,params.maxage,params.nage)
       montegrid.age = randomu(seed,params.nage,params.nmonte)*(params.maxage-params.minage)+params.minage
       for ii = 0L, params.nmonte-1 do montegrid[ii].age = montegrid[ii].age[sort(montegrid[ii].age)]

; reddening, if any; Gamma distribution is the default, unless FLATAV==1
       if (params.av[1] gt 0) then begin
          if params.flatav then begin
             montegrid.av = randomu(seed,params.nmonte)*(params.av[1]-params.av[0])+params.av[0] 
          endif else begin
             montegrid.av = params.av[0]*randomu(seed,params.nmonte,gamma=params.av[1])
;            montegrid.av = ((10.0^(randomn(seed,params.nmonte)*params.av[1]+alog10(params.av[0]))>0.0
;            montegrid.av = randomu(seed,params.nmonte,gamma=1.0)*params.av[1]+params.av[0]
;            montegrid.av = 10.0^(randomn(seed,params.nmonte)*params.av[1]+alog10(params.av[0]))
          endelse
       
; "mu" is the Charlot & Fall (2000) factor for evolved stellar
; populations; Gamma distribution is the default, unless FLATMU==1
          if (redcurvestring eq 'charlot') then begin
             if params.flatmu then begin
                montegrid.mu = randomu(seed,params.nmonte)*(params.mu[1]-params.mu[0])+params.mu[0] 
             endif else begin
                montegrid.mu = params.mu[0]*randomu(seed,params.nmonte,gamma=params.mu[1])
;               montegrid.mu = ((10.0^(randomn(seed,params.nmonte)*params.mu[1]+alog10(params.mu[0])))<1.0)>0.0
             endelse
          endif
       endif 

; now assign bursts; note that the bursts can occur outside
; (generally, before) the AGE vector; divide the time vector into
; NMAXBURST intervals of width INTERVAL_PBURST [Gyr]
       ntime = 100
;      ntime = params.nage
       
; type of burst: 0 (step function, default), 1 (gaussian), 2 (step
; function with exponential wings) 
       montegrid.bursttype = params.bursttype
       
       if (nmaxburst gt 0) then begin
          montegrid.mintburst = params.mintburst
          montegrid.maxtburst = params.maxtburst
          
          tburst = dblarr(nmaxburst,params.nmonte)-1.0
          ran = randomu(seed,nmaxburst,params.nmonte)
          for imod = 0L, params.nmonte-1 do begin
             for ib = 0, nmaxburst-1 do begin
                tmin = params.mintburst+ib*params.interval_pburst
                tmax = (params.mintburst+(ib+1)*params.interval_pburst)<params.maxtburst
;               tmin = params.minage+ib*params.interval_pburst
;               tmax = (params.minage+(ib+1)*params.interval_pburst)<params.maxage

                if (tmax le tmin) then message, 'This violates causality!'
                time1 = randomu(seed,ntime)*(tmax-tmin)+tmin
                time1 = time1[sort(time1)]

; compute the cumulative probability, dealing with edge effects correctly
                dtimeleft = (time1-shift(time1,1))/2.0
                dtimeleft[0] = time1[0]-tmin
                dtimeright = (shift(time1,-1)-time1)/2.0
                dtimeright[ntime-1] = tmax-time1[ntime-1]
                prob = params.pburst*total([[dtimeleft],[dtimeright]],2)/params.interval_pburst
                if (prob[0] le 0) then message, 'This should not happen!'

                this = where(total(prob,/cum) gt ran[ib,imod])
                if (this[0] ne -1) then tburst[ib,imod] = time1[this[0]]
;               splog, tmin, tmax, ran[ib,imod], tburst[ib,imod] & if ib eq 0 then print
             endfor
          endfor 
          montegrid.nburst = total(tburst gt -1.0,1) ; total number of bursts
       endif
       
; assign burst strengths and durations       
       hasburst = where(montegrid.nburst gt 0,nhasburst)
       nallburst = long(total(montegrid.nburst))
       if (nhasburst ne 0L) then begin
          if params.flatdtburst then begin
             dtburst = randomu(seed,nallburst)*(params.dtburst[1]-params.dtburst[0])+params.dtburst[0]
          endif else begin
             dtburst = 10.0^(randomu(seed,nallburst)*(alog10(params.dtburst[1])-$
               alog10(params.dtburst[0]))+alog10(params.dtburst[0]))
          endelse

          if params.flatfburst then begin
             fburst = randomu(seed,nallburst)*(params.fburst[1]-params.fburst[0])+params.fburst[0]
          endif else begin
             fburst = 10.0^(randomu(seed,nallburst)*(alog10(params.fburst[1])-$
               alog10(params.fburst[0]))+alog10(params.fburst[0]))
          endelse

          count = 0L
          for ii = 0L, nhasburst-1 do begin ; sort
             if (montegrid[hasburst[ii]].nburst gt 0) then begin
                nb = montegrid[hasburst[ii]].nburst
                good = where(tburst[*,hasburst[ii]] gt -1.0)
                montegrid[hasburst[ii]].tburst = tburst[good,hasburst[ii]]
                montegrid[hasburst[ii]].fburst = fburst[count:count+nb-1]
                montegrid[hasburst[ii]].dtburst = dtburst[count:count+nb-1]
                count = count + nb
             endif
          endfor

; allow a FRACTRUNC fraction of the models with bursts to have
; truncated bursts 
          if (params.fractrunc gt 0D) then begin
             ntrunc = long(params.fractrunc*nhasburst)
             if params.fractrunc eq 1D then trunc = lindgen(ntrunc) else $
               trunc = random_indices(nhasburst,ntrunc)
          endif else ntrunc = 0L
          if (ntrunc gt 0L) then begin
             if (montegrid[0].bursttype ne 1) then message, 'Should use truncated *Gaussian* bursts.'
;            montegrid[hasburst[trunc]].trunctau = params.trunctau ; [Gyr]
             montegrid[hasburst[trunc]].trunctau = randomu(seed,ntrunc)*(params.trunctau[1]-$
               params.trunctau[0])+params.trunctau[0] ; [Gyr]
          endif
       endif  

; write out
       im_mwrfits, montegrid, montefile, /clobber
    endif else begin
; read the Monte Carlo grid
       splog, 'Reading '+montefile+'.gz'
       montegrid = mrdfits(montefile+'.gz',1)
       if (n_elements(montegrid) ne params.nmonte) then begin
          splog, 'Dimensions of MONTEGRID and PARAMS.NMONTE '+$
            'do not match; call with /MAKE_MONTEGRID!'
          return
       endif
    endelse

; ---------------------------------------------------------------------------
; now build the actual composite stellar populations; do it in chunks
; to avoid memory issues and pass along all the parameters
    chunkinfofile = sfhgridpath+imf+'_chunkinfo.fits'
    chunksize = 100 ; number of models per output FITS table

    chunkinfo1 = struct_addtags(params,{imf: imf, $
      synthmodels: synthmodels, sfhgridstring: sfhgridstring, $
      redcurve: redcurvestring, chunksize: chunksize})
    
    nchunk = ceil(params.nmonte/float(chunksize))
    chunkfiles = imf+'_chunk_'+string(lindgen(nchunk)+1,$
      format='(I4.4)')+'.fits'
    chunkinfo = create_struct(chunkinfo1, 'nmodel', params.nmonte, $
      'nchunk', nchunk, 'chunkfiles', sfhgridpath+chunkfiles)

    t0 = systime(1)
    build_grid, montegrid, chunkinfo, ssppath=ssppath+synthmodels+'/', $
      sspinfo=sspinfo, redcurve=redcurve, params=params, debug=debug
    chunkinfo.chunkfiles = file_basename(chunkinfo.chunkfiles)+'.gz'
    splog, 'Total time (min) = ', (systime(1)-t0)/60.0
 
; write out the chunkinfofile
    im_mwrfits, chunkinfo, chunkinfofile, /clobber

return
end

;; QAplot for debugging
;             if keyword_set(debug) then begin
;                im_plotconfig, 6, pos
;                xr = [min(outage)>1,max(outage)]
;                yr = [min(outsfr)>1D-15,max(outsfr)>1D-10]
;
;                im_window, 0, yr=0.8
;                djs_plot, [0], [0], position=pos[*,0], xsty=3, ysty=3, xlog=0, /ylog, $
;                  xrange=xr, yrange=yr, ytitle='SFR (M_{\odot} yr^{-1})', $
;                  xtickname=replicate(' ',20)
;                djs_oplot, outage, outsfr, psym=-6
;                for iburst = 0, outinfo[indx1[jj]].nburst-1 do djs_oplot, $
;                  outinfo[indx1[jj]].tburst[iburst]*[1,1], 10^!y.crange, color='yellow'
;                im_legend, 'N_{burst}='+strtrim(outinfo[indx1[jj]].nburst,2), $
;                  /left, /bottom, box=0, margin=0
;                im_legend, '\tau='+string(outinfo[indx1[jj]].tau,format='(G0.0)')+$
;                  ' Gyr', /right, /top, box=0, margin=0
;
;                djs_plot, outage, outmgal, position=pos[*,1], /noerase, xsty=3, $
;                  xrange=xr, xlog=0, psym=-6, ysty=3, yr=[0,1.05], xtitle='Age (Gyr)'
;                djs_oplot, outage, outmstar, line=5, psym=-6
;                for iburst = 0, outinfo[indx1[jj]].nburst-1 do djs_oplot, $
;                  outinfo[indx1[jj]].tburst[iburst]*[1,1], !y.crange, color='yellow'
;                cc = get_kbrd(1)
;             endif 
