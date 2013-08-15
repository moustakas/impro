;+
; NAME:
;   ISEDFIT_MONTEGRIDS
;
; PURPOSE:
;   Build the Monte Carlo star formation history grids.
;
; INPUTS: 
;   isedfit_paramfile - iSEDfit parameter file
;
; OPTIONAL INPUTS: 
;   params - data structure with the same information contained in
;     ISEDFIT_PARAMFILE (over-rides ISEDFIT_PARAMFILE)
;   thissfhgrid - if ISEDFIT_PARAMFILE contains multiple grids then
;     build this SFHgrid (may be a vector)
;   montegrids_dir - full directory path where the Monte Carlo grids
;     should be written (default 'montegrids' subdirectory of the
;     PWD=present working directory) 
;
;   chunksize, minichunksize - to deal with disk space issues this
;     routine splits the full set of models (NMODEL, as specified in 
;     WRITE_ISEDFIT_PARAMFILE) into CHUNKSIZE (default 5000) sized
;     chunks; furthemore, to circumvent possible memory issues, each
;     chunk is further divided into MINICHUNKSIZE (default 500) sized
;     chunks; all these tricks are transparent to the user and are
;     only documented here for completeness
;
; KEYWORD PARAMETERS: 
;   clobber - delete old files from a previous call to this routine,
;     including all the Monte Carlo distributions of parameter values
;     (the user is prompted to confirm the clobber to avoid
;     inadvertently deleting a large grid!)
;   debug - make some debugging plots and wait for a keystroke
;     (currently not a very useful switch)
;
; OUTPUTS: 
;   The grids (spectra, parameters, etc.) are written out in a data
;   model that iSEDfit understands, and which should be transparent to
;   the user.  However, this routine also generates a QAplot (written
;   to the appropriate subdirectory of MONTEGRIDS_DIR) which should be
;   inspected to ensure that proper parameter priors have been
;   chosen. 
;
; COMMENTS:
;   Some ToDo items:
;     * output the UV slope, beta, and maybe D(4000)
;     * better debugging diagnostic plots 
;
; TODO:
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
;   jm13aug09siena - significant update and rewrite to conform to a
;     new and much simpler data model
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

function init_montegrid, nmodel, nmaxburst=nmaxburst
; internal support routine: initialize the structure containing the
; parameters of the Monte Carlo grid 

    burstarray1 = -1.0
    if (nmaxburst gt 1) then burstarray1 = fltarr(nmaxburst)-1.0
    
    montegrid = {$
      chunkindx:                 -1,$
      modelindx:                 -1,$
      age:                     -1.0,$
      tau:                     -1.0,$
      Zmetal:                  -1.0,$
      AV:                       0.0,$ ; always initialize with zero to accommodate the dust-free models!
      mu:                       1.0,$ ; always default to 1.0!
      oiiihb:                  -1.0,$ ; [OIII] 5007/H-beta
      nlyc:                    -1.0,$ ; number of Lyman-continuum photons
      mstar:                   -1.0,$
      sfr:                     -1.0,$
      sfr100:                  -1.0,$
      sfrage:                  -1.0,$
      nburst:                     0,$
      trunctau:                -1.0,$ ; burst truncation time scale
      tburst:           burstarray1,$
      dtburst:          burstarray1,$ ; 
      fburst:           burstarray1}  ; burst mass fraction
    montegrid = replicate(montegrid,nmodel)

return, montegrid
end

function get_bracket_files, info, sspinfo=sspinfo
; internal support routine: get the two SSPs that bracket the desired
; metallicity, to allow for interpolation
    ninfo = n_elements(info)
    bracket_files = strarr(2,ninfo)

    for ii = 0L, ninfo-1 do begin
       indx = findex(sspinfo.Zmetal,info[ii].Zmetal)
       below = fix(indx)
       above = ceil(indx)
       bracket_files[0,ii] = sspinfo.sspfile[below]
       bracket_files[1,ii] = sspinfo.sspfile[above]
    endfor
    
return, bracket_files
end

function read_and_interpolate, files, info=info, $
  ssppath=ssppath, fits_grid=fits_grid
; internal support routine: read and interpolate the base model file
; onto the desired metallicity grid 

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
       indx = findex(fits_grid.Zmetal,info.Zmetal)
; this interpolate can be very slow because we're building an
; NPIX,NAGE,NOBJ elements array
;      t0 = systime(1)
       fits.flux = interpolate(fits_grid.flux,indx)
;      splog, 'Total time (sec) = ', (systime(1)-t0)
       fits.mstar = interpolate(fits_grid.mstar,indx)
       fits.Zmetal = interpolate(fits_grid.Zmetal,indx)
    endelse       

return, fits
end

function build_modelgrid, montegrid, params=params, debug=debug, $
  sspinfo=sspinfo, ssppath=ssppath, minichunksize=minichunksize, $
  ichunk=ichunk, nchunk=nchunk
; internal support routine: workhorse code which builds the spectra  

    if (n_elements(tbc) eq 0) then tbc = 0.01D9 ; dust dispersal timescale [10 Myr]

; define some reddening curve variables
    case strlowcase(strtrim(params.redcurve,2)) of
       'none': 
       'calzetti': calzetti = 1
       'charlot': charlot = 1
       'odonnell': odonnell = 1
       'smc': smc = 1
       else: message, 'Unrecognized reddening curve'
    endcase

; loop on each mini-chunk of models
    nmodel = n_elements(montegrid)
    nmini = ceil(nmodel/float(minichunksize))
    
    for imini = 0, nmini-1 do begin
       i1 = imini*minichunksize
       i2 = ((imini*minichunksize+minichunksize)<nmodel)-1L
       nthese = i2-i1+1L & these = lindgen(nthese)+i1

       modelgrid1 = montegrid[these]
       
; get the two SSPs that bracket the desired metallicity, to allow for
; interpolation; do this here so that we're not reading in the same
; set of models multiple times; jm13aug07siena: interpolating a large
; number of files is very slow and memory-intensive, so subdivide the
; SSPs into mini-chunks
       bracket_sspfile = get_bracket_files(modelgrid1,sspinfo=sspinfo)
       uindx = uniq(bracket_sspfile[0,*],sort(bracket_sspfile[0,*]))
       ubracket_sspfile = bracket_sspfile[*,uindx]
       nssp = n_elements(uindx)

       for issp = 0, nssp-1 do begin ; loop on SSP metallicity
          sspindx = where(ubracket_sspfile[0,issp] eq bracket_sspfile[0,*],nsspindx)
          modelgrid1[sspindx].modelindx = sspindx+i1

          sspfits = read_and_interpolate(ubracket_sspfile[*,issp],$
            info=modelgrid1[sspindx],ssppath=ssppath)

; build an emission-line friendly wavelength array; this should be
; identical for all SSPs
          if params.nebular then begin
             junk = isedfit_nebular(1D,wave=emwave,vsigma=vsigma,$
               inst_vsigma=sspinfo.inst_vsigma)

             wave = [emwave,sspfits[0].wave]
             wave = float(wave[sort(wave)])
             if nsspindx eq 1L then begin
                sspflux = interpolate(sspfits.flux,findex(sspfits[0].wave,wave),$
                  findgen(n_elements(sspfits[0].age)),/grid)
             endif else begin
                sspflux = interpolate(sspfits.flux,findex(sspfits[0].wave,wave),$
                  findgen(n_elements(sspfits[0].age)),findgen(nsspindx),/grid)
             endelse
          endif else begin
             wave = sspfits[0].wave
             sspflux = sspfits.flux
          endelse
          npix = n_elements(wave)

; initialize the output SED data structure; we waste a little disk
; space by replicating the wavelength array for each model, but
; it's fine, disk is cheap
          if issp eq 0 then begin
             spectags = replicate({ewoii: -1.0, ewoiiihb: -1.0, ewniiha: -1.0, $
               wave: float(wave), flux: fltarr(npix)},nthese)
             modelgrid1 = struct_addtags(temporary(modelgrid1),spectags)
          endif

          delvarx, rv
          klam = k_lambda(wave,r_v=rv,charlot=charlot,calzetti=calzetti,$
            odonnell=odonnell,smc=smc,/silent)
          klam = rebin(reform(klam,npix,1),npix,n_elements(sspfits[0].age)) ; [NPIX,NAGE]

; convolve each model with the specified SFH
          for jj = 0, nsspindx-1 do begin
              print, format='("Chunk=",I4.4,"/",I4.4,", '+$
                'Minichunk=",I4.4,"/",I4.4,", SSP(Z)=",I3.3,"/",I4.4,", '+$
                'Model=",I4.4,"/",I4.4,"   ",A5,$)', ichunk+1, nchunk, $
                imini+1, nmini, issp+1, nssp, jj+1, nsspindx, string(13b)
; attenuation
             alam = klam*(modelgrid1[sspindx[jj]].av/rv)
             if keyword_set(charlot) then begin
                old = where(sspfits[jj].age gt tbc,nold)
                if (nold ne 0) then alam[*,old] = modelgrid1[sspindx[jj]].mu*alam[*,old]
             endif
             sspflux[*,*,jj] = sspflux[*,*,jj]*10.0^(-0.4*alam)

; do the convolution; tau=0 and no bursts is a special case
             outage = modelgrid1[sspindx[jj]].age
             if (modelgrid1[sspindx[jj]].tau eq 0.0) and $
               (modelgrid1[sspindx[jj]].nburst eq 0) then begin ; special case
                ageindx = findex(sspfits[jj].age,outage*1D9)
                outflux = interpolate(sspflux[*,*,jj],lindgen(npix),ageindx,/grid)
                outmstar = interpolate(sspfits[jj].mstar,ageindx)
                outnlyc = interpolate(sspfits[jj].nlyc,ageindx)
                outmgal = outmstar*0+1
             endif else begin
                outflux = isedfit_convolve_sfh(sspflux[*,*,jj],age=sspfits[jj].age,$
                  info=modelgrid1[sspindx[jj]],time=im_double(outage),$
                  mstar=sspfits[jj].mstar,nlyc=sspfits[jj].nlyc,$
                  cspmstar=outmstar,cspnlyc=outnlyc,nsamp=1.0,debug=debug,$
                  delayed=params.delayed,bursttype=params.bursttype)
;               test = isedfit_sfh(modelgrid1[sspindx[jj]],outage=outage,$
;                 delayed=params.delayed,bursttype=params.bursttype,/debug)
             endelse
             inf = where(finite(outflux) eq 0)
             if (inf[0] ne -1) then message, 'Bad bad bad'

; get the instantaneous and 100-Myr averaged SFR
             outsfr = isedfit_sfh(modelgrid1[sspindx[jj]],outage=outage,$
               sfr100=outsfr100,sfrage=outsfrage,b100=b100,mgalaxy=mgal,$
               delayed=params.delayed,bursttype=params.bursttype,$
               debug=debug)

;splog, 'check that the number of nlyc and sfr matches             '
             
; generate the emission-line spectrum, if desired
             if params.nebular then begin
                oiiihb = 10D^modelgrid1[sspindx[jj]].oiiihb
                nebflux = isedfit_nebular(10D^outnlyc,inst_vsigma=sspinfo.inst_vsigma,$
                  vsigma=vsigma,oiiihb=oiiihb,wave=wave,line=line,flam_line=flam_line,$
                  flam_cont=flam_cont)
; attenuate
                alam = klam[*,0]*(modelgrid1[sspindx[jj]].av/rv) ; never decrease by MU, if /CHARLOT
;               djs_plot, wave/1D4, outflux, xr=[0.1,1], xsty=3, ysty=3
;               djs_oplot, wave/1D4, outflux+nebflux, color='green'
;               djs_oplot, wave/1D4, outflux+nebflux*10.0^(-0.4*alam), color='red'
                nebflux = nebflux*10.0^(-0.4*alam)
                flam_line = flam_line*10.0^(-0.4*alam)
                flam_cont = flam_cont*10.0^(-0.4*alam)

; get the EW of some emission lines; the wavelengths correspond to
; [OII], Hb+[OIII], and Ha+[NII]; should probably push this to its own
; function 
                lmin = [3727.42,4861.325,6548.043]
                lmax = [3727.42,5006.842,6730.815]
                lwave = total([[lmin],[lmax]],2)/2.0
                cflux = lwave*0.0
                for cc = 0, n_elements(lwave)-1 do begin
                   factor = sqrt(vsigma^2+sspinfo.inst_vsigma^2)/im_light(/kms)*lwave[cc]
                   cpix = where((wave lt (lmin[cc]-4*factor) and wave gt (lmin[cc]-12*factor)) or $
                     (wave gt (lmax[cc]+4*factor) and wave lt (lmax[cc]+12*factor)))
                   cflux[cc] = djs_median(outflux[cpix]+flam_cont[cpix])
                   if keyword_set(debug) then begin
                      djs_plot, wave/1D4, outflux+nebflux, xsty=3, ysty=3, $
                        xr=[lmin[cc]-30*factor,lmax[cc]+30*factor]/1D4
;                     djs_oplot, wave/1D4, outflux+flam_cont, color='red'
                      djs_oplot, wave[cpix]/1D4, outflux[cpix]+nebflux[cpix], color='green', psym=7
                      djs_oplot, [lwave[cc]]/1D4, [cflux[cc]], psym=8, color='blue', symsize=2
                      djs_oplot, !x.crange, cflux[cc]*[1,1], color='orange', line=5
                      kk = get_kbrd(1)
                   endif
                endfor
                if total(cflux le 0.0) ne 0 then message, 'This should never happen!'

                isoii = where(line.name eq '[OII]_3726' or line.name eq '[OII]_3729')
                isoiiihb = where(line.name eq '[OIII]_4959' or $
                  line.name eq '[OIII]_5007' or line.name eq 'Hbeta')
                isniiha = where(line.name eq '[NII]_6548' or $
                  line.name eq '[NII]_6584' or line.name eq 'Halpha')

                lwave = djs_mean(line[isoii].wave)
                ldust = 10.0^(-0.4*interpol(alam,wave,lwave))
                modelgrid1[sspindx[jj]].ewoii = total(line[isoii].flux*lwave)*ldust/cflux[0]

                lwave = djs_mean(line[isoiiihb].wave)
                ldust = 10.0^(-0.4*interpol(alam,wave,lwave))
                modelgrid1[sspindx[jj]].ewoiiihb = total(line[isoiiihb].flux*lwave)*ldust/cflux[1]

                lwave = djs_mean(line[isniiha].wave)
                ldust = 10.0^(-0.4*interpol(alam,wave,lwave))
                modelgrid1[sspindx[jj]].ewniiha = total(line[isniiha].flux*lwave)*ldust/cflux[2]

;               splog, modelgrid1[sspindx[jj]].age, modelgrid1[sspindx[jj]].ewoii, $
;                 modelgrid1[sspindx[jj]].ewoiiihb, modelgrid1[sspindx[jj]].ewniiha
             endif else nebflux = outflux*0.0

; ToDo: add additional interesting outputs like beta (UV slope)

; pack it in             
             modelgrid1[sspindx[jj]].age = outage
             modelgrid1[sspindx[jj]].mstar = outmstar
             modelgrid1[sspindx[jj]].sfr = outsfr
             modelgrid1[sspindx[jj]].sfr100 = outsfr100
             modelgrid1[sspindx[jj]].sfrage = outsfrage
             modelgrid1[sspindx[jj]].nlyc = outnlyc
             modelgrid1[sspindx[jj]].flux = outflux+nebflux
          endfor 
       endfor    ; close SSP loop 
       if imini eq 0 then modelgrid = temporary(modelgrid1) else $
         modelgrid = [temporary(modelgrid),temporary(modelgrid1)]
    endfor        ; close miniChunk loop   

return, modelgrid
end    

pro isedfit_montegrids, isedfit_paramfile, params=params, thissfhgrid=thissfhgrid, $
  montegrids_dir=montegrids_dir, chunksize=chunksize, minichunksize=minichunksize, $
  clobber=clobber, debug=debug

; read the SFHGRID parameter file
    if n_elements(isedfit_paramfile) eq 0 and n_elements(params) eq 0 then begin
       doc_library, 'isedfit_montegrids'
       return
    endif

    if n_elements(montegrids_dir) eq 0 then montegrids_dir = get_pwd()+'montegrids/'

; read the parameter file and then optionally call this routine
; recursively     
    if (n_elements(params) eq 0) then params = $
      read_isedfit_paramfile(isedfit_paramfile,thissfhgrid=thissfhgrid)

    ngrid = n_elements(params)
    if ngrid gt 1 then begin
       for ii = 0, ngrid-1 do begin
          isedfit_montegrids, params=params[ii], montegrids_dir=montegrids_dir, $
            chunksize=chunksize, minichunksize=minichunksize, clobber=clobber, $
            debug=debug
       endfor
       return
    endif

; number of models per output FITS table    
    if n_elements(chunksize) eq 0 then chunksize = 5000L
    if n_elements(minichunksize) eq 0 then minichunksize = 500L
    
; read the SSP information structure    
    ssppath = getenv('ISEDFIT_SSP_DIR')
    if file_test(ssppath,/dir) eq 0 then begin
       splog, 'Verify that ${ISEDFIT_SSP_DIR} environment variable is defined!'
       return
    endif
    ssppath=ssppath+'/'
    sspinfofile = ssppath+'info_'+strtrim(params.synthmodels,2)+'_'+strtrim(params.imf,2)+'.fits.gz'
    if file_test(sspinfofile) eq 0 then begin
       splog, 'SSP info file '+sspinfofile+' not found!'
       splog, 'Run the appropriate BUILD_*_SSP code!'
       return
    endif
    sspinfo = mrdfits(sspinfofile,1,/silent)

; make directories and delete old files
    sfhgridstring = 'sfhgrid'+string(params.sfhgrid,format='(I2.2)')
    sfhgridpath = montegrids_dir+sfhgridstring+$
      '/'+strtrim(params.synthmodels,2)+'/'+strtrim(params.redcurve,2)+'/'

    if file_test(sfhgridpath,/dir) eq 0 then begin
       splog, 'Making directory '+sfhgridpath
       file_mkdir, sfhgridpath
    endif

; ---------------------------------------------------------------------------
; first major step: build the Monte Carlo grid, if it doesn't
; already exist 
    montefile = sfhgridpath+strtrim(params.imf,2)+'_montegrid.fits'
    if (file_test(montefile+'.gz') eq 0) or keyword_set(clobber) then begin
       cc = 'N'
       if (keyword_set(clobber) eq 0) then begin
          splog, 'Delete all *'+strtrim(params.imf,2)+'* files from '+sfhgridpath+' [Y/N]?'
          cc = get_kbrd(1)
       endif
       if keyword_set(clobber) or (strupcase(cc) eq 'Y') then begin
          delfiles = file_search(sfhgridpath+'*'+strtrim(params.imf,2)+'*.fits.gz',count=ndel)
          if ndel ne 0 then file_delete, delfiles, /quiet
       endif
       
       splog, 'Building SFHGRID='+sfhgridstring+' REDCURVE='+strtrim(params.redcurve,2)+$
         ' NMODEL='+string(params.nmodel,format='(I0)')

       montegrid = init_montegrid(params.nmodel,nmaxburst=params.nmaxburst)

; draw uniformly from linear TAU, or 1/TAU?
       tau = randomu(seed,params.nmodel)*(params.tau[1]-params.tau[0])+params.tau[0]
       if params.oneovertau eq 1 and params.delayed eq 1 then $
         message, 'DELAYED and ONEOVERTAU may not work well together.'
       if params.oneovertau eq 1 then tau = 1D/tau
       montegrid.tau = tau

; metallicity; check to make sure that the prior boundaries do not
; exceed the metallicity range available from the chosen SYNTHMODELS
       if (params.Zmetal[0] lt min(sspinfo.Zmetal)) then begin
          splog, 'Adjusting minimum prior metallicity!'
          params.Zmetal[0] = min(sspinfo.Zmetal)
       endif
       if (params.Zmetal[1] gt max(sspinfo.Zmetal)) then begin
          splog, 'Adjusting maximum prior metallicity!'
          params.Zmetal[1] = max(sspinfo.Zmetal)
       endif
       montegrid.Zmetal = randomu(seed,params.nmodel)*(params.Zmetal[1]-params.Zmetal[0])+params.Zmetal[0]

; age; unfortunately I think we have to loop to sort
       montegrid.age = randomu(seed,params.nmodel)*(params.age[1]-params.age[0])+params.age[0]

; reddening, if any; Gamma distribution is the default, unless FLATAV==1
       if (params.av[1] gt 0) then begin
          if params.flatav then begin
             montegrid.av = randomu(seed,params.nmodel)*(params.av[1]-params.av[0])+params.av[0] 
          endif else begin
             montegrid.av = params.av[0]*randomu(seed,params.nmodel,gamma=params.av[1])
;            montegrid.av = ((10.0^(randomn(seed,params.nmodel)*params.av[1]+alog10(params.av[0]))>0.0
;            montegrid.av = randomu(seed,params.nmodel,gamma=1.0)*params.av[1]+params.av[0]
;            montegrid.av = 10.0^(randomn(seed,params.nmodel)*params.av[1]+alog10(params.av[0]))
          endelse
       
; "mu" is the Charlot & Fall (2000) factor for evolved stellar
; populations; Gamma distribution is the default, unless FLATMU==1
          if (strtrim(params.redcurve,2) eq 'charlot') then begin
             if params.flatmu then begin
                montegrid.mu = randomu(seed,params.nmodel)*(params.mu[1]-params.mu[0])+params.mu[0] 
             endif else begin
                montegrid.mu = params.mu[0]*randomu(seed,params.nmodel,gamma=params.mu[1])
;               montegrid.mu = ((10.0^(randomn(seed,params.nmodel)*params.mu[1]+alog10(params.mu[0])))<1.0)>0.0
             endelse
          endif
       endif

; add emission lines; draw [OIII]/H-beta from a uniform logarithmic
; distribution 
       if params.nebular then begin
          montegrid.oiiihb = randomu(seed,params.nmodel)*$
            (params.oiiihb[1]-params.oiiihb[0])+params.oiiihb[0] 
       endif

; now assign bursts; note that the bursts can occur outside
; (generally, before) the AGE vector; divide the time vector into
; NMAXBURST intervals of width INTERVAL_PBURST [Gyr]
       ntime = 100
;      ntime = params.nage
       
; type of burst: 0 (step function, default), 1 (gaussian), 2 (step
; function with exponential wings) 
       if (params.nmaxburst gt 0) then begin
          tburst = dblarr(params.nmaxburst,params.nmodel)-1.0
          ran = randomu(seed,params.nmaxburst,params.nmodel)
          for imod = 0L, params.nmodel-1 do begin
             for ib = 0, params.nmaxburst-1 do begin
                tmin = params.tburst[0]+ib*params.interval_pburst
                tmax = (params.tburst[0]+(ib+1)*params.interval_pburst)<params.tburst[1]
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
                prob = params.pburst*total([[dtimeleft],[dtimeright]],2)/$
                  params.interval_pburst
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
             dtburst = randomu(seed,nallburst)*(params.dtburst[1]-$
               params.dtburst[0])+params.dtburst[0]
          endif else begin
             dtburst = 10.0^(randomu(seed,nallburst)*(alog10(params.dtburst[1])-$
               alog10(params.dtburst[0]))+alog10(params.dtburst[0]))
          endelse

          if params.flatfburst then begin
             fburst = randomu(seed,nallburst)*(params.fburst[1]-$
               params.fburst[0])+params.fburst[0]
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
             if (params.bursttype ne 1) then message, 'Should use truncated *Gaussian* bursts.'
;            montegrid[hasburst[trunc]].trunctau = params.trunctau ; [Gyr]
             montegrid[hasburst[trunc]].trunctau = $
               randomu(seed,ntrunc)*(params.trunctau[1]-$
               params.trunctau[0])+params.trunctau[0] ; [Gyr]
          endif
       endif  

; write out
       im_mwrfits, montegrid, montefile, /clobber
    endif else begin
; read the Monte Carlo grid
       splog, 'Reading '+montefile+'.gz'
       montegrid = mrdfits(montefile+'.gz',1)
       if (n_elements(montegrid) ne params.nmodel) then begin
          splog, 'Dimensions of MONTEGRID and PARAMS.NMODEL '+$
            'do not match; call with /CLOBBER!'
          return
       endif
    endelse 

; ---------------------------------------------------------------------------
; now build the actual composite stellar populations; do it in chunks
; to avoid memory issues and pass along all the parameters
    chunkinfofile = sfhgridpath+strtrim(params.imf,2)+'_chunkinfo.fits'
    nchunk = ceil(params.nmodel/float(chunksize))

    chunkfiles = strtrim(params.imf,2)+'_chunk_'+string(lindgen(nchunk)+1,$
      format='(I4.4)')+'.fits'
    chunkinfo = create_struct(params, 'chunksize', chunksize, $
      'nchunk', nchunk, 'chunkfiles', sfhgridpath+chunkfiles)

; loop on each chunk of models
    t0 = systime(1)
    for ichunk = 0, nchunk-1 do begin
       i1 = ichunk*chunksize
       i2 = ((ichunk*chunksize+chunksize)<params.nmodel)-1L
       nthese = i2-i1+1L
       these = lindgen(nthese)+i1

       modelgrid = build_modelgrid(montegrid[these],params=params,debug=debug,$
         sspinfo=sspinfo,ssppath=ssppath+strtrim(params.synthmodels,2)+'/',$
         minichunksize=minichunksize,ichunk=ichunk,nchunk=nchunk)
       modelgrid.chunkindx = ichunk

       im_mwrfits, modelgrid, chunkinfo.chunkfiles[ichunk], /clobber
    endfor
    splog, 'Total time (min) = ', (systime(1)-t0)/60.0
 
; write out the chunkinfofile
    chunkinfo.chunkfiles = file_basename(chunkinfo.chunkfiles)+'.gz'
    im_mwrfits, chunkinfo, chunkinfofile, /clobber

;; build a QAplot: age, tau, zmetal, av, mu, oiiihb, sfr, t/tau, sfrage
;    qafile = sfhgridpath+'qa_'+strtrim(params.imf,2)+'_montegrid.ps'
;    im_plotconfig, 0, pos, psfile=qafile
;    im_plotconfig, psfile=qafile, /psclose, /pdf
;stop       
    
    
return
end
