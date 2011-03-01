function measure_ires, spec, wave, skyspec, header
; jm02jan31uofa
;; measure the instrumental resolution from the sky spectrum
; measure the instrumental resolution of the data as a function of
; wavelength

    cspeed = 2.99792458D5

   if not keyword_set(datapath) then datapath = cwd()

;   if (n_elements(background) EQ 0) then background = 0
    
;   pushd, '/home/ioannis/kennicutt/data/98mar/'

    arcfile = 'wdra.0104.fits'
    arc = rd2dspec(arcfile,datapath=datapath)

    ncols = arc.naxis1
    nrows = arc.naxis2
    midrow = nrows/2L

    header = *arc.header
    wave = make_wave(header,cd1_1=dwave,crval1=minwave)
    maxwave = max(wave)    

; fit the middle row

;   flux = arc.image[*,midrow]
    nmed = 10L
    flux = djs_median(arc.image[*,midrow-nmed/2L:midrow+nmed/2L],2)
    ferr = arc.sigmamap[*,midrow]
    norm = max(flux)
    flux = flux/max(flux) & ferr = ferr/norm
    invvar = 1.0/ferr^2.0
    zobj = 0.0
    
;   nmed = 10L
;   medspec = djs_median(arc.image[*,midrow-nmed/2L:midrow+nmed/2L],2) ; median filter NMED rows centered on MIDROW
;   medspec = medspec/max(medspec)
;   plot, wave, medspec, xsty=3, ysty=3

    lampwave = [3850.581,3888.648,4158.591,4510.0,4545.052,4657.901,$
                4764.865,4806.021,4847.810,4921.931,4965.080,5015.680,$
                5187.746,5495.874,5739.520,6416.307];,6752.834,6871.289]
    nline = n_elements(lampwave)
    
; ----------------------------------------------------------------------
;; read in the arc lines to fit
;    
;    lampfile = filepath('',root_dir=getenv('RKSPEC_DIR'),subdirectory='etc')+'lamphear.dat'
;    readcol, lampfile, lampwave, intensity, format='F,F', /silent
;
;    goodwave = where((lampwave gt minwave+100.0) and (lampwave lt maxwave-100.0),nline)
;    lampwave = lampwave[goodwave]
;    intensity = intensity[goodwave]
;
;; sort by line strength, choose the NLINE strongest lines, and sort
;; again
;
;    srtinten = sort(intensity)
;    lampwave = lampwave[srtinten]
;
;    nline = 10L
;    lampwave = lampwave[0:nline-1L]
;    lampwave = lampwave[sort(lampwave)]    
; ----------------------------------------------------------------------
    
; LINEBACKFIT wants vacuum wavelengths
    
    linevacwave = lampwave
    airtovac, linevacwave
    vacwave = wave
    airtovac, vacwave    
    
; background terms (constant)
    
    nback = 2L
    background = poly_array(ncols,nback)

; constraints: same redshift, positive Gaussians

    zindex = lonarr(nline)
    findex = lindgen(nline)
    fvalue = replicate(1.0,nline)    
    
    result = linebackfit(linevacwave,alog10(vacwave),flux,invvar=invvar, $
                         zguess=zobj,background=background,zindex=zindex,$
                         findex=findex,fvalue=fvalue,bterms=bterms,$
                         yfit=yfit,bfit=bfit)

    djs_plot, vacwave, flux, xsty=3, ysty=3, ps=10
    djs_oplot, vacwave, yfit, ps=10, color='green'

    print_struct, result

    good = where(result.linenpix ne 0L)
    result = result[good]
    lampwave = lampwave[good]
    
    dv = result.linesigma
    dv_err = result.linesigma_err
    dlambda = lampwave*dv/cspeed
    dlambda_err = lampwave*dv_err/cspeed

    res = lampwave/dlambda ; resolution
    
;   yrange = minmax(dlambda)
    yrange = minmax(dv)
    plot, [minwave,maxwave], yrange, xsty=3, ysty=3, yrange=yrange, /nodata
;   oploterror, lampwave, dlambda, dlambda_err, ps=4
    oploterror, lampwave, dv, dv_err, ps=4


    coeff = linfit(lampwave,dlambda)
;   coeff = linfit(lampwave,dv)
    resfit = poly(wave,coeff)
    oplot, wave, resfit, thick=2.0, line=0

;   dl_5500_nominal = interpol(5500.0/res,lampwave,5500.0)

    info = create_struct('dv', dv, 'dv_err', dv_err, 'dlambda', dlambda, $
                         'dlambda_err', dlambda_err, 'wave', wave, 'flux', $
                         flux, 'fit', yfit)

stop

;; fit the continuum
;
;    cont = dkboxstats(medspec,xwidth=50,boxstat='min')
;    contcoeff = poly_fit(wave,cont,2)
;    contfit = poly(wave,contcoeff)
;    medspec = medspec-contfit
;    
;; initialize the output structure
;
;    ires = create_struct(name='Instrumental Resolution', $
;                         'linewave'	, 0.0, $
;                         'linesigma'	, 0.0, $
;                         'lineflux'	, 0.0)
;    ires = replicate(ires,nline)
;
;; initialize the fitting parameters
;
;    parinfo = {value: 0.0D,     $
;               fixed: 0L,       $
;               limited: [0,0],  $
;               tied: ' ',       $
;               limits: [0.0D,0D]}
;    parinfo = replicate(parinfo,3*nline+nback)
;
;    functargs = {nline: nline, nback: nback, $
;                 loglam: wave, background: background}
;
;    sigguess = 4*dwave/2.35 ; Gaussian sigma guess [Angstrom]
;    
;    for i = 0L, nline-1L do begin
;
;       parinfo[0+i*3L].value = 0.5*sqrt(2.0*!pi)*sigguess ; peak flux 
;       parinfo[1+i*3L].value = lampwave[i]                 ; wavelength
;       parinfo[2+i*3L].value = sigguess                    ; sigma
;       
;    endfor
;
;    for j = 0L, nback-1L do parinfo[nline*3+j].value = 1.0 
;    
;    xindx = lindgen(ncols)
;    
;;   test = manygauss(xindx,parinfo.value,nline=nline,loglam=wave,nback=0,background=background)
;
;stop
;    
;    lfit = mpfitfun('manygauss', xindx, medspec, parinfo=parinfo, functargs=functargs, $
;                    perror=perror, yfit=yfit, nfev=nfev, covar=covar, $
;                    niter=niter, status=status)
;
;stop
;    
;    splog, 'MPFIT number of function evaluations=', nfev
;    splog, 'MPFIT number of iterations=', niter
;    splog, 'MPFIT exit status=', status
;    if (status EQ 5) then $
;      splog, 'Warning: Maximum number of iterations reached: ', niter
;    yfit[igood] = yfit1
;
;
;    
;    
;    cspeed = 2.99792458e5 ; [km/s]
;
;    skywaves = [$
;                 4046.563, $    ;   Hg I
;;                4358.3277, $   ;  Hg I
;;                5199.0, $      ; NI
;;                5460.7348, $   ; Hg I
;                 5577.339, $    ; O I
;;                5682.633, $    ; NaI
;;                5889.950, $    ; Na D (broadened) 5889.950 + 5895.924
;                 6300.32 $     ; O I blended with OH
;;                6360.2272, $   ; OI 5-0 P1
;;                6360.7747 $    ; OI 5-0 P1
;               ]              
;
;    path = '/home/ioannis/kennicutt/atlas1d/'
;    spec = 'n4736_nuc_5.ms.fits'
;
;    scube = rd1dspec(spec,datapath=path)
;    header = *scube.header
;
;    flux = scube.spec
;    wave = scube.wave
;    sky = scube.sky
;
;    goodwaves = where((skywaves gt min(wave)) and (skywaves lt max(wave)),nskywaves)
;    skywaves = skywaves[goodwaves]
;    
;    box = 20.0 ; [Angstrom]
;
;    skylam = fltarr(nskywaves)
;    siglam = fltarr(nskywaves)
;    stay = ''
;    for i = 0L, nskywaves-1L do begin
;    
;       local = where((wave gt skywaves[i]-box) and (wave lt skywaves[i]+box),nlocal)
;       skyfit = mpfitpeak(wave[local],sky[local],a,nterm=4,/gaussian,/positive)
;
;       skylam[i] = a[1]
;       siglam[i] = a[2]
;       
;;      djs_plot, wave[local], sky[local], ps=10, xsty=3, ysty=3
;;      djs_oplot, wave[local], skyfit, color='green'
;;      read, stay
;       
;    endfor
;
;    sigvel = cspeed*siglam/skylam
;
;    lamcoeff = linfit(skylam,siglam)
;    siglamfit = poly(wave,lamcoeff)
;    
;    window, 2, xs=450, ys=450
;    plot, skylam, siglam, ps=4, xsty=3, ysty=3
;    oplot, wave, siglamfit, line=0, thick=2.0
    
return, info
end
