pro test_galex_starfinder, starfinder=starfinder, compare=compare, $
  completeness=completeness
; jm10apr24ucsd - test the starfinder point-source photometry code on
; a GALEX/DIS image of the Growth strip

    aegispath = getenv('RESEARCHPATH')+'/projects/aegis/'
    datapath = getenv('RESEARCHPATH')+'/data/galex/AEGIS/'

    band = ['nuv','fuv']
    prefix = ['nd','fd']
    zpt = [20.08,18.82]
    
; ---------------------------------------------------------------------------
; photometer using starfinder
    if keyword_set(starfinder) then begin
       for ii = 0, 0 do begin
;      for ii = 0, 1 do begin
          splog, 'Processing band '+band[ii]
          psf = galex_read_psf(fuv=ii eq 1)

          img = mrdfits(datapath+'GROTH_00-'+prefix[ii]+'-int.fits.gz',0,hdr)
          sky =   mrdfits(datapath+'GROTH_00-'+prefix[ii]+'-skybg.fits.gz')

;         threshold = 0.01
          cat = galex_starfinder(img,sky,psf,hdr=hdr,threshold=threshold,$
            min_correlation=min_correlation,zpt=zpt[ii],modelimg=modelimg,$
            params=params)
          help, params, /str

          suffix = ''
;         suffix = '_test'
;         suffix = '_nodeblend'
;         suffix = '_deblost'
          catfile = datapath+band[ii]+'_catalog'+suffix+'.fits'
          modelfile = datapath+band[ii]+'_modelimg'+suffix+'.fits'
          paramsfile = datapath+band[ii]+'_params'+suffix+'.fits'
          im_mwrfits, params, paramsfile, /clobber
          im_mwrfits, cat, catfile, /clobber
          im_mwrfits, modelimg, modelfile, hdr, /clobber
       endfor
    endif

; ---------------------------------------------------------------------------
; derive the photometric completeness and mean magnitude error;
; consider a uniform set of input magnitudes in the interval [15,30];
; we want a ~noiseless estimate of the completeness in 0.25 mag wide
; bins, or 61 bins total; so, to get good statistics in each bin, say
; 100 stars or S/N=10, we need 6100 total simulated stars;
; however, we don't want the fake stars to significantly alter
; the source density, so add the stars ~1000 at a time in the NUV, and
; ~500 at a time in the FUV
    if keyword_set(completeness) then begin
       minmag = 15.0
       maxmag = 30.0
       nstars = 5000L ; total number of fake stars
       npersim = 1000L ; number of stars per simulation
       nsim = ceil(nstars/float(npersim)) ; total number of simulations

; loop on each NUV,FUV channel
       for ii = 0, 0 do begin
          splog, 'Computing completeness in band '+band[ii]
          psf = galex_read_psf(fuv=ii eq 1)

          base = datapath+'GROTH_00-'+prefix[ii]
          img = mrdfits(base+'-int.fits.gz',0,hdr)
          sky = mrdfits(base+'-skybg.fits.gz')  
          sz = size(img,/dim)

; build up the coordinates     
          nfake = 0
          delvarx, xfake, yfake
          while nfake lt nstars do begin
             xyfake = randomu(seed,2*nstars)
             xfake1 = xyfake[0:nstars-1L]*sz[0]
             yfake1 = xyfake[nstars:2*nstars-1L]*sz[1]
             keep = where(interpolate(img gt 0,xfake1,yfake1) gt 0.0,nkeep)
             if (nkeep ne 0L) then begin
                if (n_elements(xfake) eq 0) then begin
                   xfake = xfake1[keep]
                   yfake = yfake1[keep]
                endif else begin
                   xfake = [xfake,xfake1[keep]]
                   yfake = [yfake,yfake1[keep]]
                endelse
             endif
             nfake = n_elements(xfake)
          endwhile

; final coordinates and magnitudes/fluxes
          xfake = xfake[0L:nstars-1L]
          yfake = yfake[0L:nstars-1L]
          magfake = randomu(seed,nstars)*(maxmag-minmag)+minmag
          fakecat = {id: 0L, ximage: -1.0, yimage: -1.0, $
            ra: -1.0D, dec: -1.0D, flux: -1.0, mag: -1.0, $
            ximage_sim: -1.0, yimage_sim: -1.0, flux_sim: -1.0, $
            mag_sim: -1.0}
          fakecat = replicate(fakecat,nstars)
          fakecat.id = lindgen(nstars)
          fakecat.ximage = xfake
          fakecat.yimage = yfake
          fakecat.flux = 10.0^(-0.4*(magfake-zpt[ii]))
          fakecat.mag = magfake
          extast, hdr, astr
          xy2ad, xfake, yfake, astr, ra, dec
          fakecat.ra = ra
          fakecat.dec = dec

          fakefile = aegispath+'fakecat_original_'+band[ii]+'.fits'
          im_mwrfits, fakecat, fakefile, /clobber
          
; do each simulation          
          for isim = 0, nsim-1 do begin
             indx1 = isim*npersim
             indx2 = ((1+isim)*npersim-1L)<(nstars-1L)
             fakecat1 = fakecat[indx1:indx2]
             sim = image_model(fakecat1.ximage,fakecat1.yimage,$
               fakecat1.flux,sz[0],sz[1],psf,data)

             simimg = img + sim ; simulated image
;            threshold = 0.00228780
             cat = galex_starfinder(simimg,sky,psf,hdr=hdr,zpt=zpt[ii])

             catfile = aegispath+'cat'+strtrim(isim,2)+'_'+band[ii]+'.fits'
             im_mwrfits, cat, catfile, /clobber

; compare the input and output             
             compare_lists, fakecat1.ximage, fakecat1.yimage, $
               cat.ximage, cat.yimage, subscripts_1=subc1, $
               subscripts_2=subc2, max_distance=1.0
             niceprint, fakecat1[subc1].ximage, cat[subc2].ximage, $
               fakecat1[subc1].mag, cat[subc2].mag
             fakecat1[subc1].ximage_sim = cat[subc2].ximage
             fakecat1[subc1].yimage_sim = cat[subc2].yimage
             fakecat1[subc1].flux_sim = cat[subc2].flux
             fakecat1[subc1].mag_sim = cat[subc2].mag
             fakecat[indx1:indx2] = fakecat1

;            ww = where((fakecat1.mag-fakecat1.mag_sim) lt 0.5)
;            im_plothist, fakecat1.mag, bin=0.5, min=minmag, max=magmag, x1, y1, /noplot
;            im_plothist, fakecat1[ww].mag, bin=0.5, min=minmag, max=maxmag, x2, y2, /noplot
;            djs_plot, x1, y2/(y1+(y1 eq 0))*(y1 ne 0), psym=10, $
;              xsty=3, ysty=3, yr=[0,1], xrange[minmag,maxmag], $
;              xtitle='Magnitude', ytitle='Completeness'
          endfor 
; write out
          fakefile = aegispath+'fakedata_'+band[ii]+'.fits'
          im_mwrfits, fakecat, fakefile, /clobber
       endfor ; close bandpass
    endif 
       
; ---------------------------------------------------------------------------
; compare photometric catalogs
    if keyword_set(compare) then begin
       cat1 = mrdfits(aegispath+'nuv_catalog.fits.gz',1)
       ref1 = mrdfits(aegispath+'matchgalexcfhtdeep2_good.fits.gz',1)
;      ref = mrdfits(aegispath+'matchgalexcfhtdeep2.fits.gz',1)
;      good = where(ref.alpha_j2000 gt 0.0 and finite(ref.nuv_mag) and $
;        finite(ref.nuv_db))
;      im_mwrfits, ref[good], 'matchgalexcfhtdeep2_good.fits', /clobber
;      write_ds9_regionfile, ref.alpha_j2000, ref.delta_j2000, $
;        file=aegispath+'ref.reg', symbol='circle', color='green'

       spherematch, ref1.alpha_j2000, ref1.delta_j2000, $
         cat1.ra, cat1.dec, 2.0/3600.0, m1, m2
       ref = ref1[m1]
       cat = cat1[m2]

       iso = where(ref.nn10[0] gt 6.0,comp=crowd)
;      bad = where(ref.nuv_db_err gt 0.3,comp=good)

       psfile = aegispath+'qaplot_galex_aegis.ps'
       im_plotconfig, 6, pos, psfile=psfile, charsize=1.6
       magrange = [15,27]
       levels = [0.5,0.75,0.9]
       slevels = ['0.5','0.75','0.9']
       
; deblended vs pipeline
       hogg_scatterplot, ref.nuv_mag, ref.nuv_db, position=pos[*,0], $
         xrange=magrange, yrange=magrange, xsty=1, ysty=1, $
         xtickname=replicate(' ',10), ytitle='NUV Deblended', $
         cannotation=slevels, levels=levels
       djs_oplot, ref[iso].nuv_mag, ref[iso].nuv_db, psym=6, sym=0.2
       djs_oplot, !x.crange, !y.crange, line=0, color='red'
       
       hogg_scatterplot, ref.nuv_mag, ref.nuv_db-ref.nuv_mag, $
         /noerase, position=pos[*,1], $
         xrange=magrange, yrange=[-1.1,1.1], xsty=1, ysty=1, $
         xtitle='NUV SExtractor', ytitle='Residuals', $
         cannotation=slevels, levels=levels
       djs_oplot, ref[iso].nuv_mag, ref[iso].nuv_db-ref[iso].nuv_mag, $
         psym=6, sym=0.2
       djs_oplot, !x.crange, [0,0], line=0, color='red'
       !p.multi = 0

; deblended vs starfinder
       hogg_scatterplot, cat.mag, ref.nuv_db, position=pos[*,0], $
         xrange=magrange, yrange=magrange, xsty=1, ysty=1, $
         xtickname=replicate(' ',10), ytitle='NUV Deblended', $
         cannotation=slevels, levels=levels
       djs_oplot, cat[iso].mag, ref[iso].nuv_db, psym=6, sym=0.2
       djs_oplot, !x.crange, !y.crange, line=0, color='red'
       
       hogg_scatterplot, cat.mag, ref.nuv_db-cat.mag, $
         /noerase, position=pos[*,1], $
         xrange=magrange, yrange=[-1.1,1.1], xsty=1, ysty=1, $
         xtitle='NUV Starfinder', ytitle='Residuals', $
         cannotation=slevels, levels=levels
       djs_oplot, cat[iso].mag, ref[iso].nuv_db-cat[iso].mag, $
         psym=6, sym=0.2
       djs_oplot, !x.crange, [0,0], line=0, color='red'
       !p.multi = 0

; starfinder vs pipeline
       hogg_scatterplot, ref.nuv_mag, cat.mag, position=pos[*,0], $
         xrange=magrange, yrange=magrange, xsty=1, ysty=1, $
         xtickname=replicate(' ',10), ytitle='NUV Starfinder', $
         cannotation=slevels, levels=levels
       djs_oplot, ref[iso].nuv_mag, cat[iso].mag, psym=6, sym=0.2
       djs_oplot, !x.crange, !y.crange, line=0, color='red'
       
       hogg_scatterplot, ref.nuv_mag, cat.mag-ref.nuv_mag, $
         /noerase, position=pos[*,1], $
         xrange=magrange, yrange=[-1.1,1.1], xsty=1, ysty=1, $
         xtitle='NUV SExtractor', ytitle='Residuals', $
         cannotation=slevels, levels=levels
       djs_oplot, ref[iso].nuv_mag, cat[iso].mag-ref[iso].nuv_mag, $
         psym=6, sym=0.2
       djs_oplot, !x.crange, [0,0], line=0, color='red'

       im_plotconfig, psfile=psfile, /gzip, /psclose
       
    endif
    
stop    
       
return
end
