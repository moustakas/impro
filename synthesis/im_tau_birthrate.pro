;+
; NAME:
;       IM_TAU_BIRTHRATE
;
; PURPOSE:
;       Compute the birthrate parameter for various continuous (tau)
;       SFH models. 
;
; CALLING SEQUENCE:
;       im_tau_birthrate, /compute
;
; INPUTS:
;       None required.
;
; OPTIONAL INPUTS:
;       None required.
;
; KEYWORD PARAMETERS:
;       compute - if the TAGE or TAU parameters are changed then set
;                 COMPUTE=1 to recompute the models and write the
;                 results out; otherwise, read in the results from the
;                 IMPRO_DIR+BC03 directory
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;
; PROCEDURES USED:
;       MRDFITS, MWRFITS, IM_READ_BC03(), DJS_PLOT, DJS_OPLOT
;
; COMMENTS:
;
; EXAMPLES:
;
; MODIFICATION HISTORY:
;       J. Moustakas, 2005 Apr 11, U of A - written
;-

pro im_tau_birthrate, compute=compute

    rootpath = getenv('CATALOGS_DIR')+'/sfhgrid/' ; parameter I/O path
    modelspath = rootpath+'basemodels/'           ; base models output path
    outpath = getenv('IMPRO_DIR')+'/bc03/'        ; output path
    outfile = 'bc03_tau_birthrate.fits'

    sfrhaconst = 7.9D-42        ; [(M_sun/yr)/(erg/s)]
    lsun = 3.826D33             ; [erg/s]
    hahb = 2.86
;   hahb = return_tbalmer()

    if keyword_set(compute) then begin
       
; read the basemodels information structure    
       
       imfstr = 'salp'          ; IMF
       Zstr = 'm62'             ; metallicity

       sfhfile = 'sfh_base_models_'+Zstr+'_'+imfstr+'.fits.gz'
       sfhgrid = mrdfits(rootpath+sfhfile,1,/silent)

; select the TAU values of interest

       if (n_elements(mintau) eq 0L) then mintau = 2.0
       if (n_elements(maxtau) eq 0L) then maxtau = 10.0
       if (n_elements(dtau) eq 0L) then dtau = 2.0

       ntau = fix((maxtau-mintau)/dtau)+1L
;      dtau = (maxtau-mintau)/(ntau-1.0)
       tau = findgen(ntau)*dtau+mintau

       tau = [1.0,2,3,4,8,20]
       ntau = n_elements(tau)
       
       reftau = sfhgrid[0].tau
       match, reftau, tau, taukeep, tauindx
       if (n_elements(tauindx) ne n_elements(tau)) then begin
          splog, 'There was a problem with the TAU values requested.'
          return
       endif

       taufiles = sfhgrid.contfile[taukeep,*]
       
; define a logarithmic age distribution of the models

       if (n_elements(mintage) eq 0L) then mintage = 0.1D ; [Gyr]
       if (n_elements(maxtage) eq 0L) then maxtage = 13.0D ; [Gyr]
       if (n_elements(dlogtage) eq 0L) then dlogtage = 0.1D ; logarithmic age spacing

       ntage = round(1.0D + (alog10(maxtage)-alog10(mintage))/dlogtage)
       logtage = alog10(mintage) + dlogtage*findgen(ntage)
       tage = 10.0D^logtage

; initialize the output data structure

       result = {$
         tau:         0.0,$
         ntau:       ntau,$
         tage:        0.0,$
         ntage:     ntage,$
         sfr:         0.0,$
         sfr_average: 0.0,$
         fraction:    0.0,$
         birthrate:   0.0,$
         ml_b:        0.0,$
         ewha:        0.0,$
         ewhb:        0.0,$
         LHa:         0.0,$
         LHb:         0.0 $
         }
       result = replicate(result,ntau,ntage)

       result.tau  = rebin(tau,ntau,ntage)
       result.tage = transpose(rebin(tage,ntage,ntau))
       
       for itau = 0L, ntau-1L do begin

          print, format='("Tau=",I0,"/",I0,".   ",A5,$)', itau+1L, ntau, string(13b)
          
; read the SSP       

          ssp = im_read_bc03(isedfile=taufiles[itau],isedpath=modelspath,$
            age=tage,minwave=minwave,maxwave=maxwave,bc=sspextras,/silent,$
            /salpeter)

; compute the birthrate parameter

          massnorm = sspextras.mgalaxy
;         sspextras.mgalaxy = sspextras.mgalaxy/massnorm
stop
          sspextras.mgas = sspextras.mgalaxy-sspextras.m_
          
          fraction = sspextras.mgas/sspextras.mgalaxy
          sfr_average = sspextras.m_/tage/1E9

          result[itau,*].sfr = reform(sspextras.sfr_yr,1,ntage)
          result[itau,*].sfr_average = reform(sfr_average,1,ntage)
          result[itau,*].fraction = reform(fraction,1,ntage)

          result[itau,*].birthrate = reform(sspextras.sfr_yr*(1-fraction)/sfr_average,1,ntage)
          result[itau,*].ml_b = reform(sspextras.m__lb,1,ntage)

; H-alpha

          result[itau,*].LHa = result[itau,*].sfr/sfrhaconst
          result[itau,*].LHb = result[itau,*].LHa/hahb

          for iage = 0L, ntage-1L do begin
          
             balmerabs = ibalmerabs(ssp.wave,ssp.flux[*,iage],$
               ferr=ssp.flux[*,iage]*0.1,balmersigma=balmersigma,$
               balmerwaves=balmerwaves,debug=debug)

             result[itau,iage].ewhb = result[itau,iage].LHb/(lsun*balmerabs[2].babs_continuum) ; [Angstrom]
             result[itau,iage].ewha = result[itau,iage].LHa/(lsun*balmerabs[3].babs_continuum) ; [Angstrom]

          endfor
          
       endfor

       struct_print, result[*,20]
       print
       struct_print, result[5,*]

stop
       
; write out
       
       splog, 'Writing '+outpath+outfile+'.'
       mwrfits, result, outpath+outfile, /create
       return
       
    endif else begin

; read in        
       
       splog, 'Reading '+outpath+outfile+'.'
       result = mrdfits(outpath+outfile,1,/silent)
       result = reform(result,result[0,0].ntau,result[0,0].ntage)

       tau = reform(result[*,0].tau)
       tage = reform(result[0,*].tage)
       
       ntau = n_elements(tau)
       ntage = n_elements(tage)
       
    endelse

    im_window, 0, /square, xratio=0.6
    
    ncols = 2.0
    pagemaker, nx=ncols, ny=ceil(ntau/ncols), xspace=0.0, yspace=0.0, $
      xmargin=[1.2,0.2], ymargin=[0.3,1.2], position=pos, /normal

; birthrate versus age for various tau models
    
    xrange = minmax(tage)
    yrange = minmax(result.birthrate)*[1.0,1.5]
    
    for itau = 0L, ntau-1L do begin

       if odd(itau) then begin
          ytickname = replicate(' ',10)
          ytitle = ''
       endif else begin
          delvarx, ytickname
          ytitle = 'Birthrate'
          ytitle2 = 'Fraction'
       endelse

       if (itau lt ntau-2L) then begin
          xtickname = replicate(' ',10)
          xtitle = ''
       endif else begin
          delvarx, xtickname
          xtitle = 'log Age [Gyr]'
       endelse
       
       djs_plot, [0], [0], /nodata, xrange=xrange, yrange=yrange, $
         xsty=3, ysty=3, xthick=2.0, ythick=2.0, charthick=2.0, $
         charsize=1.5, noerase=(itau gt 0L), position=pos[*,itau], $
         ytickname=ytickname, xtickname=xtickname, xtitle=xtitle, $
         ytitle=ytitle, /ylog
       djs_oplot, tage, result[itau,*].birthrate, thick=2.0, ps=-4

    endfor
       
    cc = get_kbrd(1)
       
    yrange = minmax(1.0-result.fraction)*[1.0,1.1]
    for itau = 0L, ntau-1L do begin

       if odd(itau) then begin
          ytickname = replicate(' ',10)
          ytitle = ''
       endif else begin
          delvarx, ytickname
          ytitle = '(1-R)'
       endelse

       if (itau lt ntau-2L) then begin
          xtickname = replicate(' ',10)
          xtitle = ''
       endif else begin
          delvarx, xtickname
          xtitle = 'log Age [Gyr]'
       endelse
       
       djs_plot, [0], [0], /nodata, xrange=xrange, yrange=yrange, $
         xsty=3, ysty=3, xthick=2.0, ythick=2.0, charthick=2.0, $
         charsize=1.5, noerase=(itau gt 0L), position=pos[*,itau], $
         ytickname=ytickname, xtickname=xtickname, xtitle=xtitle, $
         ytitle=ytitle, ylog=0
       djs_oplot, tage, 1.0-result[itau,*].fraction, thick=2.0, ps=-4

    endfor    
    cc = get_kbrd(1)
       
    yrange = minmax(result.sfr)*[1.0,1.1]
    for itau = 0L, ntau-1L do begin

       if odd(itau) then begin
          ytickname = replicate(' ',10)
          ytitle = ''
       endif else begin
          delvarx, ytickname
          ytitle = 'SFR'
       endelse

       if (itau lt ntau-2L) then begin
          xtickname = replicate(' ',10)
          xtitle = ''
       endif else begin
          delvarx, xtickname
          xtitle = 'log Age [Gyr]'
       endelse
       
       djs_plot, [0], [0], /nodata, xrange=xrange, yrange=yrange, $
         xsty=3, ysty=3, xthick=2.0, ythick=2.0, charthick=2.0, $
         charsize=1.5, noerase=(itau gt 0L), position=pos[*,itau], $
         ytickname=ytickname, xtickname=xtickname, xtitle=xtitle, $
         ytitle=ytitle, ylog=1
       djs_oplot, tage, result[itau,*].sfr, thick=2.0, ps=-4

    endfor    
    cc = get_kbrd(1)
       
stop
    
return
end
    
