;+
; NAME:
;       PLOT_BC03_TEMPLATES
;
; PURPOSE:
;       Inspect (plot) the BC03 template spectra.
;
; CALLING SEQUENCE:
;       plot_bc03_templates, /tremonti, /postscript
;
; INPUTS:
;       None
;
; OPTIONAL INPUTS:
;       nperplot   - number of templates per plot
;
; KEYWORD PARAMETERS:
;       postscript - generate multi-page postscript output 
;
; OUTPUTS:
;
; PROCEDURES USED:
;       REPSTR(), MRDFITS(), IM_OPENCLOSE, DJS_PLOT, OPLOT,
;       MINMAX(), LEGEND
;
; COMMENTS:
;
; EXAMPLE:
;
; MODIFICATION HISTORY:
;       J. Moustakas, 2003 Nov 24, U of A
;-

pro plot_bc03_templates, nperplot=nperplot, salpeter=salpeter, postscript=postscript

    if keyword_set(postscript) then begin
       postthick = 5.0 
       postthick2 = 2.0 
    endif else begin
       postthick = 2.0
       postthick2 = 2.0
    endelse
    
    Z = ['Z004','Z02','Z05']
    Zstr = 'Z = '+['0.004','0.02','0.05']
    nZ = n_elements(Z)

    if keyword_set(salpeter) then begin
       splog, 'Selecting SALPETER IMF.'
       imf = 'salpeter'
    endif else begin
       splog, 'Selecting CHABRIER IMF.'
       imf = 'chabrier'
    endelse

    eigendir = filepath('',root_dir=getenv('ISPEC_DIR'),subdirectory='templates')
    eigenfile = 'BC03_'+Z+'_'+imf+'_templates.fits'

    pspath = filepath('',root_dir=getenv('bc03_dir'))
    psname = repstr(eigenfile,'.fits','.eps')

    normwave = 5500.0
    
    loadct, 22, /silent ; 13
    
    for i = 0L, nZ-1L do begin

       if file_test(eigendir+eigenfile[i]) eq 0B then begin

          splog, 'Template file '+eigendir+eigenfile[i]+' does not exist.'

       endif else begin
          
          info = mrdfits(eigendir+eigenfile[i],1,/silent)
          flux = mrdfits(eigendir+eigenfile[i],2,header,/silent)
          wave = make_wave(header)
          ntemplate = info.ntemplate
          age = info.template_age/1E6
          strage = string(round(info.template_age/1E6),format='(I5)')

          if n_elements(nperplot) eq 0L then nperplot = n_elements(age)
    
          newflux = fltarr(5000,ntemplate)
          for it = 0L, ntemplate-1L do newflux[*,it] = $
            frebin(im_normalize(flux[*,it],wave,normwave=normwave),5000)

          wave = frebin(wave,5000)
          flux = newflux
          
;         colors = 255*findgen(ntemplate)/(ntemplate-1)
          npage = ceil(ntemplate/float(nperplot))

          im_openclose, pspath+psname[i], postscript=postscript, /color, $
            /encapsulated, xsize=8.5, ysize=5.0; /square
          
          pagemaker, nx=1, ny=1, xspace=0, yspace=0.0, width=7.3, height=3.8, $
            xmargin=[0.9,0.3], ymargin=[0.3,0.9], xpage=8.5, ypage=5.0, $
            position=pos, /normal
          polyfill, [1,1,0,0,1], [1,0,0,1,1], /normal, color=djs_icolor('black')

          for j = 0L, npage-1L do begin

             s1 = j*nperplot
             s2 = (((j+1L)*nperplot)<ntemplate)-1L

             colors = (255*findgen(s2-s1+1)/(s2-s1))>5.0

             xrange = [3200,7300] ; xrange=minmax(wave)
             get_element, wave, xrange, xx
;            yrange = minmax(flux[xx[0]:xx[1],s1:s2])*[1,1.05]
             yrange = [0,3.9]
             
             plot, [0], [0], xsty=3, ysty=3, charsize=1.5, charthick=postthick, $
               xtitle='Wavelength [\AA]', ytitle='Relative Flux [arbitrary units]', $
               xthick=postthick, ythick=postthick, yrange=yrange, xrange=xrange, $
               position=pos, /noerase, color=djs_icolor('white')
             legend, strage[s1:s2]+' Myr', /right, box=0, charsize=1.6, charthick=postthick, $
               top=(strage[s1] lt 1E3), bottom=(strage[s1] ge 1E3), $
 ;             top=(strage[s1] lt 1E3) or (nperplot le 5L), bottom=(strage[s1] ge 1E3) or (nperplot gt 5L), $
               textcolors=colors
;            legend, Zstr[i], /left, /top, box=0, charsize=2.0, charthick=postthick
;            print, (strage[s1] lt 1E3), (strage[s1] ge 1E3)

             for k = 0L, nperplot-1L do begin

                if (k+s1 le s2) then begin
                   oplot, wave, flux[*,k+s1], ps=10, thick=postthick2, color=colors[k]
                endif
                
             endfor

;            if not keyword_set(postscript) then cc = get_kbrd(1)
             
          endfor

          im_openclose, postscript=postscript, /close

       endelse
          
    endfor 

    loadct, 0, /silent
    
return
end
