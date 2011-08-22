;+
; NAME:
;   BC03_MASSNORM_TEST
;
; PURPOSE:
;   Build a QAplot to determine whether the BC03 templates need to be
;   divided by the mstar in the .4color file.
;
; INPUTS: 
;
; KEYWORD PARAMETERS: 
;   makemodels - build the BC03 models we need
;
; OUTPUTS: 
;   The models and QAplot are written to the ${bc03_dir}/massnorm_test
;   directory. 
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2004 Nov 16, U of A
;
; Copyright (C) 2004, John Moustakas
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

pro bc03_massnorm_test, makemodels=makemodels

    bc03path = getenv('bc03_dir')+'/src/'         ; FORTRAN source path
    isedpath = getenv('bc03_dir')+'/models/Padova1994/salpeter/'
    modelspath = getenv('bc03_dir')+'/massnorm_test/'
    if (file_test(modelspath) eq 0) then spawn, 'mkdir -p '+modelspath, /sh
    
    if keyword_set(makemodels) then begin

; tau = 3 Gyr
       input = 'tau_sfr.input'
       output = 'tau_sfr.output'
       isedfile = 'tau_sfr_m62_salp.ised'
       
       openw, lun, modelspath+input, /get_lun
       printf, lun, 'bc2003_hr_m62_salp_ssp.ised'
       printf, lun, 'N'       ; include dust attenuation?
       printf, lun, '1'       ; SFH law [exponential]
       printf, lun, '3.0'     ; e-folding time
       printf, lun, 'N'       ; include gas recycling?
       printf, lun, '20.0'    ; SFR = 0 at t = 20 Gyr
       printf, lun, isedfile  ; output file
       free_lun, lun          ; close the csp_burst input file

; run csp_galaxev, write the resultant models into the output
; directory as a binary FITS file, and delete the ISED files

       pushd, isedpath
       splog, 'Generating '+isedfile
       t0 = systime(1)
       spawn, bc03path+'csp_galaxev < '+modelspath+input+' > '+modelspath+output, /sh
       splog, format='("Total time = ",G0," minutes.")', (systime(1)-t0)/60.0
       popd

       spawn, '/bin/mv -f '+isedpath+repstr(isedfile,'.ised','.*')+' '+modelspath, /sh

; continuous star formation

       input = 'const_sfr.input'
       output = 'const_sfr.output'
       isedfile = 'const_sfr_m62_salp.ised'
       
       openw, lun, modelspath+input, /get_lun
       printf, lun, 'bc2003_hr_m62_salp_ssp.ised'
       printf, lun, 'N'       ; include dust attenuation?
       printf, lun, '2'       ; SFH law [truncated burst]
       printf, lun, '20.0'    ; burst duration [Gyr]
       printf, lun, isedfile  ; output file
       free_lun, lun          ; close the csp_burst input file

; run csp_galaxev, write the resultant models into the output
; directory as a binary FITS file, and delete the ISED files

       pushd, isedpath
       splog, 'Generating '+isedfile
       t0 = systime(1)
       spawn, bc03path+'csp_galaxev < '+modelspath+input+' > '+modelspath+output, /sh
       splog, format='("Total time = ",G0," minutes.")', (systime(1)-t0)/60.0
       popd

       spawn, '/bin/mv -f '+isedpath+repstr(isedfile,'.ised','.*')+' '+modelspath, /sh

    endif

    dfpsplot, modelspath+'bc03_massnorm_test.ps', /square
    postthick = 4.0
    
    ssp = im_read_bc03(age=[0.01,0.03,0.1,0.3,1,2,5,10],bc=sspe)
    
    djs_plot, ssp.wave, ssp.flux[*,0], /xlog, /ylog, xr=[500,1E4], xsty=3, ysty=3, $
      yrange=[1D-7,1], charsize=1.8, charthick=postthick, xthick=postthick, $
      ythick=postthick, ytitle='f_{\lambda} [erg/s/M'+sunsymbol()+']', $
      xtitle='Wavelength [\AA]', xmargin=[10,5], $
      title='SSP BC03 Models'
    xyouts, 10^!x.crange[1]*1.05, interpol(ssp.flux[*,0],ssp.wave,10^!x.crange[1]), $
      strtrim(string(ssp.age[0]/1E9,format='(F5.2)'),2)+' Gyr', align=0, charsize=1.0, $
      charthick=postthick, /data
    for i = 1L, 7L do djs_oplot, ssp.wave, ssp.flux[*,i]
    for i = 1L, 7L do xyouts, 10^!x.crange[1]*1.05, interpol(ssp.flux[*,i],$
      ssp.wave,10^!x.crange[1]), strtrim(string(ssp.age[i]/1E9,format='(F5.2)'),2), $
      align=0, charsize=1.0, charthick=postthick, /data

; continuous SFR    
    isedpath = modelspath
    const = im_read_bc03(isedfile='const_sfr_m62_salp.ised',isedpath=isedpath,$
      age=[0.01,0.03,0.1,0.3,1,2,5,10],bc=conste,/salpeter)

; not normalized
    
    djs_plot, const.wave, const.flux[*,0], /xlog, /ylog, xr=[500,1E4], xsty=3, ysty=3, $
      yrange=[3E-7,6E-4], charsize=1.8, charthick=postthick, xthick=postthick, $
      ythick=postthick, ytitle='f_{\lambda} [erg/s/M'+sunsymbol()+']', $
      xtitle='Wavelength [\AA]', xmargin=[10,5], $
      title='Continuous SFR - Not normalized'
    xyouts, 10^!x.crange[1]*1.05, interpol(const.flux[*,0],const.wave,10^!x.crange[1]), $
      strtrim(string(const.age[0]/1E9,format='(F5.2)'),2)+' Gyr', align=0, charsize=1.0, $
      charthick=postthick, /data
    for i = 1L, 7L do djs_oplot, const.wave, const.flux[*,i]
    for i = 1L, 7L do xyouts, 10^!x.crange[1]*1.05, interpol(const.flux[*,i],$
      const.wave,10^!x.crange[1]), strtrim(string(const.age[i]/1E9,format='(F5.2)'),2), $
      align=0, charsize=1.0, charthick=postthick, /data

; normalized
    djs_plot, const.wave, const.flux[*,0]/conste[0].m_, /xlog, /ylog, xr=[500,1E4], xsty=3, ysty=3, $
      yrange=[1D-5,1], charsize=1.8, charthick=postthick, xthick=postthick, $
      ythick=postthick, ytitle='f_{\lambda} [erg/s/M'+sunsymbol()+']', $
      xtitle='Wavelength [\AA]', xmargin=[10,5], $
      title='Continuous SFR - Normalized'
    xyouts, 10^!x.crange[1]*1.05, interpol(const.flux[*,0]/conste[0].m_,const.wave,10^!x.crange[1]), $
      strtrim(string(const.age[0]/1E9,format='(F5.2)'),2)+' Gyr', align=0, charsize=1.0, $
      charthick=postthick, /data
    for i = 1L, 7L do djs_oplot, const.wave, const.flux[*,i]/conste[i].m_
    for i = 1L, 7L do xyouts, 10^!x.crange[1]*1.05, interpol(const.flux[*,i]/conste[i].m_,$
      const.wave,10^!x.crange[1]), strtrim(string(const.age[i]/1E9,format='(F5.2)'),2), $
      align=0, charsize=1.0, charthick=postthick, /data

; tau = 3.0 SFR    
    isedpath = modelspath
    tau = im_read_bc03(isedfile='tau_sfr_m62_salp.ised',isedpath=isedpath,$
      age=[0.01,0.03,0.1,0.3,1,2,5,10],bc=taue,/salpeter)

; not normalized
    
    djs_plot, tau.wave, tau.flux[*,0], /xlog, /ylog, xr=[500,1E4], xsty=3, ysty=3, $
      yrange=[3E-6,2E-3], charsize=1.8, charthick=postthick, xthick=postthick, $
      ythick=postthick, ytitle='f_{\lambda} [erg/s/M'+sunsymbol()+']', $
      xtitle='Wavelength [\AA]', xmargin=[10,5], $
      title='\tau = 3 Gyr - Not normalized'
    xyouts, 10^!x.crange[1]*1.05, interpol(tau.flux[*,0],tau.wave,10^!x.crange[1]), $
      strtrim(string(tau.age[0]/1E9,format='(F5.2)'),2)+' Gyr', align=0, charsize=1.0, $
      charthick=postthick, /data
    for i = 1L, 7L do djs_oplot, tau.wave, tau.flux[*,i]
    for i = 1L, 7L do xyouts, 10^!x.crange[1]*1.05, interpol(tau.flux[*,i],$
      tau.wave,10^!x.crange[1]), strtrim(string(tau.age[i]/1E9,format='(F5.2)'),2), $
      align=0, charsize=1.0, charthick=postthick, /data

; normalized
    djs_plot, tau.wave, tau.flux[*,0]/taue[0].m_, /xlog, /ylog, xr=[500,1E4], xsty=3, ysty=3, $
      yrange=[1D-5,1], charsize=1.8, charthick=postthick, xthick=postthick, $
      ythick=postthick, ytitle='f_{\lambda} [erg/s/M'+sunsymbol()+']', $
      xtitle='Wavelength [\AA]', xmargin=[10,5], $
      title='\tau = 3 Gyr - Normalized'
    xyouts, 10^!x.crange[1]*1.05, interpol(tau.flux[*,0]/taue[0].m_,tau.wave,10^!x.crange[1]), $
      strtrim(string(tau.age[0]/1E9,format='(F5.2)'),2)+' Gyr', align=0, charsize=1.0, $
      charthick=postthick, /data
    for i = 1L, 7L do djs_oplot, tau.wave, tau.flux[*,i]/taue[i].m_
    for i = 1L, 7L do xyouts, 10^!x.crange[1]*1.05, interpol(tau.flux[*,i]/taue[i].m_,$
      tau.wave,10^!x.crange[1]), strtrim(string(tau.age[i]/1E9,format='(F5.2)'),2), $
      align=0, charsize=1.0, charthick=postthick, /data

    dfpsclose

return
end
