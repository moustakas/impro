function gandalf_parse_elinefit, sol, esol, linepars, $
  vsys=vsys, fluxscale=fluxscale
; M) Add to the emission setup structure the flux, amplitude and
; kinematics of each line, and call the result the fit_results
; structure. Add to it also the A/N values, the stellar kinematics,
; and the normalised template weights.  This is to save not only the
; emission-line fitting results but also the conditions under which
; the fit was performed.

    light = 2.99792458D5        ; speed of light1 [km/s]

    if (n_elements(vsys) eq 0) then vsys = 0.0
    if (n_elements(voffset) eq 0) then voffset = 0.0
    if (n_elements(fluxscale) eq 0) then fluxscale = 1.0

; initialize the output data structure
    linename = strtrim(linepars.name,2)
    nline = n_elements(linename)
    for ii = 0, nline-1 do begin
       names = linepars[ii].name+['','_amp','_linez','_sigma',$
         '_continuum','_ew','_limit','_ew_limit','_wave']
       values = [replicate('[0.0,-1.0]',n_elements(names)-3),'-1.0','-1.0','-1.0']
       out1 = mrd_struct(names,values,1)
       if (ii eq 0) then out = out1 else out = struct_addtags(out,out1)
    endfor
    out = create_struct({linename: linename},out)

; fill in the wavelengths for all the lines    
    for ii = 0, nline-1 do begin
       tag_wave = tag_indx(out,linepars[ii].name+'_wave')
       out.(tag_wave) = linepars[ii].lambda ; wavelength
    endfor

; first loop through the lines that were ignored because they were
; outside the wavelength range and set their errors to -2.0 (as
; opposed to -1.0, for undetected lines)
    ignore = where((linepars.action eq 'i'),nignore)
    for ii = 0, nignore-1 do begin
       tag_flux = tag_indx(out,linepars[ignore[ii]].name)
       tag_amp = tag_indx(out,linepars[ignore[ii]].name+'_amp')
       tag_linez = tag_indx(out,linepars[ignore[ii]].name+'_linez')
       tag_sigma = tag_indx(out,linepars[ignore[ii]].name+'_sigma')
       tag_continuum = tag_indx(out,linepars[ignore[ii]].name+'_continuum')
       tag_ew = tag_indx(out,linepars[ignore[ii]].name+'_ew')
       tag_limit = tag_indx(out,linepars[ignore[ii]].name+'_limit')
       tag_ew_limit = tag_indx(out,linepars[ignore[ii]].name+'_ew_limit')

       out.(tag_flux) =  [0.0,-2.0]
       out.(tag_amp) = [0.0,-2.0]
       out.(tag_linez) = [0.0,-2.0]
       out.(tag_sigma) = [0.0,-2.0]
       out.(tag_continuum) = [0.0,-2.0]
       out.(tag_ew) = [0.0,-2.0]
       out.(tag_limit) = -2.0
       out.(tag_ew_limit) = -2.0
    endfor
    
; now just consider the strong lines in the multiplet that were not
; ignored 
    broad = where((linepars.kind eq 'l') and (linepars.action ne 'i') and $
      (strmatch(linepars.name,'*broad*',/fold) eq 1),nbroad)
    strong = where((linepars.kind eq 'l') and (linepars.action ne 'i'),nstrong)
    weak = where((linepars.kind ne 'l') and (linepars.action ne 'i'),nweak)

    strong_linez = fltarr(nstrong)-1.0
    strong_sigma = fltarr(nstrong)-1.0
    if (nbroad ne 0) then begin
       broad_linez = fltarr(nbroad)-1.0
       broad_sigma = fltarr(nbroad)-1.0
    endif
    
; now pack it in!
    broad_count = 0
    for ii = 0, nstrong-1 do begin
; continuum flux
       tag_continuum = tag_indx(out,linepars[strong[ii]].name+'_continuum')
       out.(tag_continuum) = fluxscale*[linepars[strong[ii]].continuum,$
         linepars[strong[ii]].continuum_err]

; if the amplitude is >0 then store everything        
       tag_flux = tag_indx(out,linepars[strong[ii]].name)
       tag_amp = tag_indx(out,linepars[strong[ii]].name+'_amp')
       tag_linez = tag_indx(out,linepars[strong[ii]].name+'_linez')
       tag_sigma = tag_indx(out,linepars[strong[ii]].name+'_sigma')
       tag_ew = tag_indx(out,linepars[strong[ii]].name+'_ew')

       if (sol[ii*4+1] gt 0.0) then begin
          out.(tag_flux) = [sol[ii*4+0],esol[ii*4+0]]*fluxscale ; flux       
          out.(tag_amp) = [sol[ii*4+1],esol[ii*4+1]]*fluxscale  ; amplitude
          out.(tag_sigma) = [sol[ii*4+3],esol[ii*4+3]] ; sigma-width [km/s]

          vv = sol[ii*4+2]      ; km/s
          vv_err = esol[ii*4+2]
          linez = exp((vv+vsys)/light)-1.0
;         linez = exp((vv+vsys-voffset)/light)-1.0
          linez_err = linez*(vv_err/light)
          out.(tag_linez) = [linez,linez_err]

          if strmatch(linepars[strong[ii]].name,'*broad*',/fold) then begin
             broad_linez[broad_count] = linez
             broad_sigma[broad_count] = sol[ii*4+3]
             broad_count++
          endif else begin
             strong_linez[ii] = linez
             strong_sigma[ii] = sol[ii*4+3]
          endelse             

; compute the EW
          ff = sol[ii*4+0]
          ff_err = esol[ii*4+0]
          cc = linepars[strong[ii]].continuum
          if (linepars[strong[ii]].continuum_npix le 0) then message, 'Problem here'
          cc_err = linepars[strong[ii]].continuum_err/sqrt(linepars[strong[ii]].continuum_npix)
          if (cc ne 0.0) and (cc_err ne 0.0) then begin
             ew = ff/cc
             ew_err = sqrt(ff_err^2.0/cc^2.0 + cc_err^2.0*(ff/cc^2.0)^2.0)
             out.(tag_ew) = [ew,ew_err]
          endif
       endif
    endfor

; store the amplitude and flux of multiplets [OIII] 4959 and [NII]
; 6548
    for ii = 0, nweak-1 do begin
       jj = where(strmid(linepars[weak[ii]].kind,1) eq linepars.i)

       tag_weak_flux = tag_indx(out,(linepars[weak[ii]].name)[0])
       tag_strong_flux = tag_indx(out,(linepars[jj].name)[0])
       if ((out.(tag_strong_flux))[1] gt 0.0) then $
         out.(tag_weak_flux) = linepars[weak[ii]].a*out.(tag_strong_flux)

       tag_weak_continuum = tag_indx(out,(linepars[weak[ii]].name)[0]+'_continuum')
       tag_strong_continuum = tag_indx(out,(linepars[jj].name)[0]+'_continuum')
       if ((out.(tag_strong_continuum))[1] gt 0.0) then $
         out.(tag_weak_continuum) = out.(tag_strong_continuum) ; the same for both lines

       tag_weak_amp = tag_indx(out,(linepars[weak[ii]].name)[0]+'_amp')
       tag_strong_amp = tag_indx(out,(linepars[jj].name)[0]+'_amp')
       if ((out.(tag_strong_amp))[1] gt 0.0) then $
         out.(tag_weak_amp) = linepars[weak[ii]].a*out.(tag_strong_amp)

       tag_weak_limit = tag_indx(out,linepars[weak[ii]].name+'_limit')
       tag_strong_limit = tag_indx(out,linepars[jj].name+'_limit')
       if (out.(tag_strong_limit) gt 0.0) then $
         out.(tag_weak_limit) = linepars[weak[ii]].a*out.(tag_strong_limit)

       tag_weak_ew = tag_indx(out,(linepars[weak[ii]].name+'_ew')[0])
       tag_strong_ew = tag_indx(out,(linepars[jj].name+'_ew')[0])
       if ((out.(tag_strong_ew))[1] gt 0.0) then $
         out.(tag_weak_ew) = linepars[weak[ii]].a*out.(tag_strong_ew)

; if the sigma *measurement* (not error) is >0 then store it
; (remember: if SIGMA hits the limit of what we allow then the error
; is set to zero by MPFIT in IM_GANDALF)       
       tag_weak_sigma = tag_indx(out,(linepars[weak[ii]].name)[0]+'_sigma')
       tag_strong_sigma = tag_indx(out,(linepars[jj].name)[0]+'_sigma')
       if ((out.(tag_strong_sigma))[0] gt 0.0) then $
         out.(tag_weak_sigma) = out.(tag_strong_sigma)

       tag_weak_linez = tag_indx(out,(linepars[weak[ii]].name)[0]+'_linez')
       tag_strong_linez = tag_indx(out,(linepars[jj].name)[0]+'_linez')
; sometimes the redshift can be defined, but the redshift *error* is
; zero; not sure why
       if ((out.(tag_strong_linez))[0] gt 0.0) then $
;      if ((out.(tag_strong_linez))[1] gt 0.0) then $
         out.(tag_weak_linez) = out.(tag_strong_linez)
    endfor

;; combine the [OII] doublet
;    names = 'oii_3727'+['','_ew','_wave']
;    values = [replicate('fltarr(2)',n_elements(names)-1),'0.0']
;    oii = mrd_struct(names,values,1)
;
;    oii.oii_3727_wave = djs_mean([out.oii_3726_wave,out.oii_3729_wave])
;    oii.oii_3727[0] = out.oii_3726[0]+out.oii_3729[0]
;    oii.oii_3727[1] = sqrt(out.oii_3726[1]^2+out.oii_3729[1]^2)
;    oii.oii_3727_ew[0] = out.oii_3726_ew[0]+out.oii_3729_ew[0]
;    oii.oii_3727_ew[1] = sqrt(out.oii_3726_ew[1]^2+out.oii_3729_ew[1]^2)
;    
;    out = struct_addtags(out,oii)

; store the average forbidden-line and Balmer-line redshift and
; velocity width
    out = create_struct(out,'zline_balmer',-1.0,'zline_forbidden',-1.0,$
      'sigma_balmer',-1.0,'sigma_forbidden',-1.0,'zline_broad',-1.0,$
      'sigma_broad',-1.0)
    good = where(strong_linez ne -1.0,ngood)
    if (ngood ne 0) then begin
       isbalm = where(strmatch(linepars[strong[good]].name,'h*',/fold),comp=isforb) ; not general!!
;      isbalm = balmerindx(linepars[strong[good]].lambda,forb=isforb,/all)
       if (isbalm[0] ne -1) then begin
          out.zline_balmer = djs_mean(strong_linez[good[isbalm]])
          out.sigma_balmer = djs_mean(strong_sigma[good[isbalm]])
       endif
       if (isforb[0] ne -1) then begin
          out.zline_forbidden = djs_mean(strong_linez[good[isforb]])
          out.sigma_forbidden = djs_mean(strong_sigma[good[isforb]])
       endif
       sigma_limit = djs_mean(strong_sigma[good])
    endif else sigma_limit = -1.0

; deal with the broad lines
    if (nbroad ne 0) then begin
       good = where(broad_linez ne -1.0,ngood)
       if (ngood ne 0) then begin
          out.zline_broad = djs_mean(broad_linez[good])
          out.sigma_broad = djs_mean(broad_sigma[good])
       endif
    endif

; finally compute 1-sigma limiting fluxes and EWs
    if (sigma_limit le 0.0) then sigma_limit = 120.0 ; [km/s]
    for ii = 0, nline-1 do begin
       tag_wave = tag_indx(out,out.linename[ii]+'_wave')
       tag_continuum = tag_indx(out,out.linename[ii]+'_continuum')
       tag_flux = tag_indx(out,out.linename[ii])
       tag_ew = tag_indx(out,out.linename[ii]+'_ew')
       tag_limit = tag_indx(out,out.linename[ii]+'_limit')
       tag_ew_limit = tag_indx(out,out.linename[ii]+'_ew_limit')

       if ((out.(tag_continuum))[1] gt 0.0) then begin
          out.(tag_limit) = sqrt(2.0*!pi)*(out.(tag_continuum))[1]*$
            out.(tag_wave)*sigma_limit/im_light()
          out.(tag_ew_limit) = out.(tag_limit)/(out.(tag_continuum))[0]
;         splog, out.(tag_limit), (out.(tag_flux))[0], out.(tag_ew_limit), (out.(tag_ew))[0]
       endif
    endfor 
    
;; continuum reddening    
;    out = create_struct('continuum_ebv',0.0,$
;      'continuum_ebv_err',0.0,out)
;    out.continuum_ebv = sol[nstrong*4]
;    out.continuum_ebv_err = esol[nstrong*4]

return, out
end
