function gandalf_clean_emission, restwave, restflux, bestfit, $
  mask, etemplates, linepars, sol, esol, snrcut=snrcut, velscale=velscale, $
  new_etemplates=new_etemplates, new_linepars=new_linepars, debug=debug, $
  fitlines=fitlines
; Given the galaxy spectrum, the best fit, the emission-line
; amplitudes, a vector with the A/N threshold for each line, and the
; array containing the spectra of each best-fitting emission-line
; templates, this function simply compute the residual-noise lvl,
; compares it the amplitude of each lines, and remouves from the
; galaxy spectrum only the best-matching emission-line templates for
; which the correspoding A/N exceed the input threshold.  This is a
; necessary step prior to measurements of the strength of the stellar
; absorption line features
;
; A list of goodpixels may be optionally input, for instance if any
; pixel was excluded by sigma-clipping during the continuum and
; emission-line fitting, or by excluding pixels on either ends of the
; spectra
;
; Also optionally outputs the computed A/N ratios.

    light = 2.99792458D5        ; speed of light1 [km/s]
    if (n_elements(snrcut) eq 0) then snrcut = 3.0

;   isline = fitlines
;   nline = n_elements(fitlines)
    isline = where((linepars.kind eq 'l') and (linepars.action ne 'i'),nline)
    nallline = n_elements(linepars.name)
    moretags = replicate({continuum: 0.0, continuum_err: 0.0, $
      continuum_npix: 0},nallline)
    new_linepars = struct_addtags(linepars,moretags)
;   new_linepars = create_struct(linepars,'continuum',$
;     fltarr(nallline),'continuum_err',fltarr(nallline))

    new_sol = sol
    new_etemplates = etemplates

    resid = restflux - bestfit
    for ii = 0, nline-1 do begin
       linewave = new_linepars[isline[ii]].lambda

       flux = sol[ii*4+0]
       amp = sol[ii*4+1]
       sigma = sol[ii*4+3]

       flux_err = esol[ii*4+0]
       amp_err = esol[ii*4+1]
       sigma_err = esol[ii*4+3]

; occasionally the amplitude *error* is zero because the line hit a
; boundary in MPFIT; set the flux to zero, too, which isn't
; quite the right thing to do, but only happens rather rarely
       if (abs(amp) gt 0) and (amp_err eq 0.0) then begin
          splog, 'Error is zero on '+new_linepars[isline[ii]].name
;         amp = 0.0
       end
       
; compute the noise in the continuum centered on each emission line in
; order to assess whether the line has been detected at S/N>snrcut;
; require a minimum sigma-width of 50 km/s (about 1 pixel); also
; compute the upper limit on every line
;      indx = where((abs(exp(restwave)-linewave) lt 60) and $
;        (abs(exp(restwave)-linewave) gt 15) and (mask eq 1),nindx)
       indx = where((abs(exp(restwave)-linewave) lt 25.0*linewave*(sigma>(3*velscale))/light) and $
         (abs(exp(restwave)-linewave) gt 5.0*linewave*(sigma>(3*velscale))/light) and $
         (mask eq 1),nindx)
       if (nindx ne 0) then begin
          resid_noise = robust_sigma(resid[indx],/zero)
          if keyword_set(debug) then begin
             plot, exp(restwave), resid, xr=linewave+[-70,70], $
               title=new_linepars[isline[ii]].name, psym=10, $
               yr=[min(resid[indx]),max(resid[indx])>amp>(1.2*snrcut*resid_noise)], $
               xsty=1, ysty=3
             djs_oplot, exp(restwave[indx]), resid[indx], psym=4, color='red'
             djs_oplot, linewave*[1,1], !y.crange, line=0, color='cyan'
             djs_oplot, !x.crange, [0,0], line=0, color='orange'
             djs_oplot, !x.crange, +resid_noise*[1,1], line=1, color='blue'
             djs_oplot, !x.crange, -resid_noise*[1,1], line=1, color='blue'
             djs_oplot, !x.crange, snrcut*resid_noise*[1,1], line=5, color='orange'
             djs_oplot, !x.crange, amp*[1,1], line=1, color='green', thick=3
             splog, amp, resid_noise, amp/resid_noise, $
               '  '+linepars[isline[ii]].name
             cc = get_kbrd(1)
          endif
; 1-sigma upper limit; for simplicity assume an intrinsic velocity
; width of 125 km/s
;         new_linepars[isline[ii]].limit = sqrt(2.0*!pi)*$
;           new_linepars[isline[ii]].lambda*125.0/im_light()*resid_noise 

; now deal with the non-detections          
          if (amp lt snrcut*resid_noise) and $
            (strmatch(linepars[isline[ii]].name,'*broad*',/fold) eq 0) then begin
;         if (amp/resid_noise lt snrcut) then begin
; set the amplitude and flux to zero
;            if (abs(amp) gt 0) and (amp_err eq 0.0) then $
             if keyword_set(debug) then $
               splog, 'Removing '+new_linepars[isline[ii]].name, $
               amp, resid_noise, amp/resid_noise
             new_sol[ii*4+0] = 0.0 
             new_sol[ii*4+1] = 0.0
             new_etemplates[*,ii] = 0.0
          endif else if (abs(amp) gt 0) and (amp_err eq 0.0) then message, 'Fix me'
       endif
    endfor

; now that crummy lines have been removed, go back and compute the
; mean continuum level and the error in the mean, so that we can
; compute EWs later
    if (size(etemplates,/n_dim) eq 2) then $
      resid_nolines = restflux - total(etemplates,2) else $
        resid_nolines = restflux - etemplates
    for ii = 0, nline-1 do begin
       amp = new_sol[ii*4+1]
       sigma = new_sol[ii*4+3]
       linewave = new_linepars[isline[ii]].lambda
       
;      indx = where((abs(exp(restwave)-linewave) lt 60) and $
;        (abs(exp(restwave)-linewave) gt 15) and (mask eq 1),nindx)
       indx = where((abs(exp(restwave)-linewave) lt 25.0*linewave*(sigma>(3*velscale))/light) and $
         (abs(exp(restwave)-linewave) gt 5.0*linewave*(sigma>(3*velscale))/light) and $
         (mask eq 1),nindx)
       if (nindx gt 3) then begin
;      if (nindx ne 0) and (amp gt 0.0) then begin
          djs_iterstat, resid_nolines[indx], mean=mn, median=md, $
            sigma=sig, sigrej=3.0, mask=msk
          new_linepars[isline[ii]].continuum = mn
          new_linepars[isline[ii]].continuum_err = sig
          new_linepars[isline[ii]].continuum_npix = total(msk)
          
          if keyword_set(debug) then begin
             plot, exp(restwave), restflux, xr=linewave+[-70,70], $
               title=new_linepars[isline[ii]].name, psym=10, $
               yrange=[0,(mn+amp)*1.01], xsty=3, ysty=3
             djs_oplot, linewave*[1,1], !y.crange, line=0, color='cyan'
             djs_oplot, exp(restwave), resid_nolines, color='blue', psym=10
             djs_oplot, exp(restwave[indx]), resid_nolines[indx], psym=4, color='red'
             djs_oplot, !x.crange, mn*[1,1], line=0, color='orange'
             djs_oplot, !x.crange, (mn+3.0*sig/sqrt(total(msk)))*[1,1], line=5, color='orange'
             djs_oplot, !x.crange, (mn-3.0*sig/sqrt(total(msk)))*[1,1], line=5, color='orange'
             djs_oplot, !x.crange, (mn+amp)*[1,1], line=1, color='green', thick=3
             cc = get_kbrd(1)
          endif
       endif
    endfor 

return, new_sol
end
