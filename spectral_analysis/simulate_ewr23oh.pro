function return_oiihb, hii, oiiihb
; jm06apr21uofa - sample the distribution of [OII]/Hb vs [OIII]/Hb
;                 space spanned by HII regions
    
;   ages = read_ages_mz_sample()
;   plothist, ages.oii_3727[0]/ages.h_beta[0], bin=0.1, xbin, ybin, /peak
;   res = arm_asymgaussfit(xbin,ybin,/plot,lpars=lpars,rpars=rpars)
;   plothist, alog10(ages.oii_3727[0]/ages.h_beta[0]), bin=0.05, xbin, ybin, /peak
;   yfit = mpfitpeak(xbin,ybin,a,nterms=4,/positive) & print, a & djs_oplot, xbin, yfit, color='red'
;   plothist, ages.zstrong_o32, bin=0.05, xbin, ybin, /peak
;   yfit = mpfitpeak(xbin,ybin,a,nterms=4,/positive) & print, a & djs_oplot, xbin, yfit, color='red'

    indx = where(hii.oiii_h_beta gt -900.0 and hii.oii_h_beta gt -900.0)
    plot, hii[indx].oiii_h_beta, hii[indx].oii_h_beta, ps=3
    running = im_medxbin(hii[indx].oiii_h_beta,hii[indx].oii_h_beta,0.1,$
      minx=-1.2,maxx=1.0,minpts=minpts,/verbose)
    oploterror, running.binctr, running.medy, replicate(running.binsz/2.0,n_elements(running.medy)), $
;     running.stddev, ps=-4, xsty=3, ysty=3
      running.sigy, ps=-4, xsty=3, ysty=3

    bin = 0.005
    oiiihb = findgen((running.maxx-running.minx)/bin+1)*bin+running.minx
    oiihb = interpol(running.medy,running.binctr,oiiihb)
    oiihb_err = interpol(running.sigy,running.binctr,oiiihb)

    oiihb = oiihb + randomn(seed,n_elements(oiiihb))*oiihb_err
    djs_oplot, oiiihb, oiihb, ps=4, color='red'
    
return, oiihb    
end

pro simulate_ewr23oh, hii, snr=snr, debug=debug, postscript=postscript, $
  noerror=noerror, test=test
; jm06apr21uofa - 

    if (n_elements(hii) eq 0L) then hii = read_hii_regions()
    if (n_elements(snr) eq 0L) then snr = 10.0

    alpha = [0.8,0.9,1.0,1.1]
    nalpha = n_elements(alpha)

; compute noiseless R23 and O32 values
    
    oiihb_nonoise = return_oiihb(hii,oiiihb_nonoise)
    nobj = n_elements(oiihb_nonoise)

    o32_nonoise = oiiihb_nonoise - oiihb_nonoise
    r23_nonoise = alog10(10.0^oiihb_nonoise + 10.0^oiiihb_nonoise)

    plot, r23_nonoise, o32_nonoise, ps=4, xsty=3, ysty=3

    psname = '/home/ioannis/idl/impro/nebular/simulate_ewr23oh.ps'
    if keyword_set(postscript) then dfpsplot, psname, /square, /color

    for i = 0L, nalpha-1L do begin

; compute EWR23 and EWO32 with no noise, for this value of alpha       
       
       ewo32_nonoise = oiiihb_nonoise - alpha[i]*oiihb_nonoise
       ewr23_nonoise = alog10(alpha[i]*10.0^oiihb_nonoise + 10.0^oiiihb_nonoise)

       plot, r23_nonoise, ewr23_nonoise, ps=4, xsty=3, ysty=3
       djs_oplot, !x.crange, !y.crange & cc = get_kbrd(1)

       plot, o32_nonoise, ewo32_nonoise, ps=4, xsty=3, ysty=3
       djs_oplot, !x.crange, !y.crange & cc = get_kbrd(1)

; now compute R23, EWR23, O32, and EWO32, with noise

       oiihb_err = oiihb_nonoise/snr
       oiihb = oiihb_nonoise + randomn(seed,nobj)*oiihb_err
       oiiihb_err = oiiihb_nonoise/snr
       oiiihb = oiiihb_nonoise + randomn(seed,nobj)*oiiihb_err
       
       o32 = oiiihb - oiihb
       r23 = alog10(10.0^oiihb + 10.0^oiiihb)
       ewo32 = oiiihb - alpha[i]*oiihb
       ewr23 = alog10(alpha[i]*10.0^oiihb + 10.0^oiiihb)

       plot, r23_nonoise, r23, ps=4, xsty=3, ysty=3
       djs_oplot, !x.crange, !y.crange & cc = get_kbrd(1)
       plot, r23, ewr23, ps=4, xsty=3, ysty=3
       djs_oplot, !x.crange, !y.crange & cc = get_kbrd(1)
       plot, o32, ewo32, ps=4, xsty=3, ysty=3
       djs_oplot, !x.crange, !y.crange & cc = get_kbrd(1)

stop       
       
       if (not keyword_set(postscript)) then cc = get_kbrd(1)
       
    endfor
    if keyword_set(postscript) then dfpsclose
       
stop    
    
return
end
    
