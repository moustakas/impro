pro simulate_r23oh, snrcut=snrcut, debug=debug, postscript=postscript, $
  noerror=noerror, test=test
; jm06apr18uofa - 

    light = 2.99792458D10 ; speed of light [cm/s]

; start with an LZ relation; assign LOGQ values uniformly between the
; allowable boundaries; convert the O/H values into R23 and O32

    dmb = 0.015
    mb_nonoise = findgen(((-14.0)-(-22.0))/dmb+1)*dmb+(-22.0)
    nobj = n_elements(mb_nonoise)
    
    logoh = 5.24 - 0.185*mb_nonoise
    logq_all = alog10([0.5,1.0,2.0,4.0,8.0,15.0,30.0]*1D7)
    logq = logq_all[fix(randomu(seed,n_elements(logoh))*n_elements(logq_all))]
    
    r23_nonoise = oh2r23(logoh,logq=logq,o32=o32_nonoise,debug=debug)

; add noise, compute abundances, assign branches, and plot

    if (n_elements(snr) eq 0L) then snr = [50.0,25.0,15.0,10.0,5.0,3.0]
    nsnr = n_elements(snr)

    psname = '/home/ioannis/idl/impro/nebular/simulate_r23oh.ps'
    if keyword_set(postscript) then dfpsplot, psname, /square, /color

    for i = 0L, nsnr-1L do begin

       r23_err = r23_nonoise/snr[i]
       o32_err = o32_nonoise/snr[i]
       
       r23 = r23_nonoise + randomn(seed,nobj)*r23_err
       o32 = o32_nonoise + randomn(seed,nobj)*o32_err

       mb = mb_nonoise + randomn(seed,nobj)*0.1
       
       abund = kk04oh(r23,r23_err,o32,o32_err)
       result = im_assign_r23branch(abund,mb=mb,branchmethod=1L,/kk04,test=test)

       good = where(result.zstrong_12oh_kk04 gt -900.0)
       up = where(result.r23branch_kk04 eq 'U')
       lo = where(result.r23branch_kk04 eq 'L')
       amb = where(result.r23branch_kk04 eq 'A')

;      yrange = [min(result[lo].zstrong_12oh_kk04-result[lo].zstrong_12oh_kk04_err),$
;        max(result[up].zstrong_12oh_kk04+result[up].zstrong_12oh_kk04_err)]
       yrange = minmax(logoh)+[-0.2,+0.2]

       djs_plot, mb_nonoise, logoh, ysty=3, xsty=3, line=0, thick=3.0, charthick=3.0, $
         xthick=3.0, ythick=3.0, charsize=1.5, xtitle='M_{B}', ytitle='12 + log (O/H)', $
         yrange=yrange
       sixlin, mb[good], result[good].zstrong_12oh_kk04, lzint, siga, lzslope, sigb
       oplot, mb_nonoise, poly(mb_nonoise,[lzint[2],lzslope[2]]), line=2, thick=5
       print, lzint[2], lzslope[2]
       if keyword_set(noerror) then begin
          djs_oplot, mb[up], result[up].zstrong_12oh_kk04, ps=4, color='red'
          djs_oplot, mb[lo], result[lo].zstrong_12oh_kk04, ps=4, color='blue'
          djs_oplot, mb[amb], result[amb].zstrong_12oh_kk04, ps=4, color='dark green'
       endif else begin
          oploterror, mb[up], result[up].zstrong_12oh_kk04, mb[up]*0.0, result[up].zstrong_12oh_kk04_err, $
            ps=4, color=djs_icolor('red'), errcolor=djs_icolor('red'), /nohat, errstyle=2
          oploterror, mb[lo], result[lo].zstrong_12oh_kk04, mb[lo]*0.0, result[lo].zstrong_12oh_kk04_err, $
            ps=4, color=djs_icolor('blue'), errcolor=djs_icolor('blue'), /nohat, errstyle=0
          oploterror, mb[amb], result[amb].zstrong_12oh_kk04, mb[amb]*0.0, result[amb].zstrong_12oh_kk04_err, $
            ps=4, color=djs_icolor('dark green'), errcolor=djs_icolor('dark green'), /nohat, errstyle=1
       endelse

       legend, 'S/N = '+string(snr[i],format='(I0)'), /right, /top, box=0, $
         charsize=1.5, charthick=3.0

       if (not keyword_set(postscript)) then cc = get_kbrd(1)
       
    endfor
    if keyword_set(postscript) then dfpsclose
       
stop    
    
return
end
    
