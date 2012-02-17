function isedfit_packit, isedfit, array, type=type, wquant=wquant
; support routine for isedfit_compute_posterior()
    
;   tagavg = tag_indx(isedfit,type+'_avg')
    tag50 = tag_indx(isedfit,type+'_50')
;   tagmode = tag_indx(isedfit,type+'_mode')

;   tagerr = tag_indx(isedfit,type+'_err')
    tagerr = tag_indx(isedfit,type+'_err')
;   tagefferr = tag_indx(isedfit,type+'_eff_err')

;   isedfit.(tagavg) = djs_mean(array)
;   isedfit.(tagerr) = djsig(array)

    quant = [1.0-gauss_pdf(2.0),0.5,gauss_pdf(2.0)]
    wquant = weighted_quantile(array,quant=quant)
    isedfit.(tag50) = wquant[1]
    isedfit.(tagerr) = (wquant[2]-wquant[0])/4.0
;   isedfit.(tagefferr) = (wquant[2]-wquant[0])/4.0

;; get the mode    
;    binsz = (wquant[2]-wquant[0])/n_elements(array)
;    if (binsz le 0.0) then binsz = (max(array)-min(array))/n_elements(array)
;    if (binsz le 0.0) then isedfit.(tagmode) = array[0] else begin
;       yhist = histogram(array,binsize=binsz,omin=omin,omax=omax)
;       xhist = dindgen(n_elements(yhist))*binsz+omin+binsz/2.0
;       mx = max(yhist,this)
;       isedfit.(tagmode) = xhist[this]
;    endelse
    
return, isedfit
end

