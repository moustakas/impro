function hlmean,data,nsamples=nsamp
;+
; NAME:
;
;  HLMEAN
;
; PURPOSE:
;
;  Calculate Hodges-Lehmann estimator of mean, using NSAMPLES
;   bootstraps from the data
;
; The Hodges-Lehmann estimator is, formally, the median value of
;  (x_i+x_j)/2 over all pairs of indices i,j .
; Here, we estimate that quantity using nsamp randomly chosen values
; of i & j, rather than using all possible values.
; 
; Although it has much of the robustness of an ordinary median, the
; H-L estimator yields much smaller errors (equivalent to the mean of
; >90% as much data, while the median has errors equivalent to the
; standard error of 64% as much data).
;
; CATEGORY:
;
; Statistics
;
; CALLING SEQUENCE:
;
;  result=HLMEAN(data [,NSAMPLES=nsamples])
;
; INPUTS:
;
;  data: array of values to calculate the H-L mean of
;
; KEYWORD PARAMETERS:
;
;  NSAMPLES= : if set, HLMEAN will use this number of bootstrap
;  samples to do the calculation.  If not set, it will use 10*the
;  number of elements of the data array.
;
; OUTPUTS:
;
;  result: bootstrap estimate of the H-L mean estimator
;
; EXAMPLE:
;    test=[1,2,0,1,2,20.]
;    print,HLMEAN(test)
;
; MODIFICATION HISTORY:
;  JAN 22 sep. 08
;-




ndata=n_elements(data)

if n_elements(nsamp) eq 0 then nsamp=10.*ndata

idx=floor(randomu(seed,nsamp,2)*ndata)
mn=total(data[idx],2)/2.

return,median(mn,/even)
end
