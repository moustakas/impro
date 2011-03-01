FUNCTION ROBUST, mags, sigs, numit
;+
; NAME:
;   ROBUST
; PURPOSE:
;   Makes a robust estimate of the mean of multiple observations of a
;   star's magnitude in an iterative way.
; CALLING SEQUENCE:
;   newmag = ROBUST(mags, sigs, numit) 
; INPUTS:
;   mags - an array of magnitudes measured for a single star
;   sigs - an array of the corresponding errors in the measurements
;   numit - the number of times to iterate 
; OUTPUTS:
;   newmag = the estimate of the mean magnitude
; PROCEDURE:
;   Called by calibrate
; MODIFICATION HISTORY:
;   Modified from Jeff Newman's MEANMAGF2 on 4/10/2000 by Bryan
;   Mendez, Astronomy Department, University of California at Berkeley
;-

nummin=0
mean=0.
sigtot=0.
oldmean=0.

; Convert to Data Numbers
magtmp = 10.^((30 - mags)/2.5)
sigtmp = (ALOG(10)/2.5)*magtmp*sigs

; Flag bad stars
tmp=WHERE(magtmp le 1., count)
IF count ne 0 THEN BEGIN
	magtmp[tmp]=1.
	sigtmp[tmp]=1.E30
ENDIF

; Sort the errors
sigsrt=sigtmp[SORT(magtmp)]

;Determine preliminary mean magnitude using the median
mean = MEDIAN(magtmp[WHERE(magtmp gt 1.)],/EVEN)

sigold = sigtmp
sigtmp = 1./(sigtmp)^2

; iteratively refine the mean determination
FOR j = 0,numit-1 DO BEGIN
	newmean = 0.
	oldmean = mean
	sigtot = 0.

	sigtmp = sigtmp/(1.+((magtmp-mean)/(2.*sigold))^2)
	sigtot = TOTAL(sigtmp)
	mean = TOTAL(magtmp*sigtmp)/sigtot
ENDFOR

newmag =  30 - 2.5*ALOG10(mean)

RETURN, newmag
END



