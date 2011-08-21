;+
; NAME:
;  SIGFIG
; PURPOSE:
;	SIGFIG is a function which returns an input value (INFIG) as
;	a string rounded to the appropriate number of significant
;	figures (NFIGS), as taught by Bonnie Kelly, my high school 
; 	physics teacher.  
;
;
; CALLING SEQUENCE:
;   result = SIGFIG(input, n_figures)
;
; INPUTS:
;   INPUT -- Vector of floats or whatever (they become floats)
;   N_FIGURES -- Number of figures to round to.
; KEYWORD PARAMETERS:
;   None.
;
; OUTPUTS:
;   RESULT -- A well rounded string array.
;
; MODIFICATION HISTORY:
;      Documented.
;       Wed Nov 21 12:25:29 2001, Erik Rosolowsky <eros@cosmic>
;-

function sigfig, inarr, nfigs

  outarr = strarr(n_elements(inarr))

  for ii = 0, n_elements(inarr)-1 do begin 
    infig = inarr[ii]
    num = strtrim(string(infig), 2)

    n = intarr(10)
    n[0] = strpos(num, '0')
    n[1] = strpos(num, '1')
    n[2] = strpos(num, '2')
    n[3] = strpos(num, '3')
    n[4] = strpos(num, '4')
    n[5] = strpos(num, '5')
    n[6] = strpos(num, '6')
    n[7] = strpos(num, '7')
    n[8] = strpos(num, '8')
    n[9] = strpos(num, '9')
    dec = strpos(num, '.')
    epos = strpos(num, 'e')
    neg = strpos(num, '-')
    errstr = 'Error: More Significant Figures than Input Digits'
    length = strlen(num)-1.*(neg eq 0)

    if length le nfigs then begin
      if dec ne -1 then begin
        print, errstr
        return, num
      endif
      if length lt nfigs then begin
        print, errstr
        outarr[ii] = strtrim(string(num, '.'), 2) 
        goto, doneit
      endif
    endif

    if epos ne -1 then begin
      rnd1 = strmid(num, 0, nfigs+1)
      rnd2 = strmid(num, epos, 6)
      rnd = string(rnd1, rnd2)
      outarr[ii] = rnd
      goto, doneit
    endif
    zerostring = '00000000000'

    if total(n(1:9)) eq -9 then begin
      if dec eq -1 then return, '0'
      rnd1 = '0.'
      rnd2 = strmid(zerostring, 0, nfigs-1)
      outarr[ii] = string(rnd1, rnd2)
      goto, doneit
    endif

    if n(0) eq 0 or ((neg eq 0) and (n(0) eq 1)) then begin
      n(where(n eq -1)) = 1000 
      firstval = min(n(1:9), subs)
      value = float(round(infig*10.^(firstval+nfigs-2)))*$
        (0.1^(firstval+nfigs-2))
      num = strtrim(string(value), 2)
      rnd = strmid(num, 0, firstval+nfigs)
    endif else begin
      if dec eq -1 then begin
        dec = strlen(num)-1.*(neg eq 0) 
        nfig1 = dec
      endif
      if dec le nfigs then nfig1 = nfigs+1 else nfig1 = dec
      
      n(where(n eq -1)) = 1000 
      firstval = min(n(1:9), subs)
      if neg eq 0 then begin
        nfig1 = nfig1+1
;		firstval=firstval-1
      endif
      exp = dec-(firstval+nfigs)
      value = float(round(infig*0.1^(exp)))*(10.^(exp))
      num = strtrim(string(value), 2)
      outarr[ii] = strmid(num, 0, nfig1)
      goto, doneit
    endelse

    outarr[ii] = rnd
doneit:
  endfor

  return, outarr
end


