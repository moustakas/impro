;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
; ROUTINE TO DRAW RANDOM NUMBERS FROM THE NEGATIVE BINOMIAL
; DISTRIBUTION. A RANDOM VARIABLE DISTRIBUTED ACCORDING TO THE
; NEGATIVE BINOMIAL DISTRIBUTION GIVES THE NUMBER OF TRIALS NEEDED TO
; OBTAIN N SUCCESSES WHEN THE PROBABILITY OF SUCCESS IS P.
;
; AUTHOR : BRANDON C. KELLY, STEWARD OBS., APRIL 2006
;
; INPUTS :
;
;    N - THE NUMBER OF 'SUCCESSES'.
;    P - THE PROBABILITY OF 'SUCCESS'.
;
; OPTIONAL INPUTS :
;
;    NRAND - THE NUMBER OF RANDOM NUMBER TO DRAW.
;    SEED - THE SEED FOR THE RANDOM NUMBER GENERATOR.
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function negbin, n, p, nrand, seed=seed

if n_params() lt 2 then begin
    print, 'Syntax- M = negbin( n, p, [nrand, seed=seed] )'
    return, 0
endif

if n lt 1 then begin
    print, 'N must be at least 1'
    return, 0
endif

if p le 0 or p ge 1 then begin
    print, 'P must be between 0 and 1.'
    return, 0
endif

q = 1 - p
n = long(n)

if n_elements(nrand) eq 0 then nrand = 1

u = randomu(seed, n, nrand, /double)

m = total( floor(alog(u) / alog(q)), 1 )

return, m
end
