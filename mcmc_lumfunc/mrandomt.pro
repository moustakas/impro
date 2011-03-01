;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
; ROUTINE TO DRAW RANDOM NUMBERS FROM THE MULTIVARIATE STUDENT'S
; T-DISTRIBUTION, WITH MEAN = 0 AND COVAR = COVAR.
;
; AUTHOR : BRANDON C. KELLY, STEWARD OBS., APRIL 2006
;
; INPUT :
;
;   NU - THE DEGREES OF FREEDOM FOR THE T-DISTRIBUTION
;   COVAR - THE COVARIANCE MATRIX OF THE T-DISTRIBUTION
;
; OPTIONAL INPUTS :
;
;   N - THE NUMBER OF RANDOM VARIABLES TO DRAW.
;   SEED - SEED FOR THE RANDOM NUMBER GENERATOR
;
; OUTPUT :
;
;   RANDOM NUMBER(S) DRAWN FROM A STUDENT'S MULTIVARIATE
;   T-DISTRIBUTION, WITH MEAN = 0, COVARIANCE MATRIX = COVAR, AND
;   DEGREES OF FREEDOM NU. TO GET A RANDOM VARIABLE DRAWN FROM A
;   T-DISTRIBUTION WITH MEAN = MU SIMPLY ADD MU TO THE OUTPUT.
;
; CALLED ROUTINES:
;
;   RANDOMCHI, MRANDOMN, CMREPLICATE
; 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function mrandomt, seed, nu, covar, nrand, A=A

if n_params() lt 3 then begin
    print, 'Syntax- t = mrandomt( seed, nu, covar[, nrand], A=A )'
    return, 0
endif

if nu lt 1 then begin
    print, 'NU must be at least 1.'
    return, 0
endif

if n_elements(nrand) eq 0 then n = 1
ndim = (size(covar, /dim))[0]

z = mrandomn(seed, covar, nrand, A=A)

x = randomchi(seed, nu, nrand)

t = z * sqrt(nu / cmreplicate(x, ndim))

return, t
end
