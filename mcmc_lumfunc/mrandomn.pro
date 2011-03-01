function mrandomn, seed, covar, nrand, a=a

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
; Function to draw NRAND random deviates from a multivariate normal
; distribution with zero mean and covariance matrix COVAR.
; 
; AUTHOR : Brandon C. Kelly, Steward Obs., Sept. 2004
;
; INPUTS : 
;
;    SEED - The random number generator seed, the default is IDL's
;           default in RANDOMN()
;    COVAR - The covariance matrix of the multivariate normal
;            distribution.    
; OPTIONAL INPUTS :
;
;    NRAND - The number of randomn deviates to draw. The default is
;            one.
; OUTPUT :
;
;    The random deviates, an [NRAND, NP] array where NP is the
;    dimension of the covariance matrix, i.e., the number of
;    parameters.
;
; ROUTINES CALLED :
;
;    POSITIVE, DIAG
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

if n_params() lt 2 then begin
    print, 'Syntax- Result = mrandomn( seed, covar, [nrand], A=A )'
    return, 0
endif

;check inputs and set up defaults
if n_elements(nrand) eq 0 then nrand = 1
if size(covar, /n_dim) ne 2 then begin
    print, 'COVAR must be a matrix.'
    return, 0
endif

np = (size(covar))[1]
if (size(covar))[2] ne np then begin
    print, 'COVAR must be a square matrix.'
    return, 0
endif

diag = lindgen(np) * (np + 1L)
epsilon = randomn(seed, nrand, np) ;standard normal random deviates (NP x NRAND matrix)

if n_elements(A) eq 0 then begin

    A = covar ;store covariance into dummy variable for input into TRIRED
    
    choldc, A, P, /double       ;do Cholesky decomposition
    
    for j = 0, np - 1 do for i = j, np - 1 do A[i,j] = 0d
    
    A[diag] = P

endif

;transform standard normal deviates so they have covariance matrix COVAR
epsilon = A ## epsilon

return, epsilon
end
