;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
; FIT A SCHECHTER LUMINOSITY FUNCTION TO SET OF OBSERVED
; LUMINOSITIES. USES THE METROPOLIS-HASTINGS ALGORITHM TO SIMULATE
; RANDOM DRAWS FROM THE POSTERIOR DISTRIBUTION.
;
; AUTHOR : BRANDON C. KELLY, STEWARD OBSERVATORY, DEC 2007
;
; INPUTS :
;
;    LUM - THE SET OF LUMINOSITIES
;    LLIM - DETECTION LIMIT, IN THE SAME UNITS AS LUM.
;    SKYFRAC - THE FRACTION OF THE SKY COVERED BY THE SURVEY OF
;              OBJECTS CORRESPONDING TO THE SET OF LUMINOSITIES IN
;              LUM.
;
; OPTIONAL INPUTS : 
;
;    NITER - THE TOTAL NUMBER OF ITERATIONS TO PERFORM IN THE
;            METROPOLIS-HASTINGS ALGORITHM. THE DEFAULT IS 20000.
;
;    SILENT - BE QUIET
;
; OUTPUTS :
;
;    POST - A STRUCTURE CONTAINING RANDOM DRAWS FROM THE POSTERIOR
;           DISTRIBUTION OF THE PARAMETERS. THE TAGS OF POST ARE:
;
;       LSTAR : THE APPROXIMATE LOCATION OF THE TURNOVER BETWEEN THE
;               FAINT POWER-LAW SECTION OF THE SCHECHTER FUNCTION, AND
;               THE BRIGHT EXPONENTIAL CUT-OFF SECTION. LSTAR IS IN
;               UNITS OF LOG-LUMINOSITY.
;
;           K : THE SLOPE OF THE FAINT POWER-LAW END OF THE SCHECHTER
;               LUMINOSITY FUNCTION.
;
;           N : THE TOTAL NUMBER OF OBJECTS, INCLUDING THOSE BELOW THE
;               DETECTION LIMIT. THIS IS THE NORMALIZATION FOR THE
;               LUMINOSITY FUNCTION.
;
; REFERENCE:
;
;     "A FLEXIBLE METHOD OF ESTIMATING LUMINOSITY FUNCTIONS", KELLY,
;        B.C., FAN, X., & VESTERGAARD, M. 2008, ACCEPTED BY ApJ,
;        (arXiv:0805.2946)
;
; ROUTINES CALLED :
;
;  SCHECHTER_MLE, MRANDOMN, NEGBIN, CRAMER_INVERT
;;

; REVISION HISTORY:
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function schechter_metro_update, logpost_new, logpost_old, seed, log_jrat

lograt = logpost_new - logpost_old

if n_elements(log_jrat) gt 0 then lograt = lograt + log_jrat

accept = 0

if lograt gt 0 then accept = 1 else begin
    
    u = randomu(seed)

    if alog(u) le lograt then accept = 1
    
endelse

return, accept
end

pro schechter_mha, lum, llim, skyfrac, post, niter=niter, silent=silent

if n_params() lt 4 then begin
    print, 'Syntax- schechter_mha, lum, llim, skyfrac, post, niter=niter, /silent'
    return
endif

nobs = n_elements(lum)
nchains = 5
nburn = 2000L                   ;do 2000 iterations of burn-in
if n_elements(niter) eq 0 then niter = 20000L

if n_elements(silent) eq 0 then silent = 0

kdraw = dblarr(niter, nchains)
lsdraw = dblarr(niter, nchains)

burnin = 1
checkiter = 500L

if niter le nburn then niter = nburn + 1
naccept = 0L

;first get initial guesses for MHA algorithm

schechter_mle, lum, alog10(llim), skyfrac, theta0, covar=jumpvar, /silent

theta = alog(10) * theta0[1:2] # replicate(1, nchains)
jumpvar = jumpvar[1:2,1:2]

if determ(jumpvar) le 1d-5 then begin

    jumpvar0 = jumpvar
    jumpvar = dblarr(2,2)
    jumpvar[0,0] = abs(jumpvar0[0,0])
    jumpvar[1,1] = abs(jumpvar0[1,1])

endif

for i = 0, nchains - 1 do theta[*,i] = theta[*,i] + mrandomn(seed, jumpvar)

iter = 0L

logpost_old = dblarr(nchains)

for i = 0, nchains - 1 do begin
    
    theta0 = theta[*,i]

    alpha = exp(theta0[0]) - 1
    lstar = theta0[1]
    
    det_loglik = alpha * (alog(lum) - lstar) - lstar - $
      lngamma(alpha + 1) - lum / exp(lstar)
        
    det_prob = $
      (1d - igamma(alpha + 1, llim / exp(lstar), /double, $
                   itmax=1d7, eps=1d-3)) * skyfrac
    
    logpost_old[i] = total(det_loglik) - nobs * alog(det_prob)
       
endfor

repeat begin

;    if iter lt 100 then print, iter

    for i = 0, nchains - 1 do begin

        tprop = theta[*,i] + mrandomn(seed, jumpvar, A=A_jump)

        aprop = exp(tprop[0]) - 1
        lsprop = tprop[1]

        det_loglik = aprop * (alog(lum) - lsprop) - lsprop - $
          lngamma(aprop + 1) - lum / exp(lsprop)
        
        det_prob = $
          (1d - igamma(aprop + 1, llim / exp(lsprop), /double, $
                       itmax=1d7, eps=1d-3)) * skyfrac
        
        logpost_new = total(det_loglik) - nobs * alog(det_prob)
        
        log_jrat = tprop[0] - theta[0,i]

        accept = $
          schechter_metro_update( logpost_new, logpost_old[i], seed, log_jrat )

        if accept then begin

            theta[*,i] = tprop
            logpost_old[i] = logpost_new
            naccept = naccept + 1L

        endif

    endfor

    iter = iter + 1

    if burnin then begin
    
        if n_elements(thetadraw0) eq 0 then thetadraw0 = theta $
        else thetadraw0 = [[thetadraw0], [theta]]
        
    endif else begin

        kdraw[iter - 1,*] = exp(theta[0,*])
        lsdraw[iter - 1,*] = theta[1,*]

    endelse

    if iter eq checkiter then begin

        if burnin then begin
                                ;update jumping covariance matrix
            thetadraw = reform(thetadraw0, 2, nchains, iter)

            jumpvar1 = 0d
            for i = 0, nchains - 1 do jumpvar1 = jumpvar1 + $
              cramer_invert( correlate(reform(thetadraw[*,i,*]), /cov) )
            jumpvar2 = nchains * cramer_invert(jumpvar1)

            if total(finite(jumpvar2)) eq 4 and $
              (jumpvar2[0,0] * jumpvar2[1,1] - jumpvar2[0,1]^2) gt 1d-5 then $
              jumpvar = jumpvar2 $
            else jumpvar = jumpvar / 36d

            delvarx, A_jump

        endif

        if not silent then begin

            print, ''
            print, 'Iteration:'
            print, iter
            print, ''
            print, 'Acceptance Rate since last report:'
            print, naccept / 500d / 5d

        endif

        naccept = 0L
        checkiter = checkiter + 500L

    endif

    if iter eq nburn and burnin then begin

        if not silent then begin

            print, ''
            print, '--------------------- Burn-in complete ----------------------'
            print, ''

        endif

        checkiter = 500L
        iter = 0L
        burnin = 0
        delvarx, thetadraw0

    endif

endrep until iter eq niter

post = {lstar:0d, k:0d, n:0d}

post = replicate(post, n_elements(kdraw))

post.lstar = lsdraw[*] / alog(10)
post.k = kdraw[*]

ndraw = n_elements(kdraw)
if not silent then print, 'Getting Posterior Draws for N...'
for i = 0L, ndraw - 1 do begin

    if not silent and i mod 10000 eq 0 then print, i

    detprob = $
      (1d - igamma(kdraw[i], llim / exp(lsdraw[i]), /double, $
                   itmax=1d7, eps=1d-3)) * skyfrac > 1d-8
    post[i].n = nobs + negbin(nobs, detprob, seed=seed)

endfor

return
end
