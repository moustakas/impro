FUNCTION fitfunc_pxf_optimal_template, pars, $
    BESTFIT=bestFit, BIAS=bias, CLEAN=clean, DEGREE=degree, $
    FACTOR=factor, GALAXY=galaxy, GOODPIXELS=goodPixels, MDEGREE=mdegree, $
    NOISE=noise, QUIET=quiet, SKY=sky, STAR=star, VSYST=vsyst, WEIGHTS=weights, $
    REGUL=regul, REG_DIM=reg_dim, LAMBDA=lambda, losvd=losvd, $
  weff=weff, maggies=maggies, errmaggies=errmaggies, modelmaggies=modelmaggies
compile_opt idl2, hidden

; jm09nov18ucsd - added support to fit for the reddening when MOMENTS=0

s = size(galaxy)
nspec = s[0]
npix = s[1]
npars = n_elements(pars) - mdegree*nspec  ; Parameters of the LOSVD only
if n_elements(lambda) gt 1 then npars = npars - 1 ; Fitting reddening

; jm09nov13ucsd - photometry
nphoto = n_elements(weff) 
if (nphoto gt 0) then goodphoto = lindgen(nphoto)+npix ; convenient index array

; pars = [vel,sigma,h3,h4,...,m1,m2,...]    ; Velocities are in pixels
;
dx = ceil(abs(vsyst)+abs(pars[0])+5d*pars[1]) ; Sample the Gaussian and GH at least to vel+5*sigma
n = 2*dx*factor + 1
x = range(dx,-dx,n)   ; Evaluate the Gaussian using steps of 1/factor pixel
losvd = dblarr(n,nspec,/NOZERO)
for k=0,nspec-1 do begin    ; nspec=2 for two-sided fitting, otherwise nspec=1
    s = (k eq 0) ? 1d : -1d ; s=+1 for left spectrum, s=-1 for right one
    vel = vsyst + s*pars[0]
    w = (x - vel)/pars[1]
    w2 = w^2
    losvd[*,k] = exp(-0.5d*w2)/(sqrt(2d*!dpi)*pars[1]) ; Normalized total(Gaussian)=1

    ; Hermite polynomials normalized as in Appendix A of van der Marel & Franx (1993).
    ; These coefficients are given e.g. in Appendix C of Cappellari et al. (2002)
    ;
    if npars gt 2 then begin
        poly = 1d + s*pars[2]/Sqrt(3d)*(w*(2d*w2-3d)) $     ; H3
                  + pars[3]/Sqrt(24d)*(w2*(4d*w2-12d)+3d)   ; H4
        if npars eq 6 then $
            poly = poly + s*pars[4]/Sqrt(60d)*(w*(w2*(4d*w2-20d)+15d)) $    ; H5
                        + pars[5]/Sqrt(720d)*(w2*(w2*(8d*w2-60d)+90d)-15d)  ; H6
        losvd[*,k] = losvd[*,k]*poly
    endif
endfor

; The zeroth order multiplicative term is already included in the
; linear fit of the templates. The polynomial below has mean of 1.
;
x = range(-1d,1d,npix) ; X needs to be within [-1,1] for Legendre Polynomials
mpoly = 1d  ; The loop below can be null if mdegree < 1
for j=1,mdegree do $
    if nspec eq 2 then $ ; Different multiplicative poly for left and right spectra
        mpoly = mpoly + legendre(x,j) # pars[npars+2*j-[2,1]] $
    else mpoly = mpoly + legendre(x,j) * pars[npars+j-1]

; Multiplicative polynomials do not make sense when fitting reddening.
; In that case one has to assume the spectrum is well calibrated.
;
; jm09nov13ucsd - photometry
if n_elements(lambda) gt 1 then begin
   mpoly = 10D^(-0.4*pars[npars]*k_lambda(lambda,/calzetti,/silent))
   if (nphoto gt 0) then mpoly = [mpoly,10D^(-0.4*pars[npars]*k_lambda(weff,/calzetti,/silent))]
endif
;if n_elements(lambda) gt 1 then mpoly = reddening_curve(lambda, pars[npars]) 

; Fill the columns of the design matrix of the least-squares problem
;
s = size(star)
ss = size(sky)
if ss[0] lt 2 then nsky = ss[0] else nsky = ss[2] ; Number of sky spectra
if s[0] eq 2 then ntemp = s[2] else ntemp = 1     ; Number of template spectra
nrows = (degree + 1 + nsky)*nspec + ntemp
ncols = npix*nspec+nphoto ; jm09nov16ucsd 
if regul gt 0 then begin
    nr = n_elements(reg_dim)
    nreg = product(reg_dim-2)*nr + 2*total(reg_dim-2)*(nr-1) 
    ncols = ncols + nreg
endif     
c = dblarr(ncols,nrows)  ; This array is used for estimating predictions

for j=0,degree do $ ; Fill first columns of the Design Matrix
    if nspec eq 2 then begin
        leg = legendre(x,j)
        c[0,2*j] = [leg,leg*0d]   ; Additive polynomials for left spectrum
        c[0,2*j+1] = [leg*0d,leg] ; Additive polynomials for right spectrum
    endif else c[0:npix*nspec-1,j] = legendre(x,j) 

if factor gt 1 then pix = range(0d,s[1]-1d,s[1]*factor) ; Oversampled pixels range
tmp = dblarr(s[1],nspec,/NOZERO)
for j=0,ntemp-1 do begin
    if factor eq 1 then $ ; No oversampling of the template spectrum
        for k=0,nspec-1 do tmp[*,k] = convol(star[*,j],losvd[*,k],/EDGE_TRUNCATE) $
    else begin             ; Oversample the template spectrum before convolution
        st = interpolate(star[*,j],pix,CUBIC=-0.5)   ; Sinc-like interpolation
        for k=0,nspec-1 do tmp[*,k] = rebin(convol(st,losvd[*,k],/EDGE_TRUNCATE),s[1])
    endelse
; jm09nov16ucsd 
    if (n_elements(lambda) gt 1) then begin
       c[0:npix*nspec-1,(degree+1)*nspec+j] = mpoly[0:npix-1]*tmp[0:npix-1,*] ; reform into a vector
       if (nphoto gt 0) then c[npix*nspec:npix*nspec+nphoto-1,(degree+1)*nspec+j] = $
         mpoly[npix:npix+nphoto-1]*modelmaggies[*,j]
    endif else begin
       c[0:npix*nspec-1,(degree+1)*nspec+j] = mpoly*tmp[0:npix-1,*] ; reform into a vector
       if (nphoto gt 0) then c[npix*nspec:npix*nspec+nphoto-1,(degree+1)*nspec+j] = $
         mpoly*modelmaggies[*,j]
    endelse
;   c[0,(degree+1)*nspec+j] = (mpoly*tmp[0:npix-1,*])[*] ; reform into a vector
endfor

; Add second-degree 1D or 2D linear regularization
; Numerical Recipes 2nd ed. equation (18.5.11)
;
if regul gt 0 then begin
    dim = size(reg_dim,/DIM) 
    if dim le 1 then $            ; 1D regularization
        for j=0,reg_dim-3 do $
            c[npix*nspec+j,(degree+1)*nspec+1+j+[-1,0,1]] = [-1d,2d,-1d]*regul $
    else if dim eq 2 then begin  ; 2D regularization
        p = npix*nspec
        for j=0,reg_dim[1]-1 do $  
            for k=0,reg_dim[0]-1 do begin ; First index moves faster in templates
                i = (degree+1)*nspec + j*reg_dim[0] + k
                if k ne 0 && k ne reg_dim[0]-1 then c[p++,i+[-1,0,1]] = [-1d,2d,-1d]*regul
                if j ne 0 && j ne reg_dim[1]-1 then c[p++,i+[-1,0,1]*reg_dim[0]] = [-1d,2d,-1d]*regul
            endfor
    endif else message, 'Unsupported regularization dimension'
endif     

for j=0,nsky-1 do begin
    skyj = sky[*,j]
    k = (degree+1)*nspec + ntemp
    if nspec eq 2 then begin
        c[0,k+2*j] = [skyj,skyj*0d]   ; Sky for left spectrum
        c[0,k+2*j+1] = [skyj*0d,skyj] ; Sky for right spectrum
    endif else c[0,k+j] = skyj
endfor

a = c                     ; This array is used for the actual solution of the system
; jm09nov16ucsd 
for j=0,nrows-1 do begin
   a[0:npix*nspec-1,j] = c[0:npix*nspec-1,j]/noise ; Weight all columns with errors
   if (nphoto gt 0) then $
     a[npix*nspec:npix*nspec+nphoto-1,j] = c[npix*nspec:npix*nspec+nphoto-1,j]/errmaggies ; Weight all columns with errors
endfor
;for j=0,nrows-1 do a[0,j] = c[0:npix*nspec-1,j]/noise ; Weight all columns with errors

if regul gt 0 then begin
    aa = a[[goodPixels,range(npix*nspec,ncols-1)],*]
    bb = [galaxy[goodPixels]/noise[goodPixels],replicate(0d,nreg)] 
endif else begin
; jm09nov16ucsd
   if (nphoto gt 0) then begin
      aa = a[[goodPixels,goodphoto],*]
      bb = [galaxy[goodPixels],maggies]/[noise[goodPixels],errmaggies]
   endif else begin
      aa = a[goodPixels,*]
      bb = galaxy[goodPixels]/noise[goodPixels]
   endelse
endelse

; Select the spectral region to fit and solve the overconditioned system
; using SVD/BVLS. Use unweighted array for estimating bestfit predictions.
; Iterate to exclude pixels deviating more than 3*sigma if /CLEAN keyword is set.
;
npoly = (degree+1)*nspec ; Number of additive polynomials in the fit 

repeat begin
   weights = BVLS_Solve_pxf(aa,bb,npoly)
    bestfit = c[0:npix*nspec+nphoto-1,*] # weights
; jm09nov16ucsd
    if (nphoto gt 0) then begin
       err = ([galaxy[goodPixels],maggies]-bestfit[[goodPixels,goodphoto]])/$
         [noise[goodPixels],errmaggies]
    endif else begin
       err = (galaxy[goodPixels]-bestfit[goodPixels])/noise[goodPixels]
    endelse
    if keyword_set(clean) then begin
        tmp = where(abs(err) gt 3, m, COMPLEM=w,ncomp=nw) ; select errors larger than 3*sigma
        if (m ne 0) then begin
            if ~keyword_set(quiet) then print, 'Outliers:', m
            if (nw eq 0) then message, 'Problem here!'
            goodPixels = goodPixels[w]
        endif
    endif else break
endrep until (m eq 0)

; Penalize the solution towards (h3,h4,...)=0 if the inclusion of
; these additional terms does not significantly decrease the error.
;
if npars gt 2 && bias ne 0 then $
    err = err + bias*robust_sigma(err, /ZERO)*sqrt(total(pars[2:*]^2))

return, err
END
;----------------------------------------------------------------------------
