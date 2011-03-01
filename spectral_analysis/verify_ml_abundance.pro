pro verify_ml_abundance, debug=debug
; jm04jun27uofa
; verify the abundances returned by ML_ABUNDANCE().  do this by
; drawing "HII regions" from the models themselves, adding Gaussian
; errors, and attempting to recover the true abundance

    grids = read_kewley_grids(Z=Z,U=U,/nolog)
    grids = grids[*,3]
    logU = U[3]
    nZ = n_elements(Z)

    nmonte = 100L

    rindx = fix(randomu(seed,nmonte)*nZ)
    Zinput = alog10(Z[rindx]) + 8.69
    montegrids = grids[rindx]

    ratios = ['nii_6584_oii','nii_6584_h_alpha','oiii_5007_oii']
    
    data = {$
      Z_Te:                 0.0D, $
;     R23:                  0.0D, $
;     R23_err:              0.0D, $
;     O32:                  0.0D, $
;     O32_err:              0.0D, $
      nii_6584_oii:         0.0D, $
      nii_6584_h_alpha:     0.0D, $
      oiii_5007_oii:        0.0D, $
      nii_6584_oii_err:     0.0D, $
      nii_6584_h_alpha_err: 0.0D, $
      oiii_5007_oii_err:    0.0D}
    data = replicate(data,nmonte)

    data.Z_Te = Zinput
    
    dfrac = 0.01
    
;   data.R23_err = montegrids.R23*dfrac
;   data.R23 = montegrids.R23 + randomn(seed,nmonte)*data.R23_err
;   data.O32_err = montegrids.O32*dfrac
;   data.O32 = montegrids.O32 + randomn(seed,nmonte)*data.O32_err
    data.nii_6584_oii_err = montegrids.nii_6584_oii*dfrac
    data.nii_6584_oii = montegrids.nii_6584_oii + randomn(seed,nmonte)*data.nii_6584_oii_err
    data.nii_6584_h_alpha_err = montegrids.nii_6584_h_alpha*dfrac
    data.nii_6584_h_alpha = montegrids.nii_6584_h_alpha + randomn(seed,nmonte)*data.nii_6584_h_alpha_err
    data.oiii_5007_oii_err = montegrids.oiii_5007_oii*dfrac
    data.oiii_5007_oii = montegrids.oiii_5007_oii + randomn(seed,nmonte)*data.oiii_5007_oii_err
    
    ml = ml_abundance(data,ratios=ratios,refmodel=2,debug=debug)

    zout = ml.ML_Z
    zout_err = ml.ML_Z_err

    ploterror, zinput, zout, zout_err, ps=3, xsty=3, ysty=3, $
      xrange=[7.2,9.5], yrange=[7.2,9.5], syms=2.0, errthick=2.0
;   plot, zinput, zout, ps=4, xsty=3, ysty=3, xrange=[7.2,9.5], yrange=[7.2,9.5], syms=2.0
    djs_oplot, !x.crange, !y.crange, line=0, thick=2.0

stop    
    
return
end
    
