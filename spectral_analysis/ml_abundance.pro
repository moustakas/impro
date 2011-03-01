;+
; NAME:
;       ML_ABUNDANCE()
;
; PURPOSE:
;       Compute the gas-phase abundance of a galaxy or HII region
;       using a maximum likelihood technique.
;
; CALLING SEQUENCE:
;
;
; INPUTS:
;       data - linear or logarithmic emission-line fluxes
;              [NRATIO,NOBJECT] 
;
; OPTIONAL INPUTS:
;       ratios   - emission-line ratios to use when constructing the  
;                  maximum likelihood surface
;
; KEYWORD PARAMETERS:
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;
; PROCEDURES USED:
;
; COMMENTS:
;
; EXAMPLES:
;
; MODIFICATION HISTORY:
;       J. Moustakas, 2004 Sep 07, U of A - non-functioning form 
;       jm04jun28uofa - developed, debuged, and documented
;-

function ml_abundance, fluxratios, invvar, rationames, Z=Z, U=U, $
  nmonte=nmonte, log=log, debug=debug

    Zsun = 8.70
    
    ndim = size(fluxratios,/n_dimension)
    ndimvar = size(invvar,/n_dimension)
    nratios = n_elements(rationames)

    if ((ndim ne 1L) and (ndim ne 2L)) or ((ndimvar ne 1L) and (ndimvar ne 2L)) or $
      (nratios eq 0L) then begin
       print, 'Syntax - result = ml_abundance(fluxratios,invvar,rationames,Z=,U=,/log,,/debug)'
       return, -1L
    endif

    dims = size(fluxratios,/dimension)
    dimsvar = size(invvar,/dimension)

    if (ndim eq 1L) then begin
       if (dimsvar[0] ne dims[0]) then begin
          print, 'FLUXRATIOS and INVVAR must have the same dimensions.'
          return, -1L
       endif
       nobj = dims[0]
    endif else begin
       if (dimsvar[0] ne dims[0]) and (dimsvar[1] ne dims[1]) then begin
          print, 'FLUXRATIOS and INVVAR must have the same dimensions.'
          return, -1L
       endif
       nobj = dims[1]
    endelse

    nlog = n_elements(log)
    if (nlog eq 0L) then nolog = 1L else nolog = 0L

; read the theoretical flux-ratio grid [NMODEL,NZ,NU]; 

    ngrids = 8L ; <-- hard-wired by what Lisa has done

    for igrid = 0L, ngrids-1L do begin
       model1 = read_kewley_grids(model=igrid+1L,Z=Zmodel,U=Umodel,$
         /oversamp,nolog=nolog)
       if (igrid eq 0L) then model = model1 else model = [ [ [model] ], [ [model1] ] ]
    endfor

    OHmodel = alog10(Zmodel)+Zsun
    
    nZmodel = n_elements(Zmodel)
    nUmodel = n_elements(Umodel)
    nmodel = nZmodel*nUmodel*ngrids
    
; subscript the model data to the flux ratios of interest

    modelratios = make_array(nratios,nmodel)
    for iratio = 0L, nratios-1L do modelratios[iratio,*] = $
      (struct_trimtags(model,select=rationames)).(iratio)
    
;   modelratios = make_array(nratios,nZmodel,nUmodel,ngrids)
;   for iratio = 0L, nratios-1L do modelratios[iratio,*,*,*] = $
;     (struct_trimtags(model,select=rationames)).(iratio)

; initialize variables that will speed up the Monte Carlo calculation...

; initialize the output chi2 surface arrays and results

    result = {$
      Zml_object:          '', $
      Zml_chi2surface:     fltarr(nmodel), $
      Zml_12logoh:         0.0, $
      Zml_12logoh_median:  0.0, $
      Zml_12logoh_err:     0.0, $
      Zml_12logoh_chi2min: 0.0, $
      Zml_U:               0.0, $
      Zml_U_median:        0.0, $
      Zml_U_err:           0.0, $
      Zml_U_chi2min:       0.0, $
      Zml_grid:            0L,  $
      Zml_grid_median:     0.0, $
      Zml_grid_err:        0.0, $
      Zml_grid_chi2min:    0.0, $
      Zml_chi2min:         0.0, $
      Zml_chi2indx:         0L}
;     Zml_errcode:         -1.0}
    result = replicate(result,nobj)

    multisurface = make_array(nUmodel,nZmodel,ngrids,/float,/nozero)
    
; loop on each object    
    
    for iobj = 100L, 150-1L do begin
;   for iobj = 0L, nobj-1L do begin

       goodfluxes = where(fluxratios[*,iobj] gt 0.0,ngoodfluxes)
       dof = ngoodfluxes - 1L ; degrees of freedom
       
       for imodel = 0L, nmodel-1L do begin

          if (((imodel+1L) mod 5000L) eq 0L) then print, format='("Object = ",I4,"/",I4,", '+$
            'Model = ",I5,"/",I5,".",A4,$)', iobj+1L, nobj, imodel+1L, nmodel, string(13b)

          dchi2 = total(invvar[*,iobj]*(fluxratios[*,iobj]-modelratios[*,imodel])^2)
          result[iobj].Zml_chi2surface[imodel] = dchi2

       endfor

; minimize by marginilizing over each dimension

       result[iobj].Zml_chi2min = min(result[iobj].Zml_chi2surface,chi2indx)
       result[iobj].Zml_chi2indx = chi2indx

       bestindx = array_indices(multisurface,chi2indx)

;      result[iobj].Zml_12logoh = alog10(Zmodel[bestindx[1]])+Zsun
;      result[iobj].Zml_U       = Umodel[bestindx[0]]
;      result[iobj].Zml_grid    = bestindx[2]

       OHchi2 = total(total(reform(result[iobj].Zml_chi2surface,nUmodel,nZmodel,ngrids),1),2)
       OHbest = find_nminima(OHchi2,OHmodel,dof=dof,nfind=1L,ypeak=OHchi2min,doplot=debug)
       if keyword_set(debug) then cc = get_kbrd(1)
       
       Uchi2 = total(total(reform(result[iobj].Zml_chi2surface,nUmodel,nZmodel,ngrids),2),2)
       Ubest = find_nminima(Uchi2,Umodel,dof=dof,nfind=1L,ypeak=Uchi2min,doplot=debug)
       if keyword_set(debug) then cc = get_kbrd(1)

       gridchi2 = total(total(reform(result[iobj].Zml_chi2surface,nUmodel,nZmodel,ngrids),1),1)
       gridbest = find_nminima(gridchi2,lindgen(ngrids),dof=dof,nfind=1L,ypeak=gridchi2min,doplot=debug)
       if keyword_set(debug) then cc = get_kbrd(1)

; store the results       
       
       result[iobj].Zml_12logoh         = OHbest
       result[iobj].Zml_12logoh_chi2min = OHchi2min
       result[iobj].Zml_U               = Ubest
       result[iobj].Zml_U_chi2min       = Uchi2min
       result[iobj].Zml_grid            = gridbest
       result[iobj].Zml_grid_chi2min    = gridchi2min
       
;      plot, OHmodel, OHchi2, xsty=3, ysty=3
;      plot, Umodel, Uchi2, xsty=3, ysty=3
;      plot, lindgen(ngrids), gridschi2, xsty=3, ysty=3
       
    endfor

return, result
end    
