;+
; NAME:
;   ISEDFIT_CHI2GRID_QAPLOT
;
; PURPOSE:
;   Generate quality-assurance plots from ISEDFIT output.
;
; INPUTS:
;   None required.
;
; OPTIONAL INPUTS:
;   datapath      - should match ISEDFIT [default CWD()] 
;   isedfitprefix - should match ISEDFIT (default 'isedfit')
;
; KEYWORD PARAMETERS:
;   maxold     - see ISEDFIT
;   make_png   - generate PNG output
;   postscript - generate a single postscript output file
;                (stronger than MAKE_PNG)
;
; INPUTS/OUTPUTS:
;   result      - see ISEDFIT
;   result_info - see ISEDFIT
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   The cosmological parameters are hard-wired to match
;   ISEDFIT_MODELS.  If POSTSCRIPT=1 then a single, large
;   postscript file is generated for all the objects.  If
;   MAKE_PNG=1 then individual PNG files are generated.
;
; EXAMPLES:
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2005 Feb 12, U of A
;   jm05aug04uofa - some small changes incorporated
;   jm06mar06uofa - documented; major changes to support the
;                   latest version of ISEDFIT
;   jm06mar23uofa - added FNU, FLAMBDA, MEDIANPLOT, and MAXOLD
;                   keywords  
;   jm07mar20nyu  - further developments
;   jm07jun27nyu  - significantly streamlined; the best-fitting
;                   SED models are now read by ISEDFIT_RESTORE() 
;
; Copyright (C) 2005-2007, John Moustakas
; 
; This program is free software; you can redistribute it and/or modify 
; it under the terms of the GNU General Public License as published by 
; the Free Software Foundation; either version 2 of the License, or
; (at your option) any later version. 
; 
; This program is distributed in the hope that it will be useful, but 
; WITHOUT ANY WARRANTY; without even the implied warranty of
; MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
; General Public License for more details. 
;-

pro render_ageplot, chi21, age1, ebv=ebv1, birthrate=birthrate1, $
  isedfit=isedfit, galaxy=galaxy

    title = repstr(strtrim(galaxy,2),'_',' ')+', z = '+$
      string(isedfit.zobj,format='(F6.4)')
    
    good = where(chi21 lt 1E6,ngood)
    chi2 = chi21[good]
    age = age1[good]
    ebv = ebv1[good]
    birthrate = birthrate1[good]
    
;   mass = chi2grid[*,*,igal].mass
    
    dustsf = where((ebv gt 0.0) and (birthrate gt 0.0),ndustsf)
    dustnosf = where((ebv gt 0.0) and (birthrate eq 0.0),ndustnosf)
    nodustsf = where((ebv eq 0.0) and (birthrate gt 0.0),nnodustsf)
    nodustnosf = where((ebv eq 0.0) and (birthrate eq 0.0),nnodustnosf)
    
;   nodust = where((ebv eq 0.0),nnodust)
;   dustysfr = where((ebv ge 0.0) and (birthrate gt 0.1),ndustysfr)
;   dustylowsfr = where((ebv gt 0.0) and (birthrate gt 0.0)
;     (birthrate lt 0.1),ndustylowsfr)

    if (min(chi2) gt 1E3) then chi2max = 1E4 else $
      chi2max = 2E3<max(chi2)
    djs_plot, [0], [0], /nodata, position=pos, xsty=3, ysty=3, $
      /xlog, /ylog, xtitle='Age (Gyr)', ytitle='\chi^{2}', $
      xrange=[0.1,max(age1)], yrange=[min(chi2),chi2max], title=title

    if (dustnosf[0] ne -1) then djs_oplot, age[dustnosf[sort(age[dustnosf])]], $
      chi2[dustnosf[sort(age[dustnosf])]], psym=symcat(7), color='blue'
    if (dustsf[0] ne -1) then djs_oplot, age[dustsf[sort(age[dustsf])]], $
      chi2[dustsf[sort(age[dustsf])]], psym=symcat(6,thick=!p.thick), color='grey'
    if (nodustsf[0] ne -1) then djs_oplot, age[nodustsf[sort(age[nodustsf])]], $
      chi2[nodustsf[sort(age[nodustsf])]], psym=symcat(15), color='orange'
    if (nodustnosf[0] ne -1) then djs_oplot, age[nodustnosf[sort(age[nodustnosf])]], $
      chi2[nodustnosf[sort(age[nodustnosf])]], psym=symcat(4,thick=!p.thick), $
      color='dark green', symsize=1.5
;   if (ndustynosfr gt 0L) then djs_oplot, age[dustynosfr[sort(age[dustynosfr])]], $
;     chi2[dustynosfr[sort(age[dustynosfr])]], psym=7, color='blue'
;   if (ndustylowsfr gt 0L) then djs_oplot, age[dustylowsfr[sort(age[dustylowsfr])]], $
;     chi2[dustylowsfr[sort(age[dustylowsfr])]], psym=-6, color='dark green', symsize=2
    im_legend, ['E(B-V)>0,b=0','E(B-V)>0,b>0','E(B-V)=0,b>0','E(B-V)=0,b=0'], $
      /left, /top, box=0, psym=[7,6,15,4], /clear, symthick=6, $
      color=['blue','grey','orange','dark green'], charsize=1.4, $
      spacing=1.5
    plots, isedfit.age, isedfit.chi2, psym=symcat(9,thick=6), $
      symsize=2.5, color=djs_icolor('red')
;   oploterror, isedfit.age, isedfit.chi2, isedfit.age_err, isedfit.chi2*0.0, $
;     psym=6, symsize=2.0, color=djs_icolor('red'), errcolor=djs_icolor('red')

return
end

pro isedfit_chi2grid_qaplot, paramfile, isedfit, params=params, iopath=iopath, $
  galaxy=galaxy1, outprefix=outprefix, isedfit_sfhgrid_dir=isedfit_sfhgrid_dir, $
  psfile=psfile1, index=index, clobber=clobber

    if (n_elements(paramfile) eq 0L) and (n_elements(params) eq 0) then begin
       doc_library, 'isedfit_chi2grid_qaplot'
       return
    endif

; read the parameter file and parse to get the relevant path and
; filenames  
    if (n_elements(iopath) eq 0) then iopath = './'
    if (n_elements(params) eq 0) then params = $
      read_isedfit_paramfile(paramfile)
    nsfhgrid = n_elements(params.sfhgrid)
    nredcurve = n_elements(params.redcurve)
    if (nsfhgrid gt 1) or (nredcurve gt 1) then begin
       for ii = 0, nsfhgrid-1 do begin
          newparams1 = struct_trimtags(params,except='sfhgrid')
          newparams1 = struct_addtags(newparams1,{sfhgrid: params.sfhgrid[ii]})
          for jj = 0, nredcurve-1 do begin
             newparams2 = struct_trimtags(newparams1,except='redcurve')
             newparams2 = struct_addtags(newparams2,{redcurve: params.redcurve[jj]})
             isedfit_chi2grid_qaplot, params=newparams2, iopath=iopath, galaxy=galaxy1, $
               outprefix=outprefix, psfile=psfile1, index=index, clobber=clobber
          endfor
       endfor 
       return
    endif

; read the isedfit output    
    junk = isedfit_restore(paramfile,isedfit,params=params,$
      iopath=iopath,index=index,outprefix=outprefix,silent=silent,$
      /nomodels)
    ngal = n_elements(isedfit)

; allow the user to overwrite PSFILE
    fp = isedfit_filepaths(params,outprefix=outprefix,iopath=iopath,$
      ngalaxy=ngal,ngalchunk=ngalchunk,galchunksize=galchunksize,$
      isedfit_sfhgrid_dir=isedfit_sfhgrid_dir)
    if (n_elements(psfile1) eq 0) then $
      psfile = iopath+strtrim(fp.qaplot_chi2grid_psfile,2) else $
        psfile = psfile1
    if file_test(psfile+'.gz',/regular) and $
      (keyword_set(clobber) eq 0) then begin
       splog, 'Output file '+psfile+' exists; use /CLOBBER'
       return
    endif

; generic galaxy names, unless provided
    if (n_elements(galaxy1) eq 0L) then begin
       fmt = '(I'+string(5L,format='(I0)')+'.'+string(5L,format='(I0)')+')'
       galaxy = 'Galaxy_'+string(lindgen(ngal),format=fmt)
    endif else begin
       if (n_elements(index) eq 0L) then begin
          galaxy = galaxy1
          if (n_elements(galaxy) ne ngal) then begin
             splog, 'GALAXY and ISEDFIT output must have '+$
               'the same number of elements'
             return
          endif
       endif else galaxy = galaxy1[index]
    endelse

; read all the models    
    nchunk = n_elements(fp.isedfit_models_chunkfiles)
    for ichunk = 0, nchunk-1 do begin
       chunkfile = fp.modelspath+fp.isedfit_models_chunkfiles[ichunk]+'.gz'
       if (keyword_set(silent) eq 0) then splog, 'Reading '+chunkfile
       modelgrid1 = mrdfits(chunkfile,1,/silent)
       modelgrid1 = struct_trimtags(temporary(modelgrid1),$
         except=['MODELMAGGIES'])
       if (n_elements(modelgrid) eq 0) then $
         modelgrid = temporary(modelgrid1) else $
           modelgrid = [temporary(modelgrid),temporary(modelgrid1)]
    endfor
    nmodel = n_elements(modelgrid)
    nage = n_elements(modelgrid[0].modelage)

; some plotting variables
    age = modelgrid.modelage
    birthrate = modelgrid.modelbirthrate
    tau = rebin(reform(modelgrid.tau,1,nmodel),nage,nmodel)
    ebv = rebin(reform(modelgrid.ebv,1,nmodel),nage,nmodel)
    Z = rebin(reform(modelgrid.Z,1,nmodel),nage,nmodel)

;   lnlike = 
;   contour, exp(-(chi2grid[*,*,0].chi2-isedfit[0].chi2)/2), age, chi2grid[*,*,0].mass, xsty=3, ysty=3        
;   contour, exp(-(chi2grid[*,*,0].chi2-isedfit[0].chi2)/2), age, chi2grid[*,*,0].mass, xsty=3, ysty=3        

; now loop on each galaxy chunk and read the relevant chi2 grid
    im_plotconfig, 8, pos, psfile=psfile, charsize=1.6, ymargin=[0.9,1.2]
;   for gchunk = 0, 0 do begin
    for gchunk = 0, ngalchunk-1 do begin
       g1 = gchunk*galchunksize
       g2 = ((gchunk*galchunksize+galchunksize)<ngal)-1L
       gnthese = g2-g1+1L
       gthese = lindgen(gnthese)+g1

       chi2gridfile = fp.modelspath+fp.chi2grid_gchunkfiles[gchunk]+'.gz'
       splog, 'Reading '+chi2gridfile
       chi2grid = mrdfits(chi2gridfile,1,/silent,columns=['mass','chi2'])
       chi2grid = reform(chi2grid,nage,nmodel,gnthese)
;      for igal = 55, 55 do begin            
       for igal = 0L, gnthese-1L do begin            
          render_ageplot, chi2grid[*,*,igal].chi2, age, ebv=ebv, $
            birthrate=birthrate, isedfit=isedfit[gthese[igal]], $
            galaxy=galaxy[gthese[igal]]
       endfor
    endfor
    im_plotconfig, psfile=psfile, /psclose, /gzip
    
return
end


;; now read the model chunks we need
;       allchunks = isedfit[gthese].chunkindx
;       chunks = allchunks[uniq(allchunks,sort(allchunks))]
;       nchunk = n_elements(chunks)
;       for ichunk = 0, nchunk-1 do begin
;          these = where(chunks[ichunk] eq allchunks,nthese)
;          if (nthese ne 0L) and (chunks[ichunk] ge 0) then begin
;             chunkfile = strtrim(fp.modelspath+fp.isedfit_models_chunkfiles[chunks[ichunk]],2)+'.gz'
;             if (not keyword_set(silent)) then $
;               splog, 'Reading '+chunkfile
;             modelgrid1 = mrdfits(chunkfile,1,/silent)
;             modelgrid1 = struct_trimtags(temporary(modelgrid1),$
;               except=['MODELMAGGIES'])
;             if (n_elements(modelgrid) eq 0) then $
;               modelgrid = temporary(modelgrid1) else $
;               modelgrid = [temporary(modelgrid),temporary(modelgrid1)]
;          endif 
;       endfor 
