;+
; NAME:
;   VERIFY_SPECTRAL_INDICES
; PURPOSE:
;   Compare the output of SPECTRAL_INDICES() against BC03. 
; OUTPUTS: 
;   Generates a QAplot.
; MODIFICATION HISTORY:
;   J. Moustakas, 2003 Oct 3, U of A
;-

pro verify_spectral_indices

    common verify, res
    
    bc = im_read_bc03(bc03_extras=bce);,age=findgen(10)+1.0)
    nage = n_elements(bc.age)

    if (n_elements(res) eq 0) then begin
       for ii = 0L, nage-1 do begin
          print, format='("Spectral index ",I3,"/",I3,".",A1,$)', $
            ii+1, nage, string(13b)
          res1 = spectral_indices(bc.wave,bc.flux[*,ii],$
            ivar=bc.flux[*,ii]*0.0+1.0,debug=debug,/silent) ;,/nolick)
          res1 = struct_addtags({age: bc.age[ii]},res1)
          if (ii eq 0) then res = res1 else res = [res,res1]
       endfor
       print
    endif

    result = res[1:nage-1]

    mynames = [$
      'D4000',       $
      'D4000_narrow',$
      'Lick_G4300'    ,$
      'Lick_FE4531'   ,$
      'Lick_HB'       ,$
      'Lick_FE5015'   ,$
      'Lick_MG1'      ,$
      'Lick_MG2'      ,$
      'Lick_MGB'      ,$
      'Lick_FE5270'   ,$
      'Lick_FE5335'   ,$
      'Lick_FE5406'   ,$
      'Lick_FE5709'   ,$
      'Lick_FE5782'   ,$
      'Lick_NaD'      ,$
      'Lick_Hd_A'     ,$
      'Lick_Hg_A'      $
      ]

    bcnames = [$
      'B_4000_'  ,$
      'B4_VN'    ,$
      'G4300'    ,$
      'FE4531'   ,$
      'H_BETA'   ,$
      'FE5015'   ,$
      'MG_1'     ,$
      'MG_2'     ,$
      'MG_B'     ,$
      'FE5270'   ,$
      'FE5335'   ,$
      'FE5406'   ,$
      'FE5709'   ,$
      'FE5782'   ,$
      'NA_D'     ,$
      'H_DELTA_A',$
      'H_GAMMA_A' $
      ]

    plotsym, 0, 1, fill=0, thick=3.0
    im_plotconfig, 6, pos, psfile=getenv('IMPRO_DIR')+$
      '/etc/verify_indices.ps', $
      xmargin=[1.4,0.4], yspace=1.0, width=6.7
    
    for jj = 0, n_elements(mynames)-1 do begin

       myindx = tag_indx(result,mynames[jj])
       bcindx = tag_indx(bce,bcnames[jj])

       xx = reform((result.(myindx))[0,*])
       yy = bce.(bcindx)
       resid = 100.0*(yy-xx)/(0.5*(yy+xx))

       xrange = [min(xx)<min(yy),max(xx)>max(yy)]
       plot, xx, yy, ps=8, xsty=3, ysty=3, position=pos[*,0], $
         xrange=xrange, yrange=xrange, xtitle=strupcase(mynames[jj])+' (Moustakas)', $
         ytitle=strupcase(mynames[jj])+' (BC03)'
       djs_oplot, !x.crange, !y.crange, line=0, thick=4.0, color='red'

       plot, result.age, resid, /xlog, ps=8, /noerase, xsty=3, ysty=3, $
         position=pos[*,1], yrange=im_max(abs(resid),sigrej=3.0)*[-1.5,1.5], $
         ytitle='Residuals (%)', xtitle='Age (yr)'
       djs_oplot, 10^!x.crange, [0,0], line=0, thick=4.0, color='red'
       
    endfor
    im_plotconfig, /psclose

return
end    
