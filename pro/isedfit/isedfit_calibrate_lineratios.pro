;+
; NAME:
;   ISEDFIT_CALIBRATE_LINERATIOS
; PURPOSE:
;   Empirically calibrate the forbidden emission-line ratio variations
;   for use in iSEDfit. 
; COMMENTS:
; MODIFICATION HISTORY:
;   J. Moustakas, 2013 Jul 11, Siena
;-

function read_pegase_hii
; HII.dat data file from Pegase
    datafile = getenv('IM_RESEARCH_DIR')+'/synthesis/PEGASE-HR/data/external/HII.dat'
    txt = djs_readilines(datafile,nhead=100)
    txt = txt[where(strcompress(txt,/remove) ne '')]
    nline = n_elements(txt)
    data = replicate({name: '', wave: 0.0, ratio: 0.0},nline)
    for ii = 0, nline-1 do begin
       line = strsplit(txt[ii],' ',/extract)
       data[ii].name = strtrim(line[2],2)
       data[ii].wave = float(line[0])
       data[ii].ratio = float(line[1])
    endfor
return, data
end

pro oplot_peg, linewave, xx=xx, yy=yy
    common isedfit_calibrate, oiiihbaxis, iz06, g99, peg, sdss, hii, atlas
    color = 'dodger blue' & psym = 16
    get_element, peg.wave, 5007.0, is5007
    get_element, peg.wave, linewave, thisline
    delvarx, xx, yy
    xx = alog10([peg[is5007].ratio])
    yy = alog10([peg[thisline].ratio])
    djs_oplot, xx, yy, symsize=2.5, psym=symcat(psym), $
      color=im_color(color)
return
end
pro oplot_hii, oii=oii, nii=nii, sii6716=sii6716, $
  sii6731=sii6731, xx=xx, yy=yy
    common isedfit_calibrate
    color = 'forest green' & psym = 6
    xx = hii.oiiihb
    if keyword_set(oii) then yy = hii.oiihb
    if keyword_set(nii) then yy = hii.niihb
    if keyword_set(sii6716) then yy = hii.sii6716hb
    if keyword_set(sii6731) then yy = hii.sii6731hb
    djs_oplot, xx, yy, symsize=0.5, psym=symcat(psym), $
      color=im_color(color)
return
end
pro oplot_atlas, oii=oii, nii=nii, sii6716=sii6716, $
  sii6731=sii6731, xx=xx, yy=yy
    common isedfit_calibrate
    color = 'tomato' & psym = 16
    xx = atlas.oiiihb
    if keyword_set(oii) then yy = atlas.oiihb
    if keyword_set(nii) then yy = atlas.niihb
    if keyword_set(sii6716) then yy = atlas.sii6716hb
    if keyword_set(sii6731) then yy = atlas.sii6731hb
    djs_oplot, xx, yy, symsize=1.0, psym=symcat(psym), $
      color=im_color(color)
return
end
pro oplot_g99, linewave, xx=xx, yy=yy
    common isedfit_calibrate
    color = 'orange' & psym = 15
    allreg = strtrim(g99.name,2)
    reg = allreg[uniq(allreg,sort(allreg))]
    delvarx, xx, yy
    for ii = 0, n_elements(reg)-1 do begin
       ww = where(reg[ii] eq allreg)
       get_element, g99[ww].wave*10, 5007.0, is5007
       get_element, g99[ww].wave*10, linewave, thisline
       if strtrim(g99[ww[is5007]].l_f3_1,2) eq '' and $
         strtrim(g99[ww[thisline]].l_f3_1,2) eq '' then begin
          xx1 = alog10(g99[ww[is5007]].f3_1)
          yy1 = alog10(g99[ww[thisline]].f3_1)
          oploterror, xx1, yy1, g99[ww[is5007]].e_f3_1/g99[ww[is5007]].f3_1/alog(10), $
            g99[ww[thisline]].e_f3_1/g99[ww[thisline]].f3_1/alog(10), $
            psym=symcat(psym), color=im_color(color), symsize=2.5
          if n_elements(xx) eq 0 then xx = xx1 else xx = [xx,xx1]
          if n_elements(yy) eq 0 then yy = yy1 else yy = [yy,yy1]
       endif
    endfor
return
end
pro oplot_iz06, linewave, xx=xx, yy=yy
    common isedfit_calibrate
    color = 'purple' & psym = 14

    tag = tag_indx(iz06,'f'+linewave)
    etag = tag_indx(iz06,'e_f'+linewave)
    reftag = tag_indx(iz06,'f4959')
    refetag = tag_indx(iz06,'e_f4959')

    ww = where(iz06.(tag) gt 0.0 and iz06.(reftag) gt 0.0,nww)
    if nww ne 0 then begin
       xx = alog10(2.98*iz06[ww].(reftag))
       yy = alog10(iz06[ww].(tag))
       xxerr = iz06[ww].(refetag)/iz06[ww].(reftag)/alog(10)
       yyerr = iz06[ww].(etag)/iz06[ww].(tag)/alog(10)
       oploterror, xx, yy, xxerr, yyerr, psym=symcat(psym), $
         color=im_color(color), symsize=0.7
    endif
return
end

pro oplot_calibration, out
; overplot the calibration
    common isedfit_calibrate
    djs_oplot, oiiihbaxis, poly(oiiihbaxis,out.coeff), line=0, thick=8
    djs_oplot, oiiihbaxis, poly(oiiihbaxis,out.coeff)+out.scatter, line=5, thick=6
    djs_oplot, oiiihbaxis, poly(oiiihbaxis,out.coeff)-out.scatter, line=5, thick=6
return
end
function read_mysdss
; read and parse my SDSS catalog
    common isedfit_calibrate
    snrcut = 5.0
    allsdss = read_vagc_garching(/ispecline)

; remove AGN and require S/N>5 on all emission lines of interest 
    sclass = iclassification(allsdss,snrcut=snrcut,ratios=sratios)
    keep = where($
      strmatch(sratios.final_class,'*AGN*') eq 0 and $
;     strmatch(sratios.final_class,'*SF*') and $
      allsdss.oiii_5007[0]/allsdss.oiii_5007[1] gt snrcut and $
      allsdss.oii_3727[0]/allsdss.oii_3727[1] gt snrcut and $
      allsdss.h_beta[0]/allsdss.h_beta[1] gt snrcut and $
      allsdss.h_alpha[0]/allsdss.h_alpha[1] gt snrcut and $
      allsdss.nii_6584[0]/allsdss.nii_6584[1] gt snrcut and $
      allsdss.sii_6716[0]/allsdss.sii_6716[1] gt snrcut and $
      allsdss.sii_6731[0]/allsdss.sii_6731[1] gt snrcut,nsdss)
    cat = allsdss[keep]

    lineratio, cat, 'OIII_5007', 'H_BETA', 'OII_3727', $
      'H_BETA', oiiihb, oiiihb_err, oiihb, oiihb_err, snrcut=snrcut
    lineratio, cat, 'OIII_5007', 'H_BETA', 'NII_6584', $
      'H_ALPHA', oiiihb, oiiihb_err, niiha, niiha_err, snrcut=snrcut
;   lineratio, cat, 'OIII_5007', 'H_BETA', ['SII_6716','SII_6731'], $
;     'H_ALPHA', oiiihb, oiiihb_err, siiha, siiha_err, snrcut=snrcut
    lineratio, cat, 'OIII_5007', 'H_BETA', 'SII_6716', $
      'H_ALPHA', oiiihb, oiiihb_err, sii6716ha, sii6716ha_err, snrcut=snrcut
    lineratio, cat, 'OIII_5007', 'H_BETA', 'SII_6731', $
      'H_ALPHA', oiiihb, oiiihb_err, sii6731ha, sii6731ha_err, snrcut=snrcut

    niihb = niiha + alog10(2.86) ; case B, no dust
    sii6716hb = sii6716ha + alog10(2.86)
    sii6731hb = sii6731ha + alog10(2.86)
    
    out = {oiiihb: oiiihb, oiiihb_err: oiiihb_err, $
      oiihb: oiihb, oiihb_err: oiihb_err, niihb: niihb, $
      niihb_err: niiha_err, sii6716hb: sii6716hb, $
      sii6716hb_err: sii6716ha_err, sii6731hb: sii6731hb, $
      sii6731hb_err: sii6731ha_err}

return, out
end    
function read_myhii
    common isedfit_calibrate
    snrcut = 5.0
    allhiiinfo = read_hiiregions(linefit=allhii)

    keep = where($
      allhii.oiii_5007[0]/allhii.oiii_5007[1] gt snrcut and $
      allhii.oii_3727[0]/allhii.oii_3727[1] gt snrcut and $
      allhii.h_beta[0]/allhii.h_beta[1] gt snrcut and $
      allhii.h_alpha[0]/allhii.h_alpha[1] gt snrcut and $
      allhii.nii_6584[0]/allhii.nii_6584[1] gt snrcut and $
      allhii.sii_6716[0]/allhii.sii_6716[1] gt snrcut and $
      allhii.sii_6731[0]/allhii.sii_6731[1] gt snrcut,nhii)
    cat = allhii[keep]
    hiiinfo = allhiiinfo[keep]

    lineratio, cat, 'OIII_5007', 'H_BETA', 'OII_3727', $
      'H_BETA', oiiihb, oiiihb_err, oiihb, oiihb_err, snrcut=snrcut
    lineratio, cat, 'OIII_5007', 'H_BETA', 'NII_6584', $
      'H_ALPHA', oiiihb, oiiihb_err, niiha, niiha_err, snrcut=snrcut
;   lineratio, cat, 'OIII_5007', 'H_BETA', ['SII_6716','SII_6731'], $
;     'H_ALPHA', oiiihb, oiiihb_err, siiha, siiha_err, snrcut=snrcut
    lineratio, cat, 'OIII_5007', 'H_BETA', 'SII_6716', $
      'H_ALPHA', oiiihb, oiiihb_err, sii6716ha, sii6716ha_err, snrcut=snrcut
    lineratio, cat, 'OIII_5007', 'H_BETA', 'SII_6731', $
      'H_ALPHA', oiiihb, oiiihb_err, sii6731ha, sii6731ha_err, snrcut=snrcut

    niihb = niiha + alog10(2.86) ; case B, no dust
    sii6716hb = sii6716ha + alog10(2.86)
    sii6731hb = sii6731ha + alog10(2.86)
    
    out = {oiiihb: oiiihb, oiiihb_err: oiiihb_err, $
      oiihb: oiihb, oiihb_err: oiihb_err, niihb: niihb, $
      niihb_err: niiha_err, sii6716hb: sii6716hb, $
      sii6716hb_err: sii6716ha_err, sii6731hb: sii6731hb, $
      sii6731hb_err: sii6731ha_err}

return, out
end

function read_myatlas
    common isedfit_calibrate
    snrcut = 5.0
    allatlas = mrdfits(getenv('IM_PROJECTS_DIR')+'/isedfit/'+$
      'integrated_atlas_speclinefit_v3.1.fits.gz',1)

    aclass = iclassification(allatlas,snrcut=snrcut,ratios=aratios)
    keep = where($
      strmatch(aratios.final_class,'*AGN*') eq 0 and $
;     strmatch(aratios.final_class,'*SF*') and $
      allatlas.oiii_5007[0]/allatlas.oiii_5007[1] gt snrcut and $
      allatlas.oii_3727[0]/allatlas.oii_3727[1] gt snrcut and $
      allatlas.h_beta[0]/allatlas.h_beta[1] gt snrcut and $
      allatlas.h_alpha[0]/allatlas.h_alpha[1] gt snrcut and $
      allatlas.nii_6584[0]/allatlas.nii_6584[1] gt snrcut and $
      allatlas.sii_6716[0]/allatlas.sii_6716[1] gt snrcut and $
      allatlas.sii_6731[0]/allatlas.sii_6731[1] gt snrcut,natlas)
    splog, 'ATLAS ', natlas
    cat = allatlas[keep]

    lineratio, cat, 'OIII_5007', 'H_BETA', 'OII_3727', $
      'H_BETA', oiiihb, oiiihb_err, oiihb, oiihb_err, snrcut=snrcut
    lineratio, cat, 'OIII_5007', 'H_BETA', 'NII_6584', $
      'H_ALPHA', oiiihb, oiiihb_err, niiha, niiha_err, snrcut=snrcut
;   lineratio, cat, 'OIII_5007', 'H_BETA', ['SII_6716','SII_6731'], $
;     'H_ALPHA', oiiihb, oiiihb_err, siiha, siiha_err, snrcut=snrcut
    lineratio, cat, 'OIII_5007', 'H_BETA', 'SII_6716', $
      'H_ALPHA', oiiihb, oiiihb_err, sii6716ha, sii6716ha_err, snrcut=snrcut
    lineratio, cat, 'OIII_5007', 'H_BETA', 'SII_6731', $
      'H_ALPHA', oiiihb, oiiihb_err, sii6731ha, sii6731ha_err, snrcut=snrcut

    niihb = niiha + alog10(2.86) ; case B, no dust
    sii6716hb = sii6716ha + alog10(2.86)
    sii6731hb = sii6731ha + alog10(2.86)
    
    out = {oiiihb: oiiihb, oiiihb_err: oiiihb_err, $
      oiihb: oiihb, oiihb_err: oiihb_err, niihb: niihb, $
      niihb_err: niiha_err, sii6716hb: sii6716hb, $
      sii6716hb_err: sii6716ha_err, sii6731hb: sii6731hb, $
      sii6731hb_err: sii6731ha_err}

return, out
end
pro calibrate_ratio, xx, yy, weight=weight, name=name, wave=wave, $
  ncoeff=ncoeff, fitmedian=fitmedian, indx=indx, out=out
; initialize the output data structure
    out1 = {id: 0, name: name, wave: wave, coeff: fltarr(4), scatter: 0.0}
    if ncoeff eq 1 then begin
       out1.coeff[0] = djs_mean(yy)
       out1.scatter = (djsig(yy)>0.1)<0.3 ; dex
    endif else begin
       if keyword_set(fitmedian) then begin
          md = im_medxbin(xx,yy,0.3,weight=weight)
          out1.coeff[0:ncoeff-1] = poly_fit(md.medx,md.medy,ncoeff-1,yerror=sig)
          out1.scatter = (sig>0.1)<0.3
       endif else begin
          out1.coeff[0:ncoeff-1] = poly_fit(xx,yy,ncoeff-1,yerror=sig)
          out1.scatter = (sig>0.1)<0.3
       endelse
    endelse
; build the output structure
    if n_elements(out) eq 0 then begin
       out = out1
       indx = 0
       out[0].id = 1
    endif else begin
       out = [out,out1]
       indx = n_elements(out)-1
       out[indx].id = n_elements(out)
    endelse
return
end

pro isedfit_calibrate_lineratios 
; jm13aug01siena -

    common isedfit_calibrate
    path = getenv('IM_PROJECTS_DIR')+'/isedfit/'

; ---------------------------------------------------------------------------
; write out the nebular *continuum* spectrum for use with
; ISEDFIT_NEBULAR by stealing shamelessly from Starburst99; the
; spectrum varies just with the number of Lyman-continuum photons, so
; I can safely choose any age; the units of the continuum spectrum is
; erg/s/A per unit N(Lyc)
    sb99 = im_read_starburst99('SB99_PadovaAGB_Kroupa_Z0.020',path=path)
    this = 10 ; age index

; build an oversampling wavelength array near the spectral breaks, to
; ensure adequate sampling
    outwave = range(min(sb99.wave),max(sb99.wave),n_elements(sb99.wave)*2)/1D4 ; oversample
    sbreak = [0.0912,0.342,0.365,0.816,1.455,2.277,3.277,4.460]
    for ii = 0, n_elements(sbreak)-1 do outwave = [outwave,$
      10D^(range(alog10(sbreak[ii])-0.01,alog10(sbreak[ii])+0.01,100))]
    outwave = outwave[sort(outwave)]
    outwave = outwave*1D4
    cneb = {wave_neb: outwave, flam_neb: interpolate(sb99.flux_gas[*,this]/10D^sb99.nlyc[this],$
      findex(sb99.wave,outwave))}
;   zero = where(out.wave lt 912.0)
;   out.flam_neb[zero] = 0.0

;   djs_plot, out.wave/1D4, out.flam_neb, psym=-8, ysty=3, xsty=3, /ylog, xr=[0.0,1.0]
;   djs_plot, out.wave/1D4, out.flam_neb, psym=-8, ysty=3, xsty=3, /ylog, xr=[1.0,4.5]
;   djs_oplot, sb99.wave/1D4, sb99.flux_gas[*,this]/10D^sb99.nlyc[this], color='orange'
;   djs_oplot, sb99.wave/1D4, sb99.flux_gas[*,this+20]/10D^sb99.nlyc[this+20], color='green'

    outfile = getenv('IMPRO_DIR')+'/etc/isedfit_nebular_continuum.fits'
    im_mwrfits, cneb, outfile, /clobber

;; below are some notes from Pegase and Koleva+09 for computing the
;; nebular *continuum*
;    
;; total recombination coefficient alpha2(T)
;    alpha2 = 2.575D-13 ; Aller+84 [cm^3 s^-1]
;;   alpha2 = 2.616D-13 ; Pegase.2
;
;; assumed densities of He+ and He++ relative to H+, taken from
;; Koleva+09 
;    f_HeI = 0.0897    ; Pegase uses 0.095
;    f_HeII = 1.667D-4 ; Pegase assumes zero
;
;    light = im_light(/ang) ; [Angstrom/s]
;
;; eventually replace this bit of code with my own computations of the
;; total emission coefficient (see for example Koleva+09), but for now
;; just use Pegase's calculations
;;   peg = pegase_hii()
;;   fneb(i)=1.d-40*(g1+g2+HeI*g3+HeII*g4)*c/alphaTe        /lambda(i)**2
;    
;; total emission coefficient
;    gamma_tot = gamma_b + gamma_2p + gamma_HI + $
;      f_HeI*gamma_HeI + f_HeII*gamma_HeII
;
;; see equation (1) in Koleva+09    
;    flam_neb = gamma_tot*nlyc*light/wave^2/alpha2
    
; ---------------------------------------------------------------------------    
; now do the emission lines

; for reference, here are the emission lines in Pegase
;   10691.    0.021          CI       36
;    5199.    0.027          NI       40
;    2141.    0.013          NII      41
;    5755.    0.011          NII      42
;    6300.    0.126          OI       44
;    7330.    0.045          OII      46
;    1663.    0.009          OIII     47
;    4363.    0.012          OIII     48
;   12800.    0.082          NeII     50
;    2798.    0.07           MgII     52
;    4070.    0.043          SII      53
;   10330.    0.163          SII      56
;    6312.    0.026          SIII     57
;    9069.    2.715          SIII     58
;    7135.    0.241          ArIII    59 
;   21230.    0.01355        H2       60
;   16440.    0.03522        FeII     61

; the detailed nebular wavelengths (in air!) are derived from
; IM_IONSPEC, e.g., 
;   IDL> im_ionspec, 'C_III'
; ANNOT INDEX       VACWAVE       AIRWAVE         AIJ RELATIVE_EMISSIVITY ABSOLUTE_EMISSIVITY
; ----- ----- ------------- ------------- ----------- ------------------- -------------------
; 5-->1 [4,0]      977.0195      977.0195  1.7900e+09              22.342          2.8106e-24
; 4-->1 [3,0]      1908.702      1908.702   0.0051900              2374.6          2.9872e-22
; 3-->1 [2,0]      1910.735      1910.735      12.100              1581.2          1.9891e-22
; 4-->2 [3,1]      1265638.      1265293.  1.3700e-13          9.4530e-11          1.1892e-35
; 4-->3 [3,2]      1793767.      1793278.  2.3300e-06           0.0011344          1.4270e-28
; 3-->2 [2,1]      4298686.      4297515.  2.2700e-07          1.3185e-08          1.6587e-33

; See also http://www.pa.uky.edu/~peter/atomic, Peimbert et al. 2004,
; ApJS, 150, 431, and Stasinska 2007 for additional useful data.
    
; --------------------------------------------------
; read the emission-line data we need
    if n_elements(peg) eq 0 then peg = read_pegase_hii() ; Pegase
    if n_elements(g99) eq 0 then g99 = mrdfits(getenv('CATALOGS_DIR')+$ ; Garnett+99    
      '/99garnett/99garnett_table2.fits.gz',1)
    if n_elements(iz06) eq 0 then iz06 = im_read_vizier_tsv($ ; Izotov+06
      getenv('CATALOGS_DIR')+'/hiiregions/06izotov/2006_izotov.dat')
    if n_elements(sdss) eq 0 then sdss = read_mysdss()
    if n_elements(hii) eq 0 then hii = read_myhii()
    if n_elements(atlas) eq 0 then atlas = read_myatlas()

;   galev = rsex(path+'galev_emlines.dat')
;   ratio = [[galev.ratio_z0],[galev.ratio_z1],[galev.ratio_z2],[galev.ratio_z3]]>(1E-6)
;   ratio = alog10(ratio)>(-3)

; some plotting variables    
    oiiihbtitle = 'log ([OIII] \lambda5007 / H\beta)'
    oiiihbrange = [-1.3,1.1]
    oiiihbaxis = range(-1.0,1.0,100)
    
    psfile = path+'qa_calibrate_lineratios.ps'
    im_plotconfig, 0, pos, psfile=psfile, height=5.0, width=6.6, $
      xmargin=[1.4,0.5]

; --------------------
; CIII] 1909 (is this a doublet?)
    ytitle = 'log (CIII] \lambda1909 / H\beta)'
    yrange = [-2,1]
    djs_plot, [0], [0], /nodata, position=pos, xsty=1, ysty=1, $
      xrange=oiiihbrange, yrange=yrange, ytitle=textoidl(ytitle), $
      xtitle=textoidl(oiiihbtitle)
    oplot_peg, 1909, xx=pegxx, yy=pegyy
    oplot_g99, 1909, xx=g99xx, yy=g99yy
    calibrate_ratio, [pegxx,g99xx], [pegyy,g99yy], name='CIII]_1909',$
      wave=1908.702D, ncoeff=1, out=out, indx=indx
    oplot_calibration, out[indx]
    
; --------------------
; CII] 1335
    ytitle = 'log (CII] \lambda1335 / H\beta)'
    yrange = [-2,1]
    djs_plot, [0], [0], /nodata, position=pos, xsty=1, ysty=1, $
      xrange=oiiihbrange, yrange=yrange, ytitle=textoidl(ytitle), $
      xtitle=textoidl(oiiihbtitle)
    oplot_peg, 1335, xx=pegxx, yy=pegyy
    calibrate_ratio, pegxx, pegyy, name='CII]_1335',$
      wave=1335.707D, ncoeff=1, out=out, indx=indx
    oplot_calibration, out[indx]

; --------------------
; CII] 2325
    ytitle = 'log (CII] \lambda2325 / H\beta)'
    yrange = [-2,1]
    djs_plot, [0], [0], /nodata, position=pos, xsty=1, ysty=1, $
      xrange=oiiihbrange, yrange=yrange, ytitle=textoidl(ytitle), $
      xtitle=textoidl(oiiihbtitle)
    oplot_peg, 2325, xx=pegxx, yy=pegyy
    oplot_g99, 2325, xx=g99xx, yy=g99yy
    calibrate_ratio, [pegxx,g99xx], [pegyy,g99yy], name='CII]_2325',$
      wave=2324.689D, ncoeff=1, out=out, indx=indx
    oplot_calibration, out[indx]

; --------------------
; [OII] 2470
    ytitle = 'log (OII] \lambda2470 / H\beta)'
    yrange = [-2,1]
    djs_plot, [0], [0], /nodata, position=pos, xsty=1, ysty=1, $
      xrange=oiiihbrange, yrange=yrange, ytitle=textoidl(ytitle), $
      xtitle=textoidl(oiiihbtitle)
    oplot_peg, 2470, xx=pegxx, yy=pegyy
    oplot_g99, 2470, xx=g99xx, yy=g99yy
    calibrate_ratio, [pegxx,g99xx], [pegyy,g99yy], name='[OII]_2470',$
      wave=2470.341D, ncoeff=1, out=out, indx=indx
    oplot_calibration, out[indx]

; --------------------
; [OII] 3727; *doublet*
    ytitle = 'log (OII] \lambda3727 / H\beta)'
    yrange = [-2,1]
    sxx = sdss.oiiihb
    syy = sdss.oiihb
    im_hogg_scatterplot, sxx, syy, position=pos, xsty=1, ysty=1, $
      xrange=oiiihbrange, yrange=yrange, ytitle=textoidl(ytitle), $
      xtitle=textoidl(oiiihbtitle)
    oplot_peg, 3727, xx=pegxx, yy=pegyy
    oplot_hii, xx=hxx, yy=hyy, /oii
    oplot_atlas, xx=axx, yy=ayy, /oii
    weight = [axx*0+1.0/n_elements(axx),hxx*0+1.0/n_elements(hxx),$ ; equal weight
      sxx*0+1.0/n_elements(sxx)]
    calibrate_ratio, [axx,hxx,sxx], [ayy,hyy,syy], name='[OII]_3727',$
      wave=3727.4230D, ncoeff=4, out=out, indx=indx, /fitmedian
    oplot_calibration, out[indx]

; --------------------
; [NII] 6584; *doublet*
    ytitle = 'log (NII] \lambda6584 / H\beta)'
    yrange = [-2,1]
    sxx = sdss.oiiihb
    syy = sdss.niihb
    im_hogg_scatterplot, sxx, syy, position=pos, xsty=1, ysty=1, $
      xrange=oiiihbrange, yrange=yrange, ytitle=textoidl(ytitle), $
      xtitle=textoidl(oiiihbtitle)
    oplot_peg, 6584, xx=pegxx, yy=pegyy
    oplot_hii, xx=hxx, yy=hyy, /nii
    oplot_atlas, xx=axx, yy=ayy, /nii
    weight = [axx*0+1.0/n_elements(axx),hxx*0+1.0/n_elements(hxx),$ ; equal weight
      sxx*0+1.0/n_elements(sxx)]
    calibrate_ratio, [axx,hxx,sxx], [ayy,hyy,syy], name='[NII]_6584',$
      wave=6583.458D, ncoeff=3, out=out, indx=indx, /fitmedian
    oplot_calibration, out[indx]

; --------------------
; [SII] 6716; *doublet*
    ytitle = 'log (SII] \lambda6716 / H\beta)'
    yrange = [-2,1]
    sxx = sdss.oiiihb
    syy = sdss.sii6716hb
    im_hogg_scatterplot, sxx, syy, position=pos, xsty=1, ysty=1, $
      xrange=oiiihbrange, yrange=yrange, ytitle=textoidl(ytitle), $
      xtitle=textoidl(oiiihbtitle)
    oplot_peg, 6716, xx=pegxx, yy=pegyy
    oplot_hii, xx=hxx, yy=hyy, /sii6716
    oplot_atlas, xx=axx, yy=ayy, /sii6716
    weight = [axx*0+1.0/n_elements(axx),hxx*0+1.0/n_elements(hxx),$ ; equal weight
      sxx*0+1.0/n_elements(sxx)]
    calibrate_ratio, [axx,hxx,sxx], [ayy,hyy,syy], name='[SII]_6716',$
      wave=6716.440D, ncoeff=3, out=out, indx=indx, /fitmedian
    oplot_calibration, out[indx]

;; --------------------
;; [SII] 6731; *doublet*; in the low-density limit, [SII] 6716 is about
;; 30% stronger than [SII] 6731; therefore, to make sure that [SII]
;; vary together in ISEDFIT_NEBULAR, only store the [SII] 6716 line
;; here and build the doublet in ISEDFIT_NEBULAR
;    ytitle = 'log (SII] \lambda6731 / H\beta)'
;    yrange = [-2,1]
;    sxx = sdss.oiiihb
;    syy = sdss.sii6731hb
;    im_hogg_scatterplot, sxx, syy, position=pos, xsty=1, ysty=1, $
;      xrange=oiiihbrange, yrange=yrange, ytitle=textoidl(ytitle), $
;      xtitle=textoidl(oiiihbtitle)
;    oplot_peg, 6731, xx=pegxx, yy=pegyy
;    oplot_hii, xx=hxx, yy=hyy, /sii6731
;    oplot_atlas, xx=axx, yy=ayy, /sii6731
;    weight = [axx*0+1.0/n_elements(axx),hxx*0+1.0/n_elements(hxx),$ ; equal weight
;      sxx*0+1.0/n_elements(sxx)]
;    calibrate_ratio, [axx,hxx,sxx], [ayy,hyy,syy], name='[SII]_6731',$
;      wave=6730.815D, ncoeff=3, out=out, indx=indx, /fitmedian
;    oplot_calibration, out[indx]

;; check that we get ~1.35 for the [SII] ratio, i.e., low-density 
;    get_element, out.wave, 6716, is16
;    get_element, out.wave, 6731, is31
;    print, median(10^(poly(oiiihbaxis,out[is16].coeff)-$
;      poly(oiiihbaxis,out[is31].coeff)))    
    
; --------------------
; [Ne III] 3869
    ytitle = 'log ([Ne III] \lambda3869 / H\beta)'
    yrange = [-2,1]
    djs_plot, [0], [0], /nodata, position=pos, xsty=1, ysty=1, $
      xrange=oiiihbrange, yrange=yrange, ytitle=textoidl(ytitle), $
      xtitle=textoidl(oiiihbtitle)
    oplot_iz06, '3868', xx=iz06xx, yy=iz06yy
    calibrate_ratio, [iz06xx], [iz06yy], name='[NeIII]_3869',$
      wave=3868.752D, ncoeff=2, out=out, indx=indx
    oplot_calibration, out[indx]

;; below are some helium lines which I'm ignoring because the
;; line-ratios match pretty well the theoretical case B values
;; tabulated in hydrogen_helium_emissivities.dat    
;    
;; --------------------
;; He I 5876
;    ytitle = 'log (He I \lambda5876 / H\beta)'
;    yrange = [-2,1]
;    djs_plot, [0], [0], /nodata, position=pos, xsty=1, ysty=1, $
;      xrange=oiiihbrange, yrange=yrange, ytitle=textoidl(ytitle), $
;      xtitle=textoidl(oiiihbtitle)
;    oplot_iz06, '5876', xx=iz06xx, yy=iz06yy
;    oplot_g99, 5876, xx=g99xx, yy=g99yy
;    calibrate_ratio, [g99xx,iz06xx], [g99yy,iz06yy], name='HeI_5876',$
;      wave=5876D, ncoeff=1, out=out, indx=indx
;    oplot_calibration, out[indx]
;
;; --------------------
;; He I 4471
;    ytitle = 'log (He I \lambda4471 / H\beta)'
;    yrange = [-2,1]
;    djs_plot, [0], [0], /nodata, position=pos, xsty=1, ysty=1, $
;      xrange=oiiihbrange, yrange=yrange, ytitle=textoidl(ytitle), $
;      xtitle=textoidl(oiiihbtitle)
;    oplot_g99, 4471, xx=g99xx, yy=g99yy
;    calibrate_ratio, [g99xx], [g99yy], name='HeI_4471',$
;      wave=4471D, ncoeff=1, out=out, indx=indx
;    oplot_calibration, out[indx]

    im_plotconfig, /psclose, psfile=psfile, /pdf

; write out    
    struct_print, out
    
    outfile = getenv('IMPRO_DIR')+'/etc/isedfit_forbidden_lineratios.fits'
    im_mwrfits, out, outfile, /clobber

return
end
    
