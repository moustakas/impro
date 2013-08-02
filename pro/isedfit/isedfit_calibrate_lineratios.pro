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

pro isedfit_calibrate_lineratios 

    path = getenv('IM_PROJECTS_DIR')+'/isedfit/'

    allatlas = mrdfits(path+'integrated_atlas_speclinefit_v3.1.fits.gz',1)
    allsdss = read_vagc_garching(/ispecline)
    allhiiinfo = read_hiiregions(linefit=allhii)

    snrcut = 5.0
    
; remove AGN from the galaxy samples and also require S/N>5 on all
; emission lines of interest 
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
    splog, 'SDSS ', nsdss
    sdss = allsdss[keep]
    
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
    atlas = allatlas[keep]

    keep = where($
      allhii.oiii_5007[0]/allhii.oiii_5007[1] gt snrcut and $
      allhii.oii_3727[0]/allhii.oii_3727[1] gt snrcut and $
      allhii.h_beta[0]/allhii.h_beta[1] gt snrcut and $
      allhii.h_alpha[0]/allhii.h_alpha[1] gt snrcut and $
      allhii.nii_6584[0]/allhii.nii_6584[1] gt snrcut and $
      allhii.sii_6716[0]/allhii.sii_6716[1] gt snrcut and $
      allhii.sii_6731[0]/allhii.sii_6731[1] gt snrcut,nhii)
    splog, 'HII ', nhii
    hii = allhii[keep]
    hiiinfo = allhiiinfo[keep]
    
; build the line-ratios
    lineratio, sdss, 'OIII_5007', 'H_BETA', 'OII_3727', $
      'H_BETA', soiiihb, soiiihb_err, soiihb, soiihb_err, snrcut=snrcut
    lineratio, atlas, 'OIII_5007', 'H_BETA', 'OII_3727', $
      'H_BETA', aoiiihb, aoiiihb_err, aoiihb, aoiihb_err, snrcut=snrcut
    lineratio, hii, 'OIII_5007', 'H_BETA', 'OII_3727', $
      'H_BETA', hoiiihb, hoiiihb_err, hoiihb, hoiihb_err, snrcut=snrcut

    lineratio, sdss, 'OIII_5007', 'H_BETA', 'NII_6584', $
      'H_ALPHA', soiiihb, soiiihb_err, sniiha, sniiha_err, snrcut=snrcut
    lineratio, atlas, 'OIII_5007', 'H_BETA', 'NII_6584', $
      'H_ALPHA', aoiiihb, aoiiihb_err, aniiha, aniiha_err, snrcut=snrcut
    lineratio, hii, 'OIII_5007', 'H_BETA', 'NII_6584', $
      'H_ALPHA', hoiiihb, hoiiihb_err, hniiha, hniiha_err, snrcut=snrcut

    lineratio, sdss, 'OIII_5007', 'H_BETA', ['SII_6716','SII_6731'], $
      'H_ALPHA', soiiihb, soiiihb_err, ssiiha, ssiiha_err, snrcut=snrcut
    lineratio, atlas, 'OIII_5007', 'H_BETA', ['SII_6716','SII_6731'], $
      'H_ALPHA', aoiiihb, aoiiihb_err, asiiha, asiiha_err, snrcut=snrcut
    lineratio, hii, 'OIII_5007', 'H_BETA', ['SII_6716','SII_6731'], $
      'H_ALPHA', hoiiihb, hoiiihb_err, hsiiha, hsiiha_err, snrcut=snrcut

; fit the median relations, giving equal weight to every sample
    xx = [soiiihb,aoiiihb,hoiiihb]
    weight = [soiiihb*0.0+1.0/nsdss,aoiiihb*0.0+1.0/natlas,hoiiihb*0.0+1.0/nhii]

    binsz = 0.3
    oiihbmed = im_medxbin(xx,[soiihb,aoiihb,hoiihb],binsz,weight=weight)
    niihamed = im_medxbin(xx,[sniiha,aniiha,hniiha],binsz,weight=weight)
    siihamed = im_medxbin(xx,[ssiiha,asiiha,hsiiha],binsz,weight=weight)

    ncoeff = 4
    oiihbcoeff = poly_fit(oiihbmed.medx,oiihbmed.medy,ncoeff-1)
    niihacoeff = poly_fit(niihamed.medx,niihamed.medy,ncoeff-1)
    siihacoeff = poly_fit(siihamed.medx,siihamed.medy,ncoeff-1)

    splog, '[OII]/Hb coeff ', oiihbcoeff
    splog, '[NII]/Ha coeff ', niihacoeff
    splog, '[SII]/Ha coeff ', siihacoeff
    
; pack into a structure and write out; these coefficients allow one to
; produce the [OII]/Hbeta, [NII]/Ha, and [SII]/Ha line-ratios given
; the [OIII]/Hbeta line-ratio
    out = {$
      oiihb_coeff: reform(oiihbcoeff),$
      niiha_coeff: reform(niihacoeff),$
      siiha_coeff: reform(siihacoeff),$
      oiihb_scatter: 0.1,$ ; dex
      niiha_scatter: 0.1,$
      siiha_scatter: 0.1}
    outfile = getenv('IMPRO_DIR')+'/etc/isedfit_lineratios_v1.0.fits'
    im_mwrfits, out, outfile, /clobber

stop    
; build a QAplot
    nogrey = 0
    oiiihbrange = [-1.4,1.0]
    oiiihbaxis = range(oiiihbrange[0]+0.1,oiiihbrange[1]-0.1,200)
    
    psfile = path+'qa_calibrate_lineratios.ps'
    im_plotconfig, 4, pos, psfile=psfile, xmargin=[1.3,0.4], $
      width=6.7, charsize=1.5
    
; [OII]/Hb vs [OIII]/Hb
    im_hogg_scatterplot, soiiihb, soiihb, position=pos[*,0], nogrey=nogrey, $
      xsty=1, ysty=1, xrange=oiiihbrange, yrange=[-1,1], $
      ytitle=textoidl('log ([OII] \lambda3727 / H\beta)'), $
      xtickname=replicate(' ',10)
    djs_oplot, aoiiihb, aoiihb, psym=8, color='blue'
    djs_oplot, hoiiihb, hoiihb, psym=6, color='orange', symsize=0.5
;   oploterror, oiihbmed.medx, oiihbmed.medy, oiihbmed.sigy, $
;     psym=symcat(9,thick=8), symsize=3.0, errthick=6.0
    djs_oplot, oiihbmed.medx, oiihbmed.medy, line=0, thick=10
    djs_oplot, oiihbmed.medx, oiihbmed.medy+0.2, line=5, thick=10
    djs_oplot, oiihbmed.medx, oiihbmed.medy-0.2, line=5, thick=10
;   djs_oplot, oiihbmed.medx, oiihbmed.quant75, line=5, thick=10
;   djs_oplot, oiihbmed.medx, oiihbmed.quant25, line=5, thick=10
;   djs_oplot, oiiihbaxis, poly(oiiihbaxis,oiihbcoeff), color='green'
    
; [NII]/Ha vs [OIII]/Hb
    im_hogg_scatterplot, soiiihb, sniiha, position=pos[*,1], /noerase, nogrey=nogrey, $
      xsty=1, ysty=1, xrange=oiiihbrange, yrange=[-2.2,0.1], $
      ytitle=textoidl('log ([N II] \lambda6584 / H\alpha)'), $
      xtickname=replicate(' ',10)
    djs_oplot, aoiiihb, aniiha, psym=8, color='blue'
    djs_oplot, hoiiihb, hniiha, psym=6, color='orange', symsize=0.5
;   oploterror, niihamed.medx, niihamed.medy, niihamed.sigy, $
;     psym=symcat(9,thick=8), symsize=3.0, errthick=6.0
    djs_oplot, niihamed.medx, niihamed.medy, line=0, thick=10
    djs_oplot, niihamed.medx, niihamed.medy+0.2, line=5, thick=10
    djs_oplot, niihamed.medx, niihamed.medy-0.2, line=5, thick=10
;   djs_oplot, oiiihbaxis, poly(oiiihbaxis,niihacoeff), color='green'

; [SII]/Ha vs [OIII]/Hb
    im_hogg_scatterplot, soiiihb, ssiiha, position=pos[*,2], /noerase, nogrey=nogrey, $
      xsty=1, ysty=1, xrange=oiiihbrange, yrange=[-2.0,0.2], $
      xtitle=textoidl('log ([O III] \lambda5007 / H\beta)'), $
      ytitle=textoidl('log ([S II] \lambda\lambda6716,31 / H\alpha)')
    djs_oplot, aoiiihb, asiiha, psym=8, color='blue'
    djs_oplot, hoiiihb, hsiiha, psym=6, color='orange', symsize=0.5
;   oploterror, siihamed.medx, siihamed.medy, siihamed.sigy, $
;     psym=symcat(9,thick=8), symsize=3.0, errthick=6.0
    djs_oplot, siihamed.medx, siihamed.medy, line=0, thick=10
    djs_oplot, siihamed.medx, siihamed.medy+0.2, line=5, thick=10
    djs_oplot, siihamed.medx, siihamed.medy-0.2, line=5, thick=10
;   djs_oplot, oiiihbaxis, poly(oiiihbaxis,siihacoeff), color='green'

;   im_legend, [
    
    im_plotconfig, /psclose, psfile=psfile, /pdf
    
; ---------------------------------------------------------------------------
; now show the     
    
    
    
stop    
    
return
end
    
