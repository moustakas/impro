stop

; SEE TALK_08APR_UCSC - OBSOLETE!!

pro plot_06erb, postscript=postscript, pdf=pdf
; jm08apr13nyu - plots from Erb et al. (2006); note that the Chabrier
;                IMF is assumed through-out!

; default plotting variables - note that POSTSCRIPT over-rides PDF

    postthick1 = 2.0
    postthick2 = 2.0
    postthick3 = 2.0
    charsize = 2.0
    textcolor1 = 'white'
    pspath = getenv('PAPERSPATH')+'/literature/'

    datapath = pspath+'data/'
    erbcolor = 'navy'
    liucolor_loz = 'orange'
    liucolor_hiz = 'red'

    if keyword_set(pdf) then begin
       pspath = pspath+'keynote/' ; for keynote presentations
       textcolor1 = 'white'
       charsize = 2.2
       postthick1 = 8.0
       postthick2 = 6.0
       postthick3 = 10.0
       erbcolor = 'forest green'
       liucolor_loz = 'orange'
       liucolor_hiz = 'red'
    endif
    if keyword_set(postscript) then begin
       textcolor1 = 'black'
       postthick1 = 4.0
       postthick2 = 3.0
       postthick3 = 6.0
    endif

; read the data    
    
    sohdust = read_sdss_mz_sample(/mzhii_log12oh)
    salldust = read_sdss_main(/ispec)
    skcorr = read_sdss_mz_sample(/mzhii_kcorr)
    
; the first mass bin is an upper limit    
    
    erb_all = rsex(datapath+'06erb.sex')
    ul = where(erb_all.fnii lt 0.0,comp=good)
    erb_ul = erb_all[ul]
    erb = erb_all[good]
    
    liu = rsex(datapath+'08liu_table2.sex')
    liu_n2 = rsex(datapath+'08liu_table3.sex')
    liu_o3n2 = rsex(datapath+'08liu_table4.sex')
    
; ---------------------------------------------------------------------------    
; BPT diagram
; ---------------------------------------------------------------------------    
    
    xsize = 8.5 & ysize = 7.2
    pagemaker, nx=1, ny=1, xspace=0.0, yspace=0.0, width=6.7, height=5.7, $
      xmargin=[1.5,0.3], ymargin=[0.4,1.1], xpage=xsize, ypage=ysize, $
      position=pos, /normal

    psname = 'bpt_06erb' ; bad name!!!
    im_openclose, pspath+psname, postscript=(keyword_set(postscript) or $
      keyword_set(pdf)), xsize=xsize, ysize=ysize, silent=keyword_set(pdf)

    xtitle = 'log (N II]/H\alpha)'
    ytitle = 'log (O III]/H\beta)'
    
    xrange = [-1.6,0.7]
    yrange = [-1.0,1.2]

; plot the sdss galaxies

    sdssindx = where((salldust.nii_6584[0]/salldust.nii_6584[1] gt 5.0) and $
      (salldust.oiii_5007[0]/salldust.oiii_5007[1] gt 5.0) and $
      (salldust.h_alpha[0]/salldust.h_alpha[1] gt 5.0) and $
      (salldust.h_beta[0]/salldust.h_beta[1] gt 5.0))
    
    sniiha = alog10(salldust[sdssindx].nii_6584[0]/salldust[sdssindx].h_alpha[0])
    soiiihb = alog10(salldust[sdssindx].oiii_5007[0]/salldust[sdssindx].h_beta[0])
    
    im_symbols, 108, psize=0.1, qfill=1, color=fsc_color(textcolor1,100)
    im_hogg_scatterplot, sniiha, soiiihb, outliers=1, xsty=1, ysty=1, $
      label=0, outpsym=8, outsymsize=1.0, outcolor=fsc_color(textcolor1,100), $
      levels=errorf(0.5*[1.0,2.0,3.0]), $
      xrange=xrange, yrange=yrange, xtitle=textoidl(xtitle), ytitle=textoidl(ytitle), $
      charsize=charsize, charthick=postthick2, xthick=postthick1, ythick=postthick1, position=pos[*,0], $
      color=fsc_color(textcolor1,100), axis_color=fsc_color(textcolor1,100), $
      contour_color=fsc_color(textcolor1,100), /nogrey, c_linestyle=0, cthick=postthick2

; now plot Liu et al. (2008); need to treat the upper limit
; correctly; need to also include Shapley et al. 

; z~1.0    
    
    loz = where((liu.z lt 1.1) and (liu.foiii gt 0.0))

    niiha_loz = liu[loz].fnii/liu[loz].fha
    niiha_loz_err = im_compute_error(liu[loz].fnii,liu[loz].fnii_err,$
      liu[loz].fha,liu[loz].fha_err,/quotient)

    oiiihb_loz = liu[loz].foiii/liu[loz].fhb
    oiiihb_loz_err = im_compute_error(liu[loz].foiii,liu[loz].foiii_err,$
      liu[loz].fhb,liu[loz].fhb_err,/quotient)
    
    im_symbols, 108, psize=1.5, thick=postthick3, fill=1, color=fsc_color(liucolor_loz,150)
    oploterror, alog10(niiha_loz), alog10(oiiihb_loz), niiha_loz_err/niiha_loz/alog(10.0), $
      oiiihb_loz_err/oiiihb_loz/alog(10.0), $
      thick=postthick3, errthick=postthick3, psym=8, color=fsc_color(liucolor_loz,150), $
      errcolor=fsc_color(liucolor_loz,150)

; z~1.4

    hiz = where((liu.z gt 1.1) and (liu.foiii gt 0.0))

    niiha_hiz = liu[hiz].fnii/liu[hiz].fha
    niiha_hiz_err = im_compute_error(liu[hiz].fnii,liu[hiz].fnii_err,$
      liu[hiz].fha,liu[hiz].fha_err,/quotient)

    oiiihb_hiz = liu[hiz].foiii/liu[hiz].fhb
    oiiihb_hiz_err = im_compute_error(liu[hiz].foiii,liu[hiz].foiii_err,$
      liu[hiz].fhb,liu[hiz].fhb_err,/quotient)
    
    im_symbols, 108, psize=1.5, thick=postthick3, fill=1, color=fsc_color(liucolor_hiz,150)
    oploterror, alog10(niiha_hiz), alog10(oiiihb_hiz), niiha_hiz_err/niiha_hiz/alog(10.0), $
      oiiihb_hiz_err/oiiihb_hiz/alog(10.0), $
      thick=postthick3, errthick=postthick3, psym=8, color=fsc_color(liucolor_hiz,150), $
      errcolor=fsc_color(liucolor_hiz,150)
    
    im_openclose, postscript=(keyword_set(postscript) or keyword_set(pdf)), /close    
    if keyword_set(pdf) then begin
       splog, 'Writing '+pspath+psname+'.pdf'
       spawn, 'ps2pdf13 '+pspath+psname+'.ps '+pspath+psname+'.pdf', /sh
       rmfile, pspath+psname+'.ps' ; delete the junk ps file
    endif

; ---------------------------------------------------------------------------    
; mass-metallicity relation    
; ---------------------------------------------------------------------------    
    
    xsize = 8.5 & ysize = 7.0
    pagemaker, nx=1, ny=1, xspace=0.0, yspace=0.0, width=7.0, height=5.5, $
      xmargin=[1.2,0.3], ymargin=[0.4,1.1], xpage=xsize, ypage=ysize, $
      position=pos, /normal

    psname = 'mz_06erb'
    im_openclose, pspath+psname, postscript=(keyword_set(postscript) or $
      keyword_set(pdf)), xsize=xsize, ysize=ysize, silent=keyword_set(pdf)

    xtitle = 'log (M_{*}/M'+sunsymbol()+')'
    ytitle = '12 + log (O/H)_{N2}'
    
    xrange = [8.6,11.7]
    yrange = [8.1,8.85]

; plot the sdss galaxies

    sindx = where((sohdust.zstrong_niiha gt -900))

    smass = skcorr[sindx].mass ; Chabrier IMF!
    soh = 8.90+0.57*sohdust[sindx].zstrong_niiha
    
    im_hogg_scatterplot, smass, soh, outliers=plot_outliers, xsty=1, ysty=1, $
      label=0, outpsym=8, outsymsize=1.0, outcolor='default', levels=errorf(0.5*[1.0,2.0,3.0]), $
      xrange=xrange, yrange=yrange, xtitle=textoidl(xtitle), ytitle=textoidl(ytitle), $
      charsize=charsize, charthick=postthick2, xthick=postthick1, ythick=postthick1, position=pos[*,0], $
      color=fsc_color(textcolor1,100), axis_color=fsc_color(textcolor1,100), $
      contour_color=fsc_color(textcolor1,100), /nogrey, c_linestyle=2, cthick=postthick3

;   plot, [0], [0], /nodata, xthick=postthick1, ythick=postthick1, $
;     xtitle=textoidl(xtitle), ytitle=textoidl(ytitle), $
;     charsize=charsize, charthick=postthick2, xsty=1, ysty=1, xrange=xrange, $
;      position=pos[*,0], yrange=yrange, color=fsc_color(textcolor1,150)

    im_symbols, 106, psize=1.5, thick=postthick3, fill=1, color=fsc_color(erbcolor,150)
    oploterror, alog10(1D10*erb.mass), erb.log12oh, erb.mass_err/erb.mass/alog(10.0), erb.log12oh_err, $
      thick=postthick3, errthick=postthick3, psym=8, color=fsc_color(erbcolor,150), $
      errcolor=fsc_color(erbcolor,150)
;   im_symbols, 111, psize=1.5, thick=postthick3, fill=1, color=fsc_color(erbcolor,150)
    oploterror, alog10(1D10*erb_ul.mass), -erb_ul.log12oh, erb_ul.mass_err/erb_ul.mass/alog(10.0), $
      erb_ul.log12oh_err, thick=postthick3, errthick=postthick3, psym=8, color=fsc_color(erbcolor,150), $
      errcolor=fsc_color(erbcolor,150)
    arrow, alog10(1D10*erb_ul.mass), -erb_ul.log12oh, alog10(1D10*erb_ul.mass), $
      -erb_ul.log12oh-0.05, /data, hsize=-0.6, color=fsc_color(erbcolor,150), $
      hthick=postthick3, thick=postthick3

; now plot Liu et al. (2008); use the N2 calibrator and differentiate
; between galaxies at z~1 and z~1.4

    loz = where(liu_n2.z lt 1.1,comp=hiz)
    
    im_symbols, 108, psize=1.5, thick=postthick3, fill=1, color=fsc_color(liucolor_loz,150)
    oploterror, liu_n2[loz].mass, liu_n2[loz].log12oh_n2, liu_n2[loz].mass_err, liu_n2[loz].log12oh_n2_err, $
      thick=postthick3, errthick=postthick3, psym=8, color=fsc_color(liucolor_loz,150), $
      errcolor=fsc_color(liucolor_loz,150)
    im_symbols, 108, psize=1.5, thick=postthick3, fill=1, color=fsc_color(liucolor_hiz,150)
    oploterror, liu_n2[hiz].mass, liu_n2[hiz].log12oh_n2, liu_n2[hiz].mass_err, liu_n2[hiz].log12oh_n2_err, $
      thick=postthick3, errthick=postthick3, psym=8, color=fsc_color(liucolor_hiz,150), $
      errcolor=fsc_color(liucolor_hiz,150)

;   im_symbols, 108, psize=1.5, thick=postthick3, fill=1, color=fsc_color(liucolor_o3n2,150)
;   oploterror, liu_o3n2.mass, liu_o3n2.log12oh_o3n2, liu_o3n2.mass_err, liu_o3n2.log12oh_o3n2_err, $
;     thick=postthick3, errthick=postthick3, psym=8, color=fsc_color(liucolor_o3n2,150), $
;     errcolor=fsc_color(liucolor_o3n2,150)

    im_openclose, postscript=(keyword_set(postscript) or keyword_set(pdf)), /close    
    if keyword_set(pdf) then begin
       splog, 'Writing '+pspath+psname+'.pdf'
       spawn, 'ps2pdf13 '+pspath+psname+'.ps '+pspath+psname+'.pdf', /sh
       rmfile, pspath+psname+'.ps' ; delete the junk ps file
    endif

stop    
    
return
end
