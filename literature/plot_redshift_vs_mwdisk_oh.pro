pro plot_redshift_vs_mwdisk_oh, postscript=postscript, pdf=pdf
; jm08apr09nyu - oxygen abundance vs redshift for Milky Way thin disk
;                stars from Reddy et al. (2003)

; default plotting variables - note that POSTSCRIPT over-rides PDF

    postthick1 = 2.0
    postthick2 = 2.0
    charsize = 2.0
    textcolor1 = 'white'
    pspath = getenv('PAPERSPATH')+'/literature/'
    datapath = pspath+'data/'

    if keyword_set(pdf) then begin
       pspath = pspath+'keynote/' ; for keynote presentations
       textcolor1 = 'white'
       charsize = 2.2
       postthick1 = 8.0
       postthick2 = 6.0
    endif
    if keyword_set(postscript) then begin
       textcolor1 = 'black'
       postthick1 = 4.0
       postthick2 = 3.0
    endif

    ohsun = 8.7

    xsize = 8.5 & ysize = 7.5
    pagemaker, nx=1, ny=1, xspace=0, yspace=0.0, width=7.0, height=5.5, $
      xmargin=[1.2,0.3], ymargin=[1.0,1.0], xpage=xsize, ypage=ysize, $
      position=pos, /normal

    psname = 'redshift_vs_mwdisk_oh'
    im_openclose, pspath+psname, postscript=(keyword_set(postscript) or $
      keyword_set(pdf)), xsize=xsize, ysize=ysize, silent=keyword_set(pdf)

    data03 = mrdfits(datapath+'03reddy.fits',1,silent=0)

    data06_table2 = rsex(datapath+'06reddy_table2.dat')
    data06_table45 = rsex(datapath+'06reddy_table45.dat')
    data06 = struct_addtags(struct_trimtags(data06_table2,select=['AGE','FEH']),$
      struct_trimtags(data06_table45,select=['O']))

    xtitle = 'Formation Redshift'
    xtitle2 = 'Stellar Age (Gyr)'
    ytitle = '12 + log (O/H)'
;   ytitle = textoidl('\Delta'+'log(O/H)')

    xlog = 0
    
    xrange = [-0.02,0.92]
;   xrange = [0.05,1]
;   yrange = [-0.8,0.8]
    yrange = [8.2,9.2]

    zaxis = findgen(((1.0)-(-0.05))/0.01+1)*0.01+(-0.05)
    timelabel1 = [1.0,3.0,5.0,7.0] ; [Gyr]
    
    plot, [0], [0], /nodata, xthick=postthick1, ythick=postthick1, charsize=charsize, xlog=xlog, $
      charthick=postthick2, thick=postthick1, xtitle=xtitle, ytitle=ytitle, color=fsc_color(textcolor1,100), $
      ymargin=[4,3], xstyle=11, ystyle=3, xrange=xrange, yrange=yrange, position=pos[*,0]
;   djs_oplot, !x.crange, ohsun*[1,1], line=0, thick=postthick1, color=fsc_color(textcolor1,100)
    axis, /xaxis, xthick=postthick1, xsty=1, xrange=xrange, charsize=charsize, $
      charthick=postthick2, xtitle='Stellar Age (Gyr)', color=fsc_color(textcolor1,100), $
      xtickv=getredshift(getage(0.0)-timelabel1), xticks=n_elements(timelabel1)-1L, $
      xtickname=string(timelabel1,format='(I0)')
;   im_xyouts_title, xtitle=xtitle2, charsize=charsize, charthick=postthick2, color=fsc_color(textcolor1,100), $
;     yspacing=3.5

; 2003 data; compute and plot the median values

    indx = lindgen(n_elements(data03))
    x2 = 10.0^data03[indx].logt9      ; age [Gyr]
    x = getredshift(getage(0.0)-x2) ; redshift [GETREDSHIFT requires the lookback time, Gyr]
;   y = data03[indx]._o_fe_
;   y = data03[indx]._fe_h_ + data03[indx]._o_fe_
    y = data03[indx]._fe_h_ + data03[indx]._o_fe_ + ohsun

    medcolor = 'forest green' ; 'dodger blue'
    ccolor = 'orange' ; 
    im_symbols, 108, psize=1.2, /fill, color=fsc_color(ccolor,120)
    djs_oplot, x, y, ps=8, color=fsc_color(ccolor,120)

    medbin = 0.1
    minx = 0.05
    minpts = 5
    
    running = im_medxbin(x,y,medbin,minx=minx,minpts=minpts,/verbose)
    im_symbols, 106, psize=3.0, thick=postthick1, color=fsc_color(medcolor,130)
    oploterror, running.binctr, running.meany, running.binsz/2.0+fltarr(running.nbins), $
      running.sigy, ps=8, errthick=postthick1, color=fsc_color(medcolor,130), $
      errcolor=fsc_color(medcolor,130)
;   oploterror, running.binctr, running.meany, running.binsz/2.0+fltarr(running.nbins), $
;     running.sigy95-running.meany, ps=8, errthick=(postthick1)>2L, /hibar, $
;     color=fsc_color(medcolor), errcolor=fsc_color(medcolor)
;   oploterror, running.binctr, running.meany, running.binsz/2.0+fltarr(running.nbins), $
;     running.meany-running.sigy05, ps=8, errthick=(postthick1)>2L, /lobar, $
;     color=fsc_color(medcolor), errcolor=fsc_color(medcolor)
;   oploterror, running.binctr, running.medy, running.binsz/2.0+fltarr(running.nbins), $
;     running.sigy75-running.medy, ps=8, errthick=(postthick1)>2L, /hibar, $
;     color=fsc_color(medcolor), errcolor=fsc_color(medcolor)
;   oploterror, running.binctr, running.medy, running.binsz/2.0+fltarr(running.nbins), $
;     running.medy-running.sigy25, ps=8, errthick=(postthick1)>2L, /lobar, $
;     color=fsc_color(medcolor), errcolor=fsc_color(medcolor)

; fit a line to the mean points

    coeff = linfit(running.binctr-0.5,running.meany,measure_errors=running.sigy)
;   coeff = linfit(running.binctr-0.5,running.meany)
    oplot, zaxis, poly(zaxis-0.5,coeff), line=5, thick=postthick1, color=fsc_color(textcolor1,100)
    coeff[0] = coeff[0] + 0.5
    splog, '12 + log (O/H) = '+strtrim(string(coeff[0],format='(F12.2)'),2)+' + '+$
      strtrim(string(coeff[1],format='(F12.2)'),2)+' z'
    
; 2006 data - thick disk stars

;   indx = where((data06.age gt -900.0) and (data06.feh gt -900.0) and (data06.o gt -900.0))
;
;   x2 = data06[indx].age       ; [Gyr]
;   x = getredshift(getage(0.0)-x2)
;   y = data06[indx].feh + data06[indx].o

;   im_symbols, 106, psize=1.5, fill=0, color=fsc_color('blue',!d.table_size-74)
;   djs_oplot, x, y, ps=8, color=fsc_color('blue',!d.table_size-74)

; legend    
    
;   legend, 'Reddy et al. (2003)', /right, /bottom, charsize=1.5, charthick=postthick2, $
;     box=0, textcolor=fsc_color(textcolor1,100)

    im_openclose, postscript=(keyword_set(postscript) or keyword_set(pdf)), /close    
    if keyword_set(pdf) then begin
       splog, 'Writing '+pspath+psname+'.pdf'
       spawn, 'ps2pdf13 '+pspath+psname+'.ps '+pspath+psname+'.pdf', /sh
       rmfile, pspath+psname+'.ps' ; delete the junk ps file
    endif

stop    
    
return
end
