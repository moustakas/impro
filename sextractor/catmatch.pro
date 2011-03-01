pro rcrit, x1, y1, x2, y2
; routine to help decide on a critical matching radius.  basically
; generates a curve of growth

    rad = findgen(8)+1.0
    nrad = n_elements(rad)

    nmatch = lonarr(nrad)
    
    for k = 0L, nrad-1L do begin
    
       srcor, x1, y1, x2, y2, rad[k], nm1, nm2, /option
       nmatch[k] = n_elements(nm1)
       
    endfor

    plotsym, 0, 1, /fill
    plot, rad, nmatch, ps=8, xsty=3, ysty=3, xtit='Radius (pixels)', $
      ytit='Number of Sources'

retall
end

pro catmatch, bfinal, rfinal, ifinal, check=check, crit=crit, plotcheck=plotcheck
; jm01may21uofa
; match sources based on (x,y) pixel positions (astrometry is good)

    colortable2

    if not keyword_set(crit) then crit = 5.0 ; critical matching radius (pixels)
    
; read in the three catalogs

    fnames = ['Bw_catalog.txt','R_catalog.txt','I_catalog.txt']

    readfast, fnames[0], bdata, nlines=nb & bx = bdata[1,*] & by = bdata[2,*]
    readfast, fnames[1], rdata, nlines=ni & rx = rdata[1,*] & ry = rdata[2,*]
    readfast, fnames[2], idata, nlines=nr & ix = idata[1,*] & iy = idata[2,*]

    if keyword_set(check) then rcrit, bx, by, ix, iy
    
    srcor, bx, by, rx, ry, crit, brmatch, rbmatch, /option ; match B to R
    srcor, bx, by, ix, iy, crit, bimatch, ibmatch, /option ; match B to I

    srcor, bx[brmatch], by[brmatch], bx[bimatch], by[bimatch], crit, bgood, dum, /option

    srcor, bx, by, bx[brmatch[bgood]], by[brmatch[bgood]], crit, bcat, dum, /option ; final B-catalog
    srcor, rx, ry, bx[brmatch[bgood]], by[brmatch[bgood]], crit, rcat, dum, /option ; final R-catalog
    srcor, ix, iy, bx[brmatch[bgood]], by[brmatch[bgood]], crit, icat, dum, /option ; final I-catalog

    bfinal = bdata[*,bcat]
    rfinal = rdata[*,rcat]
    ifinal = idata[*,icat]    

    if keyword_set(plotcheck) then begin
    
; overplot the final coordinates

       window, 0, xs=450, ys=450
       plotsym, 0, 0.8
       plot, bx[bcat], by[bcat], ps=8, color=3, xsty=3, ysty=3, xtit='x', ytit='y'
       plotsym, 0, 0.8, /fill
       oplot, rx[rcat], ry[rcat], ps=8, color=4
       oplot, ix[icat], iy[icat], ps=2, color=5
       
    endif
       
return
end
