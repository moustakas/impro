pro jhist, bin=bin, half=half, _extra=_extra

    if n_elements(bin)  eq 0 then bin  = 0.2
    if n_elements(half) eq 0 then half = 1

; generate a couple shifted histograms
    nx=1000
    a=randomn(seed,nx)
    b=randomn(seed,nx)+2.0

; get the binning
    plothist,a,ax,ay,bin=bin,half=half,/noplot
    plothist,b,bx,by,bin=bin,half=half,/noplot
    na = n_elements(ax)
    nb = n_elements(bx)

; show the plots
    plot, [ax[0],ax,ax[na-1]],[0,ay,0],$
          psym=10,xr=[min(ax),max(bx)],yr=[0,max([ay,by])]
    oplot,[bx[0],bx,bx[nb-1]],[0,by,0],$
          psym=10,thick=2

; calculate the overlap histogram
; i haven't done any sanity checks for half=0 ....
    xar=ax
    yar=fltarr(na)
    for i=0,na-1 do begin
       ix = where(bx gt ax[i]-bin/2 and bx lt ax[i]+bin/2)
       if ix[0] ne -1 then yar[i]=min([ay[i],by[ix]])
    endfor

; now remove the leading zeros
    zix=where(yar eq 0.0)
    shx=zix-shift(zix,-1)
    thx=where(shx ne -1)
    xar=xar[zix[thx[0]]+1:na-1]
    yar=yar[zix[thx[0]]+1:na-1]

; do a sanity-check overplot...
;   oplot,xar,yar,psym=4,symsize=2

; shamelessly crip from plothist to work out the polyfill region
    Xfill = transpose([[Xar-bin/2.0],[Xar+bin/2.0]])
    Xfill = reform(Xfill, n_elements(Xfill))
    Xfill = [Xfill[0], Xfill, Xfill[n_elements(Xfill)-1]]
    Yfill = transpose([[Yar],[Yar]])
    Yfill = reform(Yfill, n_elements(Yfill))
    Yfill = [0, Yfill, 0]
    Xfill = Xfill > !X.CRANGE[0] < !X.CRANGE[1]
    Yfill = Yfill > !Y.CRANGE[0] < !Y.CRANGE[1]

; done deal.  change polyfill any way desired...
    polyfill,xfill,yfill,_extra=_extra

 end

