pro test_bisector

    seed = -1L
    
    nn = 20L
    xx = randomu(seed,nn)
    yy = xx + 0.3*randomu(seed,nn)
    ww = xx*0.0+1.0
    ww[0] = 5.0

    yydup = [yy,replicate(yy[0],4)]
    xxdup = [xx,replicate(xx[0],4)]
    
    sixlin, xx, yy, a, siga, b, sigb
    sixlin, xxdup, yydup, adup, sigadup, bdup, sigbdup
    sixlin_weighted, xx, yy, aw, sigaw, bw, sigbw, weight=ww

;   print, a[2], b[2], siga[2], sigb[2]
;   print, adup[2], bdup[2], sigadup[2], sigbdup[2]
;   print, aw[2], bw[2], sigaw[2], sigbw[2]

    print, adup-aw, bdup-bw, sigadup-sigaw, sigbdup-sigbw

    plot, xx, yy, ps=4, xsty=3, ysty=3, sym=2.0
    plots, xx[0], yy[0], ps=7, thick=3, sym=3.0
    oplot, !x.crange, poly(!x.crange,[a[2],b[2]]), line=0, thick=3
    oplot, !x.crange, poly(!x.crange,[adup[2],bdup[2]]), line=2, thick=6
    oplot, !x.crange, poly(!x.crange,[aw[2],bw[2]]), line=0, thick=1

return
end
    
