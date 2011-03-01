function stern_irac_agn, ch1_ab, ch2_ab, ch3_ab, ch4_ab
; jm10mar16ucsd - simple function to determine whether an object
; should be classified as an AGN using the Stern+05 mid-IR color
; criteria; ch[1-4] are the irac AB *magnitudes*

    v2ab = k_vega2ab(filterlist=irac_filterlist(),/kurucz,/silent)
    ch12 = (ch1_ab-ch2_ab)-(v2ab[0]-v2ab[1])
    ch34 = (ch3_ab-ch4_ab)-(v2ab[2]-v2ab[3])

    isagn = (ch34 gt 0.6) and (ch12 gt (0.2*ch34+0.18)) and $
      (ch12 gt (2.5*ch34-3.5))
    
return, isagn
end
