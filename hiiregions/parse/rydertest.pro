pro rydertest

    data = rsex('/home/ioannis/catalogs/hiiregions/1995_ryder.sex')
    w = where(strmatch(data.hii_galaxy,'*3621*'),nhii)

    incl = replicate(67.0,nhii)
    pa = replicate(71.0+90.0,nhii)
;   pa = replicate(71.0,nhii)
    
    ryder = rsex('/home/ioannis/catalogs/hiiregions/95ryder/raw_ryder95.sex')
    w2 = where(strmatch(ryder.hii_galaxy,'*3621*'))

;   niceprint, data[w].hii_galaxy, ryder[w2].hii_galaxy, data[w].hii_region, ryder[w2].hii_region
    
    radius = im_hiiregion_deproject(incl,pa,data[w].raoffset,data[w].deoffset,hii_phi=phi)
    niceprint, data[w].hii_galaxy, ryder[w2].hii_galaxy, data[w].hii_region, ryder[w2].hii_region, $
      radius, ryder[w2].radius
    
stop    
    
return
end
    
