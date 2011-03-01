function sersicb_func, bb
    common sersicb, nn
return, gamma(2.0*nn)-2D*igamma(2*nn,bb)*gamma(2*nn)
end

function get_sersicb, sersicn
    common sersicb, nn
    nn = sersicn
return, zbrent(0D,20D,func_name='sersicb_func',max_iterations=50)
end
