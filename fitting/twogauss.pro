function twogauss, xval, pp
; jm02jul6uofa
    
    term1 = exp( - (xval - pp[1])^2 / (2. * pp[2]^2) )
    term2 = exp( - (xval - pp[4])^2 / (2. * pp[5]^2) )

    yval = pp[0] * term1 + pp[3] * term2
    
return, yval
end
