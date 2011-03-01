;invert a 2 x 2 matrix using cramer's rule

function cramer_invert, A

determ = A[0,0] * A[1,1] - A[0,1] * A[1,0]

A_inv = [[A[1,1], -A[0,1]], [-A[1,0], A[0,0]]]

return, A_inv / determ
end
