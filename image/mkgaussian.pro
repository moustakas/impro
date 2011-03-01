FUNCTION mkgaussian, dim, sigma
;lamucb
   IF NOT keyword_set(dim) THEN dim = 100
   IF NOT keyword_set(sigma) THEN sigma = 20
   x = dist(dim)
   top = -x^2/(2*(sigma)^2)
   return,shift(exp(top), dim/2, dim/2)
END 
