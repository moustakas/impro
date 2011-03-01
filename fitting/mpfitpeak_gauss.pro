; Gaussian Function
function mpfitpeak_gauss, x, p, _extra=extra
  sz = size(x)
  if sz(sz(0)+1) EQ 5 then smax = 26D else smax = 13.
  u = mpfitpeak_u(x, p)
  mask = u LT (smax^2)  ;; Prevents floating underflow
  if n_elements(p) GE 4 then f = p(3) else f = 0
  if n_elements(p) GE 5 then f = f + p(4)*x
  return,  f + p(0) * mask * exp(-0.5 * temporary(u) * mask)
end
