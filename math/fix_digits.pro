function fix_digits, num, digits
; Martin Downing
    expon = 1.0
    while num/expon gt 1 do expon = expon * 10
    fix_val = num/expon
    fstring = string(digits+2, digits, format = '("(F",I0.0,".",I0.0,")")')
    reads, string(fix_val, format = fstring), fix_val

return, fix_val * expon
end
