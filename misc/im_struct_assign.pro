function im_struct_assign, in, out, nozero=nozero
; jm09mar01nyu - simple wrapper on struct_assign to get around the
; pass-by-reference pain-in-the-neck
    struct_assign, in, out, nozero=nozero
return, out
end
