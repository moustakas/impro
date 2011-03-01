; jm02jan21uofa
; returns current working directory
function cwd

    spawn, ['pwd'], d
    d = d[0]+'/'
    
return, d
end
