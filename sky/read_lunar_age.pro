function read_lunar_age, age=age
; jm04mar04uofa
; compute the surface brightness of the Moon [mag/arcsec2] as a
; function of lunar age

; http://www.astro.utoronto.ca/~patton/astro/mags.html

    age = [0,3,7,10,14]    ; lunar age [days]
    nage = n_elements(age)
    
    lunar = {$
      U: 0.0, $
      B: 0.0, $
      V: 0.0, $
      R: 0.0, $
      I: 0.0}
    lunar = replicate(lunar,nage)

    lunar.U = [22.0,21.5,19.9,18.5,17.0]
    lunar.B = [22.7,22.4,21.6,20.7,19.5]
    lunar.V = [21.8,21.7,21.4,20.7,20.0]
    lunar.R = [20.9,20.8,20.6,20.3,19.9]
    lunar.I = [19.9,19.9,19.7,19.5,19.2]

return, lunar
end
    
