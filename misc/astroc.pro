FUNCTION astroc,help=help
   as = {name:    '', $
         pi:      dblarr(1), $
         EE:      dblarr(1), $
         Msun:    dblarr(1), $
         Lsun:    dblarr(1), $
         Rsun:    dblarr(1), $
         c:       dblarr(1), $
         G:       dblarr(1), $
         e:       dblarr(1), $
         h:       dblarr(1), $
         me:      dblarr(1), $
         mp:      dblarr(1), $
         alpha:   dblarr(1), $
         sigmaT:  dblarr(1), $
         k:       dblarr(1), $
         NA:      dblarr(1), $
         sigma:   dblarr(1), $
         keV:     dblarr(1), $
         omega:   dblarr(1), $
         pc:      dblarr(1)}
   
   as.pi    = 3.14159265358979  ;
   as.EE    = 2.718281828459045 ; exp(1)
   as.Msun  = 1.9891e+33        ; g
   as.Lsun  = 3.826e+33         ; erg/s
   as.Rsun  = 6.9599e+10        ; cm
   as.c     = 2.9979245620e+10  ; cm/s ; speed of light
   as.G     = 6.67259e-08       ; cgs ; gravitational constant
   as.e     = 4.8032068e-10     ; charge of electron
   as.h     = 6.6260755e-27     ; erg*s ; Planck's constant
   as.me    = 9.1093897e-28     ; g ; mass of electron
   as.mp    = 1.6726231e-24     ; g ; mass of proton
   as.alpha = 1/(as.h*as.c/(2*as.pi*as.e^2)) ; fine structure constant
   as.sigmaT= 6.6524616e-25     ; cm^2 ; Thomson cross section=8pi/3*re^2
   as.k     = 1.380658e-16      ; erg/K ; Boltzman constant
   as.NA    = 6.0221367e+23     ; mol^-1 ; Avogadro constant
   as.sigma = 5.67051e-5        ; Stefan-Boltzmann constant
   as.keV = 1.602192e-9         ; eV  = 1.602192e-12 ; erg
   as.omega = 3.05d-4           ; steradians per square degree
   as.pc  = 3.085678e+18        ; cm ; kpc = 1000*pc; Mpc = 1000*kpc;

   IF n_elements(help) NE 0 THEN BEGIN 
      print,'as.pi    = 3.14159265358979'
      print,'as.EE    = 2.718281828459045 ; exp(1) '
      print,'as.Msun  = 1.9891e+33        ; g '
      print,'as.Lsun  = 3.826e+33         ; erg/s '
      print,'as.Rsun  = 6.9599e+10        ; cm '
      print,'as.c     = 2.9979245620e+10  ; cm/s ; speed of light '
      print,'as.G     = 6.67259e-08       ; cgs ; gravitational constant'
      print,'as.e     = 4.8032068e-10     ; charge of electron'
      print,'as.h     = 6.6260755e-27     ; erg*s ; Planck constant'
      print,'as.me    = 9.1093897e-28     ; g ; mass of electron'
      print,'as.mp    = 1.6726231e-24     ; g ; mass of proton'
      print,'as.alpha = 1/(as.h*as.c/(2*as.pi*as.e^2)) ; fine str. constant'
      print,'as.sigmaT= 6.6524616e-25     ; cm^2 ; Thomson sigma = 8pi/3*re^2'
      print,'as.k     = 1.380658e-16      ; erg/K ; Boltzman constant'
      print,'as.NA    = 6.0221367e+23     ; mol^-1 ; Avogadro constant'
      print,'as.sigma = 5.67051e-5        ; Stefan-Boltzmann constant'
      print,'as.keV = 1.602192e-9         ; eV  = 1.602192e-12 ; erg'
      print,'as.omega = 3.05d-4           ; steradians per square degree'
      print,'as.pc  = 3.085678e+18        ; cm ; kpc = 1000*pc; Mpc = 1000*kpc'
   ENDIF 

   return, as
END 

