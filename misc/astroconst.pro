;+
; NAME:
;   ASTROCONST()
;
; PURPOSE:
;   Return a large number of astrophysics-related constants. 
;
; KEYWORD PARAMETERS: 
;   mks - default is to return cgs units unless /MKS
;
; OUTPUTS: 
;   c - data structure with all the constants you need
;
; COMMENTS:
;   A_c	speed of light in a vacuum
;   A_k	Boltzmann's constant
;   A_G	universal gravitational constant
;   A_me	mass of electron
;   A_mp	mass of proton
;   A_amu	atomic mass unit
;   A_eV	electron volt
;   A_h	Planck's constant
;   A_hbar	Planck's constant /2pi
;   A_e	electron charge (defined positive)
;   A_mn	neutron mass
;   A_mH	mass of Hydrogen atom
;   A_pc	parsec
;   A_msun	solar mass
;   A_rsun	solar radius
;   A_AU	astronomical unit
;   A_ly	light year
;   A_Lsun	solar luminosity
;   A_mearth	mass of Earth
;   A_rearth	radius of Earth
;   A_AA	Angstrom
;   A_Jy	Jansky
;   A_rhoc	critical density /h^2
;   A_sigma	Stefan-Boltzmann constant
;   A_arad	radiation density constant
;   A_sigmaT	Thomson cross-section
;   A_re	classical electron radius
;   A_a0	Bohr radius
;   A_Wien	Wien displacement law constant
;   A_alpha	fine structure constant
;   A_NA	Avogadro's number
;   A_yr	year
;   A_tH	Hubble time *h
;   A_Tsun	surface temperature of sun
;   A_TCMB	cosmic microwave background temperature
;   A_pi	trig constant pi
;   A_exp	base of natural logarithm
;   A_dsun	mean solar density
;   A_dearth	mean Earth density
;   A_gsun	solar surface gravity
;   A_gearth	Earth surface gravity
;   A_Vsun	solar V magnitude
;   A_MVsun	solar absolute V magnitude
;   A_Z0	characteristic impedance of vacuum
;   A_eps0	permittivity of free space
;   A_mu0	permeability of free space
;   A_rmoon	lunar radius
;   A_mmoon	lunar mass
;   A_amoon	lunar orbital semi-major axis
;   A_emoon	lunar orbital eccentricity
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2002 Aug 1, U of A - based almost entirely on
;     Jeremy Bailin's astroconst.idl
;-

function astroconst, mks=mks

    C = create_struct($
      'c', 2.997925D10,        $
      'k', 1.380658D-16,       $
      'G', 6.67259D-8,         $
      'me', 9.1093897D-28,     $
      'mp', 1.672631D-24,      $
      'amu', 1.6605402D-24,    $
      'eV', 1.6021772D-12,     $
      'h', 6.6260755D-27,      $
      'hbar', 1.0545727D-27,   $
      'e', 4.8032D-10,         $
      'mn', 1.6749286D-24,     $
      'mH', 1.6733D-24,        $
      'pc', 3.08567758D18,     $
      'msun', 1.99D33,         $
      'rsun', 6.96D10,         $
      'AU', 1.496D13,          $
      'ly', 9.463D17,          $
      'Lsun', 3.9D33,          $
      'mearth', 5.976D27,      $
      'rearth', 6.378D8,       $
      'AA', 1D-8,              $
      'Jy', 1D-23,             $
      'rhoc', 1.8791D-29,      $
      'sigma', 5.67040D-5,     $
      'arad', 7.5646D-15,      $
      'sigmaT', 6.6524585D-25, $
      're', 2.817941D-13,      $
      'a0', 5.29177208D-9,     $
      'Wien', 2.897769D-1,     $
      'dsun', 1.408D,          $
      'dearth', 5.515D,        $
      'gsun', 27400D,          $
      'gearth', 978D,          $
      'eps0', 1D,              $
      'mu0', 1D,               $
      'rmoon', 1738D5,         $
      'mmoon', 7.348D25,       $
      'amoon', 3.844D10        $
      )

    if keyword_set(MKS) then begin
       C.c = 2.997925D8
       C.k = 1.380658D-23
       C.G = 6.67259D-11
       C.me = 9.1093897D-31
       C.mp = 1.672631D-27
       C.amu = 1.6605402D-27
       C.eV = 1.6021772D-19
       C.h = 6.6260755D-34
       C.hbar = 1.0545727D-34
       C.e = 1.602D-19
       C.mn = 1.6749286D-27
       C.mH = 1.6733D-27
       C.pc = 3.08567758D16
       C.msun = 1.99D30
       C.rsun = 6.96D8
       C.AU = 1.496D11
       C.ly = 9.463D15
       C.Lsun = 3.9D26
       C.mearth = 5.976D24
       C.rearth = 6.378D6
       C.AA = 1D-10
       C.Jy = 1D-26
       C.rhoc = 1.8791D-26
       C.sigma = 5.67040D-8
       C.arad = 7.5646D-16
       C.sigmaT = 6.6524585D-29
       C.re = 2.817941D-15
       C.a0 = 5.29177208D-11
       C.Wien = 2.897769D-3
       C.dsun = 1408
       C.dearth = 5515
       C.gsun = 274
       C.gearth = 9.78
       C.eps0 = 8.854187817D-12
       C.mu0 = 1.256637061D-6
       C.rmoon = 1738D3
       C.mmoon = 7.348D22
       C.amoon = 3.844D8
    endif

; append dimensionless numbers
    C = create_struct(C, $
      'alpha', (1.0D/137.03599) , $
      'NA', 6.0221367D23        , $
      'yr', 3.156D7             , $
      'tH', 3.0853056D17        , $
      'Tsun', 5780D             , $
      'TCMB', 2.726D            , $
      'pi', 3.14159265358979324D, $
      'exp', 2.71828182846D     , $
      'Vsun', -26.74D           , $
      'MVsun', 4.83D            , $
      'Z0', 376.730313461D      , $
      'emoon', 0.0549D            $
      )

return, C
end
