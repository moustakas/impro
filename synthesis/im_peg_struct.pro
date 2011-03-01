function im_peg_struct, ntime, ncont, nlines
; jm05oct18uofa - compatible with IM_READ_PEG
; jm10nov04ucsd - replaced MALLSTARS with MSTARLIVE
    
    peg = {$
      file:          '',            $
      nage:          ntime,         $
      ncont:         ncont,         $
      nlines:        nlines,        $
      age:           0.0D,          $ ; need double precision
      mgalaxy:       0.0,           $
      mstar:         0.0,           $ ; M_star+M_wd+M_nsbh+M_substellar
      mwd:           0.0,           $ ; white dwarfs
      mnsbh:         0.0,           $ ; neutron stars, black holes
      msubstellar:   0.0,           $
      mstarlive:     0.0,           $ 
      mgas:          0.0,           $
      Zgas:          0.0,           $
      Zstar_mass:    0.0,           $
      Zstar_lbol:    0.0,           $
      lbol:          0.0D,          $ ; double
      tauv:          0.0,           $
      ldust_lbol:    0.0,           $
      sfr:           0.0,           $
      nlyc:          0.0D,          $ ; double
      sniirate:      0.0,           $
      sniarate:      0.0,           $
      starage_mass:  0.0,           $ ; Myr
      starage_lbol:  0.0,           $ ; Myr
      wave:          fltarr(ncont), $
      flux:          dblarr(ncont), $ ; need double precision
      linewave:      fltarr(nlines),$
      lineflux:      dblarr(nlines)}  ; need double precision
    peg = replicate(peg,ntime)

return, peg
end
    
