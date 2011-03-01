pro telluric_bands_vs_redshift
; jm08mar26nyu - print the redshifts where the emission lines of
;                interest fall in the telluric bands

    zmax = 1.2 & zmin = 0.0 & dz = 0.05
    zgrid = findgen((zmax-zmin)/dz+1)*dz+zmin
    lamgrid = findgen((1E4-3500.0)/1.0+1)*1.0+3500.0

    lamrest = [3727.0,4861.0,5007.0,6563.0]
    linename = ['[O II]','H-BETA','[O III]','H-ALPHA']
    tellbands1 = {twave1: 6850., twave2: 6960.0}
    tellbands2 = {twave1: 7150., twave2: 7350.0}
    tellbands3 = {twave1: 7560., twave2: 7720.0}
    tellbands4 = {twave1: 8105., twave2: 8240.0}
    tellbands = [tellbands1,tellbands2,tellbands3,tellbands4]

    for iline = 0L, n_elements(lamrest)-1L do begin
       print, linename[iline]+':'
       print, ' Telluric Band 1: '+strtrim(string(interpol(zgrid,lamrest[iline]*(1+zgrid),tellbands1.twave1),format='(F12.4)'),2)+'-'+$
         strtrim(string(interpol(zgrid,lamrest[iline]*(1+zgrid),tellbands1.twave2),format='(F12.4)'),2)
       print, ' Telluric Band 2: '+strtrim(string(interpol(zgrid,lamrest[iline]*(1+zgrid),tellbands2.twave1),format='(F12.4)'),2)+'-'+$
         strtrim(string(interpol(zgrid,lamrest[iline]*(1+zgrid),tellbands2.twave2),format='(F12.4)'),2)
       print, ' Telluric Band 3: '+strtrim(string(interpol(zgrid,lamrest[iline]*(1+zgrid),tellbands3.twave1),format='(F12.4)'),2)+'-'+$
         strtrim(string(interpol(zgrid,lamrest[iline]*(1+zgrid),tellbands3.twave2),format='(F12.4)'),2)
       print, ' Telluric Band 4: '+strtrim(string(interpol(zgrid,lamrest[iline]*(1+zgrid),tellbands4.twave1),format='(F12.4)'),2)+'-'+$
         strtrim(string(interpol(zgrid,lamrest[iline]*(1+zgrid),tellbands4.twave2),format='(F12.4)'),2)
       print
    endfor
    
return
end
    
