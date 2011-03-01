pro im_lineztelluric
; jm05jan19uofa
; at what redshifts do emission lines of interest fall within the
; prominent telluric bands?

    linename = ['[O II]','H-beta','[O III]','H-alpha','[N II]']
    linewave = [3727.0,4861.0,5007.0,6563.0,6584.0]
    nline = n_elements(linewave)
    
    tellbands1 = { TELLBAND, $
      twave1: 6850., $
      twave2: 6960., $
      cwave1: [6600., 6950., 0], $
      cwave2: [6860., 7200., 0] }
    tellbands2 = { TELLBAND, $
      twave1: 7150., $
      twave2: 7350., $
      cwave1: [7050., 7115., 7340.], $
      cwave2: [7160., 7130., 7440.] }
    tellbands3 = { TELLBAND, $
      twave1: 7560., $
      twave2: 7720., $
      cwave1: [7400., 7700., 0], $
      cwave2: [7580., 8000., 0] }
    tellbands4 = { TELLBAND, $
      twave1: 8105., $
      twave2: 8240., $
      cwave1: [8000., 8225., 0], $
      cwave2: [8105., 8325., 0] }
    tellbands5 = { TELLBAND, $
      twave1: 8530., $
      twave2: 8865., $
      cwave1: [8490., 8865., 0], $
      cwave2: [8530., 8905., 0] }
    tellbands6 = { TELLBAND, $
      twave1: 8644., $
      twave2: 8697., $
      cwave1: [8604., 8697., 0], $
      cwave2: [8644., 8737., 0] }

    result = {$
      linename:          '', $
      linewave:         0.0, $
      z_band1:    [0.0,0.0], $
      z_band2:    [0.0,0.0], $
      z_band3:    [0.0,0.0], $
      z_band4:    [0.0,0.0], $
      z_band5:    [0.0,0.0], $
      z_band6:    [0.0,0.0]}
    result = replicate(result,nline)

    result.linename = linename
    result.linewave = linewave
    
    for iline = 0L, nline-1L do begin

       result[iline].z_band1 = [tellbands1.twave1,tellbands1.twave2]/result[iline].linewave-1.0
       result[iline].z_band2 = [tellbands2.twave1,tellbands2.twave2]/result[iline].linewave-1.0
       result[iline].z_band3 = [tellbands3.twave1,tellbands3.twave2]/result[iline].linewave-1.0
       result[iline].z_band4 = [tellbands4.twave1,tellbands4.twave2]/result[iline].linewave-1.0
       result[iline].z_band5 = [tellbands5.twave1,tellbands5.twave2]/result[iline].linewave-1.0
       result[iline].z_band6 = [tellbands6.twave1,tellbands6.twave2]/result[iline].linewave-1.0
       
    endfor

    struct_print, result

stop    
    
return
end
    
