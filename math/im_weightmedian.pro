; taken from http://physics.open.ac.uk/~crpowell/idl.html
FUNCTION weightmedian, Y, dY, mean= MEDIAN
  Case N_params(0) Of
    0:Return, 0b
    1:Begin
      If ( Size( Y, /n_dimensions ) Eq 2 ) And ( Min( Size( Y, /dimensions ) ) Gt 1 ) Then Begin
        m = Size( Y, /dimensions )
        If m[ 0 ] Lt m[ 1 ] Then Begin
          yy     = Double( Y[ 0, * ] )
          dyy    = Double( Y[ 1, * ] )
         Endif Else Begin
          yy     = Double( Y[ *, 0 ] )
          dyy    = Double( Y[ *, 1 ] )
         Endelse
       Endif Else Begin
        yy     = Double( Y )
        dyy    = Double( Finite( Y ) )
       Endelse
     End
    Else:Begin
      yy     = Double( Y )
      dyy    = Double( dY )
     End
   Endcase
  If N_elements( dyy ) Eq 1 Then dyy = yy * 0 + Total( dyy )

  If Total( dyy Eq 0 ) Ne 0 Then Begin
    Print,"Some errors are zero; setting them to 1 to enable output."
    dyy[ Where( dyy Eq 0 ) ] = 1.
   Endif

  If Total( Finite( yy ) ) + Total( Finite( dyy ) ) Ne 2*N_elements( yy ) Then Begin
    filter = Where( Finite( yy ) And Finite( dyy ) )
    If filter[ 0 ] Eq -1 Then Begin
      Print,"No finite elements available for median."
      Return, 0b
     Endif
    yy = yy[ filter ]
    dyy = dyy[ filter ]
   Endif

  order  = Sort( yy )
  yy     = yy[  order ]
  dyy    = dyy[ order ]
stop
  weight = dyy ^ (-2)
  s      = Total( weight )
  For _0 = 0, N_elements( yy ) - 1 Do Begin
    If Total( weight[ 0:_0 ] ) Ge 0.5 * s Then Break
   Endfor

  Median = Total( weight[ 0:_0 ] ) Eq 0.5 * s ? 0.5 * ( yy[ _0 ] + yy[ _0+1 ] ) : yy[ _0 ]

  Return, MEDIAN
 End

PRO im_weightmedian, Y, dy, Median

  If N_params(0) Eq 0 Then Return
  dyy = N_params(0) Eq 1 ? Finite( Y ) : dY
  Median = Weightmedian( Y, dyy)

return
End
