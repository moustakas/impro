pro which,proname
;Prints full filenames in IDL !path search order for a particular routine.
; proname (input string) procedure name (.pro will be appended) to find
;24-Aug-92 JAV	Create.
;10-Mar-93 JAV	Fixed bug; last directory in !path ignored; pad with ': '

if n_params() lt 1 then begin
  print,'syntax: which,proname(.pro assumed)'
  retall
endif

  pathlist = '.:' + !path + ': '		;build IDL path list
  fcount = 0					;reset file counter
  il = strlen(pathlist) - 1			;length of path string
  ib = 0					;begining substring index
  ie = strpos(pathlist,':',ib)			;ending substring index
  repeat begin					;true: found path separator
    path = strmid(pathlist,ib,ie-ib)		;extract path element
    fullname = path + '/' + proname + '.pro'	;build full filename
    openr,unit,fullname,error=eno,/get_lun	;try to open file
    if eno eq 0 then begin			;true: found file
      fcount = fcount + 1			;increment file counter
      if path eq '.' then begin			;true: in current directory
	spawn,'pwd',dot				;get current working directory
	dot = dot(0)				;convert to scalar
	print,fullname + ' (. = ' + dot + ')'	;print filename + current dir
      endif else begin				;else: not in current directory
	print,fullname				;print full name
      endelse
      free_lun,unit				;close file
    endif
    ib = ie + 1					;point beyond separator
    ie = strpos(pathlist,':',ib)		;ending substring index
    if ie eq -1 then ie = il			;point at end of path string
  endrep until ie eq il				;until end of path reached
  if fcount eq 0 then begin			;true: routine not found
    print,'which: ' + proname + '.pro not found on IDL !path.'
  endif
end
