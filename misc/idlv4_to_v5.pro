pro idlv4_to_v5,infiles,outdir
;+
; NAME:
;	IDLV4_TO_V5
; PURPOSE:
;	Modify an IDL V4.0 (or earlier) procedure such that variables are 
;	indexed using square brackets, as allowed (and suggested) 
;	within IDL V5.0 and later
;
; CALLING SEQUENCE:
;	IDLV4_TO_V5, infiles, outdir 
;
; INPUTS:
;	infiles - scalar string specifying IDL procedure name(s), wild card 
;		values allowed
;	outdir - scalar string giving directory to contain output file.
;
; EXAMPLES:
;	Convert the procedure curvefit.pro in the current directory to a
;	IDL V5 version in the (Unix) idlv5 directory
;
;	IDL> idlv4_to_v5,'curvefit.pro','idlv5/'
;
;	Convert all the procedures in the current directory to IDL V5 versions
;	in the /share/idlv5 directory
;
;	IDL> idlv4_to_v5, '*.pro', '/share/idlv5/'
;
; METHOD:
;	ISFUNCTION() is used to determine all the routine names in the file,
;	and then ROUTINE_INFO() is used to determine the names of all variables
;	in the procedure.    Each (non-commented) line is scanned for
;	parentheses, and converted to square brackets if the token to the left
;	of the left parenthesis matches a variable name.
; 
; NOTES:
;	(1) Only runs under IDL V5.0 (since it calls ROUTINE_INFO())
;	(2) May possibly get confused by parenthesis within strings.
;	(3) May get confused by IDL statements that extend over multiple lines
;	    idlv4_to_v5 will supply a warning when it becomes confused by
;	    unmatched parenthesis.
;	(4) Do not include this procedure 'idlv4_to_v5' in the directory that 
;	    you are trying to convert (since it will compile the procedure 
;	    while executing it, and do a retall.)
;	(5) Conversions cannot be performed unless specified procedure(s) 
;	    already compile properly
;	(6) Will not work on IDL main programs
;	(7) May get confused by gaps between array name and parenthesis
;
; PROCEDURES CALLED:
;	FDECOMP, MATCH, REMOVE, ISFUNCTION()
; REVISION HISTORY:
;	Written  W. Landsman   Hughes STX     June 1997 
;	Variable names can have numerals      August 1997
;	Never change an intrinsic IDL function to square brackets, even if it
;	is also a variable name.
;-

 if N_params() LT 2 then begin
	print,'Syntax  - idlv4_to_v5, infiles, outdir'
	return
 endif

 a = findfile(infiles,count=n)
 if n EQ 0 then message,'No files found ' + infiles
 get_lun,inlun
 get_lun,outlun

funcnames = routine_names(S_functions=-1)

 line = ''

;Loop variables
;i - loop over all filenames (if wildcard value of 'infile' supplied)
;k - loop over all routines in the current filename
;kk -loop over all lines in the current routine
;j - loop over all left parentheses in the current line
;jj- loop over all right parentheses in the current line

 for i=0,n-1 do begin            ;loop over each procedure name

	fdecomp,a[i],disk,dir,name,ext    ;Decompose file name
        status = isfunction(a[i], outnames,numline) 
;Resolve main procedure first, even if it mean compiling twice
	g = where(outnames EQ strtrim(strlowcase(name),2),Ng)
	if Ng GT 0 then begin
		g = g[0] 
		if status[g] then resolve_routine,outnames[g],/is_function $
		else resolve_routine, outnames[g]
	endif
	

	for k = 0, N_elements(status)-1 do begin

	case status[k] of
	 1: begin
		resolve_routine,outnames[k],/is_function
		variables = routine_info(/variables,outnames[k],/functions)
	    end
	 0: begin
		resolve_routine,outnames[k]
		variables = routine_info(/variables,outnames[k])
		end
	-1: begin
	    message,a[i] + ' will not be modified',/INF
	    goto, Done_pro 
	    end
	endcase
        
 match, variables, funcnames, subv, Count = Nfunc
 if Nfunc GT 0 then remove,subv,variables

 if k EQ 0 then begin
	openr,inlun,a[i]
	openw,outlun, outdir + name + '.pro'
 endif

	for kk=0,numline[k]-1 do begin
		readf,inlun,line
		len = strlen(line)
		pos = strpos( line, ';')
		if pos EQ -1 then begin
			goodline = line
			comment = ''
		endif else begin
			goodline = strmid(line,0, pos)
			comment = strmid(line,pos,len-pos)
		endelse
		if goodline EQ '' then goto,Done_line
	    	
		bchar = byte(goodline)

		leftparen = where(bchar EQ 40b, Nparen)
		if Nparen EQ 0 then goto, Done_line

;Variable names can contain letters, digits, underscore or a dollar sign.
;To allow structure tags and system variables, we include a period and a !
		n = strlen( goodline )

		mask = bytarr(n)
		ii = WHERE( ((bchar GE 65B) and (bchar LE 90b)) OR $
	        ((bchar GE 97B) and (bchar LE 122B)) OR $
                ((bchar GE 48B) and (bchar LE 57B)) OR $
              (bchar EQ 46B) or (bchar EQ 36B) OR $
              (bchar EQ 41B) OR $
              (bchar EQ 95B) or (bchar EQ 33B), count)
		if count GT 0 then mask[ii] = 1b else goto, Done_Line
 		pconvert = bytarr(Nparen)  ;Keep track of which paren to convert

; Now we step backward from the left parenthesis until we find the first 
; character that cannot be part of a variable name

		for j = 0, Nparen - 1 do begin
			mark =  leftparen[j] - 1
			if mark EQ -1 then goto,Done_paren
			while mask[mark] do begin
				mark = mark - 1b
				if mark EQ -1 then goto, Done_search
			endwhile
		done_search:

			if mark EQ leftparen[j]-1 then goto, Done_paren
			varname = strtrim(bchar[mark+1:leftparen[j]-1],2)
			if varname EQ '' then goto, Done_paren
; Test for structure name.   Note that for a structure x, that x.tag[3] is 
; legal in V5.0 but x.[3] is not (although x.(3) is).

			dot = strpos(varname,'.')    ;Structure name
			if dot EQ strlen(varname)-1 then goto,Done_paren
			if dot GT 0 then varname = strmid(varname,0,dot)
			g = where(variables EQ strupcase(varname), Ng)
			if Ng GT 0 then pconvert[j] = 1b
			if strmid(strtrim(varname,2),0,1) EQ '!' then $
				pconvert[j] = 1b				
		Done_paren:
		endfor
		convert = where(pconvert, Nconvert)

		if Nconvert GT 0 then begin
			bchar[leftparen[convert]] = 91b    ;byte('[')=91b
			rparen = where(bchar EQ 41b, Nrparen)
			if Nrparen EQ 0 then begin 
				message, 'Warning - no right parenthesis',/INF
				print,goodline
				goto,done_line
			endif

			for jj = 0, Nrparen - 1 do begin
			g = where(leftparen LT rparen[jj], Ng)
			if Ng EQ 0 then begin
				message,'Warning - missing left parenthesis',/INF
				print,goodline
				goto, done_line
			endif 
			leftindex = max(g)
			if pconvert[leftindex] then  bchar[rparen[jj]] = 93b 
			if N_elements(leftparen) GT 1 then $ 
			remove,leftindex,leftparen,pconvert $
			else goto, Done_rparen
			endfor
		endif
done_rparen:
	goodline = string(bchar)
   Done_line:  

		printf,outlun,goodline + comment
   endfor
 endfor
   close,inlun
   close,outlun
Done_pro:
   endfor
   free_lun,inlun,outlun
   return
   end
