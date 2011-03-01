;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;
;;;  NAME:
;;;     
;;;     sort2d.pro (written by Andy Marble, Steward Observatory, July 2001)
;;;
;;;  DESCRIPTION:
;;;
;;;     sorts entire table according to data in specified column(s)
;;;
;;;  INPUT:  
;;;
;;;     table   - 2 dimensional table of numerical data
;;;     columns - array of column numbers to sort (must be 1 or greater)
;;;
;;;  KEYWORDS:
;;;
;;;     row     - sorts table by row instead of column
;;;               (also appropriate to set this parameter if the table
;;;                is indexed as [row#,col#] instead of [col#,row#])
;;;     reverse - if set, sorts table in descending order
;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro sort2d,table,columns,REVERSE=backwards,ROW=row,order=order

;   CHECK CALLING SEQUENCE
    if not keyword_set(table) then begin
       print & print,'ERROR: invalid function call'
       print & print,'       sort2d, table, column, /REVERSE, /ROW' & print
       return
    endif

;   CHECK THAT TABLE HAS TWO DIMENSIONS
    if n_elements(table[*,0]) lt 2 and n_elements(table[0,*]) lt 2 then begin
       print & print,'ERROR: table must be 2 dimensional' & print
       return
    endif

;   DETERMINE NUMBER OF ROWS OR COLUMNS
    n=[n_elements(table[*,0]),n_elements(table[0,*])]
    if keyword_set(row) then n=reverse(n)

;   CHECK COLUMNS, IF UNDEFINED SORT BY INCREASING COLUMN NUMBER
    if not keyword_set(columns) then columns=indgen(n[0])+1 else begin
       wh=where(columns lt 1 or columns gt n[0])
       if wh[0] ne -1 then begin
          print & print,'ERROR: invalid column/row number' & print
          return
       endif
    endelse
    
;   DETERMINE SORTING ORDER
    table2=table
    order=indgen(n[1])
    for j=n_elements(columns)-1,0,-1 do begin
       if not keyword_set(row) then begin
          order=order[sort(table2[columns[j]-1,order])]
       endif else begin
          order=order[sort(table2[order,columns[j]-1])]
       endelse
    endfor
    
;   REVERSE SORTING ORDER IF APPROPRIATE
    if keyword_set(backwards) then order=reverse(order)

;   SORT THE RETURNED TABLE
    if not keyword_set(row) then for i=0,n[0]-1 do table2[i,*]=table2[i,order] $
                            else for i=0,n[0]-1 do table2[*,i]=table2[order,i]
    table=table2
    return

 end








