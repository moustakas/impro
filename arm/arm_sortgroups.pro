;+
; NAME: ARM_SORTGROUPS()
;       
; CATEGORY: array manipulation
;
; PURPOSE: sort array elements while preserving desired groupings
;
; CALLING SEQUENCE: result = ARM_SORTGROUPS()
;
; INPUTS:
;   array - 1D array of N elements
;       
; OPTIONAL INPUTS:
;   group - 1D array of N elements where elements of equal value
;           indicate members of a particular grouping
;
; KEYWORDS:
;   fixed   - sort groups, but do not sort elements within groups
;   extract - return sorted array rather than array of sorted indices
;
; OUTPUTS: returns array indices in sorted order (see EXTRACT keyword)
; 
; OPTIONAL OUTPUTS:
;
; EXAMPLE:  IDL> a = [9, 8, 7, 6, 5, 4, 3, 2, 1]
;           IDL> b = ['c', 'c', 'a', 'a', 'a', 'a', 'b', 'b', 'c']
;           IDL> c = ARM_SORTGROUPS(a, b, /extract) 
;           IDL> d = ARM_SORTGROUPS(a, b, /extract, /fixed)
;            
;                c -> [1, 8, 9, 2, 3, 4, 5, 6, 7]
;                d -> [9, 8, 1, 3, 2, 7, 6, 5, 4]
;
; PROCEDURES USED:
;
; COMMENTS: Useful, for example, if you have a list of paired entries
;           which you want to sort but also keep the pairs together.
;
; MODIFICATION HISTORY:
;    written by A.R.Marble, Steward Obs., June 9, 2004
;-

function arm_sortgroups, array, group, fixed=fixed, extract=extract
    
    if N_ELEMENTS(group) eq 0L then group = INDGEN(N_ELEMENTS(array))
    
    groups = group[SORT(group)]

    groups = groups[UNIQ(groups)]
    
    n = N_ELEMENTS(groups)
    
    grouprep = array[0:n-1]
    
    for i=0,n-1 do begin
       
       wh = WHERE(group eq groups[i])
       
       grouprep[i] = array[wh[(SORT(array[wh]))[0]]]
       
    endfor
    
    grouporder = SORT(grouprep)
    
    for i=0,n-1 do begin
       
       wh = WHERE(group eq groups[grouporder[i]], count)
       
       if not KEYWORD_SET(fixed) then suborder = SORT(array[wh]) $
       else suborder = INDGEN(count)
       
       if i eq 0 then order = wh[suborder] else order = [order, wh[suborder]]
       
    endfor
    
    if KEYWORD_SET(extract) then return, array[order] else return, order
    
 end
       







