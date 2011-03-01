;+
; NAME: ARM_YESORNO
;       
; CATEGORY: miscellaneous
;
; PURPOSE: prompt user for negative or affirmative response
;
; CALLING SEQUENCE: ARM_YESORNO, query, [choice1=, choice2=]
;
; INPUTS: query - message to be displayed when prompting user
;       
; OPTIONAL INPUTS:
;   choice1 - first of two valid case-insensitive responses
;             (default='Y') 
;   choice2 - second of two valid case-insensitive responses
;             (default='N')
;
; OUTPUTS: returns 1 (true) for affirmative resonse (choice1), 
;          returns 0 (false) for negative resonse (choice2)
;
; MODIFICATION HISTORY:
;    written by A.R.Marble, Steward Obs., 2003
;    default query value added by ARM, October 23, 2003
;-

function arm_yesorno, query, choice1=choice1, choice2=choice2
    
; defaults

    if not KEYWORD_SET(choice1) then choice1 = 'Y' else $
      choice1 = STRTRIM(STRING(choice1), 2)
    if not KEYWORD_SET(choice2) then choice2 = 'N' else $
      choice2 = STRTRIM(STRING(choice2), 2)
    if not KEYWORD_set(query) then query = "Please select 'yes' or 'no'."

; print message prompting user for input

    PRINT
    PRINT, query+' ('+choice1+'/'+choice2+')'

    answer = ''                 ; initialize input
 
; prompt for input until valid response is given

    while STRUPCASE(answer) ne STRUPCASE(choice1) $ 
      and STRUPCASE(answer) ne STRUPCASE(choice2) do begin

       READ, answer             ; accept input from user
       answer = STRTRIM(STRING(answer), 2)

       if STRUPCASE(answer) eq STRUPCASE(choice1) then return, 1L
       if STRUPCASE(answer) eq STRUPCASE(choice2) then return, 0L
 
    endwhile
    
 end

