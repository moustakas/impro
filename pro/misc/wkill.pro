;+
; NAME:
;   WKILL
; PURPOSE:
;   Delete all open windows.
; MODIFICATION HISTORY:
;   J. Moustakas, 2000 Sep 08, U of A
;-

pro wkill
while !d.window ne -1L do wdelete, !d.window ; delete all windows
return
end
