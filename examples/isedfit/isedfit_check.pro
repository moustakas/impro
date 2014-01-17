;+
; ISEDFIT_CHECK -- Check the current installed versions for debug
;
; e.g. IDL> check_isedfit
; -

PRO _req, cmd, name, version
  ;; a little wrapper that runs a command to ensure that the installed
  ;; version matches the required version passed into the function.
  value = call_function(cmd)
  status = (value EQ version) ? 'ok' : 'Failed!'
  print, 'iSEDfit requires '+version+' version of '+name+'. Found: '+value+'. '+status
  
  IF value NE version THEN BEGIN
    message, ' !!! Version did not match the required version, please update'
  ENDIF
END


PRO isedfit_check
  ;; might also want to check directories ($ISEDFIT_DIR,
  ;; $KCORRECT_DIR) as well.
  _req, 'idlutils_version', 'IDLUtils', 'trunk'
  _req, 'k_version', 'kCorrect', 'NOCVS:kcorrect'
END

