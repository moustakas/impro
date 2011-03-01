pro menu_template

;+
; NAME:
;       MENU_TEMPLATE
;
; PURPOSE:
;       is a template for menu widgets
;
; CATEGORY:
;       XX Package; widgets
;
; CALLING SEQUENCE:
;
;
; INPUTS:
;
;
; OUTPUTS:
;
;
; COMMON BLOCKS:
;
;
; PROCEDURE:
;
;
; MODIFICATION HISTORY:
;       2-April-1997, Herve Dole, IAS
;
;-

@fonts

; NOTE: the name in the field "UVALUE" is the name of the procedure called !!

; CREATION OF THE MENU STRUCTURE
;-------------------------------
   wMainWindow = WIDGET_BASE(TITLE="IRAS Map Maker", MBAR=wMenuBar, TLB_FRAME_ATTR=1, XSIZE = 500)
   wMainMenu = WIDGET_BASE(wMainWindow, /COLUMN)

; FILE MENU
;-----------
   wFileMenu = WIDGET_BUTTON(wMenuBar, VALUE='File', /MENU, FONT = HELB18)
        wFileItem = WIDGET_BUTTON(wFileMenu, VALUE='File_ssmenu1', UVALUE='File_ssmenu1', FONT = HELB18I )
        wFileItem = WIDGET_BUTTON(wFileMenu, VALUE='File_ssmenu2', UVALUE='File_ssmenu2', FONT = HELB18I )
        wFileItem = WIDGET_BUTTON(wFileMenu, VALUE='File_ssmenu3', UVALUE='File_ssmenu3', FONT = HELB18I )
        wFileItem = WIDGET_BUTTON(wFileMenu, VALUE='File_ssmenu4', UVALUE='File_ssmenu4', FONT = HELB18I )
        wFileItem = WIDGET_BUTTON(wFileMenu, VALUE='File_ssmenu5', UVALUE='File_ssmenu5', FONT = HELB18I )
        wExitItem = WIDGET_BUTTON(wFileMenu, VALUE='Exit', UVALUE='EXIT', FONT = HELB18I )

; FIRST MENU
;-----------
   wMenu1Menu = WIDGET_BUTTON(wMenuBar, VALUE='Selection', /MENU, FONT = HELB18)
        wMenu1Item = WIDGET_BUTTON(wMenu1Menu, VALUE='menu1_ssmenu1...', UVALUE='menu1_ssmenu1', FONT = HELB18I )
        wMenu1Item = WIDGET_BUTTON(wMenu1Menu, VALUE='menu1_ssmenu2...', UVALUE='menu1_ssmenu2', FONT = HELB18I )
        wMenu1Item = WIDGET_BUTTON(wMenu1Menu, VALUE='menu1_ssmenu3...', UVALUE='menu1_ssmenu3', FONT = HELB18I )
        wMenu1Item = WIDGET_BUTTON(wMenu1Menu, VALUE='menu1_ssmenu4...', UVALUE='menu1_ssmenu4', FONT = HELB18I )
        wMenu1Item = WIDGET_BUTTON(wMenu1Menu, VALUE='menu1_ssmenu5...', UVALUE='menu1_ssmenu5', FONT = HELB18I )

; SECOND MENU
;-----------
   wMenu2Menu = WIDGET_BUTTON(wMenuBar, VALUE='Create Map', /MENU, FONT = HELB18)
        wMenu2Item = WIDGET_BUTTON(wMenu2Menu, VALUE='Create Map ...', UVALUE='menu2_ssmenu1', FONT = HELB18I )
        wMenu2Item = WIDGET_BUTTON(wMenu2Menu, VALUE='menu2_ssmenu2...', UVALUE='menu2_ssmenu2', FONT = HELB18I )
        wMenu2Item = WIDGET_BUTTON(wMenu2Menu, VALUE='menu2_ssmenu3...', UVALUE='menu2_ssmenu3', FONT = HELB18I )
        wMenu2Item = WIDGET_BUTTON(wMenu2Menu, VALUE='menu2_ssmenu4...', UVALUE='menu2_ssmenu4', FONT = HELB18I )
        wMenu2Item = WIDGET_BUTTON(wMenu2Menu, VALUE='menu2_ssmenu5...', UVALUE='menu2_ssmenu5', FONT = HELB18I )

; THIRD MENU
;-----------
   wMenu3Menu = WIDGET_BUTTON(wMenuBar, VALUE='Display', /MENU, FONT = HELB18)
        wMenu3Item = WIDGET_BUTTON(wMenu3Menu, VALUE='menu3_ssmenu1...', UVALUE='menu3_ssmenu1', FONT = HELB18I )
        wMenu3Item = WIDGET_BUTTON(wMenu3Menu, VALUE='menu3_ssmenu2...', UVALUE='menu3_ssmenu2', FONT = HELB18I )
        wMenu3Item = WIDGET_BUTTON(wMenu3Menu, VALUE='menu3_ssmenu3...', UVALUE='menu3_ssmenu3', FONT = HELB18I )
        wMenu3Item = WIDGET_BUTTON(wMenu3Menu, VALUE='menu3_ssmenu4...', UVALUE='menu3_ssmenu4', FONT = HELB18I )
        wMenu3Item = WIDGET_BUTTON(wMenu3Menu, VALUE='menu3_ssmenu5...', UVALUE='menu3_ssmenu5', FONT = HELB18I )

; FOURTH MENU
;-----------
   wMenu4Menu = WIDGET_BUTTON(wMenuBar, VALUE='Menu4', /MENU, FONT = HELB18)
        wMenu4Item = WIDGET_BUTTON(wMenu4Menu, VALUE='menu4_ssmenu1...', UVALUE='menu4_ssmenu1', FONT = HELB18I )
        wMenu4Item = WIDGET_BUTTON(wMenu4Menu, VALUE='menu4_ssmenu2...', UVALUE='menu4_ssmenu2', FONT = HELB18I )
        wMenu4Item = WIDGET_BUTTON(wMenu4Menu, VALUE='menu4_ssmenu3...', UVALUE='menu4_ssmenu3', FONT = HELB18I )
        wMenu4Item = WIDGET_BUTTON(wMenu4Menu, VALUE='menu4_ssmenu4...', UVALUE='menu4_ssmenu4', FONT = HELB18I )
        wMenu4Item = WIDGET_BUTTON(wMenu4Menu, VALUE='menu4_ssmenu5...', UVALUE='menu4_ssmenu5', FONT = HELB18I )

; HELP MENU
;----------
   wHelpMenu = WIDGET_BUTTON(wMenuBar, VALUE='Help', /HELP, /MENU, FONT = HELB18)
   wInfoItem = WIDGET_BUTTON(wHelpMenu, VALUE='Info...', UVALUE='INFO', FONT = HELB18I )
   wInfoItem = WIDGET_BUTTON(wHelpMenu, VALUE='About...', UVALUE='ABOUT', FONT = HELB18I )

; DISPLAY WAIT CURSOR
;--------------------
   WIDGET_CONTROL, wMainWindow, /HOURGLASS

; MAKE THE WINDOW VISIBLE
;------------------------
   WIDGET_CONTROL, /REALIZE, wMainWindow
   
; Set cursor to arrow cursor
;   DEVICE, /CURSOR_ORIGINAL

; REGISTER THIS APPLICATION WITH THE XMANAGER
;--------------------------------------------
   XMANAGER, "menu_template", wMainWindow, EVENT_HANDLER="menu_template_event"

END

;------------------------------------------------------------------------------------------------------

PRO menu_template_event, event

@fonts

WIDGET_CONTROL, GET_UVALUE=control, event.id

CASE control OF

     "INFO": BEGIN
         infoText = [ $
          "In Construction ..."+ $
          "Thanks ... "+ $
          "HD"]
         ShowInfo, TITLE="Information", GROUP=event.top, $
         WIDTH=80, HEIGHT=24, INFOTEXT=infoText, FONT=HEL14
     ENDCASE

     "ABOUT": BEGIN
         aboutText = [ $
          "This template is written by Herve Dole, IAS"+ $
          " dole@ias.fr"+ $
          " on 2-April-1997"]

         ShowInfo, TITLE="Information", GROUP=event.top, $
           WIDTH=80, HEIGHT=1, INFOTEXT=aboutText, FONT=HEL14

     ENDCASE

     "EXIT": WIDGET_CONTROL, event.top, /DESTROY
    
        ; User selected an item from the menu - 
        ; execute the proper procedure that was stored in the user value 
        ; for each menu item.
      ELSE:  returnValue = EXECUTE(control + ', GROUP=event.top')

ENDCASE


END
