
pro setplot, device, filename=file, $
  portrait=portrait, landscape=landscape, square=square, $
  banner=banner, eps=eps, color=color, guider=guider, h2guider=h2guider, talk=talk

;+
; NAME:
;	SETPLOT
;	setplot, device
;	   where Device can be 'hirez','modgr','ps','sun','tek','X',
;		               'NULL' or 'HP'
;          if Device is omitted, current information will be printed.
; INPUTS:
; OPTIONAL INPUT PARAMETERS:
; KEYWORD PARAMETERS:
; OUTPUTS:
; SIDE EFFECTS:
;	If the current graphics device is 'PS', a 'device, /close' is done.
; PROCEDURE:
; MODIFICATION HISTORY:
;       Modified extensively by Tom Bida.	

;   if not(keyword_set(filename)) then file='idl.ps'
;   if ((device eq 'ps') and (keyword_set(banner))) then begin
;     spawn,'date', date
;     xyouts, 0, 0, date+'  '+ file, /dev, size=0.5
;   endif

    on_error, 2
    ttype= strlowcase(strmid(Device, 0, 2))
    if !d.name eq 'PS' then device, /close
    com0 = ""
    case ttype of
      'ps':  begin
             com0 = "device, bits_per_pixel=8, /helvetica" 
             kw2 = ""
	     kw3 = ""
             sinit = "set_plot,'ps' & !p.font=0"
	       if keyword_set(file) and keyword_set(eps) then begin
	        kw1= ", /encapsulated, filename='" + file + "'"
	       endif
	       if keyword_set(file) and not(keyword_set(eps)) then begin
	        kw1 = ", filename='" + file + "'"
	       endif
	       if not(keyword_set(file)) then begin
	        file = 'idl.ps'
	        kw1 = ", filename='" + file + "'"
	       endif

	       if keyword_set(portrait) then begin
	        kw2 = ", /portrait, /inches, yoffset=0.6, ysize=9.75"
               endif
	       if keyword_set(landscape) then begin
	        kw2 = ", /landscape, /inches"
               endif
	       if keyword_set(square) then begin
	        kw2 = ", /portrait, /inches, yoffset=3.0, ysize=6.68"
               endif
	       if keyword_set(guider) then begin
                width = 6.5
		height = 6.5*288./384.*7/8.
	        kw2=",/portrait,/inches,xsize=width,ysize=height,yoffset=4
               endif
	       if keyword_set(h2guider) then begin
                width = 6.5
		height = 6.5*7/8.
	        kw2=",/portrait,/inches,xsize=width,ysize=height,yoffset=4
               endif
	       if keyword_set(color) then kw3=",/color"
	     end
      else:  begin
	       sinit = "set_plot, '" + Device + "'"
	     end
    endcase

    if (keyword_set(talk)) then print,' ', sinit
    jerr = execute(sinit)
    if strlen(com0) gt 0 then begin
        dcom0 = com0 + kw1(0)
        dcom1 = 'device' + kw2 + kw3
        if (keyword_set(talk)) then print,' ', dcom0
        if (keyword_set(talk)) then print, ' ', dcom1
        err = execute(dcom0)
        err = execute(dcom1)
    endif

    return
end
