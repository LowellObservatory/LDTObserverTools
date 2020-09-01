
function dvtrim,filename,noread=noread,trim=trim,sub=sub,nofit=nofit, $
oscana=afit,odataa=ava,jrange=jrange,cgain=cgain,hot=hot,oscanb=bfit, $
odatab=avb,gaincor=gaincor,postpix=postpix
	
; OPTIONAL INPUT KEYWORDS:
;	noread		set to designate filename as an input image;
;			note lack of header though
;	trim		Output a 21-row strip of the image.
;	sub		Subtract overscan.
;	nofit		Do not fit overscan to polynomial; smooth 21 pixels.
;	oscana,b	Return overscan vectors.
;	jrange		Delimit number of rows.
;	gain		Apply gain correction.
;	hot		Median filter hot columns.
;	gaincor		Apply gain boost to right (a) side, value.
;

; Parameters for DeVeny:
; 2015 DD:
  darkcol = 0
  nxpix = 2048
  prepix = 50
  postpix = 50
; e2v CCD 4210:
; darkcol = 0
; nxpix = 2040
; prepix = 58
; Up to 20 March 2006:
;postpix = 50
; After 20 March 2006:
; if not(keyword_set(postpix)) then postpix = 102

if keyword_set (noread) then begin
 image=filename
endif else begin
 header=headfits(filename)
 image=readfits(filename)
 image = float(image)
 fimage=image
endelse

;if keyword_set (jrange) then begin
 image=image(*,12:511)
;endif
sz=size(image)
ny=sz(2)

; set parameters
hotcol=[2007,2008,1128,1129]-1
numamps = sxpar(header,'NUMAMP')
;postpix = fix(sxpar(header,'POSTSCAN')*numamps)
;prepix = fix(sxpar(header,'PRESCAN')*numamps) + darkcol

;if sz(1) gt 2303 then begin
; numamps=2 & prepix=42 & postpix=468 	; HIRES sys2 dual-amp
;endif else begin
; numamps=1 & prepix=21 & postpix=234	; HIRES sys2 single-amp
;endelse

; amplifier b is left side of image; amp a right side
if keyword_set(cgain) then begin
  gaina=1.0 & gainb=1.025
endif else begin
 gaina=1. & gainb=1.	; no gain amplification
endelse

jcnt=sz(2)/2			; half length
  
if keyword_set (gaincor) then begin
 gaina=gaincor & gainb=1.00	; yes gain amplification
endif 

xx=findgen(sz(2))		; dummy array for bias fitting

; BEGIN OVERSCAN SUBTRACTION
if keyword_set(sub) then begin
  if numamps eq 2 then begin

; set columns for bias averages to be at right of postscan

  ib1=prepix+nxpix+darkcol+postpix*0.25
  ib2=prepix+nxpix+darkcol+postpix*0.75

; make average postscan columns and fit with polynomials
  avb=fltarr(sz(2))
  ava=avb

  for iy=0,sz(2)-1 do begin
    avb(iy)=(avg(image(ib1-5:ib1+5,iy)))*gainb
    ava(iy)=(avg(image(ib2-5:ib2+5,iy)))*gaina
  endfor

  fitb=poly_fit(xx,avb,3,bfit)
  fita=poly_fit(xx,ava,3,afit)

  if keyword_set(nofit) then begin
    afit = smooth(ava,21)
    bfit = smooth(avb,21)
  endif

   ixb=sz(1)-postpix-darkcol-nxpix		; left amplifier
   ixa=sz(1)-postpix-darkcol-nxpix/2		; right amplifier

; subtract the fits to the overscan

  for iy=0,sz(2)-1 do begin
   image(ixb:ixb+nxpix/2-1,iy)=(image(ixb:ixb+nxpix/2-1,iy)*gainb-bfit(iy))
   image(ixa:ixa+nxpix/2-1,iy)=(image(ixa:ixa+nxpix/2-1,iy)*gaina-afit(iy))
  endfor

  endif else if numamps eq 1 then begin

; set columns for bias averages to be at right of postscan

  ib2=prepix+nxpix+postpix*0.5

; make average postscan columns and fit with polynomials
  ava=fltarr(sz(2))
  for iy=0,sz(2)-1 do begin
   ava(iy)=(avg(image(ib2-5:ib2+5,iy)))
  endfor

  fita=poly_fit(xx,ava,3,afit)

   ixa=sz(1)-postpix-nxpix		; single amplifier

; subtract the fits to the overscan

  for iy=0,sz(2)-1 do begin
   image(ixa:ixa+nxpix-1,iy)=(image(ixa:ixa+nxpix-1,iy)-afit(iy))
  endfor

  endif   ; end of numamps section
endif	; end of subtraction section
; END OF OVERSCAN SUBTRACTION

; fix hot columns
if keyword_set (hot) then begin
  for ihot=0,n_elements(hotcol)-1 do begin
   for iy=0,ny-1 do begin
    image(hotcol(ihot)+prepix,iy) = $
     median(image(prepix+hotcol(ihot)-2:prepix+hotcol(ihot)+2,iy))
   endfor
  endfor
endif


;stop
; setting trim=1 will make a 1-amp size gain-corrected image out of a two amp
if keyword_set(trim) then begin
 ixb=sz(1)-postpix-nxpix		; left amplifier
 ixa=sz(1)-postpix-nxpix/2		; right amplifier
  for iy=0,sz(2)-1 do begin
   image(ixb:ixb+nxpix/2-1,iy)=(image(ixb:ixb+nxpix/2-1,iy)+bfit(iy))
   image(ixa:ixa+nxpix/2-1,iy)=(image(ixa:ixa+nxpix/2-1,iy)+bfit(iy))
  endfor
 imagetmp=image(prepix/2:prepix/2+nxpix+prepix/2+postpix/2-1,*)
 return,imagetmp
endif else begin
 imagetmp=image(prepix:prepix+nxpix-1,0:sz(2)-1,*)
 return,imagetmp
endelse

end
