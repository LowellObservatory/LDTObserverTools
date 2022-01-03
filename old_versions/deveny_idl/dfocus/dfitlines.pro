
pro dfitlines, spectrum, cnt, fwhm=fwhm, cplot=cplot

;
; Procedure to fit line profiles in a focus image
;

  if keyword_set(cplot) then begin
    !p.multi(1)=4
    !p.multi(2)=7
  endif
  nc = n_elements(cnt)
  fwhm = fltarr(nc)
  boxi = 17 & boxf=11
  xx = findgen(2*boxf+1)

  for j=0,nc-1 do begin
    if (cnt(j)-boxi) lt 0 or (cnt(j)+boxi) gt 2030 then goto, skip
    box = spectrum(cnt(j)-boxi:cnt(j)+boxi)
    ccnt = com1d(box)+cnt(j)-boxi
    if (ccnt-boxi) lt 0 or (ccnt+boxi) gt 2030 then goto, skip
    box = spectrum(ccnt-boxf:ccnt+boxf)
    tmp = gaussfit(xx, box, a)
    fwhm(j) = a(2)*1.177*2.0
    if keyword_set(cplot) then begin
      plot, xx, box, psym=-8, title='Line center: ' + strtrim(cnt(j), 2), $
        color=50
      oplot, xx, tmp, color=128
    endif
    skip:
  endfor

  !p.multi=0
  return
end
