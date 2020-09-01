
pro dfitcurves, fwhm, filelist, fits=fits, focus=focus, $
  fwidth=fwidth, minfoc=minfoc, fnom=fnom

;
; Fit line/star focus curves.
;
  thresh = 2.0
  norder = 2
  if not(keyword_set(fnom)) then fnom=2.7
  sz = size(fwhm) & nc = sz(1) & nf = sz(2)
  fits = fltarr(norder+1,nc)
  minfoc = fltarr(nc)
  x = findgen(nf)
  focus=fltarr(nc) & fwidth=fltarr(nc)
  xfine = findgen((nf-1)*10.0+1.0)/10.0
  wts = make_array(nf,value=1.0)

  for i = 0,nc-1 do begin
    wts(*) = 1.0
    data = fwhm(i,*)
    out = where ((data lt 1.0) or (data gt 15.0),count)
    if count ge 3 then goto,skip
    if count ne 0 then begin
      wts(out)=0.0
      data(out) = 50.0
    endif
    fit = goodpoly(x,data,norder,thresh,yfit,newx,newy)
;   fit = poly_fit(x,data,norder,measure_errors=wts)
    fits(*,i) = fit
;  Curve minimum:
    fitfine=poly(xfine,fit) & minfine=min(fitfine) & focus(i)=xfine(!c)
    minfoc(i)=minfine
;  Nominal focus position:
    fpix = quad(fnom, reverse(reform(fit)))
    fwidth(i) = fpix(0)
    skip:
  endfor

; Plot up focus curves:
; dplotfocus, x, fwhm, fits, filelist

  return
end
