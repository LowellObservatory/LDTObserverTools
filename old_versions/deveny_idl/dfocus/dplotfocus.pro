
pro dplotfocus, x, fwhm, fits, filelist, fnom=fnom

  common curves, fx, centers, fpos, fwid, fwmin

  sz = size(fwhm) & nc = sz(1)
  !p.multi(1)=6
  !p.multi(2)=nc/6+1
  if (nc/4+1) gt 5 then !p.multi(2)=5
  for i=0,nc-1 do begin

    plot, fx, fwhm(i,*), psym=4, yra=[0,max(fwhm(i,*))], $
      ticklen=1,ysty=1, ymargin=[4, 4], $
      xtitle='Collimator Position (mm)', $
      ytitle='FWHM', $
      title='LC: '+strtrim(centers(i),1)+'  Fnom: '+strmid(fnom,6,5)+' pixels'
    oplot, fx, poly(x,fits(*,i))
    plots, [fpos(i), fpos(i)], [0, fwmin(i)], /data, color=90, thick=3
    plots, [fwid(i), fwid(i)], [0, fnom], /data, color=60, thick=3

  endfor

  !p.multi=0

  return
end
