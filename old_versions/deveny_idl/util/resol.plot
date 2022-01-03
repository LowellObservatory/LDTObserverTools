
plot,wave,resol(*,0),xra=[3000,11000],yra=[0,8000], $
 xtitle='Wavelength('+'!3'+string("305b)+')',ytitle='R ('+'!7k'+'!3/'+'!7Dk'+' !3(2" slit))', $
 title='!3WhiteCam Spectral Resolution vs. Wavelength and Grating; Slit Width 1"', $
 ticklen=1.,xsty=1
for i=1,2 do begin
 oplot,wave,resol(*,i)
endfor
for i=0,2 do begin
 xyouts,1.06e4,resol(700,i),strtrim(gpmm(i),1)
endfor
end
