
function specavg, spectrum, trace, wsize

; extract an average spectrum along trace of size wsize

sz=size(spectrum) & nx = sz(1)
speca=fltarr(nx)

for i=0,nx-1 do begin
  speca(i)=avg(spectrum(i, trace(i)-wsize/2:trace(i)+wsize/2))
endfor

return, speca

end

