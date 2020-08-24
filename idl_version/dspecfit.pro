
function dspecfit,spec,trace,bwidth=bwidth,raw=raw,extwidth=extwidth

; fit the background to a spectrum and subtract it

tval=fix(trace)

if not (keyword_set (bwidth)) then bwidth=101 ; 
if not (keyword_set (extwidth)) then extwidth=19 

nfit=bwidth & win=[bwidth/2-extwidth/2,bwidth/2+extwidth/2] ;  w=0 over win
sz=size(spec)
fits=fltarr(nfit,sz(1)) & datas=fits & subs=fits

for i=0,sz(1)-1 do begin

ww=fltarr(nfit) & ww(*)=1. & ww(win(0):win(1))=0.
data=spec(i,tval(i)-nfit/2:tval(i)+nfit/2)
coef=polyfitw(findgen(nfit),data,ww,1)	; go with linear fit

fit=poly(findgen(nfit),coef)
fits(*,i)=fit
datas(*,i)=data
subs(*,i)=data-fit

;plot,data & oplot,fit
;stop

endfor

gplott=fltarr(sz(1))

for i=0,sz(1)-1 do gplott(i)=total(subs(win(0):win(1),i))

return,gplott

end

