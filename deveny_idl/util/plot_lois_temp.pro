
file = ''
print,''
print,' Enter file name:'
read, file
openr, lun, file, /get_lun
spawn, 'wc ' + file, text
nt = fix(strmid(text, 5, 3))
nt = nt(0)
time = fltarr(nt)
temp1 = fltarr(nt)

t0=0.
line = ''
for i=0, nt-1 do begin
  readf, lun, line
  ut = strmid(line, 0, 8)
  ti = raparse(ut)*12./!pi
  if ((t0 gt 23) and (ti lt 23.0)) then ti=ti+24.0
  tmp = float(strmid(line, 26, 6))
  time(i) = ti
  temp1(i) = tmp
  t0 = ti
endfor

!p.multi(2) = 3
setplot, 'ps', filename = file + '.ps', /port
  plot, time, temp1, ticklen=1, title='DeVeny Temp log ' + file, $
    xtitle='Hours', ytitle='T (C)', psym=-8, $
    yra=[-130, 30], ysty=1, xra=[time(0)-1, time(0)+5], xsty=1

  tgrad = (temp1(1:*) - temp1(0:nt-2))/(time(1:*) - time(0:nt-2))/60.
  tacc = (tgrad(1:*) - tgrad(0:nt-2))/(time(1:*) - time(0:nt-2))/60.

  plot, time, tgrad, ticklen=1, $
    title='Temperature Gradient', $
    xtitle='Hours', ytitle='dT/dt (deg C/min)', $
    yra=[-5, 5], ysty=1, xra=[time(0)-1, time(0)+5], xsty=1

  plot, time, tacc, ticklen=1, $
    title='Temp rate', $
    xtitle='Hours', ytitle='dG/dt (deg C/min)^2', $
    yra=[-1, 1], ysty=1, xra=[time(0)-1, time(0)+5], xsty=1

!p.multi=0
setplot,'x'
free_lun, lun
end
 
