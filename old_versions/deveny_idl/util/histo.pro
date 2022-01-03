
pro histo, data, n, hrange=hrange, hmedian=hmedian, havg=havg, hmax=mx, $
 title=title, subtitle=subtitle, xtitle=xtitle, plot=plot, $
 printval=printval, overplot=overplot, charsize=charsize, yrange=yrange, $
 xhist=xax, hist=hdata, top=top, sigrej=sigrej, xrange=xrange,  $
 ndigits=ndigits, norm=norm, charthick=charthick, thick=thick, offset=offset, $
 noerase=noerase, norej=norej

; This routine will calculate and plot histograms for floating point or
;  integer data.

if n_params() eq 0 then begin
 print,''
 print,' Usage:
 print,' histo,data,n, hrange=hrange, hmedian=hmedian, havg=havg, hmax=mx, $'
 print,' title=title, subtitle=subtitle, xtitle=xtitle, plot=plot, $'
 print,' printval=printval, overplot=overplot, charsize=charsize, yrange=yrange, $'
 print,' xhist=xax, hist=hdata, top=top, sigrej=sigrej, xrange=xrange, ndigits=ndigits'
 print,''
 print,'where:'
 print,''
 print,' n = binsize 
 print,' hrange = [minimum,maximum] value for histogram.  
 print,'  default: minmax of supplied data
 print,' hmedian = value: return median value of the data 
 print,' havg = value: return average value of the data
 print,' hmax = value: return data value of the histogram maximum
 print,''
 print,' Set /plot for plot, with optional titles'
 print,' Set /printval to print out median,avg,max on the plot
 print,' Set /overplot to overplot input data on existing plot
 return
endif

if keyword_set(subtitle) eq 0 then subtitle=''
if keyword_set(xtitle) eq 0 then xtitle=''
if keyword_set(title) eq 0 then title=''
if not (keyword_set(charsize)) then charsize = 0.75
if not (keyword_set(charthick)) then charthick = 1.0
if not (keyword_set(thick)) then thick = 1.0
if not (keyword_set(top)) then top = 1.0
if not (keyword_set(sigrej)) then sigrej = 6.0
if not (keyword_set(ndigits)) then ndigits = 2
if not (keyword_set(offset)) then offset = 0
if not (keyword_set(hrange)) then begin
  hran = [min(data), max(data)]
  x1 = hran[0] & x2 = hran[1]
  rdata = float(data)
endif
if keyword_set(hrange) then begin
  x1 = hrange[0] & x2 = hrange[1]
  hxd = where((data ge x1) and (data le x2))
  rdata = float(data(hxd))
endif

hdata=histogram(data, binsize=n, min=(x1), max=(x2))
xax=findgen((x2-x1)/n+1)*n+x1+float(n)/2.
nax = n_elements(xax)

; Return median value of all data, and average out to +/- sigrej-sigma

ndata = float(n_elements(data))
nr = float(n_elements(rdata))
if (keyword_set(norm)) then hdata=hdata/ndata
hmax=max(hdata) & mx=xax(!c)
htop = where(hdata gt (1.-top)*hmax)
hmedian = median(rdata)
hmediantop = median(xax(htop))
dmax = max(rdata)
dmin = min(rdata)
if not(keyword_set(norej)) then begin
  robomean, rdata, sigrej, 0.1, avdata, stddata
  xplot = [avdata-sigrej*stddata,avdata+sigrej*stddata]
endif else begin
  result = moment(rdata, mean=avdata)
  result = moment(rdata, sdev=stddata)
  xplot = [dmin,dmax]
endelse

nd = strtrim(ndigits,1)
print,''
print,title
print,' Median:',string(hmedian,format='(f10.'+nd+')')
;print,' Median, top '+strtrim(top,1)+':',string(hmediantop,format='(f10.'+nd+')')
print,' Mean: ',string(avdata,format='(f10.'+nd+')'),' +/- ',string(stddata,format='(f10.'+nd+')')
print,' Max bin:',string(mx,format='(f10.'+nd+')')
print,''

if keyword_set (plot) then begin
if not(keyword_set (yrange)) then yrange=[0,max(hdata)]
if not(keyword_set (xrange)) then xrange=xplot
plot, xax, hdata, psym=10,xrange=xrange, $
  title=title,subtitle=subtitle,xtitle=xtitle,ytitle='Count', $
  yrange=yrange,charsize=charsize,noerase=noerase, xsty=1
plots,[xax(0),xax(0)],[0,hdata(0)]
plots,[xax(nax-1),xax(nax-1)],[0,hdata(nax-2)]
endif

; Overplot option (includes printing values on plot):
if keyword_set(overplot) then begin
 oplot,xax,hdata,psym=10,linestyle=4,thick=1.5
 plots,[xax(0),xax(0)],[0,hdata(0)],linestyle=4,thick=1.5
 plots,[xax(nax-1),xax(nax-1)],[0,hdata(nax-1)],linestyle=4,thick=1.5

 if (keyword_set(printval)) then begin
; plots,[x2-0.2*(x2-x1),x2-0.05*(x2-x1)],[hmax,hmax],linestyle=4,thick=1.5
  xyouts,x2-0.3*(x2-x1),hmax-0.1*hmax, $
	'Median: '+string(hmedian,format='(f12.'+nd+')'), $
	size=0.75,/data
  xyouts,x2-0.3*(x2-x1),hmax-0.2*hmax, $
	'Mean: '+string(avdata,format='(f12.'+nd+')')+ $
	' +/- '+string(stddata,format='(f12.'+nd+')'), $
   	size=0.75,/data
  xyouts,x2-0.3*(x2-x1),hmax-0.3*hmax, $
	'Max value: '+string(dmax,format='(f12.'+nd+')'), $
   	size=0.75,/data
  xyouts,x2-0.3*(x2-x1),hmax-0.4*hmax, $
	'Max bin: '+string(mx,format='(f12.'+nd+')'), $
   	size=0.75,/data
  xyouts,x2-0.3*(x2-x1),hmax-0.4*hmax, $
	'Min value: '+string(dmin,format='(f12.'+nd+')'), $
   	size=0.75,/data
 endif

endif 

!p.font=0
; Print values on plot.
if (keyword_set(printval) and not (keyword_set(overplot))) then begin
 xyouts,xplot(1)-0.35*(xplot(1)-xplot(0))+offset,hmax-0.1*hmax+offset, $
	'Median: '+string(hmedian,format='(f10.'+nd+')'), $
	size=charsize,/data,charthick=charthick
 xyouts,xplot(1)-0.35*(xplot(1)-xplot(0))+offset,hmax-0.17*hmax+offset, $
	'Mean: ' +string(avdata,format='(f10.'+nd+')')+' +/- '+ $
	string(stddata,format='(f10.'+nd+')'), $
	size=charsize,/data,charthick=charthick
 xyouts,xplot(1)-0.35*(xplot(1)-xplot(0))+offset,hmax-0.24*hmax+offset, $
	'Max value: '+string(dmax,format='(f10.'+nd+')'), $  
	size=charsize,/data,charthick=charthick
 xyouts,xplot(1)-0.35*(xplot(1)-xplot(0))+offset,hmax-0.31*hmax+offset, $
	'Max bin: '+string(mx,format='(f10.'+nd+')'), $  
	size=charsize,/data,charthick=charthick
 xyouts,xplot(1)-0.35*(xplot(1)-xplot(0))+offset,hmax-0.38*hmax+offset, $
	'Min value: '+string(dmin,format='(f10.'+nd+')'), $  
	size=charsize,/data,charthick=charthick
 xyouts,xplot(1)-0.35*(xplot(1)-xplot(0))+offset,hmax-0.45*hmax+offset, $
	'Fraction of points: '+string(nr/ndata,format='(f10.'+nd+')'), $  
	size=charsize,/data,charthick=charthick

endif

return
end

