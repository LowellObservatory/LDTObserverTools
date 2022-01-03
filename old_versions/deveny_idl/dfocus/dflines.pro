
function dflines, image, thresh=thresh, mark=mark, fwhm=fwhm, title=title

; Procedure to automatically find and centroid lines in 1-row-image:

sz=size(image) & ny=sz(1) 
avgj = image

if not keyword_set(mark) then mark=0
if not keyword_set(thresh) then thresh=20.
if not keyword_set(title) then title=''
peaks=fltarr(1)
fwhm=fltarr(1)

; Create background from median value of the image:
bkgd=median(image)
print,'  Background level:',bkgd
print,''

; Step through the cut and identify peaks:
fmax=11		; minimum separation
fwin=15
findmax=50
j0=0

if (keyword_set (mark)) then begin
  centers=[0]
  goto,reselect
endif

; !p.multi(1)=5
; !p.multi(2)=9
for j=0,ny-1 do begin
  if j gt (ny-fmax) then goto,nextpeak
	if avgj(j) gt (bkgd+thresh) then begin
	 j1=j
	  if abs (j1-j0) lt fmax then goto, nextpeak
		for jf=0,findmax do begin
			itmp0=avgj(jf+j)
			itmp1=avgj(jf+j+1)
			if itmp1 lt itmp0 then begin 
				icntr=jf+j
				goto, getout
			endif
		endfor
			getout:
			begin
			if (icntr lt fmax/2 or icntr gt (ny-fmax/2-1)) then goto,nextpeak
 			xx=findgen(fwin)+float(icntr-fwin/2)
 			temp=avgj((icntr)-fwin/2:(icntr)+fwin/2)
			temp = median(temp,3)
 			tempfit=gaussfit(xx,temp,aa,nterms=5)
			;plot,xx,temp,psym=1
			;oplot,xx,tempfit
			;plots,aa(1),aa(0),psym=4
 			center=aa(1)
 			fw=aa(2)*1.177*2.
                        pmax=aa(0)
			if fw gt 1.00 then begin
			  peaks=[peaks,center]
			  fwhm=[fwhm,fw]
		        endif
			end
         j0=jf+j
	endif
 nextpeak:
endfor

centers=peaks(1:*)
fwhm=fwhm(1:*)

reselect:
badcnt = where((centers lt 0) or (centers gt 2100), complement=cc, /null)
centers=centers(cc)
!p.multi=0
szc=size(centers)
print,' Number of lines:',szc
print,''
; print,' Enter ymax for plot'
; read,ymax

; Plot the lines and mark the ones that were found:
;plot, avgj, xra=[0,ny+2], xsty=1, $
;  title='
plot,avgj,xra=[0,ny+2],xsty=1,yra=[0,max(avgj)+0.2*max(avgj)], $
  title=title, xtitle='CCD column', ytitle='I (DN)'
for id=0,szc(1)-1 do begin
 plots,centers(id),avgj(centers(id))+20,psym=1
 xyouts,centers(id),avgj(centers(id))+30,strtrim(centers(id),2), $
 /data,orientation=0.
endfor

goto,out

print, 'Select orders by clicking left on the plot;'
print, 'De-Select orders by clicking middle;'
print, ' Exit and replot with right'

!err=0
cr = string("15b)       ;this codes a newline
form="($,'x=',F10.3,', y=',F10.3,a)"

while !err ne 4 do begin
 cursor,x,y,3,/data,/down
;print,form = form, x, y, cr

 if (!err eq 1 and mark eq 0) then begin
  xx=findgen(fmax)+fix(x)-fmax/2
  temp=avgj(fix(x)-fmax/2:fix(x)+fmax/2)
  temp = median(temp,3)
  tempfit=gaussfit(xx,temp,aa,nterms=5)
  centers=[centers,aa(1)] 
  plots,aa(1),avgj(aa(1))+10,psym=1
  xyouts,aa(1)+10,avgj(aa(1))+10,string(aa(1)), $
    /data,orientation=90.
  print,' New center at', aa(1),' value = ',avgj(aa(1))
  print,''

 endif else if (!err eq 1 and mark) then begin
  centers=[centers,x]
  plots,x,y+10,psym=1
  xyouts,x,y+10,string(x),/data,orientation=90
  print,' New center at', x
  print,''

 endif else if (!err eq 2) then begin
  out=where(abs(centers-x) lt 4.,count)
  print,' Deleted center number', out, centers(out)
  print,''
  if count ne 0 then centers(out)=0.
 endif
; print,form = form, x,y,cr

endwhile

; sort and replot
nhits=where(centers eq 0,count)
if count ne 0 then begin
 sortc=sort(centers)
 centers=centers(sortc(count:*))
endif

szc=size(centers)
; Plot the lines and mark the ones that were found:
plot,avgj,xra=[0,ny+2],xsty=1,yra=[0,max(avgj)+0.2*max(avgj)], $
  title=title, xtitle='CCD column', ytitle='I (DN)'
for id=0,szc(1)-1 do begin
 plots,centers(id),avgj(centers(id))+10,psym=1
 xyouts,centers(id)+5,avgj(centers(id))+10,strtrim(centers(id),0), $
 /data,orientation=0.
endfor
 
; Reloop option:
yy=''
print,' Enter y to select/deselect again'
read,yy
if yy eq 'y' then goto,reselect

out:

return,centers(sort(centers))

end

