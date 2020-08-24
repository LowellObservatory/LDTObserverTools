
function com2d, box, signal=tbox

; returns the center-of-mass of box
; background subtraction from annular average
; all less than zero are set to 0


sz=size(box)

nx=sz(1)
ny=sz(2)

sky=median([box(*,0),box(*,ny-1),transpose(box(0,1:ny-2)),transpose(box(nx-1,1:ny-2))])

box2=box-sky
ltz=where(box2 lt 0.)
box2(ltz)=0.

xx=findgen(nx) & yy=findgen(ny)

cx=0. & cy=0.

tbox = total(box2)
for j=0,ny-1 do cx=cx+xx*box2(*,j)
for i=0,nx-1 do cy=cy+yy*box2(i,*)

dx=total(cx)/tbox
dy=total(cy)/tbox

return, [dx,dy,tbox]
end

