
function com1d, line

; returns the center-of-mass of line

  ltz=where(line lt 0,count)
  if count ne 0 then line(ltz)=0.
 
  sz=size(line)
  ny=sz(1)
  yy=findgen(ny)

  cx=0. & cy=0.
  cy=yy*line
  dy=total(cy)/total(line)

return, dy
end

