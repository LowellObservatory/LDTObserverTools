
pro rot1d, x, y, phi, xnew, ynew

; this procedure performs a rotation on a coordinate (x,y)
; angle in degrees

rphi = rad(phi) 
cphi = cos(rphi) & sphi=sin(rphi)

xnew =  cphi*x + sphi*y
ynew = -sphi*x + cphi*y

return

end
