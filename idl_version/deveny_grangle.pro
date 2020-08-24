
pro deveny_grangle

; Calculate the DeVeny grating angle tilt given the following variables:
;	gpmm	Grating ruling (grooves/mm)
;	wavelen	Central wavelength in imaged spectral format (Angstroms)
  common grvar, gpmm, wavelen

; Mechanical offset:
  tgoffset = 0.0

  print,' Enter grating resolution (g/mm):'
  read, gpmm
  print,' Enter central wavelength (A):'
  read, wavelen

; Call IDL function 'newton' to solve the grating equation:
  theta=rad(20.)
  grangle = newton(theta, 'grangle_eqn')
  grangle = deg(grangle) 
  amag = deveny_amag(grangle)

  print,''
  print,' Grating: ', gpmm, ' g/mm'
  print,' Central wavelength: ', wavelen, ' A'
  print,' DeVeny grating tilt = ', strtrim(grangle+tgoffset, 1),' deg', format='(a,f6.2,a)'
  print,' Slit demagnification (pixels/arcsec, 0.34 arcsec/pixel):  ', $
    2.94*amag, format='(a,f6.2,a)'
  print,''

  return
end


function grangle_eqn, theta

  common grvar, gpmm, wavelen
; DeVeny optical angles:
  camcol=rad(55.00)
  coll=rad(10.00)

  gx = (sin((coll+theta)) + sin(coll+theta-camcol))*1.e7/gpmm - wavelen

  return, gx
end


function deveny_amag, grangle

  collang = 10.0
  camcollang = 55.0
  alpha = rad(grangle + collang)
  mbeta = rad(camcollang - deg(alpha))

  return,(cos(alpha)/cos(mbeta))

end

