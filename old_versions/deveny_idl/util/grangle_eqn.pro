function grangle_eqn, theta

  common grvar, gpmm, wavelen
; DeVeny optical angles:
  camcol=rad(55.00)
  coll=rad(10.00)

  gx = (sin((coll+theta)) + sin(coll+theta-camcol))*1.e7/gpmm - wavelen

  return, gx
end
