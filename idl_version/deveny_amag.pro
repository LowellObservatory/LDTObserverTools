function deveny_amag, grangle

  collang = 10.0
  camcollang = 55.0
  alpha = rad(grangle + collang)
  mbeta = rad(camcollang - deg(alpha))

  return,(cos(alpha)/cos(mbeta))

end

