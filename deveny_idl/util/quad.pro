
function quad, y, abc

; find x using the quadratic formula
; abc: [a, b, c]
; y: polynomial solution

  a = abc[0]
  b = abc[1]
  c = abc[2] - y

  xp = (-b + sqrt(b^2 - 4*a*c))/(2.*a)
  xm = (-b - sqrt(b^2 - 4*a*c))/(2.*a)

  return, [xp, xm]

end

