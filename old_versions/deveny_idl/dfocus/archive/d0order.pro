
; Main level routine to fit a set of DeVeny focus spectral images 

pro d0order, cent, file=file
 
  if not(keyword_set(file)) then file='filelist'
  files = readstr1(file)
  nfiles = n_elements(files)
  xi=cent[0] & yi=cent[1]
  bi=100 & ci=20
  loadct, 22
  setplot,'ps',filename='thruslit.'+file+'.ps',/port, /color
  !p.multi(1) = 3
  !p.multi(2) = 5

; Run through files:
  for i=0,nfiles-1 do begin
    print,''
    print,' Processing image'+files(i)+'...'
    spectrum = readfits(files(i))
    box = spectrum(xi-bi:xi+bi, yi-bi:yi+bi)
    cb = com2d(box)
    cx = xi-bi+cb[0]
    cy = yi-bi+cb[1]
    cbox = spectrum(cx-ci:cx+ci, cy-ci:cy+ci)
    print,' image: ', files(i), cx, cy
    image_cont, cbox
  endfor

!p.multi=0
setplot, 'x'

return
end

