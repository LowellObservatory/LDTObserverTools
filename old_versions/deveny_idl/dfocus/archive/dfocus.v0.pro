
; Main level routine to fit a set of DeVeny focus spectral images 

pro dfocus, nf, df, f0, file=file, thresh=thresh
 
  common curves, fx, centers, fpos, fwid

  if not(keyword_set(file)) then file='filelist'
  if not(keyword_set(thresh)) then thresh=300.0
  files = readstr1(file)
  nfiles = n_elements(files)
  mfile = files(nfiles/2)
  x = findgen(nfiles)
  fx = findgen(nfiles)*df + f0
  loadct, 5
  setplot,'ps',filename=file+'.ps',/port, /color, /land

; Examine middle image:
  swext=4
  win = 11
  print,''
  print,' Processing object image'+mfile+'...'
  spectrum = dvtrim(mfile,/sub,postpix=postpix)
  sz = size(spectrum)
  print,sz
  traces = make_array(sz(1),value=sz(2)/2,/float)
  dextract,spectrum,traces,win,spectra=mspectra,swext=swext
; Find the lines:
  centers=dflines(mspectra, fwhm=fwhm, thresh=thresh)  
  nc = n_elements(centers)
  fw = fltarr(nc, nfiles)

; Run through files:
  for i=0,nfiles-1 do begin
    print,''
    print,' Processing arc image'+files(i)+'...'
    spectrum = dvtrim(files(i),/sub,postpix=postpix)
    print,''
    print,'   Extracting spectra from image '+files(i)+'...'
    dextract,spectrum,traces,win,spectra=spectra,swext=swext

; Find fwhm of lines:
    dfitlines,spectra,centers,fwhm=fwhm
    fw(*,i) = fwhm
  endfor

; Fit the lines:
  dfitcurves, fw, file, fits=fits, focus=focus, fwidth=fwidth
  fpos = focus*df + f0
  fwid = fwidth*df + f0
  dplotfocus, x, fw, fits, file
  plot, centers, fpos, /ynoz, psym=-5
  plot, centers, fwid, /ynoz, psym=-5
  setplot,'x'
  loadct, 0

return
end

