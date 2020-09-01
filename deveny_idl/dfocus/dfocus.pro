
; Main level routine to fit a set of DeVeny focus spectral images 

pro dfocus, flog=flog, thresh=thresh, $
      grating=grating, grangle=grangle, fnom=fnom
 
  common curves, fx, centers, fpos, fwid, fwmin

  if not(keyword_set(thresh)) then thresh=100.0
; if not(keyword_set(fnom)) then fnom=2.7

; Parse the log file to obtain file list;
;   if log is unspecified, process the last sequence:
  if not(keyword_set(flog)) then flog='last'
  flogparse, flog=flog, flist=flist
  files = readstr1(flist)
  nfiles = n_elements(files)

; Query the setup and initialize the analysis:
  file0 = '../' + files(0)
  hdr0 = headfits(file0)
  slitasec = sxpar(hdr0, 'SLITASEC')
  grating = sxpar(hdr0, 'GRATING')
  grangle = sxpar(hdr0, 'GRANGLE')
  lampcal = sxpar(hdr0, 'LAMPCAL')
  filter = sxpar(hdr0, 'FILTREAR')
  fnom = 2.94*slitasec*deveny_amag(grangle)

  f0 = sxpar(hdr0, 'COLLFOC')
  file1 = '../' + files(nfiles-1)
  hdr1 = headfits(file1)
  f1 = sxpar(hdr1, 'COLLFOC')
  df = (f1 - f0)/(nfiles-1)
  x = findgen(nfiles)
  fx = findgen(nfiles)*df + f0

  loadct, 5
; if mode eq 'manual' then goto, nocopy
  setplot,'ps',filename=flist+'.ps',/port, /color, /land
  nocopy:
; Examine middle image:
  mfile = '../' + files(nfiles/2)
  swext=4
  win = 11
  print,''
  print,' Processing object image'+mfile+'...'
  spectrum = dvtrim(mfile,/sub,postpix=postpix)
  sz = size(spectrum)
  print,sz
  traces = make_array(sz(1),value=sz(2)/2,/float)
  dextract, spectrum, traces, win, spectra=mspectra, swext=swext
; Find the lines:
  centers=dflines(mspectra, fwhm=fwhm, thresh=thresh, $
    title = mfile + '   Grating:  ' + grating + '   GRANGLE:  ' + $
    strtrim(grangle,2) + '   ' + lampcal)
  nc = n_elements(centers)
  fw = fltarr(nc, nfiles)

; Run through files:
  for i=0,nfiles-1 do begin
    print,' Processing arc image'+files(i)+'...'
    spectrum = dvtrim('../'+files(i),/sub,postpix=postpix)
    print,''
    print,'   Extracting spectra from image '+files(i)+'...'
    dextract, spectrum, traces, win, spectra=spectra, swext=swext

; Find fwhm of lines:
    dfitlines, spectra, centers, fwhm=fwhm
    fw(*,i) = fwhm
  endfor

; Fit the lines:
  dfitcurves, fw, flist, fits=fits, focus=focus, $
    fwidth=fwidth, minfoc=minfoc, fnom=fnom
  fpos = focus*df + f0
  fwid = fwidth*df + f0
  fwmin = minfoc
  mfpos = median(fpos)
  mfwid = median(fwid)
  dplotfocus, x, fw, fits, flist, fnom=fnom
  plot, centers, fpos, xra=[0,2050], yra=[f0, f1], xsty=1, $
    /ynoz, ticklen=1, psym=5, $
    title='Minimum focus position vs. line position, median = ' + $
      strmid(mfpos, 3, 8), xtitle='CCD column', ytitle='Focus (mm)'
  plot, centers, fwid, xra=[0,2050], yra=[f0, f1], xsty=1, $
    /ynoz, ticklen=1, psym=5, $
    title='Optimal focus position vs. line position, median = ' + $
      strmid(mfwid, 3, 8), xtitle='CCD column', ytitle='Optimal Focus (mm)', $
      subtitle='Grating: ' + grating + $
      '   Slit width:' + strmid(slitasec, 4, 6) + ' arcsec' + $
      '    Nominal line width:' + strmid(fnom,4,7) + ' pixels'
  plots, [0, 2050], [mfwid, mfwid], color=60, thick=3, /data 
  setplot,'x'

  window, 0, xsize=750, ysize=450, xpos=50, ypos=725
  centers=dflines(mspectra, fwhm=fwhm, thresh=thresh, $
    title = mfile + '   Grating:  ' + grating + '   GRANGLE:  ' + $
    strtrim(grangle,2) + '   ' + lampcal)
  window, 1, xsize=1500, ysize=650, xpos=50, ypos=50
  dplotfocus, x, fw, fits, flist, fnom=fnom
  window, 2, xsize=750, ysize=450, xpos=805, ypos=725
  plot, centers, fwid, xra=[0,2050], yra=[f0, f1], xsty=1, $
    /ynoz, ticklen=1, psym=5, $
    title='Optimal focus position vs. line position, median = ' + $
      strmid(mfwid, 3, 8), xtitle='CCD column', ytitle='Optimal Focus (mm)', $
      subtitle='Grating: ' + grating + $
      '   Slit width:' + strmid(slitasec, 4, 6) + ' arcsec' + $
      '    Nominal line width:' + strmid(fnom,4,7) + ' pixels'
  plots, [0, 2050], [mfwid, mfwid], color=60, thick=3, /data 
  loadct, 0

return
end

