pro dextract,spectrum,traces,nspix,swext=swext,spectra=spectra,npixavg=npixavg

if not (keyword_set(swext)) then swext=2

; Object spectral extraction:
;  Options:
;   swext = 1: extract point source spectra by averaging over window
;   swext = 2: extract full spatial orders with spatial interpolation
;   swext = 3: extract full orders without interpolation 
;   swext = 4: extract spectra by averaging over window
;
;  Output:
;   spectra: 2 or 3-d array of spectra of individual orders

sz = size(traces)
nx = sz(1)
norders = sz(2)
if sz(0) eq 1 then norders = 1

if swext eq 0 then goto,out

if swext eq 4 then begin
 if not(keyword_set(npixavg)) then npixavg=nspix
 spectra = fltarr(nx,norders)
 for io=0,norders-1 do begin
  spectra(*,io)=specavg(spectrum,traces(*,io),npixavg)
 endfor
endif

if swext eq 1 then begin
 if not(keyword_set(npixavg)) then npixavg=19
 spectra = fltarr(nx,norders)
 for io=0,norders-1 do begin
  spectra(*,io)=dspecfit(spectrum,traces(*,io),bwidth=nspix,extwidth=npixavg)
 endfor
endif

if swext eq 2 then begin
 spectra = fltarr(nx,nspix,norders)
 fnspix = findgen(nspix)-nspix/2
  for io = 0,norders-1 do begin
   for ix = 0,nx-1 do begin
; Interpolation:
    xt = fix(traces(ix,io))+fnspix
    ut = traces(ix,io)+fnspix
    vector = spectrum(ix,xt)
    tvector = interpol(vector,xt,ut)
    spectra(ix,*,io) = tvector
   endfor
  endfor
endif

if swext eq 3 then begin
 spectra = fltarr(nx,nspix,norders)
 inspix = indgen(nspix)-nspix/2
  for io = 0,norders-1 do begin
   for ix = 0,nx-1 do begin
; Interpolation:
    xt = fix(traces(ix,io))+inspix
    spectra(ix,*,io) = spectrum(ix,xt)
   endfor
  endfor
endif

out:

return

end
