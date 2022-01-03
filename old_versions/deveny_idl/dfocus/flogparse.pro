
pro flogparse, flog=flog, flist=flist

  if (flog eq 'last') then begin
    spawn, 'ls -1 deveny_focus* > focus.files'
    spawn, 'tail -1 focus.files', flog
  endif
  id = strmid(flog, 13, 15)
  flist = 'dfocus.'+id
  openr, lun, flog, /get_lun
  openw, lun1, flist, /get_lun

  line = ''
; Read header:
  readf, lun, line 
  while not(eof(lun)) do begin
    readf, lun, line 
    fr1 = strmid(line, 2, 18)
    printf, lun1, fr1
  endwhile

  free_lun, lun
  free_lun, lun1
 
return
end
