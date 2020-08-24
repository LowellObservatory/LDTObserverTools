
function readstr1,file

; function to read in 1 column strings

close,/all
openr,99,file
data = strarr(1000000)
n1 = ''

i=long(0)
while not eof(99) do begin

readf,99,n1

data(i) = [n1]
i=i+1

endwhile

close,99
return,data(0:i-1)

end

