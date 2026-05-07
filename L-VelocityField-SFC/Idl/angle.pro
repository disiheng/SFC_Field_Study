
PRO angle, a, b, n, rad

if n eq 1 then begin
  norma=sqrt(total(a^2)+1d-12)
  normb=sqrt(total(b^2)+1d-12)
  rad = total(a * b) / (norma * normb) 
endif else if n gt 1 then begin
  i=0L
  repeat begin
    norma=sqrt(total(a[*,i]^2)+1d-12)
    normb=sqrt(total(b[*,i]^2)+1d-12)
    rad[i] = total(a[*,i] * b[*,i]) / (norma * normb) 

    i=i+1L
  endrep until i ge n
endif

return
END
