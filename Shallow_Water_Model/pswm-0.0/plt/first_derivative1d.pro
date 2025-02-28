PRO first_derivative1d, input, $
                  dx, $
                  dvar

;==================================================
; A routine which calculates the partial derivative 
; of a 2-D field using centered differencing 
; (one-sided on boundaries).
;
; INPUT
;    input = FLTARR(nx) = a 1-D field
;    
;    dx = grid spacing
;
; 
; OUTPUT
;    dvar = FLTARR(nx)
;                   
;==================================================

s = SIZE(input)
nx = s[1]

dvar = FINDGEN(nx) * 0.0  ; zero the array dvar

   FOR i = 1, nx - 2 DO BEGIN
      
      dvar[i] = ( 1./(2.*dx) ) * ( input[i+1] - input[i-1] )
      
   ENDFOR
      
      dvar[0] = ( 1./dx ) * ( input[1] - input[0] )
      dvar[nx-1] = ( 1./dx ) * ( input[nx-1] - input[nx-2] )
      
END
