PRO first_derivative2d, input, $
                  dx, dy, $
                  dim, dvar

;==================================================
; A routine which calculates the partial derivative 
; of a 2-D field using centered differencing 
; (one-sided on boundaries).
;
; INPUT
;    input = FLTARR(nx, ny) = a 2-D field
;    
;    dx, dy = grid spacing
;
;    dim = dimension over which differentiation is
;          applied (1 or 2)
; 
; OUTPUT
;    dvar = FLTARR(nx, ny)
;                   
;==================================================

s = SIZE(input)
nx = s[1]
ny = s[2]

dvar = FINDGEN(nx, ny) * 0.0  ; zero the array dvar


IF (dim EQ 1) THEN BEGIN

FOR j = 0, ny - 1 DO BEGIN
   FOR i = 1, nx - 2 DO BEGIN
      
      dvar[i,j] = ( 1./(2.*dx) ) * ( input[i+1,j] - input[i-1,j] )
      
   ENDFOR
      
      dvar[0,j] = ( 1./dx ) * ( input[1,j] - input[0,j] )
      dvar[nx-1,j] = ( 1./dx ) * ( input[nx-1,j] - input[nx-2,j] )
      
ENDFOR      ; end j loop

ENDIF ELSE BEGIN

IF (dim EQ 2) THEN BEGIN

FOR i = 0, nx - 1 DO BEGIN
   FOR j = 1, ny - 2 DO BEGIN
      
      dvar[i,j] = ( 1./(2.*dy) ) * ( input[i,j+1] - input[i,j-1] )
      
   ENDFOR
      
      dvar[i,0] = ( 1./dy ) * ( input[i,1] - input[i,0] )
      dvar[i,ny-1] = ( 1./dy ) * ( input[i,ny-1] - input[i,ny-2] )
      
ENDFOR      ; end i loop

ENDIF 

ENDELSE

END
