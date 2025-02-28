;==================================================

PRO integral2d, input, $
                x_0, x_n, dim, $
                n_1, n_2, dx, $
                var_integrated
			
;==================================================
; A routine which calculates the definite integral
; of a 2_D variable using the trapezoidal rule.
;
; INPUT
;    input = FLTARR(n_1, n_2) = a 2-D field
;    
;    dx = grid spacing used in the integral sum
;
;    x_0 = index value of lower bound 
;
;    x_n = index value of upper bound
;
;    dim = dimension over which to integrate (1, 2, or 3)
; 
; OUTPUT
;    var_integrated = FLTARR(n_2, n_3), if dim = 1
;                   = FLTARR(n_1, n_3), if dim = 2
;                   = FLTARR(n_1, n_2), if dim = 3
;
;==================================================

IF (dim EQ 1) THEN BEGIN

   var_integrated = FLTARR(n_2)

   FOR j = 0, n_2 - 1 DO BEGIN
     

      var_integrated[j, k] = dx * ( input[x_0, j, k] * 0.5 + $
                       TOTAL(input[x_0+1:x_n-1, j, k], dim) + $
                       input[x_n, j, k] * 0.5 )
     
     
   ENDFOR
   
ENDIF ELSE BEGIN

IF (dim EQ 2) THEN BEGIN

   var_integrated = FLTARR(n_1)
   
   FOR i = 0, n_1 - 1 DO BEGIN
     
      
      var_integrated[i] = dx * ( input[i, x_0] * 0.5 + $
                       TOTAL(input[i, x_0+1:x_n-1], dim) + $
                       input[i, x_n] * 0.5 )
      
     
   ENDFOR
   
ENDIF 
ENDELSE

END






