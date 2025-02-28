;==================================================

PRO integral1d, input, $
                x_0, x_n, $
                dx, $
                var_integrated
			
;==================================================
; A routine which calculates the definite integral
; of a 1_D variable using the trapezoidal rule.
;
; INPUT
;    input = FLTARR(n_1) = a 1-D field
;    
;    dx = grid spacing used in the integral sum
;
;    x_0 = index value of lower bound 
;
;    x_n = index value of upper bound
;
; 
; OUTPUT
;    var_integrated = a constant value

;==================================================
  
var_integrated   = dx * ( input[x_0] * 0.5 + $
                       TOTAL(input[x_0+1:x_n-1]) + $
                       input[x_n] * 0.5 )
         
END






