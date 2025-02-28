;==================================================

FUNCTION azimuthalmean2d, input, $
                        perturbation = perturbation
			
;==================================================
; A function to return azimuthal averaged quantities
; of a field defined on a cylindrical grid.  
; Averages are calculated using the domain center as
; the center of the vortex.
;
; INPUT
;    input = FLTARR(n_r, n_theta) = a 2-D field
; 
; OUTPUT
;    <azimuthalmean> = FLTARR(n_r)
;
;    perturbation = FLTARR(n_r, n_theta) = departure 
;                   from azimuthal average
;==================================================


s = SIZE(input)
ret_val = TOTAL( input[*, 1:*], 2 ) / (s[2] - 1)

perturbation = input * 0.0   ; zero the perturbation array 


FOR i = 0, s[2] - 1 DO BEGIN

   perturbation[*, i] = input[*, i] - ret_val
   
ENDFOR

RETURN, ret_val

END



