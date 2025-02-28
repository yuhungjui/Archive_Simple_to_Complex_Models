;===================================================

PRO CARTESIANTOCYLIND2D, u, v, vr, vt,dx, dy, nazm, $
		       xCenter, $
		       yCenter, $
		       xcoord, $
		       ycoord

;===================================================
; A program which converts Cartesian wind components 
; (u & v) into cylindrical wind components (vr & vt)
; at each Cartesian grid point.  
;
; INPUT
;    u = FLTARR(nx, ny, nz) = Cartesian x-velocity
;
;    v = FLTARR(nx, ny, nz) = Cartesian y-velocity
;
;    nazm = number of azimuthal points to be used
;           on cylindrical grid
;
;    xCenter = index value of x-grid center
;
;    yCenter = index value of y-grid center
;
; OUTPUT
;    vr = FLTARR(nx, ny, nz) = radial wind velocity
;
;    vt = FLTARR(nx, ny, nz) = tangential wind 
;                              velocity
;===================================================


s = SIZE(u)
t = SIZE(v)

vr = FLTARR( s[1], s[2] )
vt = FLTARR( s[1], s[2] )

xcoord = FLTARR( s[1] )
ycoord = FLTARR( s[2] )

xcoord = FINDGEN( s[1] ) * dx - xCenter * dx   ;  Cartesian x-grid locations 
                                                   ;     relative to origin
ycoord = FINDGEN( s[2] ) * dy - yCenter * dy   ;  Cartesian y-grid locations 
                                                   ;     relative to origin
                                               
del_lambda = 2. * !PI / (nazm - 1)

FOR i = 0, s[1] - 1 DO BEGIN
   FOR j = 0, s[2] - 1 DO BEGIN    
      
      IF (xcoord[i] EQ 0.0) AND (ycoord[j] EQ 0.0) THEN BEGIN
      
         lambda = FINDGEN(nazm) * del_lambda
         
         vr_sum = 0.0
         vt_sum = 0.0
         
         FOR azm = 0, nazm - 2 DO BEGIN
         
         vr_sum = vr_sum + u[i, j] * COS( lambda(azm) ) $
                         + v[i, j] * SIN( lambda(azm) )
                         
         vt_sum = vt_sum - u[i, j] * SIN( lambda(azm) ) $
                         + v[i, j] * COS( lambda(azm) )
         
         ENDFOR
         
         vr[i, j] = vr_sum / (nazm - 1)
         vt[i, j] = vt_sum / (nazm - 1)
         
         GOTO, jump1
         
      ENDIF ELSE BEGIN
      
      IF (xcoord[i] EQ 0.0) AND (ycoord[j] GT 0.0) THEN BEGIN
         lambda = !PI / 2.0
      ENDIF ELSE BEGIN
      
      IF (xcoord[i] EQ 0.0) AND (ycoord[j] LT 0.0) THEN BEGIN
         lambda = 3.0 * !PI / 2.0
      ENDIF ELSE BEGIN
      
      IF (xcoord[i] LT 0.0) THEN BEGIN
         lambda = ATAN( ycoord[j] / xcoord[i] ) + !PI
      ENDIF ELSE BEGIN
      
      IF (ycoord[j] GE 0.0) THEN BEGIN
         lambda = ATAN( ycoord[j] / xcoord[i] )
      ENDIF ELSE BEGIN
      
         lambda = ATAN( ycoord[j] / xcoord[i] ) + 2.0 * !PI
	 
      ENDELSE
      ENDELSE
      ENDELSE
      ENDELSE
      ENDELSE
      
      vr[i, j] =   u[i, j] * COS(lambda) + v[i, j] * SIN(lambda)
      vt[i, j] = - u[i, j] * SIN(lambda) + v[i, j] * COS(lambda)
      
      jump1:   ;continue
      
     
   ENDFOR
ENDFOR

END




