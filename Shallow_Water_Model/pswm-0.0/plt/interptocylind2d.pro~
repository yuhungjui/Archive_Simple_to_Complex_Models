;===================================================
;  Procedure to interpolate Cartesian grid data to 
;  cylindrical grid data
;===================================================

PRO INTERPTOCYLIND2D, squareArray, $
                    cylindricalArray, $
	            numberRCords, $
	            numberThetaCords, $
		    xCenter, $
		    yCenter, $
                    xPositions, $
                    yPositions
 
  deltaTheta = 2. * !PI / (numberThetaCords - 1)
  
  cylindricalArray = FLTARR(numberRCords, numberThetaCords)
  xPositions = FLTARR(numberRCords, numberThetaCords)
  yPositions = FLTARR(numberRCords, numberThetaCords)

  r = FLOAT( LINDGEN(numberRCords, numberThetaCords) MOD numberRCords )
  theta = FLOAT( ( LINDGEN(numberRCords, numberThetaCords) / numberRCords ) * $
                    deltaTheta )

  xPositions = r * COS(theta) + xCenter  
  yPositions = r * SIN(theta) + yCenter
        
  cylindricalArray[*, *] = INTERPOLATE( squareArray[*, *], $
                                                 xPositions, yPositions )

RETURN

END
