pro makecolorbar, startx, $
                  starty, $
                  colors, $
		  xwidth, $
		  yheight, $
		  labels, $
		  borderWidth = borderWidth, $
                  borderColor = borderColor, $
	          labelColor = labelColor
		  

;procedure makeColorBar.pro

  if n_elements(borderWidth) ne 1 then borderwidth = xwidth / 10.0

  if n_elements(borderColor) ne 1 then begin

    borderColor = 255L + 255L * 256L + 256L * 256L * 255L 

  endif

  if n_elements(labelColor) ne 1 then labelColor = borderColor

  numberOfBars = n_elements(labels)
  ywidth = yheight / FLOAT(numberOfBars)

  xcords = [ startx - borderWidth, $
             startx + xwidth + borderWidth, $
             startx + xwidth + borderWidth, $
             startx - borderWidth ]
  
  ycords = [ starty + borderWidth, $
             starty + borderWidth, $
             starty - yheight - borderWidth, $
             starty - yheight - borderWidth ]

  polyfill, xcords, ycords, $
            color = borderColor, $
            /NORMAL

  xcords = [ startx, startx + xwidth, startx + xwidth, startx ]
  ycords = [ starty, starty, starty - ywidth, starty - ywidth ]

  for i = 0, (numberOfBars - 1) do begin

    polyfill, xcords, ycords - (i * ywidth), $
              color = colors((numberOfBars - (i + 1)) mod n_elements(colors)),$
              /NORMAL

    xyouts, startx + xwidth + borderWidth + .01, $
            starty - ((i + 1) * ywidth), $
            labels(numberOfBars - (i + 1)), $
            color = 0, $
            /NORMAL, $
            charsize = 0.7, $
            charthick = 2.0

  endfor

  return

end

