;===================================

FUNCTION c32, r, g, b, po

;===================================
; A function to return a 32 bit color
;
; INPUT
;    r, g, b = BYTE or INT from 0 to 255
;    po = INT = If set, will do an 8 bit color
;
; OUTPUT
;    <c32> = LONG 
;==========

IF NOT KEYWORD_SET(po) THEN BEGIN
   RETURN, LONG(r) + 256L * LONG(g) + 256L * 256L * LONG(b)
ENDIF ELSE BEGIN
   ret_val = r * 0
   
   where_zero = WHERE(r + g + b EQ 0, count_zero)
   IF count_zero GT 0 THEN ret_val[where_zero] = 255
   where_white = WHERE(r + g + b EQ 3 * 255, count_white)
   IF count_white GT 0 THEN ret_val[where_white] = 0
   
   where_not_yet = WHERE(r + g + b GT 0 AND r + g + b LT 3 * 255, count_not_yet)
   FOR j = 0, count_not_yet  - 1 DO BEGIN
      TVLCT, rtt, gtt, btt, /GET
      FOR i = 1, 255 DO BEGIN
         ix = where_not_yet[j]
         IF rtt[i] EQ r[ix] AND gtt[i] EQ g[ix] AND btt[i] EQ b[ix] THEN BEGIN
            ret_val[ix] = i
            GOTO, JUMP
         ENDIF ELSE IF rtt[i] + gtt[i] + btt[i] EQ 0 THEN BEGIN
            rtt[i] = r[ix]
            gtt[i] = g[ix]
            btt[i] = b[ix]
            TVLCT, rtt, gtt, btt
            ret_val[ix] = i
            GOTO, JUMP
         ENDIF
      ENDFOR

      PRINT, 'Ran out of colors'
      STOP
JUMP:
   ENDFOR
   
   RETURN, ret_val
ENDELSE

END

