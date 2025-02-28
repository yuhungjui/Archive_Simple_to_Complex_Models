;============================================

PRO new_color_table

;============================================
; A procedure to wipe out the color table for use with c32.pro in 8 bit mode.
;
;==========

TVLCT, rtt, gtt, btt, /GET

where_grey = $
      WHERE(rtt + gtt + btt GT 0 AND rtt + gtt + btt LT 3 * 255, count_grey)

IF count_grey GT 0 THEN BEGIN
   rtt[where_grey] = 0
   gtt[where_grey] = 0
   btt[where_grey] = 0
ENDIF

TVLCT, rtt, gtt, btt

END
