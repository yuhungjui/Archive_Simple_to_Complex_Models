Hey Eric,

Alrighy... PS color plotting routines.  Attached are a couple different 
ones that I've used.  I'll chat about those in a bit.

First, how to invoke them.

So, somewhere in your main program, it helps to set up your plotting 
colors.  Almost all of my programs use a 17-color palette.  Below are 
the two color schemes I use most: rainbow and cold2warm.  The numbers 
are values of the red, blue, green hues from 0 to 255 (pretty standard 
RGB convention).  Included in the routines I sent along are two things: 
new_color_table.pro, which makes sure that previous color tables are 
deleted before you set up yours; and c32.pro, which converts the 0-255 
values to the proper LONG integer values that the PS-color routines use.

;######################################################################
; Set up plotting colors
;######################################################################

new_color_table

;$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
; Rainbow Colors.
;rbw_ELM_LIST = [ 1    2    3    4    5    6    7    8    9    10   11   
12   13   14   15   16   17]
rbw_red_color = [153, 125, 077, 000, 070, 070, 070, 070, 050, 200, 255, 
238, 238, 220, 238, 220, 255]
rbw_gre_color = [050, 095, 000, 000, 070, 150, 235, 255, 220, 255, 255, 
205, 154, 105, 069, 075, 000]
rbw_blu_color = [204, 186, 186, 205, 255, 255, 255, 160, 050, 050, 000, 
000, 000, 000, 000, 075, 000]

;$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
; Cold to warm
;c2w_ELM_LIST = [ 1    2    3    4    5    6    7    8    9    10   11   
12   13   14   15   16   17]
c2w_red_color = [000, 032, 064, 096, 128, 160, 191, 223, 254, 255, 255, 
255, 255, 255, 255, 255, 255]
c2w_gre_color = [000, 032, 064, 096, 128, 160, 191, 223, 254, 213, 181, 
150, 118, 086, 054, 022, 000]
c2w_blu_color = [255, 255, 255, 255, 255, 255, 255, 255, 254, 223, 191, 
160, 128, 096, 064, 032, 000]

c2w_plot_colors = c32(c2w_red_color, c2w_gre_color, c2w_blu_color, 1)
rbw_plot_colors = c32(rbw_red_color, rbw_gre_color, rbw_blu_color, 1)
NLEV1 = FIX(SIZE(c2w_plot_colors, /DIMENSIONS))
NLEV  = NLEV1(0)
black_cons = c32(intarr(NLEV), intarr(NLEV), intarr(NLEV), 1)

Once you have these, then things are easy with my routines.  Below is a 
piece of code that invokes one of the routines:
    ;************************************************************
    ; U wind sections...
    ;************************************************************
    CON_SPACE     = 2.0
    COLORSCHEME   = rbw_plot_colors
    CONTOURLEVELS = findgen(NLEV) * CON_SPACE - 16.0
    BOUNDS        = [MIN(XPLOT2), MAX(XPLOT2), MIN(YPLOT2), MAX(YPLOT2), 
TIMEBOUND]
    for k = 0, ZMAX - 1, 3 do begin
        LEVELDIR = '/Level'+STRCOMPRESS(STRING(k), /remove_all)
        isdirthere = file_test(COMPDIR+LEVELDIR, /directory)
        if(isdirthere EQ 0) then file_mkdir, COMPDIR+LEVELDIR
        isdirthere = file_test(COMPDIR+LEVELDIR+'/U', /directory)
        if(isdirthere EQ 0) then file_mkdir, COMPDIR+LEVELDIR+'/U'
        TITLES    = ['!17U (m/s) at z = ' + 
STRMID(STRCOMPRESS(STRING(ZPLOT(k)), /remove_all), 0, 5) + '!17 km', $
                     '!17x (km)', '!17y (km)!3']
        PLOTNAME  = COMPDIR + LEVELDIR + '/U/u'+STRCOMPRESS(STRING(t), 
/remove_all)+'.ps'
        twodplotter, REFORM(U(*, *, k)), XPLOT2, YPLOT2, TITLES, 
COLORSCHEME, CONTOURLEVELS, BOUNDS, PLOTNAME
    endfor

The code invokes the twodplotter routine.  The first argument is the 2D 
array.  Second and third arguments are the plotting arrays for the x and 
y dimensions.  The fourth is a three-element string array with the plot 
title, x-title and y-title.  Fifth and sixth are two same size 1-D 
arrays with the color scheme and contour levels defined.  Seventh is a 
five-element array with the minimum and maximum x points, the minimum 
and maximum y points, and the time.  Lastly is the full path to the file 
the plot will be put into.

Most of the plotting routines I sent along use this same sort of syntax:
functionname, DATA, X-ARRAY, (Y/Z)-ARRAY, TITLES, COLORS, LEVELS, 
BOUNDS, FILE

Just check out the routines and you should be able to figure out a good 
number of them.  The routines do all the work for you, including color 
bars (if appropriate) and zero-contour overlaying.  Feel free to modify 
those plotting routines to do something else, if you like.

Let me know if you need any help with any of them.

-Wes


-- 
=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
Wesley Terwey
Ph.D., Dept. of Atmospheric Sciences
Colorado State University
=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
There are no data that cannot be plotted on a
straight line if the axis are chosen correctly.
=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
