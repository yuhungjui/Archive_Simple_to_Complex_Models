	subroutine windbar(x,y,u,v,index,size,icol)
c
c       x,y	: The X/Y coordinate position in user coordinates
c       U,V	: The (U,V) or (theta,speed) of a wind vector               

c       NOTE: THE (U,V) OR (THETA,SPEED) IS THE X-Y COMPONENT OR
c	      WIND DIRECTION RELATIVE TO THE X-Y COORDINATE,
c             NOT THE EAST-NORTH DIRECTION.

c	index	: 1: for (U,V) (default)
c		  2: for (theta, speed)
c	size	: The size of the flag
c	icol	: The color of the flag
c
        call wmseti('wbf',0)  ! default=0, 
                              ! draw a sky cover symbol at base of the barb

        call wmseti('col',icol)
        call wmsetr('wbt',.6) ! The size of the full bar
        call wmsetr('wbr',0.) ! The radius of the circle to represent the calm 
        call wmsetr('wbd',0.22) ! The spacing of the tic bar
	DMX=.01
        CALL GETSET(VL,VR,VB,VT,UL,UR,UB,UT,LL)
        VRL = size * DMX / (VR - VL)
        call wmsetr('wbs',vrl)
c	print*,vrl


	if(index.eq.2)then
	  pi=atan(1.)*4.
	  uu=v*sin(u*pi/180.)*1.94 	! 1 m/s = 1.94 knots
	  vv=v*cos(u*pi/180.)*1.94
	else
	  uu=-1.*u*1.94
	  vv=-1.*v*1.94
	endif
	call wmbarb(x,y,uu,vv)
	return
	end
