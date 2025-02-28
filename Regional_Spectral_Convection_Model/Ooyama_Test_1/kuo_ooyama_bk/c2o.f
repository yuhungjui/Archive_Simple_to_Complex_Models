      parameter ( mx=48, mz=48, nx=mx/3*2, nz=mz/3*2 )
      parameter ( mnx=nx, mnz=nz )

      common/model1/ sigs(-1:nx,0:nz), xis(-1:nx,0:nz),
     & etas(-1:nx,0:nz), uds(-1:nx,0:nz), wds(-1:nx,0:nz),
     & imx(13), imz(13), fmx(3*mx/2+1), fmz(2*mz),
     & f(-1:nx,0:nz,1:5)

      common/pressu/ pp(-1:mx,0:mz), pxp(-1:mx,0:mz+1), 
     & pzp(-1:mx,0:mz+1), bf(-1:mx,0:mz+1),p1(-1:mx,0:mz),
     & p2(-1:mx,0:mz),p3(-1:mx,0:mz)

      common/temperature/ tt(-1:mx,0:mz)

      common/const/ xi0,t0,e0,cw,ci,xlw0,xli0,rv,ra,cvv,cpv,cva,cpa,
     & ggg
 
      common/const1/ xlamd0,seta0

      common/constemp/ dddt

      common/basic/ bar(0:mnz,8), xibar(-1:mnx,0:mnz),
     &                            ttbar(-1:mnx,0:mnz)

      common/ex/ zi,h
      common/job/ dox, doz, dzz
      common/mmm/ zj(0:mz), xxi(-1:mx)
      common/diffusion/ edx, edz
      common/ou/ x(-1:mnx), z(0:mnz)
      common/filter/ fitf(-1:nx), fitc(0:nz), filt(-1:nx,0:nz),
     & fitc2(0:mz)

      common/timeset/ timax,delt
      common/output/ ito
      common/temppert/ tpmax





