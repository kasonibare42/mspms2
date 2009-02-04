c*****************************************************************
c
c   FF for silver (12/04/08)
c 
c234567
c     implicit none
c     integer i,n,nsize,ndata
c     parameter (nsize=10,ndata=400)
c     real*8 x(nsize),y(nsize),z(nsize),energy,e,eps,
c    #       xp(nsize),yp(nsize),zp(nsize),ex,ey,ez,
c    #       f,dfx(nsize),dfy(nsize),dfz(nsize)
c     parameter (eps=1.d-7)
c     character*2 atom

c     open(10,file='ffieldag.out',status='unknown')

c     if(nsize.eq.2)open(20,file='Ag2data',status='old')
c     if(nsize.eq.3)open(20,file='Ag3data',status='old')
c     if(nsize.eq.4)open(20,file='Ag4data',status='old')
c     if(nsize.eq.5)open(20,file='Ag5data',status='old')
c     if(nsize.eq.6)open(20,file='Ag6data',status='old')
c     if(nsize.eq.7)open(20,file='Ag7data',status='old')
c     if(nsize.eq.8)open(20,file='Ag8data',status='old')
c     if(nsize.eq.9)open(20,file='Ag9data',status='old')
c     if(nsize.eq.10)open(20,file='Ag10data',status='old')
c     if(nsize.eq.11)open(20,file='Ag11data',status='old')
c     if(nsize.eq.12)open(20,file='Ag12data',status='old')
c     if(nsize.eq.13)open(20,file='Ag13data',status='old')
c     if(nsize.eq.14)open(20,file='Ag14data',status='old')
c     if(nsize.eq.15)open(20,file='Ag15data',status='old')
c     if(nsize.eq.16)open(20,file='Ag16data',status='old')
c     if(nsize.eq.17)open(20,file='Ag17data',status='old')
c     if(nsize.eq.18)open(20,file='Ag18data',status='old')
c     if(nsize.eq.19)open(20,file='Ag19data',status='old')
c     if(nsize.eq.20)open(20,file='Ag20data',status='old')
c     if(nsize.eq.21)open(20,file='Ag21data',status='old')
c     if(nsize.eq.22)open(20,file='Ag22data',status='old')
c     if(nsize.eq.23)open(20,file='Ag23data',status='old')
c     if(nsize.eq.24)open(20,file='Ag24data',status='old')
c     if(nsize.eq.25)open(20,file='Ag25data',status='old')
c     if(nsize.eq.26)open(20,file='Ag26data',status='old')
c     if(nsize.eq.27)open(20,file='Ag27data',status='old')
c     if(nsize.eq.30)open(20,file='Ag30data',status='old')
c     if(nsize.eq.31)open(20,file='Ag31data',status='old')
c     if(nsize.eq.32)open(20,file='Ag32data',status='old')
c     if(nsize.eq.37)open(20,file='Ag37data',status='old')
c     if(nsize.eq.38)open(20,file='Ag38data',status='old')
c     if(nsize.eq.39)open(20,file='Ag39data',status='old')
c     if(nsize.eq.40)open(20,file='Ag40data',status='old')
c     if(nsize.eq.55)open(20,file='Ag55data',status='old')
c     if(nsize.eq.147)open(20,file='Ag147data',status='old')
c     if(nsize.eq.256)open(20,file='Ag256data',status='old')
c     if(nsize.eq.500)open(20,file='Ag500data',status='old')

c     do 100 n=1,ndata
c     read(20,*)
c     read(20,*)e   !  DFT energy
c     do i=1,nsize
c     read(20,*)atom,x(i),y(i),z(i)
c     xp(i)=x(i)
c     yp(i)=y(i)
c     zp(i)=z(i)
c     end do
c     read(20,*)
c     call ffieldag(nsize,x,y,z,f,dfx,dfy,dfz)
c     write(10,*)nsize,n,e/nsize,f
c100  continue
c  
c     stop
c     end

c*****************************************************************
c
c  FF code to evaluate energy and derivatives for Ag-clusters
c
c  Computes derivatives using analytic formula for large-N
c                       and finite differences for small-N 
c
c234567
      subroutine ffieldag(nm,x,y,z,energy,dfx,dfy,dfz)
      implicit none
      integer nm,nmax,i,j,ip,nlarge
      parameter (nmax=2000,nlarge=55)
      real*8 x(nm),y(nm),z(nm),dfx(nm),dfy(nm),dfz(nm),
     #       xp(nmax),yp(nmax),zp(nmax),energy,
     #       r(nmax,nmax),rho(nmax),gc(nmax,nmax),nc(nmax),
     #       rijx(nmax,nmax),rijy(nmax,nmax),rijz(nmax,nmax),
     #       lc1(nmax,nmax),mc1(nmax),lc2(nmax,nmax),mc2(nmax),
     #       lc3(nmax,nmax),mc3(nmax),lc4(nmax,nmax),mc4(nmax),
     #       lc5(nmax,nmax),mc5(nmax),lc6(nmax,nmax),mc6(nmax),
     #       lc7(nmax,nmax),mc7(nmax),lc8(nmax,nmax),mc8(nmax),
     #       lc1x(nmax,nmax),mc1x(nmax),lc2x(nmax,nmax),mc2x(nmax),
     #       lc3x(nmax,nmax),mc3x(nmax),lc4x(nmax,nmax),mc4x(nmax),
     #       lc5x(nmax,nmax),mc5x(nmax),lc6x(nmax,nmax),mc6x(nmax),
     #       lc7x(nmax,nmax),mc7x(nmax),lc8x(nmax,nmax),mc8x(nmax),
     #       lc1y(nmax,nmax),mc1y(nmax),lc2y(nmax,nmax),mc2y(nmax),
     #       lc3y(nmax,nmax),mc3y(nmax),lc4y(nmax,nmax),mc4y(nmax),
     #       lc5y(nmax,nmax),mc5y(nmax),lc6y(nmax,nmax),mc6y(nmax),
     #       lc7y(nmax,nmax),mc7y(nmax),lc8y(nmax,nmax),mc8y(nmax),
     #       lc1z(nmax,nmax),mc1z(nmax),lc2z(nmax,nmax),mc2z(nmax),
     #       lc3z(nmax,nmax),mc3z(nmax),lc4z(nmax,nmax),mc4z(nmax),
     #       lc5z(nmax,nmax),mc5z(nmax),lc6z(nmax,nmax),mc6z(nmax),
     #       lc7z(nmax,nmax),mc7z(nmax),lc8z(nmax,nmax),mc8z(nmax),
     #       sum1x,sum1y,sum1z,sum2x,sum2y,sum2z,sum3x,sum3y,sum3z,
     #       sum4x,sum4y,sum4z,sum5x,sum5y,sum5z,sum6x,sum6y,sum6z,
     #       sum7x,sum7y,sum7z,sum8x,sum8y,sum8z,sum9x,sum9y,sum9z,
     #       sum,sum0,sum1,sum2,sum3,sum4,sum5,sum6,sum7,sum8,sum9,
     #       xij,yij,zij,rij,vr,va,dvax,dvay,dvaz,dvrx,dvry,dvrz,
     #       rsq,fcutoffag,eij0,cij0,aij0,eij,cij,aij,eps,
     #       arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,
     #       e1,e2,e3,e4,e5,e6,e7,e8,e9,e10,e11,e12,
     #       c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,
     #       a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,
     #       sum10,sum11,sum12,arg10,arg11,arg12,
     #       rmin1,rmax1,rmin2,rmax2,rmin3,rmax3,
     #       rmin4,rmax4,rmin5,rmax5,rmin6,rmax6,
     #       rmin7,rmax7,rmin8,rmax8,rmin9,rmax9
      data e1,c1,a1 /-0.165501,  0.178600,  0.044569/
      data e2,c2,a2 /-0.164846,  1.460703,  0.105401/
      data e3,c3,a3 / 0.171883, -1.488742, -0.112043/
      data e4,c4,a4 / 0.818292, -0.649727, -0.194487/
      data e5,c5,a5 /-0.702180,  0.543692,  0.159137/
      data e6,c6,a6 / 0.021457,  0.011655, -0.007669/
      data e7,c7,a7 /-0.023476,  0.017641, -0.000518/
      data e8,c8,a8 /-0.459144,  1.920249,  0.282291/
      data e9,c9,a9 /-0.064645,  0.998784,  0.124938/
      data e10,c10,a10 /-0.037492,  0.998988,  0.092216/
      data e11,c11,a11 /-0.025497,  1.000000,  0.066137/
      data e12,c12,a12 /-0.019024,  1.000000,  0.013612/
      data rmin1,rmax1 /-0.408828,  0.217188/
      data rmin2,rmax2 / 0.422726,  1.186462/
      data rmin3,rmax3 /-0.087055,  3.226801/
      data rmin4,rmax4 /-0.486428,  4.427632/
      data rmin5,rmax5 / 0.520539, 11.181866/
      data rmin6,rmax6 / 1.281913,  4.813950/
      data rmin7,rmax7 / 4.716074,  6.072361/
      data rmin8,rmax8 / 0.164882,  6.621761/

      if(nm.gt.nmax)write(*,*)'dimension error in dfield'

      eij0=0.796671027d0
      cij0=1.52968004d0
      aij0=2.51007759d0

      do i=1,nm
      sum=0.d0
      do j=1,nm
        rsq=(x(i)-x(j))**2+(y(i)-y(j))**2+(z(i)-z(j))**2
        r(i,j)=sqrt(rsq)
        gc(i,j)=fcutoffag(r(i,j),7.d0,10.d0)
        sum=sum+gc(i,j)
      end do
      nc(i)=sum-1
      xp(i)=x(i)
      yp(i)=y(i)
      zp(i)=z(i)
      end do

      arg1=2.d0                            !N=3
      arg2=2.d0                            !N=4
      arg3=2.d0                            !N=5
      arg4=2.d0                            !N=6
      arg5=2.d0                            !N=7
      arg6=2.d0                            !N=8
      arg7=2.d0                            !N=9
      arg8=1.d0                            !N=10
      do i=1,nm
      arg9=min(1,max(0,nc(i)-13.d0))       !N=15
      arg10=min(1,max(0,nc(i)-14.d0))      !N=16
      arg11=min(1,max(0,nc(i)-16.d0))      !N=18
      arg12=min(1,max(0,nc(i)-53.d0))      !N=55
      aij=aij0+a1*arg1+a2*arg2+a3*arg3+a4*arg4+a5*arg5
     #        +a6*arg6+a7*arg7+a8*arg8+a9*arg9+a10*arg10
     #        +a11*arg11+a12*arg12
      sum=0.d0
      do j=1,nm
        rij=r(i,j)
        if(j.ne.i)sum=sum+(aij/rij)**5
      end do
      rho(i)=dsqrt(sum)
      end do

      sum0=0.d0
      sum1=0.d0
      sum2=0.d0
      sum3=0.d0
      do i=1,nm
      if(nc(i).gt.nlarge)then  ! use analytic formula
      arg9=min(1,max(0,nc(i)-13.d0))       !N=15
      arg10=min(1,max(0,nc(i)-14.d0))      !N=16
      arg11=min(1,max(0,nc(i)-16.d0))      !N=18
      arg12=min(1,max(0,nc(i)-53.d0))      !N=55
      eij=eij0+e1*arg1+e2*arg2+e3*arg3+e4*arg4+e5*arg5
     #        +e6*arg6+e7*arg7+e8*arg8+e9*arg9+e10*arg10
     #        +e11*arg11+e12*arg12
      cij=cij0+c1*arg1+c2*arg2+c3*arg3+c4*arg4+c5*arg5
     #        +c6*arg6+c7*arg7+c8*arg8+c9*arg9+c10*arg10
     #        +c11*arg11+c12*arg12
      aij=aij0+a1*arg1+a2*arg2+a3*arg3+a4*arg4+a5*arg5
     #        +a6*arg6+a7*arg7+a8*arg8+a9*arg9+a10*arg10
     #        +a11*arg11+a12*arg12
      vr=0.d0
      va=0.d0
      dvrx=0.d0
      dvax=0.d0
      dvry=0.d0
      dvay=0.d0
      dvrz=0.d0
      dvaz=0.d0
      do j=1,nm
        if(j.ne.i)then
        rij=r(i,j)
        xij=x(i)-x(j)
        yij=y(i)-y(j)
        zij=z(i)-z(j)
        vr=vr+eij*(aij/rij)**10.d0
        va=va+(cij*eij)**2*(aij/rij)**5.d0
        dvrx=dvrx-10.d0*aij**10.d0*xij/rij**12.d0
        dvax=dvax-5.d0*aij**5.d0*xij/rij**7.d0/2.d0/rho(i)
     #           -5.d0*aij**5.d0*xij/rij**7.d0/2.d0/rho(j)
        dvry=dvry-10.d0*aij**10.d0*yij/rij**12.d0
        dvay=dvay-5.d0*aij**5.d0*yij/rij**7.d0/2.d0/rho(i)
     #           -5.d0*aij**5.d0*yij/rij**7.d0/2.d0/rho(j)
        dvrz=dvrz-10.d0*aij**10.d0*zij/rij**12.d0
        dvaz=dvaz-5.d0*aij**5.d0*zij/rij**7.d0/2.d0/rho(i)
     #           -5.d0*aij**5.d0*zij/rij**7.d0/2.d0/rho(j)
        endif
      end do
      sum0=sum0+vr/2.d0-sqrt(va)
      sum1=eij*dvrx-cij*eij*dvax
      sum2=eij*dvry-cij*eij*dvay
      sum3=eij*dvrz-cij*eij*dvaz
      energy=sum0/nm
      dfx(i)=sum1/nm
      dfy(i)=sum2/nm
      dfz(i)=sum3/nm

      else  ! use finite differences

      eps=1.d-7

      xp(i)=x(i)+eps
      yp(i)=y(i)+eps
      zp(i)=z(i)+eps

      do ip=1,nm
      sum1=0.d0
      sum2=0.d0
      sum3=0.d0
      sum4=0.d0
      sum5=0.d0
      sum6=0.d0
      sum7=0.d0
      sum8=0.d0
      sum1x=0.d0
      sum1y=0.d0
      sum1z=0.d0
      sum2x=0.d0
      sum2y=0.d0
      sum2z=0.d0
      sum3x=0.d0
      sum3y=0.d0
      sum3z=0.d0
      sum4x=0.d0
      sum4y=0.d0
      sum4z=0.d0
      sum5x=0.d0
      sum5y=0.d0
      sum5z=0.d0
      sum6x=0.d0
      sum6y=0.d0
      sum6z=0.d0
      sum7x=0.d0
      sum7y=0.d0
      sum7z=0.d0
      sum8x=0.d0
      sum8y=0.d0
      sum8z=0.d0
      do j=1,nm
      lc1(ip,j)=fcutoffag(r(ip,j),rmin1,rmax1)
      lc2(ip,j)=fcutoffag(r(ip,j),rmin2,rmax2)
      lc3(ip,j)=fcutoffag(r(ip,j),rmin3,rmax3)
      lc4(ip,j)=fcutoffag(r(ip,j),rmin4,rmax4)
      lc5(ip,j)=fcutoffag(r(ip,j),rmin5,rmax5)
      lc6(ip,j)=fcutoffag(r(ip,j),rmin6,rmax6)
      lc7(ip,j)=fcutoffag(r(ip,j),rmin7,rmax7)
      lc8(ip,j)=fcutoffag(r(ip,j),rmin8,rmax8)
      sum1=sum1+lc1(ip,j)
      sum2=sum2+lc2(ip,j)
      sum3=sum3+lc3(ip,j)
      sum4=sum4+lc4(ip,j)
      sum5=sum5+lc5(ip,j)
      sum6=sum6+lc6(ip,j)
      sum7=sum7+lc7(ip,j)
      sum8=sum8+lc8(ip,j)
      rijx(ip,j)=sqrt((xp(ip)-xp(j))**2+(y(ip)-y(j))**2+(z(ip)-z(j))**2)
      rijy(ip,j)=sqrt((x(ip)-x(j))**2+(yp(ip)-yp(j))**2+(z(ip)-z(j))**2)
      rijz(ip,j)=sqrt((x(ip)-x(j))**2+(y(ip)-y(j))**2+(zp(ip)-zp(j))**2)
      lc1x(ip,j)=fcutoffag(rijx(ip,j),rmin1,rmax1)
      lc1y(ip,j)=fcutoffag(rijy(ip,j),rmin1,rmax1)
      lc1z(ip,j)=fcutoffag(rijz(ip,j),rmin1,rmax1)
      lc2x(ip,j)=fcutoffag(rijx(ip,j),rmin2,rmax2)
      lc2y(ip,j)=fcutoffag(rijy(ip,j),rmin2,rmax2)
      lc2z(ip,j)=fcutoffag(rijz(ip,j),rmin2,rmax2)
      lc3x(ip,j)=fcutoffag(rijx(ip,j),rmin3,rmax3)
      lc3y(ip,j)=fcutoffag(rijy(ip,j),rmin3,rmax3)
      lc3z(ip,j)=fcutoffag(rijz(ip,j),rmin3,rmax3)
      lc4x(ip,j)=fcutoffag(rijx(ip,j),rmin4,rmax4)
      lc4y(ip,j)=fcutoffag(rijy(ip,j),rmin4,rmax4)
      lc4z(ip,j)=fcutoffag(rijz(ip,j),rmin4,rmax4)
      lc5x(ip,j)=fcutoffag(rijx(ip,j),rmin5,rmax5)
      lc5y(ip,j)=fcutoffag(rijy(ip,j),rmin5,rmax5)
      lc5z(ip,j)=fcutoffag(rijz(ip,j),rmin5,rmax5)
      lc6x(ip,j)=fcutoffag(rijx(ip,j),rmin6,rmax6)
      lc6y(ip,j)=fcutoffag(rijy(ip,j),rmin6,rmax6)
      lc6z(ip,j)=fcutoffag(rijz(ip,j),rmin6,rmax6)
      lc7x(ip,j)=fcutoffag(rijx(ip,j),rmin7,rmax7)
      lc7y(ip,j)=fcutoffag(rijy(ip,j),rmin7,rmax7)
      lc7z(ip,j)=fcutoffag(rijz(ip,j),rmin7,rmax7)
      lc8x(ip,j)=fcutoffag(rijx(ip,j),rmin8,rmax8)
      lc8y(ip,j)=fcutoffag(rijy(ip,j),rmin8,rmax8)
      lc8z(ip,j)=fcutoffag(rijz(ip,j),rmin8,rmax8)
      sum1x=sum1x+lc1x(ip,j)
      sum1y=sum1y+lc1y(ip,j)
      sum1z=sum1z+lc1z(ip,j)
      sum2x=sum2x+lc2x(ip,j)
      sum2y=sum2y+lc2y(ip,j)
      sum2z=sum2z+lc2z(ip,j)
      sum3x=sum3x+lc3x(ip,j)
      sum3y=sum3y+lc3y(ip,j)
      sum3z=sum3z+lc3z(ip,j)
      sum4x=sum4x+lc4x(ip,j)
      sum4y=sum4y+lc4y(ip,j)
      sum4z=sum4z+lc4z(ip,j)
      sum5x=sum5x+lc5x(ip,j)
      sum5y=sum5y+lc5y(ip,j)
      sum5z=sum5z+lc5z(ip,j)
      sum6x=sum6x+lc6x(ip,j)
      sum6y=sum6y+lc6y(ip,j)
      sum6z=sum6z+lc6z(ip,j)
      sum7x=sum7x+lc7x(ip,j)
      sum7y=sum7y+lc7y(ip,j)
      sum7z=sum7z+lc7z(ip,j)
      sum8x=sum8x+lc8x(ip,j)
      sum8y=sum8y+lc8y(ip,j)
      sum8z=sum8z+lc8z(ip,j)
      end do
      mc1(ip)=sum1-1
      mc2(ip)=sum2-1
      mc3(ip)=sum3-1
      mc4(ip)=sum4-1
      mc5(ip)=sum5-1
      mc6(ip)=sum6-1
      mc7(ip)=sum7-1
      mc8(ip)=sum8-1
      mc1x(ip)=sum1x-1
      mc1y(ip)=sum1y-1
      mc1z(ip)=sum1z-1
      mc2x(ip)=sum2x-1
      mc2y(ip)=sum2y-1
      mc2z(ip)=sum2z-1
      mc3x(ip)=sum3x-1
      mc3y(ip)=sum3y-1
      mc3z(ip)=sum3z-1
      mc4x(ip)=sum4x-1
      mc4y(ip)=sum4y-1
      mc4z(ip)=sum4z-1
      mc5x(ip)=sum5x-1
      mc5y(ip)=sum5y-1
      mc5z(ip)=sum5z-1
      mc6x(ip)=sum6x-1
      mc6y(ip)=sum6y-1
      mc6z(ip)=sum6z-1
      mc7x(ip)=sum7x-1
      mc7y(ip)=sum7y-1
      mc7z(ip)=sum7z-1
      mc8x(ip)=sum8x-1
      mc8y(ip)=sum8y-1
      mc8z(ip)=sum8z-1
      end do

      sum0=0.d0
      sum1=0.d0
      sum2=0.d0
      sum3=0.d0

      do ip=1,nm
      vr=0.d0
      va=0.d0
      do j=1,nm
        if(j.ne.ip)then
        rij=r(ip,j)
        arg1=min(2,max(0,mc1(ip)+mc1(j)-2.d0*lc1(ip,j)))           !N=3
        arg2=min(2,max(0,mc2(ip)+mc2(j)-2.d0*lc2(ip,j)-2.d0))      !N=4
        arg3=min(2,max(0,mc3(ip)+mc3(j)-2.d0*lc3(ip,j)-4.d0))      !N=5
        arg4=min(2,max(0,mc4(ip)+mc4(j)-2.d0*lc4(ip,j)-6.d0))      !N=6
        arg5=min(2,max(0,mc5(ip)+mc5(j)-2.d0*lc5(ip,j)-8.d0))      !N=7
        arg6=min(2,max(0,mc6(ip)+mc6(j)-2.d0*lc6(ip,j)-10.d0))     !N=8
        arg7=min(2,max(0,mc7(ip)+mc7(j)-2.d0*lc7(ip,j)-12.d0))     !N=9
        arg8=min(1,max(0,mc8(ip)+mc8(j)-2.d0*lc8(ip,j)-14.d0))     !N=10
        arg9=min(1,max(0,nc(ip)-13.d0))       !N=15
        arg10=min(1,max(0,nc(ip)-14.d0))      !N=16
        arg11=min(1,max(0,nc(ip)-16.d0))      !N=18
        arg12=min(1,max(0,nc(ip)-53.d0))      !N=55
        eij=eij0+e1*arg1+e2*arg2+e3*arg3+e4*arg4+e5*arg5
     #          +e6*arg6+e7*arg7+e8*arg8+e9*arg9+e10*arg10
     #          +e11*arg11+e12*arg12
        cij=cij0+c1*arg1+c2*arg2+c3*arg3+c4*arg4+c5*arg5
     #          +c6*arg6+c7*arg7+c8*arg8+c9*arg9+c10*arg10
     #          +c11*arg11+c12*arg12
        aij=aij0+a1*arg1+a2*arg2+a3*arg3+a4*arg4+a5*arg5
     #          +a6*arg6+a7*arg7+a8*arg8+a9*arg9+a10*arg10
     #          +a11*arg11+a12*arg12
        vr=vr+eij*(aij/rij)**10.d0
        va=va+(cij*eij)**2*(aij/rij)**5.d0
        endif
      end do
      sum0=sum0+vr/2.d0-sqrt(va)
      end do

      do ip=1,nm
      vr=0.d0
      va=0.d0
      do j=1,nm
        if(j.ne.ip)then
        rij=rijx(ip,j)
        arg1=min(2,max(0,mc1x(ip)+mc1x(j)-2.d0*lc1x(ip,j)))           !N=3
        arg2=min(2,max(0,mc2x(ip)+mc2x(j)-2.d0*lc2x(ip,j)-2.d0))      !N=4
        arg3=min(2,max(0,mc3x(ip)+mc3x(j)-2.d0*lc3x(ip,j)-4.d0))      !N=5
        arg4=min(2,max(0,mc4x(ip)+mc4x(j)-2.d0*lc4x(ip,j)-6.d0))      !N=6
        arg5=min(2,max(0,mc5x(ip)+mc5x(j)-2.d0*lc5x(ip,j)-8.d0))      !N=7
        arg6=min(2,max(0,mc6x(ip)+mc6x(j)-2.d0*lc6x(ip,j)-10.d0))     !N=8
        arg7=min(2,max(0,mc7x(ip)+mc7x(j)-2.d0*lc7x(ip,j)-12.d0))     !N=9
        arg8=min(1,max(0,mc8x(ip)+mc8x(j)-2.d0*lc8x(ip,j)-14.d0))     !N=10
        arg9=min(1,max(0,nc(ip)-13.d0))       !N=15
        arg10=min(1,max(0,nc(ip)-14.d0))      !N=16
        arg11=min(1,max(0,nc(ip)-16.d0))      !N=18
        arg12=min(1,max(0,nc(ip)-53.d0))      !N=55
        eij=eij0+e1*arg1+e2*arg2+e3*arg3+e4*arg4+e5*arg5
     #          +e6*arg6+e7*arg7+e8*arg8+e9*arg9+e10*arg10
     #          +e11*arg11+e12*arg12
        cij=cij0+c1*arg1+c2*arg2+c3*arg3+c4*arg4+c5*arg5
     #          +c6*arg6+c7*arg7+c8*arg8+c9*arg9+c10*arg10
     #          +c11*arg11+c12*arg12
        aij=aij0+a1*arg1+a2*arg2+a3*arg3+a4*arg4+a5*arg5
     #          +a6*arg6+a7*arg7+a8*arg8+a9*arg9+a10*arg10
     #          +a11*arg11+a12*arg12
        vr=vr+eij*(aij/rij)**10.d0
        va=va+(cij*eij)**2*(aij/rij)**5.d0
        endif
      end do
      sum1=sum1+vr/2.d0-sqrt(va)
      end do

      do ip=1,nm
      vr=0.d0
      va=0.d0
      do j=1,nm
        if(j.ne.ip)then
        rij=rijy(ip,j)
        arg1=min(2,max(0,mc1y(ip)+mc1y(j)-2.d0*lc1y(ip,j)))           !N=3
        arg2=min(2,max(0,mc2y(ip)+mc2y(j)-2.d0*lc2y(ip,j)-2.d0))      !N=4
        arg3=min(2,max(0,mc3y(ip)+mc3y(j)-2.d0*lc3y(ip,j)-4.d0))      !N=5
        arg4=min(2,max(0,mc4y(ip)+mc4y(j)-2.d0*lc4y(ip,j)-6.d0))      !N=6
        arg5=min(2,max(0,mc5y(ip)+mc5y(j)-2.d0*lc5y(ip,j)-8.d0))      !N=7
        arg6=min(2,max(0,mc6y(ip)+mc6y(j)-2.d0*lc6y(ip,j)-10.d0))     !N=8
        arg7=min(2,max(0,mc7y(ip)+mc7y(j)-2.d0*lc7y(ip,j)-12.d0))     !N=9
        arg8=min(1,max(0,mc8y(ip)+mc8y(j)-2.d0*lc8y(ip,j)-14.d0))     !N=10
        arg9=min(1,max(0,nc(ip)-13.d0))       !N=15
        arg10=min(1,max(0,nc(ip)-14.d0))      !N=16
        arg11=min(1,max(0,nc(ip)-16.d0))      !N=18
        arg12=min(1,max(0,nc(ip)-53.d0))      !N=55
        eij=eij0+e1*arg1+e2*arg2+e3*arg3+e4*arg4+e5*arg5
     #          +e6*arg6+e7*arg7+e8*arg8+e9*arg9+e10*arg10
     #          +e11*arg11+e12*arg12
        cij=cij0+c1*arg1+c2*arg2+c3*arg3+c4*arg4+c5*arg5
     #          +c6*arg6+c7*arg7+c8*arg8+c9*arg9+c10*arg10
     #          +c11*arg11+c12*arg12
        aij=aij0+a1*arg1+a2*arg2+a3*arg3+a4*arg4+a5*arg5
     #          +a6*arg6+a7*arg7+a8*arg8+a9*arg9+a10*arg10
     #          +a11*arg11+a12*arg12
        vr=vr+eij*(aij/rij)**10.d0
        va=va+(cij*eij)**2*(aij/rij)**5.d0
        endif
      end do
      sum2=sum2+vr/2.d0-sqrt(va)
      end do

      do ip=1,nm
      vr=0.d0
      va=0.d0
      do j=1,nm
        if(j.ne.ip)then
        rij=rijz(ip,j)
        arg1=min(2,max(0,mc1z(ip)+mc1z(j)-2.d0*lc1z(ip,j)))           !N=3
        arg2=min(2,max(0,mc2z(ip)+mc2z(j)-2.d0*lc2z(ip,j)-2.d0))      !N=4
        arg3=min(2,max(0,mc3z(ip)+mc3z(j)-2.d0*lc3z(ip,j)-4.d0))      !N=5
        arg4=min(2,max(0,mc4z(ip)+mc4z(j)-2.d0*lc4z(ip,j)-6.d0))      !N=6
        arg5=min(2,max(0,mc5z(ip)+mc5z(j)-2.d0*lc5z(ip,j)-8.d0))      !N=7
        arg6=min(2,max(0,mc6z(ip)+mc6z(j)-2.d0*lc6z(ip,j)-10.d0))     !N=8
        arg7=min(2,max(0,mc7z(ip)+mc7z(j)-2.d0*lc7z(ip,j)-12.d0))     !N=9
        arg8=min(1,max(0,mc8z(ip)+mc8z(j)-2.d0*lc8z(ip,j)-14.d0))     !N=10
        arg9=min(1,max(0,nc(ip)-13.d0))       !N=15
        arg10=min(1,max(0,nc(ip)-14.d0))      !N=16
        arg11=min(1,max(0,nc(ip)-16.d0))      !N=18
        arg12=min(1,max(0,nc(ip)-53.d0))      !N=55
        eij=eij0+e1*arg1+e2*arg2+e3*arg3+e4*arg4+e5*arg5
     #          +e6*arg6+e7*arg7+e8*arg8+e9*arg9+e10*arg10
     #          +e11*arg11+e12*arg12
        cij=cij0+c1*arg1+c2*arg2+c3*arg3+c4*arg4+c5*arg5
     #          +c6*arg6+c7*arg7+c8*arg8+c9*arg9+c10*arg10
     #          +c11*arg11+c12*arg12
        aij=aij0+a1*arg1+a2*arg2+a3*arg3+a4*arg4+a5*arg5
     #          +a6*arg6+a7*arg7+a8*arg8+a9*arg9+a10*arg10
     #          +a11*arg11+a12*arg12
        vr=vr+eij*(aij/rij)**10.d0
        va=va+(cij*eij)**2*(aij/rij)**5.d0
        endif
      end do
      sum3=sum3+vr/2.d0-sqrt(va)
      end do

      energy=sum0/nm
      dfx(i)=(sum1-sum0)/(nm*eps)
      dfy(i)=(sum2-sum0)/(nm*eps)
      dfz(i)=(sum3-sum0)/(nm*eps)
      xp(i)=x(i)
      yp(i)=y(i)
      zp(i)=z(i)
      endif
      end do

      return
      end

c************************************************************
c
c234567
      real*8 function fcutoffag(r,rminshift,rmaxshift)
      implicit none
      integer index
      real*8 r,rmin,rmax,pi,arg,re,rminshift,rmaxshift

      pi=dacos(-1.d0)

      rmin=3.0d0+rminshift
      rmax=5.0d0+rmaxshift

      if(r.lt.rmin)then
        fcutoffag=1.d0
      elseif(r.lt.rmax)then
        arg=pi*(r-rmin)/(rmax-rmin)
        fcutoffag=(1.d0+cos(arg))/2.d0
      else
        fcutoffag=0.d0
      endif

      return
      end
