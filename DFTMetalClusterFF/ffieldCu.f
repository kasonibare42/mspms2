c*****************************************************************
c
c  call forcefield (FF) to evaluate the energy of a copper cluster
c
c  nsize - size of cluster
c
c  ndata - number of clusters with nsize 
c
c  x,y,z - coordinates of atoms in the cluster
c
c  revised (7/30/08) to include analytical derivatives for nsize>25
c
c234567
c     implicit none
c     integer i,n,nsize,ndata
c     parameter (nsize=32,ndata=1)
c     real*8 x(nsize),y(nsize),z(nsize),energy,e,eps,
c    #       xp(nsize),yp(nsize),zp(nsize),ex,ey,ez,
c    #       dfx(nsize),dfy(nsize),dfz(nsize)
c     parameter (eps=1.d-7)
c     character*2 atom

c     open(10,file='ffield.out',status='unknown')
c     open(30,file='ffield.dump',status='unknown')

c     if(nsize.eq.2)open(20,file='Cu2data',status='old')
c     if(nsize.eq.3)open(20,file='Cu3data',status='old')
c     if(nsize.eq.4)open(20,file='Cu4data',status='old')
c     if(nsize.eq.5)open(20,file='Cu5data',status='old')
c     if(nsize.eq.6)open(20,file='Cu6data',status='old')
c     if(nsize.eq.7)open(20,file='Cu7data',status='old')
c     if(nsize.eq.8)open(20,file='Cu8data',status='old')
c     if(nsize.eq.9)open(20,file='Cu9data',status='old')
c     if(nsize.eq.10)open(20,file='Cu10data',status='old')
c     if(nsize.eq.11)open(20,file='Cu11data',status='old')
c     if(nsize.eq.12)open(20,file='Cu12data',status='old')
c     if(nsize.eq.13)open(20,file='Cu13data',status='old')
c     if(nsize.eq.14)open(20,file='Cu14data',status='old')
c     if(nsize.eq.15)open(20,file='Cu15data',status='old')
c     if(nsize.eq.16)open(20,file='Cu16data',status='old')
c     if(nsize.eq.17)open(20,file='Cu17data',status='old')
c     if(nsize.eq.18)open(20,file='Cu18data',status='old')
c     if(nsize.eq.19)open(20,file='Cu19data',status='old')
c     if(nsize.eq.20)open(20,file='Cu20data',status='old')
c     if(nsize.eq.21)open(20,file='Cu21data',status='old')
c     if(nsize.eq.22)open(20,file='Cu22data',status='old')
c     if(nsize.eq.23)open(20,file='Cu23data',status='old')
c     if(nsize.eq.24)open(20,file='Cu24data',status='old')
c     if(nsize.eq.25)open(20,file='Cu25data',status='old')
c     if(nsize.eq.32)open(20,file='Cu32data',status='old')
c     if(nsize.eq.53)open(20,file='Cu53data',status='old')
c     if(nsize.eq.55)open(20,file='Cu55data',status='old')
c     if(nsize.eq.63)open(20,file='Cu63data',status='old')
c     if(nsize.eq.67)open(20,file='Cu67data',status='old')
c     if(nsize.eq.71)open(20,file='Cu71data',status='old')
c     if(nsize.eq.87)open(20,file='Cu87data',status='old')
c     if(nsize.eq.88)open(20,file='Cu88data',status='old')
c     if(nsize.eq.108)open(20,file='Cu108data',status='old')
c     if(nsize.eq.145)open(20,file='Cu145data',status='old')
c     if(nsize.eq.147)open(20,file='Cu147data',status='old')
c     if(nsize.eq.158)open(20,file='Cu158data',status='old')
c     if(nsize.eq.172)open(20,file='Cu172data',status='old')
c     if(nsize.eq.365)open(20,file='Cu365data',status='old')
c     if(nsize.eq.864)open(20,file='Cu864data',status='old')
c     if(nsize.eq.1372)open(20,file='Cu1372data',status='old')
c     if(nsize.eq.2457)open(20,file='Cu2457data',status='old')
c     if(nsize.eq.4000)open(20,file='Cu4000data',status='old')

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
c     call ffield(nsize,x,y,z,energy)
c     write(10,*)nsize,n,e/nsize,energy
c     if(nsize.gt.25)then
c     call dfield(nsize,x,y,z,dfx,dfy,dfz)
c     do i=1,nsize             
c     xp(i)=x(i)+eps   ! check against finite difference
c     yp(i)=y(i)+eps
c     zp(i)=z(i)+eps
c     call ffield(nsize,xp,y,z,ex)
c     call ffield(nsize,x,yp,z,ey)
c     call ffield(nsize,x,y,zp,ez)
c     write(30,'(2(i2,1x),6(e10.3,1x))')n,i,(ex-energy)/eps,dfx(i),
c    #               (ey-energy)/eps,dfy(i),(ez-energy)/eps,dfz(i)
c     xp(i)=x(i)
c     yp(i)=y(i)
c     zp(i)=z(i)
c     end do
c     endif
c100  continue
c  
c     stop
c     end


c*****************************************************************
c
c  FF code to evaluate analytical derivatives of Q-SC model
c
c234567
      subroutine dfield(nm,x,y,z,dfx,dfy,dfz)
      implicit none
      integer nm,nmax,i,j
      parameter (nmax=4000)
      real*8 x(nm),y(nm),z(nm),dfx(nm),dfy(nm),dfz(nm),
     #       r(nmax,nmax),rho(nmax),eij,cij,aij,rij,rsq,
     #       xij,yij,zij,sum1,sum2,sum3,dvax,dvay,dvaz,
     #       dvrx,dvry,dvrz

      if(nm.gt.nmax)write(*,*)'dimension error in dfield'

      eij=5.35d-3
      cij=84.843d0
      aij=3.603d0

      do i=1,nm
      sum1=0.d0
      do j=1,nm
        rsq=(x(i)-x(j))**2+(y(i)-y(j))**2+(z(i)-z(j))**2
        r(i,j)=sqrt(rsq)
        rij=r(i,j)
        if(j.ne.i)sum1=sum1+(aij/rij)**5
      end do
      rho(i)=dsqrt(sum1)
      end do

      sum1=0.d0
      sum2=0.d0
      sum3=0.d0
      do i=1,nm
      dvrx=0.d0
      dvax=0.d0
      dvry=0.d0
      dvay=0.d0
      dvrz=0.d0
      dvaz=0.d0
      do j=1,nm
        rij=r(i,j)
        xij=x(i)-x(j)
        yij=y(i)-y(j)
        zij=z(i)-z(j)
        if(j.ne.i)then
          dvrx=dvrx-10.d0*aij**10.d0*xij/rij**12.d0
          dvax=dvax-5.d0*aij**5.d0*xij/rij**7.d0/2.d0/rho(i)
     #             -5.d0*aij**5.d0*xij/rij**7.d0/2.d0/rho(j)
          dvry=dvry-10.d0*aij**10.d0*yij/rij**12.d0
          dvay=dvay-5.d0*aij**5.d0*yij/rij**7.d0/2.d0/rho(i)
     #             -5.d0*aij**5.d0*yij/rij**7.d0/2.d0/rho(j)
          dvrz=dvrz-10.d0*aij**10.d0*zij/rij**12.d0
          dvaz=dvaz-5.d0*aij**5.d0*zij/rij**7.d0/2.d0/rho(i)
     #             -5.d0*aij**5.d0*zij/rij**7.d0/2.d0/rho(j)
        endif
      end do
      sum1=eij*dvrx-cij*eij*dvax
      sum2=eij*dvry-cij*eij*dvay
      sum3=eij*dvrz-cij*eij*dvaz
      dfx(i)=sum1/nm
      dfy(i)=sum2/nm
      dfz(i)=sum3/nm
      end do

      return
      end

c*****************************************************************
c
c  FF code to evaluate the energy of a copper cluster
c
c  revised (7/17/08) to include Q-SC model for large clusters
c
c234567
      subroutine ffield(nm,x,y,z,ff)
      implicit none
      integer i,j,k,l,m,n,p,nmax,nm,ierror,lcount,lmax,ncoord
      parameter (nmax=4000)
      real*8 nc(nmax),mc(nmax),lc(nmax)
      real*8 f,x(nm),y(nm),z(nm),r(nmax,nmax),fc(nmax,nmax),
     #       vr,va,de,re,se,beta,rsq,e,fcutoff,rmin,rmax,
     #       costheta,sum,term,arg,pi,a0,a,b,nu,g,h,num,den,
     #       alpha1,delta1,theta1,theta2,theta3,theta4,req,deq,
     #       a1,c1,d1,h1,g1,a2,c2,d2,h2,g2,alpha,arg1,arg2,arg3,
     #       a3,c3,d3,h3,g3,a4,c4,d4,h4,g4,term1,term2,term3,
     #       bsym1,bsym2,b1(nmax,nmax),b2(nmax,nmax),ff,arg4,
     #       alpha2,alpha3,alpha4,delta2,delta3,delta4,term4,
     #       gamma1,gamma2,gamma3,gamma4,gamma5,gamma6,gamma7,
     #       seq,betaq,gc(nmax,nmax),sum1,sum2,sum3,delta0,
     #       sum4,sum5,sum6,sum7,sum8,anew,bnew,cnew,dnew,
     #       defect1,defect2,defect3,defect4,coeff1,fcclimit,
     #       sav1(nm),sav2(nm),sav3(nm),sav4(nm),sav5(nm),
     #       sav6(nm),sav7(nm),sav8(nm),eij,cij,aij,rij,ncutoff
      character*2 atom
      external scale
 
      if(nm.gt.nmax)write(*,*)'dimension error in ffield'

      de=2.15049644d0
      re=2.31115983d0
      se=3.07088775d0
      beta=1.34410778d0
      anew=1.d0
      bnew=1.d0
      cnew=1.d0
      dnew=1.d0

c     fcclimit=3.2d0/6.d0*(1.d0+11.d0*a1)**delta1
c     coeff1=fcclimit-de*exp(20.d0*gamma1)
c     write(*,*)coeff1

      do i=1,nm
      sum=0.d0
      do j=1,nm
        rsq=(x(i)-x(j))**2+(y(i)-y(j))**2+(z(i)-z(j))**2
        r(i,j)=sqrt(rsq)
        fc(i,j)=fcutoff(r(i,j),1)
        sum=sum+fc(i,j)
      end do
      lc(i)=sum-1
      call scale(defect1,defect2,defect3,defect4,
     #              a0,alpha1,delta0,delta1,sum)
      sav1(i)=defect1
      sav2(i)=defect2
      sav3(i)=defect3
      sav4(i)=defect4
      sav5(i)=a0
      sav6(i)=alpha1
      sav7(i)=delta0
      sav8(i)=delta1
      end do
   
      do i=1,nm
      sum1=0.d0
      do j=1,nm
        gc(i,j)=fcutoff(r(i,j),3)
        sum1=sum1+gc(i,j)
      end do
      mc(i)=sum1-1
      end do

      do i=1,nm
      sum1=0.d0
      do j=1,nm
        fc(i,j)=fcutoff(r(i,j),2)
        sum1=sum1+fc(i,j)
      end do
      nc(i)=sum1-1
c     write(55,*)i,lc(i)+1,nc(i)+1,mc(i)+1
      end do

      ncutoff=24
      do i=1,nm
      if(lc(i).gt.ncutoff)go to 100
      a1=sav5(i)
      alpha1=sav6(i)
      delta0=sav7(i)
      delta1=sav8(i)
      do j=1,nm
      if(lc(j).gt.ncutoff)write(*,*)'warning, j=',j
      if(i.eq.j .or. fc(i,j).eq.0.d0)then
      b1(i,j)=0.d0
      else
      term1=0.d0
      do k=1,nm
        if(lc(k).gt.ncutoff)write(*,*)'warning, k=',k
        if(k.ne.i .and. k.ne.j .and. fc(i,k).ne.0.d0)then
        arg1=alpha1*(r(i,j)-r(i,k))
        if(lc(i).lt.0)then
        costheta=(r(i,j)**2+r(j,k)**2-r(i,k)**2)/(2*r(i,j)*r(j,k))
        g1=1.d0+c1**2/d1**2-c1**2/(d1**2+(h1-costheta)**2)
        else
        g1=1.d0
        endif
        else
        g1=0.d0
        endif
        term1=term1+a1*g1*fc(i,k)*exp(arg1)
      end do
      arg2=max(0,min(8,mc(i)+mc(j)-2*gc(i,j)-2.d0))
      b1(i,j)=(1.d0+term1)**(-delta0-delta1*arg2)
      endif
      end do
 100  continue
      end do

      sum=0.d0
      do i=1,nm
      if(lc(i).gt.ncutoff)then
      vr=0.d0
      va=0.d0
c     eij=5.7921d-3
      eij=5.35d-3
      cij=84.843d0
      aij=3.603d0
      do j=1,nm
        rij=r(i,j)
        if(j.eq.i)then
          vr=vr
          va=va
        else
          vr=vr+(aij/rij)**10.d0
          va=va+(aij/rij)**5.d0
        endif
      end do
      sum=sum+eij*vr/2.d0-cij*eij*sqrt(va)
      go to 200
      endif
      g1=sav1(i)
      g2=sav2(i)
      g3=sav3(i)
      g4=sav4(i)
      do j=i+1,nm
        if(lc(j).gt.ncutoff)write(*,*)'warning, j=',j
        bsym1=(b1(i,j)+b1(j,i))/2.d0
        arg1=max(0,min(8,anew*(mc(i)+mc(j)-2*gc(i,j)-2.d0)))
        arg2=max(0,min(8,bnew*(mc(i)+mc(j)-2*gc(i,j)-2.d0)))
        arg3=max(0,min(8,cnew*(mc(i)+mc(j)-2*gc(i,j)-2.d0)))
        arg4=max(0,min(8,dnew*(mc(i)+mc(j)-2*gc(i,j)-2.d0)))
        deq=de+g1*arg1
        req=re+g2*arg2
        seq=se+g3*arg3
        betaq=beta+g4*arg4
        vr=deq/(seq-1.d0)*exp(-(sqrt(2*seq)*betaq*(r(i,j)-req)))
        va=deq*seq/(seq-1.d0)*exp(-(sqrt(2/seq)*betaq*(r(i,j)-req)))
        sum=sum+fc(i,j)*(vr-va)*bsym1
c       write(99,*)r(i,j),fc(i,j)*(vr-va),fc(i,j)*(vr-va)*bsym1
      end do
 200  continue
      end do

      ff=sum/nm

      return
      end
 
c************************************************************
c
c234567
      real*8 function fcutoff(r,index)
      implicit none
      integer index
      real*8 r,rmin,rmax,pi,arg

      pi=dacos(-1.d0)

      if(index.eq.1)then
        rmin=10.d0
        rmax=15.d0
      elseif(index.eq.2)then
        rmin=5.d0
        rmax=7.d0
      elseif(index.eq.3)then
        rmin=2.7d0
        rmax=3.2d0
      else
        write(*,*)'index error in fcutoff'
      endif

      if(r.lt.rmin)then
        fcutoff=1.d0
      elseif(r.lt.rmax)then
        arg=pi*(r-rmin)/(rmax-rmin)
        fcutoff=(1.d0+cos(arg))/2.d0 
      else
        fcutoff=0.d0
      endif

      return
      end

c************************************************************
c
c234567
      subroutine scale(shift1,shift2,shift3,shift4,
     #                   a1,alpha1,delta0,delta1,n)
      implicit none
      real*8 shift1,shift2,shift3,shift4,a1,alpha1,
     #       delta0,delta1,n
 
      if(n.le.2)then
      a1=1.d0
      shift1=0.d0
      shift2=0.d0
      shift3=0.d0
      shift4=0.d0
      alpha1=0.d0
      delta1=0.d0
      delta0=0.d0
      elseif(n.gt.2 .and. n.le.3)then
      a1=1.069328d0
      shift1=0.d0
      shift2=0.d0
      shift3=0.d0
      shift4=0.d0
      alpha1=13.597996d0
      delta1=0.d0
      delta0=0.911726d0
      elseif(n.gt.3 .and. n.le.4)then
      a1=0.648487d0
      shift1=1.460445d0
      shift2=0.128394d0
      shift3=9.946546d0
      shift4=-0.274245d0
      alpha1=1.015724d0
      delta1=0.645296d0
      delta0=0.911726d0
      elseif(n.gt.4 .and. n.le.5)then
      a1=0.503753d0
      shift1=0.281321d0
      shift2=0.019378d0
      shift3=2.555847d0
      shift4=0.278186d0
      alpha1=1.504305d0
      delta1=0.177818d0
      delta0=0.911726d0
      elseif(n.gt.5 .and. n.le.6)then
      a1=0.932259d0
      shift1=1.052853d0
      shift2=0.030376d0
      shift3=-0.197126d0
      shift4=0.042575d0
      alpha1=1.661175d0
      delta1=0.143196d0
      delta0=0.911726d0
      elseif(n.gt.6 .and. n.le.7)then
      a1=2.654231d0
      shift1=0.649477d0
      shift2=0.032205d0
      shift3=-0.034713d0
      shift4=-0.013672d0
      alpha1=2.294211d0
      delta1=0.049695d0
      delta0=0.55d0
      elseif(n.gt.7 .and. n.le.8)then
      a1=2.289658d0
      shift1=0.661942d0
      shift2=0.021854d0
      shift3=-0.023934d0
      shift4=0.006514d0
      alpha1=2.308261d0
      delta1=0.054622d0
      delta0=0.55d0
      elseif(n.gt.8 .and. n.le.9)then
      a1=2.307554d0
      shift1=0.459070d0
      shift2=0.034354d0
      shift3=-0.034298d0
      shift4=0.024123d0
      alpha1=1.898695d0
      delta1=0.037339d0
      delta0=0.55d0
      elseif(n.gt.9 .and. n.le.10)then
      a1=2.308585d0
      shift1=0.457422d0
      shift2=0.029284d0
      shift3=-0.034548d0
      shift4=0.030815d0
      alpha1=1.898865d0
      delta1=0.035043d0
      delta0=0.55d0
      elseif(n.gt.10 .and. n.le.11)then
      a1=1.347445d0
      shift1=0.193441d0
      shift2=0.027390d0
      shift3=-0.033655d0
      shift4=0.027374d0
      alpha1=1.648295d0
      delta1=0.031680d0
      delta0=0.55d0
      elseif(n.gt.11 .and. n.le.12)then
      a1=1.347445d0
      shift1=0.203000d0
      shift2=0.027390d0
      shift3=-0.033655d0
      shift4=0.027374d0
      alpha1=1.648295d0
      delta1=0.031680d0
      delta0=0.55d0
      elseif(n.gt.12 .and. n.le.13)then
      a1=1.052998d0
      shift1=0.240226d0
      shift2=0.026274d0
      shift3=0.037879d0
      shift4=0.022774d0
      alpha1=1.722316d0
      delta1=0.050567d0
      delta0=0.55d0
      elseif(n.gt.13 .and. n.le.14)then
      a1=0.980241d0
      shift1=0.143588d0
      shift2=0.029708d0
      shift3=0.179160d0
      shift4=0.030420d0
      alpha1=4.411814d0
      delta1=0.040847d0
      delta0=0.55d0
      elseif(n.gt.14 .and. n.le.18)then
      a1=1.d0
      shift1=-0.02d0
      shift2=0.0299d0
      shift3=0.d0
      shift4=0.d0
      alpha1=1.d0
      delta1=0.d0
      delta0=0.55d0
      elseif(n.gt.18 .and. n.le.25)then
      a1=1.d0
c     shift1=-0.03d0   ! best for 19
c     shift1=-0.029d0  ! best for 25
      shift1=-0.028d0  ! best for 20-24
      shift2=0.0299d0
      shift3=0.d0
      shift4=0.d0
      alpha1=3.d0
      delta1=0.d0
      delta0=0.55d0
      elseif(n.gt.25)then
      a1=1.d0
      shift1=-0.015d0
      shift2=0.0299d0
      shift3=0.d0
      shift4=0.d0
      alpha1=3.d0
      delta1=0.d0
      delta0=0.55d0
      endif

      return
      end

