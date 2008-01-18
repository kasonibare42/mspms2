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
c  test  - tests di-clusters to ensure proper break-up 
c
c234567
c     implicit none
c     integer i,n,nsize,ndata
c     parameter (nsize=20,ndata=217)
c     real*8 x(2*nsize),y(2*nsize),z(2*nsize),energy,e
c     character*2 atom
c     logical test
c     parameter (test=.true.)

c     open(10,file='ffield.out',status='unknown')

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
c     if(test)then
c     do i=1,nsize
c       read(20,*)atom,x(i),y(i),z(i)
c       x(nsize+i)=x(i)+40.d0
c       y(nsize+i)=y(i)
c       z(nsize+i)=z(i)
c     end do
c     read(20,*)
c     call ffield(2*nsize,ndata,x,y,z,energy)
c     write(10,*)2*nsize,n,e/nsize,energy
c     else
c     do i=1,nsize
c       read(20,*)atom,x(i),y(i),z(i)
c     end do
c     read(20,*)
c     call ffield(nsize,ndata,x,y,z,energy)
c     write(10,*)nsize,n,e/nsize,energy
c     endif
c100  continue
c  
c     stop
c     end


c*****************************************************************
c
c  FF code to evaluate the energy of one atom in a copper cluster
c
c234567
      subroutine ffieldCu_Single(nm,ndata,x,y,z,ff)
      implicit none
      integer i,j,k,l,m,n,p,nmax,nm,ierror,lcount,lmax,ndata
      parameter (nmax=4000)
      real*8 nc(nmax),mc(nmax)
      real*8 f,x(*),y(*),z(*),r(nmax,nmax),fc(nmax,nmax),
     #       vr,va,de,re,se,beta,rsq,e,fcutoff,rmin,rmax,
     #       costheta,sum,term,arg,pi,a0,a,b,nu,g,h,num,den,
     #       alpha1,delta1,theta1,theta2,theta3,theta4,req,deq,
     #       a1,c1,d1,h1,g1,a2,c2,d2,h2,g2,alpha,arg1,arg2,arg3,
     #       a3,c3,d3,h3,g3,a4,c4,d4,h4,g4,term1,term2,term3,
     #       bsym1,bsym2,b1(nmax,nmax),b2(nmax,nmax),ff,arg4,
     #       alpha2,alpha3,alpha4,delta2,delta3,delta4,term4,
     #       gamma1,gamma2,gamma3,gamma4,gamma5,gamma6,gamma7,
     #       seq,betaq,gc(nmax,nmax),gcutoff,sum1,sum2,
     #       anew,bnew,cnew,dnew,enew,shift,defect
      character*2 atom
 
      if(nm.gt.nmax)then
        write(*,*)'dimension error in ffield'
        return
      endif

      de=2.15049644d0
      re=2.31115983d0
      se=3.07088775d0
      beta=1.34410778d0
      alpha1=13.597996d0
      delta1=0.911726d0
      a1=1.069328d0
      c1=0.000000d0
      d1=0.100000d0
      h1=0.866030d0
      gamma1=0.024585d0
      gamma2=0.004134d0
      gamma3=0.008292d0
      gamma4=0.000761d0
      gamma5=-0.007840d0
      gamma6=0.496589d0
      gamma7=-0.001215d0

      sum1=0.d0
      do i=1,nm
      sum2=0.d0
      do j=1,nm
        rsq=(x(i)-x(j))**2+(y(i)-y(j))**2+(z(i)-z(j))**2
        r(i,j)=sqrt(rsq)
        fc(i,j)=fcutoff(r(i,j),1)
        sum2=sum2+fc(i,j)
      end do
      nc(i)=sum2-1
      sum1=sum1+shift(sum2)
c     write(88,*)i,sum2,shift(sum2)
      end do
      defect=sum1

      do i=1,nm
      sum1=0.d0
      do j=1,nm
        gc(i,j)=gcutoff(r(i,j),nc(i)+1)
        sum1=sum1+gc(i,j)
      end do
      mc(i)=sum1-1
c     write(99,*)i,nc(i)+1,mc(i)+1,defect/nm
      end do

      do i=1,nm
      sum1=0.d0
      do j=1,nm
        fc(i,j)=fcutoff(r(i,j),2)
        sum1=sum1+fc(i,j)
      end do
      nc(i)=sum1-1
      end do

      do i=1,nm
      if(nc(i).le.11)then
      enew=0.716359395
      else
      enew=1.d0
      endif
      do j=1,nm
      if(i.eq.j .or. fc(i,j).eq.0.d0)then
      b1(i,j)=0.d0
      else
      term1=0.d0
      do k=1,nm
        if(k.ne.i .and. k.ne.j .and. fc(i,k).ne.0.d0)then
c       arg=nc(k)-fc(i,k)-fc(j,k)
        arg=mc(k)-gc(i,k)-gc(j,k)
        arg1=gamma5*arg+(alpha1+gamma6*arg)*(r(i,j)-r(i,k))
        g1=a1*fc(i,k)*exp(arg1)
        else
        g1=0.d0
        endif
        term1=term1+g1
      end do
c     arg=min(20,enew*(nc(i)+nc(j)-2*fc(i,j)-2))
      arg=min(20,enew*(mc(i)+mc(j)-2*gc(i,j)-2))
      delta2=delta1*exp(gamma7*arg)
      b1(i,j)=(1.d0+term1)**(-delta2)
      endif
      end do
      end do

      sum=0.d0
c     for only the atom with index ndata
      do i=ndata,ndata
      if(nc(i).le.11)then
      anew=1.87897233
      bnew=0.945024341
      cnew=1.46407271
      dnew=1.39023338
      else
      anew=1.d0
      bnew=1.d0
      cnew=1.d0
      dnew=1.d0
      endif
      do j=1,nm
c     do NOT self calculation
        if (i .ne. j) then
        bsym1=(b1(i,j)+b1(j,i))/2.d0
        arg1=min(20,anew*(mc(i)+mc(j)-2*gc(i,j)-2))
        arg2=min(20,bnew*(mc(i)+mc(j)-2*gc(i,j)-2))
        arg3=min(20,cnew*(mc(i)+mc(j)-2*gc(i,j)-2))
        arg4=min(20,dnew*(mc(i)+mc(j)-2*gc(i,j)-2))
        deq=de*exp(gamma1*arg1)
        req=re*exp(gamma2*arg2)
        seq=se*exp(gamma3*arg3)
        betaq=beta*exp(gamma4*arg4)
        vr=deq/(seq-1.d0)*exp(-(sqrt(2*seq)*betaq*(r(i,j)-req)))
        va=deq*seq/(seq-1.d0)*exp(-(sqrt(2/seq)*betaq*(r(i,j)-req)))
        sum=sum+fc(i,j)*(vr-va)*bsym1 
        endif
      end do
      end do

c     return the energy for the atom
      ff=sum-defect/nm

      return
      end
 







c*****************************************************************
c
c  FF code to evaluate the energy of a copper cluster
c
c234567
      subroutine ffieldCu(nm,ndata,x,y,z,ff)
      implicit none
      integer i,j,k,l,m,n,p,nmax,nm,ierror,lcount,lmax,ndata
      parameter (nmax=4000)
      real*8 nc(nmax),mc(nmax)
      real*8 f,x(*),y(*),z(*),r(nmax,nmax),fc(nmax,nmax),
     #       vr,va,de,re,se,beta,rsq,e,fcutoff,rmin,rmax,
     #       costheta,sum,term,arg,pi,a0,a,b,nu,g,h,num,den,
     #       alpha1,delta1,theta1,theta2,theta3,theta4,req,deq,
     #       a1,c1,d1,h1,g1,a2,c2,d2,h2,g2,alpha,arg1,arg2,arg3,
     #       a3,c3,d3,h3,g3,a4,c4,d4,h4,g4,term1,term2,term3,
     #       bsym1,bsym2,b1(nmax,nmax),b2(nmax,nmax),ff,arg4,
     #       alpha2,alpha3,alpha4,delta2,delta3,delta4,term4,
     #       gamma1,gamma2,gamma3,gamma4,gamma5,gamma6,gamma7,
     #       seq,betaq,gc(nmax,nmax),gcutoff,sum1,sum2,
     #       anew,bnew,cnew,dnew,enew,shift,defect
      character*2 atom
 
      if(nm.gt.nmax)then
        write(*,*)'dimension error in ffield'
        return
      endif

      de=2.15049644d0
      re=2.31115983d0
      se=3.07088775d0
      beta=1.34410778d0
      alpha1=13.597996d0
      delta1=0.911726d0
      a1=1.069328d0
      c1=0.000000d0
      d1=0.100000d0
      h1=0.866030d0
      gamma1=0.024585d0
      gamma2=0.004134d0
      gamma3=0.008292d0
      gamma4=0.000761d0
      gamma5=-0.007840d0
      gamma6=0.496589d0
      gamma7=-0.001215d0

      sum1=0.d0
      do i=1,nm
      sum2=0.d0
      do j=1,nm
        rsq=(x(i)-x(j))**2+(y(i)-y(j))**2+(z(i)-z(j))**2
        r(i,j)=sqrt(rsq)
        fc(i,j)=fcutoff(r(i,j),1)
        sum2=sum2+fc(i,j)
      end do
      nc(i)=sum2-1
      sum1=sum1+shift(sum2)
c     write(88,*)i,sum2,shift(sum2)
      end do
      defect=sum1

      do i=1,nm
      sum1=0.d0
      do j=1,nm
        gc(i,j)=gcutoff(r(i,j),nc(i)+1)
        sum1=sum1+gc(i,j)
      end do
      mc(i)=sum1-1
c     write(99,*)i,nc(i)+1,mc(i)+1,defect/nm
      end do

      do i=1,nm
      sum1=0.d0
      do j=1,nm
        fc(i,j)=fcutoff(r(i,j),2)
        sum1=sum1+fc(i,j)
      end do
      nc(i)=sum1-1
      end do

      do i=1,nm
      if(nc(i).le.11)then
      enew=0.716359395
      else
      enew=1.d0
      endif
      do j=1,nm
      if(i.eq.j .or. fc(i,j).eq.0.d0)then
      b1(i,j)=0.d0
      else
      term1=0.d0
      do k=1,nm
        if(k.ne.i .and. k.ne.j .and. fc(i,k).ne.0.d0)then
c       arg=nc(k)-fc(i,k)-fc(j,k)
        arg=mc(k)-gc(i,k)-gc(j,k)
        arg1=gamma5*arg+(alpha1+gamma6*arg)*(r(i,j)-r(i,k))
        g1=a1*fc(i,k)*exp(arg1)
        else
        g1=0.d0
        endif
        term1=term1+g1
      end do
c     arg=min(20,enew*(nc(i)+nc(j)-2*fc(i,j)-2))
      arg=min(20,enew*(mc(i)+mc(j)-2*gc(i,j)-2))
      delta2=delta1*exp(gamma7*arg)
      b1(i,j)=(1.d0+term1)**(-delta2)
      endif
      end do
      end do

      sum=0.d0
      do i=1,nm-1
      if(nc(i).le.11)then
      anew=1.87897233
      bnew=0.945024341
      cnew=1.46407271
      dnew=1.39023338
      else
      anew=1.d0
      bnew=1.d0
      cnew=1.d0
      dnew=1.d0
      endif
      do j=i+1,nm
        bsym1=(b1(i,j)+b1(j,i))/2.d0
        arg1=min(20,anew*(mc(i)+mc(j)-2*gc(i,j)-2))
        arg2=min(20,bnew*(mc(i)+mc(j)-2*gc(i,j)-2))
        arg3=min(20,cnew*(mc(i)+mc(j)-2*gc(i,j)-2))
        arg4=min(20,dnew*(mc(i)+mc(j)-2*gc(i,j)-2))
        deq=de*exp(gamma1*arg1)
        req=re*exp(gamma2*arg2)
        seq=se*exp(gamma3*arg3)
        betaq=beta*exp(gamma4*arg4)
        vr=deq/(seq-1.d0)*exp(-(sqrt(2*seq)*betaq*(r(i,j)-req)))
        va=deq*seq/(seq-1.d0)*exp(-(sqrt(2/seq)*betaq*(r(i,j)-req)))
        sum=sum+fc(i,j)*(vr-va)*bsym1
      end do
      end do

      ff=(sum-defect)/nm

      return
      end
 
c************************************************************
c
c234567
      real*8 function fcutoff(r,index)
      implicit none
      integer index
      real*8 r,rmin,rmax,pi,arg

      if(index.eq.1)then
        rmin=10.d0
        rmax=15.d0
      else
        rmin=5.d0
        rmax=7.d0
      endif

      pi=dacos(-1.d0)

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
      real*8 function gcutoff(r,n)
      implicit none
c     integer n
      real*8 n,r,rmin,rmax,pi,arg,amin,amax,rmin19,rmax19

      if(n.le.13)then
        rmin=5.d0
        rmax=7.d0
      else
        amin=0.1d0
        amax=0.1d0
        rmin=2.7d0+2.3d0*exp(-amin*(n-13))
        rmax=3.2d0+3.8d0*exp(-amax*(n-13))
      endif

      pi=dacos(-1.d0)

      if(r.lt.rmin)then
        gcutoff=1.d0
      elseif(r.lt.rmax)then
        arg=pi*(r-rmin)/(rmax-rmin)
        gcutoff=(1.d0+cos(arg))/2.d0
      else
        gcutoff=0.d0
      endif

      return
      end

c************************************************************
c
c234567
      real*8 function shift(n)
      implicit none
c     integer n
      logical fit   
      parameter (fit=.false.)
      real*8 n,alpha

      if(n.le.13)then
        shift=0.d0
      elseif(n.gt.13 .and. n.le.14)then
        shift=0.02d0
      elseif(n.gt.14 .and. n.le.15)then
        shift=0.03d0
      elseif(n.gt.15 .and. n.le.16)then
        shift=0.20d0
      elseif(n.gt.16 .and. n.le.17)then
        shift=0.22d0
      elseif(n.gt.17 .and. n.le.18)then
        shift=0.28d0
      elseif(n.gt.18 .and. n.le.19)then
        shift=0.34d0
      elseif(n.gt.19 .and. n.le.20)then
        shift=0.38d0
      elseif(n.gt.20 .and. n.le.21)then
        shift=0.43d0
      elseif(n.gt.21 .and. n.le.22)then
        shift=0.49d0
      elseif(n.gt.22 .and. n.le.23)then
        shift=0.55d0
      elseif(n.gt.23 .and. n.le.24)then
        shift=0.60d0
      elseif(n.gt.24 .and. n.le.25)then
        shift=0.50d0
      else
        alpha=0.08d0
        shift=0.9d0*(1.d0-exp(-alpha*(n-13)))
      endif

      if(fit)then
      if(n.le.13)then
        shift=0.d0
      else
        alpha=0.08d0
        shift=0.9d0*(1.d0-exp(-alpha*(n-13)))
      endif
      endif

      return
      end
