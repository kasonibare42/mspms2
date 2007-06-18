!########################################
!### PERIODIC BOUNDARY CONDITIONS      ##
!### AND MINIMUM IMAGE CONVENTION      ##
!########################################
! PBCPAR = 1 ORTHORHOMBIC/CUBIC
! PBCPAR = 3 MONOCLINIC B
! PBCPAR = 6 TRICLINIC

Subroutine PBC(xo,yo,zo,xt,yt,zt)
  Use Global
  Implicit None
  Real(8) :: xo,yo,zo,xt,yt,zt
  Integer :: k,l,m
  !----------------------------
  Select Case(PBCType)
  Case(1)
     !# ORTHORHOMBIC/CUBIC
     !# pbcpar(1...3) = (Lx,Ly,Lz)
     xt = xo - floor(xo/pbcpar(1))*pbcpar(1)
     yt = yo - floor(yo/pbcpar(2))*pbcpar(2)
     zt = zo - floor(zo/pbcpar(3))*pbcpar(3)
  Case(2)
     !# NANOTUBE
     !# pbcpar(1...4) = (Lx,Ly,Lz,Tube Radius
     !# periodic boundary conditions only along the tube
     xt = xo
     yt = yo
     zt = zo - floor(zo/pbcpar(3))*pbcpar(3)
  Case(3)
     !# MONOCLINIC B - Angle on XZ plane - Y Direction orthorhombic PBC 
     !# C Directed along Z
     !------------------
     !  Lx -> pbcpar(1)
     !  Ly -> pbcpar(2)
     !  Lz -> pbcpar(3)
     ! DRx -> pbcpar(4)
     ! using floor because my box is from 0 to 1.
     k=FLOOR(xo/pbcpar(1))
     xt = xo  - k*pbcpar(1)
     yt = yo  - FLOOR(yo/pbcpar(2))*pbcpar(2)
     zt = zo  - k*pbcpar(4)
     zt = zt - FLOOR(zt/pbcpar(3))*pbcpar(3)
  Case(6)
     !# TRICLINIC - C PARALLEL TO Z , BC ON XY PLANE
     ! Auc_new -> 1 
     ! Buc_new -> 2
     ! Cuc_new -> 3
     ! ax -> 4
     ! ay -> 5
     ! az -> 6          
     ! ay/ax -> 7
     ! wz -> 8
     ! Auc*ay -> 9
     ! Auc*az -> 10
     ! Buc*Cos[phi_a] -> 11
     ! Cot[Phi_a] -> 12
     ! 1/Auc_new -> 13
     ! 1/Buc_new -> 14
     ! 1/Cuc_new -> 15
     k = FLOOR(PBCpar(13)*xo)
     l = FLOOR(pbcpar(14)*(yo - pbcpar(7)*xo))
     m = FLOOR(pbcpar(15)*(zo - pbcpar(12)*yo + pbcpar(8)*xo))
     xt = xo - k*pbcpar(1)
     yt = yo - k*pbcpar(9) - l*pbcpar(2)
     zt = zo - k*pbcpar(10) - l*pbcpar(11) - l*pbcpar(3)
  End Select
End Subroutine PBC

Subroutine PBCuc(xo,yo,zo,xt,yt,zt)
  Use Global
  Implicit None
  Real(8) :: xo,yo,zo,xt,yt,zt,x2,y2,z2,xn,yn,zn
  Integer :: k,l,m,n
  !----------------------------
  Select Case(PBCType)
  Case(1)
     !# ORTHORHOMBIC
     !# pbcpar_uc(1...3) = (A,B,C)
     xt = xo - floor(xo/pbcpar_uc(1))*pbcpar_uc(1)
     yt = yo - floor(yo/pbcpar_uc(2))*pbcpar_uc(2)
     zt = zo - floor(zo/pbcpar_uc(3))*pbcpar_uc(3)
  Case(2)
     !# NANOTUBE
     xt = xo
     yt = yo
     zt = zo - floor(zo/pbcpar_uc(3))*pbcpar_uc(3)
  Case(3)
     !# MONOCLINIC B - Angle on XZ plane - Y Direction orthorhombic PBC 
     !------------------
     !  A -> pbcpar_uc(1)
     !  B -> pbcpar_uc(2)
     !  C -> pbcpar_uc(3)
     ! DRx -> pbcpar_uc(4)
     ! using floor2 because my box is from 0 to 1.
     k=FLOOR(xo/pbcpar_uc(1))
     xt = xo  - k*pbcpar_uc(1)
     yt = yo  - FLOOR(yo/pbcpar_uc(2))*pbcpar_uc(2)
     zt = zo  - k*pbcpar_uc(4)
     zt = zt - FLOOR(zt/pbcpar_uc(3))*pbcpar_uc(3)
  Case(6)
     !# TRICLINIC - C PARALLEL TO Z , BC ON XY PLANE
     ! Auc_new -> 1 
     ! Buc_new -> 2
     ! Cuc_new -> 3
     ! ax -> 4
     ! ay -> 5
     ! az -> 6          
     ! ay/ax -> 7
     ! wz -> 8
     ! Auc*ay -> 9
     ! Auc*az -> 10
     ! Buc*Cos[phi_a] -> 11
     ! Cot[Phi_a] -> 12
     ! 1/Auc_new -> 13
     ! 1/Buc_new -> 14
     ! 1/Cuc_new -> 15
     k = FLOOR(PBCpar_uc(13)*xo)
     l = FLOOR(pbcpar_uc(14)*(yo - pbcpar_uc(7)*xo))
     m = FLOOR(pbcpar_uc(15)*(zo - pbcpar_uc(12)*yo + pbcpar_uc(8)*xo))
     xt = xo - k*pbcpar_uc(1)
     yt = yo - k*pbcpar_uc(9) - l*pbcpar_uc(2)
     zt = zo - k*pbcpar_uc(10) - l*pbcpar_uc(11) - m*pbcpar_uc(3)
     m=floor(yt/hy+1)
     n=floor(zt/hz+1)
     if (m .le. 0 .or. n .le. 0) then
     If (M .Le. 0 .AND. N.GT.0) then
                 Call INVPBC(Xt,Yt,Zt,0,1,0,Xn,Yn,Zn)
                 If (Zn .LT. 0d0) then 
                    Call INVPBC(Xn,Yn,Zn,0,0,1,X2,Y2,Z2)  
                    Xn = X2
                    Yn = Y2
                    Zn = Z2
                 EndIf
              ElseIf (N .Le. 0 .AND. M .GT. 0) then
                 Call INVPBC(Xt,Yt,Zt,0,0,1,Xn,Yn,Zn)  
              Else
                 Call INVPBC(Xt,Yt,Zt,0,1,1,Xn,Yn,Zn)
              EndIf
	 xt = xn 
         yt = yn
         zt = zn
   endif
   if (zt .lt. 0d0) Write(*,*) 'BLAKA....',Xt,Yt,Zt
   End Select
End Subroutine PBCuc
!################  

Subroutine InvPBC(xo,yo,zo,i,j,k,xt,yt,zt)
  Use Global
  Implicit None
  Real(8) :: xo,yo,zo,xt,yt,zt
  Integer :: i,j,k
  !-------------------
  Select Case(PBCType)
  Case(1)
     !# ORTHORHOMBIC
     !# pbcpar_uc(1...3) = (Auc,Buc,Cuc)
     xt = xo + i*pbcpar_uc(1)
     yt = yo + j*pbcpar_uc(2)
     zt = zo + k*pbcpar_uc(3)
  Case(2)
     !# NANOTUBE
     xt = xo 
     yt = yo 
     zt = zo + k*pbcpar_uc(3)
  Case(3)
     !# MONOCLINIC B - Angle on XZ plane - Y Direction orthorhombic PBC 
     !# C Directed along Z
     !--------------------
     !  Auc*Dsin(phi_uc) -> pbcpar_uc(1)
     !  Buc -> pbcpar_uc(2)
     !  Cuc -> pbcpar_uc(3)
     ! DRx -> pbcpar_uc(4)
     xt = xo + i*pbcpar_uc(1)  
     yt = yo + j*pbcpar_uc(2)
     zt = zo + i*pbcpar_uc(4) + k*pbcpar_uc(3) 
  Case(6)             
     !# TRICLINIC - C PARALLEL TO Z , BC ON XY PLANE
     ! Auc_new -> 1 
     ! Buc_new -> 2
     ! Cuc_new -> 3
     ! ax -> 4
     ! ay -> 5
     ! az -> 6          
     ! ay/ax -> 7
     ! wz -> 8
     ! Auc*ay -> 9
     ! Auc*az -> 10
     ! Buc*Cos[phi_a] -> 11
     ! Cot[Phi_a] -> 12
     ! 1/Auc_new -> 13
     ! 1/Buc_new -> 14
     ! 1/Cuc_new -> 15
     xt = xo + i*pbcpar_uc(1)
     yt = yo + i*pbcpar_uc(9) + j*pbcpar_uc(2)
     zt = zo + i*pbcpar_uc(10) + j*pbcpar_uc(11) + k*pbcpar_uc(3)
  End Select
End Subroutine InvPBC

!################  

Subroutine MinImage(xo,yo,zo,xt,yt,zt)
  Use Global
  Implicit None
  Real(8) :: xo,yo,zo,xt,yt,zt
  Integer :: k,l,m
  !----------------------------
  Select Case(PBCType)
  Case(1)
     !# ORTHORHOMBIC/CUBIC
     !# pbcpar(1...3) = (Lx,Ly,Lz)
     xt = xo - anint(xo/pbcpar(1))*pbcpar(1)
     yt = yo - anint(yo/pbcpar(2))*pbcpar(2)
     zt = zo - anint(zo/pbcpar(3))*pbcpar(3)
  Case(2)
     !# NANOTUBE
     xt = xo 
     yt = yo 
     zt = zo - anint(zo/pbcpar(3))*pbcpar(3)
  Case(3)                     !___
     !# MONOCLINIC B    Form: /__/ (x is horizontal, z diagonal up and right) 
     !  Lx -> pbcpar(1)
     !  Ly -> pbcpar(2)
     !  Lz -> pbcpar(3)
     ! DRx -> pbcpar(4)
     k=ANINT(xo/pbcpar(1))
     xt = xo  - k*pbcpar(1)
     zt = zo  - k*pbcpar(4)
     zt = zt - ANINT(zt/pbcpar(3))*pbcpar(3)
     yt = yo  - ANINT(yo/pbcpar(2))*pbcpar(2)     
  Case(6)     
     !# TRICLINIC - C PARALLEL TO Z , BC ON XY PLANE
     ! Auc_new -> 1 
     ! Buc_new -> 2
     ! Cuc_new -> 3
     ! ax -> 4
     ! ay -> 5
     ! az -> 6          
     ! ay/ax -> 7
     ! wz -> 8
     ! Auc*ay -> 9
     ! Auc*az -> 10
     ! Buc*Cos[phi_a] -> 11
     ! Cot[Phi_a] -> 12
     ! 1/Auc_new -> 13
     ! 1/Buc_new -> 14
     ! 1/Cuc_new -> 15
     k = ANINT(pbcpar(13)*xo)
     l = ANINT(pbcpar(14)*(yo - pbcpar(7)*xo))
     m = ANINT(pbcpar(15)*(zo - pbcpar(12)*yo + pbcpar(8)*xo))
     xt = xo - k*pbcpar(1)
     yt = yo - k*pbcpar(9) - l*pbcpar(2)
     zt = zo - k*pbcpar(10) - l*pbcpar(11) - m*pbcpar(3)
  End Select
End Subroutine MinImage

Subroutine MinImage_uc(xo,yo,zo,xt,yt,zt)
  Use Global
  Implicit None
  Real(8) :: xo,yo,zo,xt,yt,zt
  Integer :: k,l,m
  !----------------------------
  Select Case(PBCType)
  Case(1)
     !# ORTHORHOMBIC/CUBIC
     !# pbcpar(1...3) = (Lx,Ly,Lz)
     xt = xo - anint(xo/pbcpar_uc(1))*pbcpar_uc(1)
     yt = yo - anint(yo/pbcpar_uc(2))*pbcpar_uc(2)
     zt = zo - anint(zo/pbcpar_uc(3))*pbcpar_uc(3)
  Case(2)
     !# NANOTUBE
     xt = xo
     yt = yo
     zt = zo - anint(zo/pbcpar_uc(3))*pbcpar_uc(3)
  Case(3)                     !___
     !# MONOCLINIC B    Form: /__/ (x is horizontal, z diagonal up and right) 
     !  Lx -> pbcpar(1)
     !  Ly -> pbcpar(2)
     !  Lz -> pbcpar(3)
     ! DRx -> pbcpar(4)
     k=ANINT(xo/pbcpar_uc(1))
     xt = xo  - k*pbcpar_uc(1)
     zt = zo  - k*pbcpar_uc(4)
     zt = zt - ANINT(zt/pbcpar_uc(3))*pbcpar_uc(3)
     yt = yo  - ANINT(yo/pbcpar_uc(2))*pbcpar_uc(2)     
  Case(6)     
     !# TRICLINIC - C PARALLEL TO Z , BC ON XY PLANE
     ! Auc_new -> 1 
     ! Buc_new -> 2
     ! Cuc_new -> 3
     ! ax -> 4
     ! ay -> 5
     ! az -> 6          
     ! ay/ax -> 7
     ! wz -> 8
     ! Auc*ay -> 9
     ! Auc*az -> 10
     ! Buc*Cos[phi_a] -> 11
     ! Cot[Phi_a] -> 12
     ! 1/Auc_new -> 13
     ! 1/Buc_new -> 14
     ! 1/Cuc_new -> 15
     k = ANINT(pbcpar_uc(13)*xo)
     l = ANINT(pbcpar_uc(14)*(yo - pbcpar_uc(7)*xo))
     m = ANINT(pbcpar_uc(15)*(zo - pbcpar_uc(12)*yo + pbcpar_uc(8)*xo))
     xt = xo - k*pbcpar_uc(1)
     yt = yo - k*pbcpar_uc(9) - l*pbcpar_uc(2)
     zt = zo - k*pbcpar_uc(10) - l*pbcpar_uc(11) - m*pbcpar_uc(3)
  End Select
End Subroutine MinImage_uc

Subroutine MinImage_Crystal(xo,yo,zo,xt,yt,zt)
  Use Global
  Implicit None
  Real(8) ,INTENT(IN) :: xo,yo,zo
  Real(8) ,INTENT(OUT) :: xt,yt,zt
  Integer :: k,l,m
  !----------------------------
  Select Case(PBCType)
  Case(1)
     !# ORTHORHOMBIC
     !# pbcpar_crystal(1...3) = (Lxt,Lyt,Lzt)
     xt = xo - anint(xo/pbcpar_crystal(1))*pbcpar_crystal(1)
     yt = yo - anint(yo/pbcpar_crystal(2))*pbcpar_crystal(2)
     zt = zo - anint(zo/pbcpar_crystal(3))*pbcpar_crystal(3)
  Case(2)
     !# NANOTUBE
     xt = xo 
     yt = yo 
     zt = zo - anint(zo/pbcpar_crystal(3))*pbcpar_crystal(3)
  Case(3)                     !___
     !# MONOCLINIC B    Form: /__/ (x is horizontal, z diagonal up and right) 
     !  Lxt -> pbcpar_crystal(1)
     !  Lyt -> pbcpar_crystal(2)
     !  Lzt -> pbcpar_crystal(3)
     ! DRxt -> pbcpar_crystal(4)
     k=ANINT(xo/pbcpar_crystal(1))
     xt = xo  - k*pbcpar_crystal(1)
     zt = zo  - k*pbcpar_crystal(4)
     zt = zt - ANINT(zt/pbcpar_crystal(3))*pbcpar_crystal(3)
     yt = yo  - ANINT(yo/pbcpar_crystal(2))*pbcpar_crystal(2)
  Case(6)
     !# TRICLINIC - C PARALLEL TO Z , BC ON XY PLANE
     ! Auc_new -> 1 
     ! Buc_new -> 2
     ! Cuc_new -> 3
     ! ax -> 4
     ! ay -> 5
     ! az -> 6          
     ! ay/ax -> 7
     ! wz -> 8
     ! Auc*ay -> 9
     ! Auc*az -> 10
     ! Buc*Cos[phi_a] -> 11
     ! Cot[Phi_a] -> 12
     ! 1/Auc_new -> 13
     ! 1/Buc_new -> 14
     ! 1/Cuc_new -> 15
     k = ANINT(pbcpar_crystal(13)*xo)
     l = ANINT(pbcpar_crystal(14)*(yo - pbcpar_crystal(7)*xo))
     m = ANINT(pbcpar_crystal(15)*(zo - pbcpar_crystal(12)*yo + pbcpar_crystal(8)*xo))
     xt = xo - k*pbcpar_crystal(1)
     yt = yo - k*pbcpar_crystal(9) - l*pbcpar_crystal(2)
     zt = zo - k*pbcpar_crystal(10) - l*pbcpar_crystal(11) - m*pbcpar_crystal(3)
  End Select
End Subroutine MinImage_Crystal

Subroutine ABCtoXYZ(At,Bt,Ct,Xt,Yt,Zt)
  Use Global
  Implicit none
 Real(8) :: At,Bt,Ct,Xt,Yt,Zt
 Select Case(PBCType)
  Case(1)
     !# ORTHORHOMBIC
     Xt = At
     Yt = Bt
     Zt = Ct
  Case(2)
     !# NANOTUBE
     Xt = At
     Yt = Bt
     Zt = Ct
  Case(3)                     !___
     !# MONOCLINIC B    Form: /__/ (x is horizontal, z diagonal up and right) 
     Xt = At*Sin(Phi_uc)
     Yt = Bt
     Zt = Ct + At*Dcos(Phi_uc)
  Case(6)
     !# TRICLINIC - C PARALLEL TO Z , BC ON XY PLANE
     Xt = At*pbcpar_uc(4)
     Yt = At*pbcpar_uc(5) + Bt*Dsin(Phi_a)
     Zt = At*Dcos(Phi_b) + Bt*Dcos(Phi_a) + Ct
  End Select
End Subroutine ABCtoXYZ
  
 
Subroutine XYZtoABC(Xt,Yt,Zt,At,Bt,Ct)
  Use Global
  Implicit none
 Real(8) :: At,Bt,Ct,Xt,Yt,Zt
 Select Case(PBCType)
  Case(1)
     !# ORTHORHOMBIC
     At = Xt
     Bt = Yt
     Ct = Zt
  Case(2)
     !# NANOTUBE
     At = Xt
     Bt = Yt
     Ct = Zt
  Case(3)                     !___
     !# MONOCLINIC B    Form: /__/ (x is horizontal, z diagonal up and right) 
     At = Xt/Dsin(Phi_uc)
     Bt = Yt
     Ct = Zt - Xt*DCOS(Phi_Uc)/DSIN(Phi_uc)
  Case(6)
     !# TRICLINIC - C PARALLEL TO Z , BC ON XY PLANE
     At = Xt/pbcpar_uc(4)
     Bt = Yt/DSIN(phi_a)-Xt/DSIN(Phi_a)*pbcpar_uc(5)/pbcpar_uc(4)
     Ct = Zt - Yt*DCOS(Phi_a)/DSIN(Phi_a)+Xt/Dsin(phi_a)*pbcpar_uc(5)/pbcpar_uc(4)
  End Select
End Subroutine XYZtoABC
