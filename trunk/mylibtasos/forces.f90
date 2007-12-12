!#############################################################
!## Updated Subroutine for the Calculation of the Forces!!
!## July 2003
!#############################################################

Subroutine CForce_atom(AtomID,X,Y,Z,PE,Force) !### Forces Energy Due to the Adsorbent
  !--------------------------------------------------------
  !-- Sorbate - Sorbent interaction Forces and Potential --
  !-- Interpolated version                               --
  !--------------------------------------------------------
  Use Global
  !----------------------------------------------
  Implicit none
!  Integer ,INTENT(IN) :: AtomID
  Integer :: AtomID
  Real(8) ,INTENT(IN) :: X,Y,Z
  Real(8) ,DIMENSION(2) :: Xbox,Ybox,Zbox
  Real(8) :: PE,Force(3)
  Real(8) :: Xp,Yp,Zp
  Real(8) :: fintp,Forceint(3),Force2(3),PE2
  Real(8) :: F(8),DF(24),DDF(24),DDDF(8),XYZ(6),maxLJPE,minLJPE
  Integer :: I1,J1,K1,I2,J2,K2
  Integer :: FS(3)
  
  !-----------------------------------------------
  Call PBCuc(X,Y,Z,Xp,Yp,Zp)
  !------------------
  FS = 1
  !############ USING SYMMETRY OPERATORS FOR GRID REDUCTION ##
  !###
  If (USESYMMETRY) then
     SELECT Case (SYMMETRYTYPE)
     CASE(1)
        !### CUBIC LATTICE - 3 REFLECTIONS
        If (2*Xp - PBCpar_uc(1) .GT. 0.0001d0) then 
           Xp = PBCpar_uc(1) - Xp
           FS(1) = -1
        EndIf
        If (2*Yp - PBCpar_uc(2) .GT. 0.0001d0) then 
           Yp = PBCpar_uc(2) - Yp
           FS(2) = -1
        EndIf
        If (2*Zp - PBCpar_uc(3) .GT. 0.0001d0) then 
           Zp = PBCpar_uc(3) - Zp
           FS(3) = -1
        EndIf
     End SELECT
  EndIf

  !###
  !###########################################################
  I1=floor(Xp/hx+1)               ! Position in grid      
  J1=floor(Yp/hy+1)               ! Updated from floor(xx+1)
  K1=floor(Zp/hz+1)
  I2=I1+1
  J2=J1+1
  K2=K1+1
  !------------------
  If ((I1 .LT. 0).OR.(J1 .LT. 0).OR.(K1 .LT. 0)) then
     write(*,*) 'BLAKA..<0..',I1,J1,K1
     Write(*,*) 'ID ',AtomID
     write(*,*) ngx,ngy,ngz
     Write(*,*) X,Y,Z
     Write(*,*) Xp,Yp,Zp
    AtomID = -9999
!     STOP
     RETURN
  EndIf

 If ((I1 .GT. ngx+1) .OR. (J1 .GT. ngy+1) .OR. (K1 .GT. ngz+1)) then
    write(*,*) 'BLAKA..>ngx..',I1,J1,K1
    Write(*,*) 'ID ',AtomID
    write(*,*) ngx,ngy,ngz
    Write(*,*) X,Y,Z
    Write(*,*) Xp,Yp,Zp
    AtomID = -9999
!    STOP
    RETURN
 endif

  !------------------
  XYZ(1)=(xp-Xg(I1))/hx
  XYZ(2)=(yp-Yg(J1))/hy
  XYZ(3)=(zp-Zg(K1))/hz
  XYZ(4)=hx
  XYZ(5)=hy
  XYZ(6)=hz
!------------------
  !-- Function
  F(1)=FG(I1,J1,K1,AtomID)
  F(2)=FG(I2,J1,K1,AtomID)
  F(3)=FG(I1,J2,K1,AtomID)
  F(4)=FG(I1,J1,K2,AtomID)
  F(5)=FG(I2,J2,K1,AtomID)
  F(6)=FG(I2,J1,K2,AtomID)
  F(7)=FG(I1,J2,K2,AtomID)
  F(8)=FG(I2,J2,K2,AtomID)
!@  minLJPE=min(F(1),F(2),F(3),F(4),F(5),F(6),F(7),F(8))
!@  maxLJPE=max(F(1),F(2),F(3),F(4),F(5),F(6),F(7),F(8))
  !-- Dx
  DF(1)=DxFG(I1,J1,K1,AtomID)
  DF(2)=DxFG(I1,J2,K1,AtomID)
  DF(3)=DxFG(I1,J1,K2,AtomID)
  DF(4)=DxFG(I1,J2,K2,AtomID)
  DF(5)=DxFG(I2,J1,K1,AtomID)
  DF(6)=DxFG(I2,J2,K1,AtomID)
  DF(7)=DxFG(I2,J1,K2,AtomID)
  DF(8)=DxFG(I2,J2,K2,AtomID)
  !-- Dy
  DF(9)=DyFG(I1,J1,K1,AtomID)
  DF(10)=DyFG(I1,J2,K1,AtomID)
  DF(11)=DyFG(I1,J1,K2,AtomID)
  DF(12)=DyFG(I1,J2,K2,AtomID)
  DF(13)=DyFG(I2,J1,K1,AtomID)
  DF(14)=DyFG(I2,J2,K1,AtomID)
  DF(15)=DyFG(I2,J1,K2,AtomID)
  DF(16)=DyFG(I2,J2,K2,AtomID)
  !-- Dz
  DF(17)=DzFG(I1,J1,K1,AtomID)
  DF(18)=DzFG(I1,J2,K1,AtomID)
  DF(19)=DzFG(I1,J1,K2,AtomID)
  DF(20)=DzFG(I1,J2,K2,AtomID)
  DF(21)=DzFG(I2,J1,K1,AtomID)
  DF(22)=DzFG(I2,J2,K1,AtomID)
  DF(23)=DzFG(I2,J1,K2,AtomID)
  DF(24)=DzFG(I2,J2,K2,AtomID)
  !-- DxDy
  DDF(1)=DxyFG(I1,J1,K1,AtomID)
  DDF(2)=DxyFG(I1,J1,K2,AtomID)
  DDF(3)=DxyFG(I2,J1,K1,AtomID)
  DDF(4)=DxyFG(I2,J1,K2,AtomID)
  DDF(5)=DxyFG(I1,J2,K1,AtomID)
  DDF(6)=DxyFG(I1,J2,K2,AtomID)
  DDF(7)=DxyFG(I2,J2,K1,AtomID)
  DDF(8)=DxyFG(I2,J2,K2,AtomID)
  !-- DxDz
  DDF(9)=DxzFG(I1,J1,K1,AtomID)
  DDF(10)=DxzFG(I1,J1,K2,AtomID)
  DDF(11)=DxzFG(I2,J1,K1,AtomID)
  DDF(12)=DxzFG(I2,J1,K2,AtomID)
  DDF(13)=DxzFG(I1,J2,K1,AtomID)
  DDF(14)=DxzFG(I1,J2,K2,AtomID)
  DDF(15)=DxzFG(I2,J2,K1,AtomID)
  DDF(16)=DxzFG(I2,J2,K2,AtomID)
  !-- DyDz
  DDF(17)=DyzFG(I1,J1,K1,AtomID)
  DDF(18)=DyzFG(I1,J1,K2,AtomID)
  DDF(19)=DyzFG(I2,J1,K1,AtomID)
  DDF(20)=DyzFG(I2,J1,K2,AtomID)
  DDF(21)=DyzFG(I1,J2,K1,AtomID)
  DDF(22)=DyzFG(I1,J2,K2,AtomID)
  DDF(23)=DyzFG(I2,J2,K1,AtomID)
  DDF(24)=DyzFG(I2,J2,K2,AtomID)
  !-- DxDyDz
  DDDF(1)=DxyzFG(I1,J1,K1,AtomID)
  DDDF(2)=DxyzFG(I2,J1,K1,AtomID)
  DDDF(3)=DxyzFG(I1,J2,K1,AtomID)
  DDDF(4)=DxyzFG(I1,J1,K2,AtomID)
  DDDF(5)=DxyzFG(I2,J2,K1,AtomID)
  DDDF(6)=DxyzFG(I2,J1,K2,AtomID)
  DDDF(7)=DxyzFG(I1,J2,K2,AtomID)
  DDDF(8)=DxyzFG(I2,J2,K2,AtomID)
  !---------------------------------------------
  Call CHermInt3D(XYZ,F,DF,DDF,DDDF,Fintp,Forceint)
  !---------------------------------------------
  Force = -Forceint*dble(FS)
  if (fintp .LT. 1.5d0*minPE(AtomID)) then 
     PE = -fintp
  Else
     PE = fintp
  EndIf
  RETURN
END Subroutine CForce_atom

Subroutine CElectrostaticForce_atom(AtomID,X,Y,Z,Fintp,Force) !### Forces Energy Due to the Adsorbent
  !--------------------------------------------------------
  !-- Sorbate - Sorbent interaction Forces and Potential --
  !-- Interpolated version                               --
  !--------------------------------------------------------
  Use Global
  !----------------------------------------------
  Implicit none
  Integer ,INTENT(IN) :: AtomID
  Real(8) ,INTENT(IN) :: X,Y,Z
  Real(8) ,DIMENSION(2) :: Xbox,Ybox,Zbox
  Real(8) :: PE,Force(3)
  Real(8) :: Xp,Yp,Zp
  Real(8) :: fintp,Forceint(3),Force2(3),PE2,minElPE,maxElPE
  Real(8) :: F(8),DF(24),DDF(24),DDDF(8),XYZ(6)
  Integer :: I1,J1,K1,I2,J2,K2
  Integer :: FS(3)
  
  !-----------------------------------------------
  Call PBCuc(X,Y,Z,Xp,Yp,Zp)
  !------------------
  !------------------
  FS = 1 
  !############ USING SYMMETRY OPERATORS FOR GRID REDUCTION ##
  !###
  If (USESYMMETRY) then
     SELECT Case (SYMMETRYTYPE)
     CASE(1)
        !### CUBIC LATTICE - 3 REFLECTIONS
        If (2*Xp - PBCpar_uc(1) .GT. 0.0001d0) then 
           Xp = PBCpar_uc(1) - Xp
           FS(1) = -1
        EndIf
        If (2*Yp - PBCpar_uc(2) .GT. 0.0001d0) then 
           Yp = PBCpar_uc(2) - Yp
           FS(2) = -1
        EndIf
        If (2*Zp - PBCpar_uc(3) .GT. 0.0001d0) then 
           Zp = PBCpar_uc(3) - Zp
           FS(3) = -1
        EndIf
     End SELECT
  EndIf

  !###
  !###########################################################
  
  I1=floor(Xp/hx+1)               ! Position in grid      
  J1=floor(Yp/hy+1)               ! Updated from floor(xx+1)
  K1=floor(Zp/hz+1)
  I2=I1+1
  J2=J1+1
  K2=K1+1
  !------------------

 !------------------
  If ((I1 .LT. 0).OR.(J1 .LT. 0).OR.(K1 .LT. 0)) then
     write(*,*) 'Electrostatic - BLAKA..<0..',I1,J1,K1
     Write(*,*) 'ID ',AtomID
     write(*,*) ngx,ngy,ngz
     Write(*,*) X,Y,Z
     Write(*,*) Xp,Yp,Zp
     STOP
  EndIf

 If ((I1 .GT. ngx+1) .OR. (J1 .GT. ngy+1) .OR. (K1 .GT. ngz+1)) then
    write(*,*) 'Electrostatic - BLAKA..>ngx..',I1,J1,K1
    Write(*,*) 'ID ',AtomID
    write(*,*) ngx,ngy,ngz
    Write(*,*) X,Y,Z
    Write(*,*) Xp,Yp,Zp
    STOP
 endif

  !------------------
  XYZ(1)=(xp-Xg(I1))/hx
  XYZ(2)=(yp-Yg(J1))/hy
  XYZ(3)=(zp-Zg(K1))/hz
  XYZ(4)=hx
  XYZ(5)=hy
  XYZ(6)=hz
  !------------------
  !-- Function
  F(1)=FG(I1,J1,K1,TNUATOMS+1)
  F(2)=FG(I2,J1,K1,TNUATOMS+1)
  F(3)=FG(I1,J2,K1,TNUATOMS+1)
  F(4)=FG(I1,J1,K2,TNUATOMS+1)
  F(5)=FG(I2,J2,K1,TNUATOMS+1)
  F(6)=FG(I2,J1,K2,TNUATOMS+1)
  F(7)=FG(I1,J2,K2,TNUATOMS+1)
  F(8)=FG(I2,J2,K2,TNUATOMS+1)
  minElPE=min(F(1),F(2),F(3),F(4),F(5),F(6),F(7),F(8))
  maxElPE=max(F(1),F(2),F(3),F(4),F(5),F(6),F(7),F(8))
  !-- Dx
  DF(1)=DxFG(I1,J1,K1,TNUATOMS+1)
  DF(2)=DxFG(I1,J2,K1,TNUATOMS+1)
  DF(3)=DxFG(I1,J1,K2,TNUATOMS+1)
  DF(4)=DxFG(I1,J2,K2,TNUATOMS+1)
  DF(5)=DxFG(I2,J1,K1,TNUATOMS+1)
  DF(6)=DxFG(I2,J2,K1,TNUATOMS+1)
  DF(7)=DxFG(I2,J1,K2,TNUATOMS+1)
  DF(8)=DxFG(I2,J2,K2,TNUATOMS+1)
  !-- Dy
  DF(9)=DyFG(I1,J1,K1,TNUATOMS+1)
  DF(10)=DyFG(I1,J2,K1,TNUATOMS+1)
  DF(11)=DyFG(I1,J1,K2,TNUATOMS+1)
  DF(12)=DyFG(I1,J2,K2,TNUATOMS+1)
  DF(13)=DyFG(I2,J1,K1,TNUATOMS+1)
  DF(14)=DyFG(I2,J2,K1,TNUATOMS+1)
  DF(15)=DyFG(I2,J1,K2,TNUATOMS+1)
  DF(16)=DyFG(I2,J2,K2,TNUATOMS+1)
  !-- Dz
  DF(17)=DzFG(I1,J1,K1,TNUATOMS+1)
  DF(18)=DzFG(I1,J2,K1,TNUATOMS+1)
  DF(19)=DzFG(I1,J1,K2,TNUATOMS+1)
  DF(20)=DzFG(I1,J2,K2,TNUATOMS+1)
  DF(21)=DzFG(I2,J1,K1,TNUATOMS+1)
  DF(22)=DzFG(I2,J2,K1,TNUATOMS+1)
  DF(23)=DzFG(I2,J1,K2,TNUATOMS+1)
  DF(24)=DzFG(I2,J2,K2,TNUATOMS+1)
  !-- DxDy
  DDF(1)=DxyFG(I1,J1,K1,TNUATOMS+1)
  DDF(2)=DxyFG(I1,J1,K2,TNUATOMS+1)
  DDF(3)=DxyFG(I2,J1,K1,TNUATOMS+1)
  DDF(4)=DxyFG(I2,J1,K2,TNUATOMS+1)
  DDF(5)=DxyFG(I1,J2,K1,TNUATOMS+1)
  DDF(6)=DxyFG(I1,J2,K2,TNUATOMS+1)
  DDF(7)=DxyFG(I2,J2,K1,TNUATOMS+1)
  DDF(8)=DxyFG(I2,J2,K2,TNUATOMS+1)
  !-- DxDz
  DDF(9)=DxzFG(I1,J1,K1,TNUATOMS+1)
  DDF(10)=DxzFG(I1,J1,K2,TNUATOMS+1)
  DDF(11)=DxzFG(I2,J1,K1,TNUATOMS+1)
  DDF(12)=DxzFG(I2,J1,K2,TNUATOMS+1)
  DDF(13)=DxzFG(I1,J2,K1,TNUATOMS+1)
  DDF(14)=DxzFG(I1,J2,K2,TNUATOMS+1)
  DDF(15)=DxzFG(I2,J2,K1,TNUATOMS+1)
  DDF(16)=DxzFG(I2,J2,K2,TNUATOMS+1)
  !-- DyDz
  DDF(17)=DyzFG(I1,J1,K1,TNUATOMS+1)
  DDF(18)=DyzFG(I1,J1,K2,TNUATOMS+1)
  DDF(19)=DyzFG(I2,J1,K1,TNUATOMS+1)
  DDF(20)=DyzFG(I2,J1,K2,TNUATOMS+1)
  DDF(21)=DyzFG(I1,J2,K1,TNUATOMS+1)
  DDF(22)=DyzFG(I1,J2,K2,TNUATOMS+1)
  DDF(23)=DyzFG(I2,J2,K1,TNUATOMS+1)
  DDF(24)=DyzFG(I2,J2,K2,TNUATOMS+1)
  !-- DxDyDz
  DDDF(1)=DxyzFG(I1,J1,K1,TNUATOMS+1)
  DDDF(2)=DxyzFG(I2,J1,K1,TNUATOMS+1)
  DDDF(3)=DxyzFG(I1,J2,K1,TNUATOMS+1)
  DDDF(4)=DxyzFG(I1,J1,K2,TNUATOMS+1)
  DDDF(5)=DxyzFG(I2,J2,K1,TNUATOMS+1)
  DDDF(6)=DxyzFG(I2,J1,K2,TNUATOMS+1)
  DDDF(7)=DxyzFG(I1,J2,K2,TNUATOMS+1)
  DDDF(8)=DxyzFG(I2,J2,K2,TNUATOMS+1)

  !---------------------------------------------
  Call CHermInt3D(XYZ,F,DF,DDF,DDDF,Fintp,Forceint)
  !---------------------------------------------
  Force = -Forceint*aPC(AtomID)*dble(FS)
  Fintp = Fintp*aPC(AtomID)

  RETURN
END Subroutine CElectrostaticForce_atom
!=======================================================================================

Subroutine CForce_atom_full(AtomID,Xp,Yp,Zp,PE,Force)
!cccccccccccccccccccccccccccccccc
! Potential on a molecule due  cc
! to the zeolite structure     cc
!cccccccccccccccccccccccccccccccc
  Use Global
  !--------------------
  Implicit none  
  Integer ,INTENT(IN):: AtomID
  Real(8) ,INTENT(IN):: Xp,Yp,Zp
  Real(8) ,INTENT(OUT) :: PE,Force(3)
  Real(8) :: rx,ry,rz
  Real(8) :: rxt,ryt,rzt
  Real(8) :: R2,r2i,r6i,r2t
  Real(8) :: PEt,fr,Ft(3)
  Integer :: I,J   
  !--------------------
  PE = 0d0
  Force = 0d0
  !-------------------- Due to Zeolite
  Do J=1,NCSPECIES
     PEt = 0d0
     Ft = 0d0
     Do I=1,NCSPARTS(J)
        rx=Xp-Xz(I,J)
        ry=Yp-Yz(I,J)
        rz=Zp-Zz(I,J)       
        Call MinImage_Crystal(Rx,Ry,Rz,Rxt,Ryt,Rzt)
        R2=rxt**2 + ryt**2 + rzt**2
        If (R2 .LE. RCac2(AtomID,J)) then
           r2i=sigac2(AtomID,J)/R2
           r6i=r2i*r2i*r2i                    
           PEt=PEt+(r6i*r6i-r6i)-Ecutac(AtomID,J)
           fr=e48ac(AtomID,J)*(r6i*r6i-0.5d0*r6i)/R2 ! fr = Virij / R2 ### VIRIAL
           Ft(1)=Ft(1)+fr*Rxt
           Ft(2)=Ft(2)+fr*Ryt ! force on ith particle of jth species
           Ft(3)=Ft(3)+fr*Rzt
        EndIf
     EndDo
     !---------------------
     PE = PE + PEt*e4ac(AtomID,J)
     Force = Force + Ft
  EndDo
  RETURN
END Subroutine CForce_atom_full

!=======================================================================================

Subroutine TotalForce(PEaa,PEac,PEaacoul,PEaccoul,PEtotal)
!cccccccccccccccccccccccccc   
! TOTAL POTENTIAL ENERGY cc
!cccccccccccccccccccccccccc        
  Use Global
  !----------------------
  Implicit none
  Real(8) :: PEaa,PEac,PEtotal,PEaacoul,PEaccoul
  Real(8) :: PEaaself,PEaaother,PEaacoulself,PEaacoulother
  Real(8) :: Ft(3),PEt
  Real(8) :: R2,r6i,r2i,r3i,ri
  Real(8) :: Xi,Yi,Zi,Fxt,Fyt,Fzt
  Real(8) :: Fxprime(2),Fyprime(2),Fzprime(2)
  Real(8) :: Rx,Ry,Rz,fr,Rxt,Ryt,Rzt,mass(3)
  Real(8) :: coeff(3),coeff1,coeff3,Mcom
  Real(8), DIMENSION(MAXATOMSPECIES) :: X1,Y1,Z1,X2,Y2,Z2
  Integer :: ISPECIES,IMOL,JSPECIES,JMOL,IPART,JPART,ID1,ID2,INDEX,JNDEX,I

  !-----------------------
  PEaaself=0d0
  PEaaother=0d0
  PEaacoulself = 0d0
  PEaacoulother = 0d0
  PEaa = 0d0
  PEac=0d0
  PEtotal=0d0 
  PEt = 0d0
  PEaacoul = 0d0
  PEaccoul = 0d0
  !-----------------------      
  If (TNMOLECULES .EQ. 0) Return
  !----------------------------------------------------------
  !###### First calculate the forces due to the adsorbent
  !----------------------------------------------------------  
  If (ADSORBENTFIELD) then    
     SpeciesLoop: Do ISPECIES=1,NSPECIES
        MoleculesLoop: Do IMOL = 1,NMOLECULES(ISPECIES)
           ID1 = (IMOL-1)*NATOMS(ISPECIES) 
           AtomsLoop : Do IPART = 1,NATOMS(ISPECIES)
              Xi = Xm(ID1+IPART,ISPECIES)
              Yi = Ym(ID1+IPART,ISPECIES)
              Zi = Zm(ID1+IPART,ISPECIES)
              Call CForce_atom(ATOMTYPE(IPART,ISPECIES),Xi,Yi,Zi,PEt,Ft)
              PEac = PEac + PEt
              Fxm(ID1+IPART,ISPECIES) = Ft(1)
              Fym(ID1+IPART,ISPECIES) = Ft(2)
              Fzm(ID1+IPART,ISPECIES) = Ft(3)
           EndDo AtomsLoop
        EndDo MoleculesLoop
     EndDo SpeciesLoop
     If (ChargedCrystal) then 
        !##### UPDATE THIS SO THAT ONLY CHARGEDSPECIES ENTER THE LOOP!!!! THE SAME WITH THE POTENTIAL
        Do ISPECIES=1,NSPECIES
           Do IMOL = 1,NMOLECULES(ISPECIES)
              ID1 = (IMOL-1)*NATOMS(ISPECIES) 
              Do IPART = 1,NATOMS(ISPECIES)
                 Xi = Xm(ID1+IPART,ISPECIES)
                 Yi = Ym(ID1+IPART,ISPECIES)
                 Zi = Zm(ID1+IPART,ISPECIES)
                 Call CElectrostaticForce_atom(ATOMTYPE(IPART,ISPECIES),Xi,Yi,Zi,PEt,Ft)
                 PEaccoul = PEaccoul + PEt
                 Fxm(ID1+IPART,ISPECIES) = Fxm(ID1+IPART,ISPECIES) + Ft(1)
                 Fym(ID1+IPART,ISPECIES) = Fym(ID1+IPART,ISPECIES) + Ft(2)
                 Fzm(ID1+IPART,ISPECIES) = Fzm(ID1+IPART,ISPECIES) + Ft(3)
              EndDo
           EndDo
        EndDo
     EndIf
  Else
     Do ISPECIES=1,NSPECIES
        PEaccoul = 0d0
        Do IPART = 1,NATOMS(ISPECIES)*NMOLECULES(ISPECIES)
           Fxm(IPART,ISPECIES) = 0d0
           Fym(IPART,ISPECIES) = 0d0
           Fzm(IPART,ISPECIES) = 0d0
        EndDo
     EndDo    
  EndIf

  !-------------------------------------------------------------
  !###### Now calculate the forces between the adsorbates
  !-------------------------------------------------------------
  Infinite_Dilution_Check: If (.NOT. INFINITEDILUTION) then 
     OuterSpeciesLoop: Do ISPECIES = 1,NSPECIES
        
        OuterSelfLoop: Do IMOL=1,NMOLECULES(ISPECIES)-1
           INDEX = (IMOL-1)*NATOMS(ISPECIES) 
           Do I = 1,NATOMS(ISPECIES)
              X1(I) = Xm(INDEX+I,ISPECIES)
              Y1(I) = Ym(INDEX+I,ISPECIES)
              Z1(I) = Zm(INDEX+I,ISPECIES)
           EndDo
           InnerSelfLoop: Do JMOL=IMOL+1,NMOLECULES(ISPECIES)
              JNDEX = (JMOL-1)*NATOMS(ISPECIES) 
              Do I = 1,NATOMS(ISPECIES)
                 X2(I) = Xm(JNDEX+I,ISPECIES)
                 Y2(I) = Ym(JNDEX+I,ISPECIES)
                 Z2(I) = Zm(JNDEX+I,ISPECIES)
              EndDo
              OuterAtomLoop: Do IPART = 1,NATOMS(ISPECIES)
                 ID1 = AtomType(IPART,ISPECIES)
                 InnerAtomLoop: Do JPART = 1,NATOMS(ISPECIES)
                    ID2 = AtomType(JPART,ISPECIES)
                    If (aPT(ID1,ID2) .EQ. 0) cycle
                    Rx = X1(IPART)-X2(JPART)
                    Ry = Y1(IPART)-Y2(JPART)
                    Rz = Z1(IPART)-Z2(JPART)
                    Call MinImage(Rx,Ry,Rz,Rxt,Ryt,Rzt)
                    R2 = rxt*rxt + rzt*rzt + ryt*ryt 
                    Select Case(aPT(ID1,ID2)) 
                    Case(1) !### ONLY LJ
                       If (R2 .LE. RCaa2(ID1,ID2)) then  
                          r2i=sigaa2(ID1,ID2)/R2
                          r6i=r2i*r2i*r2i                    
                          PEaaself = PEaaself + e4aa(ID1,ID2)*((r6i*r6i-r6i)-Ecutaa(ID1,ID2))
                          fr=e48aa(ID1,ID2)*(r6i*r6i-0.5d0*r6i)/R2 ! fr = Virij / R2
                          Fxt = fr*Rxt
                          Fyt = fr*Ryt
                          Fzt = fr*Rzt
                          Fxm(INDEX+IPART,ISPECIES)=Fxm(INDEX+IPART,ISPECIES)+Fxt
                          Fym(INDEX+IPART,ISPECIES)=Fym(INDEX+IPART,ISPECIES)+Fyt ! force on Ith particle
                          Fzm(INDEX+IPART,ISPECIES)=Fzm(INDEX+IPART,ISPECIES)+Fzt
                          Fxm(JNDEX+JPART,ISPECIES)=Fxm(JNDEX+JPART,ISPECIES)-Fxt
                          Fym(JNDEX+JPART,ISPECIES)=Fym(JNDEX+JPART,ISPECIES)-Fyt ! force on Jth Particle (using newton's 3rd law)
                          Fzm(JNDEX+JPART,ISPECIES)=Fzm(JNDEX+JPART,ISPECIES)-Fzt
                       EndIf
                    Case(2) !### ONLY COULOMB
                       If (R2 .LE. coulomb_cutoff2) then 
                          ri = 1D0/Dsqrt(R2)
                          r3i = ri*ri*ri 
PEaacoulself = PEaacoulself + aPC(ID1)*aPC(ID2)*(ri - icoulomb_cutoff + icoulomb_cutoff**2*(Dsqrt(R2) - coulomb_cutoff))
                          fr = coulomb_constant*aPC(ID1)*aPc(ID2)*(r3i - icoulomb_cutoff**2*ri)
                          Fxt = fr*Rxt
                          Fyt = fr*Ryt
                          Fzt = fr*Rzt
                          Fxm(INDEX+IPART,ISPECIES)=Fxm(INDEX+IPART,ISPECIES)+Fxt
                          Fym(INDEX+IPART,ISPECIES)=Fym(INDEX+IPART,ISPECIES)+Fyt ! force on Ith particle
                          Fzm(INDEX+IPART,ISPECIES)=Fzm(INDEX+IPART,ISPECIES)+Fzt
                          Fxm(JNDEX+JPART,ISPECIES)=Fxm(JNDEX+JPART,ISPECIES)-Fxt
                          Fym(JNDEX+JPART,ISPECIES)=Fym(JNDEX+JPART,ISPECIES)-Fyt ! force on Jth Particle (using newton's 3rd law)
                          Fzm(JNDEX+JPART,ISPECIES)=Fzm(JNDEX+JPART,ISPECIES)-Fzt
                       EndIf
                    Case(3) !### BOTH LJ AND COULOMB
                       If (R2 .LE. RCaa2(ID1,ID2)) then  
                          r2i=sigaa2(ID1,ID2)/R2
                          r6i=r2i*r2i*r2i                    
                          PEaaself = PEaaself + e4aa(ID1,ID2)*((r6i*r6i-r6i)-Ecutaa(ID1,ID2))
                          fr=e48aa(ID1,ID2)*(r6i*r6i-0.5d0*r6i)/R2 ! fr = Virij / R2
                          Fxt = fr*Rxt
                          Fyt = fr*Ryt
                          Fzt = fr*Rzt
                          Fxm(INDEX+IPART,ISPECIES)=Fxm(INDEX+IPART,ISPECIES)+Fxt
                          Fym(INDEX+IPART,ISPECIES)=Fym(INDEX+IPART,ISPECIES)+Fyt ! force on Ith particle
                          Fzm(INDEX+IPART,ISPECIES)=Fzm(INDEX+IPART,ISPECIES)+Fzt
                          Fxm(JNDEX+JPART,ISPECIES)=Fxm(JNDEX+JPART,ISPECIES)-Fxt
                          Fym(JNDEX+JPART,ISPECIES)=Fym(JNDEX+JPART,ISPECIES)-Fyt ! force on Jth Particle (using newton's 3rd law)
                           Fzm(JNDEX+JPART,ISPECIES)=Fzm(JNDEX+JPART,ISPECIES)-Fzt
                       EndIf
                       If (R2 .LE. coulomb_cutoff2) then 
                          ri = 1D0/Dsqrt(R2)
                          r3i = ri*ri*ri
PEaacoulself = PEaacoulself + aPC(ID1)*aPC(ID2)*(ri - icoulomb_cutoff + icoulomb_cutoff**2*(Dsqrt(R2) - coulomb_cutoff))
                          fr = coulomb_constant*aPC(ID1)*aPc(ID2)*(r3i - icoulomb_cutoff**2*ri)
                          Fxt = fr*Rxt
                          Fyt = fr*Ryt
                          Fzt = fr*Rzt
                          Fxm(INDEX+IPART,ISPECIES)=Fxm(INDEX+IPART,ISPECIES)+Fxt
                          Fym(INDEX+IPART,ISPECIES)=Fym(INDEX+IPART,ISPECIES)+Fyt ! force on Ith particle
                          Fzm(INDEX+IPART,ISPECIES)=Fzm(INDEX+IPART,ISPECIES)+Fzt
                          Fxm(JNDEX+JPART,ISPECIES)=Fxm(JNDEX+JPART,ISPECIES)-Fxt
                          Fym(JNDEX+JPART,ISPECIES)=Fym(JNDEX+JPART,ISPECIES)-Fyt ! force on Jth Particle (using newton's 3rd law)
                          Fzm(JNDEX+JPART,ISPECIES)=Fzm(JNDEX+JPART,ISPECIES)-Fzt
                       EndIf
                    End Select
                 EndDo InnerAtomLoop
              EndDo OuterAtomLoop
           EndDo InnerSelfLoop
        EndDo OuterSelfLoop
        
        If (ISPECIES .LT. NSPECIES) then !### SO THAT NSPECIES DOESN'T TRIGGER A LOOP WITH ITSELF
           InnerSpeciesLoop: Do JSPECIES = ISPECIES+1,NSPECIES
              OuterMolLoop: Do IMOL = 1,NMOLECULES(ISPECIES)
                 INDEX = (IMOL - 1)*NATOMS(ISPECIES) 
                 Do I = 1,NATOMS(ISPECIES)
                    X1(I) = Xm(INDEX+I,ISPECIES)
                    Y1(I) = Ym(INDEX+I,ISPECIES)
                    Z1(I) = Zm(INDEX+I,ISPECIES)
                 EndDo
                 InnerMolLoop: Do JMOL = 1,NMOLECULES(JSPECIES)
                    JNDEX = (JMOL - 1)*NATOMS(JSPECIES)
                    Do I = 1,NATOMS(JSPECIES)
                       X2(I) = Xm(JNDEX+I,JSPECIES)
                       Y2(I) = Ym(JNDEX+I,JSPECIES)
                       Z2(I) = Zm(JNDEX+I,JSPECIES)
                    EndDo
                    OuterAtomLoop2: Do IPART = 1,NATOMS(ISPECIES) 
                       ID1 = AtomType(IPART,ISPECIES)
                       InnerAtomLoop2: Do JPART = 1,NATOMS(JSPECIES)
                          ID2 = AtomType(JPART,JSPECIES)                       
                          If (aPT(ID1,ID2) .EQ. 0) cycle
                          Rx = X1(IPART)-X2(JPART)
                          Ry = Y1(IPART)-Y2(JPART)
                          Rz = Z1(IPART)-Z2(JPART)
                          Call MinImage(Rx,Ry,Rz,Rxt,Ryt,Rzt)
                          R2 = rxt*rxt + rzt*rzt + ryt*ryt 
                          Select Case(aPT(ID1,ID2)) 
                          Case(1) !### ONLY LJ
                             If (R2 .LE. RCaa2(ID1,ID2)) then  
                                r2i=sigaa2(ID1,ID2)/R2
                                r6i=r2i*r2i*r2i                    
                                PEaaother = PEaaother + e4aa(ID1,ID2)*((r6i*r6i-r6i)-Ecutaa(ID1,ID2))
                                fr=e48aa(ID1,ID2)*(r6i*r6i-0.5d0*r6i)/R2 ! fr = Virij / R2
                                Fxt = fr*Rxt
                                Fyt = fr*Ryt
                                Fzt = fr*Rzt
                                Fxm(INDEX+IPART,ISPECIES)=Fxm(INDEX+IPART,ISPECIES)+Fxt
                                Fym(INDEX+IPART,ISPECIES)=Fym(INDEX+IPART,ISPECIES)+Fyt ! force on Ith particle
                                Fzm(INDEX+IPART,ISPECIES)=Fzm(INDEX+IPART,ISPECIES)+Fzt
                                Fxm(JNDEX+JPART,JSPECIES)=Fxm(JNDEX+JPART,JSPECIES)-Fxt
                                Fym(JNDEX+JPART,JSPECIES)=Fym(JNDEX+JPART,JSPECIES)-Fyt ! force on Jth Particle (using newton's 3rd law)
                                Fzm(JNDEX+JPART,JSPECIES)=Fzm(JNDEX+JPART,JSPECIES)-Fzt
                             EndIf
                          Case(2) !### ONLY COULOMB
                             If (R2 .LE. coulomb_cutoff2) then 
                                ri = 1D0/Dsqrt(R2)
                                r3i = ri*ri*ri 
PEaacoulother = PEaacoulother + aPC(ID1)*aPC(ID2)*(ri - icoulomb_cutoff + icoulomb_cutoff**2*(Dsqrt(R2) - coulomb_cutoff))
                                fr = coulomb_constant*aPC(ID1)*aPc(ID2)*(r3i- icoulomb_cutoff**2*ri)
                                Fxt = fr*Rxt
                                Fyt = fr*Ryt
                                Fzt = fr*Rzt
                                Fxm(INDEX+IPART,ISPECIES)=Fxm(INDEX+IPART,ISPECIES)+Fxt
                                Fym(INDEX+IPART,ISPECIES)=Fym(INDEX+IPART,ISPECIES)+Fyt ! force on Ith particle
                                Fzm(INDEX+IPART,ISPECIES)=Fzm(INDEX+IPART,ISPECIES)+Fzt
                                Fxm(JNDEX+JPART,JSPECIES)=Fxm(JNDEX+JPART,JSPECIES)-Fxt
                                Fym(JNDEX+JPART,JSPECIES)=Fym(JNDEX+JPART,JSPECIES)-Fyt ! force on Jth Particle (using newton's 3rd law)
                                Fzm(JNDEX+JPART,JSPECIES)=Fzm(JNDEX+JPART,JSPECIES)-Fzt
                             EndIf
                          Case(3) !### BOTH LJ AND COULOMB
                             If (R2 .LE. RCaa2(ID1,ID2)) then  
                                r2i=sigaa2(ID1,ID2)/R2
                                r6i=r2i*r2i*r2i                    
                                PEaaother = PEaaother + e4aa(ID1,ID2)*((r6i*r6i-r6i)-Ecutaa(ID1,ID2))
                                fr=e48aa(ID1,ID2)*(r6i*r6i-0.5d0*r6i)/R2 ! fr = Virij / R2
                                Fxt = fr*Rxt
                                Fyt = fr*Ryt
                                Fzt = fr*Rzt
                                Fxm(INDEX+IPART,ISPECIES)=Fxm(INDEX+IPART,ISPECIES)+Fxt
                                Fym(INDEX+IPART,ISPECIES)=Fym(INDEX+IPART,ISPECIES)+Fyt ! force on Ith particle
                                Fzm(INDEX+IPART,ISPECIES)=Fzm(INDEX+IPART,ISPECIES)+Fzt
                                Fxm(JNDEX+JPART,JSPECIES)=Fxm(JNDEX+JPART,JSPECIES)-Fxt
                                Fym(JNDEX+JPART,JSPECIES)=Fym(JNDEX+JPART,JSPECIES)-Fyt ! force on Jth Particle (using newton's 3rd law)
                                Fzm(JNDEX+JPART,JSPECIES)=Fzm(JNDEX+JPART,JSPECIES)-Fzt
                             EndIf
                             If (R2 .LE. coulomb_cutoff2) then 
                                ri = 1D0/Dsqrt(R2)
                                r3i = ri*ri*ri
PEaacoulother = PEaacoulother + aPC(ID1)*aPC(ID2)*(ri - icoulomb_cutoff + icoulomb_cutoff**2*(Dsqrt(R2) - coulomb_cutoff))
                                fr = coulomb_constant*aPC(ID1)*aPc(ID2)*(r3i- icoulomb_cutoff**2*ri)
                                Fxt = fr*Rxt
                                Fyt = fr*Ryt
                                Fzt = fr*Rzt
                                Fxm(INDEX+IPART,ISPECIES)=Fxm(INDEX+IPART,ISPECIES)+Fxt
                                Fym(INDEX+IPART,ISPECIES)=Fym(INDEX+IPART,ISPECIES)+Fyt ! force on Ith particle
                                Fzm(INDEX+IPART,ISPECIES)=Fzm(INDEX+IPART,ISPECIES)+Fzt
                                Fxm(JNDEX+JPART,JSPECIES)=Fxm(JNDEX+JPART,JSPECIES)-Fxt
                                Fym(JNDEX+JPART,JSPECIES)=Fym(JNDEX+JPART,JSPECIES)-Fyt ! force on Jth Particle (using newton's 3rd law)
                                Fzm(JNDEX+JPART,JSPECIES)=Fzm(JNDEX+JPART,JSPECIES)-Fzt
                             EndIf
                          End Select
                       EndDo InnerAtomLoop2
                    EndDo OuterAtomLoop2
                 EndDo InnerMolLoop
              EndDo OuterMolLoop
           EndDo InnerSpeciesLoop
        EndIf
     EndDo OuterSpeciesLoop
  EndIf Infinite_Dilution_Check

  PEaa = PEaaself + PEaaother 
  PEaacoul = coulomb_constant*(PEaacoulself + PEaacoulother)
  PETOTAL = PEaa + PEac + PEaacoul +  PEaccoul 
  !####
  !#### Now Distribute Forces for the linear molecules if they exist.
  If (LinearSpecies .GT. 0) then 
     TreatLinearSpecies: Do I = 1,LinearSpecies
        ISPECIES = LinearSpeciesID(I)
        mass = AtomicMass(:,ISPECIES)
        Mcom = molecularmass(ISPECIES)
        coeff(1)=(1d0-mass(2)/2d0/Mcom)
        coeff(3)=-mass(2)/2d0/Mcom 
        coeff1=mass(1)/Mcom
        coeff3=mass(3)/Mcom
        Do IMOL = 1,NMOLECULES(ISPECIES)
           INDEX = (IMOL-1)*NATOMS(ISPECIES)
           Fxprime(1) = coeff(1)*Fxm(INDEX+1,ISPECIES) + coeff1*Fxm(INDEX+2,ISPECIES) + coeff(3)*Fxm(INDEX+3,ISPECIES)
           Fxprime(2) = coeff(1)*Fxm(INDEX+3,ISPECIES) + coeff3*Fxm(INDEX+2,ISPECIES) + coeff(3)*Fxm(INDEX+1,ISPECIES)
           !---
           Fyprime(1) = coeff(1)*Fym(INDEX+1,ISPECIES) + coeff1*Fym(INDEX+2,ISPECIES) + coeff(3)*Fym(INDEX+3,ISPECIES)
           Fyprime(2) = coeff(1)*Fym(INDEX+3,ISPECIES) + coeff3*Fym(INDEX+2,ISPECIES) + coeff(3)*Fym(INDEX+1,ISPECIES)
           !---
           Fzprime(1) = coeff(1)*Fzm(INDEX+1,ISPECIES) + coeff1*Fzm(INDEX+2,ISPECIES) + coeff(3)*Fzm(INDEX+3,ISPECIES)
           Fzprime(2) = coeff(1)*Fzm(INDEX+3,ISPECIES) + coeff3*Fzm(INDEX+2,ISPECIES) + coeff(3)*Fzm(INDEX+1,ISPECIES)                    
           !### UPDATE BASIS PARTICLES' FORCES
           Fxm(INDEX+1,ISPECIES) = Fxprime(1)
           Fxm(INDEX+3,ISPECIES) = Fxprime(2)
           !---
           Fym(INDEX+1,ISPECIES) = Fyprime(1)
           Fym(INDEX+3,ISPECIES) = Fyprime(2)
           !---
           Fzm(INDEX+1,ISPECIES) = Fzprime(1)
           Fzm(INDEX+3,ISPECIES) = Fzprime(2)                 
        EndDo
     EndDo TreatLinearSpecies
  EndIf
  Return
END Subroutine TotalForce


