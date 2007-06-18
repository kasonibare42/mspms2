!##########################################################
!## Subroutine for the Calculation and Initialization of ##
!## the potential grids                                  ##
!##  - Calculate_Hermite_Grid                            ##
!##  - Calculate_Electrostatic_Hermite_Grid              ##       
!##  - Write_Hermite_Grid                                ##
!##  - Write_Electrostatic_Grid                          ##
!##  - Read_Grids                                        ##
!##  - Save_PES                                          ##
!##  - Calculate_Channel_Volume                          ##
!##########################################################
!## Author : Anastasios Skoulidas                        ##
!## Last Update : January, 16 2004                       ##
!##########################################################


Subroutine Write_Electrostatic_Grid()
  Use Global
  Use Functions
  Implicit none
  Integer:: I,IU,J,K
  open(Unit=3,file=electrostatic_grid_file,form='unformatted')
  Do K = 1,ngz+2
     Do J = 1,ngy+2
        Do I = 1,ngx+2
           Write(3) &
                FG(I,J,K,TNUATOMS+1),&
                DxFG(I,J,K,TNUATOMS+1),&
                DyFG(I,J,K,TNUATOMS+1),&
                DzFG(I,J,K,TNUATOMS+1),&
                DxyFG(I,J,K,TNUATOMS+1),&
                DxzFG(I,J,K,TNUATOMS+1),&
                DyzFG(I,J,K,TNUATOMS+1),&
                DxyzFG(I,J,K,TNUATOMS+1)
        EndDo
     EndDo
  EndDo
End Subroutine Write_Electrostatic_Grid

!########################################

Subroutine pass_grid_file_name(grid_file_name_cpp, idx)
  Use Global
  Implicit none
  Character(200) grid_file_name_cpp
  Integer:: idx

  grid_file(idx) = grid_file_name_cpp

!  print *, "grid file = ", grid_file(idx)

End Subroutine pass_grid_file_name

Subroutine Read_Grids(NSPECIES_yang)
  Use Global
  Use Functions
  Implicit none
  Integer:: NSPECIES_yang
  Integer:: I,IU,J,K,IATOM,ISPECIES,length

  minPE = 10d10
  Do K=1,ngz+2
     Zg(K)=(K-1)*hz     
  EndDo
  Do J=1,ngy+2   
     Yg(J)=(J-1)*hy
  EndDo
  Do I=1,ngx+2
     Xg(I)=(I-1)*hx
  EndDo
  
  Do ISPECIES = 1,NSPECIES_yang
     Open(Unit=3,file=grid_file(ISPECIES),form = 'formatted')
!     Do IATOM = 1,NUATOMS_yang(ISPECIES)
!     By default, one specie has only one unique atoms

!     One species has only one unique atoms
!     So, IU equals the id of the every speices
!  print *, "read grid file = ", grid_file(ISPECIES)
        IU = ISPECIES
        Do K = 1,ngz+2
           Do J = 1,ngy+2
              Do I = 1,ngx+2
                 Read(3,*) &
                      FG(I,J,K,IU),&
                      DxFG(I,J,K,IU),&
                      DyFG(I,J,K,IU),&
                      DzFG(I,J,K,IU),&
                      DxyFG(I,J,K,IU),&
                      DxzFG(I,J,K,IU),&
                      DyzFG(I,J,K,IU),&
                      DxyzFG(I,J,K,IU)
!    if (K==ngz+2)  then 
!      if (J==ngy+2)  then 
!        if (I==ngx+2) then 
!    write(*,*) FG(I,J,K,IU), DxFG(I,J,K,IU), DyFG(I,J,K,IU), DzFG(I,J,K,IU), DxyFG(I,J,K,IU),&
!               DxzFG(I,J,K,IU), DyzFG(I,J,K,IU), DxyzFG(I,J,K,IU)
!    endif
!    endif
!    endif

                 minPE(IU)=min(FG(I,J,K,IU),minPE(IU))
              EndDo
           EndDo
        EndDo
!     EndDo
     close(3)
  EndDo
  If (chargedcrystal) then
     open(Unit=3,file=electrostatic_grid_file,form='unformatted')
     Do K = 1,ngz+2
        Do J = 1,ngy+2
           Do I = 1,ngx+2
              Read(3) &
                   FG(I,J,K,TNUATOMS+1),&
                   DxFG(I,J,K,TNUATOMS+1),&
                   DyFG(I,J,K,TNUATOMS+1),&
                   DzFG(I,J,K,TNUATOMS+1),&
                   DxyFG(I,J,K,TNUATOMS+1),&
                   DxzFG(I,J,K,TNUATOMS+1),&
                   DyzFG(I,J,K,TNUATOMS+1),&
                   DxyzFG(I,J,K,TNUATOMS+1)
           EndDo
        EndDo
     EndDo
     close(3)
  EndIf
  Return 
End Subroutine Read_Grids

!########################################

Subroutine SavePES()
  Use Global 
  Use Functions
  Implicit none
  Integer:: I,IU,J,K,IATOM
  Open(unit=99,File='grid-data.dat',form='formatted')
  Write(99,*) ngx,ngy,ngz
  Close(99)
!@  open(unit=99,file='PES.dat',form='unformatted')
 !@ Do IU =1,TNUATOMS+1
  !@   write(99+IU) (((real(FG(I,J,K,IU)),I=1,ngx),J=1,ngy),K=1,ngz)
  !@   !@  close(99+I)
  !@EndDO
   write(100) (((real( FG(I,J,K,TNUATOMS+1) - 0.936*FG(I,J,K,TNUATOMS) ),I=1,ngx),J=1,ngy),K=1,ngz)
  !@EndDO
  Return
End Subroutine SavePES

!########################################

Subroutine CalculateChannelVolume(ATOMID)
  Use Global
  Use Functions
  Implicit None
  Real(8) Xp,Yp,Zp
  Integer I,J,K,aTOMID
  
  Ncv=0d0
  Do I=1,ngx
     Do J=1,ngy
        Do K=1,ngz
           Xp=(I-1)*hx+0.5d0*hx
           Yp=(J-1)*hy+0.5d0*hy
           Zp=(K-1)*hz+0.5d0*hz
           If (CheckIfinChannel(ATOMID,Xp,Yp,Zp)) Ncv(1)=Ncv(1)+1
        EndDo
     EndDo
  EndDo
  Vcv(1)=Ncv(1)*hx*hy*hz
  Open(Unit=55,File='channel_volume.dat',form='formatted')
  Write(55,*) 'Energy Cutoff = ',energy_cutoff
  Write(55,*) 'Total Channel Volume = ',Vcv(1)*na*nb*nc
  Write(55,*) 'Channel Volume per unit cell = ',Vcv(1)
  Write(55,*) 'Free Volume Calculated  = ',Vcv(1)/UnitCellVolume*100,'%'
  Close(55)
  Return
End Subroutine CalculateChannelVolume


Logical Function CheckIfinChannel(ATOMID,X,Y,Z)
  Use Global
  Implicit none
  
  Real*8 Xp,Yp,Zp,X,Y,Z
  Real*8 Xbox(2),Ybox(2),Zbox(2)
  Real*8 fintp,Dfintp(3)
  Real*8 F(8),DF(24),DDF(24),DDDF(8),XYZ(6)
  Integer I1,J1,K1,I2,J2,K2,FS(3),It,Jt,Kt,ATOMID
  
  Call PBCuc(X,Y,Z,XP,YP,ZP) ! PERIODIC BOUNDARY CONDITIONS TO THE UNIT CEL
  If (USESYMMETRY) then
     SELECT Case (SYMMETRYTYPE)
     CASE(1)
        !### CUBIC LATTICE - 3 REFLECTIONS
        If (2*Xp - PBCpar_uc(1) .GT. hx*0.5d0) Xp = PBCpar_uc(1) - Xp
        If (2*Yp - PBCpar_uc(2) .GT. hy*0.5d0) Yp = PBCpar_uc(2) - Yp
        If (2*Zp - PBCpar_uc(3) .GT. hz*0.5d0) Zp = PBCpar_uc(3) - Zp
     End SELECT
  EndIf
  !------------------
  I1=floor(Xp/hx+1)               ! Position in grid      
  J1=floor(Yp/hy+1)               ! Updated from floor(xx+1)
  K1=floor(Zp/hz+1)

  !------------------
  If ((I1 .LT. 0).OR.(J1 .LT. 0).OR.(K1 .LT. 0)) then
     write(*,*) 'BLAKA....',I1,J1,K1
     STOP
  EndIf
 If ((I1 .GT. ngx+1) .OR. (J1 .GT. ngy+1) .OR. (K1 .GT. ngz+1)) then
    write(*,*) 'BLAKA....',I1,J1,K1
    STOP
 endif

  I2=I1+1
  J2=J1+1
  K2=K1+1
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

  If (F(1) <= energy_cutoff .OR. F(2) <= energy_cutoff .OR. F(3) <= energy_cutoff .OR. F(4) <= energy_cutoff .OR.  &
       F(5) <= energy_cutoff .OR. F(6) <= energy_cutoff .OR. F(7) <= energy_cutoff .OR. F(8) <= energy_cutoff ) then
     CheckIfinChannel=.TRUE.
  Else
     CheckIfinChannel=.FALSE.
  EndIf
  Return
End Function CheckIfinChannel

!##################################################################################################
!==================================================================================================
