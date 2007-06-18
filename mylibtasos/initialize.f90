!################################################################################
!### AUXILIARY INITIALIZATION FUNCTIONS #########################################
!################################################################################

!Subroutine InitMoleculeArrays
!  Use Global
!  Implicit none
!  Logical :: Done
!  if (allocated(Xm)) deallocate(Xm)
!  if (allocated(Ym)) deallocate(Ym)
!  if (allocated(Zm)) deallocate(Zm)
!  if (allocated(Vxm)) deallocate(Vxm)
!  if (allocated(Vym)) deallocate(Vym)
!  if (allocated(Vzm)) deallocate(Vzm)
!  if (allocated(Fxm)) deallocate(Fxm)
!  if (allocated(Fym)) deallocate(Fym)
!  if (allocated(Fzm)) deallocate(Fzm)
!  if (allocated(atomtype)) deallocate(atomtype)
!  if (allocated(molarray)) deallocate(molarray)
!  allocate(Xm(maxatoms,nspecies))
!  allocate(Ym(maxatoms,nspecies))
!  allocate(Zm(maxatoms,nspecies))
!  allocate(Vxm(maxatoms,nspecies))
!  allocate(Vym(maxatoms,nspecies))
!  allocate(Vzm(maxatoms,nspecies))
!  allocate(Fxm(maxatoms,nspecies))
!  allocate(Fym(maxatoms,nspecies))
!  allocate(Fzm(maxatoms,nspecies))
!  allocate(AtomType(maxatoms,nspecies))
!  allocate(MolArray(maxatoms,2))
!  Done =  allocated(Xm).AND.allocated(Ym).AND.allocated(Zm).AND.allocated(Vxm).AND.allocated(Vym).AND.allocated(Vzm)&
!  .AND.allocated(Fxm).AND.allocated(Fym).AND.allocated(Fzm).AND.allocated(atomtype).AND.allocated(molarray)
!  If (.NOT.DONE) then
!     Write(*,*) 'Could not allocate all arrays!!! STOPPING'
!     STOP
!  EndIf
!  Return
!End Subroutine InitMoleculeArrays


Subroutine InitPotentialGrid(tnuatoms_cpp, auc_cpp, buc_cpp, cuc_cpp,&
                             na_cpp, nb_cpp, nc_cpp, nanotuberadius_cpp)
  Use Global
  Implicit none
  Integer:: tnuatoms_cpp, na_cpp, nb_cpp, nc_cpp
  Real(8):: auc_cpp, buc_cpp, cuc_cpp, nanotuberadius_cpp
  Logical :: Done

  !--- Only for nanotube structure

  PBCType = 2

  tnuatoms = tnuatoms_cpp

  Auc = auc_cpp
  Buc = buc_cpp
  Cuc = cuc_cpp
  NanotubeRadius = nanotuberadius_cpp

  nA = na_cpp
  nB = nb_cpp
  nC = nc_cpp

  Lx = nA*Auc
  Ly = nB*Buc
  Lz = nC*Cuc

  Volz = Lz*Pi*NanotubeRadius**2

  PBCpar(1) = Lx
  PBCpar(2) = Ly
  PBCpar(3) = Lz
  PBCpar(4) = NanotubeRadius
  PBCpar_uc(1) = Auc
  PBCpar_uc(2) = Buc
  PBCpar_uc(3) = Cuc
  PBCpar_uc(4) = NanotubeRadius

  UnitCellVolume = Cuc*Pi*NanotubeRadius**2

  !--- ---------------------------

  !--- compute grid size
  If (USESYMMETRY) then 
     Select CASE(SYMMETRYTYPE)
     CASE(1)
        ! CUBIC LATTICE
        ngx = ANINT(pbcpar_uc(1)/2d0/grid_spacing)  ! integer spacing
        hx = pbcpar_uc(1)/ngx/2d0                   ! realspacing
        ngy = ANINT(pbcpar_uc(2)/2d0/grid_spacing)  ! integer spacing
        hy = pbcpar_uc(2)/ngy/2d0                   ! realspacing
        ngz = ANINT(pbcpar_uc(3)/2d0/grid_spacing)  ! integer spacing
        hz = pbcpar_uc(3)/ngz/2d0                   ! realspacing 
     End Select
  Else
     ngx = ANINT(pbcpar_uc(1)/grid_spacing)  ! integer spacing
     hx = pbcpar_uc(1)/ngx                   ! realspacing
     ngy = ANINT(pbcpar_uc(2)/grid_spacing)  ! integer spacing
     hy = pbcpar_uc(2)/ngy                   ! realspacing
     ngz = ANINT(pbcpar_uc(3)/grid_spacing)  ! integer spacing
     hz = pbcpar_uc(3)/ngz                   ! realspacing 
  EndIf
  !--- check allocation
  if (allocated(FG)) deallocate(FG)
  if (allocated(DxFG)) deallocate(DxFG)
  if (allocated(DyFG)) deallocate(DyFG)
  if (allocated(DzFG)) deallocate(DzFG)
  if (allocated(DxyFG)) deallocate(DxyFG)
  if (allocated(DxzFG)) deallocate(DxzFG)
  if (allocated(DyzFG)) deallocate(DyzFG)
  if (allocated(DxyzFG)) deallocate(DxyzFG)
  if (allocated(Xg)) deallocate(Xg)
  if (allocated(Yg)) deallocate(Yg)
  if (allocated(Zg)) deallocate(Zg)
  if (allocated(EBiasLogical)) deallocate(EBiasLogical)
!  !--- allocate grid
   write(*,*) "tnuatoms = ", tnuatoms
  If (chargedcrystal) then 
     allocate(FG(ngx+2,ngy+2,ngz+2,tnuatoms+1))
     allocate(DxFG(ngx+2,ngy+2,ngz+2,tnuatoms+1))
     allocate(DyFG(ngx+2,ngy+2,ngz+2,tnuatoms+1))
     allocate(DzFG(ngx+2,ngy+2,ngz+2,tnuatoms+1))
     allocate(DxyFG(ngx+2,ngy+2,ngz+2,tnuatoms+1))
     allocate(DxzFG(ngx+2,ngy+2,ngz+2,tnuatoms+1))
     allocate(DyzFG(ngx+2,ngy+2,ngz+2,tnuatoms+1))
     allocate(DxyzFG(ngx+2,ngy+2,ngz+2,tnuatoms+1)) 
  else
     allocate(FG(ngx+2,ngy+2,ngz+2,tnuatoms))
     allocate(DxFG(ngx+2,ngy+2,ngz+2,tnuatoms))
     allocate(DyFG(ngx+2,ngy+2,ngz+2,tnuatoms))
     allocate(DzFG(ngx+2,ngy+2,ngz+2,tnuatoms))
     allocate(DxyFG(ngx+2,ngy+2,ngz+2,tnuatoms))
     allocate(DxzFG(ngx+2,ngy+2,ngz+2,tnuatoms))
     allocate(DyzFG(ngx+2,ngy+2,ngz+2,tnuatoms))
     allocate(DxyzFG(ngx+2,ngy+2,ngz+2,tnuatoms))
  endif
!  allocate(FG(ngx+2,ngy+2,ngz+2,1))
!  allocate(DxFG(ngx+2,ngy+2,ngz+2,1))
!  allocate(DyFG(ngx+2,ngy+2,ngz+2,1))
!  allocate(DzFG(ngx+2,ngy+2,ngz+2,1))
!  allocate(DxyFG(ngx+2,ngy+2,ngz+2,1))
!  allocate(DxzFG(ngx+2,ngy+2,ngz+2,1))
!  allocate(DyzFG(ngx+2,ngy+2,ngz+2,1))
!  allocate(DxyzFG(ngx+2,ngy+2,ngz+2,1))
  allocate(Xg(ngx+2))
  allocate(Yg(ngy+2))
  allocate(Zg(ngz+2))
  Done = ALLOCATED(FG).AND.ALLOCATED(DxFG).AND.ALLOCATED(DyFG).AND.ALLOCATED(DzFG).AND.ALLOCATED(DxyFG).AND.ALLOCATED(DxzFG)&
        .AND.ALLOCATED(DyzFG).AND.ALLOCATED(DxyzFG).AND.ALLOCATED(Xg).AND.ALLOCATED(Yg).AND.ALLOCATED(Zg)
  If (.NOT.DONE) then
     Write(*,*) 'Not all grid arrays allocated. STOPPING!!!'
     STOP
  EndIf
  Return
End Subroutine InitPotentialGrid

Subroutine Aux_Potential_Params
  Use Global
  Implicit none
  Integer :: I,J,K
  !### ADSORBATE - ADSORBATE
  e4aa = epsaa*4d0
  e48aa = epsaa*48d0
  sigaa2 = sigaa*sigaa
  sigaa3 = sigaa2*sigaa
  invsigaa = 1d0/sigaa
  invsigaa2 = invsigaa*invsigaa
  Rcaa2 = Rcaa*Rcaa
  Rcaastar = Rcaa*invsigaa
  Rcaa2star = Rcaastar*Rcaastar
  aPC_cutoff2 = aPC_cutoff*aPC_cutoff
  !### ADSORBENT - ADSORBATE
  e4ac = epsac*4d0
  e48ac = epsac*48d0
  e96ac = epsac*96d0
  e384ac = epsac*384d0
  sigac2 = sigac*sigac
  sigac3 = sigac2*sigac
  invsigac = 1d0/sigac
  invsigac2 = invsigac*invsigac
  Rcac2 = Rcac*Rcac
  Rcacstar = Rcac*invsigac
  Rcac2star = Rcacstar*Rcacstar
  cPC_cutoff2 = cPC_cutoff*cPC_cutoff

  coulomb_cutoff = 0d0
  !### maximum cutoff distance
  Do I = 1,TNUATOMS
     DO J = I,TNUATOMS
        max_LJ_cutoff = max(max_LJ_cutoff,Rcaa(I,J))
     EndDo
     coulomb_cutoff = max(coulomb_cutoff,aPC_cutoff(I))
  EndDo
  If (ADSORBENTFIELD) then
     Do I = 1,NCSPECIES
        DO J = 1,TNUATOMS
           max_LJ_cutoff = max(max_LJ_cutoff,RCac(J,I))
        EndDo
!@        maxcutoff = max(maxcutoff,cPC_cutoff(I))
!@        coulomb_cutoff = max(coulomb_cutoff,cPC_cutoff(I))
        !@@@ NOTE - THIS HAS TO BE ACTIVATED IF BRUTE FORCE CHARGES WITH CUTOFF POTENTIALS HAVE TO BE CALCULATED WITH THE STRUCTURE
     EndDo
  EndIf
  maxcutoff = max(max_LJ_cutoff,coulomb_cutoff)
  coulomb_cutoff2 = coulomb_cutoff*coulomb_cutoff
  icoulomb_cutoff = 1d0/coulomb_cutoff
  !### Cutoff Energies When Dynamics
  Select CASE(JOBREQUEST(1)) 
  CASE(1) ! DYNAMICS
     !If (Dynamics) then
     ecutaa=(1d0/Rcaa2star**6-1d0/Rcaa2star**3)
     ecutac=(1d0/Rcac2star**6-1d0/Rcac2star**3)
     LRangeCorrection = 0d0
     !  Else
  CASE(2) ! MONTE CARLO
     ecutac = 0d0
     ecutaa = 0d0
     If (LongRangeCorrect) then        
        DO I = 1,NSPECIES
           DO J = 1,NATOMS(I)
              Do K = 1,NCSPECIES
                 If (PTYPE_AC(J,I,K) .EQ. 1 .OR. PTYPE_AC(J,I,K) .EQ. 3) then                 
                    LRangeCorrection(ATOMTYPE(J,I),K)=((1d0/3d0)*1d0/RCacstar(ATOMTYPE(J,I),K)**9-1d0/RCacstar(ATOMTYPE(J,I),K)**3)
                 Else
                    LRangeCorrection(ATOMTYPE(J,I),K)=0
                 EndIf
              EndDo
           EndDo
        EndDo
     Else
        LRangeCorrection = 0d0
     Endif
     !  EndIf
  END Select
  
  PotentialInitialized = .TRUE.
  Return
End Subroutine Aux_Potential_Params

!#########################################################################################

Subroutine InitializeCoulombList()
 Use Global
 Implicit none
 Integer :: ILEVEL,ALLOCSIZE,I,J,K
 !###########################################################
 !### The condition of including a cell in the calculation
 !### is that icell^2 + jcell^2 + kcell^2 < maxcell
 !###########################################################
 !## Subroutine corrected on January 24th, 2003
 !###########################################################
 ALLOCSIZE = 0
 Do I = -CoulombLevel,CoulombLevel
    Do j = -CoulombLevel,CoulombLevel
       Do K = -CoulombLevel,CoulombLevel
          If ( I*I + J*J + K*K  .LE. CoulombLevel*CoulombLevel) then
             ALLOCSIZE = ALLOCSIZE + 1
          EndIf
       EndDo
    EndDo
 EndDO
 If (allocated(celllist)) deallocate(celllist)
 ALLOCATE(CELLLIST(ALLOCSIZE,3))

 NCELL = 0
 Do I = -CoulombLevel,CoulombLevel
    Do j = -CoulombLevel,CoulombLevel
       Do K = -CoulombLevel,CoulombLevel
          If ( I*I + J*J + K*K  .LE. CoulombLevel*CoulombLevel) then
             NCELL = NCELL + 1
             CELLLIST(NCELL,1) = I
             CELLLIST(NCELL,2) = J
             CELLLIST(NCELL,3) = K
          EndIf
       EndDo
    EndDo
 EndDo
 If (ALLOCSIZE .NE. NCELL) then
    Write(*,*) ALLOCSIZE,NCELL
    Write(*,*) 'STH IS WRONG'
 endif
End Subroutine InitializeCoulombList

