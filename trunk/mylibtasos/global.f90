! Adapting a global module representation of the simulation 
! in order to reduce the number of edited files for
! a more accurate versioning of the code. June 2003

Module Global
  Implicit none
  
  !#############################
  !## GLOBAL SIZE PARAMETERS  ##
  !#############################
  Integer, PARAMETER :: MAXADSORBENTSPECIES = 12 ! Maximum Adsorbent Species Allowed e.g. Si and O in Silicalite
  Integer, PARAMETER :: MAXADSORBATESPECIES = 3 ! Maximum Adsorbate Species Allowed e.g. CO2 and CH4 and H2
  Integer, PARAMETER :: MAXATOMSPECIES = 3      ! Maximum number of atoms each molecular species is allowed to have
  Integer, PARAMETER :: MAXINTERNALSPECIES = 2  ! Maximum number of intramolecular species allowed e.g. C and O in CO2
  Integer, PARAMETER :: MAXLJSIZE = 6           ! = MAXINTERNALSPECIES*MAXSPECIES max size for LJ adsorbate - adsorbate parameters array
  Integer, PARAMETER :: MAXATOMS = 20000         ! Maximum number of atoms for coordinate arrays

  !##################################
  !## GLOBAL CONSTANTS DECLARATION ##
  !##################################
  
  Real(8), PARAMETER :: PI = 3.14159265358979_8
  Real(8), PARAMETER :: TWOPI = 6.283185307179_8
  Real(8), PARAMETER :: SQRT2 = 1.414213562373095_8
  Real(8), PARAMETER :: Two_Thirds = 0.6666666666667_8
  Real(8), PARAMETER :: Rgas_J = 8.314410d0 ! Gas constant in J/mol/K
  Real(8), PARAMETER :: Rgas_bar = 83.1441d0 ! Gas constant in bar*cm^3/mol/K
  Real(8), PARAMETER :: activity_conversion = 0.6022045d0 ! to convert the activity in 1/A^3 using a density of mole/cm^3
  Real(8), PARAMETER :: fugacity_to_activity = 7.242902318d-3 ! to convert the activity in 1/A^3 with z = beta*f and f in bar
  Real(8), PARAMETER :: Coulomb_Constant = 167101.700d0

  !#################################################
  !## LOGICAL SETTINGS AND SELECTIONS             ##
  !#################################################
  INTEGER :: JOBREQUEST(2)
!!  Logical :: Dynamics = .True. 
!!  Logical :: MonteCarlo = .False.
!!  Logical :: ExploreGrid = .False.
!!  Logical :: NVT = .FALSE.     ! False means NVE
!!  Logical :: GCMC = .TRUE.     ! False means Henry's Law constant calculation  
!!  Logical :: PESurface = .FALSE.     ! True means the PES will be saved
  
  !##################################
  !## GLOBAL SIMULATION PARAMETERS ##
  !##################################
  Real(8) :: T                                  ! Temperature
  Real(8) :: beta                               ! Boltzman beta = kT
  
   
  !############################################
  !## ADSORBENT CRYSTAL STRUCTURE PARAMETERS ##
  !############################################
  LOGICAL :: PotentialInitialized = .FALSE. 
  LOGICAL :: ADSORBENTFIELD = .TRUE.
  LOGICAL :: Calculategrid(MAXADSORBATESPECIES) = .TRUE.,CalculateAgrid = .TRUE.
  
  Character(80) CrystalType
  Character(200) StructureFile,grid_file(MAXADSORBENTSPECIES)
  Character(2), DIMENSION(MAXADSORBENTSPECIES) :: CSPECIES_NAME
  Integer :: NCSPECIES                                  ! Number of species in the crystal structure
  Integer :: NCSPARTS(MAXADSORBENTSPECIES)              ! Number of particles for eachspecies
  Integer :: NCSPARTuc(MAXADSORBENTSPECIES)             ! Number of particles for eachspecies in the unit cell
  Real(8) :: Auc,Buc,Cuc,Phi_uc,Phi_a,Phi_b,Phi_c       ! Unit Cell Parameters
  Integer :: PBCType                                    ! Type of crystal
                                                        ! 1 CUBIC/Orhtorhombic
                                                        ! 2 Nanotube
                                                        ! 3 Monoclinic-B
                                                        ! 6 Triclinic 
  Real(8) :: PBCpar(15),PBCpar_crystal(15),PBCpar_uc(15)! Parameters for Periodic Boundaric Conditions
  Real(8) :: CrystalSpeciesDensity(MAXADSORBENTSPECIES) ! Density of crystal species for Long Range Corrections
  Real(8) :: NanotubeRadius                             ! Radius of a nanotube

  !###############################
  !## SIMULATION BOX PARAMETERS ##
  !###############################
  Integer :: nA,nB,nC                                   ! unit cells for the simulation box along the x(A),y(B),z(C) direction
  Integer :: nAt,nBt,nCt                                ! unit cells for grid generation box
  Real(8) :: Lx,Ly,Lz                                   ! Size of the simulation box
  Real(8) :: Volz,UnitCellVolume,rhoz(MAXADSORBENTSPECIES)! volume of the simulationbox and density of adsorbent species
  Real(8) :: Lxt,Lyt,Lzt                                ! Size of the minimum simulation box used for the grid generation
                                                        ! It depends on the maximum cut-off distance for 
                                                        !the adsorbate - adsorbent interaction
  Real(8), DIMENSION(1000,MAXADSORBENTSPECIES) :: Xzuc,Yzuc,Zzuc ! Coordinates of species' particles in a unit cell
  Real(8), ALLOCATABLE :: Xz(:,:),Yz(:,:),Zz(:,:)       ! Coordinates of species' particles in grid generating simulation box
  Real(8) :: maxcutoff                                  ! maximum cut-off distance
  Real(8) :: max_lj_cutoff
  !###############################################
  !## ADSORBATE MOLECULAR SYSTEM PARAMETERS     ##
  !###############################################
  Integer :: NSPECIES                                   ! Number of Species
  Integer :: NMOLECULES(MAXADSORBATESPECIES)            ! Number of molecules in each species
  Real(8) :: BULKPRESSURE(MAXADSORBATESPECIES)          ! Bulk pressure of species for Molecular Dynamics initialization
  Integer :: NATOMS(MAXADSORBATESPECIES)                ! Number of atoms in each molecular species
  Integer :: NUATOMS(MAXADSORBATESPECIES)               ! Number of unique atom identities in each species for interactions
  Integer :: TNATOMS,TNMOLECULES,TNUATOMS               ! Sum of the number of atoms and molecules of all species
  Integer :: TNATOMS_S,TNMOLECULES_S                    ! Sum of the number of atoms and molecules of all species
  Integer :: TNATOMS_L,TNMOLECULES_L                    ! Sum of the number of atoms and molecules of all species
  Real(8) :: MOLFRACTION(MAXADSORBATESPECIES)           ! Mol Fractions for Composition Input in Monte Carl
  Integer :: LinearSpecies = 0                          ! Number of species to be integrated with the Ciccotti integrator
  Integer :: LinearSpeciesID(MAXADSORBATESPECIES)       ! The Species ID of the species that are linear
  Integer :: SphericalSpecies = 0                       ! Number of spherical species 
  Integer :: SphericalSpeciesID(MAXADSORBATESPECIES)    ! The Species ID of the species that are spherical
  Integer :: ChargedSpecies = 0
  Integer :: ChargedSpeciesID(MAXADSORBATESPECIES)
  Logical :: ChargedCrystal = .FALSE.
  !###############################################
  !## COORDINATE MATRICES                       ##
  !###############################################
  ! Coordinate matrices are 2D and have dimensions (number_of_molecules , number_of_species
  Real(8), Dimension(MAXATOMS,MAXADSORBATESPECIES) :: Xm,Ym,Zm     ! Positions of particles Xm(molecule,atom,species)
  Real(8), Dimension(MAXATOMS,MAXADSORBATESPECIES) :: Vxm,Vym,Vzm  ! Velocities of particles
  Real(8), Dimension(MAXATOMS,MAXADSORBATESPECIES) :: Fxm,Fym,Fzm  ! Forces on particles
  !#### WHEN LINEAR SPECIES ARE PRESENT
  Real(8), Dimension(MAXATOMS,MAXADSORBATESPECIES) :: Xprime,Yprime,Zprime     ! uncorrected Positions of particles Xm(molecule,atom,species)
  Real(8), Dimension(MAXATOMS,MAXADSORBATESPECIES) :: Vxprime,Vyprime,Vzprime  ! uncorrected Velocities of particles
  Real(8), Dimension(MAXADSORBATESPECIES,MAXATOMS,3) :: WV                     ! angular velocity
  !### MOLECULE PROPERTIES
  Real(8) :: MolecularMass(MAXADSORBATESPECIES)                 ! Molecular mass of adsorbates 
  Real(8) :: Mumass(MAXADSORBATESPECIES)                        ! Geometric Mean Mass for Rigid Linear Molecules
  Real(8) :: Inertia(MAXADSORBATESPECIES)                       ! inertia for rigid linear molecules
  Real(8) :: AtomicMass(MAXATOMSPECIES,MAXADSORBATESPECIES) ! Atomic mass of internal atomic identities of adsorbates
  Real(8) :: iAtomicMass(MAXATOMSPECIES,MAXADSORBATESPECIES) ! inverse Atomic mass of internal atomic identities of adsorbates
  Integer :: MOLSTRCTID(MAXATOMSPECIES,MAXADSORBATESPECIES)     ! structure id array for polyatomic molecules
                                                                ! e.g.  1 -> 2 1 2 for CO2 where 1: C and 2 : O
  Real(8) :: MOLSTRCTCOORD(MAXATOMSPECIES,MAXADSORBATESPECIES)   ! structure coordinate array for polyatomic molecules
                                                                ! e.g.  1 -> -bond 0 bond for CO2 where 1: C and 2 : O
  Real(8) :: BOND(MAXADSORBATESPECIES),Bond4_2(MAXADSORBATESPECIES)
  !#### MOLECULAR ARRAY with dimension (TNMOLECULES,2)
  !#### First column is the id of the molecule and the second the id of the species
  !  Integer, Dimension(MAXATOMS,2) :: MOLARRAY
  !#### ARRAY OF ATOMIC TYPES for interaction matrices
  Integer, Dimension(MAXATOMSPECIES,MAXADSORBATESPECIES) :: ATOMTYPE
  Integer, Dimension(MAXINTERNALSPECIES,MAXADSORBATESPECIES)    :: UNIQUEATOMTYPE
  !#### ARRAY OF MOLECULAR TYPES for subroutine selection
  !#### 1 -> spherical particle
  !#### 2 -> polyatomic molecule
  !#### 3... -> future expansion
  Integer, Dimension(MAXADSORBATESPECIES) :: MOLECULETYPE
  !#### NAMES
  Character(10), Dimension(MAXADSORBATESPECIES) :: MOLECULENAME
  Character(2), Dimension(MAXATOMSPECIES,MAXADSORBATESPECIES) :: ATOMNAME
  Character(2), Dimension(MAXINTERNALSPECIES,MAXADSORBATESPECIES) :: UNIQUEATOMNAME

  !################################################
  !## INTERACTION PROPERTIES                     ##
  !################################################
  Logical :: INFINITEDILUTION = .FALSE.
  Integer :: PType_Adsorbate(MAXINTERNALSPECIES,MAXADSORBATESPECIES)
  Integer :: PType_Adsorbent(MAXADSORBENTSPECIES)
  Integer :: PTYPE_AC(MAXINTERNALSPECIES,MAXADSORBATESPECIES,MAXADSORBENTSPECIES)
  Integer :: aPT(maxljsize,maxljsize)
  Integer :: cPT(MAXADSORBENTSPECIES)
  Integer :: acPT(maxLJSIZE,MAXADSORBENTSPECIES)
  !### Assuming Lennard Jones Interactions
  !***************MOLECULE MOLECULE****************
  !#### MAIN LJ PARAMETERS
  !### LJ epsilon FOR ADSORBAT1ES
  Real(8) :: epsaa(MAXLJSIZE,MAXLJSIZE) 
  !### LJ sigma FOR ADSORBATES
  Real(8) :: sigaa(MAXLJSIZE,MAXLJSIZE)
  !### CUTOFF DISTANCE FOR LJ POTENTIAL
  Real(8) :: Rcaa(MAXLJSIZE,MAXLJSIZE)
  !### CUTOFF LJ ENERGY FOR DYNAMICS
  Real(8) :: ecutaa(MAXLJSIZE,MAXLJSIZE)
  !### DERIVATIVE PARAMETERS FOR QUICK COMPUTATION  !
  Real(8) :: e4aa(MAXLJSIZE,MAXLJSIZE)
  Real(8) :: e48aa(MAXLJSIZE,MAXLJSIZE)
  Real(8) :: sigaa2(MAXLJSIZE,MAXLJSIZE)
  Real(8) :: sigaa3(MAXLJSIZE,MAXLJSIZE)
  Real(8) :: invsigaa(MAXLJSIZE,MAXLJSIZE)
  Real(8) :: invsigaa2(MAXLJSIZE,MAXLJSIZE)
  Real(8) :: Rcaa2(MAXLJSIZE,MAXLJSIZE)
  Real(8) :: Rcaastar(MAXLJSIZE,MAXLJSIZE)
  Real(8) :: Rcaa2star(MAXLJSIZE,MAXLJSIZE)
  !###############################################
  !### FOR COULOMB INTERACTIONS                 ##
  !###############################################
  Real(8) :: aPC(MAXLJSIZE),aPC_cutoff(MAXLJSIZE)                      ! adsorbate partial charges
  Real(8) :: cPC(MAXADSORBENTSPECIES),cPC_cutoff(MAXADSORBENTSPECIES)        ! crystal partial charges
  Real(8) :: aPC_cutoff2(MAXLJSIZE),cPC_cutoff2(MAXADSORBENTSPECIES)
  Real(8) :: coulomb_cutoff,coulomb_cutoff2,icoulomb_cutoff
  !################################################
  !********* CRYSTAL SELF PARAMETERS **************
  Real(8) :: epscc(MAXADSORBENTSPECIES)
  Real(8) :: sigcc(MAXADSORBENTSPECIES)
  Real(8) :: Rccc(MAXADSORBENTSPECIES)
  Real(8) :: ecutcc(MAXADSORBENTSPECIES)
  !***************MOLECULE CRYSTAL*****************
  Real(8) :: epsac(MAXLJSIZE,MAXADSORBENTSPECIES)
  Real(8) :: sigac(MAXLJSIZE,MAXADSORBENTSPECIES)
  Real(8) :: Rcac(MAXLJSIZE,MAXADSORBENTSPECIES)
  Real(8) :: ecutac(MAXLJSIZE,MAXADSORBENTSPECIES)
  !-- Derivative parameters
  Real(8) :: e4ac(MAXLJSIZE,MAXADSORBENTSPECIES)
  Real(8) :: e48ac(MAXLJSIZE,MAXADSORBENTSPECIES)
  Real(8) :: e96ac(MAXLJSIZE,MAXADSORBENTSPECIES)
  Real(8) :: e384ac(MAXLJSIZE,MAXADSORBENTSPECIES)
  Real(8) :: sigac2(MAXLJSIZE,MAXADSORBENTSPECIES)
  Real(8) :: sigac3(MAXLJSIZE,MAXADSORBENTSPECIES)
  Real(8) :: invsigac(MAXLJSIZE,MAXADSORBENTSPECIES)
  Real(8) :: invsigac2(MAXLJSIZE,MAXADSORBENTSPECIES)
  Real(8) :: Rcac2(MAXLJSIZE,MAXADSORBENTSPECIES)
  Real(8) :: Rcacstar(MAXLJSIZE,MAXADSORBENTSPECIES)
  Real(8) :: Rcac2star(MAXLJSIZE,MAXADSORBENTSPECIES)
  !##################
  INTEGER :: CombinationRule ! 0 No combination rules, read potential
                             ! 1 Lorenz-Berthelot rules
                             ! 2 ... future expansion

  !###################################################
  !## INTERATOMIC CRYSTAL - MOLECULE POTENTIAL GRID ##
  !################################################### 
  Logical :: INTERPOLATION
  Real(8) :: grid_spacing = 0.2d0                     
  Integer :: ngx,ngy,ngz                                        ! Grid Size
  Real(8) :: hx,hy,hz                                           ! Grid Spacing
  !-- (i,j,k,id)
!  Integer ,Parameter :: GridX=135,GridY=135,GridZ=135,SPEC = 1
!  Real(8) ,Dimension(GridX,GridY,GridZ,SPEC) ::  FG,DxFG,DyFG,DzFG
!  Real(8) ,Dimension(GridX,GridY,GridZ,SPEC) ::  DxyFG,DxzFG,DyzFG
!  Real(8),Dimension(GridX,GridY,GridZ,SPEC) ::  DxyzFG
  Real(8), Allocatable,Dimension(:,:,:,:) ::  FG,DxFG,DyFG,DzFG
  Real(8), Allocatable,Dimension(:,:,:,:) ::  DxyFG,DxzFG,DyzFG
  Real(8), Allocatable,Dimension(:,:,:,:) ::  DxyzFG
  !-- i(j(k))
!  Real(8)  :: Xg(GridX),Yg(GridY),Zg(GridZ)
  Real(8), Allocatable,Dimension(:) :: Xg,Yg,Zg
  Real(8) :: minPE(MAXLJSIZE)
  Real(8) :: energy_cutoff,Vcv(MAXADSORBATEspecies) ! energy cutoff for channel volume calculation
  Integer :: Ncv(MAXADSORBATESPECIES)   ! Number of cubicles with a vertix less than energy_cuttoff
  Logical :: channelvolume = .false.
  Logical :: LongRangeCorrect = .false. ! logical parameter for the inclusion of longrange corrections in the grid/material potential generation
                                        ! for Grand Canonical Monte Carlo
  Real(8) :: LRangeCorrection(MAXLJSIZE,MAXADSORBENTSPECIES) ! Array for the precalculation of the Long Range Corrections
  Logical :: USESYMMETRY = .FALSE.
  Integer :: SYMMETRYTYPE ! 1 -> Cubic, 3D reflections in effect
  !###################################################
  !## MOLECULAR DYNAMICS PARAMETERS                 ##
  !###################################################
  !--- Timing information
  Real(8) :: dt,tsample,tmax                                   ! time step, sampling step, simulation time
  Real(8) :: dthalf,dtsqalf                                    ! dt = time step, dthalf = dt/2, dtsqalf = dt*dt/2
  Integer :: nsample,nmaxsteps
  !--- Equlibration informationm
  Integer :: MDsteps,MCsteps                                   ! Molecular Dynamic Steps, Monte Carlo steps
    
  !###################################################
  !## INTEGRATOR PARAMETERS                         ##
  !###################################################
  ! MTK Integrator for NVT
  Real(8), DIMENSION(11,3) :: GLOGS=0d0,VLOGS=0d0,XLOGS=0d0 ! Thermostat Forces, Velocities, and Positions
  Real(8), DIMENSION(11,3) :: QMASS=0d0
  Real(8), DIMENSION(5) :: WDT,WDT2,WDT4,WDT8
  Real(8) :: GNKT(3),GKT(3),QFREQUENCY
  Integer :: NNOS,NCMT,NYS,NFREE(3)
  Logical :: ZERONNOS = .FALSE.

  !#################################################
  !## MONTE CARLO PARAMETERS AND VARIABLES        ##
  !#################################################
  !--- Move information
  Integer :: CM_accmove(MAXADSORBATESPECIES),CM_attmove(MAXADSORBATESPECIES)           ! Accepted and Attempted moves
  Integer :: EM_accmove(MAXADSORBATESPECIES),EM_attmove(MAXADSORBATESPECIES)           ! Accepted and Attempted moves
  Integer :: SM_accmove(MAXADSORBATESPECIES),SM_attmove(MAXADSORBATESPECIES)           ! Accepted and Attempted moves
  Integer :: RM_accmove(MAXADSORBATESPECIES),RM_attmove(MAXADSORBATESPECIES)           ! Accepted and Attempted moves
  Real(8) :: Dr(MAXADSORBATESPECIES),Dtheta(MAXADSORBATESPECIES)
  !--- Equation of State Information and Chemical Potentials/Activities
  Real(8) :: Za(MAXADSORBATESPECIES)                                   ! activity of species
  Real(8) :: Pressure(MAXADSORBATESPECIES)                             ! pressure of species
  !--- Simulation Ranges
  Integer :: SimulationPoints,SimulationRanges
  Integer :: RangePoints(20)
  Real(8) :: PressureRange(20,2)       !! [StartPressure,EndPressure]
  Real(8) :: ActivityPoints(200,MAXADSORBATESPECIES)
  Real(8) :: PressurePoints(200,MAXADSORBATESPECIES)
  !---
  Logical :: DoDesorption
  !---- Equilibration and Production
  Integer :: GCMC_Equilibration,GCMC_Production
  !---- Moves and probabilities
  Integer :: MoveTypes
  !### movetype 1 -> Canonical Move
  !### movetype 2 -> Addition/Deletion
  !### movetype 3 -> Swapping
  !### movetype 4 -> Rotation
  Real(8) :: MoveProbability(10)   
  !---- Simulation Results
  !#### Each simulation will be separated into Sim_Bins parts - Averages will be based on those bins
  Integer, PARAMETER :: Sim_Bins = 10
  !#### data for each simulation point
  !#### all measurement variables are prefixed with m
  !---- Simulation Data
  !- N -> total number of adsorbed molecules
  !- Qst -> Isosteric Heat of adsorption (currently only of single species)  
  !- muex -> Excess chemical potential (currently only of single species)  
  !- prefix "m" -> measured
  Real(8), DIMENSION(1000,MAXADSORBATESPECIES) :: mNMOLECULES,mNMOLECULESvar
  Real(8), DIMENSION(1000,MAXADSORBATESPECIES) :: mQst,mQstvar
  Real(8), DIMENSION(1000,MAXADSORBATESPECIES) :: mPE,mPEvar
  Real(8), DIMENSION(1000,3) :: mPE_LJ,mPE_LJvar,mPE_Coul,mPE_Coulvar
  Real(8), DIMENSION(1000,MAXADSORBATESPECIES) :: muex
  !--- Henry's Law parameters
  Integer :: HenryMoves
  !--- Rotational Matrices
  Real(8), Dimension(3,3) :: aax,aay,aaz
  !#################################################
  !#### EQUATION OF STATE SPECIFICATIONS          ##
  !#################################################
  !-- EOS = 1 -> Ideal Gas
  !-- EOS = 2 -> Virial
  !-- EOS = 3... -> Other
  Integer :: EOS(MAXADSORBATESPECIES)
  Integer, PARAMETER :: NEOSPARAM(2) = (/0,3/)
  Real(8):: EOSPARAM(MAXADSORBATESPECIES,5)
  Integer :: CurrentEOS

  
  Real(8) ,PARAMETER,Dimension(2) :: CenterX=(/0.25,0.75/),CenterY=(/0.25,0.75/),CenterZ=(/0.25,0.75/)
  Real(8) ,PARAMETER :: CenterRadius = 7d0  
  Integer ,DIMENSION(1000) :: InitialPosition,CurrentPosition,Switchovers

  !#################################################
  !#### BRUTE FORCE COULOMB PARAMETERS            ##
  !#################################################
  Integer :: COULOMBLEVEL,NCELL
  Integer, ALLOCATABLE :: CELLLIST(:,:)
  Logical :: CalculateElectrostaticGrid = .TRUE.
  Character(200) :: electrostatic_grid_file

  Logical :: CHECKPOINT = .FALSE.

  !#################################################
  !### PARAMETERS AND VARIABLES FOR BIAS MODULES  ##
  !#################################################
  Integer :: BIASTYPE
  !----- ENERGY BIAS
  Integer :: EBatoms(MAXADSORBATESPECIES)
  Logical :: EBiasOn = .FALSE.
  Logical , DIMENSION(:,:,:,:),ALLOCATABLE :: EBiasLogical
  !----- CONFIGURATION BIAS
  Logical :: CBiasOn = .FALSE.

  Logical :: HardWallAdded = .FALSE.
  
  !###########################################################
  !### PARAMETERS AND VARIABLES FOR MC TRAJECTORY ANALYSIS ###
  !###########################################################
  Logical :: EnergyAnalysis = .FALSE.  
  Logical :: SaveMCtrajectory = .FALSE.
  ! Sampling Rate
  Integer :: trajsampling 
  

  !#### ENERGY HISTOGRAM
  Real(8) , ALLOCATABLE,DIMENSION(:) :: EnergyHistogram,EnergyHistogram1
  Real(8) :: minEHist,maxEHist,EHiststep
  Integer :: HistLength
  
  !#### Density Histogram
  Real(8), Allocatable,Dimension(:,:,:) ::  Density3D
  Real(8), Allocatable,Dimension(:,:,:) ::  Energy3DLattice,Energy3DTotal,Energy3DAdsorbate

  !#### Seggregation Analysis
  Real(8) :: Region_Loading(4)=0d0,TotalLoading=0d0

  !#################################################
  !#### RADIAL DISTRIBUTION FUNCTION VARIABLES    ##
  !#################################################
  ! G(r) for adsorbate-adsorbate and adsorbate -crystal
  !
  Real(8) , DIMENSION(1000,MAXADSORBATESPECIES,MAXADSORBENTSPECIES) :: CRYSTALRADIAL
  Integer :: NHISCRYSTAL,NCRYSTALGR
  Real(8) :: DELCRYSTALG
  Real(8),  DIMENSION(1000,MAXINTERNALSPECIES*MAXADSORBATESPECIES,MAXINTERNALSPECIES*MAXADSORBATESPECIES) :: ADSORBATERADIAL
  Integer :: NHISADSORBATE,NADSORBATEGR
  Real(8) :: DELADSORBATEG
End Module Global 
 
 
