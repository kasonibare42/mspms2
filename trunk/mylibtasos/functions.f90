!##############################################
!## DECLARATION OF GLOBAL FUNCTIONS          ##
!##############################################

Module Functions
  Implicit none
  !##### POTENTIAL FUNCTIONS
  Real(8) , External :: AdsorbatePE
  Real(8) , External :: CPE_atom,CPE_atom2,CPE_molecule,CPE_atom_full,CPE_molecule_full
  Real(8) , External :: CElectrostaticPE_atom,CElectrostaticPE_molecule
  Real(8) , External :: CElectrostaticPE_atom_full,CElectrostaticPE_molecule_full
  !##### EOS FUNCTIONS
  Real(8) , External :: mu_excess_pure,density_p_pure,pvirial_pure,activity_rho_pure
  Real(8) , External :: mu_excess_mixture,density_p_mixture,pvirial_mixture,activity_rho_mixture
  !##### TIMING
  Real(8) , External :: Seconds
  !##### CHECKIFINCHANNEL
  logical , External :: checkifinchannel,energybiascheck
End Module Functions
