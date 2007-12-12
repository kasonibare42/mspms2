!#################################################
! Compute an interpolant to a given 3D function ##
! at a given position on a 3D grid using a 3D   ##
! cubic hermite interpolation scheme            ##
! Schultz, M. "Spline Analysis", Prentice-Hall: ##
! Englewood Cligs: NJ, 1973                     ##
!#################################################

Subroutine CHermInt3D(XYZ,F,DF,DDF,DDDF,Fintp,DFintp)
! Return interpolant and derivatives
  Use BasePolynomials
  !--------------------
  Implicit None
  Real(8) ,INTENT(IN) :: XYZ(6),F(8),DF(24),DDF(24),DDDF(8)
  Real(8) ,INTENT(OUT) :: Fintp,DFintp(3)
  Real(8) :: FV(4,3),FVprime(4,3)
  Real(8) :: BP0D(8),BP1D(24),BP2D(24),BP3D(8)
  Real(8) ::  f0D,f1D,f2D,f3D
  Integer :: I
  !---------------------   
  Call ComputeFV(FV,XYZ(1),XYZ(2),XYZ(3),XYZ(4),XYZ(5),XYZ(6))
  Call ComputeFVprime(FVprime,XYZ(1),XYZ(2),XYZ(3),XYZ(4),XYZ(5),XYZ(6))
  !--- Function Value
  Call NoD(BP0D,FV)
  Call OneD(BP1D,FV)
  Call TwoD(BP2D,FV)
  Call ThreeD(BP3D,FV)
  !---------------------
  f0D=0d0
  f1D=0d0
  f2D=0d0
  f3D=0d0
  !---------------------
  Do I=1,8
     f0D=f0D+F(I)*BP0D(I)
     f3D=f3D+DDDF(I)*BP3D(I)
  EndDo
  !---------------------
  Do I=1,24
     f1D=f1D+DF(I)*BP1D(I)
     f2D=f2D+DDF(I)*BP2D(I)
  EndDo
  !---------------------
  Fintp = f0D+f1D+f2D+f3D
  !---------------------
  !--- Dx Value
  Call NoDPrimeX(BP0D,FV,FVprime)
  Call OneDPrimeX(BP1D,FV,FVprime)
  Call TwoDPrimeX(BP2D,FV,FVprime)
  Call ThreeDPrimeX(BP3D,FV,FVprime)
  !---------------------
  f0D=0d0
  f1D=0d0
  f2D=0d0
  f3D=0d0
  !---------------------
  Do I=1,8
     f0D=f0D+F(I)*BP0D(I)
     f3D=f3D+DDDF(I)*BP3D(I)
  EndDo
  !---------------------
  Do I=1,24
     f1D=f1D+DF(I)*BP1D(I)
     f2D=f2D+DDF(I)*BP2D(I)
  EndDo
  !---------------------
  DFintp(1) = f0D+f1D+f2D+f3D
  !---------------------
  !--- Dy Value
  Call NoDPrimeY(BP0D,FV,FVprime)
  Call OneDPrimeY(BP1D,FV,FVprime)
  Call TwoDPrimeY(BP2D,FV,FVprime)
  Call ThreeDPrimeY(BP3D,FV,FVprime)
  !---------------------
  f0D=0d0
  f1D=0d0
  f2D=0d0
  f3D=0d0
  !---------------------
  Do I=1,8
     f0D=f0D+F(I)*BP0D(I)
     f3D=f3D+DDDF(I)*BP3D(I)
  EndDo
  !---------------------
  Do I=1,24
     f1D=f1D+DF(I)*BP1D(I)
     f2D=f2D+DDF(I)*BP2D(I)
  EndDo
  !---------------------
  DFintp(2) = f0D+f1D+f2D+f3D
  !---------------------
  !--- Dz Value
  Call NoDPrimeZ(BP0D,FV,FVprime)
  Call OneDPrimeZ(BP1D,FV,Fvprime)
  Call TwoDPrimeZ(BP2D,FV,FVprime)
  Call ThreeDPrimeZ(BP3D,FV,Fvprime)
  !---------------------
  f0D=0d0
  f1D=0d0
  f2D=0d0
  f3D=0d0
  !---------------------
  Do I=1,8
     f0D=f0D+F(I)*BP0D(I)
     f3D=f3D+DDDF(I)*BP3D(I)
  EndDo
  !---------------------
  Do I=1,24
     f1D=f1D+DF(I)*BP1D(I)
     f2D=f2D+DDF(I)*BP2D(I)
  EndDo
  !---------------------
  DFintp(3) = f0D+f1D+f2D+f3D
  !---------------------
  Return
End Subroutine CHermInt3D

!======================================================================================

Subroutine CHermInt3DNoD(XYZ,F,DF,DDF,DDDF,Fintp)
! Return interpolant without derivatives
  Use BasePolynomials
  !--------------------
  Implicit None
  Real(8) ,INTENT(IN) :: XYZ(6),F(8),DF(24),DDF(24),DDDF(8)
  Real(8) ,INTENT(OUT) :: Fintp
  Real(8) :: FV(4,3),FVprime(4,3)
  Real(8) :: BP0D(8),BP1D(24),BP2D(24),BP3D(8)
  Real(8) ::  f0D,f1D,f2D,f3D
  Integer :: I
  !---------------------      
  Call ComputeFV(FV,XYZ(1),XYZ(2),XYZ(3),XYZ(4),XYZ(5),XYZ(6))
  Call ComputeFVprime(FVprime,XYZ(1),XYZ(2),XYZ(3),XYZ(4),XYZ(5),XYZ(6))
  !---------------------      
  ! --- Function Value
  Call NoD(BP0D,FV)
  Call OneD(BP1D,FV)
  Call TwoD(BP2D,FV)
  Call ThreeD(BP3D,FV)
  !---------------------      
  f0D=0d0
  f1D=0d0
  f2D=0d0
  f3D=0d0
  !---------------------      
  Do I=1,8
     f0D=f0D+F(I)*BP0D(I)
     f3D=f3D+DDDF(I)*BP3D(I)
  EndDo
  !---------------------      
  Do I=1,24
     f1D=f1D+DF(I)*BP1D(I)
     f2D=f2D+DDF(I)*BP2D(I)
  EndDo
  !---------------------      
  Fintp = f0D+f1D+f2D+f3D
  Return
End Subroutine CHermInt3DNoD
