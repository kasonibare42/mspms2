Module CubicAuxFunctions
  Implicit none
Contains
  Real(8) Function f(x)
    Implicit none
    Real(8) ,INTENT(IN):: x      
    ! [ 3 multiplications 1 addition ]
    f=x*x*(2d0*x-3d0)      
    Return
  End Function f
  !==============================================
  Real(8) Function g(x)
    Implicit none
    Real(8) ,INTENT(IN):: x      
    ! [ 3 multiplications 2 additions ]
    g=x*x*(2d0*x-3d0)+1
    Return
  End Function g
  !==============================================
  Real(8) Function w(x)
    Implicit none
    Real(8) ,INTENT(IN):: x      
    ! [ 3 multiplications 2 additions ]
    w=x*(x-1)*(x-1)
    Return
  End Function w
  !==============================================
  Real(8) Function fprime(x)
    Implicit none
    Real(8) ,INTENT(IN):: x      
    ! [ 3 multiplications 1 addition ]
    fprime=6d0*x*(x-1d0)      
    Return
  End Function fprime
  !==============================================
  Real(8) Function gprime(x)
    Implicit none
    Real(8) ,INTENT(IN) :: x  
    ! [ 3 multiplications 1 addition ]    
    gprime=6d0*x*(x-1d0)      
    Return
  End Function gprime
  !==============================================
  Real(8) Function w1prime(x)
    Implicit none
    Real(8) ,INTENT(IN):: x            
    ! [ 2 multiplications 2 additions ]
    w1prime=x*(3d0*x-4d0)+1d0
    Return
  End Function w1prime
  !==============================================
  Real(8) Function w2prime(x)
    Implicit none
    Real(8) ,INTENT(IN):: x      
    ! [ 2 multiplications 1 addition ]
    w2prime=-x*(3d0*x-2d0)
    Return
  End Function w2prime
End Module CubicAuxFunctions

!=====================================================================================
!=====================================================================================

Module BasePolynomials
!-----------------------------------------------------------------------------
!---- BASE POLYNOMIALS e.g. BP
!     ***************************************************
!     * BP : array for the computed BP values           *
!     * FV : Function Values                            *
!     *      FV(4,3): |f(tx)      f(ty)      f(tz)     |*
!     *               |g(tx)      g(ty)      g(tz)     |*
!     *               |hx*w(tx)   hy*w(ty)   hz*w(tz)  |*
!     *               |hx*w(1-tx) hy*w(1-ty) hz*w(1-tz)|*
!     ***************************************************
Contains
  Subroutine ComputeFV(FV,tx,ty,tz,hx,hy,hz)
    Use CubicAuxFunctions
    !------------------------
    Implicit None    
    Real(8) ,INTENT(IN) :: tx,ty,tz,hx,hy,hz    
    Real(8) ,INTENT(OUT) :: FV(4,3)
    !------------------------
    ! [ 6 multiplications ]
    
    FV(1,1)=-f(tx)
    FV(1,2)=-f(ty)
    FV(1,3)=-f(tz)
    
    FV(2,1)=g(tx)
    FV(2,2)=g(ty)
    FV(2,3)=g(tz)
    
    FV(3,1)=hx*w(tx)
    FV(3,2)=hy*w(ty)
    FV(3,3)=hz*w(tz)
    
    FV(4,1)=-hx*w(1-tx)
    FV(4,2)=-hy*w(1-ty)
    FV(4,3)=-hz*w(1-tz)
    
    Return
  End Subroutine ComputeFV
  !============================================================
  Subroutine ComputeFVPrime(FVprime,tx,ty,tz,hx,hy,hz)
    Use CubicAuxFunctions
    !------------------------
    Implicit None    
    Real(8) ,INTENT(IN) :: tx,ty,tz,hx,hy,hz    
    Real(8) ,INTENT(OUT) :: FVprime(4,3)
    !------------------------
    ! [ 6 multiplications ]
    
    FVprime(1,1)=-fprime(tx)/hx
    FVprime(1,2)=-fprime(ty)/hy
    FVprime(1,3)=-fprime(tz)/hz
    
    FVprime(2,1)=gprime(tx)/hx
    FVprime(2,2)=gprime(ty)/hy
    FVprime(2,3)=gprime(tz)/hz
    
    FVprime(3,1)=w1prime(tx)
    FVprime(3,2)=w1prime(ty)
    FVprime(3,3)=w1prime(tz)
    
    FVprime(4,1)=-w2prime(tx)
    FVprime(4,2)=-w2prime(ty)
    FVprime(4,3)=-w2prime(tz)
    
    Return
  End Subroutine ComputeFVPrime
  !============================================================
  Subroutine NoD(BP,FV)
    Implicit none
    Real(8) ,INTENT(OUT) :: BP(8)
    Real(8) ,INTENT(IN) :: FV(4,3)
    ! [ 8x2 = 16 multiplications ]
    BP(1)= FV(2,1)*FV(2,2)*FV(2,3) ! i   j   k
    BP(2)= FV(1,1)*FV(2,2)*FV(2,3) ! i+1 j   k
    BP(3)= FV(2,1)*FV(1,2)*FV(2,3) ! i   j+1 k
    BP(4)= FV(2,1)*FV(2,2)*FV(1,3) ! i   j   k+1   
    BP(5)= FV(1,1)*FV(1,2)*FV(2,3) ! i+1 j+1 k
    BP(6)= FV(1,1)*FV(2,2)*FV(1,3) ! i+1 j   k+1
    BP(7)= FV(2,1)*FV(1,2)*FV(1,3) ! i   j+1 k+1
    BP(8)= FV(1,1)*FV(1,2)*FV(1,3) ! i+1 j+1 k+1
    
    Return 
  End Subroutine NoD
  !============================================================
  Subroutine OneD(BP,FV)
    Implicit none
    Real(8) ,INTENT(OUT) :: BP(24)
    Real(8) ,INTENT(IN) :: FV(4,3)
    !    [ 8x3x2 = 48 multiplications ]
    BP(1)= FV(3,1)*FV(2,2)*FV(2,3) ! i'   j   k
    BP(2)= FV(3,1)*FV(1,2)*FV(2,3) ! i'   j+1 k
    BP(3)= FV(3,1)*FV(2,2)*FV(1,3) ! i'   j   k+1
    BP(4)= FV(3,1)*FV(1,2)*FV(1,3) ! i'   j+1 k+1
    BP(5)= FV(4,1)*FV(2,2)*FV(2,3) ! i+1' j   k
    BP(6)= FV(4,1)*FV(1,2)*FV(2,3) ! i+1' j+1 k 
    BP(7)= FV(4,1)*FV(2,2)*FV(1,3) ! i+1' j   k+1
    BP(8)= FV(4,1)*FV(1,2)*FV(1,3) ! i+1' j+1 k+1
    
    BP(9)= FV(2,1)*FV(3,2)*FV(2,3) ! i   j'   k
    BP(10)= FV(2,1)*FV(4,2)*FV(2,3) ! i   j+1' k
    BP(11)= FV(2,1)*FV(3,2)*FV(1,3) ! i   j'   k+1
    BP(12)= FV(2,1)*FV(4,2)*FV(1,3) ! i   j+1' k+1
    BP(13)= FV(1,1)*FV(3,2)*FV(2,3) ! i+1 j'   k
    BP(14)= FV(1,1)*FV(4,2)*FV(2,3) ! i+1 j+1' k 
    BP(15)= FV(1,1)*FV(3,2)*FV(1,3) ! i+1 j'   k+1
    BP(16)= FV(1,1)*FV(4,2)*FV(1,3) ! i+1 j+1' k+1
    
    BP(17)= FV(2,1)*FV(2,2)*FV(3,3) ! i   j   k'
    BP(18)= FV(2,1)*FV(1,2)*FV(3,3) ! i   j+1 k'
    BP(19)= FV(2,1)*FV(2,2)*FV(4,3) ! i   j   k+1'
    BP(20)= FV(2,1)*FV(1,2)*FV(4,3) ! i   j+1 k+1'
    BP(21)= FV(1,1)*FV(2,2)*FV(3,3) ! i+1 j   k'
    BP(22)= FV(1,1)*FV(1,2)*FV(3,3) ! i+1 j+1 k' 
    BP(23)= FV(1,1)*FV(2,2)*FV(4,3) ! i+1 j   k+1'
    BP(24)= FV(1,1)*FV(1,2)*FV(4,3) ! i+1 j+1 k+1'
    
    Return 
  End Subroutine OneD
  !============================================================
  Subroutine TwoD(BP,FV)
    Implicit none
    Real(8) ,INTENT(OUT) :: BP(24)
    Real(8) ,INTENT(IN) :: FV(4,3)
    !  [ 8x3x2 = 48 multiplications ]
    BP(1)= FV(3,1)*FV(3,2)*FV(2,3) ! i'   j'   k
    BP(2)= FV(3,1)*FV(3,2)*FV(1,3) ! i'   j'   k+1
    BP(3)= FV(4,1)*FV(3,2)*FV(2,3) ! i+1' j'   k
    BP(4)= FV(4,1)*FV(3,2)*FV(1,3) ! i+1' j'   k+1
    BP(5)= FV(3,1)*FV(4,2)*FV(2,3) ! i'   j+1' k
    BP(6)= FV(3,1)*FV(4,2)*FV(1,3) ! i'   j+1' k+1
    BP(7)= FV(4,1)*FV(4,2)*FV(2,3) ! i+1' j+1' k
    BP(8)= FV(4,1)*FV(4,2)*FV(1,3) ! i+1' j+1' k+1
    
    BP(9)= FV(3,1)*FV(2,2)*FV(3,3) ! i'   j   k'
    BP(10)= FV(3,1)*FV(2,2)*FV(4,3) ! i'   j   k+1'
    BP(11)= FV(4,1)*FV(2,2)*FV(3,3) ! i+1' j   k'
    BP(12)= FV(4,1)*FV(2,2)*FV(4,3) ! i+1' j   k+1'
    BP(13)= FV(3,1)*FV(1,2)*FV(3,3) ! i'   j+1 k'
    BP(14)= FV(3,1)*FV(1,2)*FV(4,3) ! i'   j+1 k+1'
    BP(15)= FV(4,1)*FV(1,2)*FV(3,3) ! i+1' j+1 k'
    BP(16)= FV(4,1)*FV(1,2)*FV(4,3) ! i+1' j+1 k+1'
    
    BP(17)= FV(2,1)*FV(3,2)*FV(3,3) ! i   j'   k'
    BP(18)= FV(2,1)*FV(3,2)*FV(4,3) ! i   j'   k+1'
    BP(19)= FV(1,1)*FV(3,2)*FV(3,3) ! i+1 j'   k'
    BP(20)= FV(1,1)*FV(3,2)*FV(4,3) ! i+1 j'   k+1'
    BP(21)= FV(2,1)*FV(4,2)*FV(3,3) ! i   j+1' k'
    BP(22)= FV(2,1)*FV(4,2)*FV(4,3) ! i   j+1' k+1'
    BP(23)= FV(1,1)*FV(4,2)*FV(3,3) ! i+1 j+1' k'
    BP(24)= FV(1,1)*FV(4,2)*FV(4,3) ! i+1 j+1' k+1'
    
    Return 
  End Subroutine TwoD
  !============================================================
  Subroutine ThreeD(BP,FV)
    Implicit none
    Real(8) ,INTENT(OUT):: BP(8)
    Real(8) ,INTENT(IN) :: FV(4,3)
    !   [ 8x2 = 16 multiplications ]
    BP(1)= FV(3,1)*FV(3,2)*FV(3,3) ! i'   j'   k'
    BP(2)= FV(4,1)*FV(3,2)*FV(3,3) ! i+1' j'   k'
    BP(3)= FV(3,1)*FV(4,2)*FV(3,3) ! i'   j+1' k'
    BP(5)= FV(3,1)*FV(3,2)*FV(4,3) ! i'   j'   k+1'
    BP(4)= FV(4,1)*FV(4,2)*FV(3,3) ! i+1' j+1' k'
    BP(6)= FV(4,1)*FV(3,2)*FV(4,3) ! i+1' j'   k+1'
    BP(7)= FV(3,1)*FV(4,2)*FV(4,3) ! i'   j+1' k+1'
    BP(8)= FV(4,1)*FV(4,2)*FV(4,3) ! i+1' j+1' k+1'
    Return 
  End Subroutine ThreeD
  !============================================================
  !============================================================
  !============================================================
  Subroutine NoDPrimeX(BP,FV,FVprime)
    Implicit none
    Real(8) ,INTENT(OUT) :: BP(8)
    Real(8) ,INTENT(IN) :: FV(4,3)
    Real(8) ,INTENT(IN) :: FVprime(4,3)
    !   [ 8x2 = 16 multiplications ]

    BP(1)= FVprime(2,1)*FV(2,2)*FV(2,3) ! i   j   k
    BP(2)= FVprime(1,1)*FV(2,2)*FV(2,3) ! i+1 j   k
    BP(3)= FVprime(2,1)*FV(1,2)*FV(2,3) ! i   j+1 k
    BP(4)= FVprime(2,1)*FV(2,2)*FV(1,3) ! i   j   k+1   
    BP(5)= FVprime(1,1)*FV(1,2)*FV(2,3) ! i+1 j+1 k
    BP(6)= FVprime(1,1)*FV(2,2)*FV(1,3) ! i+1 j   k+1
    BP(7)= FVprime(2,1)*FV(1,2)*FV(1,3) ! i   j+1 k+1
    BP(8)= FVprime(1,1)*FV(1,2)*FV(1,3) ! i+1 j+1 k+1
    
    Return 
  End Subroutine NoDPrimeX
  !============================================================
  Subroutine OneDPrimeX(BP,FV,Fvprime)
    Implicit none
    Real(8) ,INTENT(OUT) :: BP(24)
    Real(8) ,INTENT(IN) :: FV(4,3)
    Real(8) ,INTENT(IN) :: FVprime(4,3)
    !     [ 8x3x2 = 48 multiplications ]
    
    BP(1)= FVprime(3,1)*FV(2,2)*FV(2,3) ! i'   j   k
    BP(2)= FVprime(3,1)*FV(1,2)*FV(2,3) ! i'   j+1 k
    BP(3)= FVprime(3,1)*FV(2,2)*FV(1,3) ! i'   j   k+1
    BP(4)= FVprime(3,1)*FV(1,2)*FV(1,3) ! i'   j+1 k+1
    BP(5)= FVprime(4,1)*FV(2,2)*FV(2,3) ! i+1' j   k
    BP(6)= FVprime(4,1)*FV(1,2)*FV(2,3) ! i+1' j+1 k 
    BP(7)= FVprime(4,1)*FV(2,2)*FV(1,3) ! i+1' j   k+1
    BP(8)= FVprime(4,1)*FV(1,2)*FV(1,3) ! i+1' j+1 k+1
    
    BP(9)= FVprime(2,1)*FV(3,2)*FV(2,3) ! i   j'   k
    BP(10)= FVprime(2,1)*FV(4,2)*FV(2,3) ! i   j+1' k
    BP(11)= FVprime(2,1)*FV(3,2)*FV(1,3) ! i   j'   k+1
    BP(12)= FVprime(2,1)*FV(4,2)*FV(1,3) ! i   j+1' k+1
    BP(13)= FVprime(1,1)*FV(3,2)*FV(2,3) ! i+1 j'   k
    BP(14)= FVprime(1,1)*FV(4,2)*FV(2,3) ! i+1 j+1' k 
    BP(15)= FVprime(1,1)*FV(3,2)*FV(1,3) ! i+1 j'   k+1
    BP(16)= FVprime(1,1)*FV(4,2)*FV(1,3) ! i+1 j+1' k+1
    
    BP(17)= FVprime(2,1)*FV(2,2)*FV(3,3) ! i   j   k'
    BP(18)= FVprime(2,1)*FV(1,2)*FV(3,3) ! i   j+1 k'
    BP(19)= FVprime(2,1)*FV(2,2)*FV(4,3) ! i   j   k+1'
    BP(20)= FVprime(2,1)*FV(1,2)*FV(4,3) ! i   j+1 k+1'
    BP(21)= FVprime(1,1)*FV(2,2)*FV(3,3) ! i+1 j   k'
    BP(22)= FVprime(1,1)*FV(1,2)*FV(3,3) ! i+1 j+1 k' 
    BP(23)= FVprime(1,1)*FV(2,2)*FV(4,3) ! i+1 j   k+1'
    BP(24)= FVprime(1,1)*FV(1,2)*FV(4,3) ! i+1 j+1 k+1'
    
    Return 
  End Subroutine OneDPrimeX
  !============================================================
  Subroutine TwoDPrimeX(BP,FV,FVprime)
    Implicit none
    Real(8) ,INTENT(OUT) :: BP(24)
    Real(8) ,INTENT(IN) :: FV(4,3)
    Real(8) ,INTENT(IN) :: FVprime(4,3)
    !    [ 8x3x2 = 48 multiplications ]
    
    BP(1)= FVprime(3,1)*FV(3,2)*FV(2,3) ! i'   j'   k
    BP(2)= FVprime(3,1)*FV(3,2)*FV(1,3) ! i'   j'   k+1
    BP(3)= FVprime(4,1)*FV(3,2)*FV(2,3) ! i+1' j'   k
    BP(4)= FVprime(4,1)*FV(3,2)*FV(1,3) ! i+1' j'   k+1
    BP(5)= FVprime(3,1)*FV(4,2)*FV(2,3) ! i'   j+1' k
    BP(6)= FVprime(3,1)*FV(4,2)*FV(1,3) ! i'   j+1' k+1
    BP(7)= FVprime(4,1)*FV(4,2)*FV(2,3) ! i+1' j+1' k
    BP(8)= FVprime(4,1)*FV(4,2)*FV(1,3) ! i+1' j+1' k+1
    
    BP(9)= FVprime(3,1)*FV(2,2)*FV(3,3) ! i'   j   k'
    BP(10)= FVprime(3,1)*FV(2,2)*FV(4,3) ! i'   j   k+1'
    BP(11)= FVprime(4,1)*FV(2,2)*FV(3,3) ! i+1' j   k'
    BP(12)= FVprime(4,1)*FV(2,2)*FV(4,3) ! i+1' j   k+1'
    BP(13)= FVprime(3,1)*FV(1,2)*FV(3,3) ! i'   j+1 k'
    BP(14)= FVprime(3,1)*FV(1,2)*FV(4,3) ! i'   j+1 k+1'
    BP(15)= FVprime(4,1)*FV(1,2)*FV(3,3) ! i+1' j+1 k'
    BP(16)= FVprime(4,1)*FV(1,2)*FV(4,3) ! i+1' j+1 k+1'
    
    BP(17)= FVprime(2,1)*FV(3,2)*FV(3,3) ! i   j'   k'
    BP(18)= FVprime(2,1)*FV(3,2)*FV(4,3) ! i   j'   k+1'
    BP(19)= FVprime(1,1)*FV(3,2)*FV(3,3) ! i+1 j'   k'
    BP(20)= FVprime(1,1)*FV(3,2)*FV(4,3) ! i+1 j'   k+1'
    BP(21)= FVprime(2,1)*FV(4,2)*FV(3,3) ! i   j+1' k'
    BP(22)= FVprime(2,1)*FV(4,2)*FV(4,3) ! i   j+1' k+1'
    BP(23)= FVprime(1,1)*FV(4,2)*FV(3,3) ! i+1 j+1' k'
    BP(24)= FVprime(1,1)*FV(4,2)*FV(4,3) ! i+1 j+1' k+1'
    
    Return 
  End Subroutine TwoDPrimeX
  !============================================================
  Subroutine ThreeDPrimeX(BP,FV,FVprime)
    Implicit none
    Real(8) ,INTENT(OUT) :: BP(8)
    Real(8) ,INTENT(IN) :: FV(4,3)
    Real(8) ,INTENT(IN) :: FVprime(4,3)     
    !     [ 8x2 = 16 multiplications ]
    BP(1)= FVprime(3,1)*FV(3,2)*FV(3,3) ! i'   j'   k'
    BP(2)= FVprime(4,1)*FV(3,2)*FV(3,3) ! i+1' j'   k'
    BP(3)= FVprime(3,1)*FV(4,2)*FV(3,3) ! i'   j+1' k'
    BP(5)= FVprime(3,1)*FV(3,2)*FV(4,3) ! i'   j'   k+1'
    BP(4)= FVprime(4,1)*FV(4,2)*FV(3,3) ! i+1' j+1' k'
    BP(6)= FVprime(4,1)*FV(3,2)*FV(4,3) ! i+1' j'   k+1'
    BP(7)= FVprime(3,1)*FV(4,2)*FV(4,3) ! i'   j+1' k+1'
    BP(8)= FVprime(4,1)*FV(4,2)*FV(4,3) ! i+1' j+1' k+1'
    
    Return 
  End Subroutine ThreeDPrimeX
  !============================================================
  !============================================================
  Subroutine NoDPrimeY(BP,FV,FVPrime)
    Implicit none
    Real(8) ,INTENT(OUT) :: BP(8)
    Real(8) ,INTENT(IN) :: FV(4,3)
    Real(8) ,INTENT(IN) :: FVprime(4,3)
    !   [ 8x2 = 16 multiplications ]
    
    BP(1)= FV(2,1)*FVprime(2,2)*FV(2,3) ! i   j   k
    BP(2)= FV(1,1)*FVprime(2,2)*FV(2,3) ! i+1 j   k
    BP(3)= FV(2,1)*FVprime(1,2)*FV(2,3) ! i   j+1 k
    BP(4)= FV(2,1)*FVprime(2,2)*FV(1,3) ! i   j   k+1   
    BP(5)= FV(1,1)*FVprime(1,2)*FV(2,3) ! i+1 j+1 k
    BP(6)= FV(1,1)*FVprime(2,2)*FV(1,3) ! i+1 j   k+1
    BP(7)= FV(2,1)*FVprime(1,2)*FV(1,3) ! i   j+1 k+1
    BP(8)= FV(1,1)*FVprime(1,2)*FV(1,3) ! i+1 j+1 k+1
    
    Return 
  End Subroutine NoDPrimeY
  !============================================================
  Subroutine OneDPrimeY(BP,FV,FVPrime)
    Implicit none
    Real(8) ,INTENT(OUT) :: BP(24)
    Real(8) ,INTENT(IN) :: FV(4,3)
    Real(8) ,INTENT(IN) :: FVprime(4,3)
    
    !  [ 8x3x2 = 48 multiplications ]

    BP(1)= FV(3,1)*FVprime(2,2)*FV(2,3) ! i'   j   k
    BP(2)= FV(3,1)*FVprime(1,2)*FV(2,3) ! i'   j+1 k
    BP(3)= FV(3,1)*FVprime(2,2)*FV(1,3) ! i'   j   k+1
    BP(4)= FV(3,1)*FVprime(1,2)*FV(1,3) ! i'   j+1 k+1
    BP(5)= FV(4,1)*FVprime(2,2)*FV(2,3) ! i+1' j   k
    BP(6)= FV(4,1)*FVprime(1,2)*FV(2,3) ! i+1' j+1 k 
    BP(7)= FV(4,1)*FVprime(2,2)*FV(1,3) ! i+1' j   k+1
    BP(8)= FV(4,1)*FVprime(1,2)*FV(1,3) ! i+1' j+1 k+1
    
    BP(9)= FV(2,1)*FVprime(3,2)*FV(2,3) ! i   j'   k
    BP(10)= FV(2,1)*FVprime(4,2)*FV(2,3) ! i   j+1' k
    BP(11)= FV(2,1)*FVprime(3,2)*FV(1,3) ! i   j'   k+1
    BP(12)= FV(2,1)*FVprime(4,2)*FV(1,3) ! i   j+1' k+1
    BP(13)= FV(1,1)*FVprime(3,2)*FV(2,3) ! i+1 j'   k
    BP(14)= FV(1,1)*FVprime(4,2)*FV(2,3) ! i+1 j+1' k 
    BP(15)= FV(1,1)*FVprime(3,2)*FV(1,3) ! i+1 j'   k+1
    BP(16)= FV(1,1)*FVprime(4,2)*FV(1,3) ! i+1 j+1' k+1
    
    BP(17)= FV(2,1)*FVprime(2,2)*FV(3,3) ! i   j   k'
    BP(18)= FV(2,1)*FVprime(1,2)*FV(3,3) ! i   j+1 k'
    BP(19)= FV(2,1)*FVprime(2,2)*FV(4,3) ! i   j   k+1'
    BP(20)= FV(2,1)*FVprime(1,2)*FV(4,3) ! i   j+1 k+1'
    BP(21)= FV(1,1)*FVprime(2,2)*FV(3,3) ! i+1 j   k'
    BP(22)= FV(1,1)*FVprime(1,2)*FV(3,3) ! i+1 j+1 k' 
    BP(23)= FV(1,1)*FVprime(2,2)*FV(4,3) ! i+1 j   k+1'
    BP(24)= FV(1,1)*FVprime(1,2)*FV(4,3) ! i+1 j+1 k+1'
    
    Return 
  End Subroutine OneDPrimeY
  !============================================================
  Subroutine TwoDPrimeY(BP,FV,FVPrime)
    Implicit none
    Real(8) ,INTENT(OUT) :: BP(24)
    Real(8) ,INTENT(IN) :: FV(4,3)
    Real(8) ,INTENT(IN) :: FVprime(4,3)    
    ! [ 8x3x2 = 48 multiplications ]
    
    BP(1)= FV(3,1)*FVprime(3,2)*FV(2,3) ! i'   j'   k
    BP(2)= FV(3,1)*FVprime(3,2)*FV(1,3) ! i'   j'   k+1
    BP(3)= FV(4,1)*FVprime(3,2)*FV(2,3) ! i+1' j'   k
    BP(4)= FV(4,1)*FVprime(3,2)*FV(1,3) ! i+1' j'   k+1
    BP(5)= FV(3,1)*FVprime(4,2)*FV(2,3) ! i'   j+1' k
    BP(6)= FV(3,1)*FVprime(4,2)*FV(1,3) ! i'   j+1' k+1
    BP(7)= FV(4,1)*FVprime(4,2)*FV(2,3) ! i+1' j+1' k
    BP(8)= FV(4,1)*FVprime(4,2)*FV(1,3) ! i+1' j+1' k+1
    
    BP(9)= FV(3,1)*FVprime(2,2)*FV(3,3)  ! i'   j   k'
    BP(10)= FV(3,1)*FVprime(2,2)*FV(4,3) ! i'   j   k+1'
    BP(11)= FV(4,1)*FVprime(2,2)*FV(3,3) ! i+1' j   k'
    BP(12)= FV(4,1)*FVprime(2,2)*FV(4,3) ! i+1' j   k+1'
    BP(13)= FV(3,1)*FVprime(1,2)*FV(3,3) ! i'   j+1 k'
    BP(14)= FV(3,1)*FVprime(1,2)*FV(4,3) ! i'   j+1 k+1'
    BP(15)= FV(4,1)*FVprime(1,2)*FV(3,3) ! i+1' j+1 k'
    BP(16)= FV(4,1)*FVprime(1,2)*FV(4,3) ! i+1' j+1 k+1'
    
    BP(17)= FV(2,1)*FVprime(3,2)*FV(3,3) ! i   j'   k'
    BP(18)= FV(2,1)*FVprime(3,2)*FV(4,3) ! i   j'   k+1'
    BP(19)= FV(1,1)*FVprime(3,2)*FV(3,3) ! i+1 j'   k'
    BP(20)= FV(1,1)*FVprime(3,2)*FV(4,3) ! i+1 j'   k+1'
    BP(21)= FV(2,1)*FVprime(4,2)*FV(3,3) ! i   j+1' k'
    BP(22)= FV(2,1)*FVprime(4,2)*FV(4,3) ! i   j+1' k+1'
    BP(23)= FV(1,1)*FVprime(4,2)*FV(3,3) ! i+1 j+1' k'
    BP(24)= FV(1,1)*FVprime(4,2)*FV(4,3) ! i+1 j+1' k+1'
    
    Return 
  End Subroutine TwoDPrimeY
  !============================================================
  Subroutine ThreeDPrimeY(BP,FV,FVPrime)
    Implicit none
    Real(8) ,INTENT(OUT) :: BP(8)
    Real(8) ,INTENT(IN) :: FV(4,3)
    Real(8) ,INTENT(IN) :: FVprime(4,3)
    !  [ 8x2 = 16 multiplications ]
    
    BP(1)= FV(3,1)*FVprime(3,2)*FV(3,3) ! i'   j'   k'
    BP(2)= FV(4,1)*FVprime(3,2)*FV(3,3) ! i+1' j'   k'
    BP(3)= FV(3,1)*FVprime(4,2)*FV(3,3) ! i'   j+1' k'
    BP(5)= FV(3,1)*FVprime(3,2)*FV(4,3) ! i'   j'   k+1'
    BP(4)= FV(4,1)*FVprime(4,2)*FV(3,3) ! i+1' j+1' k'
    BP(6)= FV(4,1)*FVprime(3,2)*FV(4,3) ! i+1' j'   k+1'
    BP(7)= FV(3,1)*FVprime(4,2)*FV(4,3) ! i'   j+1' k+1'
    BP(8)= FV(4,1)*FVprime(4,2)*FV(4,3) ! i+1' j+1' k+1'
    
    Return 
  End Subroutine ThreeDPrimeY
  !============================================================
  !============================================================
  Subroutine NoDPrimeZ(BP,FV,FVPrime)
    Implicit none
    Real(8) ,INTENT(OUT) :: BP(8)
    Real(8) ,INTENT(IN) :: FV(4,3)
    Real(8) ,INTENT(IN) :: FVprime(4,3)
    !  [ 8x2 = 16 multiplications ]

    BP(1)= FV(2,1)*FV(2,2)*FVprime(2,3) ! i   j   k
    BP(2)= FV(1,1)*FV(2,2)*FVprime(2,3) ! i+1 j   k
    BP(3)= FV(2,1)*FV(1,2)*FVprime(2,3) ! i   j+1 k
    BP(4)= FV(2,1)*FV(2,2)*FVprime(1,3) ! i   j   k+1   
    BP(5)= FV(1,1)*FV(1,2)*FVprime(2,3) ! i+1 j+1 k
    BP(6)= FV(1,1)*FV(2,2)*FVprime(1,3) ! i+1 j   k+1
    BP(7)= FV(2,1)*FV(1,2)*FVprime(1,3) ! i   j+1 k+1
    BP(8)= FV(1,1)*FV(1,2)*FVprime(1,3) ! i+1 j+1 k+1
    
    Return 
  End Subroutine NoDPrimeZ
  !============================================================
  Subroutine OneDPrimeZ(BP,FV,FVPrime)
    Implicit none
    Real(8) ,INTENT(OUT) :: BP(24)
    Real(8) ,INTENT(IN) :: FV(4,3)
    Real(8) ,INTENT(IN) :: FVprime(4,3)
    !  [ 8x3x2 = 48 multiplications ]
    
    BP(1)= FV(3,1)*FV(2,2)*FVprime(2,3) ! i'   j   k
    BP(2)= FV(3,1)*FV(1,2)*FVprime(2,3) ! i'   j+1 k
    BP(3)= FV(3,1)*FV(2,2)*FVprime(1,3) ! i'   j   k+1
    BP(4)= FV(3,1)*FV(1,2)*FVprime(1,3) ! i'   j+1 k+1
    BP(5)= FV(4,1)*FV(2,2)*FVprime(2,3) ! i+1' j   k
    BP(6)= FV(4,1)*FV(1,2)*FVprime(2,3) ! i+1' j+1 k 
    BP(7)= FV(4,1)*FV(2,2)*FVprime(1,3) ! i+1' j   k+1
    BP(8)= FV(4,1)*FV(1,2)*FVprime(1,3) ! i+1' j+1 k+1
    
    BP(9)= FV(2,1)*FV(3,2)*FVprime(2,3) ! i   j'   k
    BP(10)= FV(2,1)*FV(4,2)*FVprime(2,3) ! i   j+1' k
    BP(11)= FV(2,1)*FV(3,2)*FVprime(1,3) ! i   j'   k+1
    BP(12)= FV(2,1)*FV(4,2)*FVprime(1,3) ! i   j+1' k+1
    BP(13)= FV(1,1)*FV(3,2)*FVprime(2,3) ! i+1 j'   k
    BP(14)= FV(1,1)*FV(4,2)*FVprime(2,3) ! i+1 j+1' k 
    BP(15)= FV(1,1)*FV(3,2)*FVprime(1,3) ! i+1 j'   k+1
    BP(16)= FV(1,1)*FV(4,2)*FVprime(1,3) ! i+1 j+1' k+1
    
    BP(17)= FV(2,1)*FV(2,2)*FVprime(3,3) ! i   j   k'
    BP(18)= FV(2,1)*FV(1,2)*FVprime(3,3) ! i   j+1 k'
    BP(19)= FV(2,1)*FV(2,2)*FVprime(4,3) ! i   j   k+1'
    BP(20)= FV(2,1)*FV(1,2)*FVprime(4,3) ! i   j+1 k+1'
    BP(21)= FV(1,1)*FV(2,2)*FVprime(3,3) ! i+1 j   k'
    BP(22)= FV(1,1)*FV(1,2)*FVprime(3,3) ! i+1 j+1 k' 
    BP(23)= FV(1,1)*FV(2,2)*FVprime(4,3) ! i+1 j   k+1'
    BP(24)= FV(1,1)*FV(1,2)*FVprime(4,3) ! i+1 j+1 k+1'
    
    Return 
  End Subroutine OneDPrimeZ
  !============================================================
  Subroutine TwoDPrimeZ(BP,FV,FVPrime)
    Implicit none
    Real(8) ,INTENT(OUT) :: BP(24)
    Real(8) ,INTENT(IN) :: FV(4,3)
    Real(8) ,INTENT(IN) :: FVprime(4,3)
    !    [ 8x3x2 = 48 multiplications ]
    
    BP(1)= FV(3,1)*FV(3,2)*FVprime(2,3) ! i'   j'   k
    BP(2)= FV(3,1)*FV(3,2)*FVprime(1,3) ! i'   j'   k+1
    BP(3)= FV(4,1)*FV(3,2)*FVprime(2,3) ! i+1' j'   k
    BP(4)= FV(4,1)*FV(3,2)*FVprime(1,3) ! i+1' j'   k+1
    BP(5)= FV(3,1)*FV(4,2)*FVprime(2,3) ! i'   j+1' k
    BP(6)= FV(3,1)*FV(4,2)*FVprime(1,3) ! i'   j+1' k+1
    BP(7)= FV(4,1)*FV(4,2)*FVprime(2,3) ! i+1' j+1' k
    BP(8)= FV(4,1)*FV(4,2)*FVprime(1,3) ! i+1' j+1' k+1
    
    BP(9)= FV(3,1)*FV(2,2)*FVprime(3,3) ! i'   j   k'
    BP(10)= FV(3,1)*FV(2,2)*FVprime(4,3) ! i'   j   k+1'
    BP(11)= FV(4,1)*FV(2,2)*FVprime(3,3) ! i+1' j   k'
    BP(12)= FV(4,1)*FV(2,2)*FVprime(4,3) ! i+1' j   k+1'
    BP(13)= FV(3,1)*FV(1,2)*FVprime(3,3) ! i'   j+1 k'
    BP(14)= FV(3,1)*FV(1,2)*FVprime(4,3) ! i'   j+1 k+1'
    BP(15)= FV(4,1)*FV(1,2)*FVprime(3,3) ! i+1' j+1 k'
    BP(16)= FV(4,1)*FV(1,2)*FVprime(4,3) ! i+1' j+1 k+1'
    
    BP(17)= FV(2,1)*FV(3,2)*FVprime(3,3) ! i   j'   k'
    BP(18)= FV(2,1)*FV(3,2)*FVprime(4,3) ! i   j'   k+1'
    BP(19)= FV(1,1)*FV(3,2)*FVprime(3,3) ! i+1 j'   k'
    BP(20)= FV(1,1)*FV(3,2)*FVprime(4,3) ! i+1 j'   k+1'
    BP(21)= FV(2,1)*FV(4,2)*FVprime(3,3) ! i   j+1' k'
    BP(22)= FV(2,1)*FV(4,2)*FVprime(4,3) ! i   j+1' k+1'
    BP(23)= FV(1,1)*FV(4,2)*FVprime(3,3) ! i+1 j+1' k'
    BP(24)= FV(1,1)*FV(4,2)*FVprime(4,3) ! i+1 j+1' k+1'
    
    Return 
  End Subroutine TwoDPrimeZ
  !============================================================
  Subroutine ThreeDPrimeZ(BP,FV,FVprime)
    Implicit none
    Real(8) ,INTENT(OUT) :: BP(8)
    Real(8) ,INTENT(IN) :: FV(4,3)
    Real(8) ,INTENT(IN) :: FVprime(4,3)
    !     [ 8x2 = 16 multiplications ]
    
    BP(1)= FV(3,1)*FV(3,2)*FVprime(3,3) ! i'   j'   k'
    BP(2)= FV(4,1)*FV(3,2)*FVprime(3,3) ! i+1' j'   k'
    BP(3)= FV(3,1)*FV(4,2)*FVprime(3,3) ! i'   j+1' k'
    BP(5)= FV(3,1)*FV(3,2)*FVprime(4,3) ! i'   j'   k+1'
    BP(4)= FV(4,1)*FV(4,2)*FVprime(3,3) ! i+1' j+1' k'
    BP(6)= FV(4,1)*FV(3,2)*FVprime(4,3) ! i+1' j'   k+1'
    BP(7)= FV(3,1)*FV(4,2)*FVprime(4,3) ! i'   j+1' k+1'
    BP(8)= FV(4,1)*FV(4,2)*FVprime(4,3) ! i+1' j+1' k+1'
    
    Return 
  End Subroutine ThreeDPrimeZ    
End Module BasePolynomials
