module Poly_Fit
  use common
  use levenberg_marquardt_mod5V
  implicit none

  private
  public :: polyfit, get_poly
  
contains
  
  real (kind=8) function d1mach(I)
    integer,intent(IN):: I
    !
    !  DOUBLE-PRECISION MACHINE CONSTANTS
    !  D1MACH( 1) = B**(EMIN-1), THE SMALLEST POSITIVE MAGNITUDE.
    !  D1MACH( 2) = B**EMAX*(1 - B**(-T)), THE LARGEST MAGNITUDE.
    !  D1MACH( 3) = B**(-T), THE SMALLEST RELATIVE SPACING.
    !  D1MACH( 4) = B**(1-T), THE LARGEST RELATIVE SPACING.
    !  D1MACH( 5) = LOG10(B)
    !
    real(KIND=8),dimension(5),save::X
    logical,save:: ran = .false.
    
    if(.not. ran)then
       X(1) = tiny(X(1))
       X(2) = huge(X(2))
       X(3) = epsilon(X(3))/radix(X(3))
       X(4) = epsilon(X(4))
       X(5) = log10(real(radix(X(5)),KIND=8))
       ran = .true.
    end if

    if (I < 1 .or. I > 5) then
       write(*,*) 'D1MACH(I): I =',I,' is out of bounds.'
       stop
    end if
    D1MACH = X(I)
    return
  end function d1mach
  
  subroutine polyfit(xdata,ydata,wdata,coeff,verbose)
    integer :: m
    integer, optional, intent(in) :: verbose
    !real, dimension(:), intent(in), target ::xdata,ydata
    real (PR), dimension(:), intent(in) ::xdata,ydata,wdata
    real (PR), dimension(size(xdata,dim=1)) :: dummy1, dummy2
    real (PR), intent(inout), dimension(:) :: coeff
    
    integer :: N,INFO
    integer, dimension(size(coeff)) ::  IWA
    real (PR) :: TOL
    REAL (PR),dimension(size(coeff)) ::x,diag,sig
    REAL (PR),dimension(size(xdata)) ::fvec
    
    if (size(ydata).ne.size(xdata)) stop
    if (size(wdata).ne.size(xdata)) stop
    
    x = coeff 
    
    m = size(xdata)
    
    fvec = 0._PR; diag = 0._PR; sig = 0._PR
    TOL = D1MACH(4) 
    INFO=0 
    N = size(coeff)
    
    call lmdif2_mod5V(polyfit_f,m,n,x,xdata,ydata,wdata,dummy1,dummy2,fvec,tol,diag,sig,info,iwa)
    
    coeff = x
    
    !else
    !   cntr = 0.0D0 ; sgm  = -9.9D0 ;  ampl = -999.9D0 ;  back = -999.9D0
    !end if
    

    return
  end subroutine polyfit
  
  SUBROUTINE polyfit_f(m,n,x,xd,yd,wd,dummy1,dummy2,fvec,iflag)
    use common
    IMPLICIT NONE
    !INTEGER,PARAMETER::dp = SELECTED_REAL_KIND(12,60)
    INTEGER,INTENT(IN)::m,n
    REAL (PR),INTENT(IN)::x(:)
    REAL (PR),INTENT(OUT)::fvec(:)
    INTEGER,INTENT(inout)::iflag
    REAL (PR), intent(in), dimension(:) :: xd,yd,wd,dummy1,dummy2

    INTEGER :: k, t
    REAL (PR) :: pol

    do k=1,M
       pol = X(n)
       do t=n-1,1,-1
          pol = pol*xd(k)+X(t)
       end do
       FVEC(k)= (yd(k) - pol)*wd(k)
    enddo
    return
  end subroutine polyfit_f
  
  
  real (PR) function get_poly(x,coeff)
    real (PR), intent(in) :: x
    real (PR), dimension(:), intent(in) :: coeff
    integer ::n,t
    n=size(coeff)
    get_poly = coeff(n)
    do t=n-1,1,-1
       get_poly = get_poly*x+coeff(t)
    end do
    return
  end function get_poly
  
  


  
end module Poly_Fit

