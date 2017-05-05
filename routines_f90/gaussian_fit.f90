module Gaussian_Fit
  use common
  use levenberg_marquardt
  implicit none

  private
  public :: gausfit,get_gaus,gausfit_asy,get_gaus_asy
  
  real (PR), dimension(:), allocatable :: xd,yd,zd,wd
  real (PR) :: fixback_value, fixasy_value
  logical :: fixback_flag,fixasy_flag
contains
  
  real (PR) function d1mach(I)
    integer,intent(IN):: I
    !
    !  DOUBLE-PRECISION MACHINE CONSTANTS
    !  D1MACH( 1) = B**(EMIN-1), THE SMALLEST POSITIVE MAGNITUDE.
    !  D1MACH( 2) = B**EMAX*(1 - B**(-T)), THE LARGEST MAGNITUDE.
    !  D1MACH( 3) = B**(-T), THE SMALLEST RELATIVE SPACING.
    !  D1MACH( 4) = B**(1-T), THE LARGEST RELATIVE SPACING.
    !  D1MACH( 5) = LOG10(B)
    !
    real(PR),dimension(5),save::X
    logical,save:: ran = .false.
    
    if(.not. ran)then
       X(1) = tiny(X(1))
       X(2) = huge(X(2))
       X(3) = epsilon(X(3))/radix(X(3))
       X(4) = epsilon(X(4))
       X(5) = log10(real(radix(X(5)),PR))
       ran = .true.
    end if

    if (I < 1 .or. I > 5) then
       write(*,*) 'D1MACH(I): I =',I,' is out of bounds.'
       stop
    end if
    D1MACH = X(I)
    return
  end function d1mach
  
  subroutine gausfit(xdata,ydata,wdata,cntr,sgm,ampl,back,fixback,verbose)
    integer :: m
    integer, optional, intent(in) :: fixback,verbose
    !real, dimension(:), intent(in), target ::xdata,ydata
    real (PR), dimension(:), intent(in) ::xdata,ydata,wdata
    real (PR), intent(inout) ::  cntr,sgm,ampl,back
    
    integer :: N,INFO
    integer, dimension(4) ::  IWA
    real (PR) :: TOL
    REAL (PR),dimension(4) ::x,diag,sig
    REAL (PR),dimension(size(xdata)) ::fvec
    
    if (size(ydata).ne.size(xdata)) stop
    if (size(wdata).ne.size(xdata)) stop
    
    fixback_flag = .false.
    fixback_value = back
    if (present(fixback)) fixback_flag = .true.
    
    X(1)=  cntr 
    X(2)=  sgm  
    X(3)=  ampl 
    X(4)=  back 
    m = size(xdata)
    allocate(xd(m),yd(m),wd(m)) ;  xd=xdata ;  yd=ydata ; wd = wdata
    
    fvec = 0._PR; diag = 0._PR; sig = 0._PR
    TOL = D1MACH(4) 
    INFO=0 
    N = 4
    
    call lmdif2(gausfit_f,m,n,x,fvec,tol,diag,sig,info,iwa)
    
    cntr = X(1)
    sgm  = X(2)
    ampl = X(3)
    back = X(4)
    if (fixback_flag) back = fixback_value
    
    
    deallocate(xd,yd,wd)
    if(present(verbose)) write(verbose,'(A25,F9.4,A10,F9.4,A9,F15.4,A9,F15.4)') &
         '   *** FIT4PGN_s: CNTR = ',cntr,'  SIGMA = ',sgm,'  AMPL = ',ampl,'  BACK = ',back

    return
  end subroutine gausfit


  SUBROUTINE gausfit_f(m,n,x,fvec,iflag)
         use common
         IMPLICIT NONE
         !INTEGER,PARAMETER::dp = SELECTED_REAL_KIND(12,60)
         INTEGER,INTENT(IN)::m,n
         REAL (PR),INTENT(IN)::x(:)
         REAL (PR),INTENT(OUT)::fvec(:)
         INTEGER,INTENT(inout)::iflag

         INTEGER :: k
         REAL (PR) :: x1,x2,x3,x4, gaus
         
         x1 = X(1)
         x2 = X(2)
         x3 = X(3)
         x4 = X(4)
         if (fixback_flag) x4 = fixback_value
         do k=1,M
            gaus=X3*exp((-(xd(k)-X1)**2)/(2._PR*X2**2))+ X4
            FVEC(k)= (yd(k) - gaus)*wd(k)!*(abs((X(3)-yd(k))/(X(3)-X(4))))**2
         enddo
         return
  end subroutine gausfit_f
  
  real (PR) function get_gaus(x,cntr,sgm,ampl,back)
    real (PR), intent(in) :: x, cntr, sgm, ampl, back
    get_gaus = back + ampl*exp((-(x-cntr)**2)/(2._PR*sgm**2)) 
    return
  end function get_gaus

 !! ************** !!

  subroutine gausfit_asy(xdata,ydata,wdata,cntr,fwhm,ampl,back,asy,fixasy,fixback,verbose)
    !! Asymmetric Gaussian fit as stated in Nardetto et al. 2006 and Fgueira et al. 2013
    !! For consistency with previous formulations of gaussian fit, the negative sign 
    !! in front of the equation has been dropped (D, or the amplitude, will be negative)
    integer :: m
    integer, optional, intent(in) :: fixasy,fixback,verbose
    !real, dimension(:), intent(in), target ::xdata,ydata
    real (PR), dimension(:), intent(in) ::xdata,ydata,wdata
    real (PR), intent(inout) ::  cntr,fwhm,ampl,back,asy
    
    integer :: N,INFO
    integer, dimension(5) ::  IWA
    real (PR) :: TOL
    REAL (PR),dimension(5) ::x,diag,sig
    REAL (PR),dimension(size(xdata)) ::fvec
    

    if (size(ydata).ne.size(xdata)) stop
    if (size(wdata).ne.size(xdata)) stop
    
    fixback_flag = .false.
    fixasy_flag  = .false.
    fixback_value = back
    fixasy_value = asy
    if (present(fixback)) fixback_flag = .true.
    if (present(fixasy)) fixasy_flag = .true.
    
    X(1)=  cntr 
    X(2)=  fwhm  
    X(3)=  ampl 
    X(4)=  back 
    X(5)=  asy 
    m = size(xdata)
    allocate(xd(m),yd(m),wd(m)) ;  xd=xdata ;  yd=ydata ; wd = wdata
    
    fvec = 0._PR; diag = 0._PR; sig = 0._PR
    TOL = D1MACH(4) 
    INFO=0 
    N = 5
    
    call lmdif2(gausfit_asy_f,m,n,x,fvec,tol,diag,sig,info,iwa)
    
    cntr = X(1)
    fwhm = X(2)
    ampl = X(3)
    back = X(4)
    asy  = X(5)
    if (fixback_flag) back = fixback_value
    if (fixasy_flag)  asy = fixasy_value
    
    
    deallocate(xd,yd,wd)
    if(present(verbose)) write(verbose,'(A25,F9.4,A10,F9.4,3(A9,F15.4))') &
         '   *** FIT_ASY_s: CNTR = ',cntr,'   FWHM = ',fwhm,'  AMPL = ',ampl,'  BACK = ',back,'  ASY = ',asy

    return
  end subroutine gausfit_asy

  SUBROUTINE gausfit_asy_f(m,n,x,fvec,iflag)
    use common
    IMPLICIT NONE
    !INTEGER,PARAMETER::dp = SELECTED_REAL_KIND(12,60)
    INTEGER,INTENT(IN)::m,n
    REAL (PR),INTENT(IN)::x(:)
    REAL (PR),INTENT(OUT)::fvec(:)
    INTEGER,INTENT(inout)::iflag
    
    INTEGER :: k
    REAL (PR) :: x1,x2,x3,x4,x5, gaus
    
    x1 = X(1)  !! Center, RV_{center}
    x2 = X(2)  !! FWHM
    x3 = X(3)  !! Amplitude, -D
    x4 = X(4)  !! Background level, C
    x5 = X(5)  !! Asy, A
    if (fixback_flag) x4 = fixback_value
    if (fixasy_flag)  x5 = fixasy_value
    do k=1,M
       if (xd(k).le.X1) gaus = x3*exp( 2._PR*log(2._PR)*(-(xd(k)-X1)**2)/( x2*(1._PR-x5))**2) + x4
       if (xd(k).gt.X1) gaus = x3*exp( 2._PR*log(2._PR)*(-(xd(k)-X1)**2)/( x2*(1._PR+x5))**2) + x4
       FVEC(k)= (yd(k) - gaus)*wd(k)!*(abs((X(3)-yd(k))/(X(3)-X(4))))**2
    enddo
    return
  end subroutine gausfit_asy_f
  
  real (PR) function get_gaus_asy(x,cntr,fwhm,ampl,back,asy)
    real (PR), intent(in) :: x, cntr, fwhm, ampl, back, asy
    
    if (x.le.cntr) get_gaus_asy = ampl*exp( 2._PR*log(2._PR)*(-(x-cntr)**2)/( fwhm*(1._PR-asy))**2) + back
    if (x.gt.cntr) get_gaus_asy = ampl*exp( 2._PR*log(2._PR)*(-(x-cntr)**2)/( fwhm*(1._PR+asy))**2) + back

    return
  end function get_gaus_asy
  
end module Gaussian_Fit

