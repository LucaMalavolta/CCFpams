module chebyshev_fit
  use common
  use levenberg_marquardt
  implicit none

  private
  public ::   cheby_fit, get_cheby, get_cheby_1der, get_cheby_2der, &
       cheby_2d_fit, get_cheby_2d, cheby_3d_fit, get_cheby_3d, cheby_4d_fit, get_cheby_4d

  real (PR), dimension(:), allocatable :: xd,yd,zd,td,fd,wd
  real (PR), dimension(:,:), allocatable :: Tc
  real (PR), dimension(:,:), allocatable :: TcXl,TcYl,TcZl,TcTl
  !!!!! ********** CHEBY 1D fit ********** !!!!!

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



  subroutine cheby_fit(xdata,ydata,wdata,order,coeff,verbose)
    integer, intent(in) :: order
    real (PR), dimension(:), intent(in) ::xdata,ydata,wdata
    real (PR), dimension(order+1), intent(out) :: coeff
    integer, optional, intent(in) :: verbose


    integer :: N,INFO
    integer, dimension(order+1) ::  IWA
    real (PR) :: TOL
    REAL (PR),dimension(order+1) ::x,diag,sig
    REAL (PR),dimension(size(xdata)) ::fvec

    integer :: m, nm, nk

    if (size(xdata).ne.size(ydata)) stop

    m = size(xdata)
    if (m.lt.order) stop
    allocate(xd(m),yd(m),wd(m)) ; xd=xdata ;  yd=ydata ;  wd = wdata

    allocate(Tc(order,m))

    ! determining first approximation of the coefficients
    ! X(k) = c_k = 2/N*sum_{1}^{N} (f(x_n)*T_{k-1}(x_n))
    ! N = number of points, c_k order
    ! X(n) id the coefficient for T(n-1), and it has size order+1!!!

    X = 0.00D0

    do nm = 1, m
       Tc(1,nm) = xd(nm)
       Tc(2,nm) = cheb_Tmp(xd(nm),Tc(1,nm),1.00D0)
       X(1) = X(1) + yd(nm)/real(m)
       X(2) = X(2) + yd(nm)*Tc(1,nm)/real(m)
       X(3) = X(3) + yd(nm)*Tc(2,nm)/real(m)
       do nk = 2, order-1
          Tc(nk+1,nm) = cheb_Tmp(xd(nm),Tc(nk,nm),Tc(nk-1,nm))
          X(nk+2) = X(nk+2) +  yd(nm)*Tc(nk+1,nm)/real(m)
       end do
    end do
    X = 2.0000D0 * X


    fvec = 0._PR; diag = 0._PR; sig = 0._PR
    TOL = D1MACH(4)
    INFO=0
    N = order+1

    call lmdif2(chebyshev_f,m,n,x,fvec,tol,diag,sig,info,iwa)

    coeff = X

    deallocate(xd,yd,wd,Tc)
    return
  end subroutine cheby_fit

  subroutine  chebyshev_f(M,N,X,FVEC,IFLAG)
    INTEGER,INTENT(IN)::m,n
    REAL (PR),INTENT(IN)::x(:)
    REAL (PR),INTENT(OUT)::fvec(:)
    INTEGER,INTENT(inout)::iflag
    real (PR) :: c_fit
    integer :: k

    do k=1,M
       c_fit = sum (Tc(:,k)*X(2:n)) - 0.500000D0*X(1)
       FVEC(k)= (yd(k) - c_fit)*wd(k)
    enddo
    return
  end subroutine chebyshev_f


  real (PR) function get_cheby(x,coeff)
    real (PR), intent(in) :: x
    real (PR), dimension(:), intent(in) :: coeff
    real (PR), dimension(size(coeff,1)-1) :: T
    integer :: n

    get_cheby = 0.00D0
    T(1) = x
    T(2) = cheb_Tmp(x,T(1),1.000D0)
    do n=2,size(coeff,1)-2
       T(n+1) =  cheb_Tmp(x,T(n),T(n-1))
    end do
    get_cheby = sum (T(:)*coeff(2:n)) - 0.500000D0*coeff(1)
    return
  end function get_cheby

  real (PR) function cheb_Tmp_1der(x,T0_n,T1_n,T1_nm)
    real (PR) :: x, T0_n, T1_n, T1_nm
    cheb_Tmp_1der = 2.000000D0*(T0_n + x*T1_n) + T1_nm
    return
  end function cheb_Tmp_1der

  real (PR) function get_cheby_1der(x,coeff)
    real (PR), intent(in) :: x
    real (PR), dimension(:), intent(in) :: coeff
    real (PR), dimension(size(coeff,1)-1) :: T0,T1
    integer :: n
    get_cheby_1der = 0.0d0
    T0(1) = x
    T1(1) = 1.000D0
    T0(2) = cheb_Tmp(x,T0(1),1.000D0)
    T1(2) = cheb_Tmp_1der(x,T0(1),T1(1),0.000D0)
    do n=2,size(coeff,1)-2
       T0(n+1) =  cheb_Tmp(x,T0(n),T0(n-1))
       T1(n+1) =  cheb_Tmp_1der(x,T0(n),T1(n),T1(n-1))
    end do
    get_cheby_1der = sum (T1(:)*coeff(2:n))
    return
  end function get_cheby_1der

  real (PR) function cheb_Tmp_2der(x,T1_n,T2_n,T2_nm)
    real (PR) :: x, T1_n, T2_n, T2_nm
    cheb_Tmp_2der = 2.000000D0*x*T2_n + 4.000000D0*T1_n + T2_nm
    return
  end function cheb_Tmp_2der

  real (PR) function get_cheby_2der(x,coeff)
    real (PR), intent(in) :: x
    real (PR), dimension(:), intent(in) :: coeff
    real (PR), dimension(size(coeff,1)-1) :: T0,T1,T2
    integer :: n
    get_cheby_2der = 0.0d0
    T0(1) = x
    T1(1) = 1.000D0
    T2(1) = 0.000D0
    T0(2) = cheb_Tmp(x,T0(1),1.000D0)
    T1(2) = cheb_Tmp_1der(x,T0(1),T1(1),0.000D0)
    T2(2) = cheb_Tmp_2der(x,T1(1),T2(1),0.000D0)
    do n=2,size(coeff,1)-2
       T0(n+1) =  cheb_Tmp(x,T0(n),T0(n-1))
       T1(n+1) =  cheb_Tmp_1der(x,T0(n),T1(n),T1(n-1))
       T2(n+1) =  cheb_Tmp_2der(x,T1(n),T2(n),T2(n-1))
    end do
    get_cheby_2der = sum (T2(:)*coeff(2:n))
    return
  end function get_cheby_2der

  !!!!! ********** CHEBY 2D fit ********** !!!!!


  subroutine cheby_2d_fit(xdata,ydata,zdata,wdata,order_x,order_y,coeff,verbose)
    integer, intent(in) :: order_x, order_y
    real (PR), dimension(:), intent(in) ::xdata,ydata,zdata,wdata
    real (PR), dimension(order_x+1,order_y+1), intent(out) :: coeff
    integer, optional, intent(in) :: verbose

    integer :: N, INFO
    integer, dimension((order_x+1)*(order_y+1)) :: IWA
    real (PR) :: TOL
    real (PR), dimension(size(xdata)) :: FVEC
    real (PR), dimension((order_x+1)*(order_y+1)) :: X,DIAG,SIG

    integer :: m, nm, nk, nx, ny, tot_order
    real (PR), dimension(:,:), allocatable :: TcX,TcY


    if (size(xdata).ne.size(ydata)) stop
    if (size(xdata).ne.size(zdata)) stop
    tot_order = (order_x+1)*(order_y+1)
    m = size(xdata)

    if (verbose.gt.0) write(verbose,'(A,I10)') ' CHEBY_3D_FIT: M         = ',m
    if (verbose.gt.0) write(verbose,'(A,I10)') ' CHEBY_3D_FIT: order_x   = ',order_x
    if (verbose.gt.0) write(verbose,'(A,I10)') ' CHEBY_3D_FIT: order_y   = ',order_y
    if (verbose.gt.0) write(verbose,'(A,I10)') ' CHEBY_3D_FIT: tot_order = ',tot_order

    if (m.lt.tot_order) stop
    allocate(xd(m),yd(m),zd(m),wd(m)) ; xd=xdata ;  yd=ydata ;  zd = zdata ;  wd = wdata

    if (verbose.gt.0) write(verbose,*) 'CHEBY_3D_FIT: VARIABLE allocation'
    allocate(TcX(order_x+1,m),TcY(order_y+1,m))
    allocate(TcXl(tot_order,m),TcYl(tot_order,m))

    ! N = number of points, c_k order
    ! Here I choose  a more consistent representation, since the other one
    ! was computationally hard to implement (a lot of ifs...)
    ! now a_ij, T_i and T_j have the same suffix meaning....

    X = 0.00D0
    if (verbose.gt.0) write(verbose,*) 'CHEBY_3D_FIT: Computing first approximation of coefficients'

    do nm = 1, m

       TcX(1,nm) = 1.0000D0
       TcX(2,nm) = xd(nm)
       do nk = 2, order_x
          TcX(nk+1,nm) = cheb_Tmp(xd(nm),TcX(nk,nm),TcX(nk-1,nm))
       end do

       TcY(1,nm) = 1.0000d0
       TcY(2,nm) = yd(nm)
       do nk = 2, order_y
          TcY(nk+1,nm) = cheb_Tmp(yd(nm),TcY(nk,nm),TcY(nk-1,nm))

       end do

       nk = 0
       do nx = 1, order_x+1
          do ny = 1, order_y+1
             nk = nk + 1; TcXl(nk,nm) = TcX(nx,nm) ; TcYl(nk,nm) = TcY(ny,nm)

          end do
       end do
       X = X + zd(nm)*(TcXl(:,nm)*TcYl(:,nm)/real(m)*real(4))
    end do

    if (verbose.gt.0) write(verbose,*) 'CHEBY_3D_FIT: Deallocation of temporary variables'
    deallocate(TcX,TcY)

!       X(1) = X(1) + yd(nm)/real(m)
!       X(2) = X(2) + yd(nm)*Tc(1,nm)/real(m)
!       X(3) = X(3) + yd(nm)*Tc(2,nm)/real(m)
!       ! write(*,*) X(1),X(2),X(3)
!       do nk = 2, order_x-1
!          X(nk+2) = X(nk+2) +  yd(nm)*Tc(nk+1,nm)/real(m)
!       end do

    if (verbose.gt.0) write(verbose,*) 'CHEBY_3D_FIT: CALLING lmdif subroutine'


    fvec = 0._PR; diag = 0._PR; sig = 0._PR
    TOL = D1MACH(4)
    INFO=0
    N = (order_x+1)*(order_y+1)

    call lmdif2(cheby_2d_f,m,n,x,fvec,tol,diag,sig,info,iwa)


    if (verbose.gt.0) write(verbose,*) 'CHEBY_3D_FIT: Saving the coefficients'

    nk = 0
    do nx = 1, order_x+1
       do ny = 1, order_y+1
          nk = nk + 1; coeff(nx,ny) = X(nk)
       end do
    end do
    if (verbose.gt.0) write(verbose,*) 'CHEBY_3D_FIT: Variable deallocation'

    deallocate(xd,yd,zd,wd,TcXl,TcYl)
    return
  end subroutine cheby_2d_fit

  subroutine  cheby_2d_f(M,N,X,FVEC,IFLAG)
    INTEGER,INTENT(IN)::m,n
    REAL (PR),INTENT(IN)::x(:)
    REAL (PR),INTENT(OUT)::fvec(:)
    INTEGER,INTENT(inout)::iflag
    real (PR) :: c_fit
    integer :: k

    do k=1,M
       c_fit = sum(X(:)*TcXl(:,k)*TcYl(:,k))
       FVEC(k)= (zd(k) - c_fit)*wd(k)
    enddo
    return
  end subroutine cheby_2d_f


  real (PR) function get_cheby_2d(x,y,coeff)
    real (PR), intent(in) :: x,y
    real (PR), dimension(:,:), intent(in) :: coeff
    real (PR), dimension(size(coeff,1)) :: Tx
    real (PR), dimension(size(coeff,2)) :: Ty
    integer :: n, nx, ny

    get_cheby_2d = 0.00D0

    Tx(1) = 1
    Tx(2) = x
    Ty(1) = 1
    Ty(2) = y

    do n=3,size(coeff,1)
       Tx(n) =  cheb_Tmp(x,Tx(n-1),Tx(n-2))
    end do
    do n=3,size(coeff,2)
       Ty(n) =  cheb_Tmp(y,Ty(n-1),Ty(n-2))
    end do

    do nx = 1,size(coeff,1)
       do ny = 1, size(coeff,2)
          get_cheby_2d = get_cheby_2d + coeff(nx,ny)*Tx(nx)*Ty(ny)
       end do
    end do

    return
  end function get_cheby_2d


  !!!!! ********** CHEBY 3D fit ********** !!!!!


  subroutine cheby_3d_fit(xdata,ydata,zdata,fdata,wdata,order_x,order_y,order_z,coeff,verbose)
    integer, intent(in) :: order_x, order_y, order_z
    real (PR), dimension(:), intent(in) ::xdata,ydata,zdata,fdata,wdata
    real (PR), dimension(order_x+1,order_y+1,order_z+1), intent(out) :: coeff
    integer, optional, intent(in) :: verbose

    integer :: N, INFO
    real (PR) :: TOL
    integer , dimension((order_x+1)*(order_y+1)*(order_z+1)) :: IWA
    real (PR), dimension(size(xdata)) :: FVEC
    real (PR), dimension((order_x+1)*(order_y+1)*(order_z+1)) :: X,DIAG,SIG

    integer :: m, nm, nk, nx, ny, nz, tot_order
    real (PR), dimension(:,:), allocatable :: TcX,TcY,TcZ


    if (size(xdata).ne.size(ydata)) stop
    if (size(xdata).ne.size(zdata)) stop
    if (size(xdata).ne.size(fdata)) stop
    tot_order = (order_x+1)*(order_y+1)*(order_z+1)
    m = size(xdata)
    if (verbose.gt.0) write(verbose,'(A,I10)') ' CHEBY_3D_FIT: M         = ',m
    if (verbose.gt.0) write(verbose,'(A,I10)') ' CHEBY_3D_FIT: order_x   = ',order_x
    if (verbose.gt.0) write(verbose,'(A,I10)') ' CHEBY_3D_FIT: order_y   = ',order_y
    if (verbose.gt.0) write(verbose,'(A,I10)') ' CHEBY_3D_FIT: order_z   = ',order_z
    if (verbose.gt.0) write(verbose,'(A,I10)') ' CHEBY_3D_FIT: tot_order = ',tot_order

    if (m.lt.tot_order) stop
    allocate(xd(m),yd(m),zd(m),fd(m),wd(m))
    xd=xdata ;  yd=ydata ;  zd = zdata ;  fd = fdata ;  wd = wdata

    if (verbose.gt.0) write(verbose,*) 'CHEBY_3D_FIT: VARIABLE allocation'
    allocate(TcX(order_x+1,m),TcY(order_y+1,m),TcZ(order_z+1,m))
    allocate(TcXl(tot_order,m),TcYl(tot_order,m),TcZl(tot_order,m))

    ! N = number of points, c_k order
    ! Here I choose  a more consistent representation, since the other one
    ! was computationally hard to implement (a lot of ifs...)
    ! now a_ij, T_i and T_j have the same suffix meaning....

    X = 0.00D0
    if (verbose.gt.0) write(verbose,*) 'CHEBY_3D_FIT: Computing first approximation of coefficients'

    do nm = 1, m

       TcX(1,nm) = 1.0000D0
       TcX(2,nm) = xd(nm)
       do nk = 2, order_x
          TcX(nk+1,nm) = cheb_Tmp(xd(nm),TcX(nk,nm),TcX(nk-1,nm))
       end do

       TcY(1,nm) = 1.0000d0
       TcY(2,nm) = yd(nm)
       do nk = 2, order_y
          TcY(nk+1,nm) = cheb_Tmp(yd(nm),TcY(nk,nm),TcY(nk-1,nm))
       end do

       TcZ(1,nm) = 1.0000d0
       TcZ(2,nm) = zd(nm)
       do nk = 2, order_z
          TcZ(nk+1,nm) = cheb_Tmp(zd(nm),TcZ(nk,nm),TcZ(nk-1,nm))
       end do


       nk = 0
       do nx = 1, order_x+1
          do ny = 1, order_y+1
             do nz = 1, order_z+1
                nk = nk + 1; TcXl(nk,nm) = TcX(nx,nm) ; TcYl(nk,nm) = TcY(ny,nm)
                TcZl(nk,nm) = TcZ(nz,nm)
             end do
          end do
       end do
       X = X + fd(nm)*(TcXl(:,nm)*TcYl(:,nm)*TcZl(:,nm)/real(m)*real(8)) !???? real(8)
    end do

    if (verbose.gt.0) write(verbose,*) 'CHEBY_3D_FIT: Deallocation of temporary variables'
    deallocate(TcX,TcY,TcZ)

    if (verbose.gt.0) write(verbose,*) 'CHEBY_3D_FIT: CALLING lmdif subroutine'

    fvec = 0._PR; diag = 0._PR; sig = 0._PR
    TOL = D1MACH(4)
    N=(order_x+1)*(order_y+1)*(order_z+1) ;  INFO=0

    call lmdif2(cheby_3d_f,m,n,x,fvec,tol,diag,sig,info,iwa)

    if (verbose.gt.0) write(verbose,*) 'CHEBY_3D_FIT: Saving the coefficients'

    nk = 0
    do nx = 1, order_x+1
       do ny = 1, order_y+1
          do nz = 1, order_z+1
             nk = nk + 1; coeff(nx,ny,nz) = X(nk)
          end do
       end do
    end do
    if (verbose.gt.0) write(verbose,*) 'CHEBY_3D_FIT: Variable deallocation'

    deallocate(xd,yd,zd,fd,wd,TcXl,TcYl,TcZl)
    return
  end subroutine cheby_3d_fit

  subroutine  cheby_3d_f(M,N,X,FVEC,IFLAG)
    INTEGER,INTENT(IN)::m,n
    REAL (PR),INTENT(IN)::x(:)
    REAL (PR),INTENT(OUT)::fvec(:)
    INTEGER,INTENT(inout)::iflag
    real (PR) :: c_fit
    integer :: k

    do k=1,M
       c_fit = sum(X(:)*TcXl(:,k)*TcYl(:,k)*TcZl(:,k))
       FVEC(k)= (fd(k) - c_fit)*wd(k)
    enddo
    return
  end subroutine cheby_3d_f


  real (PR) function get_cheby_3d(x,y,z,coeff)
    real (PR), intent(in) :: x,y,z
    real (PR), dimension(:,:,:), intent(in) :: coeff
    real (PR), dimension(size(coeff,1)) :: Tx
    real (PR), dimension(size(coeff,2)) :: Ty
    real (PR), dimension(size(coeff,3)) :: Tz
    integer :: n, nx, ny, nz

    get_cheby_3d = 0.00D0

    Tx(1) = 1
    Tx(2) = x
    Ty(1) = 1
    Ty(2) = y
    Tz(1) = 1
    Tz(2) = z

    do n=3,size(coeff,1)
       Tx(n) =  cheb_Tmp(x,Tx(n-1),Tx(n-2))
    end do
    do n=3,size(coeff,2)
       Ty(n) =  cheb_Tmp(y,Ty(n-1),Ty(n-2))
    end do
    do n=3,size(coeff,3)
       Tz(n) =  cheb_Tmp(z,Tz(n-1),Tz(n-2))
    end do

    do nx = 1,size(coeff,1)
       do ny = 1, size(coeff,2)
          do nz = 1, size(coeff,3)
             get_cheby_3d = get_cheby_3d + coeff(nx,ny,nz)*Tx(nx)*Ty(ny)*Tz(nz)
          end do
       end do
    end do

    return
  end function get_cheby_3d

  !!!!! ********** CHEBY 4D fit ********** !!!!!

  subroutine cheby_4d_fit(xdata,ydata,zdata,tdata,fdata,wdata,order_x,order_y,order_z,order_t,coeff,verbose)
    integer, intent(in) :: order_x, order_y, order_z, order_t
    real (PR), dimension(:), intent(in) ::xdata,ydata,zdata,tdata,fdata,wdata
    real (PR), dimension(order_x+1,order_y+1,order_z+1,order_t+1), intent(out) :: coeff
    integer, optional, intent(in) :: verbose

    integer :: N,INFO
    integer, dimension((order_x+1)*(order_y+1)*(order_z+1)*(order_t+1)) :: IWA
    real (PR) :: TOL
    real (PR), dimension(size(xdata)) :: FVEC
    real (PR), dimension((order_x+1)*(order_y+1)*(order_z+1)*(order_t+1)) :: X,DIAG,SIG

    integer :: m, nm, nk, nx, ny, nz, nt, tot_order
    real (PR), dimension(:,:), allocatable :: TcX,TcY,TcZ,TcT


    if (size(xdata).ne.size(ydata)) stop
    if (size(xdata).ne.size(zdata)) stop
    if (size(xdata).ne.size(tdata)) stop
    if (size(xdata).ne.size(fdata)) stop
    if (size(xdata).ne.size(wdata)) stop
    tot_order = (order_x+1)*(order_y+1)*(order_z+1)*(order_t+1)
    m = size(xdata)
    if (verbose.gt.0) write(verbose,'(A,I10)') ' CHEBY_3D_FIT: M         = ',m
    if (verbose.gt.0) write(verbose,'(A,I10)') ' CHEBY_3D_FIT: order_x   = ',order_x
    if (verbose.gt.0) write(verbose,'(A,I10)') ' CHEBY_3D_FIT: order_y   = ',order_y
    if (verbose.gt.0) write(verbose,'(A,I10)') ' CHEBY_3D_FIT: order_z   = ',order_z
    if (verbose.gt.0) write(verbose,'(A,I10)') ' CHEBY_3D_FIT: order_t   = ',order_t
    if (verbose.gt.0) write(verbose,'(A,I10)') ' CHEBY_3D_FIT: tot_order = ',tot_order
    if (verbose.gt.0) write(verbose,'(A,I10)') ' M=        ',m
    if (verbose.gt.0) write(verbose,'(A,I10)') ' tot_order=',tot_order

    if (m.lt.tot_order) stop
    allocate(xd(m),yd(m),zd(m),td(m),fd(m),wd(m))
    xd=xdata ;  yd=ydata ;  zd = zdata ; td = tdata ; fd = fdata ; wd = wdata

    if (verbose.gt.0) write(verbose,*) 'CHEBY_3D_FIT: VARIABLE allocation'
    allocate(TcX(order_x+1,m),TcY(order_y+1,m),TcZ(order_z+1,m),TcT(order_t+1,m))
    allocate(TcXl(tot_order,m),TcYl(tot_order,m),TcZl(tot_order,m),TcTl(tot_order,m))

    ! N = number of points, c_k order
    ! Here I choose  a more consistent representation, since the other one
    ! was computationally hard to implement (a lot of ifs...)
    ! now a_ij, T_i and T_j have the same suffix meaning....

    X = 0.00D0
    if (verbose.gt.0) write(verbose,*) 'CHEBY_3D_FIT: Computing first approximation of coefficients'

    do nm = 1, m

       TcX(1,nm) = 1.0000D0
       TcX(2,nm) = xd(nm)
       do nk = 2, order_x
          TcX(nk+1,nm) = cheb_Tmp(xd(nm),TcX(nk,nm),TcX(nk-1,nm))
       end do

       TcY(1,nm) = 1.0000d0
       TcY(2,nm) = yd(nm)
       do nk = 2, order_y
          TcY(nk+1,nm) = cheb_Tmp(yd(nm),TcY(nk,nm),TcY(nk-1,nm))
       end do

       TcZ(1,nm) = 1.0000d0
       TcZ(2,nm) = zd(nm)
       do nk = 2, order_z
          TcZ(nk+1,nm) = cheb_Tmp(zd(nm),TcZ(nk,nm),TcZ(nk-1,nm))
       end do

       TcT(1,nm) = 1.0000d0
       TcT(2,nm) = td(nm)
       do nk = 2, order_t
          TcT(nk+1,nm) = cheb_Tmp(td(nm),TcT(nk,nm),TcT(nk-1,nm))
       end do

       nk = 0
       do nx = 1, order_x+1
          do ny = 1, order_y+1
             do nz = 1, order_z+1
                do nt = 1, order_t+1
                   nk = nk + 1; TcXl(nk,nm) = TcX(nx,nm) ; TcYl(nk,nm) = TcY(ny,nm)
                   TcZl(nk,nm) = TcZ(nz,nm) ;  TcTl(nk,nm) = TcT(nt,nm)
                end do
             end do
          end do
       end do
       X = X + fd(nm)*(TcXl(:,nm)*TcYl(:,nm)*TcZl(:,nm)*TcTl(:,nm)/real(m)*real(16)) !???? real(8)
    end do

    if (verbose.gt.0) write(verbose,*) 'CHEBY_3D_FIT: Deallocation of temporary variables'
    deallocate(TcX,TcY,TcZ,TcT)

!       X(1) = X(1) + yd(nm)/real(m)
!       X(2) = X(2) + yd(nm)*Tc(1,nm)/real(m)
!       X(3) = X(3) + yd(nm)*Tc(2,nm)/real(m)
!       ! write(*,*) X(1),X(2),X(3)
!       do nk = 2, order_x-1
!          X(nk+2) = X(nk+2) +  yd(nm)*Tc(nk+1,nm)/real(m)
!       end do

    if (verbose.gt.0) write(verbose,*) 'CHEBY_3D_FIT: CALLING lmdif subroutine'

    fvec = 0._PR; diag = 0._PR; sig = 0._PR
    TOL = D1MACH(4)
    INFO=0
    N=(order_x+1)*(order_y+1)*(order_z+1)*(order_t+1)

    call lmdif2(cheby_4d_f,m,n,x,fvec,tol,diag,sig,info,iwa)


!!$    fvec = 0._PR; diag = 0._PR; sig = 0._PR
!!$    TOL = D1MACH(4)
!!$    INFO=0
!!$    N = 4
!!$
!!$    call lmdif2(cheby_4d_f,m,n,x,fvec,tol,diag,sig,info,iwa)


    if (verbose.gt.0) write(verbose,*) 'CHEBY_3D_FIT: Saving the coefficients'

    nk = 0
    do nx = 1, order_x+1
       do ny = 1, order_y+1
          do nz = 1, order_z+1
             do nt = 1, order_t+1
                nk = nk + 1; coeff(nx,ny,nz,nt) = X(nk)
             end do
          end do
       end do
    end do
    if (verbose.gt.0) write(verbose,*) 'CHEBY_3D_FIT: Variable deallocation'

    deallocate(xd,yd,zd,td,fd,wd,TcXl,TcYl,TcZl,TcTl)
    return
  end subroutine cheby_4d_fit

  subroutine  cheby_4d_f(M,N,X,FVEC,IFLAG)
    INTEGER,INTENT(IN)::m,n
    REAL (PR),INTENT(IN)::x(:)
    REAL (PR),INTENT(OUT)::fvec(:)
    INTEGER,INTENT(inout)::iflag
    real (PR) :: c_fit
    integer :: k

    do k=1,M
       c_fit = sum(X(:)*TcXl(:,k)*TcYl(:,k)*TcZl(:,k)*TcTl(:,k))
       FVEC(k)= (fd(k) - c_fit)*wd(k)
    enddo
    return
  end subroutine cheby_4d_f


  real (PR) function get_cheby_4d(x,y,z,t,coeff)
    real (PR), intent(in) :: x,y,z,t
    real (PR), dimension(:,:,:,:), intent(in) :: coeff
    real (PR), dimension(size(coeff,1)) :: Tx
    real (PR), dimension(size(coeff,2)) :: Ty
    real (PR), dimension(size(coeff,3)) :: Tz
    real (PR), dimension(size(coeff,4)) :: Tt
    integer :: n, nx, ny, nz, nt

    get_cheby_4d = 0.00D0

    Tx(1) = 1
    Tx(2) = x
    Ty(1) = 1
    Ty(2) = y
    Tz(1) = 1
    Tz(2) = z
    Tt(1) = 1
    Tt(2) = y

    do n=3,size(coeff,1)
       Tx(n) =  cheb_Tmp(x,Tx(n-1),Tx(n-2))
    end do
    do n=3,size(coeff,2)
       Ty(n) =  cheb_Tmp(y,Ty(n-1),Ty(n-2))
    end do
    do n=3,size(coeff,3)
       Tz(n) =  cheb_Tmp(z,Tz(n-1),Tz(n-2))
    end do
    do n=3,size(coeff,4)
       Tt(n) =  cheb_Tmp(t,Tt(n-1),Tt(n-2))
    end do

    do nx = 1,size(coeff,1)
       do ny = 1, size(coeff,2)
          do nz = 1, size(coeff,3)
             do nt = 1, size(coeff,4)
                get_cheby_4d = get_cheby_4d + coeff(nx,ny,nz,nt)*Tx(nx)*Ty(ny)*Tz(nz)*Tt(nt)
             end do
          end do
       end do
    end do
    return
  end function get_cheby_4d

  !!!!! ********** CHEBY General subroutines ********** !!!!!


  real (PR) function cheb_Tmp(x,T_n,T_nm)
    ! How to use it:
    ! T(1) = x
    ! T(2) = cheb_Tmp(x,T(1),1)
    ! do n=2,order-1
    !   T(n+1) =  cheb_Tmp(x,T(n),T(n-1))
    ! end do
    ! since T_0 is ALWAYS equal to 1, the function will not give that back
    ! instead, there will be the corrispondence T(n)=T_n
    ! T_n = Tn(x) , T_nm = T_(n-1) (x), Tnp = T_(n+1) (x)
    ! T_n and T_nm must be initialized to 1 and x if this is the first run
    real (PR) :: x, T_n, T_nm
    cheb_Tmp = 2.000000D0*x*T_n - T_nm
    return
  end function cheb_Tmp


end module chebyshev_fit
