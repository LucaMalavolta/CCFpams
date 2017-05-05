MODULE Levenberg_Marquardt_mod5V
  use common
  ! MINPACK routines which are used by both LMDIF & LMDER

  IMPLICIT NONE
  !INTEGER,PARAMETER::dp = SELECTED_REAL_KIND(12,60)

  PRIVATE
  !PUBLIC::dp,lmdif1,lmdif,lmder1,lmder,enorm
  PUBLIC::lmdif1_mod5V,lmdif2_mod5V,lmdif_mod5V,lmder1_mod5V,lmder_mod5V,enorm_mod5V

CONTAINS


  SUBROUTINE lmdif1_mod5V(fcn_mod5V,m,n,x,xd1,xd2,xd3,xd4,xd5,fvec,tol,info,iwa)

    ! Code converted using TO_F90 by Alan Miller
    ! Date: 1999-12-11  Time: 00:51:44

    ! N.B. Arguments WA & LWA have been removed.

    INTEGER,INTENT(IN)::m
    INTEGER,INTENT(IN)::n
    REAL (PR),INTENT(INOUT)::x(:)
    REAL (PR),INTENT(OUT)::fvec(:)
    REAL (PR),INTENT(IN)::tol
    INTEGER,INTENT(OUT)::info
    INTEGER,INTENT(OUT)::iwa(:)
    REAL (PR),INTENT(IN)::xd1(:),xd2(:),xd3(:),xd4(:),xd5(:)

    ! EXTERNAL fcn

    INTERFACE
       SUBROUTINE fcn_mod5V(m,n,x,xd1,xd2,xd3,xd4,xd5,fvec,iflag)
         use common
         IMPLICIT NONE
         !INTEGER,PARAMETER::dp = SELECTED_REAL_KIND(12,60)
         INTEGER,INTENT(IN)::m,n
         REAL (PR),INTENT(IN)::x(:)
         REAL (PR),INTENT(OUT)::fvec(:)
         INTEGER,INTENT(inout)::iflag
         REAL (PR),INTENT(IN)::xd1(:),xd2(:),xd3(:),xd4(:),xd5(:)
       END SUBROUTINE fcn_mod5V
    END INTERFACE

    !     **********

    !     subroutine lmdif1_mod5V

    !     the purpose of lmdif1 is to minimize the sum of the squares of
    !     m nonlinear functions in n variables by a modification of the
    !     levenberg-marquardt algorithm. this is done by using the more
    !     general least-squares solver lmdif. the user must provide a
    !     subroutine which calculates the functions. the jacobian is
    !     then calculated by a forward-difference approximation.

    !     the subroutine statement is

    !       subroutine lmdif1(fcn,m,n,x,fvec,tol,info,iwa)

    !     where

    !       fcn is the name of the user-supplied subroutine which
    !         calculates the functions. fcn must be declared
    !         in an external statement in the user calling
    !         program,and should be written as follows.

    !         subroutine fcn(m,n,x,fvec,iflag)
    !         integer m,n,iflag
    !         REAL (PR) x(n),fvec(m)
    !         ----------
    !         calculate the functions at x and
    !         return this vector in fvec.
    !         ----------
    !         return
    !         end

    !         the value of iflag should not be changed by fcn unless
    !         the user wants to terminate execution of lmdif1.
    !         in this case set iflag to a negative integer.

    !       m is a positive integer input variable set to the number of functions.

    !       n is a positive integer input variable set to the number
    !         of variables.  n must not exceed m.

    !       x is an array of length n.  on input x must contain
    !         an initial estimate of the solution vector. on output x
    !         contains the final estimate of the solution vector.

    !       fvec is an output array of length m which contains
    !         the functions evaluated at the output x.

    !       tol is a nonnegative input variable. termination occurs
    !         when the algorithm estimates either that the relative
    !         error in the sum of squares is at most tol or that
    !         the relative error between x and the solution is at most tol.

    !       info is an integer output variable. if the user has
    !         terminated execution,info is set to the (negative)
    !         value of iflag. see description of fcn. otherwise,
    !         info is set as follows.

    !         info = 0  improper input parameters.

    !         info = 1  algorithm estimates that the relative error
    !                   in the sum of squares is at most tol.

    !         info = 2  algorithm estimates that the relative error
    !                   between x and the solution is at most tol.

    !         info = 3  conditions for info = 1 and info = 2 both hold.

    !         info = 4  fvec is orthogonal to the columns of the
    !                   jacobian to machine precision.

    !         info = 5  number of calls to fcn has reached or
    !                   exceeded 200*(n+1).

    !         info = 6  tol is too small. no further reduction in
    !                   the sum of squares is possible.

    !         info = 7  tol is too small. no further improvement in
    !                   the approximate solution x is possible.

    !       iwa is an integer work array of length n.

    !       wa is a work array of length lwa.

    !       lwa is a positive integer input variable not less than
    !         m*n+5*n+m.

    !     subprograms called

    !       user-supplied ...... fcn

    !       minpack-supplied ... lmdif

    !     argonne national laboratory. minpack project. march 1980.
    !     burton s. garbow,kenneth e. hillstrom,jorge j. more

    !     **********
    INTEGER::maxfev,mode,nfev,nprint
    REAL (PR)::epsfcn,ftol,gtol,xtol,wa(2*n),fjac(m,n)
    REAL (PR),PARAMETER::factor = 100._PR,zero = 0.0_PR

    info = 0

    !     check the input parameters for errors.

    IF (n <= 0 .OR. m < n .OR. tol < zero) GO TO 10

    !     call lmdif.

    maxfev = 200*(n + 1)
    ftol = tol
    xtol = tol
    gtol = zero
    epsfcn = zero
    mode = 1
    nprint = 0
    CALL lmdif_mod5V(fcn_mod5V,m,n,x,xd1,xd2,xd3,xd4,xd5,fvec,ftol,xtol,gtol,maxfev,epsfcn,wa,  &
         mode,factor,nprint,info,nfev,fjac,iwa,wa(n+1:))
    IF (info == 8) info = 4

10  RETURN

    !     last card of subroutine lmdif1_mod5V.

  END SUBROUTINE lmdif1_mod5V


  SUBROUTINE lmdif2_mod5V(fcn_mod5V,m,n,x,xd1,xd2,xd3,xd4,xd5,fvec,tol,diag,sig,info,iwa)
    !  use functions,only:newunit
    INTEGER,INTENT(IN)::m
    INTEGER,INTENT(IN)::n
    REAL (PR),INTENT(INOUT)::x(:)
    REAL (PR),INTENT(OUT)::fvec(:),diag(:),sig(:)
    REAL (PR),INTENT(IN)::tol
    INTEGER,INTENT(OUT)::info
    INTEGER,INTENT(OUT)::iwa(:)
    REAL (PR),INTENT(IN)::xd1(:),xd2(:),xd3(:),xd4(:),xd5(:)

    INTERFACE
       SUBROUTINE fcn_mod5V(m,n,x,xd1,xd2,xd3,xd4,xd5,fvec,iflag)
         use common
         IMPLICIT NONE
         !INTEGER,PARAMETER::dp = SELECTED_REAL_KIND(12,60)
         INTEGER,INTENT(IN)::m,n
         REAL (PR),INTENT(IN)::x(:)
         REAL (PR),INTENT(OUT)::fvec(:)
         INTEGER,INTENT(inout)::iflag
         REAL (PR),INTENT(IN)::xd1(:),xd2(:),xd3(:),xd4(:),xd5(:)
       END SUBROUTINE fcn_mod5V
    END INTERFACE

    INTEGER::maxfev,mode,nfev,nprint,uread
    REAL (PR)::epsfcn,ftol,gtol,xtol,wa(2*n),fjac(m,n)
    REAL (PR),PARAMETER::factor = 100._PR,zero = 0.0_PR

    ! added for the covar subroutine
    real(PR),dimension(m,n)::r
    integer,dimension(size(iwa))::ipvt
    real(PR),dimension(m)::wa2
    integer::ldr,j

    info = 0

    !     check the input parameters for errors.
    IF (n <= 0 .OR. m < n .OR. tol < zero) GO TO 10
!!$    !     call lmdif.
!!$    uread=newunit(1)
!!$    open(uread,file="lm.opt",status='OLD')
!!$    read(uread,*)maxfev
!!$    if(maxfev.le.0) maxfev = 200*(n + 1)
!!$    read(uread,*)ftol
!!$    if(ftol.le.0._PR) ftol = tol
!!$    read(uread,*)xtol
!!$    if(xtol.le.0._PR) xtol = tol
!!$    read(uread,*)gtol
!!$    if(gtol.lt.0._PR) gtol = zero
!!$    read(uread,*)epsfcn
!!$    if(epsfcn.lt.0._PR) epsfcn=sqrt((10._PR)**(-precision(tol)))
!!$    write(*,*)" lmdif2_mod5V setted tolerances: "
!!$    write(*,*)" ftol = ",ftol," xtol = ",xtol," gtol = ",gtol," epsfcn = ",epsfcn
!!$    mode = 1
!!$    read(uread,*)nprint
!!$    if(nprint.lt.0) nprint = 0
!!$    close(uread)

    maxfev = 200*(n + 1)
    ftol = tol
    xtol = tol
    gtol = zero
    !epsfcn=sqrt((10._PR)**(-precision(tol)))
    epsfcn=epsilon(1._PR)
!!$ write(*,*)" lmdif2_mod5V setted tolerances: "
!!$ write(*,*)" ftol = ",ftol," xtol = ",xtol," gtol = ",gtol," epsfcn = ",epsfcn
    mode = 1
    nprint = 0

!!$ write(*,*)" nprint = ",nprint
    CALL lmdif_mod5V(fcn_mod5V,m,n,x,xd1,xd2,xd3,xd4,xd5,fvec,ftol,xtol,gtol,maxfev,epsfcn,wa,&
         mode,factor,nprint,info,nfev,fjac,iwa,wa(n+1:))
    if(info == 8)then
!!$ write(*,*)"** info = 8 will be set to 4 **"
       info = 4
    end if
!!$ write(*,*)""
!!$ write(*,*)" - - - - - - - - - - - - - - - - - "
!!$ write(*,*)""
    ! add by Luca Borsato
    r=fjac
    ipvt=iwa
    ldr=m
    call covar(n,r,ldr,ipvt,tol,wa2)
    do j=1,n
       diag(j)=r(j,j)
       sig(j)=sqrt(diag(j))
    end do

10  RETURN

    !     last card of subroutine lmdif2_mod5V

  END SUBROUTINE lmdif2_mod5V




  SUBROUTINE lmder1_mod5V(fcn_mod5V,m,n,x,xd1,xd2,xd3,xd4,xd5,fvec,fjac,tol,info,ipvt)

    ! Code converted using TO_F90 by Alan Miller
    ! Date: 1999-12-09  Time: 12:45:54

    ! N.B. Arguments LDFJAC,WA & LWA have been removed.

    INTEGER,INTENT(IN)::m
    INTEGER,INTENT(IN)::n
    REAL (PR),INTENT(INOUT)::x(:)
    REAL (PR),INTENT(INOUT)::fvec(:)
    REAL (PR),INTENT(INOUT)::fjac(:,:)    ! fjac(ldfjac,n)
    REAL (PR),INTENT(IN)::tol
    INTEGER,INTENT(OUT)::info
    INTEGER,INTENT(INOUT)::ipvt(:)
    REAL (PR),INTENT(IN)::xd1(:),xd2(:),xd3(:),xd4(:),xd5(:)


    ! EXTERNAL fcn

    INTERFACE
       SUBROUTINE fcn_mod5V(m,n,x,xd1,xd2,xd3,xd4,xd5,fvec,fjac,iflag)
         use common
         IMPLICIT NONE
         !INTEGER,PARAMETER::PR = SELECTED_REAL_KIND(12,60)
         INTEGER,INTENT(IN)::m,n
         REAL (PR),INTENT(IN)::x(:)
         REAL (PR),INTENT(OUT)::fvec(:)
         REAL (PR),INTENT(OUT)::fjac(:,:)
         INTEGER,INTENT(inout)::iflag
         REAL (PR),INTENT(IN)::xd1(:),xd2(:),xd3(:),xd4(:),xd5(:)
       END SUBROUTINE fcn_mod5V
    END INTERFACE

    !     **********

    !     subroutine lmder1_mod5V

    !     The purpose of lmder1_mod5V is to minimize the sum of the squares of
    !     m nonlinear functions in n variables by a modification of the
    !     levenberg-marquardt algorithm.  This is done by using the more
    !     general least-squares solver lmder.  The user must provide a
    !     subroutine which calculates the functions and the jacobian.

    !     the subroutine statement is

    !       subroutine lmder1_mod5V(fcn,m,n,x,fvec,fjac,tol,info,ipvt)

    !     where

    !       fcn is the name of the user-supplied subroutine which
    !         calculates the functions and the jacobian.  fcn must
    !         be declared in an interface statement in the user
    !         calling program,and should be written as follows.

    !         subroutine fcn(m,n,x,fvec,fjac,iflag)
    !         integer::m,n,ldfjac,iflag
    !         REAL (PR)::x(:),fvec(:),fjac(:,:)
    !         ----------
    !         if iflag = 1 calculate the functions at x and
    !         return this vector in fvec. do not alter fjac.
    !         if iflag = 2 calculate the jacobian at x and
    !         return this matrix in fjac. do not alter fvec.
    !         ----------
    !         return
    !         end

    !         the value of iflag should not be changed by fcn unless
    !         the user wants to terminate execution of lmder1_mod5V.
    !         in this case set iflag to a negative integer.

    !       m is a positive integer input variable set to the number of functions.

    !       n is a positive integer input variable set to the number
    !         of variables.  n must not exceed m.

    !       x is an array of length n. on input x must contain
    !         an initial estimate of the solution vector. on output x
    !         contains the final estimate of the solution vector.

    !       fvec is an output array of length m which contains
    !         the functions evaluated at the output x.

    !       fjac is an output m by n array. the upper n by n submatrix
    !         of fjac contains an upper triangular matrix r with
    !         diagonal elements of nonincreasing magnitude such that

    !                t     t           t
    !               p *(jac *jac)*p = r *r,

    !         where p is a permutation matrix and jac is the final calculated
    !         Jacobian.  Column j of p is column ipvt(j) (see below) of the
    !         identity matrix.  The lower trapezoidal part of fjac contains
    !         information generated during the computation of r.

    !       ldfjac is a positive integer input variable not less than m
    !         which specifies the leading dimension of the array fjac.

    !       tol is a nonnegative input variable. termination occurs
    !         when the algorithm estimates either that the relative
    !         error in the sum of squares is at most tol or that
    !         the relative error between x and the solution is at most tol.

    !       info is an integer output variable.  If the user has terminated
    !         execution,info is set to the (negative) value of iflag.
    !         See description of fcn.  Otherwise,info is set as follows.

    !         info = 0  improper input parameters.

    !         info = 1  algorithm estimates that the relative error
    !                   in the sum of squares is at most tol.

    !         info = 2  algorithm estimates that the relative error
    !                   between x and the solution is at most tol.

    !         info = 3  conditions for info = 1 and info = 2 both hold.

    !         info = 4  fvec is orthogonal to the columns of the
    !                   jacobian to machine precision.

    !         info = 5  number of calls to fcn with iflag = 1 has reached 100*(n+1).

    !         info = 6  tol is too small.  No further reduction in
    !                   the sum of squares is possible.

    !         info = 7  tol is too small.  No further improvement in
    !                   the approximate solution x is possible.

    !       ipvt is an integer output array of length n. ipvt
    !         defines a permutation matrix p such that jac*p = q*r,
    !         where jac is the final calculated jacobian,q is
    !         orthogonal (not stored),and r is upper triangular
    !         with diagonal elements of nonincreasing magnitude.
    !         column j of p is column ipvt(j) of the identity matrix.

    !       wa is a work array of length lwa.

    !       lwa is a positive integer input variable not less than 5*n+m.

    !     subprograms called

    !       user-supplied ...... fcn

    !       minpack-supplied ... lmder

    !     argonne national laboratory. minpack project. march 1980.
    !     burton s. garbow,kenneth e. hillstrom,jorge j. more

    !     **********
    INTEGER::maxfev,mode,nfev,njev,nprint
    REAL (PR)::ftol,gtol,xtol,wa(2*n)
    REAL (PR),PARAMETER::factor = 100._PR,zero = 0.0_PR

    info = 0

    !     check the input parameters for errors.

    IF ( n <= 0 .OR. m < n .OR. tol < zero ) GO TO 10

    !     call lmder_mod5V.

    maxfev = 100*(n + 1)
    ftol = tol
    xtol = tol
    gtol = zero
    mode = 1
    nprint = 0
    CALL lmder_mod5V(fcn_mod5V,m,n,x,xd1,xd2,xd3,xd4,xd5,fvec,fjac,ftol,xtol,gtol,maxfev, &
         wa,mode,factor,nprint,info,nfev,njev,ipvt,wa(n+1:) )
    IF (info == 8) info = 4

10  RETURN

    !     last card of subroutine lmder1_mod5V.

  END SUBROUTINE lmder1_mod5V


  SUBROUTINE lmder_mod5V(fcn_mod5V,m,n,x,xd1,xd2,xd3,xd4,xd5,fvec,fjac,ftol,xtol,gtol,maxfev,&
       diag,mode,factor,nprint,info,nfev,njev,ipvt,qtf)

    ! Code converted using TO_F90 by Alan Miller
    ! Date: 1999-12-09  Time: 12:45:50

    ! N.B. Arguments LDFJAC,WA1,WA2,WA3 & WA4 have been removed.

    INTEGER,INTENT(IN)::m
    INTEGER,INTENT(IN)::n
    REAL (PR),INTENT(INOUT)::x(:)
    REAL (PR),INTENT(OUT)::fvec(m)
    REAL (PR),INTENT(OUT)::fjac(:,:)    ! fjac(ldfjac,n)
    REAL (PR),INTENT(IN)::ftol
    REAL (PR),INTENT(IN)::xtol
    REAL (PR),INTENT(INOUT)::gtol
    INTEGER,INTENT(INOUT)::maxfev
    REAL (PR),INTENT(OUT)::diag(:)
    INTEGER,INTENT(IN)::mode
    REAL (PR),INTENT(IN)::factor
    INTEGER,INTENT(IN)::nprint
    INTEGER,INTENT(OUT)::info
    INTEGER,INTENT(OUT)::nfev
    INTEGER,INTENT(OUT)::njev
    INTEGER,INTENT(OUT)::ipvt(:)
    REAL (PR),INTENT(OUT)::qtf(:)
    REAL (PR),INTENT(IN)::xd1(:),xd2(:),xd3(:),xd4(:),xd5(:)

    INTERFACE
       SUBROUTINE fcn_mod5V(m,n,x,xd1,xd2,xd3,xd4,xd5,fvec,fjac,iflag)
         use common
         IMPLICIT NONE
         !INTEGER,PARAMETER::PR = SELECTED_REAL_KIND(12,60)
         INTEGER,INTENT(IN)::m,n
         REAL (PR),INTENT(IN)::x(:)
         REAL (PR),INTENT(OUT)::fvec(:)
         REAL (PR),INTENT(OUT)::fjac(:,:)
         INTEGER,INTENT(inout)::iflag
         REAL (PR),INTENT(IN)::xd1(:),xd2(:),xd3(:),xd4(:),xd5(:)
       END SUBROUTINE fcn_mod5V
    END INTERFACE


    !     **********

    !     subroutine lmder

    !     the purpose of lmder is to minimize the sum of the squares of
    !     m nonlinear functions in n variables by a modification of
    !     the levenberg-marquardt algorithm. the user must provide a
    !     subroutine which calculates the functions and the jacobian.

    !     the subroutine statement is

    !       subroutine lmder(fcn,m,n,x,fvec,fjac,ldfjac,ftol,xtol,gtol,
    !                        maxfev,diag,mode,factor,nprint,info,nfev,
    !                        njev,ipvt,qtf,wa1,wa2,wa3,wa4)

    !     where

    !       fcn is the name of the user-supplied subroutine which
    !         calculates the functions and the jacobian. fcn must
    !         be declared in an external statement in the user
    !         calling program,and should be written as follows.

    !         subroutine fcn(m,n,x,fvec,fjac,ldfjac,iflag)
    !         integer m,n,ldfjac,iflag
    !         REAL (PR) x(:),fvec(m),fjac(ldfjac,n)
    !         ----------
    !         if iflag = 1 calculate the functions at x and
    !         return this vector in fvec. do not alter fjac.
    !         if iflag = 2 calculate the jacobian at x and
    !         return this matrix in fjac. do not alter fvec.
    !         ----------
    !         return
    !         end

    !         the value of iflag should not be changed by fcn unless
    !         the user wants to terminate execution of lmder.
    !         in this case set iflag to a negative integer.

    !       m is a positive integer input variable set to the number
    !         of functions.

    !       n is a positive integer input variable set to the number
    !         of variables. n must not exceed m.

    !       x is an array of length n. on input x must contain
    !         an initial estimate of the solution vector. on output x
    !         contains the final estimate of the solution vector.

    !       fvec is an output array of length m which contains
    !         the functions evaluated at the output x.

    !       fjac is an output m by n array. the upper n by n submatrix
    !         of fjac contains an upper triangular matrix r with
    !         diagonal elements of nonincreasing magnitude such that

    !                t     t           t
    !               p *(jac *jac)*p = r *r

    !         where p is a permutation matrix and jac is the final calculated
    !         jacobian.  Column j of p is column ipvt(j) (see below) of the
    !         identity matrix.  The lower trapezoidal part of fjac contains
    !         information generated during the computation of r.

    !       ldfjac is a positive integer input variable not less than m
    !         which specifies the leading dimension of the array fjac.

    !       ftol is a nonnegative input variable.  Termination occurs when both
    !         the actual and predicted relative reductions in the sum of squares
    !         are at most ftol.   Therefore,ftol measures the relative error
    !         desired in the sum of squares.

    !       xtol is a nonnegative input variable. termination
    !         occurs when the relative error between two consecutive
    !         iterates is at most xtol. therefore,xtol measures the
    !         relative error desired in the approximate solution.

    !       gtol is a nonnegative input variable.  Termination occurs when the
    !         cosine of the angle between fvec and any column of the jacobian is
    !         at most gtol in absolute value.  Therefore,gtol measures the
    !         orthogonality desired between the function vector and the columns
    !         of the jacobian.

    !       maxfev is a positive integer input variable.  Termination occurs when
    !         the number of calls to fcn with iflag = 1 has reached maxfev.

    !       diag is an array of length n.  If mode = 1 (see below),diag is
    !         internally set.  If mode = 2,diag must contain positive entries
    !         that serve as multiplicative scale factors for the variables.

    !       mode is an integer input variable.  if mode = 1,the
    !         variables will be scaled internally.  if mode = 2,
    !         the scaling is specified by the input diag.  other
    !         values of mode are equivalent to mode = 1.

    !       factor is a positive input variable used in determining the
    !         initial step bound. this bound is set to the product of
    !         factor and the euclidean norm of diag*x if nonzero,or else
    !         to factor itself. in most cases factor should lie in the
    !         interval (.1,100.).100. is a generally recommended value.

    !       nprint is an integer input variable that enables controlled printing
    !         of iterates if it is positive.  In this case,fcn is called with
    !         iflag = 0 at the beginning of the first iteration and every nprint
    !         iterations thereafter and immediately prior to return,with x,fvec,
    !         and fjac available for printing.  fvec and fjac should not be
    !         altered.  If nprint is not positive,no special calls of fcn with
    !         iflag = 0 are made.

    !       info is an integer output variable.  If the user has terminated
    !         execution,info is set to the (negative) value of iflag.
    !         See description of fcn.  Otherwise,info is set as follows.

    !         info = 0  improper input parameters.

    !         info = 1  both actual and predicted relative reductions
    !                   in the sum of squares are at most ftol.

    !         info = 2  relative error between two consecutive iterates
    !                   is at most xtol.

    !         info = 3  conditions for info = 1 and info = 2 both hold.

    !         info = 4  the cosine of the angle between fvec and any column of
    !                   the jacobian is at most gtol in absolute value.

    !         info = 5  number of calls to fcn with iflag = 1 has reached maxfev.

    !         info = 6  ftol is too small.  No further reduction in
    !                   the sum of squares is possible.

    !         info = 7  xtol is too small.  No further improvement in
    !                   the approximate solution x is possible.

    !         info = 8  gtol is too small.  fvec is orthogonal to the
    !                   columns of the jacobian to machine precision.

    !       nfev is an integer output variable set to the number of
    !         calls to fcn with iflag = 1.

    !       njev is an integer output variable set to the number of
    !         calls to fcn with iflag = 2.

    !       ipvt is an integer output array of length n.  ipvt
    !         defines a permutation matrix p such that jac*p = q*r,
    !         where jac is the final calculated jacobian,q is
    !         orthogonal (not stored),and r is upper triangular
    !         with diagonal elements of nonincreasing magnitude.
    !         column j of p is column ipvt(j) of the identity matrix.

    !       qtf is an output array of length n which contains
    !         the first n elements of the vector (q transpose)*fvec.

    !       wa1,wa2,and wa3 are work arrays of length n.

    !       wa4 is a work array of length m.

    !     subprograms called

    !       user-supplied ...... fcn

    !       minpack-supplied ... PRmpar,enorm,lmpar,qrfac

    !       fortran-supplied ... ABS,MAX,MIN,SQRT,mod

    !     argonne national laboratory. minpack project. march 1980.
    !     burton s. garbow,kenneth e. hillstrom,jorge j. more

    !     **********
    INTEGER::i,iflag,iter,j,l
    REAL (PR)::actred,delta,dirder,epsmch,fnorm,fnorm1,gnorm, &
         par,pnorm,prered,ratio,sum,temp,temp1,temp2,xnorm
    REAL (PR)::wa1(n),wa2(n),wa3(n),wa4(m)
    REAL (PR),PARAMETER::one = 1.0_PR,p1 = 0.1_PR,p5 = 0.5_PR, &
         p25 = 0.25_PR,p75 = 0.75_PR,p0001 = 0.0001_PR,&
         zero = 0.0_PR

    !     epsmch is the machine precision.

    epsmch = EPSILON(zero)

    info = 0
    iflag = 0
    nfev = 0
    njev = 0

    !     check the input parameters for errors.

    IF (n <= 0 .OR. m < n .OR. ftol < zero .OR. xtol < zero .OR. gtol < zero  &
         .OR. maxfev <= 0 .OR. factor <= zero) GO TO 300
    IF (mode /= 2) GO TO 20
    DO  j = 1,n
       IF (diag(j) <= zero) GO TO 300
    END DO

    !     evaluate the function at the starting point and calculate its norm.

20  iflag = 1
    CALL fcn_mod5V(m,n,x,xd1,xd2,xd3,xd4,xd5,fvec,fjac,iflag)
    nfev = 1
    IF (iflag < 0) GO TO 300
    fnorm = enorm_mod5V(m,fvec)

    !     initialize levenberg-marquardt parameter and iteration counter.

    par = zero
    iter = 1

    !     beginning of the outer loop.

    !        calculate the jacobian matrix.

30  iflag = 2
    CALL fcn_mod5V(m,n,x,xd1,xd2,xd3,xd4,xd5,fvec,fjac,iflag)
    njev = njev + 1
    IF (iflag < 0) GO TO 300

    !        if requested,call fcn to enable printing of iterates.

    IF (nprint <= 0) GO TO 40
    iflag = 0
    IF (MOD(iter-1,nprint) == 0) CALL fcn_mod5V(m,n,x,xd1,xd2,xd3,xd4,xd5,fvec,fjac,iflag)
    IF (iflag < 0) GO TO 300

    !        compute the qr factorization of the jacobian.

40  CALL qrfac(m,n,fjac,.true.,ipvt,wa1,wa2)

    !        on the first iteration and if mode is 1,scale according
    !        to the norms of the columns of the initial jacobian.

    IF (iter /= 1) GO TO 80
    IF (mode == 2) GO TO 60
    DO  j = 1,n
       diag(j) = wa2(j)
       IF (wa2(j) == zero) diag(j) = one
    END DO

    !        on the first iteration,calculate the norm of the scaled x
    !        and initialize the step bound delta.

60  wa3(1:n) = diag(1:n)*x(1:n)
    xnorm = enorm_mod5V(n,wa3)
    delta = factor*xnorm
    IF (delta == zero) delta = factor

    !        form (q transpose)*fvec and store the first n components in qtf.

80  wa4(1:m) = fvec(1:m)
    DO  j = 1,n
       IF (fjac(j,j) == zero) GO TO 120
       sum = DOT_PRODUCT( fjac(j:m,j),wa4(j:m) )
       temp = -sum/fjac(j,j)
       DO  i = j,m
          wa4(i) = wa4(i) + fjac(i,j)*temp
       END DO
120    fjac(j,j) = wa1(j)
       qtf(j) = wa4(j)
    END DO

    !        compute the norm of the scaled gradient.

    gnorm = zero
    IF (fnorm == zero) GO TO 170
    DO  j = 1,n
       l = ipvt(j)
       IF (wa2(l) == zero) CYCLE
       sum = zero
       DO  i = 1,j
          sum = sum + fjac(i,j)*(qtf(i)/fnorm)
       END DO
       gnorm = MAX(gnorm,ABS(sum/wa2(l)))
    END DO

    !        test for convergence of the gradient norm.

170 IF (gnorm <= gtol) info = 4
    IF (info /= 0) GO TO 300

    !        rescale if necessary.

    IF (mode == 2) GO TO 200
    DO  j = 1,n
       diag(j) = MAX(diag(j),wa2(j))
    END DO

    !        beginning of the inner loop.

    !           determine the levenberg-marquardt parameter.

200 CALL lmpar(n,fjac,ipvt,diag,qtf,delta,par,wa1,wa2)

    !           store the direction p and x + p. calculate the norm of p.

    DO  j = 1,n
       wa1(j) = -wa1(j)
       wa2(j) = x(j) + wa1(j)
       wa3(j) = diag(j)*wa1(j)
    END DO
    pnorm = enorm_mod5V(n,wa3)

    !           on the first iteration,adjust the initial step bound.

    IF (iter == 1) delta = MIN(delta,pnorm)

    !           evaluate the function at x + p and calculate its norm.

    iflag = 1
    CALL fcn_mod5V(m,n,wa2,xd1,xd2,xd3,xd4,xd5,wa4,fjac,iflag)
    nfev = nfev + 1
    IF (iflag < 0) GO TO 300
    fnorm1 = enorm_mod5V(m,wa4)

    !           compute the scaled actual reduction.

    actred = -one
    IF (p1*fnorm1 < fnorm) actred = one - (fnorm1/fnorm)**2

    !           compute the scaled predicted reduction and
    !           the scaled directional derivative.

    DO  j = 1,n
       wa3(j) = zero
       l = ipvt(j)
       temp = wa1(l)
       DO  i = 1,j
          wa3(i) = wa3(i) + fjac(i,j)*temp
       END DO
    END DO
    temp1 = enorm_mod5V(n,wa3)/fnorm
    temp2 = (SQRT(par)*pnorm)/fnorm
    prered = temp1**2 + temp2**2/p5
    dirder = -(temp1**2 + temp2**2)

    !           compute the ratio of the actual to the predicted reduction.

    ratio = zero
    IF (prered /= zero) ratio = actred/prered

    !           update the step bound.

    IF (ratio <= p25) THEN
       IF (actred >= zero) temp = p5
       IF (actred < zero) temp = p5*dirder/(dirder + p5*actred)
       IF (p1*fnorm1 >= fnorm .OR. temp < p1) temp = p1
       delta = temp*MIN(delta,pnorm/p1)
       par = par/temp
    ELSE
       IF (par /= zero .AND. ratio < p75) GO TO 260
       delta = pnorm/p5
       par = p5*par
    END IF

    !           test for successful iteration.

260 IF (ratio < p0001) GO TO 290

    !           successful iteration. update x,fvec,and their norms.

    DO  j = 1,n
       x(j) = wa2(j)
       wa2(j) = diag(j)*x(j)
    END DO
    fvec(1:m) = wa4(1:m)
    xnorm = enorm_mod5V(n,wa2)
    fnorm = fnorm1
    iter = iter + 1

    !           tests for convergence.

290 IF (ABS(actred) <= ftol .AND. prered <= ftol .AND. p5*ratio <= one) info = 1
    IF (delta <= xtol*xnorm) info = 2
    IF (ABS(actred) <= ftol .AND. prered <= ftol  &
         .AND. p5*ratio <= one .AND. info == 2) info = 3
    IF (info /= 0) GO TO 300

    !           tests for termination and stringent tolerances.

    IF (nfev >= maxfev) info = 5
    IF (ABS(actred) <= epsmch .AND. prered <= epsmch  &
         .AND. p5*ratio <= one) info = 6
    IF (delta <= epsmch*xnorm) info = 7
    IF (gnorm <= epsmch) info = 8
    IF (info /= 0) GO TO 300

    !           end of the inner loop. repeat if iteration unsuccessful.

    IF (ratio < p0001) GO TO 200

    !        end of the outer loop.

    GO TO 30

    !     termination,either normal or user imposed.

300 IF (iflag < 0) info = iflag
    iflag = 0
    IF (nprint > 0) CALL fcn_mod5V(m,n,x,xd1,xd2,xd3,xd4,xd5,fvec,fjac,iflag)
    RETURN

    !     last card of subroutine lmder.

  END SUBROUTINE lmder_mod5V


  SUBROUTINE lmpar(n,r,ipvt,diag,qtb,delta,par,x,sdiag)

    ! Code converted using TO_F90 by Alan Miller
    ! Date: 1999-12-09  Time: 12:46:12

    ! N.B. Arguments LDR,WA1 & WA2 have been removed.

    INTEGER,INTENT(IN)::n
    REAL (PR),INTENT(INOUT)::r(:,:)
    INTEGER,INTENT(IN)::ipvt(:)
    REAL (PR),INTENT(IN)::diag(:)
    REAL (PR),INTENT(IN)::qtb(:)
    REAL (PR),INTENT(IN)::delta
    REAL (PR),INTENT(OUT)::par
    REAL (PR),INTENT(OUT)::x(:)
    REAL (PR),INTENT(OUT)::sdiag(:)


    !     **********

    !     subroutine lmpar

    !     given an m by n matrix a,an n by n nonsingular diagonal
    !     matrix d,an m-vector b,and a positive number delta,
    !     the problem is to determine a value for the parameter
    !     par such that if x solves the system

    !           a*x = b ,    sqrt(par)*d*x = 0 ,

    !     in the least squares sense,and dxnorm is the euclidean
    !     norm of d*x,then either par is zero and

    !           (dxnorm-delta) <= 0.1*delta ,

    !     or par is positive and

    !           abs(dxnorm-delta) <= 0.1*delta .

    !     this subroutine completes the solution of the problem
    !     if it is provided with the necessary information from the
    !     qr factorization,with column pivoting,of a. that is,if
    !     a*p = q*r,where p is a permutation matrix,q has orthogonal
    !     columns,and r is an upper triangular matrix with diagonal
    !     elements of nonincreasing magnitude,then lmpar expects
    !     the full upper triangle of r,the permutation matrix p,
    !     and the first n components of (q transpose)*b. on output
    !     lmpar also provides an upper triangular matrix s such that

    !            t   t                   t
    !           p *(a *a + par*d*d)*p = s *s .

    !     s is employed within lmpar and may be of separate interest.

    !     only a few iterations are generally needed for convergence
    !     of the algorithm. if,however,the limit of 10 iterations
    !     is reached,then the output par will contain the best
    !     value obtained so far.

    !     the subroutine statement is

    !       subroutine lmpar(n,r,ldr,ipvt,diag,qtb,delta,par,x,sdiag,wa1,wa2)

    !     where

    !       n is a positive integer input variable set to the order of r.

    !       r is an n by n array. on input the full upper triangle
    !         must contain the full upper triangle of the matrix r.
    !         on output the full upper triangle is unaltered,and the
    !         strict lower triangle contains the strict upper triangle
    !         (transposed) of the upper triangular matrix s.

    !       ldr is a positive integer input variable not less than n
    !         which specifies the leading dimension of the array r.

    !       ipvt is an integer input array of length n which defines the
    !         permutation matrix p such that a*p = q*r. column j of p
    !         is column ipvt(j) of the identity matrix.

    !       diag is an input array of length n which must contain the
    !         diagonal elements of the matrix d.

    !       qtb is an input array of length n which must contain the first
    !         n elements of the vector (q transpose)*b.

    !       delta is a positive input variable which specifies an upper
    !         bound on the euclidean norm of d*x.

    !       par is a nonnegative variable. on input par contains an
    !         initial estimate of the levenberg-marquardt parameter.
    !         on output par contains the final estimate.

    !       x is an output array of length n which contains the least
    !         squares solution of the system a*x = b,sqrt(par)*d*x = 0,
    !         for the output par.

    !       sdiag is an output array of length n which contains the
    !         diagonal elements of the upper triangular matrix s.

    !       wa1 and wa2 are work arrays of length n.

    !     subprograms called

    !       minpack-supplied ... PRmpar,enorm_mod5V,qrsolv

    !       fortran-supplied ... ABS,MAX,MIN,SQRT

    !     argonne national laboratory. minpack project. march 1980.
    !     burton s. garbow,kenneth e. hillstrom,jorge j. more

    !     **********
    INTEGER::i,iter,j,jm1,jp1,k,l,nsing
    REAL (PR)::dxnorm,dwarf,fp,gnorm,parc,parl,paru,sum,temp
    REAL (PR)::wa1(n),wa2(n)
    REAL (PR),PARAMETER::p1 = 0.1_PR,p001 = 0.001_PR,zero = 0.0_PR

    !     dwarf is the smallest positive magnitude.

    dwarf = TINY(zero)

    !     compute and store in x the gauss-newton direction. if the
    !     jacobian is rank-deficient,obtain a least squares solution.

    nsing = n
    DO  j = 1,n
       wa1(j) = qtb(j)
       IF (r(j,j) == zero .AND. nsing == n) nsing = j - 1
       IF (nsing < n) wa1(j) = zero
    END DO

    DO  k = 1,nsing
       j = nsing - k + 1
       wa1(j) = wa1(j)/r(j,j)
       temp = wa1(j)
       jm1 = j - 1
       DO  i = 1,jm1
          wa1(i) = wa1(i) - r(i,j)*temp
       END DO
    END DO

    DO  j = 1,n
       l = ipvt(j)
       x(l) = wa1(j)
    END DO

    !     initialize the iteration counter.
    !     evaluate the function at the origin,and test
    !     for acceptance of the gauss-newton direction.

    iter = 0
    DO  j = 1,n
       wa2(j) = diag(j)*x(j)
    END DO
    dxnorm = enorm_mod5V(n,wa2)
    fp = dxnorm - delta
    IF (fp <= p1*delta) GO TO 220

    !     if the jacobian is not rank deficient,the newton
    !     step provides a lower bound,parl,for the zero of
    !     the function.  Otherwise set this bound to zero.

    parl = zero
    IF (nsing < n) GO TO 120
    DO  j = 1,n
       l = ipvt(j)
       wa1(j) = diag(l)*(wa2(l)/dxnorm)
    END DO
    DO  j = 1,n
       sum = zero
       jm1 = j - 1
       DO  i = 1,jm1
          sum = sum + r(i,j)*wa1(i)
       END DO
       wa1(j) = (wa1(j) - sum)/r(j,j)
    END DO
    temp = enorm_mod5V(n,wa1)
    parl = ((fp/delta)/temp)/temp

    !     calculate an upper bound,paru,for the zero of the function.

120 DO  j = 1,n
       sum = zero
       DO  i = 1,j
          sum = sum + r(i,j)*qtb(i)
       END DO
       l = ipvt(j)
       wa1(j) = sum/diag(l)
    END DO
    gnorm = enorm_mod5V(n,wa1)
    paru = gnorm/delta
    IF (paru == zero) paru = dwarf/MIN(delta,p1)

    !     if the input par lies outside of the interval (parl,paru),
    !     set par to the closer enPRoint.

    par = MAX(par,parl)
    par = MIN(par,paru)
    IF (par == zero) par = gnorm/dxnorm

    !     beginning of an iteration.

150 iter = iter + 1

    !        evaluate the function at the current value of par.

    IF (par == zero) par = MAX(dwarf,p001*paru)
    temp = SQRT(par)
    wa1(1:n) = temp*diag(1:n)
    CALL qrsolv(n,r,ipvt,wa1,qtb,x,sdiag)
    wa2(1:n) = diag(1:n)*x(1:n)
    dxnorm = enorm_mod5V(n,wa2)
    temp = fp
    fp = dxnorm - delta

    !        if the function is small enough,accept the current value
    !        of par. also test for the exceptional cases where parl
    !        is zero or the number of iterations has reached 10.

    IF (ABS(fp) <= p1*delta .OR. parl == zero .AND. fp <= temp  &
         .AND. temp < zero .OR. iter == 10) GO TO 220

    !        compute the newton correction.

    DO  j = 1,n
       l = ipvt(j)
       wa1(j) = diag(l)*(wa2(l)/dxnorm)
    END DO
    DO  j = 1,n
       wa1(j) = wa1(j)/sdiag(j)
       temp = wa1(j)
       jp1 = j + 1
       DO  i = jp1,n
          wa1(i) = wa1(i) - r(i,j)*temp
       END DO
    END DO
    temp = enorm_mod5V(n,wa1)
    parc = ((fp/delta)/temp)/temp

    !        depending on the sign of the function,update parl or paru.

    IF (fp > zero) parl = MAX(parl,par)
    IF (fp < zero) paru = MIN(paru,par)

    !        compute an improved estimate for par.

    par = MAX(parl,par+parc)

    !        end of an iteration.

    GO TO 150

    !     termination.

220 IF (iter == 0) par = zero
    RETURN

    !     last card of subroutine lmpar.

  END SUBROUTINE lmpar



  SUBROUTINE qrfac(m,n,a,pivot,ipvt,rdiag,acnorm)

    ! Code converted using TO_F90 by Alan Miller
    ! Date: 1999-12-09  Time: 12:46:17

    ! N.B. Arguments LDA,LIPVT & WA has been removed.

    INTEGER,INTENT(IN)::m
    INTEGER,INTENT(IN)::n
    REAL (PR),INTENT(INOUT)::a(:,:)
    LOGICAL,INTENT(IN)::pivot
    INTEGER,INTENT(OUT)::ipvt(:)
    REAL (PR),INTENT(OUT)::rdiag(:)
    REAL (PR),INTENT(OUT)::acnorm(:)


    !     **********

    !     subroutine qrfac

    !     this subroutine uses householder transformations with column
    !     pivoting (optional) to compute a qr factorization of the
    !     m by n matrix a. that is,qrfac determines an orthogonal
    !     matrix q,a permutation matrix p,and an upper trapezoidal
    !     matrix r with diagonal elements of nonincreasing magnitude,
    !     such that a*p = q*r. the householder transformation for
    !     column k,k = 1,2,...,min(m,n),is of the form

    !                           t
    !           i - (1/u(k))*u*u

    !     where u has zeros in the first k-1 positions. the form of
    !     this transformation and the method of pivoting first
    !     appeared in the corresponding linpack subroutine.

    !     the subroutine statement is

    !       subroutine qrfac(m,n,a,lda,pivot,ipvt,lipvt,rdiag,acnorm,wa)

    !     where

    !       m is a positive integer input variable set to the number of rows of a.

    !       n is a positive integer input variable set to the number
    !         of columns of a.

    !       a is an m by n array. on input a contains the matrix for
    !         which the qr factorization is to be computed.  on output
    !         the strict upper trapezoidal part of a contains the strict
    !         upper trapezoidal part of r,and the lower trapezoidal
    !         part of a contains a factored form of q (the non-trivial
    !         elements of the u vectors described above).

    !       lda is a positive integer input variable not less than m
    !         which specifies the leading dimension of the array a.

    !       pivot is a logical input variable.  if pivot is set true,
    !         then column pivoting is enforced.  if pivot is set false,
    !         then no column pivoting is done.

    !       ipvt is an integer output array of length lipvt.  ipvt
    !         defines the permutation matrix p such that a*p = q*r.
    !         column j of p is column ipvt(j) of the identity matrix.
    !         if pivot is false,ipvt is not referenced.

    !       lipvt is a positive integer input variable. if pivot is false,
    !         then lipvt may be as small as 1. if pivot is true,then
    !         lipvt must be at least n.

    !       rdiag is an output array of length n which contains the
    !         diagonal elements of r.

    !       acnorm is an output array of length n which contains the
    !         norms of the corresponding columns of the input matrix a.
    !         If this information is not needed,then acnorm can coincide
    !         with rdiag.

    !       wa is a work array of length n. if pivot is false,then wa
    !         can coincide with rdiag.

    !     subprograms called

    !       minpack-supplied ... PRmpar,enorm_mod5V

    !       fortran-supplied ... MAX,SQRT,MIN

    !     argonne national laboratory. minpack project. march 1980.
    !     burton s. garbow,kenneth e. hillstrom,jorge j. more

    !     **********
    INTEGER::i,j,jp1,k,kmax,minmn
    REAL (PR)::ajnorm,epsmch,sum,temp,wa(n)
    REAL (PR),PARAMETER::one = 1.0_PR,p05 = 0.05_PR,zero = 0.0_PR

    !     epsmch is the machine precision.

    epsmch = EPSILON(zero)

    !     compute the initial column norms and initialize several arrays.

    DO  j = 1,n
       acnorm(j) = enorm_mod5V(m,a(1:,j))
       rdiag(j) = acnorm(j)
       wa(j) = rdiag(j)
       IF (pivot) ipvt(j) = j
    END DO

    !     reduce a to r with householder transformations.

    minmn = MIN(m,n)
    DO  j = 1,minmn
       IF (.NOT.pivot) GO TO 40

       !        bring the column of largest norm into the pivot position.

       kmax = j
       DO  k = j,n
          IF (rdiag(k) > rdiag(kmax)) kmax = k
       END DO
       IF (kmax == j) GO TO 40
       DO  i = 1,m
          temp = a(i,j)
          a(i,j) = a(i,kmax)
          a(i,kmax) = temp
       END DO
       rdiag(kmax) = rdiag(j)
       wa(kmax) = wa(j)
       k = ipvt(j)
       ipvt(j) = ipvt(kmax)
       ipvt(kmax) = k

       !        compute the householder transformation to reduce the
       !        j-th column of a to a multiple of the j-th unit vector.

40     ajnorm = enorm_mod5V(m-j+1,a(j:,j))
       IF (ajnorm == zero) CYCLE
       IF (a(j,j) < zero) ajnorm = -ajnorm
       DO  i = j,m
          a(i,j) = a(i,j)/ajnorm
       END DO
       a(j,j) = a(j,j) + one

       !        apply the transformation to the remaining columns
       !        and update the norms.

       jp1 = j + 1
       DO  k = jp1,n
          sum = zero
          DO  i = j,m
             sum = sum + a(i,j)*a(i,k)
          END DO
          temp = sum/a(j,j)
          DO  i = j,m
             a(i,k) = a(i,k) - temp*a(i,j)
          END DO
          IF (.NOT.pivot .OR. rdiag(k) == zero) CYCLE
          temp = a(j,k)/rdiag(k)
          rdiag(k) = rdiag(k)*SQRT(MAX(zero,one-temp**2))
          IF (p05*(rdiag(k)/wa(k))**2 > epsmch) CYCLE
          rdiag(k) = enorm_mod5V(m-j,a(jp1:,k))
          wa(k) = rdiag(k)
       END DO
       rdiag(j) = -ajnorm
    END DO
    RETURN

    !     last card of subroutine qrfac.

  END SUBROUTINE qrfac



  SUBROUTINE qrsolv(n,r,ipvt,diag,qtb,x,sdiag)

    ! Code converted using TO_F90 by Alan Miller
    ! Date: 1999-12-09  Time: 12:46:21

    ! N.B. Arguments LDR & WA has been removed.

    INTEGER,INTENT(IN)::n
    REAL (PR),INTENT(INOUT)::r(:,:)
    INTEGER,INTENT(IN)::ipvt(:)
    REAL (PR),INTENT(IN)::diag(:)
    REAL (PR),INTENT(IN)::qtb(:)
    REAL (PR),INTENT(OUT)::x(:)
    REAL (PR),INTENT(OUT)::sdiag(:)


    !     **********

    !     subroutine qrsolv

    !     given an m by n matrix a,an n by n diagonal matrix d,
    !     and an m-vector b,the problem is to determine an x which
    !     solves the system

    !           a*x = b ,    d*x = 0 ,

    !     in the least squares sense.

    !     this subroutine completes the solution of the problem
    !     if it is provided with the necessary information from the
    !     qr factorization,with column pivoting,of a. that is,if
    !     a*p = q*r,where p is a permutation matrix,q has orthogonal
    !     columns,and r is an upper triangular matrix with diagonal
    !     elements of nonincreasing magnitude,then qrsolv expects
    !     the full upper triangle of r,the permutation matrix p,
    !     and the first n components of (q transpose)*b. the system
    !     a*x = b,d*x = 0,is then equivalent to

    !                  t       t
    !           r*z = q *b , p *d*p*z = 0 ,

    !     where x = p*z. if this system does not have full rank,
    !     then a least squares solution is obtained. on output qrsolv
    !     also provides an upper triangular matrix s such that

    !            t   t               t
    !           p *(a *a + d*d)*p = s *s .

    !     s is computed within qrsolv and may be of separate interest.

    !     the subroutine statement is

    !       subroutine qrsolv(n,r,ldr,ipvt,diag,qtb,x,sdiag,wa)

    !     where

    !       n is a positive integer input variable set to the order of r.

    !       r is an n by n array. on input the full upper triangle
    !         must contain the full upper triangle of the matrix r.
    !         on output the full upper triangle is unaltered,and the
    !         strict lower triangle contains the strict upper triangle
    !         (transposed) of the upper triangular matrix s.

    !       ldr is a positive integer input variable not less than n
    !         which specifies the leading dimension of the array r.

    !       ipvt is an integer input array of length n which defines the
    !         permutation matrix p such that a*p = q*r. column j of p
    !         is column ipvt(j) of the identity matrix.

    !       diag is an input array of length n which must contain the
    !         diagonal elements of the matrix d.

    !       qtb is an input array of length n which must contain the first
    !         n elements of the vector (q transpose)*b.

    !       x is an output array of length n which contains the least
    !         squares solution of the system a*x = b,d*x = 0.

    !       sdiag is an output array of length n which contains the
    !         diagonal elements of the upper triangular matrix s.

    !       wa is a work array of length n.

    !     subprograms called

    !       fortran-supplied ... ABS,SQRT

    !     argonne national laboratory. minpack project. march 1980.
    !     burton s. garbow,kenneth e. hillstrom,jorge j. more

    !     **********
    INTEGER::i,j,jp1,k,kp1,l,nsing
    REAL (PR)::COS,cotan,qtbpj,SIN,sum,TAN,temp,wa(n)
    REAL (PR),PARAMETER::p5 = 0.5_PR,p25 = 0.25_PR,zero = 0.0_PR

    !     copy r and (q transpose)*b to preserve input and initialize s.
    !     in particular,save the diagonal elements of r in x.

    DO  j = 1,n
       DO  i = j,n
          r(i,j) = r(j,i)
       END DO
       x(j) = r(j,j)
       wa(j) = qtb(j)
    END DO

    !     eliminate the diagonal matrix d using a givens rotation.

    DO  j = 1,n

       !        prepare the row of d to be eliminated,locating the
       !        diagonal element using p from the qr factorization.

       l = ipvt(j)
       IF (diag(l) == zero) CYCLE
       sdiag(j:n) = zero
       sdiag(j) = diag(l)

       !        the transformations to eliminate the row of d
       !        modify only a single element of (q transpose)*b
       !        beyond the first n,which is initially zero.

       qtbpj = zero
       DO  k = j,n

          !           determine a givens rotation which eliminates the
          !           appropriate element in the current row of d.

          IF (sdiag(k) == zero) CYCLE
          IF (ABS(r(k,k)) < ABS(sdiag(k))) THEN
             cotan = r(k,k)/sdiag(k)
             SIN = p5/SQRT(p25 + p25*cotan**2)
             COS = SIN*cotan
          ELSE
             TAN = sdiag(k)/r(k,k)
             COS = p5/SQRT(p25 + p25*TAN**2)
             SIN = COS*TAN
          END IF

          !           compute the modified diagonal element of r and
          !           the modified element of ((q transpose)*b,0).

          r(k,k) = COS*r(k,k) + SIN*sdiag(k)
          temp = COS*wa(k) + SIN*qtbpj
          qtbpj = -SIN*wa(k) + COS*qtbpj
          wa(k) = temp

          !           accumulate the tranformation in the row of s.

          kp1 = k + 1
          DO  i = kp1,n
             temp = COS*r(i,k) + SIN*sdiag(i)
             sdiag(i) = -SIN*r(i,k) + COS*sdiag(i)
             r(i,k) = temp
          END DO
       END DO

       !        store the diagonal element of s and restore
       !        the corresponding diagonal element of r.

       sdiag(j) = r(j,j)
       r(j,j) = x(j)
    END DO

    !     solve the triangular system for z. if the system is
    !     singular,then obtain a least squares solution.

    nsing = n
    DO  j = 1,n
       IF (sdiag(j) == zero .AND. nsing == n) nsing = j - 1
       IF (nsing < n) wa(j) = zero
    END DO

    DO  k = 1,nsing
       j = nsing - k + 1
       sum = zero
       jp1 = j + 1
       DO  i = jp1,nsing
          sum = sum + r(i,j)*wa(i)
       END DO
       wa(j) = (wa(j) - sum)/sdiag(j)
    END DO

    !     permute the components of z back to components of x.

    DO  j = 1,n
       l = ipvt(j)
       x(l) = wa(j)
    END DO
    RETURN

    !     last card of subroutine qrsolv.

  END SUBROUTINE qrsolv



  FUNCTION enorm_mod5V(n,x) RESULT(fn_val)

    ! Code converted using TO_F90 by Alan Miller
    ! Date: 1999-12-09  Time: 12:45:34

    INTEGER,INTENT(IN)::n
    REAL (PR),INTENT(IN)::x(:)
    REAL (PR)::fn_val


    !     **********

    !     function enorm_mod5V

    !     given an n-vector x,this function calculates the euclidean norm of x.

    !     the euclidean norm is computed by accumulating the sum of squares in
    !     three different sums.  The sums of squares for the small and large
    !     components are scaled so that no overflows occur.  Non-destructive
    !     underflows are permitted.  Underflows and overflows do not occur in the
    !     computation of the unscaled sum of squares for the intermediate
    !     components.  The definitions of small,intermediate and large components
    !     depend on two constants,rdwarf and rgiant.  The main restrictions on
    !     these constants are that rdwarf**2 not underflow and rgiant**2 not
    !     overflow.  The constants given here are suitable for every known computer.

    !     the function statement is

    !       REAL (PR) function enorm_mod5V(n,x)

    !     where

    !       n is a positive integer input variable.

    !       x is an input array of length n.

    !     subprograms called

    !       fortran-supplied ... ABS,SQRT

    !     argonne national laboratory. minpack project. march 1980.
    !     burton s. garbow,kenneth e. hillstrom,jorge j. more

    !     **********
    INTEGER::i
    REAL (PR)::agiant,floatn,s1,s2,s3,xabs,x1max,x3max
    REAL (PR),PARAMETER::one = 1.0_PR,zero = 0.0_PR,rdwarf = 3.834E-20_PR, &
         rgiant = 1.304E+19_PR

    s1 = zero
    s2 = zero
    s3 = zero
    x1max = zero
    x3max = zero
    floatn = n
    agiant = rgiant/floatn
    DO  i = 1,n
       xabs = ABS(x(i))
       IF (xabs > rdwarf .AND. xabs < agiant) GO TO 70
       IF (xabs <= rdwarf) GO TO 30

       !              sum for large components.

       IF (xabs <= x1max) GO TO 10
       s1 = one + s1*(x1max/xabs)**2
       x1max = xabs
       GO TO 20

10     s1 = s1 + (xabs/x1max)**2

20     GO TO 60

       !              sum for small components.

30     IF (xabs <= x3max) GO TO 40
       s3 = one + s3*(x3max/xabs)**2
       x3max = xabs
       GO TO 60

40     IF (xabs /= zero) s3 = s3 + (xabs/x3max)**2

60     CYCLE

       !           sum for intermediate components.

70     s2 = s2 + xabs**2
    END DO

    !     calculation of norm.

    IF (s1 == zero) GO TO 100
    fn_val = x1max*SQRT(s1 + (s2/x1max)/x1max)
    GO TO 120

100 IF (s2 == zero) GO TO 110
    IF (s2 >= x3max) fn_val = SQRT(s2*(one + (x3max/s2)*(x3max*s3)))
    IF (s2 < x3max) fn_val = SQRT(x3max*((s2/x3max) + (x3max*s3)))
    GO TO 120

110 fn_val = x3max*SQRT(s3)

120 RETURN

    !     last card of function enorm_mod5V.

  END FUNCTION enorm_mod5V



  SUBROUTINE fdjac2(fcn_mod5V,m,n,x,xd1,xd2,xd3,xd4,xd5,fvec,fjac,iflag,epsfcn)

    ! Code converted using TO_F90 by Alan Miller
    ! Date: 1999-12-09  Time: 12:45:44

    ! N.B. Arguments LDFJAC & WA have been removed.

    INTEGER,INTENT(IN)::m
    INTEGER,INTENT(IN)::n
    REAL (PR),INTENT(INOUT)::x(n)
    REAL (PR),INTENT(IN)::fvec(m)
    REAL (PR),INTENT(OUT)::fjac(:,:)    ! fjac(ldfjac,n)
    INTEGER,INTENT(inout)::iflag
    REAL (PR),INTENT(IN)::epsfcn
    REAL (PR),INTENT(IN)::xd1(:),xd2(:),xd3(:),xd4(:),xd5(:)

    INTERFACE
       SUBROUTINE fcn_mod5V(m,n,x,xd1,xd2,xd3,xd4,xd5,fvec,iflag)
         use common
         IMPLICIT NONE
         !INTEGER,PARAMETER::PR = SELECTED_REAL_KIND(12,60)
         INTEGER,INTENT(IN)::m,n
         REAL (PR),INTENT(IN)::x(:)
         REAL (PR),INTENT(OUT)::fvec(:)
         INTEGER,INTENT(inout)::iflag
         REAL (PR),INTENT(IN)::xd1(:),xd2(:),xd3(:),xd4(:),xd5(:)

       END SUBROUTINE fcn_mod5V
    END INTERFACE

    !     **********

    !     subroutine fdjac2

    !     this subroutine computes a forward-difference approximation
    !     to the m by n jacobian matrix associated with a specified
    !     problem of m functions in n variables.

    !     the subroutine statement is

    !       subroutine fdjac2(fcn,m,n,x,fvec,fjac,ldfjac,iflag,epsfcn,wa)

    !     where

    !       fcn is the name of the user-supplied subroutine which
    !         calculates the functions. fcn must be declared
    !         in an external statement in the user calling
    !         program,and should be written as follows.

    !         subroutine fcn(m,n,x,fvec,iflag)
    !         integer m,n,iflag
    !         REAL (PR) x(n),fvec(m)
    !         ----------
    !         calculate the functions at x and
    !         return this vector in fvec.
    !         ----------
    !         return
    !         end

    !         the value of iflag should not be changed by fcn unless
    !         the user wants to terminate execution of fdjac2.
    !         in this case set iflag to a negative integer.

    !       m is a positive integer input variable set to the number of functions.

    !       n is a positive integer input variable set to the number
    !         of variables. n must not exceed m.

    !       x is an input array of length n.

    !       fvec is an input array of length m which must contain the
    !         functions evaluated at x.

    !       fjac is an output m by n array which contains the
    !         approximation to the jacobian matrix evaluated at x.

    !       ldfjac is a positive integer input variable not less than m
    !         which specifies the leading dimension of the array fjac.

    !       iflag is an integer variable which can be used to terminate
    !         the execution of fdjac2.  see description of fcn.

    !       epsfcn is an input variable used in determining a suitable
    !         step length for the forward-difference approximation. this
    !         approximation assumes that the relative errors in the
    !         functions are of the order of epsfcn. if epsfcn is less
    !         than the machine precision,it is assumed that the relative
    !         errors in the functions are of the order of the machine precision.

    !       wa is a work array of length m.

    !     subprograms called

    !       user-supplied ...... fcn

    !       minpack-supplied ... PRmpar

    !       fortran-supplied ... ABS,MAX,SQRT

    !     argonne national laboratory. minpack project. march 1980.
    !     burton s. garbow,kenneth e. hillstrom,jorge j. more

    !     **********
    INTEGER::i,j
    REAL (PR)::eps,epsmch,h,temp,wa(m)
    REAL (PR),PARAMETER::zero = 0.0_PR

    !     epsmch is the machine precision.

    epsmch = EPSILON(zero)

    eps = SQRT(MAX(epsfcn,epsmch))
    DO  j = 1,n
       temp = x(j)
       h = eps*ABS(temp)
       IF (h == zero) h = eps
       x(j) = temp + h
       CALL fcn_mod5V(m,n,x,xd1,xd2,xd3,xd4,xd5,wa,iflag)
       IF (iflag < 0) EXIT
       x(j) = temp
       DO  i = 1,m
          fjac(i,j) = (wa(i) - fvec(i))/h
       END DO
    END DO

    RETURN

    !     last card of subroutine fdjac2.

  END SUBROUTINE fdjac2


  SUBROUTINE covar(n,r,ldr,ipvt,tol,wa)

    ! Code converted using TO_F90 by Alan Miller
    ! Date: 2013-01-12  Time: 18:09:08

    INTEGER,intent(in)::n
    REAL(PR),intent(out)::r(ldr,n)
    INTEGER,intent(inout)::ldr
    INTEGER,intent(in)::ipvt(n)
    REAL(PR),intent(in)::tol
    REAL(PR),intent(out)::wa(n)




    !     **********

    !     subroutine covar

    !     given an m by n matrix a,the problem is to determine
    !     the covariance matrix corresponding to a,defined as

    !              t
    !     inverse(a *a) .

    !     this subroutine completes the solution of the problem
    !     if it is provided with the necessary information from the
    !     qr factorization,with column pivoting,of a. that is,if
    !     a*p = q*r,where p is a permutation matrix,q has orthogonal
    !     columns,and r is an upper triangular matrix with diagonal
    !     elements of nonincreasing magnitude,then covar expects
    !     the full upper triangle of r and the permutation matrix p.
    !     the covariance matrix is then computed as

    !                t     t
    !     p*inverse(r *r)*p  .

    !     if a is nearly rank deficient,it may be desirable to compute
    !     the covariance matrix corresponding to the linearly independent
    !     columns of a. to define the numerical rank of a,covar uses
    !     the tolerance tol. if l is the largest integer such that

    !     abs(r(l,l)) .gt. tol*abs(r(1,1)) ,

    !     then covar computes the covariance matrix corresponding to
    !     the first l columns of r. for k greater than l,column
    !     and row ipvt(k) of the covariance matrix are set to zero.

    !     the subroutine statement is

    !     subroutine covar(n,r,ldr,ipvt,tol,wa)

    !     where

    !     n is a positive integer input variable set to the order of r.

    !     r is an n by n array. on input the full upper triangle must
    !     contain the full upper triangle of the matrix r. on output
    !     r contains the square symmetric covariance matrix.

    !     ldr is a positive integer input variable not less than n
    !     which specifies the leading dimension of the array r.

    !     ipvt is an integer input array of length n which defines the
    !     permutation matrix p such that a*p = q*r. column j of p
    !     is column ipvt(j) of the identity matrix.

    !     tol is a nonnegative input variable used to define the
    !     numerical rank of a in the manner described above.

    !     wa is a work array of length n.

    !     subprograms called

    !     fortran-supplied ... dabs

    !     argonne national laboratory. minpack project. august 1980.
    !     burton s. garbow,kenneth e. hillstrom,jorge j. more

    !     **********
    INTEGER::i,ii,j,jj,k,km1,l
    LOGICAL::sing
    REAL(PR)::one,temp,tolr,zero
    DATA one,zero /1.0D0,0.0D0/

    !     form the inverse of r in the full upper triangle of r.

    tolr = tol*DABS(r(1,1))
    l = 0
    DO  k = 1,n
       IF (DABS(r(k,k)) <= tolr) EXIT
       r(k,k) = one/r(k,k)
       km1 = k - 1
       IF (km1 < 1) GO TO 30
       DO  j = 1,km1
          temp = r(k,k)*r(j,k)
          r(j,k) = zero
          DO  i = 1,j
             r(i,k) = r(i,k) - temp*r(i,j)
          END DO
       END DO
30     CONTINUE
       l = k
    END DO
    ! 50 CONTINUE

    !     form the full upper triangle of the inverse of (r transpose)*r
    !     in the full upper triangle of r.

    IF (l < 1) GO TO 110
    DO  k = 1,l
       km1 = k - 1
       IF (km1 < 1) GO TO 80
       DO  j = 1,km1
          temp = r(j,k)
          DO  i = 1,j
             r(i,j) = r(i,j) + temp*r(i,k)
          END DO
       END DO
80     CONTINUE
       temp = r(k,k)
       DO  i = 1,k
          r(i,k) = temp*r(i,k)
       END DO
    END DO
110 CONTINUE

    !     form the full lower triangle of the covariance matrix
    !     in the strict lower triangle of r and in wa.

    DO  j = 1,n
       jj = ipvt(j)
       sing = j > l
       DO  i = 1,j
          IF (sing) r(i,j) = zero
          ii = ipvt(i)
          IF (ii > jj) r(ii,jj) = r(i,j)
          IF (ii < jj) r(jj,ii) = r(i,j)
       END DO
       wa(jj) = r(j,j)
    END DO

    !     symmetrize the covariance matrix in r.

    DO  j = 1,n
       DO  i = 1,j
          r(i,j) = r(j,i)
       END DO
       r(j,j) = wa(j)
    END DO
    RETURN

    !     last card of subroutine covar.

  END SUBROUTINE covar



  SUBROUTINE lmdif_mod5V(fcn_mod5v,m,n,x,xd1,xd2,xd3,xd4,xd5,fvec,ftol,xtol,gtol,maxfev,epsfcn, &
       diag,mode,factor,nprint,info,nfev,fjac,ipvt,qtf)
    !This modification of LMDIF can take up to five data arrays (define din DOUBLE PRECISION)
    !to be passed to FCN for residuals computation

    ! Code converted using TO_F90 by Alan Miller
    ! Date: 1999-12-09  Time: 12:45:59

    ! N.B. Arguments LDFJAC,WA1,WA2,WA3 & WA4 have been removed.

    INTEGER,INTENT(IN)::m
    INTEGER,INTENT(IN)::n
    REAL (PR),INTENT(INOUT)::x(:)
    REAL (PR),INTENT(OUT)::fvec(:)
    REAL (PR),INTENT(IN)::ftol
    REAL (PR),INTENT(IN)::xtol
    REAL (PR),INTENT(INOUT)::gtol
    INTEGER,INTENT(INOUT)::maxfev
    REAL (PR),INTENT(INOUT)::epsfcn
    REAL (PR),INTENT(OUT)::diag(:)
    INTEGER,INTENT(IN)::mode
    REAL (PR),INTENT(IN)::factor
    INTEGER,INTENT(IN)::nprint
    INTEGER,INTENT(OUT)::info
    INTEGER,INTENT(OUT)::nfev
    REAL (PR),INTENT(OUT)::fjac(:,:)    ! fjac(ldfjac,n)
    INTEGER,INTENT(OUT)::ipvt(:)
    REAL (PR),INTENT(OUT)::qtf(:)

    REAL (PR),INTENT(IN)::xd1(:),xd2(:),xd3(:),xd4(:),xd5(:)

    ! EXTERNAL fcn_mod5V

    INTERFACE
       SUBROUTINE fcn_mod5V(m,n,x,xd1,xd2,xd3,xd4,xd5,fvec,iflag)
         use common
         IMPLICIT NONE
         !INTEGER,PARAMETER::PR = SELECTED_REAL_KIND(12,60)
         INTEGER,INTENT(IN)::m,n
         REAL (PR),INTENT(IN)::x(:)
         REAL (PR),INTENT(OUT)::fvec(:)
         INTEGER,INTENT(inout)::iflag
         REAL (PR),INTENT(IN)::xd1(:),xd2(:),xd3(:),xd4(:),xd5(:)
       END SUBROUTINE fcn_mod5V
    END INTERFACE

    !     **********

    !     subroutine lmdif

    !     the purpose of lmdif is to minimize the sum of the squares of
    !     m nonlinear functions in n variables by a modification of
    !     the levenberg-marquardt algorithm. The user must provide a
    !     subroutine which calculates the functions. The jacobian is
    !     then calculated by a forward-difference approximation.

    !     the subroutine statement is

    !       subroutine lmdif(fcn_mod5V,m,n,x,xd1,xd2,xd3,xd4,xd5,fvec,ftol,xtol,gtol,maxfev,epsfcn,
    !                        diag,mode,factor,nprint,info,nfev,fjac,
    !                        ldfjac,ipvt,qtf,wa1,wa2,wa3,wa4)

    !     where

    !       fcn is the name of the user-supplied subroutine which
    !         calculates the functions. fcn must be declared
    !         in an external statement in the user calling
    !         program,and should be written as follows.

    !         subroutine fcn_mod5V(m,n,x,xd1,xd2,xd3,xd4,xd5,fvec,iflag)
    !         integer m,n,iflag
    !         REAL (PR) x(:),fvec(m),xd1(:),xd2(:),xd3(:),xd4(:),xd5(:)
    !         ----------
    !         calculate the functions at x and
    !         return this vector in fvec.
    !         ----------
    !         return
    !         end

    !         the value of iflag should not be changed by fcn unless
    !         the user wants to terminate execution of lmdif.
    !         in this case set iflag to a negative integer.

    !       m is a positive integer input variable set to the number
    !         of functions.

    !       n is a positive integer input variable set to the number
    !         of variables. n must not exceed m.

    !       x is an array of length n. on input x must contain
    !         an initial estimate of the solution vector. on output x
    !         contains the final estimate of the solution vector.

    !       fvec is an output array of length m which contains
    !         the functions evaluated at the output x.

    !       ftol is a nonnegative input variable. termination
    !         occurs when both the actual and predicted relative
    !         reductions in the sum of squares are at most ftol.
    !         therefore,ftol measures the relative error desired
    !         in the sum of squares.

    !       xtol is a nonnegative input variable. termination
    !         occurs when the relative error between two consecutive
    !         iterates is at most xtol. therefore,xtol measures the
    !         relative error desired in the approximate solution.

    !       gtol is a nonnegative input variable. termination
    !         occurs when the cosine of the angle between fvec and
    !         any column of the jacobian is at most gtol in absolute
    !         value. therefore,gtol measures the orthogonality
    !         desired between the function vector and the columns
    !         of the jacobian.

    !       maxfev is a positive integer input variable. termination
    !         occurs when the number of calls to fcn is at least
    !         maxfev by the end of an iteration.

    !       epsfcn is an input variable used in determining a suitable
    !         step length for the forward-difference approximation. this
    !         approximation assumes that the relative errors in the
    !         functions are of the order of epsfcn. if epsfcn is less
    !         than the machine precision,it is assumed that the relative
    !         errors in the functions are of the order of the machine
    !         precision.

    !       diag is an array of length n. if mode = 1 (see
    !         below),diag is internally set. if mode = 2,diag
    !         must contain positive entries that serve as
    !         multiplicative scale factors for the variables.

    !       mode is an integer input variable. if mode = 1,the
    !         variables will be scaled internally. if mode = 2,
    !         the scaling is specified by the input diag. other
    !         values of mode are equivalent to mode = 1.

    !       factor is a positive input variable used in determining the
    !         initial step bound. this bound is set to the product of
    !         factor and the euclidean norm of diag*x if nonzero,or else
    !         to factor itself. in most cases factor should lie in the
    !         interval (.1,100.). 100. is a generally recommended value.

    !       nprint is an integer input variable that enables controlled
    !         printing of iterates if it is positive. in this case,
    !         fcn is called with iflag = 0 at the beginning of the first
    !         iteration and every nprint iterations thereafter and
    !         immediately prior to return,with x and fvec available
    !         for printing. if nprint is not positive,no special calls
    !         of fcn with iflag = 0 are made.

    !       info is an integer output variable.  If the user has terminated
    !         execution,info is set to the (negative) value of iflag.
    !         See description of fcn.  Otherwise,info is set as follows.

    !         info = 0  improper input parameters.

    !         info = 1  both actual and predicted relative reductions
    !                   in the sum of squares are at most ftol.

    !         info = 2  relative error between two consecutive iterates
    !                   is at most xtol.

    !         info = 3  conditions for info = 1 and info = 2 both hold.

    !         info = 4  the cosine of the angle between fvec and any column of
    !                   the Jacobian is at most gtol in absolute value.

    !         info = 5  number of calls to fcn has reached or exceeded maxfev.

    !         info = 6  ftol is too small. no further reduction in
    !                   the sum of squares is possible.

    !         info = 7  xtol is too small. no further improvement in
    !                   the approximate solution x is possible.

    !         info = 8  gtol is too small. fvec is orthogonal to the
    !                   columns of the jacobian to machine precision.

    !       nfev is an integer output variable set to the number of calls to fcn.

    !       fjac is an output m by n array. the upper n by n submatrix
    !         of fjac contains an upper triangular matrix r with
    !         diagonal elements of nonincreasing magnitude such that

    !                t     t           t
    !               p *(jac *jac)*p = r *r,

    !         where p is a permutation matrix and jac is the final calculated
    !         Jacobian.  Column j of p is column ipvt(j) (see below) of the
    !         identity matrix. the lower trapezoidal part of fjac contains
    !         information generated during the computation of r.

    !       ldfjac is a positive integer input variable not less than m
    !         which specifies the leading dimension of the array fjac.

    !       ipvt is an integer output array of length n. ipvt
    !         defines a permutation matrix p such that jac*p = q*r,
    !         where jac is the final calculated jacobian,q is
    !         orthogonal (not stored),and r is upper triangular
    !         with diagonal elements of nonincreasing magnitude.
    !         column j of p is column ipvt(j) of the identity matrix.

    !       qtf is an output array of length n which contains
    !         the first n elements of the vector (q transpose)*fvec.

    !       wa1,wa2,and wa3 are work arrays of length n.

    !       wa4 is a work array of length m.

    !     subprograms called

    !       user-supplied ...... fcn_mod5V

    !       minpack-supplied ... PRmpar,enorm_mod5V,fdjac2,lmpar,qrfac

    !       fortran-supplied ... dabs,dmax1,dmin1,dsqrt,mod

    !     argonne national laboratory. minpack project. march 1980.
    !     burton s. garbow,kenneth e. hillstrom,jorge j. more

    !     **********
    INTEGER::i,iflag,iter,j,l
    REAL (PR)::actred,delta,dirder,epsmch,fnorm,fnorm1,gnorm, &
         par,pnorm,prered,ratio,sum,temp,temp1,temp2,xnorm
    REAL (PR)::wa1(n),wa2(n),wa3(n),wa4(m)
    REAL (PR),PARAMETER::one = 1.0_PR,p1 = 0.1_PR,p5 = 0.5_PR, &
         p25 = 0.25_PR,p75 = 0.75_PR,p0001 = 0.0001_PR,&
         zero = 0.0_PR

    !     epsmch is the machine precision.

    epsmch = EPSILON(zero)

    info = 0
    iflag = 0
    nfev = 0

    !     check the input parameters for errors.

    IF (n <= 0 .OR. m < n .OR. ftol < zero .OR. xtol < zero .OR. gtol < zero  &
         .OR. maxfev <= 0 .OR. factor <= zero) GO TO 300
    IF (mode /= 2) GO TO 20
    DO  j = 1,n
       IF (diag(j) <= zero) GO TO 300
    END DO

    !     evaluate the function at the starting point and calculate its norm.

20  iflag = 1
    CALL fcn_mod5V(m,n,x,xd1,xd2,xd3,xd4,xd5,fvec,iflag)
    nfev = 1
    IF (iflag < 0) GO TO 300
    fnorm = enorm_mod5V(m,fvec)

    !     initialize levenberg-marquardt parameter and iteration counter.

    par = zero
    iter = 1

    !     beginning of the outer loop.


    !        calculate the jacobian matrix.

30  iflag = 2
    CALL fdjac2(fcn_mod5V,m,n,x,xd1,xd2,xd3,xd4,xd5,fvec,fjac,iflag,epsfcn)
    nfev = nfev + n
    IF (iflag < 0) GO TO 300

    !        if requested,call fcn to enable printing of iterates.

    IF (nprint <= 0) GO TO 40
    iflag = 0
    IF (MOD(iter-1,nprint) == 0) CALL fcn_mod5V(m,n,x,xd1,xd2,xd3,xd4,xd5,fvec,iflag)
    IF (iflag < 0) GO TO 300

    !        compute the qr factorization of the jacobian.

40  CALL qrfac(m,n,fjac,.true.,ipvt,wa1,wa2)

    !        on the first iteration and if mode is 1,scale according
    !        to the norms of the columns of the initial jacobian.

    IF (iter /= 1) GO TO 80
    IF (mode == 2) GO TO 60
    DO  j = 1,n
       diag(j) = wa2(j)
       IF (wa2(j) == zero) diag(j) = one
    END DO

    !        on the first iteration,calculate the norm of the scaled x
    !        and initialize the step bound delta.

60  wa3(1:n) = diag(1:n)*x(1:n)
    xnorm = enorm_mod5V(n,wa3)
    delta = factor*xnorm
    IF (delta == zero) delta = factor

    !        form (q transpose)*fvec and store the first n components in qtf.

80  wa4(1:m) = fvec(1:m)
    DO  j = 1,n
       IF (fjac(j,j) == zero) GO TO 120
       sum = DOT_PRODUCT( fjac(j:m,j),wa4(j:m) )
       temp = -sum/fjac(j,j)
       DO  i = j,m
          wa4(i) = wa4(i) + fjac(i,j)*temp
       END DO
120    fjac(j,j) = wa1(j)
       qtf(j) = wa4(j)
    END DO

    !        compute the norm of the scaled gradient.

    gnorm = zero
    IF (fnorm == zero) GO TO 170
    DO  j = 1,n
       l = ipvt(j)
       IF (wa2(l) == zero) CYCLE
       sum = zero
       DO  i = 1,j
          sum = sum + fjac(i,j)*(qtf(i)/fnorm)
       END DO
       gnorm = MAX(gnorm,ABS(sum/wa2(l)))
    END DO

    !        test for convergence of the gradient norm.

170 IF (gnorm <= gtol) info = 4
    IF (info /= 0) GO TO 300

    !        rescale if necessary.

    IF (mode == 2) GO TO 200
    DO  j = 1,n
       diag(j) = MAX(diag(j),wa2(j))
    END DO

    !        beginning of the inner loop.

    !           determine the levenberg-marquardt parameter.

200 CALL lmpar(n,fjac,ipvt,diag,qtf,delta,par,wa1,wa2)

    !           store the direction p and x + p. calculate the norm of p.

    DO  j = 1,n
       wa1(j) = -wa1(j)
       wa2(j) = x(j) + wa1(j)
       wa3(j) = diag(j)*wa1(j)
    END DO
    pnorm = enorm_mod5V(n,wa3)

    !           on the first iteration,adjust the initial step bound.

    IF (iter == 1) delta = MIN(delta,pnorm)

    !           evaluate the function at x + p and calculate its norm.

    iflag = 1
    CALL fcn_mod5V(m,n,wa2,xd1,xd2,xd3,xd4,xd5,wa4,iflag)
    nfev = nfev + 1
    IF (iflag < 0) GO TO 300
    fnorm1 = enorm_mod5V(m,wa4)

    !           compute the scaled actual reduction.

    actred = -one
    IF (p1*fnorm1 < fnorm) actred = one - (fnorm1/fnorm)**2

    !           compute the scaled predicted reduction and
    !           the scaled directional derivative.

    DO  j = 1,n
       wa3(j) = zero
       l = ipvt(j)
       temp = wa1(l)
       DO  i = 1,j
          wa3(i) = wa3(i) + fjac(i,j)*temp
       END DO
    END DO
    temp1 = enorm_mod5V(n,wa3)/fnorm
    temp2 = (SQRT(par)*pnorm)/fnorm
    prered = temp1**2 + temp2**2/p5
    dirder = -(temp1**2 + temp2**2)

    !           compute the ratio of the actual to the predicted reduction.

    ratio = zero
    IF (prered /= zero) ratio = actred/prered

    !           update the step bound.

    IF (ratio <= p25) THEN
       IF (actred >= zero) temp = p5
       IF (actred < zero) temp = p5*dirder/(dirder + p5*actred)
       IF (p1*fnorm1 >= fnorm .OR. temp < p1) temp = p1
       delta = temp*MIN(delta,pnorm/p1)
       par = par/temp
    ELSE
       IF (par /= zero .AND. ratio < p75) GO TO 260
       delta = pnorm/p5
       par = p5*par
    END IF

    !           test for successful iteration.

260 IF (ratio < p0001) GO TO 290

    !           successful iteration. update x,fvec,and their norms.

    DO  j = 1,n
       x(j) = wa2(j)
       wa2(j) = diag(j)*x(j)
    END DO
    fvec(1:m) = wa4(1:m)
    xnorm = enorm_mod5V(n,wa2)
    fnorm = fnorm1
    iter = iter + 1

    !           tests for convergence.

290 IF (ABS(actred) <= ftol .AND. prered <= ftol .AND. p5*ratio <= one) info = 1
    IF (delta <= xtol*xnorm) info = 2
    IF (ABS(actred) <= ftol .AND. prered <= ftol  &
         .AND. p5*ratio <= one .AND. info == 2) info = 3
    IF (info /= 0) GO TO 300

    !           tests for termination and stringent tolerances.

    IF (nfev >= maxfev) info = 5
    IF (ABS(actred) <= epsmch .AND. prered <= epsmch  &
         .AND. p5*ratio <= one) info = 6
    IF (delta <= epsmch*xnorm) info = 7
    IF (gnorm <= epsmch) info = 8
    IF (info /= 0) GO TO 300

    !           end of the inner loop. repeat if iteration unsuccessful.

    IF (ratio < p0001) GO TO 200

    !        end of the outer loop.

    GO TO 30

    !     termination,either normal or user imposed.

300 IF (iflag < 0) info = iflag
    iflag = 0
    IF (nprint > 0) CALL fcn_mod5V(m,n,x,xd1,xd2,xd3,xd4,xd5,fvec,iflag)
    RETURN

    !     last card of subroutine lmdif.

  END SUBROUTINE lmdif_mod5V



END MODULE Levenberg_Marquardt_mod5V
