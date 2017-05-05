!!! --- FITS module ---!!!
!!! This module provides the basic routines to open, read and write the most standard FITS files
!!! Some routines for basic HEADER operations are provided too
!!!
!!! fits_checkdim : gives back the number of dimensions of a FITS file 
!!! fits_actv     : open a FITS file without performing any operation on it
!!! fits_actv_rw  : open a FITS file in read-write mode
!!! fits_size_1d  : gives back the size of a monodimensional FITS file
!!! fits_size_2d  : gives back the two size (x_axis,y_axis) of a bidimensional FITS file
!!! wcs_size      : gives back the size a given axis (FITS must have been already open)
!!! wcs_info      : gives back basic WCS informations for a given axis
!!! wcs_to_scale  : gives back all the necessary WCS information to rebuild a WCS axis
!!! ccfscale_to_wcs : convert a CCF scale into WCS and write it in the FITS file
!!! fits_open_s1d : open a 1D, double precision (PD=8) FITS file
!!! fits_open_s2d : open a 2D, double precision (PD=8) FITS file
!!! fits_write_1d : write a 1D, double precision (PD=8) FITS file
!!! fits_write_2d : write a 2D, double precision (PD=8) FITS file
!!! fits_close    : close a FITS file previously opened with fits_actv
!!! fits_check    : check the existence of a FITS file (sts_out=1 -> the file exists)
!!! fits_delete   : delete a given fits file, if it exists
!!! reset_sts     : reset the status of a CFITSIO subroutine, and provides informations in case of error

module fits
  use common
  
  implicit none

contains
  
  subroutine fits_checkdim(file,naxis,verbose)

    integer :: unit
    character (len = nch_file), intent(in) :: file
    integer, intent(in), optional :: verbose

    integer, intent(out) :: naxis
    integer :: verb, rwmod, sts, blksz
    call get_lun(unit)
    
    rwmod = 0 ; sts = 0 ;  verb = 0
    if (present(verbose)) verb=verbose
    call FTOPEN(unit,file,rwmod,blksz,sts)
    call reset_sts(sts,1,verb)
    call FTGIDM(unit,naxis,sts)
    call reset_sts(sts,1,verb)
    call FTCLOS(unit, sts)
    call reset_sts(sts,1,verb)
    
    call free_lun(unit)

    return
    
  end subroutine fits_checkdim


  subroutine fits_actv(unit,file,verbose)
    !     open a fits file $file without reading the image
    !     Useful if you only have to check his existence or only do
    !     header operations

    integer, intent(in) :: unit
    character (len=nch_file), intent(in) :: file
    integer, intent(in), optional :: verbose

    integer rwmod, blksz, grp
    integer sts
    integer nrows, ncols, maxdim, naxis, naxes(2), bitpix, verb


    rwmod = 0 ;  grp = 0 ;  maxdim = 2 ;  sts = 0 ;  verb = 0
    if (present(verbose)) verb = verbose
    !     opens fits 'file'
    if (verb.gt.0) write(verb,'(A25,A,A9,I4)') '  *** FTSACTV: Opening: ',trim(file), ' in unit ',unit
    call FTOPEN(unit,file,rwmod,blksz,sts)
    call reset_sts(sts,1,verb)
    return
  end subroutine fits_actv

  subroutine fits_actv_rw(unit,file,verbose)
    !     open a fits file $file without reading the image
    !     File is open in read-write mode, so use this routine carefully!!

    integer, intent(in) :: unit
    character (len=nch_file), intent(in) :: file
    integer, intent(in), optional :: verbose

    integer rwmod, blksz, grp
    integer sts
    integer nrows, ncols, maxdim, naxis, naxes(2), bitpix, verb

    
    rwmod = 1 ;  grp = 0 ;  maxdim = 2 ;  sts = 0 ;  verb = 0
    if (present(verbose)) verb = verbose
    !     opens fits 'file'
    if (verb.gt.0) write(verb,'(A25,A,A9,I4,A)') '  *** FTSACTV: Opening: ',trim(file), ' in unit ',unit,' in RW mode!!' 
    call FTOPEN(unit,file,rwmod,blksz,sts)
    call reset_sts(sts,1,verb)
    return
  end subroutine fits_actv_rw

  subroutine fits_size_1d(file,x_axis,verbose)

    integer :: unit
    character (len = nch_file), intent(in) :: file
    integer, intent(in), optional :: verbose

    integer, intent(out) :: x_axis
    integer :: rwmod, blksz, sts,  maxdim, naxis, bitpix, verb
    logical :: anyf

    rwmod = 0 ;  maxdim = 1 ;  sts = 0 ;  verb = 0
    if (present(verbose)) verb=verbose
    
    call get_lun(unit)
    call FTOPEN(unit,file,rwmod,blksz,sts); call reset_sts(sts,1,verb)
    !     extracts basic header info from fits
    call FTGIPR(unit, maxdim, bitpix, naxis, x_axis, sts);  call reset_sts(sts,1,verb)
    call FTCLOS(unit, sts); call reset_sts(sts,1,verb)
    call free_lun(unit)

    return
  end subroutine fits_size_1d

  subroutine fits_size_2d(file,x_axis,y_axis,verbose)

    integer :: unit
    character (len = nch_file), intent(in) :: file
    integer, intent(in), optional :: verbose

    integer, intent(out) :: x_axis, y_axis
    integer :: rwmod, blksz, sts,  maxdim, naxis, bitpix, verb
    integer, dimension(2) :: naxes
    logical :: anyf

    rwmod = 0 ;  maxdim = 2 ;  sts = 0 ;  verb = 0
    if (present(verbose)) verb=verbose
    
    call get_lun(unit)
    call FTOPEN(unit,file,rwmod,blksz,sts); call reset_sts(sts,1,verb)
    !     extracts basic header info from fits
    call FTGIPR(unit, maxdim, bitpix, naxis, naxes, sts);  call reset_sts(sts,1,verb)
    call FTCLOS(unit, sts); call reset_sts(sts,1,verb)
    x_axis = naxes(1)
    y_axis = naxes(2)
    call free_lun(unit)

    return
  end subroutine fits_size_2d
    

  subroutine wcs_size(unit,naxis,size,verbose)
    integer, intent(in)  :: unit,naxis
    integer, intent(in), optional :: verbose
    integer, intent(out) :: size
    
    character (len=1) :: nc
    integer ::  sts, verb
    character (len=64) :: comment

    size = 0
    write(nc,'(I1)') naxis
    
    call FTGKYJ(unit,'NAXIS'//nc,size  ,comment,sts); call reset_sts(sts,0,verb,'CRVAL'//nc)
    return
  end subroutine wcs_size
  
  subroutine wcs_info(file,naxis,size,crval,cdelt,verbose)
    character (len=nch_file), intent(in)  :: file
    integer, intent(in)  :: naxis
    integer, intent(in), optional :: verbose
    integer, intent(out) :: size
    real   , intent(out) :: crval, cdelt

    character (len=1) :: nc
    integer ::  sts, verb, unit
    character (len=64) :: comment
    
    integer :: rwmod, blksz

    rwmod = 0 ; sts = 0 ;  verb = 0
    size = 0 ; cdelt = 0.00; crval = 0.00
    write(nc,'(I1)') naxis

    call get_lun(unit)
    call FTOPEN(unit,file,rwmod,blksz,sts) ;  call reset_sts(sts,1,verb)
    call FTGKYJ(unit,'NAXIS'//nc,size ,comment,sts); call reset_sts(sts,0,verb,'NAXIS'//nc)
    call FTGKYE(unit,'CDELT'//nc,cdelt,comment,sts); call reset_sts(sts,0,verb,'CDELT'//nc)
    call FTGKYE(unit,'CRVAL'//nc,crval,comment,sts); call reset_sts(sts,0,verb,'CRVAL'//nc)
    call FTCLOS(unit, sts)
    call free_lun(unit)
    
    return
  end subroutine wcs_info


  subroutine wcs_to_scale(unit, naxis, scale, verbose, ctype, cunit, crval, crpix, cdelt)
    integer, intent(in) :: unit,naxis

    real (kind=8), dimension(:), intent(out) :: scale

    integer, intent(in), optional :: verbose
    character (len=nch_key), intent(out), optional :: ctype, cunit
    real, intent(out), optional :: crval, crpix, cdelt

    character (len=nch_key) :: ctype_i, cunit_i
    real :: crval_i, crpix_i, cdelt_i
    character (len=1) :: nc
    ! _i : internal variables

    integer ::  sts, verb, nsize, i
    character (len=64) :: comment
    
    !fits with more than 9 dimensions are not supported
    write(nc,'(I1)') naxis
    
    nsize = 0 ;  cunit_i = '' ;  ctype_i = '' ;  crval_i = 1.0 ;  crpix_i = 1.0 ;  cdelt_i = 1.0 
    sts = 0 ; verb = 0
    
    if (present(verbose)) verb=verbose

    call FTGKYJ(unit,'NAXIS'//nc,nsize  ,comment,sts); call reset_sts(sts,0,verb,'CRVAL'//nc)
    call FTGKYS(unit,'CTYPE'//nc,ctype_i,comment,sts); call reset_sts(sts,0,verb,'CTYPE'//nc)
    call FTGKYS(unit,'CUNIT'//nc,cunit_i,comment,sts); call reset_sts(sts,0,verb,'CUNIT'//nc)
    call FTGKYE(unit,'CDELT'//nc,cdelt_i,comment,sts); call reset_sts(sts,0,verb,'CDELT'//nc)
    call FTGKYE(unit,'CRVAL'//nc,crval_i,comment,sts); call reset_sts(sts,0,verb,'CRVAL'//nc)
    call FTGKYE(unit,'CRPIX'//nc,crpix_i,comment,sts); call reset_sts(sts,0,verb,'CRPIX'//nc)

    if (size(scale,1) .ne. nsize) stop
    do i=1,nsize
       scale(i)=dble(crval_i) + dble((real(i)-crpix_i)*cdelt_i)
    end do
    
    if (present(ctype)) ctype=ctype_i
    if (present(cunit)) cunit=cunit_i
    if (present(crval)) crval=crval_i
    if (present(crpix)) crpix=crpix_i
    if (present(cdelt)) cdelt=cdelt_i

    return    
  end subroutine wcs_to_scale

  subroutine scale_to_wcs(unit, naxis, verbose, nsize, ctype, cunit, crval, crpix, cdelt, bunit, &
       nsize_com, ctype_com, cunit_com, crval_com, crpix_com, cdelt_com, bunit_com)
    integer, intent(in) :: unit, naxis 
    integer, intent(in), optional :: verbose,  nsize
    character (len=nch_key), intent(in), optional :: ctype, cunit, bunit
    real, intent(in), optional :: crval, crpix, cdelt
    character (len=64), optional :: &
         nsize_com, ctype_com, cunit_com, crval_com, crpix_com, cdelt_com, bunit_com
    
    
    character (len=64) :: comment
    character (len=1) :: nc
    ! _i : internal variables
    
    integer ::  sts, verb, i
    
    !fits with more than 9 dimensions are not supported
    write(nc,'(I1)') naxis
    if (naxis.eq.0) nc=''
    sts = 0 ; verb = 0
    if (present(verbose)) verb=verbose
    
!!$    call FTPKYJ(unit,'NAXIS'//nc,nsize,comment,sts); call reset_sts(sts,0,verb,'CRVAL'//nc//'      ')
!!$    call FTPKYS(unit,'CTYPE'//nc,ctype,comment,sts); call reset_sts(sts,0,verb,'CTYPE'//nc//'      ')
!!$    call FTPKYS(unit,'CUNIT'//nc,cunit,comment,sts); call reset_sts(sts,0,verb,'CUNIT'//nc//'      ')
!!$    call FTPKYE(unit,'CDELT'//nc,cdelt,6,comment,sts); call reset_sts(sts,0,verb,'CDELT'//nc//'      ')
!!$    call FTPKYE(unit,'CRVAL'//nc,crval,6,comment,sts); call reset_sts(sts,0,verb,'CRVAL'//nc//'      ')
!!$    call FTPKYE(unit,'CRPIX'//nc,crpix,6,comment,sts); call reset_sts(sts,0,verb,'CRPIX'//nc//'      ')
!!$    

    if (present(nsize)) then 
       comment = ''
       if (present(nsize_com)) comment=nsize_com
       call FTPKYJ(unit,'NAXIS'//nc,nsize,comment,sts)
       call reset_sts(sts,0,verb,'CRVAL'//nc)
    end if
    
    if (present(ctype)) then 
       comment = ''
       if (present(ctype_com)) comment=ctype_com
       call FTPKYS(unit,'CTYPE'//nc,ctype,comment,sts)
       call reset_sts(sts,0,verb,'CTYPE'//nc)
    end if
    
    if (present(cunit)) then
       comment = ''
       if (present(cunit_com)) comment=cunit_com
       call FTPKYS(unit,'CUNIT'//nc,cunit,comment,sts)
       call reset_sts(sts,0,verb,'CUNIT'//nc)
    end if
    
    if (present(bunit)) then
       comment = ''
       if (present(bunit_com)) comment=bunit_com
       call FTPKYS(unit,'BUNIT',bunit,comment,sts)
       call reset_sts(sts,0,verb,'BUNIT'//nc)
    end if
    

    if (present(cdelt)) then
       comment = ''
       if (present(cdelt_com)) comment=cdelt_com
       call FTPKYE(unit,'CDELT'//nc,cdelt,6,comment,sts)
       call reset_sts(sts,0,verb,'CDELT'//nc)
    end if
    
    if (present(crval)) then
       comment = ''
       if (present(crval_com)) comment=crval_com
       call FTPKYE(unit,'CRVAL'//nc,crval,6,comment,sts)
       call reset_sts(sts,0,verb,'CRVAL'//nc)
    end if
    
    if (present(crpix)) then
       comment = ''
       if (present(crpix_com)) comment=crpix_com
       call FTPKYE(unit,'CRPIX'//nc,crpix,6,comment,sts)
       call reset_sts(sts,0,verb,'CRPIX'//nc)
    end if
    
    return    
  end subroutine scale_to_wcs


  subroutine ccfscale_to_wcs(unit, crval, crpix, cdelt, verbose)
    !this subroutine put a RV reference scale in a 2d images 
    !so it acts differently from wcs_to_scale, that reads a 1D image
    integer, intent(in) :: unit
    real, intent(in) :: crval, crpix, cdelt

    integer, intent(in), optional :: verbose
    !character (len=16), intent(out), optional :: ctype, cunit

    integer ::  sts, verb, nsize
    character (len=64) :: comment
    
    !fits with more than 9 dimensions are not supported
    
    sts = 0 ; verb = 0
    
    if (present(verbose)) verb=verbose
    comment = ''
    !call FTPKYJ(unit,'NAXIS1',nsize  ,comment,sts); call reset_sts(sts,0,verb,'CRVAL1      ')
    call FTPKYS(unit,'CTYPE1','km/s    ',comment,sts)
    call reset_sts(sts,0,verb,'CTYPE1')
    
    call FTPKYS(unit,'CUNIT1','[km/s]  ',comment,sts)
    call reset_sts(sts,0,verb,'CUNIT1')
    
    call FTPKYE(unit,'CDELT1',cdelt,4,'CCF steps [km/s]'  ,sts)
    call reset_sts(sts,0,verb,'CDELT1')
    
    call FTPKYE(unit,'CRVAL1',crval,4,'value of ref pixel',sts)
    call reset_sts(sts,0,verb,'CRVAL1')

    call FTPKYE(unit,'CRPIX1',crpix,1,'Reference pixel'   ,sts)
    call reset_sts(sts,0,verb,'CRPIX1')
    
    return    
  end subroutine ccfscale_to_wcs
  

  subroutine fits_open_s1d (unit,file,dpix,verbose)
    !!! This routine open and read a 1D fits file with unknown dimension
    !!! The input variable must be an allocatable REAL*8, even if the
    !!! fits file has been save as a REAL*4
    !!! Subroutine 'reset_sts' has been temporalily turned down
    !!! The file is not closed at the return point
    
    integer, intent(in) :: unit
    character (len = nch_file), intent(in) :: file

    integer, intent(in), optional :: verbose
    
    real (kind=8), dimension(:), intent(out) :: dpix

    integer :: rwmod, blksz, sts,  maxdim, naxis, naxes, bitpix, verb
    character (len = 30) :: errtext
    logical :: anyf

    rwmod = 0 ;  maxdim = 1 ;  sts = 0 ;  verb = 0
    if (present(verbose)) verb=verbose

    !     opens fits 'file'
    if (verb.gt.0) write(verb,'(A45,I3,A3)',advance='no') &
         '   *** FTSOAL: Opening in REAL  mode in unit ',unit,' : '
    if (verb.gt.0) write(verb,*) trim(file)

    call FTOPEN(unit,file,rwmod,blksz,sts) ;  call reset_sts(sts,1,verb)

    !     extracts basic header info from fits
    call FTGIPR(unit, maxdim, bitpix, naxis, naxes, sts) ;  call reset_sts(sts,1,verb)

    if (verb.gt.0) then
       write (verb,'(1X, A18, I9)') '     image depth: ', bitpix
       write (verb,'(1X, A18, I9)') '     axis:        ', naxis
       write (verb,'(1X, A18, I9)') '     points:      ', naxes
    end if

    if (size(dpix,1).ne.naxes) stop

    !    extracts the image array from fits and puts it in 'pix'
    call FTGPVD(unit,0,1,naxes,nullval_d,dpix,anyf,sts) ;  call reset_sts(sts,1,verb)

    if (verb.gt.0.and.anyf) write(verb,*) '   *** FTSOAL: Warning: there are undefined data values (anyf=T)'
    return
  end subroutine fits_open_s1d

  subroutine fits_open_s2d (unit,file,d2pix,verbose)
    !!! This routine open and read a 2D fits file with unknown dimensions
    !!! The input variable must be an allocatable REAL*8, even if the
    !!! fits file has been save as a REAL*4
    !!! The file is not closed at the return point
    
    !!! Pay attention: in the HARPS format, spectra are save in the horizontal
    !!! direction!
    
    integer, intent(in) :: unit
    character (len = nch_file), intent(in) :: file

    integer, intent(in), optional :: verbose
    
    real (kind=8), dimension(:,:), intent(out) :: d2pix

    integer :: rwmod, blksz, sts,  maxdim, naxis, bitpix, verb
    integer, dimension(2) :: naxes
    character (len = 30) :: errtext
    logical :: anyf

    rwmod = 0 ;  maxdim = 2 ;  sts = 0 ;  verb = 0
    if (present(verbose)) verb=verbose

    !     opens fits 'file'
    if (verb.gt.0) write(verb,'(A45,I3,A3,A)') &
         '   *** FTSOAL: Opening in REAL  mode in unit ',unit,' : ', trim(file)
    
    call FTOPEN(unit,file,rwmod,blksz,sts); call reset_sts(sts,1,verb)

    !     extracts basic header info from fits
    call FTGIPR(unit, maxdim, bitpix, naxis, naxes, sts);  call reset_sts(sts,1,verb)

    if (verb.gt.0) then
       write (verb,'(1X, A18, I9)') '     image depth: ', bitpix
       write (verb,'(1X, A18, I9)') '     axis:        ', naxis
       write (verb,'(1X, A18, I9)') '     nrows:       ', naxes(2)
       write (verb,'(1X, A18, I9)') '     ncols:       ', naxes(1)
    end if
    
    if (size(d2pix,1).ne.naxes(1) .and. size(d2pix,2).ne.naxes(2)) stop

    !    extracts the image array from fits and puts it in 'pix'
    call FTG2DD(unit,0,nullval_d,naxes(1),naxes(1),naxes(2),d2pix,anyf,sts); call reset_sts(sts,1,verb)
    
    if (verb.gt.0.and.anyf) write(verb,*) '   *** FTSOAL: Warning: there are undefined data values (anyf=T)'
    return
  end subroutine fits_open_s2d
  

  subroutine fits_write_1d(unit,filename,dpix,verbose)
    !     open a fits file $file and write all the pixel values from the 
    !     array pix(), as reals
    integer, intent(in) ::unit
    real (kind = 8), dimension(:), intent(in) :: dpix
    character (len=nch_file), intent(in) :: filename
    integer, intent(in), optional :: verbose

    integer :: sts, naxes, verb
    character (len= 20) :: hdutype

    sts = 0
    naxes = size(dpix,1)

    verb = 0 ;  if (present(verbose)) verb=verbose
    
    if (verb.gt.0) write(verb,'(A45,I3,A3,A)') &
         '   *** FTSWRI: Writing in DBLE  mode in unit ',unit,' : ' , trim(filename)

    call fits_delete(filename,sts)    ; call reset_sts(sts,0,verb,'FTSDEL')
    call ftinit(unit,filename,1,sts)  ; call reset_sts(sts,1,verb,'FTINIT')
    call ftiimg(unit,-64,1,naxes,sts) ; call reset_sts(sts,1,verb,'FTIIMG')
    call FTPPRD(unit,0,1,naxes,dpix,sts) ;  call reset_sts(sts,1,verb,'FTPPRD')
    call FTMAHD(unit,1,hdutype,sts) ;  call reset_sts(sts,1,verb,'FTMAHD')
    ! return to the first HDU to give the possibility to add more header keys
    
    return
  end subroutine fits_write_1d

  subroutine fits_write_2d(unit,filename,d2pix,verbose)
    !     open a fits file $file and write all the pixel values from the 
    !     array pix(), as reals
    integer, intent(in) ::unit
    real (kind = 8), dimension(:,:), intent(in) :: d2pix
    character (len=nch_file), intent(in) :: filename
    integer, intent(in), optional :: verbose

    integer :: sts, naxes(2),verb
    character (len= 20) :: hdutype

    sts = 0
    naxes(1) = size(d2pix,1) ; naxes(2) = size(d2pix,2)

    verb = 0 ;  if (present(verbose)) verb=verbose
    
    if (verb.gt.0) write(verb,'(A45,I3,A3,A)') &
         '   *** FTSWRI: Writing in DBLE  mode in unit ',unit,' : ' , trim(filename)

    call fits_delete(filename,sts)    ; call reset_sts(sts,0,verb,'FTSDEL')
    call ftinit(unit,filename,1,sts)  ; call reset_sts(sts,1,verb,'FTINIT')
    call ftiimg(unit,-64,2,naxes,sts) ; call reset_sts(sts,1,verb,'FTIIMG')
    call ftp2dd(unit,0,naxes(1),naxes(1),naxes(2),d2pix,sts) ;  call reset_sts(sts,1,verb,'FPT2DD')
    call FTMAHD(unit,1,hdutype,sts) ;  call reset_sts(sts,1,verb,'FTMAHD')
    ! return to the first HDU to give the possibility to add more header keys
    
    return
  end subroutine fits_write_2d

  subroutine fits_close(unit)
    integer, intent(in) :: unit
    integer :: sts

    sts = 0
    call FTCLOS(unit, sts)
    call reset_sts(sts,1,0)

    return
  end subroutine fits_close

  subroutine fits_check(filename,sts_out)
    !  A simple little routine to check the existence of a file and delete it if
    ! the fits file contains errors
    ! sts_out = 0 -> file doesn't exist
    ! sts_out = 1 -> file already exists
    integer :: sts,sts_out,unit,blocksize
    character ( len= nch_file) :: filename
    sts = 0
    call get_lun(unit) ! Get an unused Logical Unit Number to use to open the FITS file
    call ftopen (unit,filename,0,blocksize,sts) ! Try to open the file, to see if it exists
    sts_out = 0
    if (sts .eq. 0) then
       sts_out = 1 ! file was opened (it exist)
    else if (sts .ne. 103) then
       sts=0 ; call ftcmsg ; call ftdelt(unit,sts)
       !         there was some other error opening the file; delete the file anyway
    end if
    call fits_close(unit) !  Free the unit number for later reuse
    call free_lun(unit)

    return
  end subroutine fits_check


!!! delete file routine: check if the file already exist... 
  subroutine fits_delete(filename,sts)
    !  A simple little routine to delete a FITS file
    !  from the FITSIO Cookbook (readapted to F90) "deletefile"
    integer :: sts,unit,blocksize
    character ( len= nch_file) :: filename
    
    if (sts .gt. 0) return  ! Simply return if sts is greater than zero
    call ftgiou(unit,sts) ! Get an unused Logical Unit Number to use to open the FITS file
    call ftopen(unit,filename,1,blocksize,sts) ! Try to open the file, to see if it exists

    if (sts .eq. 0) then
       call ftdelt(unit,sts) ! file was opened;  so now delete it 
    else if (sts .eq. 103) then
       sts=0 ; call ftcmsg ! file doesn't exist, so just reset status to zero and clear errors
    else
       sts=0 ; call ftcmsg ; call ftdelt(unit,sts)
       !         there was some other error opening the file; delete the file anyway
    end if

    call ftfiou(unit, sts) !  Free the unit number for later reuse

  end subroutine fits_delete

  subroutine reset_sts(sts,do,verbose,key)
    integer, intent(inout) :: sts
    integer, intent(in) :: do, verbose
    character (len = *), intent(in), optional :: key
    
    character (len = 80) :: errtext
    
    if (sts.ne.0) then
       call FTGERR(sts,errtext)
       
       if (verbose.gt.0) then
!!$          key_t='UNDEFINED'
!!$          if (present(key)) key_t=key
!!$          write(verbose,'(A25,I5,A21,A)') & 
!!$               '   ****** CFITSIO error: ', sts,' while operating on: ', trim(key_t)
          if (present(key)) then
             write(verbose,'(A25,I5,A21,A)') & 
                  '   ****** CFITSIO error: ', sts,' while operating on: ', trim(key)
             else
                write(verbose,'(A25,I5,A)') & 
                  '   ****** CFITSIO error: ', sts,' while operating on: UNDEFINED'
             end if
          write(verbose,*) '  ****** CFITSIO error: ',trim(errtext)
       endif

       sts=0 
       if (do.gt.0) then 
          write(*,*) '  ****** CFITSIO error: ',trim(errtext)
          stop
       end if
       
    end if
    return
    
  end subroutine reset_sts


end module fits
