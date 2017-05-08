Program harps_input2pams
  use Common
  use Fits
  use Instrument_HARPS
  use CCF
  use Chebyshev_Fit
  use Poly_Fit
  use Gaussian_Fit

  implicit none

  character (len=*), parameter :: &
    code_path = '/Users/malavolta/Astro/CODE/CCFpams/', &
    archive_harps = '/Users/malavolta/data/HARPS/'

  !!$ hard-coded files to easily run the script
  character (len=*), parameter :: &
          calib_fdir='mask_calib/', sccf_dir ='../mask_ccf/', &
          calib_mask1 = 'TGdirect_666_direct_calib_cheb.dat', &
          calib_mask2 = 'LOGG_ewfind_33222_ewfind_logg4_calib_cheb.dat', &
          line_list = 'Sunsky_linelist.dat'

  character (len=*), parameter :: coelho2005_sun= trim(calib_fdir) // 'Coelho2005_5750_45_p00p00.ms.fits.gz'

  integer :: n_list, n_com, date_opt
  integer :: n_count

  logical :: verbose
  character (len=nch_file) :: option, option_in
  character (len=nch_file) :: file, outname_sum,  mask_out, file_blaze, object_list, &
       object_file, file_night, file_rad, fits_name, templ_name, archive_files

  character (len=2) :: mask_sel
  character (len=1) :: fiber_sel
  real (PR) :: rvc_fin, rv_drift, berv, airmass

  integer :: star_iostat,read_iostat,cheb_iostat,sts,sts_out     !these are working varibles
  integer :: fits_lun, input_lun, list_lun, chb_lun,&
       out1_lun, out2_lun, out3_lun !these are logical unit
  integer :: output_lun, extended_lun
  integer :: i, ii, n, nl, no, ns, nx, ny, nt ! these are counting variables
  integer :: nh_sel, check

  character (len=nch_object) :: object

  ! fits images informations
  real (PR), dimension(n_pix,n_order) :: d2pix,w2pix,p2pix,blaze
  real (PR), dimension(n_order) :: wl_min, wl_max

  real (PR), allocatable, dimension(:) ::  cont_val
  real (PR), allocatable, dimension(:,:) ::  dwpix

  integer :: wl_is, wl_ie, wl_ii
  real (kind=8), dimension(n_order) :: int_flux, int_wlav, int_tmpl
  real (kind=8), dimension(n_order) :: wlr_min, wlr_max, wlr_rng

  real (kind=8), dimension(:,:), allocatable :: cont_fits,cont_tmp
  real (PR) :: m_factor

  integer, parameter :: chb_order = 6
  integer :: ADC_size
  real (kind=8)   , dimension(n_order) :: ADC_x, ADC_y1, ADC_z
  real (kind=8)   , dimension(chb_order+1) :: chb_coeff1, chb_coeff2



  real :: ccd_sigdet


!!! Gaussian fitting
  real (kind=PR) :: ampl_ccf, cntr_ccf, sgm_ccf, back_ccf

  !! mask parameters
  !! CCF
  integer :: x_axis, y_axis
  real (kind=PR) :: mask_wb
  real (kind=PR), dimension(:,:), allocatable:: ccf_sum
  real (kind=PR), dimension(:,:), allocatable:: ccf_sum_weighted, ccf_sum_total
  real (kind=PR), dimension(rv_nbin) :: ccf_x, ccf_y, ccf_w, ccf_z, ccf_s

  !input linelist
  integer :: lmask_n
  real (PR), dimension(:), allocatable :: lmask_wl
  real         , dimension(:), allocatable :: lmask_ep,lmask_lgf,lmask_sew,lmask_elm
  real (PR), dimension(:), allocatable::  lmask_wg1, lmask_wd
  logical, allocatable, dimension(:,:) :: lmask_sel

  character (len=4) :: chr_rvnbin

  real (PR) ::  cntr
  real :: elm, sew, ep, lgf

!!! Calibrations variables

  integer :: tmp_int, tmp_int2
  real (PR) :: tmp_real

  !! Chebishev grid result
  integer ::  tot_niter
  integer ,parameter :: logg_iter = 14
  real   , dimension(logg_iter), parameter :: &
       logg_values = (/3.600,3.700,3.800,3.900,4.000,4.100,4.200,4.300,4.400,4.500,4.600,4.700,4.800,4.900/)

  logical, dimension(logg_iter) :: cal2_logg_sel



  integer, parameter :: n_cal2 = 1 , n_ord2 = 2
  real (kind=PR), dimension(n_ord2,logg_iter) :: cheb2_xmin, cheb2_xmax, cheb2_xavg, cheb2_xrng
  real (kind=PR) :: cheb2_x_in,cheb2_y_in
  integer ::  cheb2_xord, cheb2_yord, cheb2_polyord
  real (kind=PR), allocatable, dimension(:,:,:) :: cheb2_coeff_g


  integer, parameter :: n_cal3 = 3 , n_ord3 = 3
  real (kind=PR), dimension(n_ord3) :: cheb3_xmin, cheb3_xmax, cheb3_xavg, cheb3_xrng
  real (kind=PR) :: cheb3_x_in,cheb3_y_in,cheb3_z_in
  integer ::  cheb3_xord, cheb3_yord, cheb3_zord
  real (kind=PR), allocatable, dimension(:,:,:) :: cheb3_coeff_t, cheb3_coeff_m


  real (kind=PR) ::   gcoeff_x0, gcoeff_y0, gcoeff_y1, gcoeff_y2

  real (kind=PR), dimension(logg_iter):: logg_xf,  logg_yf, logg_wf
  real (kind=PR), dimension(:), allocatable :: poly_coeff

  integer, parameter :: n_cal = n_cal3 + n_cal2
  character (len=nch_file), dimension(n_cal) :: cal_name


  real (PR), dimension(n_cal) :: area_out, sigma_out, center_out, contrast_out, back_out

  real (PR) :: teff_val,gfeh_val,logg_val, bjd_obs, snr_pip, teff_err,gfeh_err,logg_err
  real :: snr_tmp
  real (PR) :: snr_tot, bjd_avg

  real  :: ccf_pcval

  integer :: nl_i, lun_iostat
  real (PR) :: lmask_wg1_i, cont_val_i
  logical :: it_exists
  character (len=nch_file) :: file_name


  real(PR) :: dummy_dbl1

  real (PR) :: cont_ccf_i
  real (PR), dimension(:), allocatable :: cont_ccf

  !!call getarg(1,line_list)
  !!call getarg(2,calib_mask1)
  !!call getarg(3,calib_mask2)

  n_com = iargc()

  if (n_com.lt.1) then
     write(*,*) 'Input 1: object list file'
     write(*,*) 'Input 2: output prefix (usually the object name)'
     write(*,*) 'Input 3: (optional): archive directory '
     write(*,*) 'Input 4: (optional): 0 = do not use date as prefix '
     write(*,*) 'Input 4  to be used if all the files are in the same directory (as Yabi reprocessed files are) '
     write(*,*) 'Input 5: (optional): verbosity on '
     stop
  end if


  call getarg(1,option_in);  object_list=trim(option_in)
  call getarg(2,option_in);  object=trim(option_in)
  write(*,*)  object_list

  archive_files=archive_harps

  !! check wheter the data are store in a local folder or in the default HARPS-N folder (default)
  if (n_com .ge. 3) then
    call getarg(3,option_in) !;  read (option,*) option_in
    archive_files=trim(option_in)
    write(*,*) archive_files
  end if

  date_opt = 1
  if (n_com .ge. 4) then
    call getarg(4,option_in); read (option_in,*) date_opt
    write(*,*) 'Prepend date: ', date_opt
  end if

  verbose = .False.
  if (n_com .ge. 5) then
    call getarg(5,option_in)!;  read (option,*) option_in
    if (option_in .eq. 'verbose') verbose = .True.
  end if


  !! FIRST PART: Building the CCFs...
  allocate(dwpix(n_pix,n_order))

  call get_lun(list_lun)
  call get_lun(fits_lun)
  call get_lun(chb_lun)

  call get_lun(out1_lun)
  call get_lun(out2_lun)
  call get_lun(out3_lun)
  call get_lun(input_lun)
  call get_lun(output_lun)
  call get_lun(extended_lun)

  do ii=1,rv_nbin
     ccf_x(ii) =  rv_bin*(ii - int(rv_nbin/2))
  end do
  ccf_pcval = rv_bin*(1 - int(rv_nbin/2))

  !!!! write(*,*) 'Linelist file: ', trim(calib_fdir)//trim(line_list)
  open(list_lun,file=code_path//calib_fdir//trim(line_list),status='old')
  read(list_lun,*) lmask_n

  !!!! write(*,*) 'Linelist file: ', trim(calib_fdir)//trim(line_list), '     N. of lines:',lmask_n

  allocate(lmask_wl(lmask_n),lmask_ep(lmask_n),lmask_lgf(lmask_n),lmask_sew(lmask_n),lmask_elm(lmask_n))
  allocate(lmask_wg1(lmask_n),lmask_wd(lmask_n))
  allocate(ccf_sum(rv_nbin,lmask_n))
  allocate(ccf_sum_weighted(rv_nbin,lmask_n))
  allocate(ccf_sum_total(rv_nbin,lmask_n))

  allocate(cont_val(lmask_n))

  lmask_ep = 0._4  ; lmask_lgf= 0._4  ; lmask_sew= 0._4  ; lmask_elm= 0._4
  lmask_wl = 0._PR ; lmask_wg1 = 0._PR

  star_iostat=0
  do while (star_iostat.eq.0)
     read (list_lun,*,iostat=star_iostat)   nl, elm, sew, ep, lgf, cntr
     lmask_wl (nl) = cntr ;  lmask_elm(nl) = elm ;  lmask_sew(nl) = sew
     lmask_ep (nl) = ep   ;  lmask_lgf(nl) = lgf
  end do
  close(list_lun)

  call int2chr4(rv_nbin,chr_rvnbin)
  mask_out = 'inp_' // chr_rvnbin

  where(lmask_elm.gt.0) lmask_wg1 = 1._PR

  lmask_wd = lmask_wl/cc * 0.75_PR !mask_width

  templ_name = code_path//calib_fdir//'Coelho2005_5750_45_p00p00.ms.fits.gz'

  call fits_size_2d (templ_name,x_axis,y_axis)
  !!!! write(*,*) 'TEMPLATE: ', trim(templ_name), x_axis,y_axis

  !! Original template file is stored in cont_tmp; informations that we need are saved in cont_fits
  !! cont_fits(:,1) = wavelength scale of the template
  !! cont_fits(:,2) = size of each unit pixel in wavelength space
  !! cont_fits(:,3) = flux normalized template spectra
  !! cont_fits(:,4) = flux template spectra (with absorption lines)
  !! cont_fits(:,5) = flux continuum
  !! cont_fits(:,6) = smoothed spectra (not normalized)

  allocate(cont_tmp(x_axis,y_axis), cont_fits(x_axis,6))
  cont_fits = 0.00000
  !! wcs_to_scale convert the wcs informtion in an array of wavelength values
  call fits_open_s2d(fits_lun,templ_name,cont_tmp)
  call wcs_to_scale (fits_lun, 1, cont_fits(:,1))
  call fits_close   (fits_lun)

  cont_fits(:,3) = cont_tmp(:,1)
  cont_fits(:,4) = cont_tmp(:,2)
  cont_fits(:,5) = cont_fits(:,4)/cont_fits(:,3)

  deallocate (cont_tmp)

  cont_fits(2:x_axis-1,2) = (cont_fits(3:x_axis,1) - cont_fits(1:x_axis-2,1))/2
  cont_fits(1,2) = cont_fits(2,2);  cont_fits(x_axis,2) = cont_fits(x_axis-1,2)
  cont_fits(:,6) = cont_fits(:,4)
  !!!! write(*,*) 'template reading complete'


  wl_ii = minloc(abs(cont_fits(:,1)-5500._PR),dim=1)


  cont_fits(:,4) = cont_fits(:,4) / cont_fits(wl_ii,2)
  cont_fits(:,5) = cont_fits(:,5) / cont_fits(wl_ii,2)
  cont_fits(:,6) = cont_fits(:,6) / cont_fits(wl_ii,2)


  cont_val = 0._PR
  do nl=1,lmask_n
     wl_ii = minloc(abs(cont_fits(:,1)-lmask_wl(nl)),dim=1)
     cont_val(nl) = cont_fits(wl_ii,5)/cont_fits(wl_ii,2)
  end do



  !! READING THE MASKS


  !! Calibration files are read: Direct Teff+[Fe/H] file, log(g)_EW file, and Teff+[Fe/H]+logg_EW file
  !! This one is the direct calibration for Teff+[Fe/H], called n_cal3 because inerithered from
  !! an old version of the program


  write(*,*) 'opening file:', code_path//calib_fdir//trim(calib_mask1)
  open(chb_lun,file=code_path//calib_fdir//trim(calib_mask1),status='old')
  read(chb_lun,*) tmp_int
  if (tmp_int.ne.n_cal3) then
     !!!! write(*,*) 'WRONG FILE!!!'
     stop
  end if

  do ii=1,n_cal3
     read(chb_lun,*) cal_name(ii)
  write(*,*)  trim(cal_name(ii))
  end do

  read(chb_lun,*) cheb3_xord, cheb3_yord, cheb3_zord
  allocate(cheb3_coeff_t(cheb3_xord+1,cheb3_yord+1,cheb3_zord+1), &
       cheb3_coeff_m(cheb3_xord+1,cheb3_yord+1,cheb3_zord+1) )

  cheb3_xmin=0.000d0; cheb3_xmax=0.000d0; cheb3_xavg=0.000d0;  cheb3_xrng=0.000d0
  cheb3_coeff_t=0.000d0 ; cheb3_coeff_m=0.000d0

  do ii=1,n_cal3
     read(chb_lun,*) cheb3_xmin(ii), cheb3_xmax(ii), cheb3_xavg(ii), cheb3_xrng(ii)
  end do
  do nx=1,cheb3_xord+1
     do ny=1,cheb3_yord+1
        read(chb_lun,*) cheb3_coeff_t(nx,ny,:)
     end do
  end do
  do nx=1,cheb3_xord+1
     do ny=1,cheb3_yord+1
        read(chb_lun,*) cheb3_coeff_m(nx,ny,:)
     end do
  end do
  close(chb_lun)


  !! This one is the direct calibration for logg_EWs

  write(*,*) 'opening file:', code_path//calib_fdir//trim(calib_mask2)
  open(chb_lun,file=code_path//calib_fdir//trim(calib_mask2),status='old')

  read(chb_lun,*) tmp_int, tmp_int2, tot_niter
  if (tmp_int.ne.n_cal2) then
     write(*,*) 'only ',n_cal3,'-dimensional functions are supported'
     stop
  end if

  if (tot_niter.ne.logg_iter) then
     write(*,*) 'ERROR in LOGG-derived values'
     stop
  end if

  do i=n_cal3+1,n_cal
     read(chb_lun,*) cal_name(i)
     write(*,*) trim(cal_name(i))
  end do
  read(chb_lun,*) cheb2_xord, cheb2_yord, cheb2_polyord
!!! write(*,*) cheb2_xord, cheb2_yord, cheb2_polyord

  allocate(cheb2_coeff_g(cheb2_xord+1,cheb2_yord+1,tot_niter))

  cheb2_xmin=0.000d0; cheb2_xmax=0.000d0; cheb2_xavg=0.000d0;  cheb2_xrng=0.000d0
  cheb2_coeff_g=0.000d0

  cal2_logg_sel = .false.

  cheb_iostat = 0
  nt = 0
  !!do nt=1,tot_niter
  do while (cheb_iostat.eq.0)
     read(chb_lun,*,iostat=cheb_iostat) tmp_int, tmp_real
     if (cheb_iostat.eq.0) then
        cal2_logg_sel(tmp_int)=.true.
     else
        exit
     end if
     nt = tmp_int
     do ii=1,n_ord2
        read(chb_lun,*) cheb2_xmin(ii,nt), cheb2_xmax(ii,nt), cheb2_xavg(ii,nt), cheb2_xrng(ii,nt)
     end do

     do nx=1,cheb2_xord+1
        read(chb_lun,*) cheb2_coeff_g(nx,:,nt)
     end do
     if (nt.eq.tot_niter) exit
  end do


  read (chb_lun,*)gcoeff_x0 , gcoeff_y0 , gcoeff_y1, gcoeff_y2
  close(chb_lun)

  allocate (poly_coeff(cheb2_polyord+1))




  !! select   lmask_sel
  allocate(lmask_sel(lmask_n,n_cal))  ;  lmask_sel=.false.
  do i=1,n_cal

     !! The input selection is store, useful when line parameter criteria instead of direct input are used
     write(*,*) 'OPENING linelist file: ', trim(cal_name(i))//'_linelist.dat'
     open(list_lun,file=code_path//calib_fdir//trim(cal_name(i))//'_linelist.dat',status='old')

     !! BUT "nl" must remain untouched (it identifies the correct row in the CCF file)
     star_iostat=0
     do while (star_iostat.eq.0)
        read (list_lun,*,iostat=star_iostat) n, nl, elm, cntr, ep, sew
        if (star_iostat.ne.0) exit
        if (abs(lmask_wl(nl) - cntr) .lt. 0.500d0) lmask_sel(nl,i) = .true.

     end do
     close(list_lun)

  end do

  bjd_avg = 0.0d0
  write(*,*)  object_list
  open(input_lun,file=trim(object_list),status='old')

  !! CCFs and configuration files for each exsposure are created individually
  !! This sectionis a mere copy&past from harps_list2ccfs
  call system('mkdir -p output_dir')
  object_file =  'output_dir/'// trim(object) // '_outcal.dat'
  open(output_lun,file=trim(object_file),status='new')

  object_file =  'output_dir/'// trim(object) // '_outcal_extended.dat'
  open(extended_lun,file=trim(object_file),status='new')

  n_list = 0

  nh_sel = 0
  n_count = 0
  read_iostat = 0


  ccf_sum_total = 0._PR
  snr_tot = 0._PR

  allocate(cont_ccf(lmask_n))

  !! New addition to better perform flux correction

  fiber_sel = 'A'
  do while (read_iostat.eq.0)
     read(input_lun,*,iostat=read_iostat) file_night, file_rad, mask_sel
     if (read_iostat.ne.0) exit
     write(*,*)
     write(*,*) 'Doing spectrum ', trim(file_rad)
     n_list = n_list + 1

     if (date_opt .eq. 0) file_night=''

     lmask_wg1 = 0._PR
     ccf_sum = 0.0

     outname_sum = sccf_dir//trim(object)//'/'//trim(file_rad)//'_sccf.fits.gz'
     write(*,*) 'OUTNAME_sum ', trim(outname_sum)
     call fits_check(outname_sum,check)

     if (check.eq.0) then
        !File has not been processed before, going now

        call system('mkdir -p '//sccf_dir//trim(object))

        !reading the fits file and
        file = trim(archive_files)//'/'//trim(file_night)//'/'//trim(file_rad)//'_e2ds_'//trim(fiber_sel)//'.fits'
        call fits_check(file,sts_out)
        if (sts_out.eq.0) then
           write(*,*) '     FILE '//trim(file)//' does not exist, no output file created for this exposure'
           cycle
        end if
        call fits_open_s2d(fits_lun,file,d2pix)
        call key2wave     (fits_lun,w2pix)
        call key2dwave    (fits_lun,dwpix)
        call fits_close(fits_lun)

        !we retrieve the blaze correction fits
        file = trim(archive_files)//'/'//trim(file_night)//'/'//trim(file_rad)//'_ccf_'// &
             trim(mask_sel)//'_'//trim(fiber_sel)//'.fits'

        call fits_check(file,sts_out)
        if (sts_out.eq.0) then
           write(*,*) '     FILE '//trim(file)//' does not exist, no output file created for this exposure'
           cycle
        end if

        call fits_actv(fits_lun,file)
        call singleheader_AIRMASS(fits_lun,airmass ,sts_out)
        call singleheader_BJD(fits_lun,bjd_obs ,sts_out)
        call singleheader_RVC    (fits_lun,rvc_fin ,sts_out)
        call singleheader_RVDRIFT(fits_lun,rv_drift,sts_out)
        call singleheader_BERV   (fits_lun,berv    ,sts_out)
        call singleheader_BLAZE  (fits_lun,file_blaze,sts_out)
        call singleheader_CCD_SIGDET(fits_lun,ccd_sigdet,sts_out)
        call singleheader_DRS_SNR(fits_lun,60,snr_tmp,sts)
        !! add RVs
        call fits_close(fits_lun)
        snr_pip = real(snr_tmp)

        !correcting for earth barycentric radial velocity, drift and measured radial velocity
        !! UPDATE: We are now using directly the drift-uncorrected RV
        w2pix = (1._PR)*((berv-rv_drift-rvc_fin)/cc+1._PR)*w2pix

        file = trim(archive_files)//'/'//trim(file_night)//'/'//trim(file_blaze)
        call fits_check(file,sts_out)
        if (sts_out.eq.0) then
           write(*,*) '     FILE '//trim(file)//' does not exist, no output file created for this exposure'
           cycle
        end if
        call fits_open_s2d(fits_lun,file,blaze)
        call fits_close(fits_lun)
        ! blaze = 1._PR

        !! All the files checked exist, we can go further
        !! open(out1_lun,file= trim(object)//'/'//trim(file_rad)//'_pams.dat',status='new')
        !! open(out3_lun,file= trim(object)//'/'//trim(file_rad)//'_pams_check.dat',status='new')


        nh_sel = nh_sel + 1
        d2pix = d2pix / blaze
        p2pix = d2pix / (dwpix / 0.01_PR)

        !we have already corrected the wavelength scale for the radial velocity + baryc velocity
        !we only have to take in account the range in RV over which the CCF is computed
        !rv_pix * dw (min_val) -> wl at which I have to start to compute the CCF
        wl_min=w2pix(    1,:)*(1._PR)/(1._PR-real((rv_nbin)/2)*rv_bin/cc)
        wl_max=w2pix(n_pix,:)*(1._PR)/(1._PR+real((rv_nbin)/2)*rv_bin/cc)

        ! ccf_x -> single line CCF for single order
        ! ccf_s -> single line CCF for all the orders in which is present


        do nl=1,lmask_n
           if (lmask_elm(nl).le.0) cycle
           ccf_s = 0._PR
           mask_wb = 0._PR
           ns = 0
           do no = 1,n_order

              !remember that we moved to zero-centered radial velocity!
              if (lmask_wl(nl).lt.wl_min(no)) cycle
              if (lmask_wl(nl).gt.wl_max(no)) cycle
              ccf_y = 0._PR
              ii = minloc(abs(w2pix(:,no)-lmask_wl(nl)),1)
!!$           mask_wb = lmask_wg(nl) * blaze(ii,no)

              call ccf_singleline(lmask_wl(nl),lmask_wd(nl),&
                   d2pix(:,no),w2pix(:,no),dwpix(:,no),ccf_x,ccf_y)
              ccf_s = ccf_s + ccf_y
              ns = ns + 1

!!$           write(out2_lun,'(4I6,4E17.9)') nh_sel, nl, no, ns, mask_wb, lmask_wd(nl), blaze(ii,no), dwpix(ii,no)
           end do
           if (ns.le.0) cycle

           !do i=1,rv_nbin
           !   write(out_lun(nm),'(5I5,E12.4,E17.9)') nm,nl,no,nh,i,ccf_x(i),ccf_y(i)*mask_wb(nl)
           !end do

!!$        ccf_sum(:,nl) = ccf_s/real(ns)
           ccf_sum(:,nl) = ccf_sum(:,nl) + ccf_s/real(ns)

        end do

        cont_ccf = (sum(ccf_sum(1:40,:),dim=1) + sum(ccf_sum(141:180,:),dim=1)/40.)/80.

        !do i=1,rv_nbin
        !   write(out_lun(nm),'(5I5,E12.4,E17.9)') nm,0,no,nh,i,ccf_x(i), ccf_y(i)
        !end do
        !end do

        if (verbose) then
          call system ('mkdir -p '//sccf_dir//trim(object)//'/')

          call fits_write_2d(fits_lun,outname_sum,ccf_sum)
          call ccfscale_to_wcs(fits_lun, real(ccf_x(1)),real( 1.00), real(rv_bin))
          call fits_close(fits_lun)
        end if
!!$     close(out2_lun)



        ! We determine the differential refraction point to be interpolated
        ADC_x = 0.000d0
        ADC_y1 = 0.000d0
        ADC_z = 0.000d0
        ADC_size = 0


        wlr_min(:) = w2pix(n_pix/8,:)
        wlr_max(:) = w2pix(n_pix*7/8,:)
        wlr_rng(:) = (w2pix(n_pix*7/8,:)-w2pix(n_pix/8,:))

        int_flux(:) = sum(d2pix(n_pix/8:n_pix*7/8,:),dim=1)/wlr_rng(:)
        int_wlav(:) = (w2pix(n_pix*7/8,:)+w2pix(n_pix/8,:))/2.000d0

        do no=1,n_order
           wl_is = minloc(abs(cont_fits(:,1)-wlr_min(no)),dim=1)
           wl_ie = minloc(abs(cont_fits(:,1)-wlr_max(no)),dim=1)
           int_tmpl(no) = sum(cont_fits(wl_is:wl_ie,6),dim=1)/wlr_rng(no)
        end do

        do no=1,n_order
           if (order_mask(no).eqv. .false.) cycle
           adc_size = adc_size + 1

           adc_x(adc_size) = (int_wlav(no)-4950.000d0)/ 3300.000d0
           adc_y1(adc_size) = int_tmpl(no)/int_flux(no)
           adc_z(adc_size) = 1.0000d0

           if (no.eq.1       ) adc_x(1)=wlr_min(no)

        end do

        !! The ADC function is normalized to the reference wavelength value
        !! anf the logarithmic value is taken

        adc_y1(1:adc_size) =  log10( adc_y1(1:adc_size) * (int_flux(no_ref)/int_tmpl(no_ref)) )
        m_factor = int_flux(no_ref) / int_tmpl(no_ref)

        call cheby_fit(ADC_x(1:ADC_size),ADC_y1(1:ADC_size),ADC_z(1:ADC_size), &
             chb_order,chb_coeff1,verbose=6)

        do nl=1,lmask_n
           lmask_wg1(nl) = 10.000d0**(get_cheby( (lmask_wl(nl)-4950.000d0)/3300.000d0, chb_coeff1))
        end do


        if (verbose) then
          open(out2_lun,file=sccf_dir//trim(object)//'/'//trim(file_rad)//'_cpm.dat',status='new')
          do no=1,n_order
            write(out2_lun,*) no, wlr_min(no), wlr_max(no), wlr_rng(no), int_wlav(no), int_flux(no), &
              int_tmpl(no)/int_flux(no)*m_factor, &
              10.000d0**(get_cheby( (int_wlav(no)-4950.000d0)/3300.000d0, chb_coeff1)), &
              10.000d0**(get_cheby( (int_wlav(no)-4950.000d0)/3300.000d0, chb_coeff2))
          end do
          close(out2_lun)

          open(out3_lun,file=sccf_dir//trim(object)//'/'//trim(file_rad)//'_pams_check.dat',status='new')
          do ii=0,330
            tmp_real =  3700.0d0 + real(ii)*10.0d0
            write(out3_lun,*) tmp_real, &
              10.000d0**(get_cheby( (tmp_real-4950.000d0)/3300.000d0, chb_coeff1)), &
              10.000d0**(get_cheby( (tmp_real-4950.000d0)/3300.000d0, chb_coeff2))
          end do
          close(out3_lun)


          open(out1_lun,file=sccf_dir//trim(object)//'/'//trim(file_rad)//'_pams.dat',status='new')
          write(out1_lun,*) m_factor, trim(coelho2005_sun)
          do nl=1,lmask_n
            write(out1_lun,'(I5,4E17.9)') nl, lmask_wg1(nl), 0.0d0, cont_val(nl), cont_ccf(nl)
          end do
          close(out1_lun)

        end if





     else

!!$
        !!file_name = '../mask_ccf/'//trim(object)//'/'//trim(file_rad)//'_pams.dat'
        file_name = sccf_dir//trim(object)//'/'//trim(file_rad)//'_pams.dat'
        !! Check if the required file exists. If not, this exposure is skipped
        inquire(file=trim(file_name),exist=it_exists)
        if (it_exists.eqv. .false.) cycle

        !! fits_name is the fits file with all the CCFs from a single exposure
        !!fits_name ='../mask_ccf/'//trim(object)//'/'//trim(file_rad)//'_sccf.fits'
        fits_name =sccf_dir//trim(object)//'/'//trim(file_rad)//'_sccf.fits.gz'
        !! Check if the required file exists. If not, this exposure is skipped
        write(*,*) trim(fits_name)
        call fits_check(fits_name,sts_out) ; if (sts_out.eq.0) cycle
        ccf_sum = 0._PR

        !! The fits file must be compliant with the rv_nbin and lmask_n parameters,
        call fits_size_2d (fits_name,x_axis,y_axis)
        if (x_axis.ne.rv_nbin .or. y_axis.ne. lmask_n) then
           write(*,*) 'MISMATCH between the WL line list and the size of the fits file'
           write(*,*) ' Probably fits files were generated using another list'
           stop
        end if
        call fits_open_s2d(fits_lun,fits_name,ccf_sum)
        call fits_close   (fits_lun)

        !!$
!!$     outname_sum = trim(object)//'/'//trim(file_rad)//'_sccf_OUT_'//trim(mask_out)//'.fits'
!!$
!!$     call fits_write_2d(fits_lun,outname_sum,ccf_sum)
!!$     call ccfscale_to_wcs(fits_lun, real(ccf_x(1)),real( 1.00), real(rv_bin))
!!$     call fits_close(fits_lun)
!!$

        !! LMASK READING
        open(list_lun,file=trim(file_name),status='old')
        read(list_lun,*) m_factor
        lun_iostat = 0
        do while (lun_iostat.eq.0)
           read(list_lun,*,iostat=lun_iostat) nl_i, lmask_wg1_i,dummy_dbl1, cont_val_i, cont_ccf_i
           if (lun_iostat.ne.0) exit
           lmask_wg1(nl_i) = lmask_wg1_i
           cont_ccf(nl_i) = cont_ccf_i
        end do
        close(list_lun)

     end if

     !lmask_wg = 1.0
     where(lmask_elm.le.0)
        ccf_sum(:,nl) = 0._PR
     end where

     do nl=1,lmask_n
       ccf_sum_weighted(:,nl) = ccf_sum(:,nl)*lmask_wg1(nl)
       !write(*,*) nl, cont_ccf(nl), (cont_ccf(nl)/cont_ccf(cont_reforder)), &
      !  cont_ref(nl), cont_ref(nl)/(cont_ccf(nl)/cont_ccf(cont_reforder))
     end do

     !! this part is from harps_ccfs2area
       area_out = 0._PR
       do i=1,n_cal
          ccf_y = 0._PR
          do ii=1,rv_nbin
             ccf_y(ii) = sum(ccf_sum_weighted(ii,:),dim=1,mask=lmask_sel(:,i))
          end do
          call get_gaussfit()
          if (back_ccf .le. 0._PR) cycle
          sigma_out(i) = sgm_ccf
          center_out(i) = cntr_ccf
          contrast_out(i) = ampl_ccf
          back_out(i) = back_ccf
          area_out(i) =  sqrt(2*pid) * sgm_ccf * ampl_ccf
          write(*,*) sgm_ccf, ampl_ccf, area_out(i)
       end do
       write(*,*) 'AREA_out ', area_out
       call get_pams_all()
       write(*,*) 'Parameters: ',teff_val, gfeh_val, logg_val

       if (teff_val.lt.3000 .or. teff_val.gt.7500) teff_val = 9000
       if (gfeh_val.lt.-3.0 .or. gfeh_val.gt. 1.0) gfeh_val = 3.00
       if (logg_val.lt. 3.0 .or. logg_val.gt. 5.1) logg_val = 8.00

       if ( teff_val.lt.8900 .and. gfeh_val.lt.2.90 .and. logg_val .lt. 6.00 ) then
         call get_errors(snr_pip,teff_err,gfeh_err,logg_err)

         n_count = n_count + 1
         write(output_lun,'(2I6,F18.6,2F9.2,8F12.5,F10.2,F14.7,3x,A)') &
            n_list, n_count, bjd_obs, teff_val, teff_err, gfeh_val, gfeh_err, logg_val,logg_err, &
            area_out, snr_pip, airmass, trim(file_rad)

         write(extended_lun,'(2I6,F18.6,6F9.3,20F15.5,F10.2,F14.7,3x,A)') &
                n_list, n_count, bjd_obs, teff_val, teff_err, gfeh_val, gfeh_err, logg_val,logg_err, &
                center_out(1), sigma_out(1), contrast_out(1), back_out(1), area_out(1), &
                center_out(2), sigma_out(2), contrast_out(2), back_out(2), area_out(2), &
                center_out(3), sigma_out(3), contrast_out(3), back_out(3), area_out(3), &
                center_out(4), sigma_out(4), contrast_out(4), back_out(4), area_out(4), &
                snr_pip, airmass, trim(file_rad)


          snr_tot = sqrt(snr_tot**2 + snr_pip**2)
          ccf_sum_total(:,:) = ccf_sum_total(:,:) + ccf_sum_weighted(:,:)
          bjd_avg = bjd_avg + bjd_obs

       end if




!!$     ! !this part is from harps_ccfs2sumccf
!!$     fits_name = trim(object)//'/'//trim(file_rad)//'_sccf_'//trim(mask_out)//'.fits'
!!$     call fits_write_2d(fits_lun,fits_name,ccf_sum)
!!$     call ccfscale_to_wcs(fits_lun, ccf_pcval,real( 1.00), real(rv_bin))
!!$     call fits_close(fits_lun)
!!$     !! end fo part from harps_ccfs2sumccf


     !! end fo part from harps_ccfs2sumccf


      write(*,*) 'Processing of spectra ', trim(file_rad), '  completed'


  end do
  close(input_lun)



    snr_pip = snr_tot
    bjd_obs = snr_tot / real(n_count)
    !! this part is from harps_ccfs2area
    area_out = 0._PR
    do i=1,n_cal
      ccf_y = 0._PR
      do ii=1,rv_nbin
          ccf_y(ii) = sum(ccf_sum_total(ii,:),dim=1,mask=lmask_sel(:,i))
      end do
      call get_gaussfit()
      if (back_ccf .le. 0._PR) cycle
      sigma_out(i) = sgm_ccf
      center_out(i) = cntr_ccf
      contrast_out(i) = ampl_ccf
      back_out(i) = back_ccf
      area_out(i) =  sqrt(2*pid) * sgm_ccf * ampl_ccf
     !!!! write(*,*) area_out
    end do
    call get_pams_all()


    if (teff_val.lt.3000 .or. teff_val.gt.7500) teff_val = 9000
    if (gfeh_val.lt.-3.0 .or. gfeh_val.gt. 1.0) gfeh_val = 3.00
    if (logg_val.lt. 3.0 .or. logg_val.gt. 5.1) logg_val = 8.00

    call get_errors(snr_pip,teff_err,gfeh_err,logg_err)
    bjd_obs = 0.0
    n_list = 0
    n_count = 0
    airmass = 0.0
    write(output_lun,'(2I6,F18.6,2F9.2,8F12.5,F10.2,F14.7,3x,A)') &
         n_list, n_count, bjd_obs, teff_val, teff_err, gfeh_val, gfeh_err, logg_val,logg_err, &
         area_out, snr_pip, airmass, 'COADDED'

    close(output_lun)

    write(extended_lun,'(2I6,F18.6,6F9.3,20F15.5,F10.2,F14.7,3x,A)') &
          n_list, n_count, bjd_obs, teff_val, teff_err, gfeh_val, gfeh_err, logg_val,logg_err, &
          center_out(1), sigma_out(1), contrast_out(1), back_out(1), area_out(1), &
          center_out(2), sigma_out(2), contrast_out(2), back_out(2), area_out(2), &
          center_out(3), sigma_out(3), contrast_out(3), back_out(3), area_out(3), &
          center_out(4), sigma_out(4), contrast_out(4), back_out(4), area_out(4), &
          snr_pip, airmass, 'COADDED'

    close(extended_lun)

!!$  ! !this part is from harps_ccfs2sumccf
!!$  fits_name = trim(object)//'/'//trim(file_rad)//'_sccf_'//trim(mask_out)//'.fits'
!!$  call fits_write_2d(fits_lun,fits_name,ccf_sum_total)
!!$  call ccfscale_to_wcs(fits_lun, ccf_pcval,real( 1.00), real(rv_bin))
!!$  call fits_close(fits_lun)
!!$  !! end fo part from harps_ccfs2sumccf

  call free_lun(list_lun)
  call free_lun(fits_lun)
  call free_lun(chb_lun)

  call free_lun(out1_lun)
  call free_lun(out2_lun)
  call free_lun(out3_lun)
  call free_lun(input_lun)
  call free_lun(output_lun)

  deallocate(lmask_wl,lmask_ep,lmask_lgf,lmask_sew,lmask_elm)
  deallocate(lmask_wg1,lmask_wd)
  deallocate(ccf_sum,ccf_sum_weighted,ccf_sum_total)
  deallocate(cont_val,cont_fits)

  deallocate(cheb3_coeff_t,cheb3_coeff_m )
  deallocate(cheb2_coeff_g)
  deallocate(poly_coeff)
  deallocate(lmask_sel)
  deallocate(dwpix)
  deallocate(cont_ccf)

  write(*,*)
  write(*,*) 'Program completed'
  write(*,*)
  !! Internal subroutines se  ction
contains

  subroutine get_gaussfit
    implicit none
    !! ****** Second Gaussian fit using just CCF continuum ******
    back_ccf =   maxval(ccf_y,dim=1,mask=(abs(ccf_x).lt.14))
    ccf_z = ccf_y
    ccf_w = 1._PR/sqrt(abs(ccf_x)+1._PR)
    where (abs(ccf_x).gt.14) ccf_w = 0._PR

    ampl_ccf = minval(ccf_z,dim=1,mask=(abs(ccf_x).lt.6)) - back_ccf
    cntr_ccf = 0._PR
    sgm_ccf  = (minloc(abs(ccf_z+ampl_ccf/2),dim=1,mask=(ccf_x.gt.0.and.ccf_x.lt.6))-int(rv_nbin/2))*rv_bin
    !write(*,'(A,5E17.9)') 'Input values for Second Gaussian fit: ' , &
    !     cntr_ccf,sgm_ccf,ampl_ccf,back_ccf,  abs(ampl_ccf)/back_ccf
    call gausfit(ccf_x,ccf_z,ccf_w,cntr_ccf,sgm_ccf,ampl_ccf,back_ccf)

    !write(*,'(A,5E17.9)') 'Output values for Second Gaussian fit:' , &
    !     cntr_ccf,sgm_ccf,ampl_ccf,back_ccf,  abs(ampl_ccf)/back_ccf
    ampl_ccf = abs(ampl_ccf)/back_ccf

    return
  end subroutine get_gaussfit

  subroutine get_pams_all
    implicit none
    teff_val = 0._PR
    gfeh_val = 0._PR
    logg_val = 0._PR

    call  get_pams_cal3(teff_val,gfeh_val)
    call  get_pams_cal2(teff_val,gfeh_val,logg_val)

    return
  end subroutine get_pams_all


  subroutine get_pams_cal3(teff_out,gfeh_out)
    implicit none
    real (PR), intent(out) :: teff_out, gfeh_out
!!! Temperature and Metallicity with no Gravity dependance

    teff_out = -9.9d0
    gfeh_out = -9.9d0

    cheb3_x_in = (area_out(1) - cheb3_xavg(1))/cheb3_xrng(1)
    cheb3_y_in = (area_out(2) - cheb3_xavg(2))/cheb3_xrng(2)
    cheb3_z_in = (area_out(3) - cheb3_xavg(3))/cheb3_xrng(3)

    if (cheb3_x_in.lt.-1._PR .or. cheb3_x_in.gt.1._PR) return
    if (cheb3_y_in.lt.-1._PR .or. cheb3_y_in.gt.1._PR) return
    if (cheb3_z_in.lt.-1._PR .or. cheb3_z_in.gt.1._PR) return

    teff_out = get_cheby_3d(cheb3_x_in,cheb3_y_in,cheb3_z_in,cheb3_coeff_t(:,:,:))
    gfeh_out = get_cheby_3d(cheb3_x_in,cheb3_y_in,cheb3_z_in,cheb3_coeff_m(:,:,:))

    return
  end subroutine get_pams_cal3

  subroutine get_pams_cal2(teff_in,gfeh_in,logg_out)
    implicit none
    real (kind=PR), intent(in)  :: teff_in, gfeh_in
    real (kind=PR), intent(out) :: logg_out
    real (kind=8) :: logg_tmp, area_rmd
    integer :: il, ii, ib, ie

    logg_out = -9.9

    logg_wf = 0.000
    logg_yf = logg_values

    do il=1,logg_iter
       if (.not.(cal2_logg_sel(il))) cycle
       cheb2_x_in = (teff_in - cheb2_xavg(1,il))/cheb2_xrng(1,il)
       cheb2_y_in = (gfeh_in - cheb2_xavg(2,il))/cheb2_xrng(2,il)
       if (cheb2_x_in.lt.-1._PR .or. cheb2_x_in.gt.1._PR) cycle
       if (cheb2_y_in.lt.-1._PR .or. cheb2_y_in.gt.1._PR) cycle
       logg_xf(il) = get_cheby_2d(cheb2_x_in,cheb2_y_in,cheb2_coeff_g(:,:,il))
       if( logg_xf(il).gt.0._PR .and. logg_xf(il).lt.6._PR ) logg_wf(il)=1._PR
    end do
    ii = minloc(abs(logg_xf-area_out(n_cal3+1)),dim=1,mask=(logg_wf.gt.0))
    ib = ii - int(cheb2_polyord/2._PR + 0.5_PR)
    if (ib.lt.1) ib=1
    ie = ib + 2*int(cheb2_polyord/2._PR + 0.5_PR) + 1
    if (ie.gt.logg_iter) ie=logg_iter
    do while (sum(logg_wf(ib:ie)).le.cheb2_polyord+1)
       if (ib.gt.1        ) ib = ib - 1
       if (ie.lt.logg_iter) ie = ie + 1
       if (ib.eq.1 .and. ie.eq.logg_iter) exit
    end do
    if (ib.gt.1) logg_wf(1:ib-1) = 0._PR
    if (ie.lt.logg_iter) logg_wf(ie+1:logg_iter) = 0._PR
    poly_coeff = 0._PR
    call polyfit(logg_xf,logg_yf,logg_wf,poly_coeff)
    logg_tmp = get_poly(area_out(n_cal3+1),poly_coeff)

    area_rmd =  area_out(n_cal3+1) - (gcoeff_y2*(logg_tmp-gcoeff_x0)**2+gcoeff_y1*(logg_tmp-gcoeff_x0)+ gcoeff_y0)
    logg_out = real(get_poly(area_rmd,poly_coeff))
    !!!! write(*,*) ' transformed: ', logg_tmp , ' --> ',logg_out

    return

  end subroutine get_pams_cal2

  subroutine get_errors(snr,teff_err,gfeh_err,logg_err)
    real (PR), intent(in) :: snr
    real (PR), intent(out) :: teff_err,gfeh_err,logg_err

    teff_err = sqrt((50.0d0/snr)*(50.0d0**2) + 30.0d0**2)
    gfeh_err = sqrt((50.0d0/snr)*(0.035d0**2) + 0.03d0**2)
    logg_err = sqrt((50.0d0/snr)*(0.09d0**2) + 0.06d0**2)
    return
  end

end Program harps_input2pams
