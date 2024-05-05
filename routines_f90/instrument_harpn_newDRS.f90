module instrument_harpn_newdrs
  use common
  use fits
  implicit none

  character (len=24) :: head_object
  real :: head_radvel
  real (PR) :: head_mjdobs

  ! fits images informations (HARPS)
  integer, parameter :: n_order = 69, n_pix = 4096

  !parameters for the cross-correlation
  real (PR), parameter :: mask_factor = 0.75_PR
  real (PR), parameter ::  mask_width = 0.01_PR, rv_bin = 0.25_PR
  integer, parameter ::  rv_nbin = 181

!!! ADC variables -> out-dated!
  real (PR) ,parameter :: wl_ref = 5500._PR, wl_start = 3870._PR, wl_range = 3100._PR
  integer ,parameter :: no_ref = 48 !! This is the order where the wl=5500 is located

  real (PR), dimension(n_order), parameter :: sun_fratio_wl = (/ &
       3897.526855,3922.503906,3947.803955,3973.431152,3999.395020,4025.699951,4052.352295,4079.360352,4106.730957,4134.470703,&
       4162.588867,4191.091309,4219.987305,4249.284180,4278.989746,4309.114746,4339.667480,4370.655273,4402.090332,4433.978516,&
       4466.333496,4499.165527,4532.482422,4566.296387,4600.617676,4635.460938,4670.833984,4706.752930,4743.228027,4780.272461,&
       4817.900391,4856.124512,4894.960938,4934.423828,4974.527832,5015.289062,5056.723145,5098.847656,5141.681152,5185.238770,&
       5229.541992,5274.608887,5320.457520,5367.112305,5414.591309,5462.917969,5512.115723,5562.208008,5613.217285,5665.173828,&
       5718.098633,5772.022949,5826.972656,5882.980469,5940.073730,5998.287109,6057.652832,6118.204590,6179.979004,6243.013672,&
       6307.348145,6373.022949,6440.078613,6508.561523,6578.513672,6649.988281,6723.034180,6797.702148,6874.046387 /)
  real (PR), dimension(n_order), parameter :: sun_fratio_fr = (/ &
     0.09123444,0.09520470,0.10491812,0.11562670,0.14275232,0.15516481,0.16324241,0.17960419,0.18659672,0.19208333, &
     0.21250103,0.22123569,0.23343271,0.24274538,0.25218421,0.26695126,0.27169895,0.29289892,0.29748344,0.33527717,&
     0.35688359,0.38282895,0.40076405,0.41762412,0.43370289,0.45078081,0.47687817,0.50921392,0.53813905,0.57327318,&
     0.60110962,0.60215247,0.64524764,0.68501562,0.72149301,0.75691658,0.79910159,0.81737721,0.83772814,0.89564115,&
     0.93638444,0.96261996,0.98816472,0.99493301,1.00203371,1.00938320,1.00684154,1.00083709,0.99643660,1.00469756,&
     1.00000000,0.99421966,0.99982262,0.98707432,0.97677267,0.96507144,0.96048629,0.94806033,0.93232155,0.89727080,&
     0.87964207,0.86183149,0.84010774,0.82472765,0.78342932,0.72417307,0.72469378,0.72518879,0.68109179/)

  logical, dimension(n_order), parameter :: order_mask = (/ &
    .true.,.true.,.true.,.true.,.true.,.true.,.true.,.true.,.true.,.true., &
    .true.,.true.,.true.,.true.,.true.,.false.,.false.,.true.,.true.,.true., &
    .true.,.true.,.true.,.true.,.true.,.true.,.true.,.true.,.true.,.true., &
    .true.,.true.,.true.,.true.,.true.,.true.,.true.,.true.,.true.,.true., &
    .true.,.true.,.true.,.true.,.true.,.true.,.true.,.true.,.true.,.true., &
    .true.,.true.,.false.,.false.,.true.,.true.,.true.,.true.,.true.,.true., &
    .true.,.true.,.true.,.false.,.true.,.true.,.true.,.true.,.true./)

contains
  subroutine singleheader_AIRMASS(unit,airmass,sts)
    integer :: unit, sts
    real (kind=8) ::  airmass
    character (len=nch_com) :: comment
    sts = 0
    airmass = 0._PR
    call FTGKYD(unit,'AIRMASS',airmass,comment,sts);  sts = 0
    return
  end subroutine singleheader_AIRMASS

  subroutine singleheader_BJD(unit,bjd,sts)
    integer :: unit, sts
    real (kind=8) ::  bjd
    character (len=nch_com) :: comment
    sts = 0
    bjd = 0._PR
    call FTGKYD(unit,'HIERARCH TNG QC BJD',bjd,comment,sts);  sts = 0
    return
  end subroutine singleheader_BJD

  subroutine singleheader_RVC(unit,rv_fin,sts)
    integer :: unit, sts
    real (kind=8) ::  rv_fin
    character (len=nch_com) :: comment
    sts = 0
    rv_fin = 0._PR
    call FTGKYD(unit,'HIERARCH TNG DRS CCF RV',rv_fin,comment,sts);  sts = 0
    return
  end subroutine singleheader_RVC

  subroutine singleheader_RVDRIFT(unit,rv_drift,sts)
    integer :: unit, sts
    real (kind=8) ::  rv_drift
    character (len=nch_com) :: comment
    sts = 0
    rv_drift = 0._PR
  !  call FTGKYD(unit,'HIERARCH TNG DRS CAL DRIFT RV USED',rv_drift,comment,sts);  sts = 0
  !  rv_drift = rv_drift/1000._PR
    return
  end subroutine singleheader_RVDRIFT

  subroutine singleheader_BERV(unit,berv,sts)
    integer :: unit, sts
    real (kind=8) ::  berv
    character (len=nch_com) :: comment
    sts = 0
    berv = 0._PR
    call FTGKYD(unit,'HIERARCH TNG QC BERV',berv,comment,sts);  sts = 0
    return
  end subroutine singleheader_BERV

  subroutine singleheader_BLAZE(unit,file_blaze,sts)
    integer :: unit, sts
    character (len=nch_file) :: file_blaze
    character (len=nch_com) :: comment
    sts = 0
    file_blaze = ''
    call FTGKLS(unit,'HIERARCH ESO PRO REC1 CAL27 NAME',file_blaze,comment,sts);  sts = 0
    return
  end subroutine singleheader_BLAZE

  ! Subroutine removed as the corresponding keyword cannot be found
  ! this subroutines are not used by CCFpams anyway

  !subroutine singleheader_WAVE(unit,file_wave,sts)
  !  integer :: unit, sts
  !  character (len=nch_file) :: file_wave
  !  character (len=nch_com) :: comment
  !  sts = 0
  !  file_wave = ''
  !  call FTGKLS(unit,'HIERARCH TNG DRS CAL TH FILE',file_wave,comment,sts);  sts = 0
  !  return
  !end subroutine singleheader_WAVE

  !subroutine singleheader_CCD_SIGDET(unit,ccd_sigdet,sts)
  !  integer :: unit, sts
  !  real  ::  ccd_sigdet
  !  character (len=nch_com) :: comment
  !  sts = 0
  !  ccd_sigdet = 0._4
  !  call FTGKYE(unit,'HIERARCH TNG DRS CCD SIGDET',ccd_sigdet,comment,sts);  sts = 0
  !  return
  !end subroutine singleheader_CCD_SIGDET

  !subroutine singleheader_CCD_GAIN(unit,ccd_gain,sts)
  !  integer :: unit, sts
  !  real (kind=4) :: ccd_gain, ccd_gain_4
  !  real (kind=8) :: ccd_gain_8
  !  character (len=nch_com) :: comment
  !  ccd_gain_4 = 0._4
  !  ccd_gain_8 = 0._8
  !  call FTGKYE(unit,'HIERARCH TNG DRS CCD CONAD',ccd_gain_4,comment,sts);  sts = 0
  !  call FTGKYD(unit,'HIERARCH TNG DRS CCD CONAD',ccd_gain_8,comment,sts);  sts = 0
  !  if (ccd_gain_4 .gt. 0._4) then
  !    ccd_gain = ccd_gain_4
  !  else
  !    ccd_gain = ccd_gain_8
  !  end if
  !  return
  !end subroutine singleheader_CCD_GAIN

  subroutine singleheader_DRS_SNR(unit,order,snr_key,sts)
    integer :: unit, sts
    real  ::  snr_key
    integer :: order
    character (len=nch_com) :: comment
    character (len=2) :: ord_chr
    sts = 0
    snr_key = 0._4
    write(ord_chr,'(I2)') order+1 ! new pipeline starts counting from 1...
    call FTGKYE(unit,'HIERARCH TNG QC ORDER'//trim(ord_chr)//' SNR',snr_key,comment,sts);  sts = 0
    return
  end subroutine singleheader_DRS_SNR

  !subroutine singleheader_pCOEFF(unit,p_coeff,sts)
  !  integer :: unit, sts
  !  real (PR), dimension(:), intent(out) :: p_coeff
  !  character (len=nch_com) :: comment
  !  sts = 0
  !  p_coeff = 0._PR
  !  call FTGKYD(unit,'HIERARCH TNG DRS FLUX CORR COEFF0',p_coeff(1),comment,sts);  sts = 0
  !  call FTGKYD(unit,'HIERARCH TNG DRS FLUX CORR COEFF1',p_coeff(2),comment,sts);  sts = 0
  !  call FTGKYD(unit,'HIERARCH TNG DRS FLUX CORR COEFF2',p_coeff(3),comment,sts);  sts = 0
  !  call FTGKYD(unit,'HIERARCH TNG DRS FLUX CORR COEFF3',p_coeff(4),comment,sts);  sts = 0
  !  call FTGKYD(unit,'HIERARCH TNG DRS FLUX CORR COEFF4',p_coeff(5),comment,sts);  sts = 0
  !  call FTGKYD(unit,'HIERARCH TNG DRS FLUX CORR COEFF5',p_coeff(6),comment,sts);  sts = 0
  !  return
  !end subroutine singleheader_pCOEFF

  subroutine fits_header_getinfo(unit,verbose)
    integer, intent(in) :: unit
    integer, intent(in), optional :: verbose
    integer :: sts, verb
    character (len=64) :: comment

    sts = 0 ; verb = 0
    head_radvel = 0.00
    head_mjdobs = 0.00
    head_object = ''

    if (present(verbose)) verb=verbose
    call FTGKYE(unit,'HIERARCH TNG TEL TARG RADVEL',head_radvel,comment,sts); call reset_sts(sts,0,verb,'TARG RADVEL  ')
    call FTGKYD(unit,'MJD-OBS',head_mjdobs,comment,sts); call reset_sts(sts,0,verb,'MJD-OBS     ')
    call FTGKYS(unit,'OBJECT' ,head_object,comment,sts); call reset_sts(sts,0,verb,'OBJECT      ')
    return
  end subroutine fits_header_getinfo

  subroutine fits_header_putinfo(unit,verbose)
    integer, intent(in) :: unit
    integer, intent(in), optional :: verbose
    integer :: sts, verb
    character (len=64) :: comment

    sts = 0 ; verb = 0
    comment = ''

    if (present(verbose)) verb=verbose
    call FTPKYE(unit,'HIERARCH TNG TEL TARG RADVEL',head_radvel,6,comment,sts); call reset_sts(sts,0,verb,'TARG RADVEL  ')
    call FTPKYD(unit,'MJD-OBS',head_mjdobs,6,comment,sts); call reset_sts(sts,0,verb,'MJD-OBS     ')
    call FTPKYS(unit,'OBJECT' ,head_object,comment,sts);   call reset_sts(sts,0,verb,'OBJECT      ')
    return
  end subroutine fits_header_putinfo



    subroutine key2wave(unit,wl_pix)

    integer, intent(in) :: unit
    real (kind=8), dimension(:,:), intent(out)  :: wl_pix

    character(len=3) :: int2ch
    character (len=64) :: comment

    integer :: sts, verb, d, i, n, a_sel

    real (kind=8) :: a_coeff
    real (kind=8), dimension(size(wl_pix,1)) :: x,y

    sts = 0 ; verb = 6
    call FTGKYK(unit,'HIERARCH TNG DRS CAL TH DEG LL',d,comment,sts); call reset_sts(sts,0,verb,'DEG LL      ')

    do n = 1,size(x)
       x(n) = real(n)-1
    end do

    do n=0,size(wl_pix,dim=2)-1
       do i = d,0,-1
          a_sel = i + n*(1+d)
          if (a_sel.lt.10)                    write(int2ch,'(I1)') a_sel
          if (a_sel.ge.10 .and. a_sel.lt.100) write(int2ch,'(I2)') a_sel
          if (a_sel.ge.100)                   write(int2ch,'(I3)') a_sel
          call FTGKYD(unit,'HIERARCH TNG DRS CAL TH COEFF LL'//trim(int2ch),a_coeff,comment,sts)
          call reset_sts(sts,0,verb,'COEFF LL    ')
          if (i.eq.d)then
             y = a_coeff
          else
             y = y*x + a_coeff
          end if
       end do
       wl_pix(:,n+1) = y
    end do

    return
  end subroutine key2wave

  subroutine key2dwave(unit,dw_pix)

    integer, intent(in) :: unit
    real (kind=8), dimension(:,:), intent(out)  :: dw_pix

    character(len=3) :: int2ch
    character (len=64) :: comment

    integer :: sts, verb, d, i, n, a_sel

    real (kind=8) :: a_coeff
    real (kind=8), dimension(size(dw_pix,1)) :: x,y

    sts = 0 ; verb = 6
    call FTGKYK(unit,'HIERARCH TNG DRS CAL TH DEG LL',d,comment,sts); call reset_sts(sts,0,verb,'DEG LL      ')

    do n = 1,size(x)
       x(n) = real(n)-1
    end do

    do n=0,size(dw_pix,dim=2)-1
       do i = d,1,-1
          a_sel = i + n*(1+d)
          if (a_sel.lt.10)                    write(int2ch,'(I1)') a_sel
          if (a_sel.ge.10 .and. a_sel.lt.100) write(int2ch,'(I2)') a_sel
          if (a_sel.ge.100)                   write(int2ch,'(I3)') a_sel
          call FTGKYD(unit,'HIERARCH TNG DRS CAL TH COEFF LL'//trim(int2ch),a_coeff,comment,sts)
          call reset_sts(sts,0,verb,'COEFF LL    ')
          if (i.eq.d) then
             y = i*a_coeff
          else
             y = y*x + i*a_coeff
          end if
       end do
       dw_pix(:,n+1) = y
    end do

    return
  end subroutine key2dwave

end module instrument_harpn_newdrs
