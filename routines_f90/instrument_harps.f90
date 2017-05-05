module instrument_harps
  use common
  use fits
  implicit none

  character (len=24) :: head_object
  real :: head_radvel
  real (PR) :: head_mjdobs

  ! fits images informations (HARPS)
  integer, parameter :: n_order = 72, n_pix = 4096

  !parameters for the cross-correlation
  real (PR), parameter :: mask_factor = 0.75_PR
  real (PR), parameter ::  mask_width = 0.01_PR, rv_bin = 0.25_PR
  integer, parameter ::  rv_nbin = 181

!!! ADC variables:
  real (PR) ,parameter :: wl_ref = 5500._PR, wl_start = 3880._PR, wl_range = 3100._PR
  integer ,parameter :: no_ref = 50 !! This is the order where the wl=5500 is located

  real (PR), dimension(n_order), parameter :: sun_fratio_wl = (/ &
       3802.510498,3826.271973,3850.332275,3874.697754,3899.372070,3924.363281,3949.677490,3975.319336,4001.297119,4027.616699,&
       4054.283936,4081.307373,4108.693359,4136.448242,4164.583008,4193.102051,4222.012695,4251.326660,4281.049805,4311.193359,&
       4341.762695,4372.769043,4404.221191,4436.128906,4468.502441,4501.351074,4534.688477,4568.522461,4602.863281,4637.727051,&
       4673.121582,4709.060547,4745.556152,4782.623047,4820.273438,4858.520020,4897.378418,4936.864258,4976.992188,5017.778809,&
       5059.237305,5101.386719,5144.246094,5187.830078,5232.160156,5277.253906,5369.783203,5417.288086,5465.641113,5514.863281,&
       5564.983398,5616.020508,5668.003906,5720.956543,5774.909180,5829.889160,5885.926270,5943.050781,6001.296387,6060.693848,&
       6121.278320,6183.086426,6246.154785,6310.525391,6376.234863,6443.326660,6511.845703,6581.838379,6653.352051,6726.436523,&
       6801.145508,6877.530762 /)
  real (PR), dimension(n_order), parameter :: sun_fratio_fr = (/ &
       0.18457793,0.16224839,0.22078945,0.24485369,0.27570051,0.27919459,0.29702580,0.35933146,0.42998219,0.45696911,&
       0.48872158,0.53037584,0.58232337,0.60384458,0.63682848,0.68755084,0.72581655,0.77243394,0.80721265,0.79464334,&
       0.84515959,0.88998204,0.91599000,0.96766847,0.99118561,1.02519357,1.04644203,1.05707550,1.07670105,1.07888186,&
       1.08538854,1.09213519,1.09652936,1.09957159,1.10992992,1.08735752,1.10215831,1.10551476,1.10882211,1.08303225,&
       1.08185148,1.08086920,1.09056175,1.08499134,1.08898497,1.08152997,1.02332914,1.01081014,1.01026058,1.01168680,&
       1.00000000,0.99162650,0.99851990,0.98995537,0.98053437,1.00373387,0.98397088,0.96174371,0.97236907,0.98990828,&
       1.00325775,1.01358664,1.00989652,0.98933047,0.99231029,0.98956072,0.94818705,0.94377893,0.93649596,0.89615142,&
       0.86011368,0.80748463 /)

  logical, dimension(n_order), parameter :: order_mask = (/ &
       .false., .false., .true., .true., .true., .true., .true., .true., .true., .true.,&
       .true., .true., .true., .true., .true., .true., .true., .true., .true., .false.,&
       .false., .true., .true., .true., .true., .true., .true., .true., .true., .true.,&
       .true., .true., .true., .true., .true., .true., .true., .true., .true., .true.,&
       .true., .true., .true., .true., .true., .true., .true., .true., .true., .true.,&
       .true., .true., .true., .true., .true., .false., .false., .true., .true., .true.,&
       .true., .true., .true., .true., .true., .true., .false., .true., .true., .true.,&
       .true., .true. /)

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
    call FTGKYD(unit,'HIERARCH ESO DRS BJD',bjd,comment,sts);  sts = 0
    return
  end subroutine singleheader_BJD

  subroutine singleheader_RVC(unit,rv_fin,sts)
    integer :: unit, sts
    real (kind=8) ::  rv_fin
    character (len=nch_com) :: comment
    sts = 0
    rv_fin = 0._PR
    call FTGKYD(unit,'HIERARCH ESO DRS CCF RVC',rv_fin,comment,sts);  sts = 0
    return
  end subroutine singleheader_RVC

  subroutine singleheader_RVDRIFT(unit,rv_drift,sts)
    integer :: unit, sts
    real (kind=8) ::  rv_drift
    character (len=nch_com) :: comment
    sts = 0
    rv_drift = 0._PR
    call FTGKYD(unit,'HIERARCH ESO DRS CAL DRIFT RV USED',rv_drift,comment,sts);  sts = 0
    rv_drift = rv_drift/1000._PR
    return
  end subroutine singleheader_RVDRIFT

  subroutine singleheader_BERV(unit,berv,sts)
    integer :: unit, sts
    real (kind=8) ::  berv
    character (len=nch_com) :: comment
    sts = 0
    berv = 0._PR
    call FTGKYD(unit,'HIERARCH ESO DRS BERV',berv,comment,sts);  sts = 0
    return
  end subroutine singleheader_BERV

  subroutine singleheader_BLAZE(unit,file_blaze,sts)
    integer :: unit, sts
    character (len=nch_file) :: file_blaze
    character (len=nch_com) :: comment
    sts = 0
    file_blaze = ''
    call FTGKLS(unit,'HIERARCH ESO DRS BLAZE FILE',file_blaze,comment,sts);  sts = 0
    return
  end subroutine singleheader_BLAZE

  subroutine singleheader_WAVE(unit,file_wave,sts)
    integer :: unit, sts
    character (len=nch_file) :: file_wave
    character (len=nch_com) :: comment
    sts = 0
    file_wave = ''
    call FTGKLS(unit,'HIERARCH ESO DRS CAL TH FILE',file_wave,comment,sts);  sts = 0
    return
  end subroutine singleheader_WAVE

  subroutine singleheader_CCD_SIGDET(unit,ccd_sigdet,sts)
    integer :: unit, sts
    real  ::  ccd_sigdet
    character (len=nch_com) :: comment
    sts = 0
    ccd_sigdet = 0._4
    call FTGKYE(unit,'HIERARCH ESO DRS CCD SIGDET',ccd_sigdet,comment,sts);  sts = 0
    return
  end subroutine singleheader_CCD_SIGDET

  subroutine singleheader_CCD_GAIN(unit,ccd_gain,sts)
    integer :: unit, sts
    real (kind=4) :: ccd_gain, ccd_gain_4
    real (kind=8) :: ccd_gain_8
    character (len=nch_com) :: comment
    ccd_gain_4 = 0._4
    ccd_gain_8 = 0._8
    call FTGKYE(unit,'HIERARCH ESO DRS CCD CONAD',ccd_gain_4,comment,sts);  sts = 0
    call FTGKYD(unit,'HIERARCH ESO DRS CCD CONAD',ccd_gain_8,comment,sts);  sts = 0
    if (ccd_gain_4 .gt. 0._4) then
      ccd_gain = ccd_gain_4
    else
      ccd_gain = ccd_gain_8
    end if
    return
  end subroutine singleheader_CCD_GAIN

  subroutine singleheader_DRS_SNR(unit,order,snr_key,sts)
    integer :: unit, sts
    real  ::  snr_key
    integer :: order
    character (len=nch_com) :: comment
    character (len=2) :: ord_chr
    sts = 0
    snr_key = 0._4
    write(ord_chr,'(I2)') order
    call FTGKYE(unit,'HIERARCH ESO DRS SPE EXT SN'//trim(ord_chr),snr_key,comment,sts);  sts = 0
    return
  end subroutine singleheader_DRS_SNR

  subroutine singleheader_pCOEFF(unit,p_coeff,sts)
    integer :: unit, sts
    real (PR), dimension(:), intent(out) :: p_coeff
    character (len=nch_com) :: comment
    sts = 0
    p_coeff = 0._PR
    call FTGKYD(unit,'HIERARCH ESO DRS FLUX CORR COEFF0',p_coeff(1),comment,sts);  sts = 0
    call FTGKYD(unit,'HIERARCH ESO DRS FLUX CORR COEFF1',p_coeff(2),comment,sts);  sts = 0
    call FTGKYD(unit,'HIERARCH ESO DRS FLUX CORR COEFF2',p_coeff(3),comment,sts);  sts = 0
    call FTGKYD(unit,'HIERARCH ESO DRS FLUX CORR COEFF3',p_coeff(4),comment,sts);  sts = 0
    call FTGKYD(unit,'HIERARCH ESO DRS FLUX CORR COEFF4',p_coeff(5),comment,sts);  sts = 0
    call FTGKYD(unit,'HIERARCH ESO DRS FLUX CORR COEFF5',p_coeff(6),comment,sts);  sts = 0
    return
  end subroutine singleheader_pCOEFF

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
    call FTGKYE(unit,'HIERARCH ESO TEL TARG RADVEL',head_radvel,comment,sts); call reset_sts(sts,0,verb,'TARG RADVEL  ')
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
    call FTPKYE(unit,'HIERARCH ESO TEL TARG RADVEL',head_radvel,6,comment,sts); call reset_sts(sts,0,verb,'TARG RADVEL  ')
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
    call FTGKYK(unit,'HIERARCH ESO DRS CAL TH DEG LL',d,comment,sts); call reset_sts(sts,0,verb,'DEG LL      ')

    do n = 1,size(x)
       x(n) = real(n)-1
    end do

    do n=0,size(wl_pix,dim=2)-1
       do i = d,0,-1
          a_sel = i + n*(1+d)
          if (a_sel.lt.10)                    write(int2ch,'(I1)') a_sel
          if (a_sel.ge.10 .and. a_sel.lt.100) write(int2ch,'(I2)') a_sel
          if (a_sel.ge.100)                   write(int2ch,'(I3)') a_sel
          call FTGKYD(unit,'HIERARCH ESO DRS CAL TH COEFF LL'//trim(int2ch),a_coeff,comment,sts)
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
    call FTGKYK(unit,'HIERARCH ESO DRS CAL TH DEG LL',d,comment,sts); call reset_sts(sts,0,verb,'DEG LL      ')

    do n = 1,size(x)
       x(n) = real(n)-1
    end do

    do n=0,size(dw_pix,dim=2)-1
       do i = d,1,-1
          a_sel = i + n*(1+d)
          if (a_sel.lt.10)                    write(int2ch,'(I1)') a_sel
          if (a_sel.ge.10 .and. a_sel.lt.100) write(int2ch,'(I2)') a_sel
          if (a_sel.ge.100)                   write(int2ch,'(I3)') a_sel
          call FTGKYD(unit,'HIERARCH ESO DRS CAL TH COEFF LL'//trim(int2ch),a_coeff,comment,sts)
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

end module instrument_harps
