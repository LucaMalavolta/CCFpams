module ccf
  use common
  implicit none
  
contains
  
  subroutine shift_mask(mask_wl,mask_wg,mask_wd,radvel,flux,wl_scale,dw_scale,ccf_val)
    real (kind=8), dimension(:), intent(in) :: mask_wl, mask_wg, mask_wd, wl_scale, flux,dw_scale
    real (kind=8), intent(in) ::  radvel
    
    integer ::  ns, nl, i, ii, mask_ie, mask_ib, ib, ir
    integer :: mask_size, sci_size, ns_prv
    real (kind=8) :: mask_wlb,mask_wle, pix_e, pix_s
    
    real (kind=8) ::    ccf_val

    !mask_wl_shift = (1.000d0+1.55d-8)*(radvel / cc + 1.000000D0) * mask_wl
    
    if (size(mask_wl,1 ).ne.size(mask_wg)) stop
    
    sci_size  = size(wl_scale,1)
    mask_size = size(mask_wl,1)

    mask_ib = 0
    mask_ie = 0
    do nl = 1, mask_size
       if (mask_wg(nl).gt.0.00d0) then
          mask_ib = nl
          exit
       end if
    end do
    do nl = mask_size, 1,-1
       if (mask_wg(nl).gt.0.00d0) then
          mask_ie = nl
          exit
       end if
    end do
    
    if (mask_ie.lt.mask_ib) return
        
    ccf_val  = 0.00d0
    
    ns_prv = 1
    do nl = mask_ib,mask_ie
       
       mask_wlb =  (radvel  / cc + 1.000000D0) * (mask_wl(nl)-0.5000d0*mask_wd(nl))
       mask_wle  = (radvel  / cc + 1.000000D0) * (mask_wl(nl)+0.5000d0*mask_wd(nl))
       !if (mask_wg(nl).eq.0.00) cycle
       
       ib = ns_prv ;  ir = ns_prv
       do ns=ns_prv,sci_size-1
          if (wl_scale(ib)+dw_scale(ib)*0.500d0.lt.mask_wlb) ib = ib+1
          if (wl_scale(ir)+dw_scale(ir)*0.500d0.lt.mask_wle) ir = ir+1
          if (ir.lt.ns-3) exit
       end do
       ns_prv = ib - 2
       
       if (ib.le.0 .or. ir.gt.sci_size) cycle
       if (ib.gt.ir) cycle
       
       if (ib.eq.ir)  then
          pix_s    = (mask_wle-mask_wlb ) / dw_scale(ib)
          ccf_val   = ccf_val   + pix_s*flux(ib)*mask_wg(nl)
       else if (ib+1.eq.ir) then
          pix_s     = (wl_scale(ib)+dw_scale(ib)*0.5-mask_wlb) / dw_scale(ib)
          !pix_e     = (mask_wle-wl_scale(ib)-dw_scale(ib)*0.5) / dw_scale(ir)
          pix_e     = (mask_wle-(wl_scale(ir)-dw_scale(ir)*0.5))/ dw_scale(ir)
          ccf_val   = ccf_val   + (pix_s*flux(ib) + pix_e*flux(ir))*mask_wg(nl)
       else 
          pix_s     = (wl_scale(ib)+dw_scale(ib)*0.5-mask_wlb)/dw_scale(ib)
          !pix_e     = (mask_wle-wl_scale(ir-1)-dw_scale(ir-1)*0.5)/ dw_scale(ir)
          pix_e     = (mask_wle-(wl_scale(ir)-dw_scale(ir)*0.5))/ dw_scale(ir)
          ccf_val   = ccf_val   + (pix_s*flux(ib)+pix_e*flux(ir))*mask_wg(nl)
          do i=ib+1,ir-1
             ccf_val   = ccf_val   + flux(i)*mask_wg(nl)
          end do
       end if
       !pay attention to the ir-i index, 
       !I made a terrible mess before realizing the correct indexing
       
    enddo
    return
    
  end subroutine shift_mask
  
  subroutine flux2ccf(mask_wl,mask_wg,mask_wd,flux,wl_scale,dw_scale,ccf_x,ccf_y)
    
    real (kind=8), dimension(:), intent(in) :: mask_wl,mask_wg,mask_wd,flux,wl_scale,dw_scale,ccf_x
    real (kind=8),dimension(:),intent(out) :: ccf_y
    
    integer :: ii
    real (kind=8) :: rv, ccf_val
    
    if (size(ccf_y).ne.size(ccf_x)) stop
    
    ccf_y = 0.00D0
    do ii = 1, size(ccf_x,1)
       rv = ccf_x(ii)
       call shift_mask(mask_wl,mask_wg,mask_wd,rv,flux,wl_scale, dw_scale,ccf_val )
       ccf_y(ii)=ccf_val
    end do
    return
  end subroutine flux2ccf
  
  subroutine ccf_singleline(line_wl,line_wd,flux,wl_scale,dw_scale,ccf_x,ccf_y)
  !!!--- CCF_SINGLELINE ---!!!
  ! This subroutine do the same of flux2ccf and shift_mask combined, with the difference that the input line list
  ! consists of a single line; the subroutine has been optimized to speed up processing time
  !!! LINE is NOT weighted: since we are processing a single line, the result can be weighted later
  !!! (if necessary)

    real (kind=8), dimension(:), intent(in) :: flux, wl_scale,dw_scale,ccf_x
    real (kind=8), intent(in) :: line_wl, line_wd
    real (kind=8),dimension(:),intent(out) :: ccf_y
    
    integer ::  i,ii,ns, wl_size, ib,ir
    real (kind=8) :: rv, line_wlb,line_wle, pix_e, pix_s
    
    wl_size = size(wl_scale,1)
    if (size(ccf_y).ne.size(ccf_x)) stop
    
    ccf_y = 0.000D0
    do ii = 1, size(ccf_x,1)
       rv = ccf_x(ii)
       !! Here starts the part taken from shift_mask
       line_wlb =  (rv  / cc + 1.000000D0) * (line_wl-0.5000d0*line_wd)
       line_wle  = (rv  / cc + 1.000000D0) * (line_wl+0.5000d0*line_wd)
       
       ib = 2 ;  ir = 2
       do ns=2,wl_size-1
          if (wl_scale(ib)+dw_scale(ib)*0.500d0.lt.line_wlb) ib = ib+1
          if (wl_scale(ir)+dw_scale(ir)*0.500d0.lt.line_wle) ir = ir+1
          if (ir.lt.ns-3) exit
       end do
       
       if (ib.le.0 .or. ir.gt.wl_size) cycle
       if (ib.gt.ir) cycle
       
       if (ib.eq.ir)  then
          pix_s       = (line_wle-line_wlb ) / dw_scale(ib)
          ccf_y(ii)   =  pix_s*flux(ib)
       else if (ib+1.eq.ir) then
          pix_s     = (wl_scale(ib)+dw_scale(ib)*0.5-line_wlb) / dw_scale(ib)
          pix_e     = (line_wle-(wl_scale(ir)-dw_scale(ir)*0.5))/ dw_scale(ir)
          !pix_e     = (line_wle-wl_scale(ib)-dw_scale(ib)*0.5) / dw_scale(ir)
          ccf_y(ii) = pix_s*flux(ib) + pix_e*flux(ir)
       else 
          pix_s     = (wl_scale(ib)+dw_scale(ib)*0.5-line_wlb)/dw_scale(ib)
          pix_e     = (line_wle-(wl_scale(ir)-dw_scale(ir)*0.5))/ dw_scale(ir)
          ccf_y(ii)   = pix_s*flux(ib)+pix_e*flux(ir)
          do i=ib+1,ir-1
             ccf_y(ii)   = ccf_y(ii)   + flux(i)
          end do

       end if
       !write(*,*) ii, line_wl,line_wd,(line_wlb+line_wle)/2,'->', ccf_x(ii),ccf_y(ii)
    end do
    !stop
    return
  end subroutine ccf_singleline

  real (kind=8) function delta_raw(sigma,dx)
    real (kind=8) :: sigma, dx
    delta_raw = exp(-dx**2/(2*PId*sigma**2))
    return
  end function delta_raw
  
end module ccf
