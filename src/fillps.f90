!
! SPDX-License-Identifier: MIT
!
module mod_fillps
  !
  use mod_types, only: rp
  !@cuf use cudafor
  !
  implicit none
  !
  private
  public  :: fillps
  !
  contains
  !
  subroutine fillps(nx,ny,nz,nh_d,nh_u,dxi,dyi,dzi,dzfi,dti,rho0,rho1,rho2,rho,u_grad_rho, &
         src,u,v,w,p,mom_based,mfx,mfy,mfz,nse_switch,out_rk_safety)
    !
    !u and v and w contains density as well
    use mod_param,      only: IBM_h,is_outflow
    use mod_common_mpi, only: MPI_PROC_NULL,left,right,front,back,top,bottom
    implicit none
    !
    integer , intent(in )                                     :: nx,ny,nz,out_rk_safety
    integer , intent(in )                                     :: nh_d,nh_u
    real(rp), intent(in )                                     :: dxi,dyi,dzi,rho1,rho2,nse_switch
    real(rp), intent(in ), dimension(1-nh_d:)                 :: dzfi
    real(rp), intent(in )                                     :: dti,rho0
    real(rp), intent(in ), dimension(1-nh_u:,1-nh_u:,1-nh_u:) :: u,v,w,mfx,mfy,mfz
    real(rp), intent(in ), dimension(1-nh_u:,1-nh_u:,1-nh_u:) :: src,rho
    real(rp), intent(in ), dimension(1:,1:,1:)                :: u_grad_rho
    real(rp), intent(out), dimension(0:,0:,0:)                :: p
    logical , intent(in )                                     :: mom_based
    !
    real(rp) :: dtidxi,dtidyi,rxm,rxp,rym,ryp,rzm,rzp
    integer  :: i,j,k,im,jm,km
    real(rp), allocatable, dimension(:,:,:)                   :: outflow_ind
    integer                                                             :: x_start,x_stop
    integer                                                             :: y_start,y_stop
    integer                                                             :: z_start,z_stop
    !
    !@cuf attributes(managed) :: p, u, v, w,mfx,mfy,mfz,src,rho,outflow_ind
    !
    x_start = 0
    x_stop  = 0
    y_start = 0
    y_stop  = 0
    z_start = 0
    z_stop  = 0   
    if    (is_outflow(0,1).AND.(left .eq.MPI_PROC_NULL)) then
      x_start = out_rk_safety+1
    elseif(is_outflow(1,1).AND.(right.eq.MPI_PROC_NULL)) then
      x_stop = out_rk_safety+1
    elseif(is_outflow(0,2).AND.(front.eq.MPI_PROC_NULL)) then
      y_start = out_rk_safety+1
    elseif(is_outflow(1,2).AND.(back .eq.MPI_PROC_NULL)) then
      y_stop = out_rk_safety+1
    elseif(is_outflow(0,3)) then
      z_start = out_rk_safety+1
    elseif(is_outflow(1,3)) then
      z_stop = out_rk_safety+1
    endif
    !
    allocate(outflow_ind(0:nx+1,0:ny+1,0:nz+1))
    outflow_ind(:,:,:) = 0.0_rp
    !$acc kernels
    do k=IBM_h+z_start,nz-z_stop+1
      do j=y_start,ny-y_stop+1
        do i=x_start,nx-x_stop+1
          outflow_ind(i,j,k) = 1.0_rp
        enddo
      enddo
    enddo
    !$acc end kernels 
    !
    dtidxi = dti*dxi
    dtidyi = dti*dyi
    p = 0.0_rp
    !
#if defined(_OPENACC)
    !$acc kernels
#else
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP SHARED(nx,ny,nz,rho0,p,u,v,w,dtidxi,dtidyi,dti,dzci,u_grad_rho,rho,src,rho1,rho2) &
    !$OMP PRIVATE(i,j,k,im,jm,km)
#endif
    if(.NOT.mom_based) then
      do k=1,nz
        do j=1,ny
          do i=1,nx
            !
            im = i-1
            jm = j-1
            km = k-1
            !    
            rxm = 0.5*(rho(i,j,k)+rho(i-1,j  ,k  ))
            rxp = 0.5*(rho(i,j,k)+rho(i+1,j  ,k  ))
            rym = 0.5*(rho(i,j,k)+rho(i  ,j-1,k  ))
            ryp = 0.5*(rho(i,j,k)+rho(i  ,j+1,k  ))
            rzm = 0.5*(rho(i,j,k)+rho(i  ,j  ,k-1))
            rzp = 0.5*(rho(i,j,k)+rho(i  ,j  ,k+1))
            !
            p(i,j,k) =  &
                       ((w(i,j,k)-w(i,j,km))*dti*dzfi(k) + &
                        (v(i,j,k)-v(i,jm,k))*dtidyi        &
                       +(u(i,j,k)-u(im,j,k))*dtidxi)*rho0 - &
                        src(i,j,k)*(1.d0/rho1 - 1.d0/rho2)*dti*rho0 &
     +(rho1-rho2)*((outflow_ind(i,j,k)*mfz(i,j,k)/rzp-outflow_ind(i,j,km)*mfz(i,j,km)/rzm)*dti*dzfi(k) + &
                   (outflow_ind(i,j,k)*mfy(i,j,k)/ryp-outflow_ind(i,jm,k)*mfy(i,jm,k)/rym)*dtidyi      + &
                   (outflow_ind(i,j,k)*mfx(i,j,k)/rxp-outflow_ind(im,j,k)*mfx(im,j,k)/rxm)*dtidxi)*rho0*nse_switch                  
          enddo
        enddo
      enddo
    else
      do k=1,nz
        do j=1,ny
          do i=1,nx
            !
            im = i-1
            jm = j-1
            km = k-1
            !
            p(i,j,k) =  &
                      ((w(i,j,k)-w(i,j,km))*dti*dzfi(k) + &
                       (v(i,j,k)-v(i,jm,k))*dtidyi     &
                      +(u(i,j,k)-u(im,j,k))*dtidxi) - &
                      (u_grad_rho(i,j,k))*dti           &
                      -rho(i,j,k)*src(i,j,k)*(1.d0/rho1 - 1.d0/rho2)*dti  &
          +(rho1-rho2)*((outflow_ind(i,j,k)*mfz(i,j,k)-outflow_ind(i,j,km)*mfz(i,j,km))*dti*dzfi(k) + &
                        (outflow_ind(i,j,k)*mfy(i,j,k)-outflow_ind(i,jm,k)*mfy(i,jm,k))*dtidyi      + &
                        (outflow_ind(i,j,k)*mfx(i,j,k)-outflow_ind(im,j,k)*mfx(im,j,k))*dtidxi      )*nse_switch

          enddo
        enddo
      enddo
    endif    
#if defined(_OPENACC)
    !$acc end kernels 
#else
    !$OMP END PARALLEL DO
#endif
    deallocate(outflow_ind)
    !
    return
  end subroutine fillps
  !
end module mod_fillps
