!
! SPDX-License-Identifier: MIT
!
module mod_rk
  !
  use mpi
  use mod_mom   , only: momad_xyz_tw_cen,momad_xyz_tw_fll,momad_xyz_tw_cen_test
  use mod_mom   , only: momad_xyz_sp
  use mod_sanity, only: flutas_error 
  use mod_types
  !@cuf use cudafor
  !
  implicit none
  !
  real(rp), allocatable, dimension(:,:,:) :: dudtrk, dvdtrk, dwdtrk
  !@cuf attributes(managed) :: dudtrk, dvdtrk, dwdtrk
  !
  private
  public  :: rk,cmpt_time_factors
  !
  contains
  !
  subroutine rk(space_scheme_mom,f_t1,f_t2,nx,ny,nz,dxi,dyi,dzi,nh_d,nh_u,dzci,dzfi,u,v,w,mfx,mfy,mfz,mu,rho, &
     dudtrko,dvdtrko,dwdtrko,nor,psi,eps_int,u_reg,dl,dli,dzc,dzf,halo_u,nh_p,nh_v,f_t12,src_ph,f_t12_o,ls,nornew,psinew)
    !
    !
    !inout velocities are onlu velocities but the outputs are the velocity and density product
    !
    ! subroutine to compute the velocity prediction
    !
    ! Note: --> suitable for both 2nd order Adams-Bashforth and 3rd order low-storage Runge-Kutta;
    !       --> source terms are included in the main program;
    !
    use mod_boiling,     only: correct_mom,use_mod_nse,out_rk_safety,rex,rey,rez
    use mod_param  ,     only: no_outflow,rho1,rho2,is_outflow,IBM_h
    use mod_bound  ,     only: bounduvw
    use mod_common_mpi,  only: MPI_PROC_NULL,left,right,front,back,top,bottom
    !
    implicit none
    !
    character(len=3), intent(in   )                                     :: space_scheme_mom
    real(rp)        , intent(in   )                                     :: f_t1,f_t2,f_t12,f_t12_o
    integer         , intent(in   )                                     :: nx,ny,nz
    real(rp)        , intent(in   )                                     :: dxi,dyi,dzi
    integer         , intent(in   )                                     :: nh_d,nh_u
    real(rp)        , intent(in   ), dimension(1-nh_d:)                 :: dzci,dzfi
    real(rp)        , intent(inout), dimension(1-nh_u:,1-nh_u:,1-nh_u:) :: u,v,w,mfx,mfy,mfz
    real(rp)        , intent(in   ), dimension(     0:,     0:,     0:) :: mu,rho
    real(rp)        , intent(inout), dimension(      :,      :,      :) :: dudtrko,dvdtrko,dwdtrko
    !Test
    real(rp)        , intent(in   ), dimension(0:,0:,0:,1:)             :: nor,nornew
    real(rp)        , intent(in   ), dimension(0:,0:,0:)                :: psi,src_ph,ls,psinew
    real(rp)        , intent(in   ), dimension(3)                       :: dl,dli
    real(rp)        , intent(in   )                                     :: eps_int,u_reg
    real(rp)        , intent(in   ), dimension(1-nh_d:)                 :: dzc,dzf
    integer         , intent(in   ), dimension(3)                       :: halo_u
    integer         , intent(in   )                                     :: nh_p,nh_v
    real(rp), allocatable, dimension(:,:,:)                             :: outflow_ind
    !@cuf attributes(managed) :: dzci, dzfi, u, v, w, dudtrko, dvdtrko, dwdtrko
    !@cuf attributes(managed) :: mu, rho
    !until here nor,psi,eps_int,u_reg
    !
    real(rp) :: rhox,rhoy,rhoz
    !
    integer :: i,j,k
    !
    integer                                                             :: x_start,x_stop
    integer                                                             :: y_start,y_stop
    integer                                                             :: z_start,z_stop
    x_start = 0
    x_stop  = 0
    y_start = 0
    y_stop  = 0
    z_start = 0
    z_stop  = 0   
    if    (is_outflow(0,1).AND.(left .eq.MPI_PROC_NULL)) then
      x_start = out_rk_safety
    elseif(is_outflow(1,1).AND.(right.eq.MPI_PROC_NULL)) then
      x_stop = out_rk_safety
    elseif(is_outflow(0,2).AND.(front.eq.MPI_PROC_NULL)) then
      y_start = out_rk_safety
    elseif(is_outflow(1,2).AND.(back .eq.MPI_PROC_NULL)) then
      y_stop = out_rk_safety
    elseif(is_outflow(0,3)) then
      z_start = out_rk_safety
    elseif(is_outflow(1,3)) then
      z_stop = out_rk_safety
    endif
    !
    allocate(outflow_ind(1:nx,1:ny,1:nz))
    outflow_ind(:,:,:) = 0.0_rp
    do k=IBM_h+1+z_start,nz-z_stop
      do j=1+y_start,ny-y_stop
        do i=x_start+1,nx-x_stop
          outflow_ind(i,j,k) = 1.0_rp
        enddo
      enddo
    enddo
    !
    ! we allocate the arrays dudtrko,dvdtrko and dwdtrko only once
    !
    if(.not.allocated(dudtrk)) allocate(dudtrk(nx,ny,nz))
    if(.not.allocated(dvdtrk)) allocate(dvdtrk(nx,ny,nz))
    if(.not.allocated(dwdtrk)) allocate(dwdtrk(nx,ny,nz))
    !
    ! compute the r.h.s. of the momentum equation 
    !
    !call momad_xyz_tw_cen(nx,ny,nz,dxi,dyi,dzi,nh_d,nh_u,dzci,dzfi,u,v,w,mu,rho,mfx,mfy,mfz,dudtrk,dvdtrk,dwdtrk)  
    call momad_xyz_tw_cen_test(nx,ny,nz,dxi,dyi,dzi,nh_d,nh_u,dzci,dzfi,u,v,w,mu,rho,mfx,mfy,mfz,dudtrk,dvdtrk, &
                                                                     dwdtrk,src_ph,dl,dzc,dzf,halo_u)!,psi)
    !  
    if(use_mod_nse) then
      !
      call correct_mom(nx,ny,nz,dxi,dyi,dzi,nh_d,nh_u,dzci,dzfi,u,v,w,mu,rho &
                        ,dudtrk,dvdtrk,dwdtrk,f_t12,psi,dl,dzc,dzf,halo_u,f_t12_o)
    endif
    !
    ! compute the velocity prediction 
    !
    !
#if defined(_OPENACC)
    !$acc kernels
#else
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(i,j,k) &
    !$OMP SHARED(nx,ny,nz,f_t1,f_t2,u,v,w,dudtrk,dvdtrk,dwdtrk,dudtrko,dvdtrko,dwdtrko)
#endif
    if(use_mod_nse) then
      do k=1,nz
        do j=1,ny
          do i=1,nx
            !
            rhox = 0.5*(rho(i,j,k)+rho(i+1,j,k))
            rhoy = 0.5*(rho(i,j,k)+rho(i,j+1,k))
            rhoz = 0.5*(rho(i,j,k)+rho(i,j,k+1))
            !
            u(i,j,k) = (rhox*u(i,j,k) - outflow_ind(i,j,k)*(rho1-rho2)*rex(i,j,k)) &
                                      + f_t1*dudtrk(i,j,k) + f_t2*dudtrko(i,j,k)
            v(i,j,k) = (rhoy*v(i,j,k) - outflow_ind(i,j,k)*(rho1-rho2)*rey(i,j,k)) &
                                      + f_t1*dvdtrk(i,j,k) + f_t2*dvdtrko(i,j,k)
            w(i,j,k) = (rhoz*w(i,j,k) - outflow_ind(i,j,k)*(rho1-rho2)*rez(i,j,k)) &
                                      + f_t1*dwdtrk(i,j,k) + f_t2*dwdtrko(i,j,k)
            dudtrko(i,j,k) = dudtrk(i,j,k)
            dvdtrko(i,j,k) = dvdtrk(i,j,k)
            dwdtrko(i,j,k) = dwdtrk(i,j,k)
            !
          enddo
        enddo
      enddo
    else
      do k=1,nz
        do j=1,ny
          do i=1,nx
            !
            rhox = 0.5*(rho(i,j,k)+rho(i+1,j,k))
            rhoy = 0.5*(rho(i,j,k)+rho(i,j+1,k))
            rhoz = 0.5*(rho(i,j,k)+rho(i,j,k+1))
            !
            u(i,j,k) = rhox*u(i,j,k) + f_t1*dudtrk(i,j,k) + f_t2*dudtrko(i,j,k)! &
            v(i,j,k) = rhoy*v(i,j,k) + f_t1*dvdtrk(i,j,k) + f_t2*dvdtrko(i,j,k)! &
            w(i,j,k) = rhoz*w(i,j,k) + f_t1*dwdtrk(i,j,k) + f_t2*dwdtrko(i,j,k)! &
            !
            dudtrko(i,j,k) = dudtrk(i,j,k)
            dvdtrko(i,j,k) = dvdtrk(i,j,k)
            dwdtrko(i,j,k) = dwdtrk(i,j,k)
            !
          enddo
        enddo
      enddo
    endif
#if defined(_OPENACC)
    !$acc end kernels 
#else
    !$OMP END PARALLEL DO
#endif
    !
    deallocate(outflow_ind)
    return
  end subroutine rk
  !
  subroutine cmpt_time_factors(time_scheme,restart,istep,irk_ss,rkcoeff,dt,dto,f_t1,f_t2,f_t12)
    !
    ! compute f_t1 and f_t2 according to the input time_scheme 
    !
    implicit none
    !
    character(len=3), intent(in )                 :: time_scheme
    logical         , intent(in )                 :: restart
    integer         , intent(in )                 :: istep,irk_ss
    real(rp)        , intent(in ), dimension(2,3) :: rkcoeff
    real(rp)        , intent(in )                 :: dt,dto
    real(rp)        , intent(out)                 :: f_t1,f_t2,f_t12
    !
    select case(time_scheme)
    case('ab2')
      if(istep.eq.1) then ! we use 1st order Euler in the first time-step (istep.eq.1)
        f_t1 = 1._rp*dt
        f_t2 = 0._rp*dt
      else
        f_t1 =  (1._rp+0.5_rp*(dt/dto))*dt   !1._rp*dt!
        f_t2 =  (     -0.5_rp*(dt/dto))*dt   !0._rp*dt!
      endif
    case('rk3')
      f_t1 = rkcoeff(1,irk_ss)*dt
      f_t2 = rkcoeff(2,irk_ss)*dt
    case default
      call flutas_error('Invalid time discretization: only AB2 and RK3 available for now - check dns.in')
    end select
    f_t12 = f_t1+f_t2
    !
    return
  end subroutine cmpt_time_factors
  !
end module mod_rk
