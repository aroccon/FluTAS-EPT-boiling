! SPDX-License-Identifier: MIT
!
module mod_di_acdi
  !
  ! based on:
  ! --> Mirjalili, S. et. al. "A conservative diffuse interface method for two-phase flows
  !     with provable boundedness properties." 
  !     Journal of Computational Physics 401 (2020): 109006.
  ! --> Jain, S. S. "Accurate conservative phase-field method 
  !     for simulation of two-phase flows." 
  !     Journal of Computational Physics 469 (2022): 111529.
  !
  ! --> Attraction model developed by Salar Zamani & Armin Shahmardi  
  !
  use mod_types
  use mod_common_mpi, only: myid,ierr,ijk_start,comm_cart
  use mod_sanity    , only: flutas_error
  !
  implicit none
  !
  private
  public :: update_property,rk_psi,rk_ent,psi_to_ls,ls_to_psi,cmpt_norm,initls,initpsi,cmpt_umax,inittmp, &
            rk_psi_ibm,rk_ent_ibm
  public rk_c1,initc1,cmpt_rel_vel, interfacial_cflux

  public static_contact_angle_m


  public cmpt_delta,marangoni_he
  !
  contains
    !
    !
  subroutine rk_ent(f_t1,f_t2,n,dli,nh_d,nh_u,nh_v,dzc,dzf,eps_int,u_reg, &
                      u,v,w,s,nor,tmp,cond,ent,dentdtrko)
    !
    ! low-storage 3rd-order Runge-Kutta scheme 
    ! for time integration of the Enthalpy field.
    !
    implicit none
    !
    real(rp), intent(in   )                                        :: f_t1,f_t2
    real(rp), intent(in   )                                        :: eps_int,u_reg
    integer , intent(in   ), dimension(3)                          :: n
    real(rp), intent(in   ), dimension(3)                          :: dli
    integer , intent(in   )                                        :: nh_d,nh_u,nh_v
    real(rp), intent(in   ), dimension(1-nh_d:)                    :: dzc,dzf
    real(rp), intent(in   ), dimension(1-nh_u:,1-nh_u:,1-nh_u:)    :: u,v,w
    real(rp), intent(in   ), dimension(1-nh_v:,1-nh_v:,1-nh_v:)    :: cond,s 
    real(rp), intent(inout), dimension(1-nh_v:,1-nh_v:,1-nh_v:)    :: tmp
    real(rp), intent(in   ), dimension(0:,0:,0:,1:)                :: nor
    real(rp), intent(inout), dimension(1-nh_v:,1-nh_v:,1-nh_v:)    :: ent  
    real(rp), intent(inout), dimension(     1:,     1:,     1:)    :: dentdtrko
    !
    real(rp) :: gradpxp,gradpxm,gradpyp,gradpym,gradpzp,gradpzm
    real(rp) :: psifxp,psifxm,psifyp,psifym,psifzp,psifzm
    real(rp) :: regxp,regxm,regyp,regym,regzp,regzm
    real(rp) :: tmpxp,tmpxm,tmpyp,tmpym,tmpzp,tmpzm
    !
    real(rp), dimension(n(1),n(2),n(3)) :: dentdtrk
    real(rp)  ::mxp,mxm,myp,mym,mzp,mzm
    integer  :: i,j,k,ip,im,jp,jm,kp,km,zp,zm
    !
    ! 1. compute diffusion
    !
    call diff_ent(n,dli,nh_d,nh_v,dzc,dzf,cond,tmp,dentdtrk)
    !
    !2 . Advection 
    call moms_reg(n,dli,nh_d,nh_u,nh_v,u,v,w,ent,dentdtrk)
    !
    ! 3. Advance the color function
    
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
        !
        ent(i,j,k)       = ent(i,j,k) + f_t1*dentdtrk(i,j,k) + f_t2*dentdtrko(i,j,k)
        dentdtrko(i,j,k) = dentdtrk(i,j,k)
        !
        enddo
      enddo
    enddo
    !
    return
  end subroutine rk_ent
  !
subroutine rk_ent_ibm(f_t1,f_t2,n,dli,nh_d,nh_u,nh_v,dzc,dzf,eps_int,u_reg, &
    u,v,w,s,nor,tmp,cond,ent,dentdtrko,rhocp,halo_u,dl,nh_p,psi,ls,f_t12,nh_van, &
     halo_van,ent_bound_o)
  !
  !
  use mod_boiling,     only: rex,rey,rez,boundsb_arrayo
  use mod_param  ,     only: no_outflow,rho1,rho2,cp1,cp2,cbctmp,bctmp
  use mod_bound  ,     only: bounduvw,boundp
  use mod_gradls ,     only: weno5
  ! low-storage 3rd-order Runge-Kutta scheme 
  ! for time integration of the Enthalpy field.
  !
  implicit none
  !
  real(rp), intent(in   )                                        :: f_t1,f_t2,f_t12
  real(rp), intent(in   )                                        :: eps_int,u_reg
  integer , intent(in   ), dimension(3)                          :: n
  real(rp), intent(in   ), dimension(3)                          :: dli,dl
  integer , intent(in   )                                        :: nh_d,nh_u,nh_v,nh_p,nh_van
  real(rp), intent(in   ), dimension(1-nh_d:)                    :: dzc,dzf
  real(rp), intent(in   ), dimension(1-nh_u:,1-nh_u:,1-nh_u:)    :: u,v,w
  real(rp), intent(in   ), dimension(1-nh_v:,1-nh_v:,1-nh_v:)    :: cond,s,psi,ls
  real(rp), intent(inout), dimension(1-nh_v:,1-nh_v:,1-nh_v:)    :: tmp
  real(rp), intent(in   ), dimension(0:,0:,0:,1:)                :: nor
  real(rp), intent(inout), dimension(1-nh_v:,1-nh_v:,1-nh_v:)    :: ent  
  real(rp), intent(inout), dimension(     1:,     1:,     1:)    :: dentdtrko
  real(rp), intent(in   ), dimension(0:,0:,0:)                   :: rhocp
  integer , intent(in   ), dimension(3)                          :: halo_u,halo_van
  real(rp), intent(in   ), dimension(1:,1:)                      :: ent_bound_o
  !@cuf attributes(managed) :: velxc,velyc,velzc,weno2,dzc_van,dzf,ent_test   ! dl: plain array (see param.f90)
  real(rp), dimension(n(1),n(2),n(3)) :: dentdtrk
  integer  :: i,j,k
  real(rp), dimension(:,:,:), allocatable :: ent_test,weno2
  real(rp), dimension(:,:,:), allocatable :: velxc,velyc,velzc
  real(rp), dimension(:), allocatable     :: dzc_van
  !_______________________________
  ! Allocate arrays for transprt term calculate with WENO scheme
  allocate(weno2(1:n(1),1:n(2),1:n(3)))
  allocate(dzc_van( 1-nh_van:n(3)+nh_van))
  dzc_van(:) = dl(3)
  weno2(:,:,:) = 0.0
  ! Create auxiliary array and copy temperature field on intermediate array
  allocate(ent_test(1-nh_van:n(1)+nh_van,1-nh_van:n(2)+nh_van,1-nh_van:n(3)+nh_van))
  ent_test(:,:,:) = 0.0_rp
  do k=1,n(3)
    do j=1,n(2)
      do i=1,n(1)
        ent_test(i,j,k) = ent(i,j,k)
      enddo
    enddo
  enddo
  ! MPI communicate 2 halo points for WENO scheme
  call boundp(cbctmp,n,bctmp,nh_van,nh_van,halo_van,dl,dzc_van,dzf,ent_test)
  call boundsb_arrayo(n,nh_van,halo_van,ent_test,ent_bound_o)
  !
  !_______________________________
  !
  ! 1. transport by velocities + regularisation terms
  !
  allocate(velxc(1-nh_u:n(1)+nh_u,1-nh_u:n(2)+nh_u,1-nh_u:n(3)+nh_u))
  allocate(velyc(1-nh_u:n(1)+nh_u,1-nh_u:n(2)+nh_u,1-nh_u:n(3)+nh_u))
  allocate(velzc(1-nh_u:n(1)+nh_u,1-nh_u:n(2)+nh_u,1-nh_u:n(3)+nh_u))
  velxc(:,:,:) = 0.0
  velyc(:,:,:) = 0.0
  velzc(:,:,:) = 0.0 
  do k=1,n(3)
    do j=1,n(2)
      do i=1,n(1)
        velxc(i,j,k)= rhocp(i,j,k)*0.5*(u(i,j,k)+u(i-1,j  ,k  )) - (rho1*cp1-rho2*cp2)*0.5*(rex(i,j,k)+rex(i-1,j  ,k  ))
        velyc(i,j,k)= rhocp(i,j,k)*0.5*(v(i,j,k)+v(i  ,j-1,k  )) - (rho1*cp1-rho2*cp2)*0.5*(rey(i,j,k)+rey(i  ,j-1,k  ))
        velzc(i,j,k)= rhocp(i,j,k)*0.5*(w(i,j,k)+w(i  ,j  ,k-1)) - (rho1*cp1-rho2*cp2)*0.5*(rez(i,j,k)+rez(i  ,j  ,k-1))
      enddo
    enddo
  enddo
  call weno5(n(1),n(2),n(3),dli(1),dli(2),dli(3),nh_u,.FALSE.,ent_test,velxc,velyc,velzc,weno2) 
  deallocate(velxc,velyc,velzc,ent_test,dzc_van)
  !
  !_______________________________
  ! 3. transport by diffusion and phenomenological current (anti trapping current) 
  dentdtrk(:,:,:) = 0.0
  call diff_ent_ibm(n,dli,nh_d,nh_v,dzc,dzf,cond,ent,dentdtrk,rhocp,s,psi,nor)
  !_______________________________
  !
  do k=1,n(3)
    do j=1,n(2)
      do i=1,n(1)
        !
        ! Assemble RHS
        dentdtrk(i,j,k)  = dentdtrk(i,j,k) + weno2(i,j,k)/rhocp(i,j,k)
        ! Update temperature field
        ent(i,j,k)       = ent(i,j,k) + f_t1*dentdtrk(i,j,k) + f_t2*dentdtrko(i,j,k)
        dentdtrko(i,j,k) = dentdtrk(i,j,k)
        !
      enddo
    enddo
  enddo
  deallocate(weno2)
  !
  return
end subroutine rk_ent_ibm
  !
    subroutine diff_ent(n,dli,nh_d,nh_v,dzc,dzf,cond,s,dsdt) 
    !
      use mod_param     , only: small
    implicit none
    !
    integer , intent(in ), dimension(3)               :: n
    real(rp), intent(in ), dimension(3)               :: dli
    integer , intent(in )                             :: nh_d,nh_v
    real(rp), intent(in ), dimension(1-nh_d:)         :: dzc,dzf
    real(rp), intent(in ), dimension(1-nh_v:,1-nh_v:,1-nh_v:)    :: s,cond
    real(rp), intent(out), dimension(1:,1:,1:)                   :: dsdt
    !
    real(rp), dimension(3) :: dl
    real(rp) :: sfxp,sfxm,sfyp,sfym,sfzp,sfzm
    integer  :: i,j,k,ip,im,jp,jm,kp,km
    !
    dl(:) = 1.d0/dli(:)
    dsdt(:,:,:) = 0.d0
    !
    !
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          ip = i+1
          im = i-1
          jp = j+1
          jm = j-1
          kp = k+1
          km = k-1
          !
          ! along x-dir
          !
          !
          sfxp = (cond(ip,j,k)+cond(i,j,k))*(s(ip,j,k)-s(i,j,k))*dli(1)*0.5
          sfxm = (cond(im,j,k)+cond(i,j,k))*(s(i,j,k)-s(im,j,k))*dli(1)*0.5
          !
          dsdt(i,j,k) = (sfxp-sfxm)*dli(1)
          !
          !
          ! along y-dir
          !
          sfyp = (cond(i,jp,k)+cond(i,j,k))*(s(i,jp,k)-s(i,j,k))*dli(2)*0.5
          sfym = (cond(i,jm,k)+cond(i,j,k))*(s(i,j,k)-s(i,jm,k))*dli(2)*0.5
          !
          dsdt(i,j,k) = dsdt(i,j,k) + (sfyp-sfym)*dli(2) 
          !
          !
          ! along z-dir
          !
          sfzp = (cond(i,j,kp)+cond(i,j,k))*(s(i,j,kp)-s(i,j,k))*dli(3)*0.5
          sfzm = (cond(i,j,km)+cond(i,j,k))*(s(i,j,k)-s(i,j,km))*dli(3)*0.5
          !
          !
          dsdt(i,j,k)  = dsdt(i,j,k)+ (sfzp-sfzm)*dli(3)
          !
          !
        enddo
      enddo
    enddo
    !
    return
  end subroutine diff_ent


  subroutine diff_ent_ibm(n,dli,nh_d,nh_v,dzc,dzf,cond,s,dsdt,rhocp,src,psi,nor) 
    ! corrections required for stability reasons because of jumps on the conductivity at the immersed surface
    use mod_param       ,only: cond1,cond2,eps_int,IBM_h,h01,h02,cp1,cp2,heater_power
    use mod_boiling     ,only: heater_ind
    !
    implicit none
    !
    integer , intent(in ), dimension(3)                               :: n
    real(rp), intent(in ), dimension(3)                               :: dli
    integer , intent(in )                                             :: nh_d,nh_v
    real(rp), intent(in ), dimension(1-nh_d:)                         :: dzc,dzf
    real(rp), intent(in ), dimension(1-nh_v:,1-nh_v:,1-nh_v:)         :: s,cond
    real(rp), intent(out), dimension(1:,1:,1:)                        :: dsdt
    real(rp), intent(in ), dimension(0:,0:,0:,1:)                     :: nor
    real(rp), intent(in ), dimension(0:,0:,0:)                        :: rhocp,src,psi
    !
    real(rp), dimension(3) :: dl
    real(rp) :: sfxp,sfxm,sfyp,sfym,sfzp,sfzm,dsdt_cache,T_fake,maxmax,T_wall
    real(rp) :: crxp,crxm,cryp,crym,crzp,crzm,zm_switch
    integer  :: i,j,k,ip,im,jp,jm,kp,km,skip_layer
    real(rp) :: psip,psim,srcp,srcm,condp,condm,aux1,aux2,aux3,aux4,norm,norp,tmpp,tmpm,const_part
    real(rp), parameter ::     phi_param = 0.0d0
    !
    dl(:) = 1.d0/dli(:)
    dsdt(:,:,:) = 0.d0
    maxmax = - 1e6_rp
    const_part = eps_int*((cond1/cond2-1.0))  !
    !
    !============================================================================================
    !1. Calcluate heat coduction in the IBM wall
    !============================================================================================
    !
    if(IBM_h.GE.2)then
      do k=1,IBM_h-1
        do j=1,n(2)
          do i=1,n(1)
            ip = i+1
            im = i-1
            jp = j+1
            jm = j-1
            kp = k+1
            km = k-1
            dsdt_cache = 0.0
            !
            ! along x-dir
            !
            sfxp = (cond(ip,j,k)+cond(i,j,k))*(s(ip,j,k)-s(i ,j,k))*dli(1)*0.5
            sfxm = (cond(im,j,k)+cond(i,j,k))*(s(i ,j,k)-s(im,j,k))*dli(1)*0.5
            dsdt_cache = dsdt_cache + (sfxp-sfxm)*dli(1)
            !
            !
            ! along y-dir
            !
            sfyp = (cond(i,jp,k)+cond(i,j,k))*(s(i,jp,k)-s(i,j ,k))*dli(2)*0.5
            sfym = (cond(i,jm,k)+cond(i,j,k))*(s(i,j ,k)-s(i,jm,k))*dli(2)*0.5
            dsdt_cache = dsdt_cache + (sfyp-sfym)*dli(2)
            !
            !
            ! along z-dir
            !
            sfzp = (cond(i,j,kp)+cond(i,j,k))*(s(i,j,kp)-s(i,j,k ))*dli(3)*0.5
            sfzm = (cond(i,j,km)+cond(i,j,k))*(s(i,j,k )-s(i,j,km))*dli(3)*0.5
            !
            dsdt_cache   = dsdt_cache + (sfzp-sfzm)*dli(3)
            dsdt(i,j,k)  = dsdt_cache/rhocp(i,j,k)
          enddo
        enddo
      enddo
    endif
    !
    !============================================================================================
    !2. Calcluate heat coduction in the first layer of solid points (wall temperature treatment)
    !============================================================================================
    !
    if(IBM_h.GE.1)then
      k=IBM_h
      do j=1,n(2)
        do i=1,n(1)
          ip = i+1
          im = i-1
          jp = j+1
          jm = j-1
          kp = k+1
          km = k-1
          dsdt_cache = 0.0
          T_fake = 0.0
          T_wall = 0.0
          !
          !
          ! along x-dir
          !
          sfxp = (cond(ip,j,k)+cond(i,j,k))*(s(ip,j,k)-s(i ,j,k))*dli(1)*0.5
          sfxm = (cond(im,j,k)+cond(i,j,k))*(s(i ,j,k)-s(im,j,k))*dli(1)*0.5
          dsdt_cache = dsdt_cache + (sfxp-sfxm)*dli(1)
          !
          !
          ! along y-dir
          !
          sfyp = (cond(i,jp,k)+cond(i,j,k))*(s(i,jp,k)-s(i,j ,k))*dli(2)*0.5
          sfym = (cond(i,jm,k)+cond(i,j,k))*(s(i,j ,k)-s(i,jm,k))*dli(2)*0.5
          dsdt_cache = dsdt_cache + (sfyp-sfym)*dli(2)
          !
          !
          ! along z-dir
          !
          ! Calculate the fiktive Temperature for the solid inside the fluid:
          T_wall = (cond(i,j,k)*s(i,j,k)+cond(i,j,kp)*s(i,j,kp)) &
                  /(cond(i,j,k)         +cond(i,j,kp)          ) &
                  +(heater_power*dl(3)*heater_ind(i,j))          &
                  /(2.0_rp*(cond(i,j,k)+cond(i,j,kp)))
          T_fake = 2.0*T_wall-s(i,j,k)
          ! Now replace the conductivity and T_fake:
          sfzp = (cond(i,j,k )+cond(i,j,k))*( T_fake -s(i,j,k))*dli(3)*0.5
          sfzm = (cond(i,j,km)+cond(i,j,k))*(s(i,j,k)-s(i,j,km))*dli(3)*0.5
          !
          maxmax = max(maxmax,T_wall)
          !
          dsdt_cache = dsdt_cache + (sfzp-sfzm)*dli(3)
          dsdt(i,j,k)  = dsdt_cache/rhocp(i,j,k)
        enddo
      enddo
      !
      call mpi_allreduce(MPI_IN_PLACE,maxmax ,1,MPI_REAL_RP,MPI_MAX,MPI_COMM_WORLD,ierr)
      if(myid.eq.0) print*, 'Max T_wall:    ', maxmax
      !
    endif
    !
    !============================================================================================
    !3. Calcluate heat coduction in the first layer of fluid points (wall temperature treatment)
    !============================================================================================
    !
    if(IBM_h.GE.1)then
      k=IBM_h+1
      do j=1,n(2)
        do i=1,n(1)
          ip = i+1
          im = i-1
          jp = j+1
          jm = j-1
          kp = k+1
          km = k-1
          dsdt_cache = 0.0
          T_fake = 0.0
          T_wall = 0.0
          !
          ! along x-dir
          !
          sfxp = (cond(ip,j,k)+cond(i,j,k))*(s(ip,j,k)-s(i,j,k))*dli(1)*0.5
          sfxm = (cond(im,j,k)+cond(i,j,k))*(s(i,j,k)-s(im,j,k))*dli(1)*0.5
          dsdt_cache = dsdt_cache + (sfxp-sfxm)*dli(1)
          !
          !
          ! along y-dir
          !
          sfyp = (cond(i,jp,k)+cond(i,j,k))*(s(i,jp,k)-s(i,j,k))*dli(2)*0.5
          sfym = (cond(i,jm,k)+cond(i,j,k))*(s(i,j,k)-s(i,jm,k))*dli(2)*0.5
          dsdt_cache = dsdt_cache + (sfyp-sfym)*dli(2)
          !
          !
          ! along z-dir
          !
          ! Calculate the fiktive Temperature for the fluid inside the solid:
          T_wall = (cond(i,j,k)*s(i,j,k)+cond(i,j,km)*s(i,j,km)) &
                  /(cond(i,j,k)         +cond(i,j,km)          ) &
                  +(heater_power*dl(3)*heater_ind(i,j))          &
                  /(2.0_rp*(cond(i,j,k)+cond(i,j,km)))
          T_fake = 2.0*T_wall-s(i,j,k)
          ! Now replace the conductivity and T_fake:
          sfzp = (cond(i,j,kp)+cond(i,j,k))*(s(i,j,kp)-s(i,j,k))*dli(3)*0.5
          sfzm = (cond(i,j,k )+cond(i,j,k))*(s(i,j,k) - T_fake )*dli(3)*0.5
          !
          dsdt_cache = dsdt_cache + (sfzp-sfzm)*dli(3)
          !
          dsdt(i,j,k)  = dsdt_cache/rhocp(i,j,k)
        enddo
      enddo
    endif
    !
    !======================================================================================================
    !4. Calcluate heat coduction in the rest of the fluid domain (i.e. the stadard for non IBM simulations)
    !======================================================================================================
    !
    skip_layer = 0
    if(IBM_h.GE.1) skip_layer = 1
    do k=(1+skip_layer+IBM_h),n(3)
      do j=1,n(2)
        do i=1,n(1)
          ip = i+1
          im = i-1
          jp = j+1
          jm = j-1
          kp = k+1
          km = k-1
          dsdt_cache = 0.0
          !
          ! along x-dir
          !
          sfxp = (cond(ip,j,k)+cond(i,j,k))*(s(ip,j,k)-s(i ,j,k))*dli(1)*0.5
          sfxm = (cond(im,j,k)+cond(i,j,k))*(s(i ,j,k)-s(im,j,k))*dli(1)*0.5
          dsdt_cache = dsdt_cache + (sfxp-sfxm)*dli(1)
          !
          !
          ! along y-dir
          !
          sfyp = (cond(i,jp,k)+cond(i,j,k))*(s(i,jp,k)-s(i,j ,k))*dli(2)*0.5
          sfym = (cond(i,jm,k)+cond(i,j,k))*(s(i,j ,k)-s(i,jm,k))*dli(2)*0.5
          dsdt_cache = dsdt_cache + (sfyp-sfym)*dli(2)
          !
          !
          ! along z-dir
          !
          sfzp = (cond(i,j,kp)+cond(i,j,k))*(s(i,j,kp)-s(i,j,k ))*dli(3)*0.5
          sfzm = (cond(i,j,km)+cond(i,j,k))*(s(i,j,k )-s(i,j,km))*dli(3)*0.5
          !
          dsdt_cache   = dsdt_cache + (sfzp-sfzm)*dli(3)
          dsdt(i,j,k)  = dsdt_cache/rhocp(i,j,k)
        enddo
      enddo
    enddo
    !
    !======================================================================================================
    !5. Calcluate the anti-trapping currents
    !======================================================================================================
    !
    do k=1+IBM_h,n(3)
      zm_switch = 1.0_rp
      if(k.EQ.(1+IBM_h)) zm_switch = 0.0_rp  ! Anti-trapping current is zero on wall -> should not be necessary since the interface normal vector is zero in the wall
      do j=1,n(2)
        do i=1,n(1)
          ip = i+1
          im = i-1
          jp = j+1
          jm = j-1
          kp = k+1
          km = k-1
          dsdt_cache = 0.0
          ! In X - direction
          !
          tmpp=0.5*(  s(i,j,k  )  +  s(i+1,j,k)  )
          srcp=0.5*(src(i,j,k  )  +src(i+1,j,k)  )
          norp=0.5*(nor(i,j,k  ,1)+nor(i+1,j,k,1))
          crxp= -norp*const_part*(h01-h02 + (cp1-cp2)*tmpp)*srcp
          !
          tmpm=0.5*(  s(i,j,k  )  +  s(i-1,j,k)  )
          srcm=0.5*(src(i,j,k  )  +src(i-1,j,k)  )
          norm=0.5*(nor(i,j,k  ,1)+nor(i-1,j,k,1))
          crxm= -norm*const_part*(h01-h02 + (cp1-cp2)*tmpm)*srcm
          !
          dsdt_cache = dsdt_cache + (crxp-crxm)*dli(1)
          ! In Y - direction
          !
          tmpp=0.5*(  s(i,j,k  )  +  s(i,j+1,k)  )
          srcp=0.5*(src(i,j,k  )  +src(i,j+1,k)  )
          norp=0.5*(nor(i,j,k  ,2)+nor(i,j+1,k,2))
          cryp= -norp*const_part*(h01-h02 + (cp1-cp2)*tmpp)*srcp
          !
          tmpm=0.5*(  s(i,j,k  )  +  s(i,j-1,k)  )
          srcm=0.5*(src(i,j,k  )  +src(i,j-1,k)  )
          norm=0.5*(nor(i,j,k  ,2)+nor(i,j-1,k,2))
          crym= -norm*const_part*(h01-h02 + (cp1-cp2)*tmpm)*srcm
          !
          dsdt_cache = dsdt_cache + (cryp-crym)*dli(2)
          ! In Z - direction
          tmpp=0.5*(  s(i,j,k  )  +  s(i,j,k+1)  )
          srcp=0.5*(src(i,j,k  )  +src(i,j,k+1)  )
          norp=0.5*(nor(i,j,k  ,3)+nor(i,j,k+1,3))
          crzp= -norp*const_part*(h01-h02 + (cp1-cp2)*tmpp)*srcp
          !
          tmpm=0.5*(  s(i,j,k  )  +  s(i,j,k-1)  )
          srcm=0.5*(src(i,j,k  )  +src(i,j,k-1)  )
          norm=0.5*(nor(i,j,k  ,3)+nor(i,j,k-1,3))
          crzm= -norm*const_part*(h01-h02 + (cp1-cp2)*tmpm)*srcm*zm_switch
          !
          dsdt_cache = dsdt_cache + (crzp-crzm)*dli(3)
          !
          dsdt(i,j,k)  = dsdt(i,j,k) + dsdt_cache/rhocp(i,j,k)
          !
          !
        enddo
      enddo
    enddo
    !
    return
  end subroutine diff_ent_ibm
  !
  subroutine static_contact_angle_m(n,dl,nh_v,eps,theta,ap,ap_bound)
    !                             
    ! To do: we can treat the term multiplying grad_rho implicitly
    ! 
    implicit none
    !
    integer , intent(in ), dimension(3) :: n
    real(rp), intent(in ), dimension(3) :: dl
    integer , intent(in )               :: nh_v
    real(rp), intent(in )               :: eps,theta
    real(rp), intent(in ), dimension(1-nh_v:,1-nh_v:,1-nh_v:) :: ap
    real(rp), intent(out), dimension( 1:, 1:)     :: ap_bound
    real(rp) :: f,mv_wall,ap_wall,tan_abs,tan_x,tan_y
    integer  :: i,j,k,kpp,kp,ip,im,jp,jm
    !
    do j=1,n(2)
      do i=1,n(1)
        ! written by Lorenz Weber 19th DEC 2023 (lorenzw@kth.se)
        ! Further reading in the Literature:
        ! - equation 36 and 37 in Zarareh (2021): Improving the staircase approximation for wettability implementation of phase-field model: Part 1 – Static contact angle
        ! - equation 42, 43 and 44 in Liang (2019): Lattice Boltzmann method for contact-line motion of binary fluids with high density ratio
        ! - equation 18 in Ding (2007): Wetting condition in diffuse interface simulations of contact line motion
        k   =   0
        kp  = k + 1
        kpp = k + 2
        ip  = i + 1
        im  = i - 1
        jp  = j + 1
        jm  = j - 1
        !
        tan_x         = 1.5d0*((ap(ip,j ,kp)-ap(im,j ,kp))/dl(1)/(2.0)) - 0.5*((ap(ip,j ,kpp)-ap(im,j ,kpp))/dl(1)/(2.0))
        tan_y         = 1.5d0*((ap(i ,jp,kp)-ap(i ,jm,kp))/dl(2)/(2.0)) - 0.5*((ap(i ,jp,kpp)-ap(i ,jm,kpp))/dl(2)/(2.0))
        tan_abs       = sqrt(tan_x**2+tan_y**2)
        ap_bound(i,j) = max(0.d0,min(1.d0,ap(i,j,kp) + tan_abs*TAND(90.d0 - theta)*dl(3)))
      enddo
    enddo
    return
  end subroutine static_contact_angle_m
  !
  subroutine cmpt_delta(n,dli,phi,delta)
    !
    implicit none
    !
    integer , intent(in ), dimension(3) :: n
    real(rp), intent(in ), dimension(3) :: dli
    real(rp), intent(in ), dimension(0:,0:,0:) :: phi
    real(rp), intent(out), dimension(0:,0:,0:) :: delta
    !
    real(rp), dimension(8) :: nx,ny,nz,mx,my,mz
    real(rp), dimension(3) :: dl
    real(rp) :: dphidx,dphidy,dphidz
    real(rp) :: small = 1.e-12
    integer  :: i,j,k,p
    ! 
    dl(:) = dli(:)**(-1)
    !
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          !
#ifdef TWOD 
          mx(1:8) = 0.d0
#else
          !i+1/2 j+1/2 k+1/2
          mx(1)=((phi(i+1,j  ,k  )+phi(i+1,j+1,k  )+phi(i+1,j  ,k+1)+phi(i+1,j+1,k+1))-&
                 (phi(i  ,j  ,k  )+phi(i  ,j+1,k  )+phi(i  ,j  ,k+1)+phi(i  ,j+1,k+1)))*dli(1)*0.25
          !i+1/2 j-1/2 k+1/2
          mx(2)=((phi(i+1,j  ,k  )+phi(i+1,j-1,k  )+phi(i+1,j  ,k+1)+phi(i+1,j-1,k+1))-&
                 (phi(i  ,j  ,k  )+phi(i  ,j-1,k  )+phi(i  ,j  ,k+1)+phi(i  ,j-1,k+1)))*dli(1)*0.25
          !i+1/2 j+1/2 k-1/2
          mx(3)=((phi(i+1,j  ,k  )+phi(i+1,j+1,k  )+phi(i+1,j  ,k-1)+phi(i+1,j+1,k-1))-&
                 (phi(i  ,j  ,k  )+phi(i  ,j+1,k  )+phi(i  ,j  ,k-1)+phi(i  ,j+1,k-1)))*dli(1)*0.25
          !i+1/2 j-1/2 k-1/2
          mx(4)=((phi(i+1,j  ,k  )+phi(i+1,j-1,k  )+phi(i+1,j  ,k-1)+phi(i+1,j-1,k-1))-&
                 (phi(i  ,j  ,k  )+phi(i  ,j-1,k  )+phi(i  ,j  ,k-1)+phi(i  ,j-1,k-1)))*dli(1)*0.25
          !i-1/2 j+1/2 k+1/2
          mx(5)=((phi(i  ,j  ,k  )+phi(i  ,j+1,k  )+phi(i  ,j  ,k+1)+phi(i  ,j+1,k+1))-&
                 (phi(i-1,j  ,k  )+phi(i-1,j+1,k  )+phi(i-1,j  ,k+1)+phi(i-1,j+1,k+1)))*dli(1)*0.25
          !i-1/2 j-1/2 k+1/2
          mx(6)=((phi(i  ,j  ,k  )+phi(i  ,j-1,k  )+phi(i  ,j  ,k+1)+phi(i  ,j-1,k+1))-&
                 (phi(i-1,j  ,k  )+phi(i-1,j-1,k  )+phi(i-1,j  ,k+1)+phi(i-1,j-1,k+1)))*dli(1)*0.25
          !i-1/2 j+1/2 k-1/2
          mx(7)=((phi(i  ,j  ,k  )+phi(i  ,j+1,k  )+phi(i  ,j  ,k-1)+phi(i  ,j+1,k-1))-&
                 (phi(i-1,j  ,k  )+phi(i-1,j+1,k  )+phi(i-1,j  ,k-1)+phi(i-1,j+1,k-1)))*dli(1)*0.25
          !i-1/2 j-1/2 k-1/2
          mx(8)=((phi(i  ,j  ,k  )+phi(i  ,j-1,k  )+phi(i  ,j  ,k-1)+phi(i  ,j-1,k-1))-&
                 (phi(i-1,j  ,k  )+phi(i-1,j-1,k  )+phi(i-1,j  ,k-1)+phi(i-1,j-1,k-1)))*dli(1)*0.25
          !
#endif 
         my(1)=((phi(i  ,j+1,k  )+phi(i+1,j+1,k  )+phi(i  ,j+1,k+1)+phi(i+1,j+1,k+1))-&
                 (phi(i  ,j  ,k  )+phi(i+1,j  ,k  )+phi(i  ,j  ,k+1)+phi(i+1,j  ,k+1)))*dli(2)*0.25
          !i+1/2 j-1/2 k+1/2
          my(2)=((phi(i  ,j  ,k  )+phi(i+1,j  ,k  )+phi(i  ,j  ,k+1)+phi(i+1,j  ,k+1))-&
                 (phi(i  ,j-1,k  )+phi(i+1,j-1,k  )+phi(i  ,j-1,k+1)+phi(i+1,j-1,k+1)))*dli(2)*0.25
          !i+1/2 j+1/2 k-1/2
          my(3)=((phi(i  ,j+1,k  )+phi(i+1,j+1,k  )+phi(i  ,j+1,k-1)+phi(i+1,j+1,k-1))-&
                 (phi(i  ,j  ,k  )+phi(i+1,j  ,k  )+phi(i  ,j  ,k-1)+phi(i+1,j  ,k-1)))*dli(2)*0.25
          !i+1/2 j-1/2 k-1/2
          my(4)=((phi(i  ,j  ,k  )+phi(i+1,j  ,k  )+phi(i  ,j  ,k-1)+phi(i+1,j  ,k-1))-&
                 (phi(i  ,j-1,k  )+phi(i+1,j-1,k  )+phi(i  ,j-1,k-1)+phi(i+1,j-1,k-1)))*dli(2)*0.25
          !i-1/2 j+1/2 k+1/2
          my(5)=((phi(i  ,j+1,k  )+phi(i-1,j+1,k  )+phi(i  ,j+1,k+1)+phi(i-1,j+1,k+1))-&
                 (phi(i  ,j  ,k  )+phi(i-1,j  ,k  )+phi(i  ,j  ,k+1)+phi(i-1,j  ,k+1)))*dli(2)*0.25
          !i-1/2 j-1/2 k+1/2
          my(6)=((phi(i  ,j  ,k  )+phi(i-1,j  ,k  )+phi(i  ,j  ,k+1)+phi(i-1,j  ,k+1))-&
                 (phi(i  ,j-1,k  )+phi(i-1,j-1,k  )+phi(i  ,j-1,k+1)+phi(i-1,j-1,k+1)))*dli(2)*0.25
          !i-1/2 j+1/2 k-1/2
          my(7)=((phi(i  ,j+1,k  )+phi(i-1,j+1,k  )+phi(i  ,j+1,k-1)+phi(i-1,j+1,k-1))-&
                 (phi(i  ,j  ,k  )+phi(i-1,j  ,k  )+phi(i  ,j  ,k-1)+phi(i-1,j  ,k-1)))*dli(2)*0.25
          !i-1/2 j-1/2 k-1/2
          my(8)=((phi(i  ,j  ,k  )+phi(i-1,j  ,k  )+phi(i  ,j  ,k-1)+phi(i-1,j  ,k-1))-&
                 (phi(i  ,j-1,k  )+phi(i-1,j-1,k  )+phi(i  ,j-1,k-1)+phi(i-1,j-1,k-1)))*dli(2)*0.25
          !
          !i+1/2 j+1/2 k+1/2
          mz(1)=((phi(i  ,j  ,k+1)+phi(i+1,j  ,k+1)+phi(i  ,j+1,k+1)+phi(i+1,j+1,k+1))-&
                 (phi(i  ,j  ,k  )+phi(i+1,j  ,k  )+phi(i  ,j+1,k  )+phi(i+1,j+1,k  )))*dli(3)*0.25
          !i+1/2 j-1/2 k+1/2
          mz(2)=((phi(i  ,j  ,k+1)+phi(i+1,j  ,k+1)+phi(i  ,j-1,k+1)+phi(i+1,j-1,k+1))-&
                 (phi(i  ,j  ,k  )+phi(i+1,j  ,k  )+phi(i  ,j-1,k  )+phi(i+1,j-1,k  )))*dli(3)*0.25
          !i+1/2 j+1/2 k-1/2
          mz(3)=((phi(i  ,j  ,k  )+phi(i+1,j  ,k  )+phi(i  ,j+1,k  )+phi(i+1,j+1,k  ))-&
                 (phi(i  ,j  ,k-1)+phi(i+1,j  ,k-1)+phi(i  ,j+1,k-1)+phi(i+1,j+1,k-1)))*dli(3)*0.25
          !i+1/2 j-1/2 k-1/2
          mz(4)=((phi(i  ,j  ,k  )+phi(i+1,j  ,k  )+phi(i  ,j-1,k  )+phi(i+1,j-1,k  ))-&
                 (phi(i  ,j  ,k-1)+phi(i+1,j  ,k-1)+phi(i  ,j-1,k-1)+phi(i+1,j-1,k-1)))*dli(3)*0.25
          !i-1/2 j+1/2 k+1/2
          mz(5)=((phi(i  ,j  ,k+1)+phi(i-1,j  ,k+1)+phi(i  ,j+1,k+1)+phi(i-1,j+1,k+1))-&
                 (phi(i  ,j  ,k  )+phi(i-1,j  ,k  )+phi(i  ,j+1,k  )+phi(i-1,j+1,k  )))*dli(3)*0.25
          !i-1/2 j-1/2 k+1/2
          mz(6)=((phi(i  ,j  ,k+1)+phi(i-1,j  ,k+1)+phi(i  ,j-1,k+1)+phi(i-1,j-1,k+1))-&
                 (phi(i  ,j  ,k  )+phi(i-1,j  ,k  )+phi(i  ,j-1,k  )+phi(i-1,j-1,k  )))*dli(3)*0.25
          !i-1/2 j+1/2 k-1/2
          mz(7)=((phi(i  ,j  ,k  )+phi(i-1,j  ,k  )+phi(i  ,j+1,k  )+phi(i-1,j+1,k  ))-&
                 (phi(i  ,j  ,k-1)+phi(i-1,j  ,k-1)+phi(i  ,j+1,k-1)+phi(i-1,j+1,k-1)))*dli(3)*0.25
          !i-1/2 j-1/2 k-1/2
          mz(8)=((phi(i  ,j  ,k  )+phi(i-1,j  ,k  )+phi(i  ,j-1,k  )+phi(i-1,j-1,k  ))-&
                 (phi(i  ,j  ,k-1)+phi(i-1,j  ,k-1)+phi(i  ,j-1,k-1)+phi(i-1,j-1,k-1)))*dli(3)*0.25
          !
          ! compute the delta
          !
          dphidx = 0._rp
          dphidy = 0._rp
          dphidz = 0._rp
          do p=1,8
            dphidx = dphidx + 0.125_rp*mx(p)
            dphidy = dphidy + 0.125_rp*my(p)
            dphidz = dphidz + 0.125_rp*mz(p)
          enddo
          !
          delta(i,j,k) = sqrt(dphidx**2+dphidy**2+dphidz**2+small)
          !
        enddo
      enddo
    enddo
    !
    return
  end subroutine cmpt_delta
    !
    !
  subroutine marangoni_he(n,dli,dzc,dzf,tmp,nh_t,delta_arr,kappa,nor,a1,rho,mar_force)
    !
    use mod_param, only: sigma_t,sigma,tl0,tg0
    !
    implicit none
    !
    integer , intent(in ), dimension(3)                         :: n
    real(rp), intent(in ), dimension(3)                         :: dli
    real(rp), intent(in ), dimension(-2:)                       :: dzc,dzf
    integer , intent(in   )                                     :: nh_t
    real(rp), intent(in ), dimension(1-nh_t:,1-nh_t:,1-nh_t:)   :: tmp
    real(rp), intent(in ), dimension( 0:, 0:, 0:)               :: delta_arr,kappa
    real(rp), intent(in ), dimension(0:,0:,0:,1:)               :: nor
    real(rp), intent(in ), dimension(0:,0:,0:)                  :: rho
    real(rp), intent(in ), dimension(0:,0:,0:)                  :: a1
    real(rp), intent(out), dimension(-2:,-2:,-2:,1:)            :: mar_force
    !
    real(rp) :: tmpf,kappas
    real(rp) :: norx,nory,norz
    real(rp) :: dpsidx,dpsidy,dpsidz
    real(rp) :: dtmpdx,dtmpdy,dtmpdz,dt_nof
    real(rp) :: rhox, rhoy, rhoz, rhoxi, rhoyi, rhozi
    real(rp) :: tmp_min,tmp_max,sigma_min,sigma_max
    real(rp) :: am, ap, delta, deltas
    real(rp), dimension(3) :: dl
    real(rp), dimension(n(1),n(2),n(3)) :: sigma_arr
    real(rp):: tmp1
    integer  :: i,j,k,p
    !
    tmp1=tl0
    dl(:) = 1.d0/dli(:)
    !
    do k=1,n(3)
      do j=1,n(2)
         do i=1,n(1)
           !
           !
#ifdef TWOD
           mar_force(i,j,k,1) = 0.d0
#else
           !
           ! X-direction (three contributions)
           !
           tmpf   = 0.5*(tmp(  i+1,j,k  )+tmp(  i,j,k  ))
           kappas = 0.5*(kappa(i+1,j,k  )+kappa(i,j,k  ))
           norx   = 0.5*(nor(  i+1,j,k,1)+nor(  i,j,k,1))
           nory   = 0.5*(nor(  i+1,j,k,2)+nor(  i,j,k,2))
           norz   = 0.5*(nor(  i+1,j,k,3)+nor(  i,j,k,3))
           dpsidx = (a1( i+1,j,k)-a1( i,j,k))*dli(1)
           !
           rhox   = 0.5d0*(rho(i+1,j,k)+rho(i,j,k))
           rhoxi  = 1.d0 / rhox
           deltas = 0.5d0*(delta_arr(i+1,j,k)+delta_arr(i,j,k))
           !
           dtmpdx = (tmp(i+1,j,k)-tmp(i,j,k))*dli(1)
           dtmpdy = 0.25*( (tmp(i,j,k)+tmp(i+1,j,k)+tmp(i+1,j+1,k)+tmp(i,j+1,k)) - &
                           (tmp(i,j,k)+tmp(i+1,j,k)+tmp(i+1,j-1,k)+tmp(i,j-1,k)) )*dli(2)
           dtmpdz = 0.25*( (tmp(i,j,k)+tmp(i+1,j,k)+tmp(i+1,j,k+1)+tmp(i,j,k+1)) - &
                           (tmp(i,j,k)+tmp(i+1,j,k)+tmp(i+1,j,k-1)+tmp(i,j,k-1)) )*dli(3)
           dt_nof = dtmpdx*norx+dtmpdy*nory+dtmpdz*norz
           !
           ap     =  a1(i+1,j,k)
           am     =  a1(i  ,j,k)
           delta  =  ((ap-am)*dli(1))**2
           ap     =  0.25*(a1(i,j,k)+a1(i+1,j,k)+a1(i+1,j+1,k)+a1(i,j+1,k))
           am     =  0.25*(a1(i,j,k)+a1(i+1,j,k)+a1(i+1,j-1,k)+a1(i,j-1,k))
           delta  =  delta+((ap-am)*dli(2))**2
           ap     =  0.25*(a1(i,j,k)+a1(i+1,j,k)+a1(i+1,j,k+1)+a1(i,j,k+1))
           am     =  0.25*(a1(i,j,k)+a1(i+1,j,k)+a1(i+1,j,k-1)+a1(i,j,k-1))
           delta  =  sqrt(delta+((ap-am)*dli(3))**2)
           !
           !mar_force(i,j,k,1) = sigma_t*( (tmpf-tmp1)*kappas*dpsidx + dtmpdx*delta - dt_nof*dpsidx )*rhoxi
           mar_force(i,j,k,1) = sigma_t*( (tmpf-tmp1)*kappas*dpsidx + dtmpdx*deltas - dt_nof*dpsidx )*rhoxi
#endif
           !
           !
           ! Y-direction (three contribution)
           !
           tmpf   = 0.5*(tmp(  i,j+1,k  )+tmp(  i,j,k))
           kappas = 0.5*(kappa(i,j+1,k  )+kappa(i,j,k))
           norx   = 0.5*(nor(  i,j+1,k,1)+nor(  i,j,k,1))
           nory   = 0.5*(nor(  i,j+1,k,2)+nor(  i,j,k,2))
           norz   = 0.5*(nor(  i,j+1,k,3)+nor(  i,j,k,3))
           dpsidy = (a1( i,j+1,k)-a1( i,j,k))*dli(1)
           !
           rhoy   = 0.5d0*(rho(i,j+1,k)+rho(i,j,k))
           rhoyi  = 1.d0 / rhoy
           deltas = 0.5d0*(delta_arr(i,j+1,k)+delta_arr(i,j,k))
           !
           dtmpdx = 0.25*( (tmp(i,j,k)+tmp(i+1,j,k)+tmp(i+1,j+1,k)+tmp(i,j+1,k)) - &
                           (tmp(i,j,k)+tmp(i-1,j,k)+tmp(i-1,j+1,k)+tmp(i,j+1,k)) )*dli(1)
           dtmpdy = (tmp(i,j+1,k)-tmp(i,j,k))*dli(2)
           dtmpdz = 0.25*( (tmp(i,j,k)+tmp(i,j+1,k)+tmp(i,j+1,k+1)+tmp(i,j,k+1)) - &
                           (tmp(i,j,k)+tmp(i,j+1,k)+tmp(i,j+1,k+1)+tmp(i,j,k+1)) )*dli(3)
           dt_nof = dtmpdx*norx+dtmpdy*nory+dtmpdz*norz
           !
           ap     =  0.25*(a1(i,j,k)+a1(i+1,j,k)+a1(i+1,j+1,k)+a1(i,j+1,k))
           am     =  0.25*(a1(i,j,k)+a1(i-1,j,k)+a1(i-1,j+1,k)+a1(i,j+1,k))
           delta  =  ((ap-am)*dli(1))**2
           ap     =  a1(i,j+1,k)
           am     =  a1(i,j  ,k)
           delta  =  delta+((ap-am)*dli(2))**2
           ap     =  0.25*(a1(i,j,k)+a1(i,j+1,k)+a1(i,j+1,k+1)+a1(i,j,k+1))
           am     =  0.25*(a1(i,j,k)+a1(i,j+1,k)+a1(i,j+1,k-1)+a1(i,j,k-1))
           delta  =  sqrt(delta+((ap-am)*dli(3))**2)
           !
           !
           !mar_force(i,j,k,2) = sigma_t*( (tmpf-tmp1)*kappas*dpsidy + dtmpdy*delta - dt_nof*dpsidy )*rhoyi
           mar_force(i,j,k,2) = sigma_t*( (tmpf-tmp1)*kappas*dpsidy + dtmpdz*deltas - dt_nof*dpsidy )*rhoyi
           !
           ! Z-contribution (three contribution)
           !
           tmpf   = 0.5*(tmp(  i,j,k+1  )+tmp(  i,j,k))
           kappas = 0.5*(kappa(i,j,k+1  )+kappa(i,j,k))
           norx   = 0.5*(nor(  i,j,k+1,1)+nor(  i,j,k,1))
           nory   = 0.5*(nor(  i,j,k+1,2)+nor(  i,j,k,2))
           norz   = 0.5*(nor(  i,j,k+1,3)+nor(  i,j,k,3))
           dpsidz = (a1( i,j,k+1)-a1( i,j,k))*dli(1)
           !
           rhoz   = 0.5d0*(rho(i,j,k+1)+rho(i,j,k))
           rhozi  = 1.d0 / rhoz
           deltas = 0.5d0*(delta_arr(i,j,k+1)+delta_arr(i,j,k))
           !
           dtmpdx = 0.25*( (tmp(i,j,k)+tmp(i+1,j,k)+tmp(i+1,j,k+1)+tmp(i,j,k+1)) - &
                           (tmp(i,j,k)+tmp(i-1,j,k)+tmp(i-1,j,k+1)+tmp(i,j,k+1)) )*dli(1)
           dtmpdy = 0.25*( (tmp(i,j,k)+tmp(i,j+1,k)+tmp(i,j+1,k+1)+tmp(i,j,k+1)) - &
                           (tmp(i,j,k)+tmp(i,j-1,k)+tmp(i,j-1,k+1)+tmp(i,j,k+1)) )*dli(2)
           dtmpdz = (tmp(i,j,k+1)-tmp(i,j,k))*dli(3)
           dt_nof = dtmpdx*norx+dtmpdy*nory+dtmpdz*norz
           !
           ap     =  0.25*(a1(i,j,k)+a1(i+1,j,k)+a1(i+1,j,k+1)+a1(i,j,k+1))
           am     =  0.25*(a1(i,j,k)+a1(i-1,j,k)+a1(i-1,j,k+1)+a1(i,j,k+1))
           delta  =  ((ap-am)*dli(1))**2
           ap     =  0.25*(a1(i,j,k)+a1(i,j+1,k)+a1(i,j+1,k+1)+a1(i,j,k+1))
           am     =  0.25*(a1(i,j,k)+a1(i,j-1,k)+a1(i,j-1,k+1)+a1(i,j,k+1))
           delta  =  delta+((ap-am)*dli(2))**2
           ap     =  a1(i,j,k+1)
           am     =  a1(i,j,k)
           delta  =  sqrt(delta+((ap-am)*dli(3))**2)
           !
           !mar_force(i,j,k,3) = sigma_t*( (tmpf-tmp1)*kappas*dpsidz + dtmpdz*delta - dt_nof*dpsidz )*rhozi
           mar_force(i,j,k,3) = sigma_t*( (tmpf-tmp1)*kappas*dpsidz + dtmpdz*deltas - dt_nof*dpsidz )*rhozi
           !
         enddo
      enddo
    enddo
    !
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          sigma_arr(i,j,k) = sigma + sigma_t*(tmp(i,j,k)-tmp1)
        enddo
      enddo
    enddo
    !
    sigma_min = minval(sigma_arr(1:n(1),1:n(2),1:n(3)))
    call mpi_allreduce(MPI_IN_PLACE,sigma_min,1,mpi_real8,mpi_min,comm_cart,ierr)
    sigma_max = maxval(sigma_arr(1:n(1),1:n(2),1:n(3)))
    call mpi_allreduce(MPI_IN_PLACE,sigma_max,1,mpi_real8,mpi_max,comm_cart,ierr)
    !
    if(myid.eq.0) print*, "sigma_min", sigma_min
    if(myid.eq.0) print*, "sigma_max", sigma_max
    !
    tmp_min = minval(tmp(1:n(1),1:n(2),1:n(3)))
    call mpi_allreduce(MPI_IN_PLACE,tmp_min,1,mpi_real8,mpi_min,comm_cart,ierr)
    tmp_max = maxval(tmp(1:n(1),1:n(2),1:n(3)))
    call mpi_allreduce(MPI_IN_PLACE,tmp_max,1,mpi_real8,mpi_max,comm_cart,ierr)
    !
    if(myid.eq.0) print*, "tmp_min", tmp_min
    if(myid.eq.0) print*, "tmp_max", tmp_max
    !
    return
  end subroutine marangoni_he
    !
    !
    !
    subroutine rk_c1(f_t1,f_t2,n,dli,nh_d,nh_u,nh_v, &
                      dzc,dzf,diff1,diff2,diff,u_reg,eps_int,u,v,w,ur,vr,wr,c1,c2,src,phi,nor,dc1dtrko,dc2dtrko,id)     
                                                                          !phi is order paramter
                                                                          !c1 is scalar
    !
    ! low-storage 3rd-order Runge-Kutta scheme 
    ! for time integration of the c1 field.
    !
    implicit none
    !
    real(rp), intent(in   )                                        :: f_t1,f_t2
    integer , intent(in   ), dimension(3)                          :: n
    real(rp), intent(in   ), dimension(3)                          :: dli
    integer , intent(in   )                                        :: nh_d,nh_u,nh_v
    real(rp), intent(in   ), dimension(1-nh_d:)                    :: dzc,dzf
    real(rp), intent(in   )                                        :: diff1,diff2,eps_int,u_reg
    real(rp), intent(in   ), dimension(1-nh_u:,1-nh_u:,1-nh_u:)    :: u,v,w,ur,vr,wr
    real(rp), intent(inout), dimension(1-nh_v:,1-nh_v:,1-nh_v:)    :: c1,c2  ! concentraion at current time level
    real(rp), intent(in   ), dimension(1-nh_v:,1-nh_v:,1-nh_v:)    :: phi,src,diff ! level-set,sai
    real(rp), intent(in   ), dimension(1-nh_v:,1-nh_v:,1-nh_v:,1:) :: nor
    real(rp), intent(inout), dimension(     1:,     1:,     1:)    :: dc1dtrko,dc2dtrko
    integer , intent(in   )                                        :: id
    ! 
    real(rp), dimension(n(1),n(2),n(3)) :: dc1dtrk,dc2dtrk
    real(rp) :: factor1,factor2,factor3,tot_reg,int_den,q_t,tot_lag
    integer  :: i,j,k
    !
    factor1 = f_t1!rkpar(1)*dt
    factor2 = f_t2!rkpar(2)*dt
    factor3 = f_t1+f_t2!factor1+factor2
    !
    dc1dtrk = 0.d0
    dc2dtrk = 0.d0
    !
    ! 1. compute regularization
    !
    !
    if (id.eq.1)then
    call regularization_c1(n,dli,nh_d,nh_v,dzc,dzf,diff1,u_reg,eps_int,c1,Phi,nor,dc1dtrk)
    elseif(id.eq.-1)then 
    call regularization_c1(n,dli,nh_d,nh_v,dzc,dzf,diff2,u_reg,eps_int,c2,1.d0-Phi,-1.d0*nor,dc2dtrk)
    endif 
    !
    ! 2. compute interfacial scalar transfer
    !
!    if (id.eq.1)then 
!    call interfacial_cflux(n,dli,nh_d,nh_v,dzc,dzf,diff1,diff2,c1,c2,phi,dc1dtrk,id)
!    elseif(id.eq.-1)then 
!    call interfacial_cflux(n,dli,nh_d,nh_v,dzc,dzf,diff1,diff2,c1,c2,phi,dc2dtrk,id)
!    endif
    !
    ! 3. Add advection 
    ! 
    if (id.eq.1)then 
    call moms_c1(n,dli,nh_d,nh_u,nh_v,u,v,w,ur,vr,wr,c1,Phi,nor,dc1dtrk)
    elseif(id.eq.-1)then 
    call moms_c1(n,dli,nh_d,nh_u,nh_v,u,v,w,ur,vr,wr,c2,1.d0-Phi,-1.d0*nor,dc2dtrk)
    !
   do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
        dc2dtrk(i,j,k) = dc2dtrk(i,j,k) + src(i,j,k)
      enddo
     enddo
    enddo
   endif   
    !
    !
    !
   tot_reg = 0._rp
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          tot_reg = tot_reg + dc1dtrk(i,j,k)/(dli(1)*dli(2)*dli(3))
        enddo
      enddo
    enddo
    call mpi_allreduce(MPI_IN_PLACE,tot_reg,1,MPI_REAL_RP,MPI_SUM,comm_cart,ierr)
    if(myid.eq.0) print*, "Tot. reg.", tot_reg
    !
#if defined(_CONTACT_LINE)
    int_den = 0._rp
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          int_den = int_den + 6._rp*c1(i,j,k)*(1._rp-c1(i,j,k))/(dli(1)*dli(2)*dli(3))
        enddo
      enddo
    enddo
    call mpi_allreduce(MPI_IN_PLACE,int_den,1,MPI_REAL_RP,MPI_SUM,comm_cart,ierr)
    q_t = tot_reg/int_den
#endif
    !
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(i,j,k) &
    !$OMP SHARED(n,f_t1,f_t2,psi,dpsidtrk,dpsidtrko)
    !
    !
    !
    !
    !
    !
    !
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          !
#if defined(_CONTACT_LINE)
          tot_lag          = 6._rp*c1(i,j,k)*(1._rp-c1(i,j,k))*q_t
          c1(i,j,k)        = c1(i,j,k) + f_t1*( dc1dtrk(i,j,k)-tot_lag) + f_t2*dc1dtrko(i,j,k)
          dc1dtrko(i,j,k)  = dc1dtrk(i,j,k)-tot_lag
#else
   if(id.eq.1)then
          c1(i,j,k)        = c1(i,j,k) + f_t1*dc1dtrk(i,j,k) + f_t2*dc1dtrko(i,j,k)
          dc1dtrko(i,j,k)  = dc1dtrk(i,j,k)
   elseif(id.eq.-1) then 
          c2(i,j,k)        = c2(i,j,k) + f_t1*dc2dtrk(i,j,k) + f_t2*dc2dtrko(i,j,k)
          dc2dtrko(i,j,k)  = dc2dtrk(i,j,k)
   endif
#endif
          !
        enddo
      enddo
    enddo
    !$OMP END PARALLEL DO
    return
  !
  end subroutine rk_c1
  !
  subroutine rk_psi(f_t1,f_t2,n,dli,nh_d,nh_u,nh_v, &
                    dzc,dzf,u_reg,eps_int,u,v,w,psi,src,ls,nor,kappa,dpsidtrko,tmp,ent,f_t12)
    !
    ! low-storage 3rd-order Runge-Kutta scheme 
    ! for time integration of the psi field.
    !
    use mod_param     , only: h01,h02,rho1,rho2,cp1,cp2
    implicit none
    !
    real(rp), intent(in   )                                        :: f_t1,f_t2,f_t12
    integer , intent(in   ), dimension(3)                          :: n
    real(rp), intent(in   ), dimension(3)                          :: dli
    integer , intent(in   )                                        :: nh_d,nh_u,nh_v
    real(rp), intent(in   ), dimension(1-nh_d:)                    :: dzc,dzf
    real(rp), intent(in   )                                        :: u_reg,eps_int
    real(rp), intent(in   ), dimension(1-nh_u:,1-nh_u:,1-nh_u:)    :: u,v,w
    real(rp), intent(inout), dimension(1-nh_v:,1-nh_v:,1-nh_v:)    :: psi,ent  ! volume-fraction at current time level
    real(rp), intent(in   ), dimension(1-nh_v:,1-nh_v:,1-nh_v:)    :: ls,src,tmp   ! level-set
    real(rp), intent(in   ), dimension(1-nh_v:,1-nh_v:,1-nh_v:,1:) :: nor
    real(rp), intent(in ), dimension( 0:, 0:, 0:)               :: kappa
    real(rp), intent(inout), dimension(     1:,     1:,     1:)    :: dpsidtrko
    ! 
    real(rp), dimension(n(1),n(2),n(3)) :: dpsidtrk
    real(rp) :: tot_reg,phase_change,diffent
!#if defined(_CONTACT_LINE)
    real(rp) :: int_den,q_t,tot_lag,d_psi
!#endif
    integer  :: i,j,k
    !
    ! 1. compute regularization
    !
#if defined(_USE_ACDI)    
    call regularization_psi(n,dli,nh_d,nh_v,dzc,dzf,u_reg,eps_int,psi,ls,nor,dpsidtrk)
#else    
    call regularization_psi_cdi(n,dli,nh_d,nh_v,dzc,dzf,u_reg,eps_int,psi,nor,dpsidtrk)
#endif
    !
    ! 2. Compute lag mult
    !
    tot_reg = 0._rp
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          tot_reg = tot_reg + dpsidtrk(i,j,k)/(dli(1)*dli(2)*dli(3))
        enddo
      enddo
    enddo
    call mpi_allreduce(MPI_IN_PLACE,tot_reg,1,MPI_REAL_RP,MPI_SUM,comm_cart,ierr)
    if(myid.eq.0) print*, "Tot. reg.", tot_reg
    !
!#if defined(_CONTACT_LINE)
    int_den = 0._rp
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          int_den = int_den + 6._rp*psi(i,j,k)*(1._rp-psi(i,j,k))/(dli(1)*dli(2)*dli(3))
        enddo
      enddo
    enddo
    call mpi_allreduce(MPI_IN_PLACE,int_den,1,MPI_REAL_RP,MPI_SUM,comm_cart,ierr)
    !q_t = tot_reg/int_den
    q_t = 0.0 !! no lag mult
!#endif
    !
    ! 3. Add advection 
    ! 
    call moms_psi(n,dli,nh_d,nh_u,nh_v,u,v,w,psi,dpsidtrk)
    !
     do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
        dpsidtrk(i,j,k) = dpsidtrk(i,j,k)!-src(i,j,k)
       enddo
      enddo
     enddo
    !
    ! 4. Advance the color function
    !
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(i,j,k) &
    !$OMP SHARED(n,f_t1,f_t2,psi,dpsidtrk,dpsidtrko)
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          !
          phase_change     = 0.d0
          tot_lag          = 6._rp*psi(i,j,k)*(1._rp-psi(i,j,k))*q_t
          psi(i,j,k)       = psi(i,j,k) + f_t1*( dpsidtrk(i,j,k)-tot_lag-phase_change ) + f_t2*dpsidtrko(i,j,k)
          dpsidtrko(i,j,k) = dpsidtrk(i,j,k)-tot_lag
          !
        enddo
      enddo
    enddo
    !$OMP END PARALLEL DO
    !
    return
  end subroutine rk_psi

  subroutine rk_psi_ibm(f_t1,f_t2,n,dli,dl,nh_d,nh_u,nh_v,nh_p, &
    dzc,dzf,u_reg,eps_int,u,v,w,psi,src,ls,nor,kappa,dpsidtrko,halo_u,f_t12)
    !
    use mod_boiling,     only: cmpt_regularization,cbcreg,bcreg,regularization_psi_cdi_ibm, &
                               cmpt_regularization_acdi
    use mod_param  ,     only: no_outflow,is_outflow
    use mod_bound  ,     only: bounduvw
    use mod_common_mpi,  only: MPI_PROC_NULL,left,right,front,back,top,bottom
    !
    ! low-storage 3rd-order Runge-Kutta scheme 
    ! for time integration of the psi field.
    !
    implicit none
    !
    real(rp), intent(in   )                                        :: f_t1,f_t2,f_t12
    integer , intent(in   ), dimension(3)                          :: n
    real(rp), intent(in   ), dimension(3)                          :: dli,dl
    integer , intent(in   )                                        :: nh_d,nh_u,nh_v,nh_p
    real(rp), intent(in   ), dimension(1-nh_d:)                    :: dzc,dzf
    real(rp), intent(in   )                                        :: u_reg,eps_int
    real(rp), intent(in   ), dimension(1-nh_u:,1-nh_u:,1-nh_u:)    :: u,v,w
    real(rp), intent(inout), dimension(1-nh_v:,1-nh_v:,1-nh_v:)    :: psi  ! volume-fraction at current time level
    real(rp), intent(in   ), dimension(1-nh_v:,1-nh_v:,1-nh_v:)    :: ls,src   ! level-set
    real(rp), intent(in   ), dimension(1-nh_v:,1-nh_v:,1-nh_v:,1:) :: nor
    real(rp), intent(in ), dimension( 0:, 0:, 0:)               :: kappa
    real(rp), intent(inout), dimension(     1:,     1:,     1:)    :: dpsidtrko
    integer , intent(in   ), dimension(3)                       :: halo_u
    ! 
    real(rp), dimension(n(1),n(2),n(3)) :: dpsidtrk
    real(rp) :: tot_reg,phase_change
    !#if defined(_CONTACT_LINE)
    real(rp) :: int_den,q_t,tot_lag
    !#endif
    integer  :: i,j,k
    real(rp), dimension(:,:,:), allocatable :: helper
    !
    dpsidtrk(:,:,:) = 0.0
    !
    !
    !
    ! 1. compute regularization
    !      
    allocate(helper(1:n(1),1:n(2),1:n(3)))
    helper(:,:,:) = 0.0
    call regularization_psi_cdi_ibm(n,dli,nh_d,nh_v,nh_u,helper)
    !
    !
    ! 2. Compute lag mult
    !
    tot_reg = 0._rp
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          tot_reg = tot_reg + helper(i,j,k)/(dli(1)*dli(2)*dli(3))
        enddo
      enddo
    enddo
    call mpi_allreduce(MPI_IN_PLACE,tot_reg,1,MPI_REAL_RP,MPI_SUM,comm_cart,ierr)
    if(myid.eq.0) print*, "Tot. reg.", tot_reg
    !
    ! 3. Add advection 
    ! 
    call moms_psi(n,dli,nh_d,nh_u,nh_v,u,v,w,psi,dpsidtrk)
    !
    ! 4. Advance the color function
    !
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(i,j,k) &
    !$OMP SHARED(n,f_t1,f_t2,psi,dpsidtrk,dpsidtrko)
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          !                   advection      +  regularization
          dpsidtrk(i,j,k)  = dpsidtrk(i,j,k) + helper(i,j,k)
          psi(i,j,k)       = psi(i,j,k) +  f_t1*dpsidtrk(i,j,k)  &
                                        +  f_t2*dpsidtrko(i,j,k)
          dpsidtrko(i,j,k) = dpsidtrk(i,j,k)
        enddo
      enddo
    enddo
    !
    !$OMP END PARALLEL DO
    !
    deallocate(helper)
    return
  end subroutine rk_psi_ibm
  !
  subroutine moms_psi(n,dli,nh_d,nh_u,nh_v,u,v,w,s,dsdt)
    !
    implicit none
    !
    integer , intent(in   ), dimension(3)                       :: n
    real(rp), intent(in   ), dimension(3)                       :: dli
    integer , intent(in   )                                     :: nh_d,nh_u,nh_v
    real(rp), intent(in   ), dimension(1-nh_u:,1-nh_u:,1-nh_u:) :: u,v,w
    real(rp), intent(in   ), dimension(1-nh_v:,1-nh_v:,1-nh_v:) :: s
    real(rp), intent(inout), dimension(     1:,     1:,     1:) :: dsdt
    !
    real(rp) :: s_flxp,s_flxm,s_flyp,s_flym,s_flzp,s_flzm
    integer  :: i,j,k
    !
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          !
          ! along x
          !
          s_flxp = 0.5_rp*( s(i+1,j,k)+s(i,j,k) )*u(i  ,j,k)
          s_flxm = 0.5_rp*( s(i-1,j,k)+s(i,j,k) )*u(i-1,j,k)
          !
          ! along y
          !
          s_flyp = 0.5_rp*( s(i,j+1,k)+s(i,j,k) )*v(i,j  ,k)
          s_flym = 0.5_rp*( s(i,j-1,k)+s(i,j,k) )*v(i,j-1,k)
          !
          ! along z
          !
          s_flzp = 0.5_rp*( s(i,j,k  )+s(i,j,k+1) )*w(i,j,k  )
          s_flzm = 0.5_rp*( s(i,j,k  )+s(i,j,k-1) )*w(i,j,k-1)
          !
          ! overall balance
          !
          dsdt(i,j,k) = - ( dli(1)*(s_flxp-s_flxm) + &
                            dli(2)*(s_flyp-s_flym) + &
                            dli(3)*(s_flzp-s_flzm) ) + dsdt(i,j,k)
          !
        enddo
      enddo
    enddo
    !
    return
  end subroutine moms_psi
  !
  subroutine moms_reg(n,dli,nh_d,nh_u,nh_v,u,v,w,s,dsdt)
    !
    use mod_boiling     ,only: rex,rey,rez

    !
    implicit none
    !
    integer , intent(in   ), dimension(3)                       :: n
    real(rp), intent(in   ), dimension(3)                       :: dli
    integer , intent(in   )                                     :: nh_d,nh_u,nh_v
    real(rp), intent(in   ), dimension(1-nh_u:,1-nh_u:,1-nh_u:) :: u,v,w
    real(rp), intent(in   ), dimension(1-nh_v:,1-nh_v:,1-nh_v:) :: s
    real(rp), intent(inout), dimension(     1:,     1:,     1:) :: dsdt
    !
    real(rp) :: s_flxp,s_flxm,s_flyp,s_flym,s_flzp,s_flzm
    integer  :: i,j,k
    !
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          !
          ! along x
          !
          s_flxp = 0.5_rp*( s(i+1,j,k)+s(i,j,k) )*u(i  ,j,k) -  rex(i,j,k)
          s_flxm = 0.5_rp*( s(i-1,j,k)+s(i,j,k) )*u(i-1,j,k) -  rex(i-1,j,k)
          !
          ! along y
          !
          s_flyp = 0.5_rp*( s(i,j+1,k)+s(i,j,k) )*v(i,j  ,k) -  rey(i,j,k)
          s_flym = 0.5_rp*( s(i,j-1,k)+s(i,j,k) )*v(i,j-1,k) -  rey(i,j-1,k)
          !
          ! along z
          !
          s_flzp = 0.5_rp*( s(i,j,k+1)+s(i,j,k) )*w(i,j,k  ) -  rez(i,j,k)
          s_flzm = 0.5_rp*( s(i,j,k-1)+s(i,j,k) )*w(i,j,k-1) -  rez(i,j,k-1)
          !
          ! overall balance
          !
          !
          dsdt(i,j,k) = - ( dli(1)*(s_flxp-s_flxm) + &
                            dli(2)*(s_flyp-s_flym) + &
                            dli(3)*(s_flzp-s_flzm) ) + dsdt(i,j,k)
        enddo
      enddo
    enddo
    !
    return
  end subroutine moms_reg

  !
   subroutine moms_c1(n,dli,nh_d,nh_u,nh_v,u,v,w,ur,vr,wr,s,phi,nor,dsdt)
    !
    implicit none
    !
    integer , intent(in ), dimension(3)                       :: n
    real(rp), intent(in ), dimension(3)                       :: dli
    integer , intent(in )                                     :: nh_d,nh_u,nh_v
    real(rp), intent(in ), dimension(1-nh_u:,1-nh_u:,1-nh_u:) :: u,v,w
    real(rp), intent(in ), dimension(1-nh_v:,1-nh_v:,1-nh_v:) :: s      !concentration
    real(rp), intent(in ), dimension(1-nh_v:,1-nh_v:,1-nh_v:)    :: phi !phase field 
    real(rp), intent(in ), dimension(1-nh_v:,1-nh_v:,1-nh_v:,1:) :: nor
    real(rp), intent(in ), dimension(0:,0:,0:)                   :: ur,vr,wr
    real(rp), intent(inout), dimension(      1:,     1:,     1:)    :: dsdt
    !
    integer  :: i,j,k,im,jm,km,ip,jp,kp,imm,ipp,jmm,jpp,kmm,kpp
    real(rp) :: s_flxp,s_flxm,s_flyp,s_flym,s_flzp,s_flzm      ,&
                Urx_p,Urx_m,Ury_p,Ury_m,Urz_p,Urz_m,A, &
                phi_face_xp,phi_face_xm,phi_face_yp,phi_face_ym, &
                phi_face_zp,phi_face_zm, nor_face_xp,nor_face_xm, &
                nor_face_yp,nor_face_ym,nor_face_zp,nor_face_zm, &
                c_face_xp,c_face_xm,c_face_yp,c_face_ym,c_face_zp, &
                c_face_zm,ur_face_xp,ur_face_xm,vr_face_yp,vr_face_ym, &
                wr_face_zp, wr_face_zm
    !    
    real(rp) :: abs_u_ijk,ureg_x,ureg_y,ureg_z,urmax
    !
    !
    do k=1,n(3)
      kp = k+1
      km = k-1
      do j=1,n(2)
        jp = j+1
        jm = j-1
        do i=1,n(1)
          ip  = i+1
          im  = i-1

          ! along x
          !
          ! interpolate relative velocity to face
          !
          ur_face_xp  = (ur(ip,j,k)+ur(i,j,k))*0.5
          ur_face_xm  = (ur(im,j,k)+ur(i,j,k))*0.5
          vr_face_yp  = (vr(i,jp,k)+vr(i,j,k))*0.5
          vr_face_ym  = (vr(i,jm,k)+vr(i,j,k))*0.5
          wr_face_zp  = (wr(i,j,kp)+wr(i,j,k))*0.5
          wr_face_zm  = (wr(i,j,km)+wr(i,j,k))*0.5
          !
          ! interpolate phi to face
          !
          phi_face_xp = (phi(ip,j,k)+phi(i,j,k))*0.5
          phi_face_xm = (phi(im,j,k)+phi(i,j,k))*0.5
          !
          ! interpolate normal to face
          !
          nor_face_xp = (nor(ip,j,k,1)+nor(i,j,k,1))*0.5
          nor_face_xm = (nor(im,j,k,1)+nor(i,j,k,1))*0.5
          ! 
          ! interpolate scalar to face 
          !
          c_face_xp   = (s(ip,j,k)+s(i,j,k))*0.5
          c_face_xm   = (s(im,j,k)+s(i,j,k))*0.5
          !
          ! compute u*c term at the face
          !
          s_flxp      = u(i,j,k)*c_face_xp
          s_flxm      = u(im,j,k)*c_face_xm
          !
          Urx_p       = ur_face_xp*c_face_xp
          Urx_m       = ur_face_xm*c_face_xm
          !
          !
          !
          ! along y
          !
          ! interpolate phi to face
          !
          phi_face_yp = (phi(i,jp,k)+phi(i,j,k))*0.5
          phi_face_ym = (phi(i,jm,k)+phi(i,j,k))*0.5
          !
          ! interpolate normal to face
          !
          nor_face_yp = (nor(i,jp,k,2)+nor(i,j,k,2))*0.5
          nor_face_ym = (nor(i,jm,k,2)+nor(i,j,k,2))*0.5
          ! 
          ! interpolate scalar to face 
          !
          c_face_yp   = (s(i,jp,k)+s(i,j,k))*0.5
          c_face_ym   = (s(i,jm,k)+s(i,j,k))*0.5
          !
          ! compute u*c term at the face
          !
          s_flyp      = v(i,j,k)*c_face_yp
          s_flym      = v(i,jm,k)*c_face_ym
          !
          !
          Ury_p       = vr_face_yp*c_face_yp
          Ury_m       = vr_face_ym*c_face_ym
          !        
          !
          ! along z
          !
          ! interpolate phi to face
          !
          phi_face_zp = (phi(i,j,kp)+phi(i,j,k))*0.5
          phi_face_zm = (phi(i,j,km)+phi(i,j,k))*0.5
          !
          ! interpolate normal to face
          !
          nor_face_zp = (nor(i,j,kp,3)+nor(i,j,k,3))*0.5
          nor_face_zm = (nor(i,j,km,3)+nor(i,j,k,3))*0.5
          ! 
          ! interpolate scalar to face 
          !
          c_face_zp   = (s(i,j,kp)+s(i,j,k))*0.5
          c_face_zm   = (s(i,j,km)+s(i,j,k))*0.5
          !
          ! compute u*c term at the face
          !
          s_flzp      = w(i,j,k)*c_face_zp
          s_flzm      = w(i,j,km)*c_face_zm
          !
          !        
          Urz_p        = wr_face_zp*c_face_zp
          Urz_m        = wr_face_zm*c_face_zm
          !
          ! compute the divergence between i+1/2 and i-1/2
          !
          ! overall balance
          !
#if defined(_TWOD)          
          s_flxp = 0.d0
          s_flxm = 0.d0
          Urx_p = 0.d0
          Urx_m = 0.d0
#endif
          dsdt(i,j,k) = - ( dli(1)*((s_flxp+Urx_p)-(s_flxm+Urx_m)) + &
                            dli(2)*((s_flyp+Ury_p)-(s_flym+Ury_m)) + &
                            dli(3)*((s_flzp+Urz_p)-(s_flzm+Urz_m)) )  + dsdt(i,j,k) 
        enddo
      enddo
    enddo

    return
  end subroutine moms_c1
  !
  !

   subroutine cmpt_rel_vel(n,dli,nh_d,nh_u,nh_v,u,v,w,s,phi,nor,ur,vr,wr)
    !
    !
    use mod_param, only: A_rel 
    !
    implicit none
    !
    integer , intent(in ), dimension(3)                       :: n
    real(rp), intent(in ), dimension(3)                       :: dli
    integer , intent(in )                                     :: nh_d,nh_u,nh_v
    real(rp), intent(in ), dimension(1-nh_u:,1-nh_u:,1-nh_u:) :: u,v,w
    real(rp), intent(in ), dimension(1-nh_v:,1-nh_v:,1-nh_v:) :: s      !concentration
    real(rp), intent(in ), dimension(1-nh_v:,1-nh_v:,1-nh_v:)    :: phi !phase field 
    real(rp), intent(in ), dimension(1-nh_v:,1-nh_v:,1-nh_v:,1:) :: nor
    real(rp), intent(out), dimension(      0:,      0:,      0:) :: ur,vr,wr
    !
    integer  :: i,j,k,im,jm,km,ip,jp,kp,imm,ipp,jmm,jpp,kmm,kpp
    real(rp) :: s_flxp,s_flxm,s_flyp,s_flym,s_flzp,s_flzm      ,&
                Urx_p,Urx_m,Ury_p,Ury_m,Urz_p,Urz_m,c_sat, &
                phi_face_xp,phi_face_xm,phi_face_yp,phi_face_ym, &
                phi_face_zp,phi_face_zm, nor_face_xp,nor_face_xm, &
                nor_face_yp,nor_face_ym,nor_face_zp,nor_face_zm, &
                c_face_xp,c_face_xm,c_face_yp,c_face_ym,c_face_zp, &
                c_face_zm
    !    
    real(rp) :: abs_u_ijk,ureg_x,ureg_y,ureg_z,eps
    !
    eps   = 1e-30
    !
    c_sat = 1.0d0
    !
    !
    !Attraction Coefficient  
    !
    !compute rel_vel*phi at the cell centers
  do k=1,n(3)
   do j=1,n(2)
     do i=1,n(1)
          ur(i,j,k) =  A_rel*(0.5-phi(i,j,k))*abs(phi(i,j,k)-1.d0)/(abs(phi(i,j,k)-1.d0)+eps)* &
                          (c_sat - s(i,j,k))*s(i,j,k)**2*nor(i,j,k,1)*phi(i,j,k)
          !
          vr(i,j,k) =  A_rel*(0.5-phi(i,j,k))*abs(phi(i,j,k)-1.d0)/(abs(phi(i,j,k)-1.d0)+eps)* &
                          (c_sat - s(i,j,k))*s(i,j,k)**2*nor(i,j,k,2)*phi(i,j,k)
          !        
          wr(i,j,k) =  A_rel*(0.5-phi(i,j,k))*abs(phi(i,j,k)-1.d0)/(abs(phi(i,j,k)-1.d0)+eps)* &
                          (c_sat - s(i,j,k))*s(i,j,k)**2*nor(i,j,k,3)*phi(i,j,k)
         enddo
      enddo
    enddo

    return
  end subroutine cmpt_rel_vel
  !
  !
  subroutine cmpt_umax(n,nh_u,u,v,w,gamma_v,u_reg)
    !
    use  my_is_nan_mod , only: my_is_nan_real64
    !
    implicit none
    !
    integer , intent(in ), dimension(3)                       :: n
    integer , intent(in )                                     :: nh_u
    real(rp), intent(in ), dimension(1-nh_u:,1-nh_u:,1-nh_u:) :: u,v,w
    real(rp), intent(in )                                     :: gamma_v
    real(rp), intent(out)                                     :: u_reg
    !
    real(rp) :: abs_u_ijk ,ux,vx,wx,uy,vy,wy,uz,vz,wz
    integer  :: i,j,k
    !
    u_reg = 0._rp 
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          !
          ux = abs(u(i,j,k))
          vx = 0.25_rp*abs( v(i,j,k)+v(i,j-1,k)+v(i+1,j,k)+v(i+1,j-1,k) )
          wx = 0.25_rp*abs( w(i,j,k)+w(i,j,k-1)+w(i+1,j,k)+w(i+1,j,k-1) )
          abs_u_ijk = sqrt(ux**2+vx**2+wx**2)
          u_reg     = max(u_reg,abs_u_ijk)
          !
          uy = 0.25_rp*abs( u(i,j,k)+u(i,j+1,k)+u(i-1,j+1,k)+u(i-1,j,k) )
          vy = abs(v(i,j,k))
          wy = 0.25_rp*abs( w(i,j,k)+w(i,j+1,k)+w(i,j+1,k-1)+w(i,j,k-1) )
          abs_u_ijk = sqrt(uy**2+vy**2+wy**2)
          u_reg     = max(u_reg,abs_u_ijk)
          !
          uz = 0.25_rp*abs( u(i,j,k)+u(i-1,j,k)+u(i-1,j,k+1)+u(i,j,k+1) )
          vz = 0.25_rp*abs( v(i,j,k)+v(i,j-1,k)+v(i,j-1,k+1)+v(i,j,k+1) )
          wz = abs(w(i,j,k))
          abs_u_ijk = sqrt(uz**2+vz**2+wz**2)
          u_reg     = max(u_reg,abs_u_ijk)
          !
        enddo
      enddo
    enddo
    !u_reg = 0.1
    call mpi_allreduce(MPI_IN_PLACE,u_reg,1,MPI_REAL_RP,MPI_MAX,comm_cart,ierr)
    if (u_reg.EQ.0.0_rp ) u_reg = 0.01_rp
    u_reg = gamma_v*u_reg
    !
    return
  end subroutine cmpt_umax
    !  
!    subroutine regularization_c1(n,dli,nh_d,nh_v, &
!                                 dzc,dzf,diff,u_reg,eps,s,phi,nor,dsdt) !s corresponds to c as concent.
!                                                                        !phi corresponds to volume_fraction/phi
!    !
!    implicit none
!    !
!    integer , intent(in ), dimension(3)               :: n
!    real(rp), intent(in ), dimension(3)               :: dli 
!    integer , intent(in )                             :: nh_d,nh_v
!    real(rp), intent(in ), dimension(1-nh_d:)         :: dzc,dzf
!    real(rp), intent(in )                             :: u_reg,eps
!    real(rp), intent(in ), dimension(1-nh_v:,1-nh_v:,1-nh_v:)    :: s
!    real(rp), intent(in ), dimension(1-nh_v:,1-nh_v:,1-nh_v:)    :: Phi,diff
!    real(rp), intent(in ), dimension(1-nh_v:,1-nh_v:,1-nh_v:,1:) :: nor
!    real(rp), intent(out), dimension(1:,1:,1:)                   :: dsdt
!    real(rp),              dimension(0:n(1)+1,0:n(2)+1,0:n(3)+1) :: gradsy,gradsz
!    !
!    real(rp), dimension(3) :: dl
!    real(rp) :: s_face,term_p,term_m,nor_face,s_min,s_max,rho_p,rho_m 
!    real(rp) :: Phi_face,A,D,urmax
!    integer  :: i,j,k
!    !
!    dl(:) = 1.d0/dli(:)
!    dsdt(:,:,:) = 0.d0
!    !
!    !Attraction coefficient
!    !A = 2.d0
!    !urmax = 0.074*A
!    !
!    ! tune the diffusion coefficient to satisfy the stability/boundness criteria 
!    !
!    !D = 2*dl(2)*(u_reg+urmax)
!    !
!    !D = diff
!    !
!    do k=1,n(3)
!      do j=1,n(2)
!        do i=1,n(1)
!         !
!
!         ! !
!         ! ! along x-dir
!         ! !
!         ! ! Salar/Armin Model
!         ! !
!         ! dsdt(i,j,k) = diff(i,j,k)*(s(i+1,j,k)-2.d0*s(i,j,k)+s(i-1,j,k))*dli(1)**2 - &
!         !               diff(i,j,k)/eps*((1-phi(i+1,j,k))*nor(i+1,j,k,1)*s(i+1,j,k) - &
!         !               (1-phi(i-1,j,k))*nor(i-1,j,k,1)*s(i-1,j,k))*0.5d0*dli(1)*0.d0
!         ! !
!         ! ! along y-dir
!         ! !
!         ! dsdt(i,j,k) = dsdt(i,j,k)+ &
!         !               diff(i,j,k)*(s(i,j+1,k)-2.d0*s(i,j,k)+s(i,j-1,k))*dli(2)**2 - &
!         !               D/eps*((1-phi(i,j+1,k))*nor(i,j+1,k,2)*s(i,j+1,k)  - &
!         !               (1-phi(i,j-1,k))*nor(i,j-1,k,2)*s(i,j-1,k))*0.5d0*dli(2)*0.d0
!         ! !
!         ! !      
!         ! ! along z-dir
!         ! !
!         ! dsdt(i,j,k)  = dsdt(i,j,k)+ &
!         !                diff(i,j,k)*(s(i,j,k+1)-2.d0*s(i,j,k)+s(i,j,k-1))*dli(3)**2 -&
!         !                diff(i,j,k)/eps*((1-phi(i,j,k+1))*nor(i,j,k+1,3)*s(i,j,k+1) - &
!         !                (1-phi(i,j,k-1))*nor(i,j,k-1,3)*s(i,j,k-1))*0.5d0*dli(3)*0.d0
!         ! !
!         ! !
!        enddo
!      enddo
!    enddo
!    !
!  return 
!  end subroutine regularization_c1

    subroutine regularization_c1(n,dli,nh_d,nh_v, &
                                 dzc,dzf,diff,gm,eps,s,phi,nor,dsdt) !s corresponds to c as concent.
                                                                        !phi corresponds to volume_fraction/phi
    !
    implicit none
    !
    integer , intent(in ), dimension(3)               :: n
    real(rp), intent(in ), dimension(3)               :: dli 
    integer , intent(in )                             :: nh_d,nh_v
    real(rp), intent(in ), dimension(1-nh_d:)         :: dzc,dzf
    real(rp), intent(in )                             :: diff,gm,eps
    real(rp), intent(in ), dimension(1-nh_v:,1-nh_v:,1-nh_v:)    :: s
    real(rp), intent(in ), dimension(1-nh_v:,1-nh_v:,1-nh_v:)    :: Phi
    real(rp), intent(in ), dimension(1-nh_v:,1-nh_v:,1-nh_v:,1:) :: nor
    real(rp), intent(out), dimension(1:,1:,1:)                   :: dsdt
    real(rp),              dimension(0:n(1)+1,0:n(2)+1,0:n(3)+1) :: r
    !
    real(rp), dimension(3) :: dl
    real(rp) :: s_face,term_p,term_m,nor_face,s_min,s_max,rho_p,rho_m 
    real(rp) :: Phi_face,A,D,urmax
    real(rp) :: gradrxp,gradrxm,gradryp,gradrym,gradrzp,gradrzm
    real(rp) :: ratioxp,ratioxm,ratioyp,ratioym,ratiozp,ratiozm
    real(rp) :: gradphixp,gradphixm,gradphiyp,gradphiym,gradphizp,gradphizm
    real(rp) :: phixp,phixm,phiyp,phiym,phizp,phizm
    real(rp) :: regxp,regxm,regyp,regym,regzp,regzm
    real(rp) :: diffxp,diffxm,diffyp,diffym,diffzp,diffzm
    real(rp) :: rhsxp,rhsxm,rhsyp,rhsym,rhszp,rhszm
    real(rp) :: small
    integer  :: i,j,k
    integer  :: ip,im,jp,jm,kp,km
    !
    dl(:) = 1.d0/dli(:)
    dsdt(:,:,:) = 0.d0
    !
    !
    D = diff
    !
    small = 1e-13
    !
 do k=0,n(3)+1
  do j=0,n(2)+1
   do i=0,n(1)+1
   r(i,j,k) = s(i,j,k)/(phi(i,j,k)+small)
   enddo
  enddo
 enddo

    !
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          !
          ip = i+1
          im = i-1
          jp = j+1
          jm = j-1
          kp = k+1
          km = k-1
          !
          ! along x-dir
          !
          ! Salar/Armin Model
          !
          !dsdt(i,j,k) = D*(s(i+1,j,k)-2.d0*s(i,j,k)+s(i-1,j,k))*dli(1)**2 - &
          !              D/eps*((1-phi(i+1,j,k))*nor(i+1,j,k,1)*s(i+1,j,k) - &
          !              (1-phi(i-1,j,k))*nor(i-1,j,k,1)*s(i-1,j,k))*0.5d0*dli(1)

#if defined(_TWOD)
          rhsxp = 0.d0
          rhsxm = 0.d0
#else          
          gradrxp = ( r(ip,j,k) - r(i,j,k)   )*dli(1) !gradient of r=c/phi is computed at face i+1/2
          gradrxm = ( r(i,j,k)  - r(im,j,k)  )*dli(1) !gradient of r=c/phi is computed at face i-1/2

          ratioxp = ( r(i,j,k)   + r(ip,j,k) )*0.5    !phi/c is computed at face i+1/2
          ratioxm = ( r(i,j,k)   + r(im,j,k) )*0.5    !phi/c is computed at face i-1/2

          
          phixp   = ( phi(ip,j,k) + phi(i,j,k) )*0.5  !phi computed at i+1/2
          phixm   = ( phi(im,j,k) + phi(i,j,k) )*0.5  !phi computed at i-1/2

          gradphixp = ( phi(ip,j,k) - phi(i,j,k)   )*dli(1) !gradient phi computed at i+1/2
          gradphixm = ( phi(i,j,k)  - phi(im,j,k)  )*dli(1) !gradient phi computed at i-1/2

          regxp     = gm*(eps*gradphixp - ((phi(ip,j,k)*(1.d0-phi(ip,j,k))*nor(ip,j,k,1))  + &
                                           (phi(i,j,k)*(1.d0-phi(i,j,k))*nor(i,j,k,1)))*0.5)*ratioxp         !reg*ratio term computed at i+1/2

          regxm     = gm*(eps*gradphixm - ((phi(im,j,k)*(1.d0-phi(im,j,k))*nor(im,j,k,1))  + &
                                           (phi(i,j,k)*(1.d0-phi(i,j,k))*nor(i,j,k,1)))*0.5)*ratioxm         !reg*ratio term computed at i-1/2
          
          diffxp    = D*phixp*gradrxp                                                                        !diff term computed at i+1/2
          diffxm    = D*phixm*gradrxm                                                                        !diff term computed at i-1/2
          ! 
          rhsxp     = regxp + diffxp
          rhsxm     = regxm + diffxm                                                                                                            
#endif
          !       
          !
          !
          ! along y-dir
          !
          !dsdt(i,j,k) = dsdt(i,j,k)+ &
          !              D*(s(i,j+1,k)-2.d0*s(i,j,k)+s(i,j-1,k))*dli(2)**2 - &
          !              D/eps*((1-phi(i,j+1,k))*nor(i,j+1,k,2)*s(i,j+1,k)  - &
          !              (1-phi(i,j-1,k))*nor(i,j-1,k,2)*s(i,j-1,k))*0.5d0*dli(2)
 
           gradryp = ( r(i,jp,k) - r(i,j,k)   )*dli(2) !gradient of r=c/phi is computed at face j+1/2
           gradrym = ( r(i,j,k)  - r(i,jm,k)  )*dli(2) !gradient of r=c/phi is computed at face j-1/2
                                                                                                                                                 
           ratioyp = ( r(i,j,k)   + r(i,jp,k) )*0.5    !phi/c is computed at face j+1/2
           ratioym = ( r(i,j,k)   + r(i,jm,k) )*0.5    !phi/c is computed at face j-1/2
                                                                                                                                                 
           phiyp   = ( phi(i,jp,k) + phi(i,j,k) )*0.5  !phi computed at j+1/2
           phiym   = ( phi(i,jm,k) + phi(i,j,k) )*0.5  !phi computed at j-1/2
                                                                                                                                                 
           gradphiyp = ( phi(i,jp,k) - phi(i,j,k)   )*dli(2)  !gradient phi computed at j+1/2
           gradphiym = ( phi(i,j,k)  - phi(i,jm,k)  )*dli(2)  !gradient phi computed at j-1/2
                                                                                                                                                 
           regyp     = gm*(eps*gradphiyp - ((phi(i,jp,k)*(1.d0-phi(i,jp,k))*nor(i,jp,k,2))  + &
                                            (phi(i,j,k)*(1.d0-phi(i,j,k))*nor(i,j,k,2)))*0.5)*ratioyp         !reg*ratio term computed at j+1/2
                                                                                                                                                 
           regym     = gm*(eps*gradphiym - ((phi(i,jm,k)*(1.d0-phi(i,jm,k))*nor(i,jm,k,2))  + &
                                            (phi(i,j,k)*(1.d0-phi(i,j,k))*nor(i,j,k,2)))*0.5)*ratioym         !reg*ratio term computed at j-1/2
           
           diffyp    = D*phiyp*gradryp                                                                        !diff term computed at j+1/2
           diffym    = D*phiym*gradrym                                                                        !diff term computed at j-1/2
            
           rhsyp     = regyp + diffyp
           rhsym     = regym + diffym                                                                                                                                                                                                                      
          
                
          ! along z-dir
          !
          !dsdt(i,j,k)  = dsdt(i,j,k)+ &
          !               D*(s(i,j,k+1)-2.d0*s(i,j,k)+s(i,j,k-1))*dli(3)**2 -&
          !               D/eps*((1-phi(i,j,k+1))*nor(i,j,k+1,3)*s(i,j,k+1) - &
          !               (1-phi(i,j,k-1))*nor(i,j,k-1,3)*s(i,j,k-1))*0.5d0*dli(3)
          !
          gradrzp = ( r(i,j,kp) - r(i,j,k)   )*dli(3) !gradient of r=c/phi is computed at face k+1/2
          gradrzm = ( r(i,j,k)  - r(i,j,km)  )*dli(3) !gradient of r=c/phi is computed at face k-1/2
                                                                                                                                                
          ratiozp = ( r(i,j,k)   + r(i,j,kp) )*0.5    !phi/c is computed at face k+1/2
          ratiozm = ( r(i,j,k)   + r(i,j,km) )*0.5    !phi/c is computed at face k-1/2
                                                                                                                                                
          phizp   = ( phi(i,j,kp) + phi(i,j,k) )*0.5  !phi computed at k+1/2
          phizm   = ( phi(i,j,km) + phi(i,j,k) )*0.5  !phi computed at k-1/2
                                                                                                                                                
          gradphizp = ( phi(i,j,kp) - phi(i,j,k)   )*dli(3)  !gradient phi computed at k+1/2
          gradphizm = ( phi(i,j,k)  - phi(i,j,km)  )*dli(3)  !gradient phi computed at k-1/2
                                                                                                                                                
          regzp     = gm*(eps*gradphizp - ((phi(i,j,kp)*(1.d0-phi(i,j,kp))*nor(i,j,kp,3))  + &
                                           (phi(i,j,k)*(1.d0-phi(i,j,k))*nor(i,j,k,3)))*0.5)*ratiozp         !reg*ratio term computed at k+1/2
                                                                                                                                                
          regzm     = gm*(eps*gradphizm - ((phi(i,j,km)*(1.d0-phi(i,j,km))*nor(i,j,km,3))  + &
                                           (phi(i,j,k)*(1.d0-phi(i,j,k))*nor(i,j,k,3)))*0.5)*ratiozm         !reg*ratio term computed at k-1/2
          
          diffzp    = D*phizp*gradrzp                                                           !diff term computed at k+1/2
          diffzm    = D*phizm*gradrzm                                                           !diff term computed at k-1/2
           
          rhszp     = regzp + diffzp
          rhszm     = regzm + diffzm                                                                                                              
          !
          dsdt(i,j,k)    =    (rhsxp - rhsxm) * dli(1) + &
                              (rhsyp - rhsym) * dli(2) + &
                              (rhszp - rhszm) * dli(3)    
          !
        enddo
      enddo
    enddo
    !
    return
  end subroutine regularization_c1
  !
    subroutine interfacial_cflux(n,dli,nh_d,nh_v, &
                                 dzc,dzf,diff1,diff2,c1,c2,phi,dsdt,id) 
    !
    use mod_param, only: A_m,k_eq 
    !
    implicit none
    !
    integer , intent(in ), dimension(3)               :: n
    real(rp), intent(in ), dimension(3)               :: dli 
    integer , intent(in )                             :: nh_d,nh_v
    real(rp), intent(in ), dimension(1-nh_d:)         :: dzc,dzf
    real(rp), intent(in )                             :: diff1,diff2
    real(rp), intent(in ), dimension(1-nh_v:,1-nh_v:,1-nh_v:)    :: c1,c2
    real(rp), intent(in ), dimension(1-nh_v:,1-nh_v:,1-nh_v:)    :: Phi
    real(rp), intent(inout), dimension(1:,1:,1:)                 :: dsdt
    integer , intent(in)                                         :: id
    !
    real(rp), dimension(3) :: dl
    real(rp) :: term1,term2_x,term2_y,term2_z,D_m
    integer  :: i,j,k
    !
    dl(:) = 1.d0/dli(:)
    term1 = 0.d0
    term2_x = 0.d0
    term2_y = 0.d0
    term2_z = 0.d0
    !
    !
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          !
          ! along x-dir
          !
          D_m = diff1*diff2/(k_eq*diff1*(1.d0-phi(i,j,k))+diff2*phi(i,j,k))
          !
          !First-term 
          term1 = A_m*D_m*(k_eq*c2(i,j,k)*phi(i,j,k)-c1(i,j,k)*(1.d0-phi(i,j,k)))
          !
          !Second-term along-x
          
          term2_x = -D_m*(phi(i+1,j,k)-phi(i-1,j,k))*0.5*dli(1)*(c1(i+1,j,k)-c1(i-1,j,k))*0.5*dli(1) - &
                     D_m*k_eq*(phi(i+1,j,k)-phi(i-1,j,k))*0.5*dli(1)*(c2(i+1,j,k)-c2(i-1,j,k))*0.5*dli(1) 
          !
          !Second-term along-y
          !
          term2_y = -D_m*(phi(i,j+1,k)-phi(i,j-1,k))*0.5*dli(2)*(c1(i,j+1,k)-c1(i,j-1,k))*0.5*dli(2) - &
                     D_m*k_eq*(phi(i,j+1,k)-phi(i,j-1,k))*0.5*dli(2)*(c2(i,j+1,k)-c2(i,j-1,k))*0.5*dli(2)
          !
          !      
          !Second-term along-z
          !
          term2_z = -D_m*(phi(i,j,k+1)-phi(i,j,k-1))*0.5*dli(3)*(c1(i,j,k+1)-c1(i,j,k-1))*0.5*dli(3) - &
                     D_m*k_eq*(phi(i,j,k+1)-phi(i,j,k-1))*0.5*dli(3)*(c2(i,j,k+1)-c2(i,j,k-1))*0.5*dli(3)
          !
          dsdt(i,j,k) = dsdt(i,j,k)+(term1+term2_x+term2_y+term2_z)*id
          !
        enddo
      enddo
    enddo
    !
    return
  end subroutine interfacial_cflux
  !
    subroutine regularization_psi_cdi(n,dli,nh_d,nh_v, &
                                     dzc,dzf,u_reg,eps_int,s,nor,dsdt)  !s is volume fraction/phi
    !
    implicit none
    !
    integer , intent(in ), dimension(3)               :: n
    real(rp), intent(in ), dimension(3)               :: dli 
    integer , intent(in )                             :: nh_d,nh_v
    real(rp), intent(in ), dimension(1-nh_d:)         :: dzc,dzf
    real(rp), intent(in )                             :: u_reg,eps_int
    real(rp), intent(in ), dimension(1-nh_v:,1-nh_v:,1-nh_v:)    :: s
    real(rp), intent(in ), dimension(1-nh_v:,1-nh_v:,1-nh_v:,1:) :: nor
    real(rp), intent(out), dimension(1:,1:,1:)                   :: dsdt
    !
    real(rp), dimension(3) :: dl
    real(rp) :: s_face,term_p,term_m,nor_face,s_min,s_max,rho_p,rho_m 
    real(rp) :: gradpxp,gradpxm,gradpyp,gradpym,gradpzp,gradpzm
    real(rp) :: psifxp,psifxm,psifyp,psifym,psifzp,psifzm
    real(rp) :: regxp,regxm,regyp,regym,regzp,regzm
    real(rp) :: abs_u_ijk
    real(rp) :: Psi_face
    integer  :: i,j,k,ip,im,jp,jm,kp,km
    !
    dsdt(:,:,:) = 0.d0
    !
    !
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
    !
        ip = i + 1
        im = i - 1
        jp = j + 1
        jm = j - 1
        kp = k + 1
        km = k - 1
#if defined(_TWOD)
      regxp   = 0.d0
      regxm   = 0.d0
#else
      gradpxp        =   ( s(ip,j,k) - s(i,j,k)  ) * dli(1)       !grad(phi) at i+0.5
      gradpxm        =   ( s(i,j,k)  - s(im,j,k) ) * dli(1)       !grad(phi) at i-0.5

      psifxp         =   ( s(i,j,k)*(1.d0-s(i,j,k))*nor(i,j,k,1) + s(ip,j,k)*(1.d0-s(ip,j,k))*nor(ip,j,k,1))*0.5 !phi*(1-phi)*nor at i+0.5 as p.15 Shahab.M
      psifxm         =   ( s(i,j,k)*(1.d0-s(i,j,k))*nor(i,j,k,1) + s(im,j,k)*(1.d0-s(im,j,k))*nor(im,j,k,1))*0.5 !phi*(1-phi)*nor at i-0.5 as p.15 Shahab.M

      regxp          =   u_reg*( eps_int*gradpxp - psifxp)        !regularization at i+0.5
      regxm          =   u_reg*( eps_int*gradpxm - psifxm)        !regularization at i-0.5
#endif

      gradpyp        =   ( s(i,jp,k) - s(i,j,k)  ) * dli(2)       !grad(phi) at j+0.5
      gradpym        =   ( s(i,j,k)  - s(i,jm,k) ) * dli(2)       !grad(phi) at j-0.5

      psifyp         =   ( s(i,j,k)*(1.d0-s(i,j,k))*nor(i,j,k,2) + s(i,jp,k)*(1.d0-s(i,jp,k))*nor(i,jp,k,2))*0.5 !phi*(1-phi)*nor at j+0.5 as p.15 Shahab.M
      psifym         =   ( s(i,j,k)*(1.d0-s(i,j,k))*nor(i,j,k,2) + s(i,jm,k)*(1.d0-s(i,jm,k))*nor(i,jm,k,2))*0.5 !phi*(1-phi)*nor at j-0.5 as p.15 Shahab.M

      regyp          =   u_reg*( eps_int*gradpyp - psifyp)        !regularization at j+0.5
      regym          =   u_reg*( eps_int*gradpym - psifym)        !regularization at j-0.5


      gradpzp        =   ( s(i,j,kp) - s(i,j,k)  ) * dli(3)       !grad(phi) at k+0.5
      gradpzm        =   ( s(i,j,k)  - s(i,j,km) ) * dli(3)       !grad(phi) at k-0.5

      psifzp         =   ( s(i,j,k)*(1.d0-s(i,j,k))*nor(i,j,k,3) + s(i,j,kp)*(1.d0-s(i,j,kp))*nor(i,j,kp,3))*0.5 !phi*(1-phi)*nor at k+0.5 as p.15 Shahab.M
      psifzm         =   ( s(i,j,k)*(1.d0-s(i,j,k))*nor(i,j,k,3) + s(i,j,km)*(1.d0-s(i,j,km))*nor(i,j,km,3))*0.5 !phi*(1-phi)*nor at k-0.5 as p.15 Shahab.M

     ! if (k.EQ.n(3)) then
       ! psifzp       =   ( s(i,j,k)*(1.d0-s(i,j,k))*nor(i,j,k,3))*1.0
       ! gradpzp      =   0.0
      !endif

      regzp          =   u_reg*( eps_int*gradpzp - psifzp)        !regularization at k+0.5
      regzm          =   u_reg*( eps_int*gradpzm - psifzm)        !regularization at k-0.5

      if (km.EQ.0) then
        regzm        =   0.0
      endif

      if (k.EQ.n(3)) then
        regzp        =   0.0
      endif

      dsdt(i,j,k)    =    (regxp - regxm) * dli(1) + &
                          (regyp - regym) * dli(2) + &
                          (regzp - regzm) * dli(3)




!     dsdt(i,j,k) =    ! U_reg*eps*(s(i+1,j,k)-2.d0*s(i,j,k)+s(i-1,j,k))*dli(1)**2 - &
!                      ! U_reg*    ((1-s(i+1,j,k))*nor(i+1,j,k,1)*s(i+1,j,k) - &
!                      ! (1-s(i-1,j,k))*nor(i-1,j,k,1)*s(i-1,j,k))*0.5d0*dli(1)
!          !      
!          !along y-dir
!          !
!          dsdt(i,j,k) = dsdt(i,j,k)+ &
!                        U_reg*eps*(s(i,j+1,k)-2.d0*s(i,j,k)+s(i,j-1,k))*dli(2)**2 - &
!                        U_reg*    ((1-s(i,j+1,k))*nor(i,j+1,k,2)*s(i,j+1,k) - &
!                        (1-s(i,j-1,k))*nor(i,j-1,k,2)*s(i,j-1,k))*0.5d0*dli(2)
!          !    
!          !along z-dir
!          !
!          dsdt(i,j,k)  = dsdt(i,j,k)+ &
!                         U_reg*eps*(s(i,j,k+1)-2.d0*s(i,j,k)+s(i,j,k-1))*dli(3)**2 - &
!                         U_reg*    ((1-s(i,j,k+1))*nor(i,j,k+1,3)*s(i,j,k+1) - &
!                         (1-s(i,j,k-1))*nor(i,j,k-1,3)*s(i,j,k-1))*0.5d0*dli(3)
!         !
         !
        enddo
      enddo
    enddo
    !
    return
  end subroutine regularization_psi_cdi
  !
  !
  subroutine regularization_psi(n,dli,nh_d,nh_v, &
                                dzc,dzf,u_reg,eps_int,s,ls,nor,dsdt)
    !
    implicit none
    !
    integer , intent(in ), dimension(3)                          :: n
    real(rp), intent(in ), dimension(3)                          :: dli 
    integer , intent(in )                                        :: nh_d,nh_v
    real(rp), intent(in ), dimension(1-nh_d:)                    :: dzc,dzf
    real(rp), intent(in )                                        :: u_reg,eps_int
    real(rp), intent(in ), dimension(1-nh_v:,1-nh_v:,1-nh_v:)    :: s
    real(rp), intent(in ), dimension(1-nh_v:,1-nh_v:,1-nh_v:)    :: ls
    real(rp), intent(in ), dimension(1-nh_v:,1-nh_v:,1-nh_v:,1:) :: nor
    real(rp), intent(out), dimension(     1:,     1:,     1:)    :: dsdt
    !
    real(rp) :: ls_face,term_p,term_m,nor_face 
    integer  :: i,j,k
    !
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          !
          ! along x
          !
          ls_face  = 0.5_rp*(ls(i+1,j,k)+ls(i,j,k))
          nor_face = 0.5_rp*(nor(i+1,j,k,1)+nor(i,j,k,1))
          term_p   = u_reg*(eps_int*(s(i+1,j,k)-s(i,j,k))*dli(1) - &
                            0.25_rp*( 1._rp-tanh(ls_face/(2._rp*eps_int))**2 )*nor_face)
          ls_face  = 0.5_rp*(ls(i-1,j,k)+ls(i,j,k))
          nor_face = 0.5_rp*(nor(i-1,j,k,1)+nor(i,j,k,1))
          term_m   = u_reg*(eps_int*(s(i,j,k)-s(i-1,j,k))*dli(1) - &
                            0.25_rp*( 1._rp-tanh(ls_face/(2._rp*eps_int))**2 )*nor_face)
          !
          dsdt(i,j,k) = (term_p-term_m)*dli(1)
          !
          ! along y
          !
          ls_face  = 0.5_rp*(ls(i,j+1,k)+ls(i,j,k))
          nor_face = 0.5_rp*(nor(i,j+1,k,2)+nor(i,j,k,2))
          term_p   = u_reg*(eps_int*(s(i,j+1,k)-s(i,j,k))*dli(2) - &
                            0.25_rp*( 1._rp-tanh(ls_face/(2._rp*eps_int))**2 )*nor_face)
          ls_face  = 0.5*(ls(i,j-1,k)+ls(i,j,k))
          nor_face = 0.5*(nor(i,j-1,k,2)+nor(i,j,k,2))
          term_m   = u_reg*(eps_int*(s(i,j,k)-s(i,j-1,k))*dli(2) - &
                            0.25_rp*( 1._rp-tanh(ls_face/(2._rp*eps_int))**2 )*nor_face)
          !
          dsdt(i,j,k) = dsdt(i,j,k) + (term_p-term_m)*dli(2)
          !
          ! along z
          !
          ls_face  = 0.5_rp*(ls(i,j,k+1)+ls(i,j,k))
          nor_face = 0.5_rp*(nor(i,j,k+1,3)+nor(i,j,k,3))
          term_p   = u_reg*(eps_int*(s(i,j,k+1)-s(i,j,k))*dli(3) - &
                            0.25_rp*( 1._rp-tanh(ls_face/(2._rp*eps_int))**2 )*nor_face)
          ls_face  = 0.5_rp*(ls(i,j,k-1)+ ls(i,j,k))
          nor_face = 0.5_rp*(nor(i,j,k-1,3)+nor(i,j,k,3))
          term_m   = u_reg*(eps_int*(s(i,j,k)-s(i,j,k-1))*dli(3) - &
                            0.25_rp*( 1._rp-tanh(ls_face/(2._rp*eps_int))**2 )*nor_face)
          !
          dsdt(i,j,k) = dsdt(i,j,k) + (term_p-term_m)*dli(3)
          !
        enddo
      enddo
    enddo
    !
    return
  end subroutine regularization_psi
  !
  subroutine psi_to_ls(n,dl,nh_p,psi,eps_int,ls)
    !
    implicit none
    !
    integer , intent(in ), dimension(3) :: n
    real(rp), intent(in ), dimension(3) :: dl
    integer , intent(in )               :: nh_p
    real(rp), intent(in ), dimension(1-nh_p:,1-nh_p:,1-nh_p:) :: psi
    real(rp), intent(in )                                     :: eps_int
    real(rp), intent(out), dimension(1-nh_p:,1-nh_p:,1-nh_p:) :: ls
    !
    real(rp) :: clip_psi
    integer  :: i,j,k
    !
    real(rp), parameter :: small = real(2.e-16,rp)
    !
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          !
          clip_psi  = min(max(0.0_rp,psi(i,j,k)),1.0_rp) 
          ls(i,j,k) = eps_int*log( (clip_psi+small)/(1.0_rp-clip_psi+small) )
          !
        enddo
      enddo
    enddo
    !
    return
  end subroutine psi_to_ls
  !
  subroutine ls_to_psi(n,dl,nh_p,ls,eps_int,psi)
    !
    implicit none
    !
    integer , intent(in ), dimension(3) :: n
    real(rp), intent(in ), dimension(3) :: dl
    integer , intent(in )               :: nh_p
    real(rp), intent(in ), dimension(1-nh_p:,1-nh_p:,1-nh_p:) :: ls
    real(rp), intent(in )                                     :: eps_int
    real(rp), intent(out), dimension(1-nh_p:,1-nh_p:,1-nh_p:) :: psi
    !
    integer :: i,j,k
    !
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          psi(i,j,k) = 0.5_rp*( 1._rp+tanh(ls(i,j,k)/(2._rp*eps_int)) )
        enddo
      enddo
    enddo
    !
    return
  end subroutine ls_to_psi
  !
  subroutine cmpt_norm(n,dli,nh_p,phi,kappa,nor)
    !
    use mod_boiling      , only: max_kappa
    !
    implicit none
    !
    integer , intent(in ), dimension(3) :: n
    real(rp), intent(in ), dimension(3) :: dli
    integer , intent(in )               :: nh_p
    real(rp), intent(in ), dimension(1-nh_p:,1-nh_p:,1-nh_p:)    :: phi
    real(rp), intent(out), dimension(     0:,     0:,     0:)    :: kappa
    real(rp), intent(out), dimension(1-nh_p:,1-nh_p:,1-nh_p:,1:) :: nor
    !
    real(rp), dimension(8) :: nx,ny,nz,mx,my,mz
    real(rp), dimension(3) :: dl
    real(rp) :: norm,cur_x,cur_y,cur_z
    integer  :: i,j,k,p
    !
    real(rp), parameter :: small = real(1.0e-12,rp)
    !
    !
    dl(:) = dli(:)**(-1)
    !
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          !
          !i+1/2 j+1/2 k+1/2
          mx(1)=((phi(i+1,j  ,k  )+phi(i+1,j+1,k  )+phi(i+1,j  ,k+1)+phi(i+1,j+1,k+1))-&
                 (phi(i  ,j  ,k  )+phi(i  ,j+1,k  )+phi(i  ,j  ,k+1)+phi(i  ,j+1,k+1)))*dli(1)*0.25_rp
          !i+1/2 j-1/2 k+1/2
          mx(2)=((phi(i+1,j  ,k  )+phi(i+1,j-1,k  )+phi(i+1,j  ,k+1)+phi(i+1,j-1,k+1))-&
                 (phi(i  ,j  ,k  )+phi(i  ,j-1,k  )+phi(i  ,j  ,k+1)+phi(i  ,j-1,k+1)))*dli(1)*0.25_rp
          !i+1/2 j+1/2 k-1/2
          mx(3)=((phi(i+1,j  ,k  )+phi(i+1,j+1,k  )+phi(i+1,j  ,k-1)+phi(i+1,j+1,k-1))-&
                 (phi(i  ,j  ,k  )+phi(i  ,j+1,k  )+phi(i  ,j  ,k-1)+phi(i  ,j+1,k-1)))*dli(1)*0.25_rp
          !i+1/2 j-1/2 k-1/2
          mx(4)=((phi(i+1,j  ,k  )+phi(i+1,j-1,k  )+phi(i+1,j  ,k-1)+phi(i+1,j-1,k-1))-&
                 (phi(i  ,j  ,k  )+phi(i  ,j-1,k  )+phi(i  ,j  ,k-1)+phi(i  ,j-1,k-1)))*dli(1)*0.25_rp
          !i-1/2 j+1/2 k+1/2
          mx(5)=((phi(i  ,j  ,k  )+phi(i  ,j+1,k  )+phi(i  ,j  ,k+1)+phi(i  ,j+1,k+1))-&
                 (phi(i-1,j  ,k  )+phi(i-1,j+1,k  )+phi(i-1,j  ,k+1)+phi(i-1,j+1,k+1)))*dli(1)*0.25_rp
          !i-1/2 j-1/2 k+1/2
          mx(6)=((phi(i  ,j  ,k  )+phi(i  ,j-1,k  )+phi(i  ,j  ,k+1)+phi(i  ,j-1,k+1))-&
                 (phi(i-1,j  ,k  )+phi(i-1,j-1,k  )+phi(i-1,j  ,k+1)+phi(i-1,j-1,k+1)))*dli(1)*0.25_rp
          !i-1/2 j+1/2 k-1/2
          mx(7)=((phi(i  ,j  ,k  )+phi(i  ,j+1,k  )+phi(i  ,j  ,k-1)+phi(i  ,j+1,k-1))-&
                 (phi(i-1,j  ,k  )+phi(i-1,j+1,k  )+phi(i-1,j  ,k-1)+phi(i-1,j+1,k-1)))*dli(1)*0.25_rp
          !i-1/2 j-1/2 k-1/2
          mx(8)=((phi(i  ,j  ,k  )+phi(i  ,j-1,k  )+phi(i  ,j  ,k-1)+phi(i  ,j-1,k-1))-&
                 (phi(i-1,j  ,k  )+phi(i-1,j-1,k  )+phi(i-1,j  ,k-1)+phi(i-1,j-1,k-1)))*dli(1)*0.25_rp
          !
          !i+1/2 j+1/2 k+1/2
          my(1)=((phi(i  ,j+1,k  )+phi(i+1,j+1,k  )+phi(i  ,j+1,k+1)+phi(i+1,j+1,k+1))-&
                 (phi(i  ,j  ,k  )+phi(i+1,j  ,k  )+phi(i  ,j  ,k+1)+phi(i+1,j  ,k+1)))*dli(2)*0.25_rp
          !i+1/2 j-1/2 k+1/2
          my(2)=((phi(i  ,j  ,k  )+phi(i+1,j  ,k  )+phi(i  ,j  ,k+1)+phi(i+1,j  ,k+1))-&
                 (phi(i  ,j-1,k  )+phi(i+1,j-1,k  )+phi(i  ,j-1,k+1)+phi(i+1,j-1,k+1)))*dli(2)*0.25_rp
          !i+1/2 j+1/2 k-1/2
          my(3)=((phi(i  ,j+1,k  )+phi(i+1,j+1,k  )+phi(i  ,j+1,k-1)+phi(i+1,j+1,k-1))-&
                 (phi(i  ,j  ,k  )+phi(i+1,j  ,k  )+phi(i  ,j  ,k-1)+phi(i+1,j  ,k-1)))*dli(2)*0.25_rp
          !i+1/2 j-1/2 k-1/2
          my(4)=((phi(i  ,j  ,k  )+phi(i+1,j  ,k  )+phi(i  ,j  ,k-1)+phi(i+1,j  ,k-1))-&
                 (phi(i  ,j-1,k  )+phi(i+1,j-1,k  )+phi(i  ,j-1,k-1)+phi(i+1,j-1,k-1)))*dli(2)*0.25_rp
          !i-1/2 j+1/2 k+1/2
          my(5)=((phi(i  ,j+1,k  )+phi(i-1,j+1,k  )+phi(i  ,j+1,k+1)+phi(i-1,j+1,k+1))-&
                 (phi(i  ,j  ,k  )+phi(i-1,j  ,k  )+phi(i  ,j  ,k+1)+phi(i-1,j  ,k+1)))*dli(2)*0.25_rp
          !i-1/2 j-1/2 k+1/2
          my(6)=((phi(i  ,j  ,k  )+phi(i-1,j  ,k  )+phi(i  ,j  ,k+1)+phi(i-1,j  ,k+1))-&
                 (phi(i  ,j-1,k  )+phi(i-1,j-1,k  )+phi(i  ,j-1,k+1)+phi(i-1,j-1,k+1)))*dli(2)*0.25_rp
          !i-1/2 j+1/2 k-1/2
          my(7)=((phi(i  ,j+1,k  )+phi(i-1,j+1,k  )+phi(i  ,j+1,k-1)+phi(i-1,j+1,k-1))-&
                 (phi(i  ,j  ,k  )+phi(i-1,j  ,k  )+phi(i  ,j  ,k-1)+phi(i-1,j  ,k-1)))*dli(2)*0.25_rp
          !i-1/2 j-1/2 k-1/2
          my(8)=((phi(i  ,j  ,k  )+phi(i-1,j  ,k  )+phi(i  ,j  ,k-1)+phi(i-1,j  ,k-1))-&
                 (phi(i  ,j-1,k  )+phi(i-1,j-1,k  )+phi(i  ,j-1,k-1)+phi(i-1,j-1,k-1)))*dli(2)*0.25_rp
          !
          !i+1/2 j+1/2 k+1/2
          mz(1)=((phi(i  ,j  ,k+1)+phi(i+1,j  ,k+1)+phi(i  ,j+1,k+1)+phi(i+1,j+1,k+1))-&
                 (phi(i  ,j  ,k  )+phi(i+1,j  ,k  )+phi(i  ,j+1,k  )+phi(i+1,j+1,k  )))*dli(3)*0.25_rp
          !i+1/2 j-1/2 k+1/2
          mz(2)=((phi(i  ,j  ,k+1)+phi(i+1,j  ,k+1)+phi(i  ,j-1,k+1)+phi(i+1,j-1,k+1))-&
                 (phi(i  ,j  ,k  )+phi(i+1,j  ,k  )+phi(i  ,j-1,k  )+phi(i+1,j-1,k  )))*dli(3)*0.25_rp
          !i+1/2 j+1/2 k-1/2
          mz(3)=((phi(i  ,j  ,k  )+phi(i+1,j  ,k  )+phi(i  ,j+1,k  )+phi(i+1,j+1,k  ))-&
                 (phi(i  ,j  ,k-1)+phi(i+1,j  ,k-1)+phi(i  ,j+1,k-1)+phi(i+1,j+1,k-1)))*dli(3)*0.25_rp
          !i+1/2 j-1/2 k-1/2
          mz(4)=((phi(i  ,j  ,k  )+phi(i+1,j  ,k  )+phi(i  ,j-1,k  )+phi(i+1,j-1,k  ))-&
                 (phi(i  ,j  ,k-1)+phi(i+1,j  ,k-1)+phi(i  ,j-1,k-1)+phi(i+1,j-1,k-1)))*dli(3)*0.25_rp
          !i-1/2 j+1/2 k+1/2
          mz(5)=((phi(i  ,j  ,k+1)+phi(i-1,j  ,k+1)+phi(i  ,j+1,k+1)+phi(i-1,j+1,k+1))-&
                 (phi(i  ,j  ,k  )+phi(i-1,j  ,k  )+phi(i  ,j+1,k  )+phi(i-1,j+1,k  )))*dli(3)*0.25_rp
          !i-1/2 j-1/2 k+1/2
          mz(6)=((phi(i  ,j  ,k+1)+phi(i-1,j  ,k+1)+phi(i  ,j-1,k+1)+phi(i-1,j-1,k+1))-&
                 (phi(i  ,j  ,k  )+phi(i-1,j  ,k  )+phi(i  ,j-1,k  )+phi(i-1,j-1,k  )))*dli(3)*0.25_rp
          !i-1/2 j+1/2 k-1/2
          mz(7)=((phi(i  ,j  ,k  )+phi(i-1,j  ,k  )+phi(i  ,j+1,k  )+phi(i-1,j+1,k  ))-&
                 (phi(i  ,j  ,k-1)+phi(i-1,j  ,k-1)+phi(i  ,j+1,k-1)+phi(i-1,j+1,k-1)))*dli(3)*0.25_rp
          !i-1/2 j-1/2 k-1/2
          mz(8)=((phi(i  ,j  ,k  )+phi(i-1,j  ,k  )+phi(i  ,j-1,k  )+phi(i-1,j-1,k  ))-&
                 (phi(i  ,j  ,k-1)+phi(i-1,j  ,k-1)+phi(i  ,j-1,k-1)+phi(i-1,j-1,k-1)))*dli(3)*0.25_rp
          !
          do p=1,8
            norm  = sqrt(mx(p)**2+my(p)**2+mz(p)**2) + small
            nx(p) = mx(p)/norm
            ny(p) = my(p)/norm
            nz(p) = mz(p)/norm
          enddo
          !
          ! compute the normal vector
          !
          nor(i,j,k,1) = 0._rp
          nor(i,j,k,2) = 0._rp
          nor(i,j,k,3) = 0._rp
          do p=1,8
            nor(i,j,k,1) = nor(i,j,k,1) + 0.125_rp*mx(p)
            nor(i,j,k,2) = nor(i,j,k,2) + 0.125_rp*my(p)
            nor(i,j,k,3) = nor(i,j,k,3) + 0.125_rp*mz(p)
          enddo
          !
#if defined(_TWOD)
          nor(i,j,k,1) = 0.0
#endif
          norm         = sqrt(nor(i,j,k,1)**2+nor(i,j,k,2)**2+nor(i,j,k,3)**2) + small
          nor(i,j,k,1) = nor(i,j,k,1)/norm
          nor(i,j,k,2) = nor(i,j,k,2)/norm
          nor(i,j,k,3) = nor(i,j,k,3)/norm
          ! 
          cur_x = ((nx(1)+nx(2)+nx(3)+nx(4))-(nx(5)+nx(6)+nx(7)+nx(8)))*dl(1)*0.25_rp
          cur_y = ((ny(1)+ny(3)+ny(5)+ny(7))-(ny(2)+ny(4)+ny(6)+ny(8)))*dl(2)*0.25_rp
          cur_z = ((nz(1)+nz(2)+nz(5)+nz(6))-(nz(3)+nz(4)+nz(7)+nz(8)))*dl(3)*0.25_rp
          !
          !
          !
          kappa(i,j,k) = -(cur_x*dli(1)**2+cur_y*dli(2)**2+cur_z*dli(3)**2) !max(-max_kappa, min(-(cur_x*dli(1)**2+cur_y*dli(2)**2+cur_z*dli(3)**2),max_kappa))
     
          !
        enddo
      enddo
    enddo
    
   
    
    
    !
    return
  end subroutine cmpt_norm
  !
  subroutine update_property(n,nh_p,prop12,psi,prop)
    ! Assume that psi is fully comminucated by now
    ! no commnuication for fluid property needed after this subroutine
    implicit none
    !
    integer , intent(in ), dimension(3)                       :: n
    integer , intent(in )                                     :: nh_p
    real(rp), intent(in ), dimension(2)                       :: prop12
    real(rp), intent(in ), dimension(1-nh_p:,1-nh_p:,1-nh_p:) :: psi
    real(rp), intent(out), dimension(1-nh_p:,1-nh_p:,1-nh_p:) :: prop
    !
    real(rp) :: prop12_1, prop12_2
    integer  :: i,j,k
    integer  :: n1, n2, n3
    ! 
    n1 = n(1)
    n2 = n(2)
    n3 = n(3)
    !
    prop12_1 = prop12(1)
    prop12_2 = prop12(2)
    !
    !$acc parallel loop collapse(3)
    do k=1-nh_p,n3+nh_p
      do j=1-nh_p,n2+nh_p
        do i=1-nh_p,n1+nh_p
          prop(i,j,k) = psi(i,j,k)*prop12_1+(1.0_rp-psi(i,j,k))*prop12_2
        enddo
      enddo
    enddo
    !$acc end parallel loop 
    !
    return
  end subroutine update_property
  !
  subroutine initls(inils,n,dl,l,nh_p,ls)
    !
    use mod_param, only: xc,yc,zc,r,nbub
    !
    implicit none
    !
    character(len=3), intent(in )                                     :: inils
    integer         , intent(in ), dimension(3)                       :: n
    real(rp)        , intent(in ), dimension(3)                       :: dl
    real(rp)        , intent(in ), dimension(3)                       :: l
    integer         , intent(in )                                     :: nh_p
    real(rp)        , intent(out), dimension(1-nh_p:,1-nh_p:,1-nh_p:) :: ls
    !
    real(rp) :: phi
    real(rp) :: d,x,y,z,ry
    real(rp) :: rho1,rho2
    integer  :: i,j,k
    !
    !ry = 0.5d0*d0(1)
    ry = r(1)
    !
    select case(inils)
    case('tgn')
      do k=1,n(3)
        z = (k+ijk_start(3)-0.5)*dl(3)
        do j=1,n(2)
          y = (j+ijk_start(2)-0.5)*dl(2)
          do i=1,n(1)
            x = (i+ijk_start(1)-0.5)*dl(1)
            !
#if defined(_TWOD)
            d = sqrt(                (y-yc(1))**2.0+(z-zc(1))**2.0 )
#else
            d = sqrt( (x-xc(1))**2.0+(y-yc(1))**2.0+(z-zc(1))**2.0 )
#endif
            !
            ls(i,j,k) = ry-d
            !
          enddo
        enddo
      enddo
    case default
      call flutas_error('ERROR: invalid name for initial level-set field: check di_acdi.in')
    end select
    !
    return
  end subroutine initls
  !
  subroutine initpsi(inipsi,n,dl,l,nh_v,eps_int,psi,tmp,cond,ent_yesno,ent,pres,trigger_update,time)
    !
    use mod_param     , only: xc,yc,zc,r,nbub,h01,h02,rho1,rho2,IBM_h,cp1,cp2,sigma, &
                              gacc_x,gacc_y,gacc_z,pi,heater_power,IBM_rho,IBM_cp
    use mod_boiling   , only: nuc_sites,num_sites,nuc_rad,max_kappa,heater_ind,k_max_nuc, &
                              nuc_sites_chk,nuc_sites_consider,site_activation,nuc_sites_consider_dir
    !
    implicit none
    !
    character(len=3), intent(in   )                                     :: inipsi
    integer         , intent(in   ), dimension(3)                       :: n
    real(rp)        , intent(in   ), dimension(3)                       :: dl
    real(rp)        , intent(in   ), dimension(3)                       :: l
    integer         , intent(in   )                                     :: nh_v
    real(rp)        , intent(in   )                                     :: eps_int,time
    real(rp)        , intent(in   ), dimension(1-nh_v:,1-nh_v:,1-nh_v:) :: tmp,cond
    real(rp)        , intent(inout), dimension(1-nh_v:,1-nh_v:,1-nh_v:) :: psi
    logical         , intent(in   )                                     :: ent_yesno
    real(rp)        , intent(inout), dimension(1-nh_v:,1-nh_v:,1-nh_v:) :: ent
    real(rp)        , intent(inout), dimension(1-nh_v:,1-nh_v:,1-nh_v:) :: pres
    logical         , intent(out  )                                     :: trigger_update
    !
    real(rp) :: d,x,y,z,diffent,diff_pres,A,C,p_loc,tanh_val,rgh,d_psi,tmp_var,limiter
    real(rp) :: xmin,xmax,ymin,ymax,zmin,zmax,gamma3,loctmp,locpsi,gamma1,gamma2
    real(rp) :: loctmp1,loctmp2,loctmp3,loctmp4,locls,clip_psi
    real(rp) :: ibm_switch,ent_switch,t_corr
    integer  :: i,j,k,p,k_corr,p_chk,is_active,p_consider
    real(rp), parameter :: small = real(2.e-16,rp)
    !
    !
    trigger_update = .FALSE.
    select case(inipsi)
    case('bub')
      do k=IBM_h+1,n(3) !
        z = (k+ijk_start(3)-0.5_rp)*dl(3)
        do j=1,n(2)
          y = (j+ijk_start(2)-0.5_rp)*dl(2)
          do i=1,n(1)
            x = (i+ijk_start(1)-0.5_rp)*dl(1)
             !
             psi(i,j,k) = 0._rp
             do p=1,nbub
#if defined(_TWOD)
               d = sqrt(                  (y-yc(p))**2._rp+(z-zc(p))**2._rp )
#else
               d = sqrt( (x-xc(p))**2._rp+(y-yc(p))**2._rp+(z-zc(p))**2._rp )
#endif
               psi(i,j,k) = psi(i,j,k) + 0.5_rp*(1._rp+tanh( (r(p)-d)/(2._rp*eps_int) ))
               if     (gacc_y.NE.0.0) then
                rgh = gacc_y*(rho2-rho1)*(yc(p)-y)
               elseif (gacc_z.NE.0.0) then
                rgh = gacc_z*(rho2-rho1)*(zc(p)-z)
               elseif (gacc_x.NE.0.0) then
                rgh = gacc_x*(rho2-rho1)*(xc(p)-x)
               endif
#if defined(_TWOD)
               diff_pres = sigma*min(1.0/r(p),max_kappa)
#else
               diff_pres = 2.0*sigma*min(1.0/r(p),max_kappa)
#endif
               diff_pres = diff_pres
               tanh_val = 0.5_rp*(1._rp+tanh( (r(p)-d)/(2._rp*eps_int) ))
               p_loc = tanh_val*diff_pres
               pres(i,j,k) = pres(i,j,k) +  p_loc
               !presold(i,j,k) = presold(i,j,k) +  p_loc
             enddo
             !
          enddo
        enddo
      enddo
    case('sfn')
      do k=1,n(3) !
        z = (k+ijk_start(3)-0.5_rp)*dl(3)
        do j=1,n(2)
          y = (j+ijk_start(2)-0.5_rp)*dl(2)
          do i=1,n(1)
            x = (i+ijk_start(1)-0.5_rp)*dl(1)
             !
             psi(i,j,k) = 0._rp
             do p=1,nbub
#if defined(_TWOD)
               d = (z-zc(p))
#else
               d = (z-zc(p))
#endif
               psi(i,j,k) = psi(i,j,k) + 0.5_rp*(1._rp+tanh( (r(p)-d)/(2._rp*eps_int) ))
             enddo
             !
          enddo
        enddo
      enddo
    case('one')
      do k=1,n(3) !
        z = (k+ijk_start(3)-0.5_rp)*dl(3)
        do j=1,n(2)
          y = (j+ijk_start(2)-0.5_rp)*dl(2)
          do i=1,n(1)
            x = (i+ijk_start(1)-0.5_rp)*dl(1)
             !
             psi(i,j,k) = 0._rp
             do p=1,nbub
               d = abs(z-zc(p))
               psi(i,j,k) = psi(i,j,k) + 0.5_rp*(1._rp+tanh( (r(p)-d)/(2._rp*eps_int) ))
             enddo
             !
          enddo
        enddo
      enddo
    case('flm')
      do k=1,n(3) !
        z = (k+ijk_start(3)-0.5_rp)*dl(3)
        do j=1,n(2)
          y = (j+ijk_start(2)-0.5_rp)*dl(2)
          do i=1,n(1)
            x = (i+ijk_start(1)-0.5_rp)*dl(1)
             !
             psi(i,j,k) = 0._rp
             d = 0.07868440885009152  !Lambda_0 in Sun et al 2014
             !d = 0.0023003957362390023  !Lambda_0 in NTNU_test_4
             d = (d/128.0)*(4.0+cos(2.0*pi*(y-0.5*d)/d))           ! in Sun et al 2014
             !d = (d/1.0)*(0.025+0.0025*cos(2.0*pi*(y-0.5*d)/d))    ! in NTNU_test_4
             psi(i,j,k) = psi(i,j,k) + 0.5_rp*(1._rp+tanh( (d-z)/(2._rp*eps_int) ))
             !
          enddo
        enddo
      enddo
    case('nuc')
      ent_switch = 0.0_rp
      if(ent_yesno) ent_switch = 1.0_rp
      ibm_switch = 0.0_rp
      if(IBM_h.NE.0) ibm_switch = 1.0_rp
      ! Default activation
      site_activation(:) = 0
      ! First: check which nucleation site is activated on this rank
      if (allocated(nuc_sites_chk)) then
        do p_chk=1,size(nuc_sites_chk)
          p = nuc_sites_chk(p_chk)
          is_active =   0
          loctmp    =   0.0
          locpsi    = 999.0
          ! Now the color function and the temperature on the nucleation site:
          i = floor(nuc_sites(p,1)/dl(1) - ijk_start(1) + 0.5)
          j = floor(nuc_sites(p,2)/dl(2) - ijk_start(2) + 0.5)
          k = floor(nuc_sites(p,3)/dl(3) - ijk_start(3) + 0.5)
          !
          gamma1  = ((i+ijk_start(1)+0.5_rp)*dl(1) - nuc_sites(p,1))/dl(1)
          gamma2  = ((j+ijk_start(2)+0.5_rp)*dl(2) - nuc_sites(p,2))/dl(2)
          gamma3  = ((k+ijk_start(3)+0.5_rp)*dl(3) - nuc_sites(p,3))/dl(3)
          !
          ! nucleation site has to be on cell face! check if this is true:
          !if(nint(gamma3*100000).NE.50000) print*, 'Nucleation site no. ' , p , ' is not on a cell face!!!'
          !
          ! Now calculate the wall temperature and interpolated color function
          !
          if((time.EQ.0.0_rp).OR.(IBM_h.EQ.0))then
            loctmp =  (      gamma1)*(      gamma2)*(      gamma3)*tmp(i  ,j  ,k  ) &
                     +(1.0 - gamma1)*(      gamma2)*(      gamma3)*tmp(i+1,j  ,k  ) &
                     +(      gamma1)*(1.0 - gamma2)*(      gamma3)*tmp(i  ,j+1,k  ) &
                     +(1.0 - gamma1)*(1.0 - gamma2)*(      gamma3)*tmp(i+1,j+1,k  ) &
                     +(      gamma1)*(      gamma2)*(1.0 - gamma3)*tmp(i  ,j  ,k+1) &
                     +(1.0 - gamma1)*(      gamma2)*(1.0 - gamma3)*tmp(i+1,j  ,k+1) &
                     +(      gamma1)*(1.0 - gamma2)*(1.0 - gamma3)*tmp(i  ,j+1,k+1) &
                     +(1.0 - gamma1)*(1.0 - gamma2)*(1.0 - gamma3)*tmp(i+1,j+1,k+1)
          else
            loctmp1 = (1.0/(cond(i,j,k)+cond(i,j,k+1))) &
                    *(tmp(i,j,k+1)*cond(i,j,k+1)+tmp(i,j,k)*cond(i,j,k))  &
                    +(heater_power*dl(3)*heater_ind(i,j))          &
                    /(2.0_rp*(cond(i,j,k)+cond(i,j,k+1)))
            loctmp2 = (1.0/(cond(i+1,j,k)+cond(i+1,j,k+1))) &
                    *(tmp(i+1,j,k+1)*cond(i+1,j,k+1)+tmp(i+1,j,k)*cond(i+1,j,k))  &
                    +(heater_power*dl(3)*heater_ind(i+1,j))          &
                    /(2.0_rp*(cond(i+1,j,k)+cond(i+1,j,k+1)))
            loctmp3 = (1.0/(cond(i,j+1,k)+cond(i,j+1,k+1))) &
                    *(tmp(i,j+1,k+1)*cond(i,j+1,k+1)+tmp(i,j+1,k)*cond(i,j+1,k))  &
                    +(heater_power*dl(3)*heater_ind(i,j+1))          &
                    /(2.0_rp*(cond(i,j+1,k)+cond(i,j+1,k+1)))
            loctmp4 = (1.0/(cond(i+1,j+1,k)+cond(i+1,j+1,k+1))) &
                    *(tmp(i+1,j+1,k+1)*cond(i+1,j+1,k+1)+tmp(i+1,j+1,k)*cond(i+1,j+1,k))  &
                    +(heater_power*dl(3)*heater_ind(i+1,j+1))          &
                    /(2.0_rp*(cond(i+1,j+1,k)+cond(i+1,j+1,k+1)))
            loctmp  = (      gamma1)*(      gamma2)*loctmp1 &
                     +(1.0 - gamma1)*(      gamma2)*loctmp2 &
                     +(      gamma1)*(1.0 - gamma2)*loctmp3 &
                     +(1.0 - gamma1)*(1.0 - gamma2)*loctmp4
          endif
          !
          gamma3 = 0.0
          locpsi = (      gamma1)*(      gamma2)*(      gamma3)*psi(i  ,j  ,k  ) &
                  +(1.0 - gamma1)*(      gamma2)*(      gamma3)*psi(i+1,j  ,k  ) &
                  +(      gamma1)*(1.0 - gamma2)*(      gamma3)*psi(i  ,j+1,k  ) &
                  +(1.0 - gamma1)*(1.0 - gamma2)*(      gamma3)*psi(i+1,j+1,k  ) &
                  +(      gamma1)*(      gamma2)*(1.0 - gamma3)*psi(i  ,j  ,k+1) &
                  +(1.0 - gamma1)*(      gamma2)*(1.0 - gamma3)*psi(i+1,j  ,k+1) &
                  +(      gamma1)*(1.0 - gamma2)*(1.0 - gamma3)*psi(i  ,j+1,k+1) &
                  +(1.0 - gamma1)*(1.0 - gamma2)*(1.0 - gamma3)*psi(i+1,j+1,k+1)
          clip_psi  = min(max(0.0_rp,locpsi),1.0_rp) 
          locls = eps_int*log( (clip_psi+small)/(1.0_rp-clip_psi+small) )
          ! should the nucleation cite be activated?
          if((loctmp.GE.nuc_sites(p,4)).AND. &
             !(locpsi.LT.0.0001)) then
              (locls.LT.(-6.0*dl(1)))) then
            is_active = 1
          endif
          !
          if(is_active.EQ.1) print*, 'Nucleation site no. '                , p     , &
                                                        ' has been activated at surface temp: ', loctmp, &
                                                        ' and phase field value '              , locpsi, &
                                                        ' at time '                            , time
          site_activation(p) = is_active
        enddo
      endif
      !
      ! MPI: which nucleation sites where activated else where?
      call mpi_allreduce(MPI_IN_PLACE,site_activation,num_sites,MPI_INTEGER,mpi_max,comm_cart,ierr)
      !
      if(sum(site_activation).NE.0) trigger_update = .TRUE.
      !
      ! Now inizialize all active nucleation cites
      !
      do k=IBM_h+1,k_max_nuc
        z = (k+ijk_start(3)-0.5_rp)*dl(3)
        do j=1,n(2)
          y = (j+ijk_start(2)-0.5_rp)*dl(2)
          do i=1,n(1)
            x = (i+ijk_start(1)-0.5_rp)*dl(1)
             !
             do p_consider=1,size(nuc_sites_consider)
               p = nuc_sites_consider(p_consider)
               if(site_activation(p).EQ.1) then
#if defined(_TWOD)
               d = sqrt( (y-nuc_sites(p,2))**2._rp & 
                        +(z-nuc_sites(p,3))**2._rp )
#else
               d = sqrt( (x-(nuc_sites(p,1)+nuc_sites_consider_dir(p_consider)*l(1)))**2._rp &
                        +(y- nuc_sites(p,2))**2._rp &
                        +(z- nuc_sites(p,3))**2._rp )
#endif
               d_psi = 0.5_rp*(1._rp+tanh((nuc_rad-d)/(2._rp*eps_int) ))
               if(d_psi.LT.0.0001) d_psi = 0.0  ! - important for stability - 0.0001 worked well, 0.001 not tested yet!
               ! calculate the acuta change in phi
               d_psi = max(psi(i,j,k) , d_psi) - psi(i,j,k) ! this is if thwo bubble apprear too close to each other
               ! correct phi
               psi(i,j,k) = psi(i,j,k) + d_psi
               !
#if defined(_TWOD)
               diff_pres = sigma*min(1.0/nuc_rad,max_kappa)
#else
               diff_pres = 2.0*sigma*min(1.0/nuc_rad,max_kappa)
#endif
               ! calculate local pressure difference based on nucleation diameter
               !p_loc = d_psi*diff_pres
               !  Depends on surface tension model: this following one is optimal for the truncated model with a limiter of 0.75
               limiter = 0.55
               p_loc = diff_pres*d_psi! (1.0/(1.0-limiter))*(min(max(d_psi,0.5*limiter),1.0-0.5*limiter)-0.5*limiter)
               ! correct the pressure
               pres(i,j,k) = pres(i,j,k) +  p_loc
               ! Calcluate temperature increment
               tmp_var    = - d_psi*ent_switch*tmp(i,j,k)
               ! Correct the temperature
               ent(i,j,k) = ent(i,j,k) + tmp_var
               ! Calculate IBM temperature change to remove 
               ! required latent heat from wall (from nucleation)
               t_corr = - ibm_switch*((h01-h02)*rho1*d_psi)/(IBM_rho*IBM_cp)
               ! Calculate location where to remove the heat from
               k_corr = max(1,2*IBM_h+1-k)
               ! Correct wall heat:
               ent(i,j,k_corr) = ent(i,j,k_corr) + t_corr
               endif
             enddo
             !
          enddo
        enddo
      enddo
    case default
      call flutas_error('ERROR: invalid name for initial DI field: check di_acdi.in')
    end select
    !
    return
  end subroutine initpsi
    !
subroutine inittmp(initmp,n,dl,l,nh_v,eps_int,tmp)
    !
    use mod_common_mpi, only: ijk_start
    !
    use mod_param, only: xc,yc,zc,r,tl0,tg0,nbub
    !
    implicit none
    !
    character(len=3), intent(in )                         :: initmp
    integer         , intent(in ), dimension(3)           :: n
    real(rp)        , intent(in ), dimension(3)           :: dl
    real(rp)        , intent(in ), dimension(3)           :: l
    integer         , intent(in )                         :: nh_v
    real(rp)        , intent(in )                         :: eps_int
    real(rp)        , intent(out), dimension(1-nh_v:,1-nh_v:,1-nh_v:) :: tmp
    !
    real(rp) :: phi
    real(rp) :: d,x,y,z,ry
    real(rp) :: rho1,rho2
    integer  :: i,j,k,p
    !
    select case(initmp)
    case('ins')
      do k=1,n(3)
        z = (k+ijk_start(3)-0.5_rp)*dl(3)
        do j=1,n(2)
          y = (j+ijk_start(2)-0.5_rp)*dl(2)
          do i=1,n(1)
            x = (i+ijk_start(1)-0.5_rp)*dl(1)
             !
             tmp(i,j,k) = 0._rp
             do p=1,nbub
#if defined(_TWOD)
               d = sqrt(                  (y-yc(p))**2._rp+(z-zc(p))**2._rp )
#else
               d = sqrt( (x-xc(p))**2._rp+(y-yc(p))**2._rp+(z-zc(p))**2._rp )
#endif
               tmp(i,j,k) = tmp(i,j,k) + 0.5_rp*(1._rp+tanh( (r(p)-d)/(2._rp*eps_int) ))*tl0
               
              ! tmp(i,j,k) = tg0+(tl0-tg0)*phi(i,j,k)
               
             enddo
             !
          enddo
        enddo
      enddo

    case default
      call flutas_error('ERROR: invalid name for initial TMP field: check di_acdi.in')
    end select
    !
    return
    !
!    ry = r(1) 
!    !
!    select case(iniy1)
!    case('tgh')
!    do k=1,n(3)
!      !
!      z = (k+ijk_start(3)-0.5)*dl(3) 
!      !
!      do j=1,n(2)
!      !
!      y = (j+ijk_start(2)-0.5)*dl(2)
!      !
!        do i=1,n(1)
!      !
!      x = (i+ijk_start(1)-0.5)*dl(1)
!          !
!          a1(i,j,k) = 0.5d0*(1-tanh((sqrt((x-xc(1))**2+(y-yc(1))**2)-ry)/(1.d0*eps_int)))
!          !
!        enddo
!      enddo
!    enddo
!    end select
!    !
!    return
  end subroutine inittmp
    !
    subroutine initc1(iniy1,n,dl,l,nh_v,eps_int,phi,c1,c2)
    !
    use mod_common_mpi, only: ijk_start
    !
    use mod_param, only: xc,yc,zc,r
    !
    implicit none
    !
    character(len=3), intent(in )                         :: iniy1
    integer         , intent(in ), dimension(3)           :: n
    real(rp)        , intent(in ), dimension(3)           :: dl
    real(rp)        , intent(in ), dimension(3)           :: l
    integer         , intent(in )                         :: nh_v
    real(rp)        , intent(in )                         :: eps_int
    real(rp)        , intent(in ), dimension(1-nh_v:,1-nh_v:,1-nh_v:)  :: phi
    real(rp)        , intent(out), dimension(1-nh_v:,1-nh_v:,1-nh_v:) :: c1,c2
    !
    real(rp) :: d,x,y,z,ry
    integer  :: i,j,k
    !
    ry = r(1)
    !
    select case(iniy1)
    case ('lin')
    do k=1,n(3)
      !
      z = (k+ijk_start(3)-0.5)*dl(3) 
      !
      do j=1,n(2)
      !
      y = (j+ijk_start(2)-0.5)*dl(2)
      !
        do i=1,n(1)
      !
      x = (i+ijk_start(1)-0.5)*dl(1)
      !
     ! if (phi(i,j,k).ne.1.0)then
      if (phi(i,j,k).lt.0.999d0.and.phi(i,j,k).gt.0.001d0)then
      c2(i,j,k) = 0.1        
     ! if (phi(i,j,k).lt.0.999d0.and.phi(i,j,k).gt.0.001d0)then
    !  c2(i,j,k) = 0.0!0.028*phi(i,j,k)*(1.0-phi(i,j,k))!1.d0-phi(i,j,k)!*(1-phi(i,j,k))!*0.028d0
     ! else
!      c2(i,j,k) = 0.d0
 !     endif
      else
      c2(i,j,k) = 0.d0
      endif
      
        enddo
      enddo
    enddo
    c1(:,:,:) = 0.d0
    end select
    !
    return
    !
  end subroutine initc1
    !
end module mod_di_acdi
