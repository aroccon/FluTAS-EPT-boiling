!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@(     #@@@@@@@@@@#        &@        @@                  @
!@@@@@@@@@@@@@@@@#/@@@@@      #@@@@@@@@         &@@@        @@                  @
!@@@@@@@@@@@@(//////@@@@@     #@@@@@@         @@@@@@        @@                  @
!@@@@@@@@@@@@@(//////@@@@#    #@@@&         @@@@@@@@        @@@@@@@        @@@@@@
!@@@@@@@/&@@@@@@#/////&@@@    #@,         @@@@@@@@@@        @@@@@@@        @@@@@@
!@@@@@(/////%@@@@@%////#@@@   #@        @@@@@@@@@@@@        @@@@@@@        @@@@@@
!@@@@//////////#@@@@&////@@/  #@         @@@@@@@@@@@        @@@@@@@        @@@@@@
!@@@@@@@@@////////(@@@&///@@  #@@#         @@@@@@@@@        @@@@@@@        @@@@@@
!@@@@@@@@@@@@@@%/////(@@@//@& #@@@@@         @@@@@@@        @@@@@@@        @@@@@@
!@&////////(#@@@@@@@@////@@/@ #@@@@@@@         &@@@@        @@@@@@@        @@@@@@
!@////////////////////(@@@#/@#%@@@@@@@@@&        &@@        @@@@@@@        @@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!#################################################
!###     Implementation based on:              ###
!### Akash Dhruv. A vortex damping outflow     ###
!### forcing for multiphase flows with sharp   ###
!### interfacial jumps. J. of Comp. Phys.,     ###
!### 511:113122, 2024. ISSN 0021-9991. doi:    ###
!### https://doi.org/10.1016/j.jcp.2024.113122 ###
!#################################################
!
! Implemented by Lorenz Weber
!
! Further relevant literature:
!
! [1] Lorenz Weber, Aritra Mukherjee, Andreas G. Class, and Luca Brandt. 
! A consistent and scalable framework suitable for boiling flows using the 
! conservative diffuse interface method. Journal of Computational Physics, 
! 551:114680, 2026. ISSN 0021-9991. doi: https://doi.org/10.1016/j.jcp.2026.114680
!
! [2] Lorenz Weber. Towards the Statistical Evaluation of Phase-Changing Flows: The
! Development of a Framework for Large-Scale Interface-Resolved Simulations 
! (PhD Thesis), Karlsruher Institut für Technologie (KIT), 02 2026, 
! URL https://doi.org/10.5445/IR/1000190080*
!
! There are two outflow options:
! (1)  the original outflow BC method by A. Dhruv, J. of Comp. Phys., 2024
! (2)  the modified version of (1) which is more restrictive and used for forced convection channel flows. See chapter 4.6 in [2].
!
! The switch between the two is the variable 'out_forced': True = Method (2); False = Method (1)
!
module mod_outflow
  !
  use mod_types
  use mpi
  use mod_common_mpi, only: myid,ierr,ijk_start
  use mod_sanity    , only: flutas_error
  use mod_bound     , only: updthalo
  !
  implicit none
  !
  real(rp), dimension(:,:)  , allocatable         :: psi_bound_o,ent_bound_o,mu_bound_o
  integer                                         :: out_dir,out_rk_safety
  real(rp)                                        :: u_gas,u_liquid,u_cap,out_len
  logical                                         :: out_forced
  !
  public :: init_outflow,out_dir,cmpt_avg_out_vel,corr_u_liquid,add_outflow_force, &
            outflow_psi,boundsb_arrayo,out_rk_safety,out_len,psi_bound_o,ent_bound_o,mu_bound_o
  !
  contains
  !
  !_______________________________________________________________________________
  !
  subroutine init_outflow(n)
    implicit none
    integer , intent(in   ), dimension(3)     :: n 
    !
    ! Allocate array sizes based on outflow direction
    !
    select case(out_dir)
    case(3)
      allocate(psi_bound_o(1:n(1),1:n(2)))
      allocate(mu_bound_o(1:n(1),1:n(2)))
      allocate(ent_bound_o(1:n(1),1:n(2)))
    case(-3)
      allocate(psi_bound_o(1:n(1),1:n(2)))
      allocate(mu_bound_o(1:n(1),1:n(2)))
      allocate(ent_bound_o(1:n(1),1:n(2)))
    case(2)
      allocate(psi_bound_o(1:n(1),1:n(3)))
      allocate(mu_bound_o(1:n(1),1:n(3)))
      allocate(ent_bound_o(1:n(1),1:n(3)))
    case(0)
      ! Do nothing
    case default
      call flutas_error('ERROR: This outflow direction is not implemented ... yet!')
    end select
    if(allocated(psi_bound_o)) psi_bound_o(:,:) = 0.0_rp
    if(allocated( mu_bound_o))  mu_bound_o(:,:) = 0.0_rp
    if(allocated(ent_bound_o)) ent_bound_o(:,:) = 0.0_rp
    !
    ! This variable defines the number of grid points
    ! from the outflow in which the momentum corrections 
    ! proposed in [1] in eq. 16 are switched off for 
    ! stability reasons:
    out_rk_safety = 5
    !
  end subroutine init_outflow
  !
  !_______________________________________________________________________________
  !
  subroutine cmpt_avg_out_vel(n,dl,l,nh_p,nh_u,psi,u,v,w)
    ! For Method (1): Calculates the liquid and gas velocities in the outflow region as described in A. Dhruv, J. of Comp. Phys., 2024
    ! For Method (2): No destiction between gas and liquid velocity.
    use mod_param     , only: sigma,rho1,rho2,gacc_x,gacc_y,gacc_z,IBM_h
    !
    implicit none
    integer         , intent(in ), dimension(3)                       :: n
    real(rp)        , intent(in ), dimension(3)                       :: dl,l
    integer         , intent(in )                                     :: nh_p,nh_u
    real(rp)        , intent(in ), dimension(1-nh_p:,1-nh_p:,1-nh_p:) :: psi
    real(rp)        , intent(in ), dimension(1-nh_u:,1-nh_u:,1-nh_u:) :: u,v,w
    real(rp)                                                          :: x,y,z,h_val
    real(rp)                                                          :: cache_g1,cache_g2
    real(rp)                                                          :: phiz,l_cap,phiy
    real(rp)                                                          :: cache_l1,cache_l2
    integer                                                           :: i,j,k
    !
    ! Remerks:
    ! Input velocities are supposed to be from the old timestep
    ! Normal direction of outflow is assumed to be the primary flow direction
    !
    if((sigma.EQ.0.0).OR.((max(abs(gacc_x),abs(gacc_y),abs(gacc_z))).EQ.0.0)) then
      l_cap = 0.0
    else
      l_cap = sqrt(sigma/((rho2-rho1)*max(abs(gacc_x),abs(gacc_y),abs(gacc_z))))    
    endif
    u_cap = sqrt(l_cap*max(abs(gacc_x),abs(gacc_y),abs(gacc_z)))
    !
    cache_g1 = 0.0
    cache_g2 = 0.0
    cache_l1 = 0.0
    cache_l2 = 0.0
    !
    select case(out_dir)
    !###################### Positive Y-direction ###################################
    case(2)
      ! x, z are face centered locations for y velocity!!!
      if (out_forced) then
        do k=IBM_h+1,n(3)
          do j=1,n(2)
            y = (j+ijk_start(2))*dl(2)
            do i=1,n(1)
              h_val    = (2.0/(1.0+exp(4.0*(l(2)-y)/out_len)))
              cache_g1 = cache_g1 + v(i,j,k)*h_val*dl(1)*dl(2)*dl(3)
              cache_g2 = cache_g2 +          h_val*dl(1)*dl(2)*dl(3)
              cache_l1 = cache_l1 + v(i,j,k)*h_val*dl(1)*dl(2)*dl(3)
              cache_l2 = cache_l2 +          h_val*dl(1)*dl(2)*dl(3)
            enddo
          enddo
        enddo
      else
        do k=IBM_h+1,n(3)
          do j=1,n(2)
            y = (j+ijk_start(2))*dl(2)
            do i=1,n(1)
              h_val    = (2.0/(1.0+exp(4.0*(l(2)-y)/out_len)))
              phiy     = 0.5*(psi(i,j,k)+psi(i,j+1,k))
              cache_g1 = cache_g1 + v(i,j,k)*phiy*h_val*dl(1)*dl(2)*dl(3)
              cache_g2 = cache_g2 +          phiy*h_val*dl(1)*dl(2)*dl(3)
              cache_l1 = cache_l1 + v(i,j,k)*(1.0-phiy)*h_val*dl(1)*dl(2)*dl(3)
              cache_l2 = cache_l2 +          (1.0-phiy)*h_val*dl(1)*dl(2)*dl(3)
            enddo
          enddo
        enddo
      endif
    !###################### Positive Z-direction ###################################
    case(3)
      ! x,y are face centered locations for z velocity!!!
      if (out_forced) then
        do k=IBM_h+1,n(3) !
          z = (k+ijk_start(3))*dl(3)
          do j=1,n(2)
            y = (j+ijk_start(2)-0.5_rp)*dl(2)
            do i=1,n(1)
              x        = (i+ijk_start(1)-0.5_rp)*dl(1)
              h_val    = (2.0/(1.0+exp(4.0*(l(3)-z)/out_len)))
              cache_g1 = cache_g1 + w(i,j,k)*h_val*dl(1)*dl(2)*dl(3)
              cache_g2 = cache_g2 +          h_val*dl(1)*dl(2)*dl(3)
              cache_l1 = cache_l1 + w(i,j,k)*h_val*dl(1)*dl(2)*dl(3)
              cache_l2 = cache_l2 +          h_val*dl(1)*dl(2)*dl(3)
            enddo
          enddo
        enddo
      else
        do k=IBM_h+1,n(3) !
          z = (k+ijk_start(3))*dl(3)
          do j=1,n(2)
            y = (j+ijk_start(2)-0.5_rp)*dl(2)
            do i=1,n(1)
              x        = (i+ijk_start(1)-0.5_rp)*dl(1)
              h_val    = (2.0/(1.0+exp(4.0*(l(3)-z)/out_len)))
              phiz     = 0.5*(psi(i,j,k)+psi(i,j,k+1))
              cache_g1 = cache_g1 + w(i,j,k)*phiz*h_val*dl(1)*dl(2)*dl(3)
              cache_g2 = cache_g2 +          phiz*h_val*dl(1)*dl(2)*dl(3)
              cache_l1 = cache_l1 + w(i,j,k)*(1.0-phiz)*h_val*dl(1)*dl(2)*dl(3)
              cache_l2 = cache_l2 +          (1.0-phiz)*h_val*dl(1)*dl(2)*dl(3)
            enddo
          enddo
        enddo
      endif
    !###################### Negative Z-direction ###################################
    case(-3)
      ! x,y are face centered locations for z velocity!!!
      if (out_forced) then
        do k=IBM_h+1,n(3) !
          z = (k+ijk_start(3))*dl(3)
          do j=1,n(2)
            y = (j+ijk_start(2)-0.5_rp)*dl(2)
            do i=1,n(1)
              x        = (i+ijk_start(1)-0.5_rp)*dl(1)
              h_val    = (2.0/(1.0+exp(4.0*(z)/out_len)))
              cache_g1 = cache_g1 - w(i,j,k)*h_val*dl(1)*dl(2)*dl(3)
              cache_g2 = cache_g2 -          h_val*dl(1)*dl(2)*dl(3)
              cache_l1 = cache_l1 - w(i,j,k)*h_val*dl(1)*dl(2)*dl(3)
              cache_l2 = cache_l2 -          h_val*dl(1)*dl(2)*dl(3)
            enddo
          enddo
        enddo
      else
        do k=IBM_h+1,n(3) !
          z = (k+ijk_start(3))*dl(3)
          do j=1,n(2)
            y = (j+ijk_start(2)-0.5_rp)*dl(2)
            do i=1,n(1)
              x        = (i+ijk_start(1)-0.5_rp)*dl(1)
              h_val    = (2.0/(1.0+exp(4.0*(z)/out_len)))
              phiz     = 0.5*(psi(i,j,k)+psi(i,j,k+1))
              cache_g1 = cache_g1 - w(i,j,k)*phiz*h_val*dl(1)*dl(2)*dl(3)
              cache_g2 = cache_g2 -          phiz*h_val*dl(1)*dl(2)*dl(3)
              cache_l1 = cache_l1 - w(i,j,k)*(1.0-phiz)*h_val*dl(1)*dl(2)*dl(3)
              cache_l2 = cache_l2 -          (1.0-phiz)*h_val*dl(1)*dl(2)*dl(3)
            enddo
          enddo
        enddo
      endif
    case(0)
      ! Do nothing
    case default
      call flutas_error('ERROR: This outflow direction is not implemented ... yet!')
    end select
    !
    if(out_dir.NE.0) then! If outflow treatment is activated
      call mpi_allreduce(MPI_IN_PLACE,cache_g1 ,1,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      call mpi_allreduce(MPI_IN_PLACE,cache_g2 ,1,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      call mpi_allreduce(MPI_IN_PLACE,cache_l1 ,1,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      call mpi_allreduce(MPI_IN_PLACE,cache_l2 ,1,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      u_gas    = cache_g1/cache_g2
      if (cache_g2.EQ.0.0) u_gas = 0.0 !avoid NaN
      u_liquid = cache_l1/cache_l2
      u_liquid = max(u_cap,u_liquid)
      u_gas = max(u_cap,u_gas)
    endif
    !
  end subroutine cmpt_avg_out_vel
  !
  !_______________________________________________________________________________
  !
  subroutine corr_u_liquid(n,dl,l,nh_p,nh_u,psi,f_t12,src_ph,u,v,w)
    ! Different from A. Dhruv, J. of Comp. Phys., 2024, the liquid velocity is corrected to 
    ! avoid pressure oscillation by ensuring that the outflow velocities do not try to violate
    ! incompressability constraint.
    !
    ! Assumption: all liquid flows out through the outflow
    use mod_param     , only: sigma,rho1,rho2,gacc_x,gacc_y,gacc_z,IBM_h,bcvel,is_outflow
    use mod_common_mpi, only: MPI_PROC_NULL,back
    !
    implicit none
    integer         , intent(in ), dimension(3)                       :: n
    real(rp)        , intent(in ), dimension(3)                       :: dl,l
    real(rp)        , intent(in )                                     :: f_t12
    integer         , intent(in )                                     :: nh_p,nh_u
    real(rp)        , intent(in ), dimension(1-nh_p:,1-nh_p:,1-nh_p:) :: psi,src_ph
    real(rp)        , intent(in ), dimension(1-nh_u:,1-nh_u:,1-nh_u:) :: u,v,w
    real(rp)                                                          :: integral1,integral2,integral3
    real(rp)                                                          :: vol_change,u_liquid_corr
    integer                                                           :: i,j,k
    !
    integral1 = 0.0
    integral2 = 0.0
    integral3 = 0.0
    !
    vol_change = 0.0
    ! Calculate the volume change due to phase change
    ! Only used for method (1). Not used for (2) to avoid dynamic changes of the outflow velocity (more restrictive)
    if((out_dir.NE.0).AND.(.NOT.out_forced)) then
      do k=IBM_h+1,n(3) !
        do j=1,n(2)
          do i=1,n(1)
            vol_change = vol_change + (1.d0/rho1 - 1.d0/rho2) &
                                    *dl(3)*dl(1)*dl(2)*src_ph(i,j,k)
          enddo
        enddo
      enddo
      !
      call mpi_allreduce(MPI_IN_PLACE,vol_change,1,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
    endif
    !
    select case(out_dir)
      ! Calculate:
      ! 1. Gas volume flux leaving the domian in the outflow boundary: integral1
      ! 2. Area of outflow affected by the liquid phase: integral2
      ! 3. Total area of outflow boundary: integral3
      ! 
      case(2)
        !###################### Positive Y-direction ###################################
        ! x,y,z are face centered locations for y velocity!!!
        if(is_outflow(1,2).AND.(back .eq.MPI_PROC_NULL)) then
        do k=IBM_h+1,n(3)
          j=n(2)
          do i=1,n(1)
            !
            integral1 = integral1 + dl(1)*dl(3)*(      0.5*(psi(i,j,k)+psi(i,j+1,k)))*v(i,j,k)
            integral2 = integral2 + dl(1)*dl(3)*(1.0 - 0.5*(psi(i,j,k)+psi(i,j+1,k)))
            integral3 = integral3 + dl(1)*dl(3)
            !
          enddo
        enddo
        endif
      case(3)
        !###################### Positive Z-direction ###################################
        ! x,y,z are face centered locations for z velocity!!!
        k=n(3)
        do j=1,n(2)
          do i=1,n(1)
            !
            integral1 = integral1 + dl(1)*dl(2)*(      0.5*(psi(i,j,k)+psi(i,j,k+1)))*w(i,j,k)
            integral2 = integral2 + dl(1)*dl(2)*(1.0 - 0.5*(psi(i,j,k)+psi(i,j,k+1)))
            integral3 = integral3 + dl(1)*dl(2)
            !
          enddo
        enddo
      case(-3)
        !###################### Negative Z-direction ###################################
        ! x,y,z are face centered locations for z velocity!!!
        k=0
        do j=1,n(2)
          do i=1,n(1)
            !
            integral1 = integral1 - dl(1)*dl(2)*(      0.5*(psi(i,j,k)+psi(i,j,k+1)))*w(i,j,k)
            integral2 = integral2 + dl(1)*dl(2)*(1.0 - 0.5*(psi(i,j,k)+psi(i,j,k+1)))
            integral3 = integral3 + dl(1)*dl(2)
            !
          enddo
        enddo
      case(0)
        ! Do nothing
      case default
        call flutas_error('ERROR: This outflow direction is not implemented ... yet!')
    end select
    !
    if(out_dir.NE.0) then
      if(.NOT.out_forced) then
        ! integral3 not needed for method (3)
        call mpi_allreduce(MPI_IN_PLACE,integral1 ,1,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
        call mpi_allreduce(MPI_IN_PLACE,integral2 ,1,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
        ! required liquid velocity to ensure not violating the incompressability constraint:
        !                 ,-> Safety margin: 2%
        u_liquid_corr = 1.02*abs((compute_influxes(n,nh_u,u,v,w,f_t12,dl) + vol_change - integral1)/integral2)
        !
      else
        call mpi_allreduce(MPI_IN_PLACE,integral3 ,1,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
        !                ,-> Safety margin: 40% to account for phase change (see p.52-53 in [2])
        u_liquid_corr = 1.4*abs(compute_influxes(n,nh_u,u,v,w,f_t12,dl)/integral3)
        !
      endif
      ! Now update:
      u_liquid_corr = max(u_cap,u_liquid_corr,u_liquid) ! Take max to be on the safe side
      if(myid.eq.0) print*, 'Change liquid outflow vel by:',   100*(u_liquid_corr-u_liquid)/u_liquid , " %"
      u_liquid = u_liquid_corr
      u_gas = u_liquid_corr
      !
    endif
    !
  end subroutine corr_u_liquid
  !
  !_______________________________________________________________________________
  !
  subroutine add_outflow_force(n,dl,l,nh_p,nh_u,psi,rho,u,v,w,uold,vold,wold,f_t12)
    ! Adds a forcing term in the outflow region so allow for a stable exit of bubbles through the outflow boundary
    ! For Method (1): see A. Dhruv, J. of Comp. Phys., 2024
    ! For Method (2): see [2] eq. 4.68 and eq. 4.69
    use mod_param     , only: sigma,rho1,rho2,gacc_x,gacc_y,gacc_z,IBM_h
    implicit none
    integer         , intent(in   ), dimension(3)                       :: n
    real(rp)        , intent(in   )                                     :: f_t12
    real(rp)        , intent(in   ), dimension(3)                       :: dl,l
    integer         , intent(in   )                                     :: nh_p,nh_u
    real(rp)        , intent(in   ), dimension(1-nh_p:,1-nh_p:,1-nh_p:) :: psi,rho
    real(rp)        , intent(inout), dimension(1-nh_u:,1-nh_u:,1-nh_u:) :: u,v,w,uold,vold,wold
    real(rp)                                                            :: x,y,z,h_val,u_e
    real(rp)                                                            :: phiz,phiy,phix
    real(rp)                                                            :: rhoz,rhoy,rhox
    real(rp)                                                            :: dudt,dvdt,dwdt,switch
    integer                                                             :: i,j,k
    !
    ! Remerks:
    ! Input velocities are supposed to be from the old timestep
    ! Normal direction of outflow is assumed to be the primary flow direction
    !
    !
    select case(out_dir)
      !###################### Positive y-direction ###################################
      case(2)
        do k=IBM_h+1,n(3) !
          z = (k+ijk_start(3)-0.5_rp)*dl(3)
          do j=1,n(2)
            y = (j+ijk_start(2))*dl(2)
            do i=1,n(1)
              x        = (i+ijk_start(1)-0.5_rp)*dl(1)
              phix     = 0.5*(psi(i,j,k)+psi(i+1,j,k))
              phiy     = 0.5*(psi(i,j,k)+psi(i,j+1,k))
              phiz     = 0.5*(psi(i,j,k)+psi(i,j,k+1))
              rhox     = 0.5*(rho(i,j,k)+rho(i+1,j,k))
              rhoy     = 0.5*(rho(i,j,k)+rho(i,j+1,k))
              rhoz     = 0.5*(rho(i,j,k)+rho(i,j,k+1))
              !_______________________________________
              ! Scaling in z
              h_val    = (2.0/(1.0+exp(4.0*(l(2)-y+0.5*dl(2))/out_len)))          
              u_e      = phiz*u_gas + (1.0-phiz)*u_liquid
              ! scaling part:
              dwdt = h_val*(min(1.0,abs(u_e/wold(i,j,k)))-1.0_rp)*wold(i,j,k)/f_t12
              if(out_forced) dwdt = - 0.2*h_val*wold(i,j,k)/f_t12
              ! convection 
              dwdt = dwdt - abs(u_e)*h_val*(wold(i,j,k)-wold(i,j-1,k))/dl(2)   ! first order upwind scheme was important for stability!
              !_______________________________________
              ! Scaling in y
              h_val    = (2.0/(1.0+exp(4.0*(l(2)-y)/out_len)))
              u_e      = phiy*u_gas + (1.0-phiy)*u_liquid
              ! scaling part:
              dvdt = h_val*(min(1.d0,abs(u_e/vold(i,j,k)))-1.0_rp)*vold(i,j,k)/f_t12
              if((vold(i,j,k).le.(0.1*u_e)).AND.out_forced) dvdt = - h_val*(vold(i,j,k)-(0.1*u_e))/f_t12
              ! convection 
              dvdt = dvdt - abs(u_e)*h_val*(vold(i,j,k)-vold(i,j-1,k))/dl(2) ! first order upwind scheme was important for stability!
              !_______________________________________
              ! Scaling in x
              h_val    = (2.0/(1.0+exp(4.0*(l(2)-y+0.5*dl(2))/out_len))) 
#if defined(_TWOD)
              u_e      = 0.0
#else
              u_e      = phix*u_gas + (1.0-phix)*u_liquid
#endif
              ! scaling part:
              dudt = h_val*(min(1.0,abs(u_e/uold(i,j,k)))-1.0_rp)*uold(i,j,k)/f_t12
              if(out_forced) dudt = - 0.2*h_val*uold(i,j,k)/f_t12
              ! convection 
              dudt = dudt - abs(u_e)*h_val*(uold(i,j,k)-uold(i,j-1,k))/dl(2) ! first order upwind scheme was important for stability!
              !_______________________________________
              u(i,j,k) = u(i,j,k) + rhox*dudt*f_t12
              v(i,j,k) = v(i,j,k) + rhoy*dvdt*f_t12
              w(i,j,k) = w(i,j,k) + rhoz*dwdt*f_t12
            enddo
          enddo
        enddo
      !###################### Positive Z-direction ###################################
      case(3)
        if(out_forced.AND.(myid.eq.0)) print*, 'WARNING: Forced convection is not implemented for this direction!'
        do k=IBM_h+1,n(3) !
          z = (k+ijk_start(3))*dl(3)
          do j=1,n(2)
            y = (j+ijk_start(2)-0.5_rp)*dl(2)
            do i=1,n(1)
              x        = (i+ijk_start(1)-0.5_rp)*dl(1)
              phix     = 0.5*(psi(i,j,k)+psi(i+1,j,k))
              phiy     = 0.5*(psi(i,j,k)+psi(i,j+1,k))
              phiz     = 0.5*(psi(i,j,k)+psi(i,j,k+1))
              rhox     = 0.5*(rho(i,j,k)+rho(i+1,j,k))
              rhoy     = 0.5*(rho(i,j,k)+rho(i,j+1,k))
              rhoz     = 0.5*(rho(i,j,k)+rho(i,j,k+1))
              !_______________________________________
              ! Scaling in z
              h_val    = (2.0/(1.0+exp(4.0*(l(3)-z)/out_len)))          
              u_e      = phiz*u_gas + (1.0-phiz)*u_liquid
              ! scaling part:
              dwdt = h_val*(min(1.d0,abs(u_e/wold(i,j,k)))-1.d0)*wold(i,j,k)/f_t12
              ! convection 
              dwdt = dwdt - abs(u_e)*h_val*(wold(i,j,k)-wold(i,j,k-1))/dl(3)   ! first order upwind scheme was important for stability!
              !_______________________________________
              ! Scaling in y
              h_val    = (2.0/(1.0+exp(4.0*(l(3)-z+0.5*dl(3))/out_len)))          
              u_e      = phiy*u_gas + (1.0-phiy)*u_liquid
              ! scaling part:
              dvdt = h_val*(min(1.d0,abs(u_e/vold(i,j,k)))-1.d0)*vold(i,j,k)/f_t12
              ! convection 
              dvdt = dvdt - abs(u_e)*h_val*(vold(i,j,k)-vold(i,j,k-1))/dl(3) ! first order upwind scheme was important for stability!
              !_______________________________________
              ! Scaling in x
              h_val    = (2.0/(1.0+exp(4.0*(l(3)-z+0.5*dl(3))/out_len)))
#if defined(_TWOD)
              u_e      = 0.0
#else
              u_e      = phix*u_gas + (1.0-phix)*u_liquid
#endif
              ! scaling part:
              dudt = h_val*(min(1.d0,abs(u_e/uold(i,j,k)))-1.d0)*uold(i,j,k)/f_t12
              ! convection 
              dudt = dudt - abs(u_e)*h_val*(uold(i,j,k)-uold(i,j,k-1))/dl(3) ! first order upwind scheme was important for stability!
              !_______________________________________
              u(i,j,k) = u(i,j,k) + rhox*dudt*f_t12
              v(i,j,k) = v(i,j,k) + rhoy*dvdt*f_t12
              w(i,j,k) = w(i,j,k) + rhoz*dwdt*f_t12
            enddo
          enddo
        enddo
      !###################### Negative Z-direction ###################################
      case(-3)
        if(out_forced.AND.(myid.eq.0)) print*, 'WARNING: Forced convection is not implemented for this direction!'
        do k=IBM_h+1,n(3) !
          z = (k+ijk_start(3))*dl(3)
          do j=1,n(2)
            y = (j+ijk_start(2)-0.5_rp)*dl(2)
            do i=1,n(1)
              x        = (i+ijk_start(1)-0.5_rp)*dl(1)
              phix     = 0.5*(psi(i,j,k)+psi(i+1,j,k))
              phiy     = 0.5*(psi(i,j,k)+psi(i,j+1,k))
              phiz     = 0.5*(psi(i,j,k)+psi(i,j,k+1))
              rhox     = 0.5*(rho(i,j,k)+rho(i+1,j,k))
              rhoy     = 0.5*(rho(i,j,k)+rho(i,j+1,k))
              rhoz     = 0.5*(rho(i,j,k)+rho(i,j,k+1))
              !_______________________________________
              ! Scaling in z
              h_val    = (2.0/(1.0+exp(4.0*(z)/out_len)))          
              u_e      = phiz*u_gas + (1.0-phiz)*u_liquid
              ! scaling part:
              dwdt = h_val*(min(1.d0,abs(u_e/wold(i,j,k)))-1.d0)*wold(i,j,k)/f_t12
              ! convection 
              dwdt = dwdt + abs(u_e)*h_val*(wold(i,j,k+1)-wold(i,j,k))/dl(3)   ! first order upwind scheme was important for stability!!
              !           '-> Note the sign is switched here because of negative direction
              !_______________________________________
              ! Scaling in y
              h_val    = (2.d0/(1.d0+exp(4.0*(z-0.5*dl(3))/out_len)))          
              u_e      = phiy*u_gas + (1.d0-phiy)*u_liquid
              ! scaling part:
              dvdt = h_val*(min(1.d0,abs(u_e/vold(i,j,k)))-1.d0)*vold(i,j,k)/f_t12
              ! convection 
              dvdt = dvdt + abs(u_e)*h_val*(vold(i,j,k+1)-vold(i,j,k))/dl(3)   ! first order upwind scheme was important for stability!!
              !           '-> Note the sign is switched here because of negative direction
              !_______________________________________
              ! Scaling in x
              h_val    = (2.d0/(1.d0+exp(4.d0*(z-0.d0*dl(3))/out_len)))          
#if defined(_TWOD)
              u_e      = 0.0
#else
              u_e      = phix*u_gas + (1.0-phix)*u_liquid
#endif
              ! scaling part:
              dudt = h_val*(min(1.d0,abs(u_e/uold(i,j,k)))-1.d0)*uold(i,j,k)/f_t12
              ! convection 
              dudt = dudt + abs(u_e)*h_val*(uold(i,j,k+1)-uold(i,j,k))/dl(3)   ! first order upwind scheme was important for stability!!
              !           '-> Note the sign is switched here because of negative direction
              !_______________________________________
              u(i,j,k) = u(i,j,k) + rhox*dudt*f_t12
              v(i,j,k) = v(i,j,k) + rhoy*dvdt*f_t12
              w(i,j,k) = w(i,j,k) + rhoz*dwdt*f_t12
            enddo
          enddo
        enddo
      case(0)
        ! Do nothing
      case default
        call flutas_error('ERROR: This outflow direction is not implemented ... yet!')
    end select
  end subroutine add_outflow_force
  !
  !_______________________________________________________________________________
  !
  subroutine outflow_psi(n,dl,nh_v,phiold,u,v,w,dt,bound_o,tmp,back_val)
    ! Calculates the boundary value at the outflow of the input field "phiold" by solving a simple convection equation
    use mod_common_mpi, only: MPI_PROC_NULL,left,right,front,back,top,bottom
    !                             
    ! Calculates the convected boundary condition of the input
    !
    implicit none
    !
    integer , intent(in ), dimension(3)                       :: n
    real(rp), intent(in ), dimension(3)                       :: dl
    integer , intent(in )                                     :: nh_v
    real(rp), intent(in ), dimension(1-nh_v:,1-nh_v:,1-nh_v:) :: phiold,tmp
    real(rp), intent(in ), dimension(1-nh_v:,1-nh_v:,1-nh_v:) :: u,v,w
    real(rp), intent(in )                                     :: dt,back_val
    real(rp), intent(out),dimension(1:,1:)                    :: bound_o
    real(rp)                                                  :: psi_wall
    real(rp)                                                  :: phi_gradx,phi_grady,phi_gradz
    real(rp)                                                  :: u_loc,v_loc,w_loc,aux,t_loc
    integer                                                   :: i,j,k,kpp,kp,ip,im,jp,jm,km
    !
    select case(out_dir)
      !###################### Positive Y-direction ###################################
      case(2)
        !
        bound_o(:,:) = 0.0
        ! we have to make sure that only the rank that is at the outflow is doing the following:
        if(back .eq.MPI_PROC_NULL) then 
          !
          j   =   n(2)
          jp  = j + 1
          jm  = j - 1    
          !
          do k=1,n(3)
            do i=1,n(1)
              !
              ip  = i + 1
              im  = i - 1
              km  = k - 1
              kp  = k + 1
              ! Transport velocities assumption du/dn = 0:
              u_loc             = 0.5*(u(i ,j ,k )+u(im,j ,k ))
              w_loc             = 0.5*(w(i ,j ,k )+w(i ,j ,km))
              v_loc             =     (v(i ,j ,k )            )
              ! Setting backflow value:
              if((v_loc.LT.0.0)) then
                aux             = back_val
              else
                aux             = phiold(i,jp,k)
              endif
              ! Old boundary value:
              psi_wall          = 0.5*((phiold(i ,jp,k )+phiold(i ,j ,k ))            )
              ! Gradients:
              phi_grady         =     ((aux             -phiold(i ,j ,k ))/dl(3)      )
              phi_gradz         = 0.5*((phiold(i ,j ,kp)-phiold(i ,j ,km))/dl(2)/(2.0) &
                                      +(phiold(i ,jp,kp)-phiold(i ,jp,km))/dl(2)/(2.0))
              phi_gradx         = 0.5*((phiold(ip,j ,k )-phiold(im,j ,k ))/dl(1)/(2.0) &
                                      +(phiold(ip,jp,k )-phiold(im,jp,k ))/dl(1)/(2.0))
              !  Boundary value for the next timestep by solving the transport equation:
              psi_wall          = max(0.d0,psi_wall - dt*(u_loc*phi_gradx &
                                                         +v_loc*phi_grady &
                                                         +w_loc*phi_gradz))
              !                                  
              bound_o(i,k)    =  psi_wall
              !
            enddo
          enddo
        endif
      case(3)
        !###################### Positive Z-direction ###################################
        ! if the outflow is in positiv z direction, all ranks are outflow ranks...
        bound_o(:,:) = 0.0
        k   =   n(3)
        km  = k - 1
        kp  = k + 1
        !
        do j=1,n(2)
          do i=1,n(1)
            !
            ip  = i + 1
            im  = i - 1
            jp  = j + 1
            jm  = j - 1
            ! Transport velocities assumption du/dn = 0:
            u_loc             = 0.5*(u(i ,j ,k )+u(im,j ,k ))
            v_loc             = 0.5*(v(i ,j ,k )+v(i ,jm,k ))
            w_loc             =     (w(i ,j ,k )            )   
            ! Setting backflow value:
            if((w_loc.LT.0.0)) then
              aux             = back_val
            else
              aux             = phiold(i,j,kp)
            endif
            ! Old boundary value:
            psi_wall          = 0.5*((phiold(i ,j ,k )+phiold(i ,j ,kp))            )
            ! Gradients:
            phi_gradz         =     ((aux             -phiold(i ,j ,k ))/dl(3)      )
            phi_grady         = 0.5*((phiold(i ,jp,k )-phiold(i ,jm,k ))/dl(2)/(2.0) &
                                    +(phiold(i ,jp,kp)-phiold(i ,jm,kp))/dl(2)/(2.0))
            phi_gradx         = 0.5*((phiold(ip,j ,k )-phiold(im,j ,k ))/dl(1)/(2.0) &
                                    +(phiold(ip,j ,kp)-phiold(im,j ,kp))/dl(1)/(2.0))
            ! Boundary value for the next timestep by solving the transport equation:
            psi_wall          = psi_wall - dt*(u_loc*phi_gradx &
                                              +v_loc*phi_grady &
                                              +w_loc*phi_gradz)
            !                                  
            bound_o(i,j)    =  psi_wall
            !
          enddo
        enddo
      case(-3)
        !###################### Negative Z-direction ###################################
        ! if the outflow is in negative z direction, all ranks are outflow ranks...
        bound_o(:,:) = 0.0
        k   =   1
        km  = k - 1
        kp  = k + 1
        !
        do j=1,n(2)
          do i=1,n(1)
            !
            ip  = i + 1
            im  = i - 1
            jp  = j + 1
            jm  = j - 1
            ! Transport velocities assumption du/dn = 0:
            u_loc             = 0.5*(u(i ,j ,k )+u(im,j ,k ))
            v_loc             = 0.5*(v(i ,j ,k )+v(i ,jm,k ))
            w_loc             =     (w(i ,j ,k-1)           )   
            ! Setting backflow value:
            if((w_loc.GT.0.0)) then
              aux             = back_val
            else
              aux             = phiold(i,j,km)
            endif
            ! Old boundary value:
            psi_wall          = 0.5*((phiold(i ,j ,k )+phiold(i ,j ,km))            )
            ! Gradients:
            phi_gradz         =     ((phiold(i ,j ,k )-aux             )/dl(3)      )
            phi_grady         = 0.5*((phiold(i ,jp,k )-phiold(i ,jm,k ))/dl(2)/(2.0) &
                                    +(phiold(i ,jp,km)-phiold(i ,jm,km))/dl(2)/(2.0))
            phi_gradx         = 0.5*((phiold(ip,j ,k )-phiold(im,j ,k ))/dl(1)/(2.0) &
                                    +(phiold(ip,j ,km)-phiold(im,j ,km))/dl(1)/(2.0))
            ! Boundary value for the next timestep by solving the transport equation:
            psi_wall          = psi_wall - dt*(u_loc*phi_gradx &
                                              +v_loc*phi_grady &
                                              +w_loc*phi_gradz)
            !                                  
            bound_o(i,j)    =  psi_wall
            !
          enddo
        enddo
      case(0)
        ! Do nothing
      case default
        call flutas_error('ERROR: This outflow direction is not implemented ... yet!')
    end select
    !
    return
  end subroutine outflow_psi
  !
  !_______________________________________________________________________________
  !
  subroutine boundsb_arrayo(n,nh_p,halo,s,bound_o)
    ! applies Dirichlet BC at the outflow that was calculated from outflow_psi
    use mod_common_mpi, only: MPI_PROC_NULL,left,right,front,back,top,bottom
    !   
    implicit none
    !   
    integer , intent(in   ), dimension(3)                       :: n 
    integer , intent(in   )                                     :: nh_p
    real(rp), intent(inout), dimension(1-nh_p:,1-nh_p:,1-nh_p:) :: s
    integer , intent(in   ), dimension(3)                       :: halo
    real(rp), intent(in   ), dimension(1:,1:)                   :: bound_o
    !   
    integer :: q,i,j,k
    integer :: ind1,ind2
    !
    !@cuf attributes(managed) :: s  
#if defined(_DECOMP_X)
    ind1 = 2 
    ind2 = 3 
#elif _DECOMP_Y
    ind1 = 1 
    ind2 = 3 
#else
    ind1 = 1 
    ind2 = 2 
#endif
    !
    select case(out_dir)
      !###################### Positive Y-direction ###################################
      case(2)
        ! we have to make sure that only the rank that is at the outflow is doing the following:
        if(back .eq.MPI_PROC_NULL) then 
          !
          do q=0,nh_p-1
            do k=1,n(3)
              do i=1,n(1)
                s(i,n(2)+1+q,k) = 2.0*bound_o(i,k) - s(i,n(2)-q,k)
              enddo
            enddo
          enddo
        endif
        !###################### Positive Z-direction ###################################
      case(3)
        ! if the outflow is in positiv z-direction, all ranks are outflow ranks...
        do q=0,nh_p-1
          do j=1,n(2)
            do i=1,n(1)
              s(i,j,n(3)+1+q) = 2.0*bound_o(i,j) - s(i,j,n(3)-q)
            enddo
          enddo
        enddo
        !###################### Negative Z-direction ###################################
      case(-3)
        ! if the outflow is in negative z-direction, all ranks are outflow ranks...
        do q=0,nh_p-1
          do j=1,n(2)
            do i=1,n(1)
              s(i,j,-q) = 2.0*bound_o(i,j) - s(i,j,1+q)
            enddo
          enddo
        enddo
      case(0)
        ! Do nothing
      case default
        call flutas_error('ERROR: This outflow direction is not implemented ... yet!')
    end select
    !
    call updthalo(n(1),n(2),n(3),nh_p,halo(ind1),ind1,s)
    call updthalo(n(1),n(2),n(3),nh_p,halo(ind2),ind2,s)
    !   
    return
  end subroutine boundsb_arrayo
  !
  !_______________________________________________________________________________
  !
  real(rp) function compute_influxes(n,nh_u,u,v,w,dt,dl)
    ! returns the current volume inflow through all faces that are NOT and outflow boundary
    use mod_param,      only: is_outflow
    !   
    ! imposes Dirichlet BC at the outflow
    !   
    implicit none
    !   
    integer , intent(in   ), dimension(3)                       :: n
    integer , intent(in   )                                     :: nh_u
    real(rp), intent(in   )                                     :: dt
    real(rp), intent(in   ), dimension(3)                       :: dl
    real(rp), intent(in   ), dimension(1-nh_u:,1-nh_u:,1-nh_u:) :: u,v,w
    integer                                                     :: i,j,k,q,sgn,idir,ioutflowdir
    real(rp)                                                    :: change,normal
    ! sgn is the inward facing vector
    change = 0.0
    do q = 1,3
      do idir = 0,1
        if(.NOT.is_outflow(idir,q)) then
          if(idir.eq.0) sgn = -1
          if(idir.eq.1) sgn = +1
          ioutflowdir = q*sgn
          normal = -1.0*sgn
          change = change + normal*surface_integral(n(1),n(2),n(3),ioutflowdir,nh_u,dl(1),dl(2),dl(3),u,v,w)
        endif
      enddo
    enddo
    !
    call mpi_allreduce(MPI_IN_PLACE,change,1,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
    !
    compute_influxes = change
    !
  end function compute_influxes
  !
  !_______________________________________________________________________________
  !
  real(rp) function surface_integral(nx,ny,nz,idir,nh_u,dx,dy,dz,u,v,w)
    ! calculate flux through faces with facecentered velocities
    use mod_common_mpi, only: MPI_PROC_NULL,left,right,front,back,top,bottom
    !
    implicit none
    !
    integer , intent(in   )                                     :: nx,ny,nz
    integer , intent(in   )                                     :: idir
    integer , intent(in   )                                     :: nh_u
    real(rp), intent(in   )                                     :: dx,dy,dz
    real(rp), intent(in   ), dimension(1-nh_u:,1-nh_u:,1-nh_u:) :: u,v,w
    integer                                                     :: i,j,k
    real(rp)                                                    :: integral
    !
    !
    integral = 0.0
    select case(idir)
      case(1) ! x direction, right
        if(right.eq.MPI_PROC_NULL) then
          i = nx
          do k=1,nz
            do j=1,ny
              integral = integral + u(i,j,k)*dy*dz
            enddo
          enddo
        endif
      case(2) ! y direction, back
        if(back.eq.MPI_PROC_NULL) then
          j = ny
          do k=1,nz
            do i=1,nx
              integral = integral + v(i,j,k)*dx*dz
            enddo
          enddo 
        endif
      case(3) ! z direction, top
        if(top.eq.MPI_PROC_NULL) then
          k = nz
          do j=1,ny
            do i=1,nx
              integral = integral + w(i,j,k)*dx*dy
            enddo
          enddo 
        endif
      case(-1) ! x direction, left
        if(left.eq.MPI_PROC_NULL) then
          i = 0
          do k=1,nz
            do j=1,ny
              integral = integral + u(i,j,k)*dy*dz
            enddo
          enddo 
        endif
      case(-2) ! y direction, front
        if(front.eq.MPI_PROC_NULL) then
          j = 0
          do k=1,nz
            do i=1,nx
              integral = integral + v(i,j,k)*dx*dz
            enddo
          enddo 
        endif
      case(-3) ! z direction, bottom
        if(bottom.eq.MPI_PROC_NULL) then
          k = 0
          do j=1,ny
            do i=1,nx
              integral = integral + w(i,j,k)*dx*dy
            enddo
          enddo
        endif
    end select
    !
    surface_integral = integral
    !
    return
  end function surface_integral
  !
  !_______________________________________________________________________________
  !
end module mod_outflow