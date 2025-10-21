!
! SPDX-License-Identifier: MIT
!
module mod_source
  !
  use mod_types
  use mod_common_mpi, only: myid,ierr,ijk_start
  use mod_sanity    , only: flutas_error
  !
  implicit none
  !
  private
  public :: bulk_forcing_src,surft_src_IBM
  public :: surft_src,stress_IBM
  public :: grav_tw_src,grav_tw_src_ibm
  public :: grav_sp_src,pres_sp_src
#if defined(_TURB_FORCING)
  public :: forc_src
#endif
  public mar_src
  !
  contains
  !
  !
  
  subroutine mar_src(dzfi,nx,ny,nz,nh_d,nh_u,nh_t,f_t12,dxi,dyi,dzi,dzci,delta_arr,nor,kappa,psi,tmp,rho,u,v,w)
    !
    use mod_param, only: sigma
    !   
    implicit none 
    !   
    integer , intent(in   )                                       :: nx,ny,nz
    integer , intent(in   )                                       :: nh_d,nh_u,nh_t
    real(rp), intent(in   )                                       :: f_t12
    real(rp), intent(in   )                                       :: dxi,dyi,dzi
    real(rp)        , intent(in   ), dimension(1-nh_d:)           :: dzfi
    real(rp), intent(in   ), dimension(1-nh_d:)                   :: dzci
    real(rp), intent(in   ), dimension(0:,0:,0:)                  :: psi,rho
    real(rp), intent(in   ), dimension(1-nh_t:,1-nh_t:,1-nh_t:)   :: tmp
    real(rp), intent(in   ), dimension(0:,0:,0:)                  :: kappa,delta_arr
    real(rp), intent(in   ), dimension(0:,0:,0:,1:)               :: nor
    real(rp), intent(inout), dimension(1-nh_u:,1-nh_u:,1-nh_u:)  :: u,v,w
    real(rp), dimension(-2:nx+3,-2:ny+3,-2:nz+3,1:3)                                      :: mar_force   
    real(rp) :: kappax,kappay,kappaz
    real(rp) :: rhox,rhoy,rhoz
    integer  :: i,j,k
    !@cuf attributes(managed) :: kappa, psi, rho, dzci, u, v, w
    !   
    !$acc kernels
 
   !call marangoni_he((/nx,ny,nz/),(/dxi,dyi,dzi/),1./dzci,1./dzfi,tmp,nh_t,delta_arr,kappa,nor,psi,rho,mar_force)


   do k=1,nz
      do j=1,ny
        do i=1,nx
          u(i,j,k) = u(i,j,k) + f_t12*mar_force(i,j,k,1)
          v(i,j,k) = v(i,j,k) + f_t12*mar_force(i,j,k,2)
          w(i,j,k) = w(i,j,k) + f_t12*mar_force(i,j,k,3)
        enddo
      enddo
    enddo
    !$acc end kernels
    !
    return
  end subroutine mar_src
  !
  !
  subroutine surft_src(nx,ny,nz,nh_d,nh_u,f_t12,dxi,dyi,dzi,dzci,kappa,psi,rho,u,v,w,l,dl)
    !
    !the input velocities are actually rho.u and the output as well
    !
    use mod_param,   only: sigma
    use mod_boiling, only: out_len,out_dir
    !
    implicit none 
    !
    integer , intent(in   )                                     :: nx,ny,nz
    integer , intent(in   )                                     :: nh_d,nh_u
    real(rp), intent(in   )                                     :: f_t12
    real(rp), intent(in   )                                     :: dxi,dyi,dzi
    real(rp), intent(in   ), dimension(3)                       :: l,dl
    real(rp), intent(in   ), dimension(1-nh_d:)                 :: dzci
    real(rp), intent(in   ), dimension(0:,0:,0:)                :: kappa,psi,rho
    real(rp), intent(inout), dimension(1-nh_u:,1-nh_u:,1-nh_u:) :: u,v,w
    !
    real(rp) :: kappax,kappay,kappaz
    real(rp) :: phix,phiy,phiz,z,h_val,switch
    integer  :: i,j,k,ip,jp,kp
    !@cuf attributes(managed) :: kappa, psi, rho, dzci, u, v, w
    !
    !$acc kernels
    switch = 1.0
    if(out_dir.eq.0) switch = 0.0
    do k=1,nz
      do j=1,ny
        do i=1,nx
          !
          ip = i+1   
          jp = j+1   
          kp = k+1  
          !
          phix   = 0.5_rp*(psi(ip,j,k)+psi(i,j,k))
          phiy   = 0.5_rp*(psi(i,jp,k)+psi(i,j,k))
          phiz   = 0.5_rp*(psi(i,j,kp)+psi(i,j,k))
          kappax = 0.5_rp*(kappa(ip,j,k)+kappa(i,j,k))
          kappay = 0.5_rp*(kappa(i,jp,k)+kappa(i,j,k))
          kappaz = 0.5_rp*(kappa(i,j,kp)+kappa(i,j,k))
          !
#if defined(_TWOD)
          phix = 0.0
#endif
          z = (k+ijk_start(3))*dl(3)
          h_val    = switch*(2.0/(1.0+exp(4.0*(l(3)-z+0.5*dl(3))/out_len)))
          u(i,j,k) = u(i,j,k) + (1.0-h_val)*f_t12*(6.d0*sigma*kappax*phix*(1.d0-phix)*(psi(ip,j,k)-psi(i,j,k)))*dxi
          v(i,j,k) = v(i,j,k) + (1.0-h_val)*f_t12*(6.d0*sigma*kappay*phiy*(1.d0-phiy)*(psi(i,jp,k)-psi(i,j,k)))*dyi
          h_val    = switch*(2.0/(1.0+exp(4.0*(l(3)-z)/out_len)))
          w(i,j,k) = w(i,j,k) + (1.0-h_val)*f_t12*(6.d0*sigma*kappaz*phiz*(1.d0-phiz)*(psi(i,j,kp)-psi(i,j,k)))*dzci(k)

          !u(i,j,k) = u(i,j,k) + f_t12*(6.d0*sigma*kappax*(psi(ip,j,k)-psi(i,j,k)))*dxi
          !v(i,j,k) = v(i,j,k) + f_t12*(6.d0*sigma*kappay*(psi(i,jp,k)-psi(i,j,k)))*dyi
          !w(i,j,k) = w(i,j,k) + f_t12*(6.d0*sigma*kappaz*(psi(i,j,kp)-psi(i,j,k)))*dzci(k)
          !
          !
          !u(i,j,k) = u(i,j,k) + f_t12*(sigma*kappax*(psi(ip,j,k)-psi(i,j,k))/rhox)*dxi
          !v(i,j,k) = v(i,j,k) + f_t12*(sigma*kappay*(psi(i,jp,k)-psi(i,j,k))/rhoy)*dyi
          !w(i,j,k) = w(i,j,k) + f_t12*(sigma*kappaz*(psi(i,j,kp)-psi(i,j,k))/rhoz)*dzci(k)
          !
        enddo
      enddo
    enddo
    !$acc end kernels
    !
    return
  end subroutine surft_src

  subroutine surft_src_IBM(nx,ny,nz,nh_d,nh_u,f_t12,dxi,dyi,dzi,dzci,kappa,psi,rho,u,v,w,ind_IBM,l,dl,ls,truncated,limiter)
    !
    !the input velocities are actually rho.u and the output as well
    !
    use mod_param    ,only: sigma,IBM_h,eps_int
    use mod_boiling  ,only: out_len,out_dir,max_kappa,copy_fields
    !
    implicit none 
    !
    integer , intent(in   )                                     :: nx,ny,nz
    integer , intent(in   )                                     :: nh_d,nh_u
    real(rp), intent(in   ), dimension(3)                       :: l,dl
    real(rp), intent(in   )                                     :: f_t12,limiter
    real(rp), intent(in   )                                     :: dxi,dyi,dzi
    real(rp), intent(in   ), dimension(1-nh_d:)                 :: dzci
    real(rp), intent(in   ), dimension(0:,0:,0:)                :: kappa,psi,rho,ind_IBM,ls
    real(rp), intent(inout), dimension(1-nh_u:,1-nh_u:,1-nh_u:) :: u,v,w
    logical,  intent(in   )                                     :: truncated
    !
    real(rp) :: kappax,kappay,kappaz,lsx,lsy,lsz,psifx,psify,psifz
    real(rp) :: phix,phiy,phiz,z,h_val,switch,kxm,kxp,kym,kyp,kzm,kzp,y,h_val_f,h_val_c
    integer  :: i,j,k,ip,jp,kp
    real(rp), allocatable, dimension(:,:,:)                     :: psi_tr
    !@cuf attributes(managed) :: kappa, psi, rho, dzci, u, v, w
    !
    !$acc kernels
    allocate(psi_tr(0:nx+1,0:ny+1,0:nz+1))
    psi_tr(:,:,:) = 0.0
    if(truncated) then
      do k=0,nz+1
        do j=0,ny+1
          do i=0,nx+1
            psi_tr(i,j,k) = (1.0/(1.0-limiter))*(min(max(psi(i,j,k),0.5*limiter),1.0-0.5*limiter)-0.5*limiter)
          enddo
        enddo
      enddo
    else
      do k=0,nz+1
        do j=0,ny+1
          do i=0,nx+1
            psi_tr(i,j,k) = psi(i,j,k)
          enddo
        enddo
      enddo
    endif
    !
    select case(out_dir)
    case(3)
      do k=IBM_h+1,nz
        z        = (k+ijk_start(3))*dl(3)
        h_val_f  = (2.0/(1.0+exp(4.0*(l(3)-z          )/out_len)))
        h_val_c  = (2.0/(1.0+exp(4.0*(l(3)-z+0.5*dl(3))/out_len)))
        do j=1,ny
          do i=1,nx
            !
            ip = i+1   
            jp = j+1   
            kp = k+1  
            !
            phix   = 0.5_rp*(psi(ip,j,k)+psi(i,j,k))
            phiy   = 0.5_rp*(psi(i,jp,k)+psi(i,j,k))
            phiz   = 0.5_rp*(psi(i,j,kp)+psi(i,j,k))
            kxp    = kappa(ip,j,k)
            kxm    = kappa(i ,j,k)
            kyp    = kappa(i,jp,k)
            kym    = kappa(i,j ,k)
            kzp    = kappa(i,j,kp)
            kzm    = kappa(i,j,k )
            kappax = 0.5_rp*(kxp+kxm)
            kappay = 0.5_rp*(kyp+kym)
            kappaz = 0.5_rp*(kzp+kzm)
            ! Kappa is constrained as in https://doi.org/10.1023/A:1011178417620
            kappax = max(-max_kappa, min(kappax,max_kappa))
            kappay = max(-max_kappa, min(kappay,max_kappa))
            kappaz = max(-max_kappa, min(kappaz,max_kappa))
            ! 
#if defined(_TWOD)
            kappax = 0.0
#endif
            !
            u(i,j,k) = u(i,j,k) + f_t12*(1.0-h_val_c)*(sigma*kappax*(psi_tr(ip,j,k)-psi_tr(i,j,k)))*dxi   !*6.0*phix*(1.d0-phix)  
            v(i,j,k) = v(i,j,k) + f_t12*(1.0-h_val_c)*(sigma*kappay*(psi_tr(i,jp,k)-psi_tr(i,j,k)))*dyi   !*6.0*phiy*(1.d0-phiy) 
            w(i,j,k) = w(i,j,k) + f_t12*(1.0-h_val_f)*(sigma*kappaz*(psi_tr(i,j,kp)-psi_tr(i,j,k)))*dzci(k)   !*6.0*phiz*(1.d0-phiz)
            !
          enddo
        enddo
      enddo
      !$acc end kernels
    case(2)
      do k=IBM_h+1,nz
        do j=1,ny
          y        = (j+ijk_start(2))*dl(2)
          h_val_f  = (2.0/(1.0+exp(4.0*(l(2)-y          )/out_len)))
          h_val_c  = (2.0/(1.0+exp(4.0*(l(2)-y+0.5*dl(2))/out_len)))
          do i=1,nx
            !
            ip = i+1   
            jp = j+1   
            kp = k+1  
            !
            phix   = 0.5_rp*(psi(ip,j,k)+psi(i,j,k))
            phiy   = 0.5_rp*(psi(i,jp,k)+psi(i,j,k))
            phiz   = 0.5_rp*(psi(i,j,kp)+psi(i,j,k))
            kxp    = kappa(ip,j,k)
            kxm    = kappa(i ,j,k)
            kyp    = kappa(i,jp,k)
            kym    = kappa(i,j ,k)
            kzp    = kappa(i,j,kp)
            kzm    = kappa(i,j,k )
            kappax = 0.5_rp*(kxp+kxm)
            kappay = 0.5_rp*(kyp+kym)
            kappaz = 0.5_rp*(kzp+kzm)
            ! Kappa is constrained as in https://doi.org/10.1023/A:1011178417620
            kappax = max(-max_kappa, min(kappax,max_kappa))
            kappay = max(-max_kappa, min(kappay,max_kappa))
            kappaz = max(-max_kappa, min(kappaz,max_kappa))
            ! 
#if defined(_TWOD)
            kappax = 0.0
#endif
            !
            u(i,j,k) = u(i,j,k) + f_t12*(1.0-h_val_c)*(sigma*kappax*(psi_tr(ip,j,k)-psi_tr(i,j,k)))*dxi   !*6.0*phix*(1.d0-phix)  
            v(i,j,k) = v(i,j,k) + f_t12*(1.0-h_val_f)*(sigma*kappay*(psi_tr(i,jp,k)-psi_tr(i,j,k)))*dyi   !*6.0*phiy*(1.d0-phiy) 
            w(i,j,k) = w(i,j,k) + f_t12*(1.0-h_val_c)*(sigma*kappaz*(psi_tr(i,j,kp)-psi_tr(i,j,k)))*dzci(k)   !*6.0*phiz*(1.d0-phiz)
            !
          enddo
        enddo
      enddo
      !$acc end kernels
    case(0)
      do k=IBM_h+1,nz
        do j=1,ny
          do i=1,nx
            !
            ip = i+1   
            jp = j+1   
            kp = k+1  
            !
            phix   = 0.5_rp*(psi(ip,j,k)+psi(i,j,k))
            phiy   = 0.5_rp*(psi(i,jp,k)+psi(i,j,k))
            phiz   = 0.5_rp*(psi(i,j,kp)+psi(i,j,k))
            kxp    = kappa(ip,j,k)
            kxm    = kappa(i ,j,k)
            kyp    = kappa(i,jp,k)
            kym    = kappa(i,j ,k)
            kzp    = kappa(i,j,kp)
            kzm    = kappa(i,j,k )
            kappax = 0.5_rp*(kxp+kxm)
            kappay = 0.5_rp*(kyp+kym)
            kappaz = 0.5_rp*(kzp+kzm)
            ! Kappa is constrained as in https://doi.org/10.1023/A:1011178417620
            kappax = max(-max_kappa, min(kappax,max_kappa))
            kappay = max(-max_kappa, min(kappay,max_kappa))
            kappaz = max(-max_kappa, min(kappaz,max_kappa))
            ! 
#if defined(_TWOD)
            kappax = 0.0
#endif
            !
            u(i,j,k) = u(i,j,k) + f_t12*(sigma*kappax*(psi_tr(ip,j,k)-psi_tr(i,j,k)))*dxi   !*6.0*phix*(1.d0-phix)  
            v(i,j,k) = v(i,j,k) + f_t12*(sigma*kappay*(psi_tr(i,jp,k)-psi_tr(i,j,k)))*dyi   !*6.0*phiy*(1.d0-phiy) 
            w(i,j,k) = w(i,j,k) + f_t12*(sigma*kappaz*(psi_tr(i,j,kp)-psi_tr(i,j,k)))*dzci(k)   !*6.0*phiz*(1.d0-phiz)
            !
          enddo
        enddo
      enddo
      !$acc end kernels
    case default
      call flutas_error('ERROR: This outflow direction is not implemented ... yet! (surface tension)')
    end select
    !
    deallocate(psi_tr)
    !
    return
  end subroutine surft_src_IBM

  !
  subroutine grav_tw_src(nx,ny,nz,rho,f_t12,cbcpre,dxi,dyi,dzi,nh_d,dzfi, &
#if defined(_BOUSSINESQ)
                         vof,nh_t,tmp, & 
#endif
                         nh_u,u,v,w)  
    !
    use mod_param     , only: lx,ly,lz,gacc_x,gacc_y,gacc_z
#if defined(_BOUSSINESQ)
    use mod_param     , only: tmp0,beta1_th,beta2_th,rho1,rho2
#endif
    !
    implicit none
    ! 
    integer         , intent(in   )                                     :: nx,ny,nz
    real(rp)        , intent(in   ), dimension(     0:,     0:,     0:) :: rho
    real(rp)        , intent(in   )                                     :: f_t12
    character(len=1), intent(in   ), dimension(0:1,3)                   :: cbcpre
    real(rp)        , intent(in   )                                     :: dxi,dyi,dzi
    integer         , intent(in   )                                     :: nh_d
    real(rp)        , intent(in   ), dimension(1-nh_d:)                 :: dzfi
#if defined(_BOUSSINESQ)
    real(rp)        , intent(in   ), dimension(     0:,     0:,     0:) :: vof
    integer         , intent(in   )                                     :: nh_t
    real(rp)        , intent(in   ), dimension(1-nh_t:,1-nh_t:,1-nh_t:) :: tmp
    !@cuf attributes(managed) :: vof,tmp
#endif
    integer         , intent(in   )                                     :: nh_u
    real(rp)        , intent(inout), dimension(1-nh_u:,1-nh_u:,1-nh_u:) :: u,v,w
    !
    real(rp) :: rho_av,rhox,rhoy,rhoz,vofpx,vofpy,vofpz
    real(rp) :: termx,termy,termz
    real(rp) :: tmppx,tmppy,tmppz
    integer  :: i,j,k,ip,jp,kp
    !@cuf attributes(managed) :: rho, dzfi, u, v, w
    !
    ! in case the gravity force acts along a periodic direction
    ! calculate average density (to subtract net gravitational force per unit mass)
    ! 
    rho_av = 0._rp
    if(((cbcpre(0,1)//cbcpre(1,1).eq.'PP').and.gacc_x.ne.0._rp).or. &
       ((cbcpre(0,2)//cbcpre(1,2).eq.'PP').and.gacc_y.ne.0._rp).or. &
       ((cbcpre(0,3)//cbcpre(1,3).eq.'PP').and.gacc_z.ne.0._rp)     ) then
#if defined(_OPENACC)
      !$acc parallel loop reduction(+:rho_av)
#endif
      do k=1,nz
        do j=1,ny
          do i=1,nx
            rho_av = rho_av + rho(i,j,k)/(dxi*dyi*dzfi(k))
          enddo
        enddo
      enddo
#if defined(_OPENACC)
      !$acc end parallel loop 
#endif
      call MPI_ALLREDUCE(MPI_IN_PLACE,rho_av,1,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      rho_av = rho_av/(lx*ly*lz)
    endif
    !
    !$acc kernels
    do k=1,nz
      do j=1,ny
        do i=1,nx
          !
          ip = i+1
          jp = j+1
          kp = k+1
          !
          rhox = 0.5_rp*(rho(ip,j,k)+rho(i,j,k))
          rhoy = 0.5_rp*(rho(i,jp,k)+rho(i,j,k))
          rhoz = 0.5_rp*(rho(i,j,kp)+rho(i,j,k))
          !
#if defined(_BOUSSINESQ)
          tmppx = 0.5_rp*(tmp(ip,j,k)+tmp(i,j,k))
          tmppy = 0.5_rp*(tmp(i,jp,k)+tmp(i,j,k))
          tmppz = 0.5_rp*(tmp(i,j,kp)+tmp(i,j,k))
          !
          vofpx = 0.5_rp*(vof(ip,j,k)+vof(i,j,k))
          vofpy = 0.5_rp*(vof(i,jp,k)+vof(i,j,k))
          vofpz = 0.5_rp*(vof(i,j,kp)+vof(i,j,k))
          !
          termx = rhox - (rho1*beta1_th*vofpx+rho2*beta2_th*(1._rp-vofpx))*(tmppx-tmp0)
          termy = rhoy - (rho1*beta1_th*vofpy+rho2*beta2_th*(1._rp-vofpy))*(tmppy-tmp0)
          termz = rhoz - (rho1*beta1_th*vofpz+rho2*beta2_th*(1._rp-vofpz))*(tmppz-tmp0)
#else
          termx = rhox
          termy = rhoy
          termz = rhoz
#endif
          !
          u(i,j,k) = u(i,j,k) + f_t12*gacc_x*(termx-rho_av)!/rhox
          v(i,j,k) = v(i,j,k) + f_t12*gacc_y*(termy-rho_av)!/rhoy
          w(i,j,k) = w(i,j,k) + f_t12*gacc_z*(termz-rho_av)!/rhoz
          !
        enddo
      enddo
    enddo
    !$acc end kernels
    !
    return
  end subroutine grav_tw_src

subroutine grav_tw_src_ibm(nx,ny,nz,rho,f_t12,cbcpre,dxi,dyi,dzi,nh_d,dzfi, &

#if defined(_BOUSSINESQ)
                         vof,nh_t,tmp, & 
#endif
                         nh_u,u,v,w)  
    !
    use mod_param     , only: lx,ly,lz,gacc_x,gacc_y,gacc_z,IBM_h
#if defined(_BOUSSINESQ)
    use mod_param     , only: tmp0,beta1_th,beta2_th,rho1,rho2
#endif
    !
    implicit none
    ! 
    integer         , intent(in   )                                     :: nx,ny,nz
    real(rp)        , intent(in   ), dimension(     0:,     0:,     0:) :: rho
    real(rp)        , intent(in   )                                     :: f_t12
    character(len=1), intent(in   ), dimension(0:1,3)                   :: cbcpre
    real(rp)        , intent(in   )                                     :: dxi,dyi,dzi
    integer         , intent(in   )                                     :: nh_d
    real(rp)        , intent(in   ), dimension(1-nh_d:)                 :: dzfi
#if defined(_BOUSSINESQ)
    real(rp)        , intent(in   ), dimension(     0:,     0:,     0:) :: vof
    integer         , intent(in   )                                     :: nh_t
    real(rp)        , intent(in   ), dimension(1-nh_t:,1-nh_t:,1-nh_t:) :: tmp
    !@cuf attributes(managed) :: vof,tmp
#endif
    integer         , intent(in   )                                     :: nh_u
    real(rp)        , intent(inout), dimension(1-nh_u:,1-nh_u:,1-nh_u:) :: u,v,w
    !
    real(rp) :: rho_av,rhox,rhoy,rhoz,vofpx,vofpy,vofpz
    real(rp) :: termx,termy,termz
    real(rp) :: tmppx,tmppy,tmppz
    integer  :: i,j,k,ip,jp,kp
    !@cuf attributes(managed) :: rho, dzfi, u, v, w
    !
    ! in case the gravity force acts along a periodic direction
    ! calculate average density (to subtract net gravitational force per unit mass)
    ! 
    rho_av = 0._rp
    if(((cbcpre(0,1)//cbcpre(1,1).eq.'PP').and.gacc_x.ne.0._rp).or. &
       ((cbcpre(0,2)//cbcpre(1,2).eq.'PP').and.gacc_y.ne.0._rp).or. &
       ((cbcpre(0,3)//cbcpre(1,3).eq.'PP').and.gacc_z.ne.0._rp)     ) then
#if defined(_OPENACC)
      !$acc parallel loop reduction(+:rho_av)
#endif
      do k=1+IBM_h,nz
        do j=1,ny
          do i=1,nx
            rho_av = rho_av + rho(i,j,k)/(dxi*dyi*dzfi(k))
          enddo
        enddo
      enddo
#if defined(_OPENACC)
      !$acc end parallel loop 
#endif
      call MPI_ALLREDUCE(MPI_IN_PLACE,rho_av,1,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      rho_av = rho_av/(lx*ly*lz)
    endif
    !
    !$acc kernels
    do k=1+IBM_h,nz
      do j=1,ny
        do i=1,nx
          !
          ip = i+1
          jp = j+1
          kp = k+1
          !
          rhox = 0.5_rp*(rho(ip,j,k)+rho(i,j,k))
          rhoy = 0.5_rp*(rho(i,jp,k)+rho(i,j,k))
          rhoz = 0.5_rp*(rho(i,j,kp)+rho(i,j,k))
          !
#if defined(_BOUSSINESQ)
          tmppx = 0.5_rp*(tmp(ip,j,k)+tmp(i,j,k))
          tmppy = 0.5_rp*(tmp(i,jp,k)+tmp(i,j,k))
          tmppz = 0.5_rp*(tmp(i,j,kp)+tmp(i,j,k))
          !
          vofpx = 0.5_rp*(vof(ip,j,k)+vof(i,j,k))
          vofpy = 0.5_rp*(vof(i,jp,k)+vof(i,j,k))
          vofpz = 0.5_rp*(vof(i,j,kp)+vof(i,j,k))
          !
          termx = rhox - (rho1*beta1_th*vofpx+rho2*beta2_th*(1._rp-vofpx))*(tmppx-tmp0)
          termy = rhoy - (rho1*beta1_th*vofpy+rho2*beta2_th*(1._rp-vofpy))*(tmppy-tmp0)
          termz = rhoz - (rho1*beta1_th*vofpz+rho2*beta2_th*(1._rp-vofpz))*(tmppz-tmp0)
#else
          termx = rhox
          termy = rhoy
          termz = rhoz
#endif
          !
          u(i,j,k) = u(i,j,k) + f_t12*gacc_x*(termx-rho_av)!/rhox
          v(i,j,k) = v(i,j,k) + f_t12*gacc_y*(termy-rho_av)!/rhoy
          w(i,j,k) = w(i,j,k) + f_t12*gacc_z*(termz-rho_av)!/rhoz
          !
        enddo
      enddo
    enddo
    if(IBM_h.NE.0) then   !add force on surface of immersed boudary 
      k = IBM_h
      do j=1,ny
        do i=1,nx
            kp = k+1
            !
            termz = rho(i,j,kp)
            !w(i,j,k) = w(i,j,k) + f_t12*gacc_z*(termz-rho_av)
        enddo
      enddo
    endif
    !$acc end kernels
    !
    return

end subroutine grav_tw_src_ibm

  subroutine stress_IBM(nx,ny,nz,f_t12,cbcpre,dxi,dyi,dzi,nh_d,dzfi,dzci,nh_u,u,v,w,ind_IBM,dudtrko,dvdtrko, &
                                          mu,rho,force_IBM,un,vn,wn,p,pold,f_t12_o,f_t1)  
        !
        use mod_param     , only: lx,ly,lz,gacc_x,gacc_y,gacc_z,rho0,IBM_h
        !
        !  Please check https://doi.org/10.1063/1.1835771 for details
        !  Adapted for adam bashforth scheme
        !
        implicit none
        ! 
        integer         , intent(in   )                                     :: nx,ny,nz
        real(rp)        , intent(in   )                                     :: f_t12
        character(len=1), intent(in   ), dimension(0:1,3)                   :: cbcpre
        real(rp)        , intent(in   )                                     :: dxi,dyi,dzi
        integer         , intent(in   )                                     :: nh_d
        real(rp)        , intent(in   ), dimension(1-nh_d:)                 :: dzfi,dzci
        integer         , intent(in   )                                     :: nh_u
        real(rp)        , intent(in   ), dimension(1-nh_u:,1-nh_u:,1-nh_u:) :: u,v,w,mu,rho
        real(rp)        , intent(inout), dimension(1:,1:,1:)                 :: dudtrko,dvdtrko
        real(rp)        , intent(inout), dimension(1-nh_u:,1-nh_u:,1-nh_u:) :: un,vn,wn,force_IBM
        real(rp)        , intent(in )  , dimension(0:,0:,0:)                :: ind_IBM
        real(rp)        , intent(in   ), dimension(1-nh_u:,1-nh_u:,1-nh_u:) :: p,pold
        real(rp)        , intent(in   )                                     :: f_t12_o,f_t1
        real(rp) :: rhop,rhom
        integer  :: i,j,k,ip,jp,kp,im,jm,km
        real(rp) :: wvkp,wvkm,wvkm_new,dvdzp,dvdzm,dvdzm_new,dvwdz_old,dvwdz_new,d2vdz2_old,d2vdz2_new
        real(rp) :: dwdyp,dwdym,dwdym_new
        real(rp) :: uwkp,uwkm,uwkm_new,dudzp,dudzm,dudzm_new,duwdz_old,duwdz_new,d2udz2_old,d2udz2_new
        real(rp) :: dudt,dvdt,dwdt,muzp,muzm,f1,f2,dwdxp,dwdxm,dwdxm_new,dpdz_old,dpdz_new
        f1 = 1._rp+(f_t12/f_t12_o)
        f2 =       (f_t12/f_t12_o)
        !@cuf attributes(managed) :: dzfi, u, v, w
        !
        !$acc kernels
        force_IBM(:,:,:) = 0.0    
        k = 1 + IBM_h
        do j=1,ny
          do i=1,nx
            !
            ip = i + 1
            jp = j + 1
            kp = k + 1
            im = i - 1
            jm = j - 1
            km = k - 1        
            !
            !!!!MomY!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            dudt = 0.0
            dvdt = 0.0
            dwdt = 0.0
            !if ((ind_IBM(i,j,k).eq.0.0) .and. (ind_IBM(i,j,k-1).eq.1.0)) then
              !
              rhop =     0.50*(rho(  i,jp,k )+rho(  i,jp,kp))
              rhom =     0.50*(rho(  i,j ,k )+rho(  i,j ,kp))
              wvkp  =    0.25*( rhom*w(i  ,j,k  )+rhop*w(i  ,j+1,k  ) )*( v(i,j,k+1)+v(i,j,k) )
              rhop =     0.50*(rho(  i,jp,k )+rho(  i,jp,km))
              rhom =     0.50*(rho(  i,j ,k )+rho(  i,j ,km))
              wvkm  =    0.25*( rhom*w(i  ,j,k-1)+rhop*w(i  ,j+1,k-1) )*( v(i,j,k-1)+v(i,j,k) )

              wvkm_new = 0.0
        
              dvdzp = (v(i,j,kp) -v(i,j,k  ))*dzci(k)
              dvdzm = (v(i,j,k ) -v(i,j,km ))*dzci(km)
              dvdzm_new = (2.0*v(i,j,k ))*dzci(km)
              dwdyp = (w(i,jp,k )-w(i,j,k  ))*dyi
              dwdym = (w(i,jp,km)-w(i,j,km ))*dyi
              dwdym_new = 0.0

              muzp = 0.25_rp*(mu(i,j,k)+mu(i,jp,k)+mu(i,jp,kp)+mu(i,j,kp))
              muzm = 0.25_rp*(mu(i,j,k)+mu(i,jp,k)+mu(i,jp,km)+mu(i,j,km))  

              dvwdz_old =  dzfi(k)*( -wvkp + wvkm )
              dvwdz_new =  dzfi(k)*( -wvkp + wvkm_new )
        
              d2vdz2_old = dzfi(k)*((dvdzp+dwdyp)*muzp-(dvdzm    +dwdym    )*muzm)
              d2vdz2_new = dzfi(k)*((dvdzp+dwdyp)*muzp-(dvdzm_new+dwdym_new)*muzm)

              dvdt = (- dvwdz_old - d2vdz2_old + dvwdz_new + d2vdz2_new)    !*rhoy

              force_IBM(i,j,k) = dvdt

              vn(i,j,k) = vn(i,j,k)  + f_t1*dvdt
              dvdtrko(i,j,k) = dvdtrko(i,j,k)+ dvdt
              !
            !endif
            !
            !!!!MomX!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !
            !if ((ind_IBM(i,j,k).eq.0.0) .and. (ind_IBM(i,j,k-1).eq.1.0)) then
#if defined(_TWOD)
              dudt = 0.0
#else
              rhop  = 0.50*(rho(ip,j ,k )+rho(ip,j ,kp))
              rhom  = 0.50*(rho(i ,j ,k )+rho(i ,j ,kp))
              uwkp  = 0.25_rp*( u(i,j,kp)+u(i,j,k) )*( rhop*w(i+1,j  ,k  )+rhom*w(i,j  ,k  ) )
              rhop  = 0.50*(rho(ip,j ,k )+rho(ip,j ,km))
              rhom  = 0.50*(rho(i ,j ,k )+rho(i ,j ,km)) 
              uwkm  = 0.25_rp*( u(i,j,km)+u(i,j,k) )*( rhop*w(i+1,j  ,k-1)+rhom*w(i,j  ,k-1) )

              uwkm_new = 0.0     !Should not matter
        
              dudzp = (u(i ,j,kp)-u(i ,j,k ))*dzci(k)
              dudzm = (u(i ,j,k )-u(i ,j,km ))*dzci(km)
              dudzm_new = (2.0*u(i ,j,k ))*dzci(km)   !This is what really matters
              dwdxp = (w(ip,j,k )-w(i ,j,k ))*dxi
              dwdxm = (w(ip,j,km)-w(i ,j,km))*dxi
              dwdxm_new = 0.0     !Should not matter
  
              muzp = 0.25_rp*(mu(i,j,k)+mu(i,j,kp)+mu(ip,j,kp)+mu(ip,j,k))
              muzm = 0.25_rp*(mu(i,j,k)+mu(i,j,km)+mu(ip,j,km)+mu(ip,j,k))
        
              duwdz_old = dzfi(k)*( -uwkp + uwkm )
              duwdz_new = dzfi(k)*( -uwkp + uwkm_new )
        
              d2udz2_old = dzfi(k)*((dudzp+dwdxp)*muzp-(dudzm    +dwdxm    )*muzm)
              d2udz2_new = dzfi(k)*((dudzp+dwdxp)*muzp-(dudzm_new+dwdxm_new)*muzm)

              dudt = (- duwdz_old - d2udz2_old + duwdz_new + d2udz2_new)
#endif 
              un(i,j,k)  = un(i,j,k)  + f_t1*dudt
              dudtrko(i,j,k) = dudtrko(i,j,k) + dudt
            !endif
            !
            !!!!MomZ!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! Test
            !dpdz_old = - ( p(i,j,k)-p(i,j,km) )*dzci(km)
            !dpdz_new = - ( p(i,j,k)-0.0       ) *dzci(km) !0.0_rp    ! zero gradient at wall
            !
            !dwdt = dpdz_new - dpdz_old
            !wn(i,j,km)  = wn(i,j,km)  + f_t12*dwdt

          enddo
        enddo
        !$acc end kernels
        !
        return
      end subroutine stress_IBM

  subroutine grav_sp_src(nx,ny,nz,f_t12,cbcpre,dxi,dyi,dzi,nh_d,dzfi, &
#if defined(_BOUSSINESQ)
                         nh_t,tmp, & 
#endif
                         nh_u,u,v,w)  
    !
    use mod_param     , only: lx,ly,lz,gacc_x,gacc_y,gacc_z
#if defined(_BOUSSINESQ)
    use mod_param     , only: tmp0,beta_sp_th
#endif
    !
    implicit none
    ! 
    integer         , intent(in   )                                     :: nx,ny,nz
    real(rp)        , intent(in   )                                     :: f_t12
    character(len=1), intent(in   ), dimension(0:1,3)                   :: cbcpre
    real(rp)        , intent(in   )                                     :: dxi,dyi,dzi
    integer         , intent(in   )                                     :: nh_d
    real(rp)        , intent(in   ), dimension(1-nh_d:)                 :: dzfi
#if defined(_BOUSSINESQ)
    integer         , intent(in   )                                     :: nh_t
    real(rp)        , intent(in   ), dimension(1-nh_t:,1-nh_t:,1-nh_t:) :: tmp
    !@cuf attributes(managed) :: tmp
#endif
    integer         , intent(in   )                                     :: nh_u
    real(rp)        , intent(inout), dimension(1-nh_u:,1-nh_u:,1-nh_u:) :: u,v,w
    !
    real(rp) :: termx,termy,termz
    real(rp) :: tmppx,tmppy,tmppz
    integer  :: i,j,k,ip,jp,kp
    !@cuf attributes(managed) :: dzfi, u, v, w
    !
    !$acc kernels
    do k=1,nz
      do j=1,ny
        do i=1,nx
          !
          ip = i+1
          jp = j+1
          kp = k+1
          !
#if defined(_BOUSSINESQ)
          tmppx = 0.5_rp*(tmp(ip,j,k)+tmp(i,j,k))
          tmppy = 0.5_rp*(tmp(i,jp,k)+tmp(i,j,k))
          tmppz = 0.5_rp*(tmp(i,j,kp)+tmp(i,j,k))
          !
          termx = 1._rp - beta_sp_th*(tmppx-tmp0)
          termy = 1._rp - beta_sp_th*(tmppy-tmp0)
          termz = 1._rp - beta_sp_th*(tmppz-tmp0)
#else
          termx = 1._rp
          termy = 1._rp
          termz = 1._rp
#endif
          !
          u(i,j,k) = u(i,j,k) + f_t12*gacc_x*termx
          v(i,j,k) = v(i,j,k) + f_t12*gacc_y*termy
          w(i,j,k) = w(i,j,k) + f_t12*gacc_z*termz
          !
        enddo
      enddo
    enddo
    !$acc end kernels
    !
    return
  end subroutine grav_sp_src
  !
  subroutine pres_sp_src(nx,ny,nz,f_t12,dxi,dyi,dzi,nh_d,nh_u,dzci,rho0i,pold,u,v,w)
    !
    implicit none 
    !
    integer , intent(in   )                                     :: nx,ny,nz
    real(rp), intent(in   )                                     :: f_t12
    real(rp), intent(in   )                                     :: dxi,dyi,dzi
    integer , intent(in   )                                     :: nh_d,nh_u
    real(rp), intent(in   ), dimension(1-nh_d:)                 :: dzci
    real(rp), intent(in   )                                     :: rho0i
    real(rp), intent(in   ), dimension(0:,0:,0:)                :: pold
    real(rp), intent(inout), dimension(1-nh_u:,1-nh_u:,1-nh_u:) :: u,v,w
    !
    integer :: i,j,k,ip,jp,kp
    !@cuf attributes(managed) :: dzci, pold, u, v, w
    !
    !$acc kernels
    do k=1,nz
      do j=1,ny
        do i=1,nx
          !
          ip = i+1
          jp = j+1
          kp = k+1
          !
          u(i,j,k) = u(i,j,k) + f_t12*( - ( pold(ip,j,k)-pold(i,j,k) )*dxi     )*rho0i
          v(i,j,k) = v(i,j,k) + f_t12*( - ( pold(i,jp,k)-pold(i,j,k) )*dyi     )*rho0i
          w(i,j,k) = w(i,j,k) + f_t12*( - ( pold(i,j,kp)-pold(i,j,k) )*dzci(k) )*rho0i
          !
        enddo
      enddo
    enddo
    !$acc end kernels
    !
    return
  end subroutine pres_sp_src
  !
#if defined(_TURB_FORCING)
  subroutine forc_src(nx,ny,nz,nh_d,nh_u,f_t12,dx,dy,dz,zc,rho,u,v,w)
    !
    use mod_param, only: abc_x,abc_y,abc_z,lx,ly,lz,f0_t,k0_t,pi,turb_type
    !
    implicit none
    ! 
    integer , intent(in   )                                     :: nx,ny,nz
    integer , intent(in   )                                     :: nh_d,nh_u
    real(rp), intent(in   )                                     :: f_t12
    real(rp), intent(in   )                                     :: dx,dy,dz
    real(rp), intent(in   ), dimension(1-nh_d:)                 :: zc
    real(rp), intent(in   ), dimension(0:,0:,0:)                :: rho
    real(rp), intent(inout), dimension(1-nh_u:,1-nh_u:,1-nh_u:) :: u,v,w
    !
    real(rp) :: xcc,ycc,zcc,xff,yff
    real(rp) :: rhox,rhoy,rhoz
    integer  :: i,j,k,ijk_start_x,ijk_start_y,ijk_start_z
    !@cuf attributes(managed) :: u, v, w, zc
    !
    ijk_start_x = ijk_start(1)
    ijk_start_y = ijk_start(2)
    ijk_start_z = ijk_start(3)
    !
    !
    select case(turb_type)
!    case('tgv')
!      !$acc kernels
!      do k=1,nz
!        do j=1,ny
!          do i=1,nx
!            !
!            !zcc = zc(k)/lz*2._rp*pi
!            zcc = (k+ijk_start_z-0.5_rp)*dz/lz*2._rp*pi ! in triperiodic we typically use unstretched grid
!            ycc = (j+ijk_start_y-0.5_rp)*dy/ly*2._rp*pi
!            yff = (j+ijk_start_y-0.0_rp)*dy/ly*2._rp*pi
!            xcc = (i+ijk_start_x-0.5_rp)*dx/lx*2._rp*pi
!            xff = (i+ijk_start_x-0.0_rp)*dx/lx*2._rp*pi
!            !
!            u(i,j,k) = u(i,j,k) + f_t12*( + f0_t*sin(k0_t*xff)*cos(k0_t*ycc)*cos(k0_t*zcc) )
!            v(i,j,k) = v(i,j,k) + f_t12*( - f0_t*cos(k0_t*xcc)*sin(k0_t*yff)*cos(k0_t*zcc) )
!            !
!          enddo
!        enddo
!      enddo
!      !$acc end kernels
    case('abc')
      !$acc kernels
      do k=1,nz
        do j=1,ny
          do i=1,nx
            !
            !zcc = zc(k)/lz*2._rp*pi
            zcc = (k+ijk_start_z-0.5_rp)*dz/lz*2._rp*pi ! in triperiodic we typically use unstretched grid
            ycc = (j+ijk_start_y-0.5_rp)*dy/ly*2._rp*pi
            yff = (j+ijk_start_y-0.0_rp)*dy/ly*2._rp*pi
            xcc = (i+ijk_start_x-0.5_rp)*dx/lx*2._rp*pi
            xff = (i+ijk_start_x-0.0_rp)*dx/lx*2._rp*pi
            !
            rhox = (rho(i,j,k)+rho(i+1,j,k))*0.5d0
            rhoy = (rho(i,j,k)+rho(i,j+1,k))*0.5d0
            rhoz = (rho(i,j,k)+rho(i,j,k+1))*0.5d0
            !
            u(i,j,k) = u(i,j,k) + f_t12*( + f0_t*(abc_x*sin(k0_t*zcc) + abc_z*cos(k0_t*ycc)) )*rhox
            v(i,j,k) = v(i,j,k) + f_t12*( + f0_t*(abc_y*sin(k0_t*xcc) + abc_x*cos(k0_t*zcc)) )*rhoy
            w(i,j,k) = w(i,j,k) + f_t12*( + f0_t*(abc_z*sin(k0_t*ycc) + abc_y*cos(k0_t*xcc)) )*rhoz
            !
          enddo
        enddo
      enddo
      !$acc end kernels
    case default
      call flutas_error('Invalid forcing for HIT - check forcing.in')
    end select
    !
    return
  end subroutine forc_src
#endif
  !
  subroutine bulk_forcing_src(forcing_type,is_forced,nx,ny,nz,dx,dy,dz,f_t12, &
                              nh_d,nh_u,dzc,dzf,rho,u,v,w,f)
    !
    ! to add bulk velocity forcing following three control strategies:
    !   --> constant flow rate (CFR);
    !   --> constant pressure gradient (CPR);
    !   --> constant power input (CPI) - not implemented yet.
    !
    use mod_debug, only: cmpt_mean
    use mod_param, only: lx,ly,lz,bvel_x,bvel_y,bvel_z,dpdl_x,dpdl_y,dpdl_z
#if !defined(_TWO_PHASE)
    use mod_param, only: rho_sp
#endif
    !
    implicit none
    !
    character(len=3), intent(in   )                                     :: forcing_type
    logical         , intent(in   ), dimension(3)                       :: is_forced
    integer         , intent(in   )                                     :: nx,ny,nz
    real(rp)        , intent(in   )                                     :: dx,dy,dz
    real(rp)        , intent(in   )                                     :: f_t12
    integer         , intent(in   )                                     :: nh_d,nh_u
    real(rp)        , intent(in   ), dimension(1-nh_d:)                 :: dzc,dzf
    real(rp)        , intent(in   ), dimension(     0:,     0:,     0:) :: rho
    real(rp)        , intent(inout), dimension(1-nh_u:,1-nh_u:,1-nh_u:) :: u,v,w
    real(rp)        , intent(out  ), dimension(3)                       :: f
    !
    real(rp) :: meanvel,rhof
    integer  :: i,j,k
    !@cuf attributes(managed) :: dzc, dzf, u, v, w, rho
    !
#if !defined(_TWO_PHASE)
    rhof = rho_sp
#endif
    !
    f(1:3) = 0._rp
    select case(forcing_type)
    case('cfr')
      !
      ! ensure the flow rate is exactly constant
      !  note: to be generalized when the density field is not uniform
      !  (e.g., two-phase, low-mach flows)
      !
      if(is_forced(1)) then
        call cmpt_mean(nx,ny,nz,nh_d,nh_u,dx,dy,dzf,lx,ly,lz,u,meanvel)
        !$acc kernels
        do k=1,nz
          do j=1,ny
            do i=1,nx
              u(i,j,k) = u(i,j,k) + (bvel_x-meanvel)
            enddo
          enddo
        enddo
        !$acc end kernels
        f(1) = bvel_x-meanvel
      endif
      if(is_forced(2)) then
        call cmpt_mean(nx,ny,nz,nh_d,nh_u,dx,dy,dzf,lx,ly,lz,v,meanvel)
        !$acc kernels
        do k=1,nz
          do j=1,ny
            do i=1,nx
              v(i,j,k) = v(i,j,k) + (bvel_y-meanvel)
            enddo
          enddo
        enddo
        !$acc end kernels
        f(2) = bvel_y-meanvel
      endif
      if(is_forced(3)) then
        call cmpt_mean(nx,ny,nz,nh_d,nh_u,dx,dy,dzc,lx,ly,lz,w,meanvel)
        !$acc kernels
        do k=1,nz
          do j=1,ny
            do i=1,nx
              w(i,j,k) = w(i,j,k) + (bvel_z-meanvel)
            enddo
          enddo
        enddo
        !$acc end kernels
        f(3) = bvel_z-meanvel
      endif
    case('cpr')
      !
      ! ensure the pressure gradient is exactly constant
      !
      if(is_forced(1)) then
        !$acc kernels
        do k=1,nz
          do j=1,ny
            do i=1,nx
#if defined(_TWO_PHASE)
              rhof     = 0.5_rp*(rho(i+1,j,k)+rho(i,j,k))
#endif
              u(i,j,k) = u(i,j,k) + f_t12*dpdl_x/rhof
            enddo
          enddo
        enddo
        !$acc end kernels
      endif
      if(is_forced(2)) then
        !$acc kernels
        do k=1,nz
          do j=1,ny
            do i=1,nx
#if defined(_TWO_PHASE)
              rhof     = 0.5_rp*(rho(i,j+1,k)+rho(i,j,k))
#endif
              v(i,j,k) = v(i,j,k) + f_t12*dpdl_y/rhof
            enddo
          enddo
        enddo
        !$acc end kernels
      endif
      if(is_forced(3)) then
        !$acc kernels
        do k=1,nz
          do j=1,ny
            do i=1,nx
#if defined(_TWO_PHASE)
              rhof     = 0.5_rp*(rho(i,j,k+1)+rho(i,j,k))
#endif
              w(i,j,k) = w(i,j,k) + f_t12*dpdl_z/rhof
            enddo
          enddo
        enddo
        !$acc end kernels
      endif
    case('cpi')
      !
      ! ensure the power input is exactly constant
      !
      if(is_forced(1).or.is_forced(2).or.is_forced(3)) then
        call flutas_error('CPI control strategy not implemented yet')
      endif
    case default
      call flutas_error('Invalid selected control strategy - check dns.in')
    end select
    !
    return
  end subroutine bulk_forcing_src
  !
end module mod_source
