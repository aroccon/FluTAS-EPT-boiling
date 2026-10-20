!
!
! SPDX-FileCopyrightText: Copyright (c) 2017-2022 Pedro Costa and the CaNS contributors.
! SPDX-FileCopyrightText: Modifications Copyright (c) 2023-2024 Maochao Xiao and the CaLES contributors.
! SPDX-License-Identifier: MIT
! Taken from CaNS-World GitHub on 3rd July 2025
! Modified for two-phase flow FluTAS by Lorenz Weber lorenz.weber@kit.edu
! Assumptions: only static smagorinsky and no wall model
!
module mod_sgs
  use mpi
  use mod_types
  implicit none
  private
  public cmpt_sgs
  contains
  !
  subroutine cmpt_sgs(n,dzc,dzf,dzci,dzfi,visc,u,v,w,visct,rho)
    use mod_param     , only: l,dl,dli,cbcvel,IBM_h
    use mod_boiling   , only: sgstype,c_smag,wall
    use mod_common_mpi, only: MPI_PROC_NULL,left,right,front,back,top,bottom,myid
    !
    ! compute subgrid viscosity at cell centers
    ! the LES with the dynamic model is ~2 times the cost of the LES with the static one.
    ! Acceleration can be further achieved by calling only one boundp to do the data
    ! exchange of multiple variables. Note that it does not save time to simply mask
    ! wall bc's in boundp. The dynamic version yields quite good results for Re_tau=
    ! 395,550,1000, with <=5% errors in the friction coefficient. Clipping is necessary
    ! to avoid negative values of averaged eddy viscosity (duct flow).
    !
    ! 2D filter is used in the first off-wall layer (Bae, Orlandi, LESGO and Balaras, 1995).
    ! It is difficult to do 3D filtering of Sij in the first layer, though feasible for
    ! velocity.
    !
    implicit none
    integer , intent(in ), dimension(3)                      :: n
    real(rp), intent(in ), dimension(0:)                     :: dzc,dzf,dzci,dzfi
    real(rp), intent(in ), dimension(0:,0:,0:)               :: u,v,w,rho,visc
    real(rp), intent(out), dimension(0:,0:,0:)               :: visct
    real(rp), allocatable, dimension(:,:,:)  , save          :: s0
    real(rp), target,  allocatable, dimension(:,:,:,:), save :: wk
    real(rp)                                                 :: dxi,dyi,visci,tauw(2),dw(6),dw_plus, &
                                                                 fd,del,dw_min,tauw_s
    real(rp), save                                           :: is_wall(6)
    logical, save                                            :: is_first = .true.
    integer                                                  :: i,j,k,m,loc
    real(rp), parameter                                      :: one_third = 1._rp/3._rp
    real(rp), parameter                                      :: big = huge(1.0_rp)
    logical , dimension(0:1,3), save                         :: is_bound
    real(rp)                                                 :: ibm_factor
    integer                                                  :: n1,n2,n3
    integer, save                                            :: ncall_dbg = 0
    integer                                                  :: ibmh_dev
    real(rp)                                                 :: dl1_dev
    !
    select case(trim(sgstype))
    case('none')
      if(is_first) then
        is_first = .false.
        !$acc kernels default(present) async(1)
        visct(:,:,:) = 0._rp
        !$acc end kernels
      end if
    case('smag')
      if(is_first) then
        is_first = .false.
        allocate(s0(0:n(1)+1,0:n(2)+1,0:n(3)+1  ), &
                 wk(0:n(1)+1,0:n(2)+1,0:n(3)+1,3))
        !$acc enter data create(s0,wk) async(1)
        is_wall   = 0._rp
        is_bound  = .FALSE.
        !
        if(left  .eq.MPI_PROC_NULL)    is_bound(0,1) = .TRUE.
        if(right .eq.MPI_PROC_NULL)    is_bound(1,1) = .TRUE.
        if(front .eq.MPI_PROC_NULL)    is_bound(0,2) = .TRUE.
        if(back  .eq.MPI_PROC_NULL)    is_bound(1,2) = .TRUE.
        if(bottom.eq.MPI_PROC_NULL)    is_bound(0,3) = .TRUE.
        if(top   .eq.MPI_PROC_NULL)    is_bound(1,3) = .TRUE.
        !
        if(wall(0,1).and.is_bound(0,1)) is_wall(1) = 1._rp
        if(wall(1,1).and.is_bound(1,1)) is_wall(2) = 1._rp
        if(wall(0,2).and.is_bound(0,2)) is_wall(3) = 1._rp
        if(wall(1,2).and.is_bound(1,2)) is_wall(4) = 1._rp
        if(wall(0,3).and.is_bound(0,3)) is_wall(5) = 1._rp
        if(wall(1,3).and.is_bound(1,3)) is_wall(6) = 1._rp
        !$acc enter data copyin(is_wall) async(1)
      end if
      !$acc kernels default(present) async(1)
      wk(:,:,:,1) = u(:,:,:)
      wk(:,:,:,2) = v(:,:,:)
      wk(:,:,:,3) = w(:,:,:)
      !$acc end kernels
      call extrapolate(n,is_bound,dzci,wk(:,:,:,1),1,wall)
      call extrapolate(n,is_bound,dzci,wk(:,:,:,2),2,wall)
      call extrapolate(n,is_bound,dzci,wk(:,:,:,3),3,wall)
      call strain_rate(n,dli,dzci,dzfi,wk(:,:,:,1),wk(:,:,:,2),wk(:,:,:,3),s0)
      !
      dxi = dli(1)
      dyi = dli(2)
      !
      ! scalar copies of n: the host array n is not present on the device (default(present))
      n1 = n(1); n2 = n(2); n3 = n(3)
      !
      ibm_factor = 1.0
      if(IBM_h.NE.0) ibm_factor = 2.0
      !
      !$acc parallel loop collapse(3) default(present) async(1) &
      !$acc private(dw,tauw,loc,fd,del,tauw_s,dw_min,dw_plus)
      do k=1+IBM_h,n3
        do j=1,n2
          do i=1,n1
            !
            if(sum(is_wall(1:6))==0) then
              ! triperiodic
              fd = 1._rp
            else
              ! van Driest damping
              dw(1) = dl(1)*(i-0.5)
              dw(2) = dl(1)*(n1-i+0.5)
              dw(3) = dl(2)*(j-0.5)
              dw(4) = dl(2)*(n2-j+0.5)
              dw(5) = dl(3)*(k-IBM_h-0.5)
              dw(6) = dl(3)*(n3-j+0.5)
              dw = dw*is_wall + big*(1._rp-is_wall)
              loc = minloc(dw,1)
              dw_min = dw(loc)
              !
              if(loc==1) then
                tauw(1) = v(1,j,k)-v(0,j,k)+v(1,j-1,k)-v(0,j-1,k)
                tauw(2) = w(1,j,k)-w(0,j,k)+w(1,j,k-1)-w(0,j,k-1)
                tauw_s = sqrt(tauw(1)*tauw(1)+tauw(2)*tauw(2))*dxi
                tauw_s   =  0.5_rp*visc(1,j,k)*tauw_s/rho(1,j,k)
                visci    = rho(1,j,k)/visc(1,j,k)
              else if(loc==2) then
                tauw(1) = v(n1,j,k)-v(n1+1,j,k)+v(n1,j-1,k)-v(n1+1,j-1,k)
                tauw(2) = w(n1,j,k)-w(n1+1,j,k)+w(n1,j,k-1)-w(n1+1,j,k-1)
                tauw_s = sqrt(tauw(1)*tauw(1)+tauw(2)*tauw(2))*dxi
                tauw_s   =  0.5_rp*visc(n1,j,k)*tauw_s/rho(n1,j,k)
                visci    = rho(n1,j,k)/visc(n1,j,k)
              else if(loc==3) then
                tauw(1) = u(i,1,k)-u(i,0,k)+u(i-1,1,k)-u(i-1,0,k)
                tauw(2) = w(i,1,k)-w(i,0,k)+w(i,1,k-1)-w(i,0,k-1)
                tauw_s  = sqrt(tauw(1)*tauw(1)+tauw(2)*tauw(2))*dyi
                tauw_s   =  0.5_rp*visc(i,1,k)*tauw_s/rho(i,1,k)
                visci    = rho(i,1,k)/visc(i,1,k)
              else if(loc==4) then
                tauw(1) = u(i,n2,k)-u(i,n2+1,k)+u(i-1,n2,k)-u(i-1,n2+1,k)
                tauw(2) = w(i,n2,k)-w(i,n2+1,k)+w(i,n2,k-1)-w(i,n2+1,k-1)
                tauw_s  = sqrt(tauw(1)*tauw(1)+tauw(2)*tauw(2))*dyi
                tauw_s   =  0.5_rp*visc(i,n2,k)*tauw_s/rho(i,n2,k)
                visci    = rho(i,n2,k)/visc(i,n2,k)
              else if(loc==5) then
                tauw(1) = u(i,j,1+IBM_h)-u(i,j,IBM_h)+u(i-1,j,1+IBM_h)-u(i-1,j,IBM_h)
                tauw(2) = v(i,j,1+IBM_h)-v(i,j,IBM_h)+v(i,j-1,1+IBM_h)-v(i,j-1,IBM_h)
                tauw_s  = sqrt(tauw(1)*tauw(1)+tauw(2)*tauw(2))*dzci(0)*ibm_factor
                tauw_s   =  0.5_rp*visc(i,j,1+IBM_h)*tauw_s/rho(i,j,1+IBM_h)
                visci    = rho(i,j,1+IBM_h)/visc(i,j,1+IBM_h)
              else if(loc==6) then
                tauw(1) = u(i,j,n3)-u(i,j,n3+1)+u(i-1,j,n3)-u(i-1,j,n3+1)
                tauw(2) = v(i,j,n3)-v(i,j,n3+1)+v(i,j-1,n3)-v(i,j-1,n3+1)
                tauw_s  = sqrt(tauw(1)*tauw(1)+tauw(2)*tauw(2))*dzci(n3)
                tauw_s   =  0.5_rp*visc(i,j,n3)*tauw_s/rho(i,j,n3)
                visci    = rho(i,j,n3)/visc(i,j,n3)
              end if
              dw_plus  =  dw_min*sqrt(tauw_s)*visci
              fd       =  1._rp-exp(-dw_plus/25._rp)
            end if
            !
            del = (dl(1)*dl(2)*dzf(k))**(one_third)
            visct(i,j,k) = rho(i,j,k)*(c_smag*del*fd)**2*s0(i,j,k)
          end do
        end do
      end do
    case default
      print*, 'ERROR: unknown SGS model'
    end select
    !
    ! the kernels above run on async queue 1, while the rest of the code uses the default
    ! (synchronous) queue: wait here so that visct is complete before it is used
    !$acc wait(1)
    !
    ! debug (temporary): per-rank checksums of the SGS intermediates, first 2 calls only
    if(trim(sgstype) == 'smag' .and. ncall_dbg < 2) then
      ncall_dbg = ncall_dbg + 1
      write(*,'(a,i3,a,i3,a,6es14.6)') ' SGSDBG call',ncall_dbg,' rank',myid, &
        ' sum|wk1|,|wk2|,|wk3|,|s0|,|visct|,max visct:', &
        sum(abs(wk(1:n(1),1:n(2),1:n(3),1))),sum(abs(wk(1:n(1),1:n(2),1:n(3),2))), &
        sum(abs(wk(1:n(1),1:n(2),1:n(3),3))),sum(abs(s0(1:n(1),1:n(2),1:n(3)))), &
        sum(abs(visct(1:n(1),1:n(2),1:n(3)))),maxval(visct(1:n(1),1:n(2),1:n(3)))
      ! what the GPU sees for module variables of mod_param (host value vs device value)
      ibmh_dev = -1
      dl1_dev  = -1._rp
      !$acc parallel loop reduction(max:ibmh_dev,dl1_dev)
      do i=1,1
        ibmh_dev = max(ibmh_dev,IBM_h)
        dl1_dev  = max(dl1_dev,dl(1))
      enddo
      write(*,'(a,i3,a,i3,a,2i4,a,2es12.4)') ' SGSDBG call',ncall_dbg,' rank',myid, &
        ' IBM_h host/device:',IBM_h,ibmh_dev,' dl(1) host/device:',dl(1),dl1_dev
      write(*,'(a,i3,a,i3,a,6f4.0,a,6l2)') ' SGSDBG call',ncall_dbg,' rank',myid,' is_wall:',is_wall, &
        ' is_bound:',is_bound(0,1),is_bound(1,1),is_bound(0,2),is_bound(1,2),is_bound(0,3),is_bound(1,3)
    endif
  end subroutine cmpt_sgs
  !
  !
  subroutine extrapolate(n,is_bound,dzci,p,iface,wall)
    use mod_param     , only: IBM_h
    !
    ! linear extrapolation of wall-parallel velocity
    !
    ! called before filter/strain_rate for 2D filter/one-sided difference. The
    ! extrapolation requires >= 2 off-wall layers. Wall-normal velocity is
    ! extrapolated when stored at cell centers (iface=0), but is not
    ! when stored at cell faces (iface=1,2,3). i.e., no-penetration bc.
    !
    implicit none
    integer , intent(in)   , dimension(3)        :: n
    logical , intent(in)   , dimension(0:1,3)    :: is_bound,wall
    real(rp), intent(in)   , dimension(0:)       :: dzci
    real(rp), intent(inout), dimension(0:,0:,0:) :: p
    integer , intent(in)                         :: iface
    logical, dimension(0:1,3)                    :: is_done
    real(rp)                                     :: dzc(0:n(3)+1),factor0,factor1
    integer                                      :: i,j,k
    integer                                      :: n1,n2,n3
    !
    dzc = 1._rp/dzci
    !
    ! scalar copies of n: the host array n is not present on the device (default(present))
    n1 = n(1); n2 = n(2); n3 = n(3)
    !
    factor0 = 1._rp   ! all uniform grid
    factor1 = 1._rp   ! all uniform grid
    !
    is_done(:,1) = (is_bound(:,1).and.wall(:,1).and.iface/=1)
    is_done(:,2) = (is_bound(:,2).and.wall(:,2).and.iface/=2)
    is_done(:,3) = (is_bound(:,3).and.wall(:,3).and.iface/=3)
    !
    if(is_done(0,1)) then
      !$acc parallel loop collapse(2) default(present) async(1)
      do k = 0,n3+1
        do j = 0,n2+1
          p(0,j,k) = 2._rp*p(1,j,k) - p(2,j,k)
        end do
      end do
    end if
    if(is_done(1,1)) then
      !$acc parallel loop collapse(2) default(present) async(1)
      do k = 0,n3+1
        do j = 0,n2+1
          p(n1+1,j,k) = 2._rp*p(n1,j,k) - p(n1-1,j,k)
        end do
      end do
    end if
    if(is_done(0,2)) then
      !$acc parallel loop collapse(2) default(present) async(1)
      do k = 0,n3+1
        do i = 0,n1+1
          p(i,0,k) = 2._rp*p(i,1,k) - p(i,2,k)
        end do
      end do
    end if
    if(is_done(1,2)) then
      !$acc parallel loop collapse(2) default(present) async(1)
      do k = 0,n3+1
        do i = 0,n1+1
          p(i,n2+1,k) = 2._rp*p(i,n2,k) - p(i,n2-1,k)
        end do
      end do
    end if
    if(is_done(0,3)) then
      !$acc parallel loop collapse(2) default(present) async(1)
      do j = 0,n2+1
        do i = 0,n1+1
          p(i,j,0+IBM_h) = (1._rp+factor0)*p(i,j,1+IBM_h) - factor0*p(i,j,2+IBM_h)
        end do
      end do
    end if
    if(is_done(1,3)) then
      !$acc parallel loop collapse(2) default(present) async(1)
      do j = 0,n2+1
        do i = 0,n1+1
          p(i,j,n3+1) = (1._rp+factor1)*p(i,j,n3) - factor1*p(i,j,n3-1)
        end do
      end do
    end if
  end subroutine extrapolate
  !
  !
  !
  subroutine strain_rate(n,dli,dzci,dzfi,u,v,w,s0)
    !
    ! compute the strain rate field
    !
    ! Sij is first computed at (or averaged to) cell center, then s0=sqrt(2SijSij)
    ! at cell center (Bae and Orlandi's codes). The implementation is efficient, since
    ! it avoids repetitive computation of derivatives. Costa first averages SijSij to
    ! cell center, then computes s0. This leads to larger s0, especially when Sij
    ! at the cell edges have opposite signs. The second and third loops cannot combine.
    !
    implicit none
    integer , intent(in ), dimension(3)        :: n
    real(rp), intent(in ), dimension(3)        :: dli
    real(rp), intent(in ), dimension(0:)       :: dzci,dzfi
    real(rp), intent(in ), dimension(0:,0:,0:) :: u,v,w
    real(rp), intent(out), dimension(0:,0:,0:) :: s0
    real(rp) :: s11,s22,s33,s12,s13,s23
    real(rp) :: dxi,dyi
    integer :: i,j,k
    integer :: n1,n2,n3
    real(rp) :: u_mcm,u_ccm,u_mmc,u_cmc,u_mcc,u_ccc,u_mpc,u_cpc,u_mcp,u_ccp, &
                v_cmm,v_ccm,v_mmc,v_cmc,v_pmc,v_mcc,v_ccc,v_pcc,v_cmp,v_ccp, &
                w_cmm,w_mcm,w_ccm,w_pcm,w_cpm,w_cmc,w_mcc,w_ccc,w_pcc,w_cpc
    !
    dxi = dli(1)
    dyi = dli(2)
    ! scalar copies of n: the host array n is not present on the device (default(present))
    n1 = n(1); n2 = n(2); n3 = n(3)
    !
    !$acc parallel loop collapse(3) default(present) async(1) &
    !$acc private(u_mcm,u_ccm,u_mmc,u_cmc,u_mcc,u_ccc,u_mpc,u_cpc,u_mcp,u_ccp) &
    !$acc private(v_cmm,v_ccm,v_mmc,v_cmc,v_pmc,v_mcc,v_ccc,v_pcc,v_cmp,v_ccp) &
    !$acc private(w_cmm,w_mcm,w_ccm,w_pcm,w_cpm,w_cmc,w_mcc,w_ccc,w_pcc,w_cpc) &
    !$acc private(s11,s22,s33,s12,s13,s23)
    do k=1,n3
      do j=1,n2
        do i=1,n1
          !
          u_mcm = u(i-1,j  ,k-1)
          u_ccm = u(i  ,j  ,k-1)
          u_mmc = u(i-1,j-1,k  )
          u_cmc = u(i  ,j-1,k  )
          u_mcc = u(i-1,j  ,k  )
          u_ccc = u(i  ,j  ,k  )
          u_mpc = u(i-1,j+1,k  )
          u_cpc = u(i  ,j+1,k  )
          u_mcp = u(i-1,j  ,k+1)
          u_ccp = u(i  ,j  ,k+1)
          !
          v_cmm = v(i  ,j-1,k-1)
          v_ccm = v(i  ,j  ,k-1)
          v_mmc = v(i-1,j-1,k  )
          v_cmc = v(i  ,j-1,k  )
          v_pmc = v(i+1,j-1,k  )
          v_mcc = v(i-1,j  ,k  )
          v_ccc = v(i  ,j  ,k  )
          v_pcc = v(i+1,j  ,k  )
          v_cmp = v(i  ,j-1,k+1)
          v_ccp = v(i  ,j  ,k+1)
          !
          w_cmm = w(i  ,j-1,k-1)
          w_mcm = w(i-1,j  ,k-1)
          w_ccm = w(i  ,j  ,k-1)
          w_pcm = w(i+1,j  ,k-1)
          w_cpm = w(i  ,j+1,k-1)
          w_cmc = w(i  ,j-1,k  )
          w_mcc = w(i-1,j  ,k  )
          w_ccc = w(i  ,j  ,k  )
          w_pcc = w(i+1,j  ,k  )
          w_cpc = w(i  ,j+1,k  )
          !
          s11 = (u_ccc-u_mcc)*dxi    
          s22 = (v_ccc-v_cmc)*dyi    
          s33 = (w_ccc-w_ccm)*dzfi(k)
          s12 = .125_rp*((u_cpc-u_ccc)*dyi + (v_pcc-v_ccc)*dxi + &
                         (u_ccc-u_cmc)*dyi + (v_pmc-v_cmc)*dxi + &
                         (u_mpc-u_mcc)*dyi + (v_ccc-v_mcc)*dxi + &
                         (u_mcc-u_mmc)*dyi + (v_cmc-v_mmc)*dxi)
          s13 = .125_rp*((u_ccp-u_ccc)*dzci(k  ) + (w_pcc-w_ccc)*dxi + &
                         (u_ccc-u_ccm)*dzci(k-1) + (w_pcm-w_ccm)*dxi + &
                         (u_mcp-u_mcc)*dzci(k  ) + (w_ccc-w_mcc)*dxi + &
                         (u_mcc-u_mcm)*dzci(k-1) + (w_ccm-w_mcm)*dxi)
          s23 = .125_rp*((v_ccp-v_ccc)*dzci(k  ) + (w_cpc-w_ccc)*dyi + &
                         (v_ccc-v_ccm)*dzci(k-1) + (w_cpm-w_ccm)*dyi + &
                         (v_cmp-v_cmc)*dzci(k  ) + (w_ccc-w_cmc)*dyi + &
                         (v_cmc-v_cmm)*dzci(k-1) + (w_ccm-w_cmm)*dyi)
          s0(i,j,k) = sqrt(max(0._rp,  2._rp*(                             &     !Add a max(...,0) to avoid complex numbers
                                      -2._rp*(s11*s22+s11*s33+s22*s33)     &     !This was: s11**2+s22**2+s33**2
                                      +2._rp*(s12**2 +s13**2 +s23**2 ))))
        end do
      end do
    end do
  end subroutine strain_rate
end module mod_sgs