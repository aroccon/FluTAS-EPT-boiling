!
! SPDX-License-Identifier: MIT
!
module mod_bound
  !
  use mpi
  use mod_types
  !@cuf use cudafor
  !
  implicit none
  !
  private
  public  :: bounduvw, boundp, updt_rhs_b, boundp_di, boundp_IBM, bounduvw_mom

public boundsb_array,boundsb_array_ibm,updthalo



  !
  contains
  !
  subroutine bounduvw_mom(cbc,n,bc,nh_d,nh_u,halo,isoutflow,dl,dzc,dzf,u,v,w,rho,mfx,mfy,mfz,nse_switch,reg_switch,sem_switch)
    ! same as bounduvw but for the momentum, i.e. Dirichlet BC's need to be corrected with density
    use mod_common_mpi, only: MPI_PROC_NULL,left,right,front,back,top,bottom,myid
    use mod_sem       , only: SEM_apply
    !
    ! imposes velocity boundary conditions 
    ! on an arbitrary stecil size
    !
    implicit none
    !
    character(len=1), intent(in   ), dimension(0:1,3,3)                 :: cbc
    integer         , intent(in   ), dimension(3)                       :: n 
    logical         , intent(in   )                                     :: sem_switch
    real(rp)        , intent(in   ), dimension(0:1,3,3)                 :: bc
    integer         , intent(in   )                                     :: nh_d,nh_u
    integer         , intent(in   ), dimension(3)                       :: halo
    real(rp)        , intent(in   ), dimension(3)                       :: dl
    real(rp)        , intent(in   )                                     :: nse_switch,reg_switch
    real(rp)        , intent(in   ), dimension(1-nh_d:)                 :: dzc,dzf
    real(rp)        , intent(inout), dimension(1-nh_u:,1-nh_u:,1-nh_u:) :: u,v,w
    real(rp)        , intent(in   ), dimension(1-nh_u:,1-nh_u:,1-nh_u:) :: mfx,mfy,mfz
    real(rp)        , intent(in   ), dimension(1-nh_u:,1-nh_u:,1-nh_u:) :: rho
    logical         , intent(in   ), dimension(0:1,3)                   :: isoutflow
    !@cuf attributes(managed) :: u,v,w
    !
    real(rp), dimension(0:nh_u-1) :: dr
    real(rp) :: dx,dy,dz
    integer  :: q,idir,sgn,ioutflowdir,qmin
    integer  :: ind1,ind2
    !
    dx = dl(1)
    dy = dl(2)
    dz = dl(3)
    !
    qmin = abs(1-nh_u)
    !
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
    !call updthalo(n(1),n(2),n(3),nh_u,halo(ind1),ind1,u)
    !call updthalo(n(1),n(2),n(3),nh_u,halo(ind2),ind2,u)
    !
    !call updthalo(n(1),n(2),n(3),nh_u,halo(ind1),ind1,v)
    !call updthalo(n(1),n(2),n(3),nh_u,halo(ind2),ind2,v)
    !
    !call updthalo(n(1),n(2),n(3),nh_u,halo(ind1),ind1,w)
    !call updthalo(n(1),n(2),n(3),nh_u,halo(ind2),ind2,w)
    !
    ! along x
    !
    do q = 0,qmin
      dr(q) = dl(1)
    enddo
    if(left .eq.MPI_PROC_NULL) then ! x - bottom boundary
      call set_bc_mom(n(1),n(2),n(3),cbc(0,1,1),0,1,.false.,bc(0,1,1),qmin,nh_u,dr,u,rho) ! changed
      call set_bc(n(1),n(2),n(3),cbc(0,1,2),0,1,.true. ,bc(0,1,2),qmin,nh_u,dr,v)
      call set_bc(n(1),n(2),n(3),cbc(0,1,3),0,1,.true. ,bc(0,1,3),qmin,nh_u,dr,w)
    endif
    if(right.eq.MPI_PROC_NULL) then ! x - top boundary
      call set_bc_mom(n(1),n(2),n(3),cbc(1,1,1),1,1,.false.,bc(1,1,1),qmin,nh_u,dr,u,rho) !changed
      call set_bc(n(1),n(2),n(3),cbc(1,1,2),1,1,.true. ,bc(1,1,2),qmin,nh_u,dr,v)
      call set_bc(n(1),n(2),n(3),cbc(1,1,3),1,1,.true. ,bc(1,1,3),qmin,nh_u,dr,w)
    endif
    !
    ! along y
    !
    do q = 0,qmin
      dr(q) = dl(2)
    enddo
    if(front.eq.MPI_PROC_NULL) then ! y - bottom boundary
      call set_bc(n(1),n(2),n(3),cbc(0,2,1),0,2,.true. ,bc(0,2,1),qmin,nh_u,dr,u)
      call set_bc_mom(n(1),n(2),n(3),cbc(0,2,2),0,2,.false.,bc(0,2,2),qmin,nh_u,dr,v,rho) !changed
      call set_bc(n(1),n(2),n(3),cbc(0,2,3),0,2,.true. ,bc(0,2,3),qmin,nh_u,dr,w)
    endif
    if(back .eq.MPI_PROC_NULL) then ! y - top boundary
      call set_bc(n(1),n(2),n(3),cbc(1,2,1),1,2,.true. ,bc(1,2,1),qmin,nh_u,dr,u)
      call set_bc_mom(n(1),n(2),n(3),cbc(1,2,2),1,2,.false.,bc(1,2,2),qmin,nh_u,dr,v,rho) !changed
      call set_bc(n(1),n(2),n(3),cbc(1,2,3),1,2,.true. ,bc(1,2,3),qmin,nh_u,dr,w)
    endif
    !
    ! along z
    !
    do q = 0,qmin
      dr(q) = dzc(-q)
    enddo
    if(bottom.eq.MPI_PROC_NULL) then ! z - bottom boundary
      call set_bc(n(1),n(2),n(3),cbc(0,3,1),0,3,.true. ,bc(0,3,1),qmin,nh_u,dr,u)
      call set_bc(n(1),n(2),n(3),cbc(0,3,2),0,3,.true. ,bc(0,3,2),qmin,nh_u,dr,v)
    endif
    do q = 0,qmin
      dr(q) = dzf(-q)
    enddo
    if(bottom.eq.MPI_PROC_NULL) then ! z - bottom boundary
      call set_bc_mom(n(1),n(2),n(3),cbc(0,3,3),0,3,.false.,bc(0,3,3),qmin,nh_u,dr,w,rho) !changed
    endif
    !
    do q = 0,qmin 
      dr(q) = dzc(n(3)+q)
    enddo
    if(top.eq.MPI_PROC_NULL) then ! z - top boundary 
      call set_bc(n(1),n(2),n(3),cbc(1,3,1),1,3,.true. ,bc(1,3,1),qmin,nh_u,dr,u)
      call set_bc(n(1),n(2),n(3),cbc(1,3,2),1,3,.true. ,bc(1,3,2),qmin,nh_u,dr,v)
    endif
    do q = 0,qmin 
      dr(q) = dzf(n(3)+q)
    enddo
    if(top.eq.MPI_PROC_NULL) then ! z - top boundary
      call set_bc_mom(n(1),n(2),n(3),cbc(1,3,3),1,3,.false.,bc(1,3,3),qmin,nh_u,dr,w,rho) !changed
    endif
    !
    if(sem_switch) call SEM_apply(n,nh_u,u,v,w,.TRUE.)
    !
    call updthalo(n(1),n(2),n(3),nh_u,halo(ind1),ind1,u)
    call updthalo(n(1),n(2),n(3),nh_u,halo(ind2),ind2,u)
    !
    call updthalo(n(1),n(2),n(3),nh_u,halo(ind1),ind1,v)
    call updthalo(n(1),n(2),n(3),nh_u,halo(ind2),ind2,v)
    !
    call updthalo(n(1),n(2),n(3),nh_u,halo(ind1),ind1,w)
    call updthalo(n(1),n(2),n(3),nh_u,halo(ind2),ind2,w)
    !
    ! outflow bc.
    !
    do q = 1,3
      do idir = 0,1
        if(isoutflow(idir,q)) then
          if(idir.eq.0) sgn = -1
          if(idir.eq.1) sgn = +1
          ioutflowdir = q*sgn
          call outflow_mom(n(1),n(2),n(3),ioutflowdir,nh_d,nh_u,dx,dy,dz,dzf,u,v,w,rho,mfx,mfy,mfz,nse_switch,reg_switch)
        endif
      enddo
    enddo
    !
    return
  end subroutine bounduvw_mom
  !
  subroutine bounduvw(cbc,n,bc,nh_d,nh_u,halo,isoutflow,dl,dzc,dzf,u,v,w,rho,mfx,mfy,mfz,nse_switch,reg_switch,sem_switch)
    !
    use mod_common_mpi, only: MPI_PROC_NULL,left,right,front,back,top,bottom,myid
    use mod_sem       , only: SEM_apply
    !
    ! imposes velocity boundary conditions 
    ! on an arbitrary stecil size
    !
    implicit none
    !
    character(len=1), intent(in   ), dimension(0:1,3,3)                 :: cbc
    integer         , intent(in   ), dimension(3)                       :: n 
    real(rp)        , intent(in   ), dimension(0:1,3,3)                 :: bc
    integer         , intent(in   )                                     :: nh_d,nh_u
    integer         , intent(in   ), dimension(3)                       :: halo
    logical         , intent(in   ), dimension(0:1,3)                   :: isoutflow
    logical         , intent(in   )                                     :: sem_switch
    real(rp)        , intent(in   ), dimension(3)                       :: dl
    real(rp)        , intent(in   ), dimension(1-nh_d:)                 :: dzc,dzf
    real(rp)        , intent(inout), dimension(1-nh_u:,1-nh_u:,1-nh_u:) :: u,v,w
    real(rp)        , intent(in   ), dimension(1-nh_u:,1-nh_u:,1-nh_u:) :: mfx,mfy,mfz
    real(rp)        , intent(in   )                                     :: nse_switch,reg_switch
    real(rp)        , intent(in   ), dimension(0:,0:,0:)                :: rho
    !
    real(rp), dimension(0:nh_u-1) :: dr
    real(rp) :: dx,dy,dz
    integer  :: q,idir,sgn,ioutflowdir,qmin
    integer  :: ind1,ind2
    !@cuf attributes(managed) :: u, v, w, dzc, dzf
    !
    dx = dl(1)
    dy = dl(2)
    dz = dl(3)
    !
    qmin = abs(1-nh_u)
    !
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
    !call updthalo(n(1),n(2),n(3),nh_u,halo(ind1),ind1,u)
    !call updthalo(n(1),n(2),n(3),nh_u,halo(ind2),ind2,u)
    !
    !call updthalo(n(1),n(2),n(3),nh_u,halo(ind1),ind1,v)
    !call updthalo(n(1),n(2),n(3),nh_u,halo(ind2),ind2,v)
    !
    !call updthalo(n(1),n(2),n(3),nh_u,halo(ind1),ind1,w)
    !call updthalo(n(1),n(2),n(3),nh_u,halo(ind2),ind2,w)
    !
    ! along x
    !
    do q = 0,qmin
      dr(q) = dl(1)
    enddo
    if(left .eq.MPI_PROC_NULL) then ! x - bottom boundary
      call set_bc(n(1),n(2),n(3),cbc(0,1,1),0,1,.false.,bc(0,1,1),qmin,nh_u,dr,u)
      call set_bc(n(1),n(2),n(3),cbc(0,1,2),0,1,.true. ,bc(0,1,2),qmin,nh_u,dr,v)
      call set_bc(n(1),n(2),n(3),cbc(0,1,3),0,1,.true. ,bc(0,1,3),qmin,nh_u,dr,w)
    endif
    if(right.eq.MPI_PROC_NULL) then ! x - top boundary
      call set_bc(n(1),n(2),n(3),cbc(1,1,1),1,1,.false.,bc(1,1,1),qmin,nh_u,dr,u)
      call set_bc(n(1),n(2),n(3),cbc(1,1,2),1,1,.true. ,bc(1,1,2),qmin,nh_u,dr,v)
      call set_bc(n(1),n(2),n(3),cbc(1,1,3),1,1,.true. ,bc(1,1,3),qmin,nh_u,dr,w)
    endif
    !
    ! along y
    !
    do q = 0,qmin
      dr(q) = dl(2)
    enddo
    if(front.eq.MPI_PROC_NULL) then ! y - bottom boundary
      call set_bc(n(1),n(2),n(3),cbc(0,2,1),0,2,.true. ,bc(0,2,1),qmin,nh_u,dr,u)
      call set_bc(n(1),n(2),n(3),cbc(0,2,2),0,2,.false.,bc(0,2,2),qmin,nh_u,dr,v)
      call set_bc(n(1),n(2),n(3),cbc(0,2,3),0,2,.true. ,bc(0,2,3),qmin,nh_u,dr,w)
    endif
    if(back .eq.MPI_PROC_NULL) then ! y - top boundary
      call set_bc(n(1),n(2),n(3),cbc(1,2,1),1,2,.true. ,bc(1,2,1),qmin,nh_u,dr,u)
      call set_bc(n(1),n(2),n(3),cbc(1,2,2),1,2,.false.,bc(1,2,2),qmin,nh_u,dr,v)
      call set_bc(n(1),n(2),n(3),cbc(1,2,3),1,2,.true. ,bc(1,2,3),qmin,nh_u,dr,w)
    endif
    !
    ! along z
    !
    do q = 0,qmin
      dr(q) = dzc(-q)
    enddo
    if(bottom.eq.MPI_PROC_NULL) then ! z - bottom boundary
      call set_bc(n(1),n(2),n(3),cbc(0,3,1),0,3,.true. ,bc(0,3,1),qmin,nh_u,dr,u)
      call set_bc(n(1),n(2),n(3),cbc(0,3,2),0,3,.true. ,bc(0,3,2),qmin,nh_u,dr,v)
    endif
    do q = 0,qmin
      dr(q) = dzf(-q)
    enddo
    if(bottom.eq.MPI_PROC_NULL) then ! z - bottom boundary
      call set_bc(n(1),n(2),n(3),cbc(0,3,3),0,3,.false.,bc(0,3,3),qmin,nh_u,dr,w)
    endif
    !
    do q = 0,qmin 
      dr(q) = dzc(n(3)+q)
    enddo
    if(top.eq.MPI_PROC_NULL) then ! z - top boundary 
      call set_bc(n(1),n(2),n(3),cbc(1,3,1),1,3,.true. ,bc(1,3,1),qmin,nh_u,dr,u)
      call set_bc(n(1),n(2),n(3),cbc(1,3,2),1,3,.true. ,bc(1,3,2),qmin,nh_u,dr,v)
    endif
    do q = 0,qmin 
      dr(q) = dzf(n(3)+q)
    enddo
    if(top.eq.MPI_PROC_NULL) then ! z - top boundary
      call set_bc(n(1),n(2),n(3),cbc(1,3,3),1,3,.false.,bc(1,3,3),qmin,nh_u,dr,w)
    endif
    !
    if(sem_switch) call SEM_apply(n,nh_u,u,v,w,.FALSE.)
    !
    call updthalo(n(1),n(2),n(3),nh_u,halo(ind1),ind1,u)
    call updthalo(n(1),n(2),n(3),nh_u,halo(ind2),ind2,u)
    !
    call updthalo(n(1),n(2),n(3),nh_u,halo(ind1),ind1,v)
    call updthalo(n(1),n(2),n(3),nh_u,halo(ind2),ind2,v)
    !
    call updthalo(n(1),n(2),n(3),nh_u,halo(ind1),ind1,w)
    call updthalo(n(1),n(2),n(3),nh_u,halo(ind2),ind2,w)
    !
    ! outflow bc.
    !
    do q = 1,3
      do idir = 0,1
        if(isoutflow(idir,q)) then
          if(idir.eq.0) sgn = -1
          if(idir.eq.1) sgn = +1
          ioutflowdir = q*sgn
          call outflow(n(1),n(2),n(3),ioutflowdir,nh_d,nh_u,dx,dy,dz,dzf,u,v,w,rho,mfx,mfy,mfz,nse_switch,reg_switch)
        endif
      enddo
    enddo
    !
    return
  end subroutine bounduvw
  !
  subroutine boundp(cbc,n,bc,nh_d,nh_p,halo,dl,dzc,dzf,p)
    !
    use mod_common_mpi, only: MPI_PROC_NULL,left,right,front,back,top,bottom
    !
    ! imposes pressure and scalar boundary conditions
    ! on an arbitrary stencil size
    !
    implicit none
    !
    character(len=1), intent(in   ), dimension(0:1,3)                   :: cbc
    integer         , intent(in   ), dimension(3)                       :: n
    real(rp)        , intent(in   ), dimension(0:1,3)                   :: bc
    integer         , intent(in   )                                     :: nh_d,nh_p 
    integer         , intent(in   ), dimension(3)                       :: halo 
    real(rp)        , intent(in   ), dimension(3)                       :: dl
    real(rp)        , intent(in   ), dimension(1-nh_d:)                 :: dzc,dzf
    real(rp)        , intent(inout), dimension(1-nh_p:,1-nh_p:,1-nh_p:) :: p
    !
    real(rp), dimension(0:nh_p-1) :: dr
    integer :: q,qmin
    integer :: ind1,ind2
    !@cuf attributes(managed) :: p, dzc, dzf
    !
    qmin = abs(1-nh_p)
    !
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
    call updthalo(n(1),n(2),n(3),nh_p,halo(ind1),ind1,p)
    call updthalo(n(1),n(2),n(3),nh_p,halo(ind2),ind2,p)
    !
    ! along x
    !
    do q = 0,qmin
      dr(q) = dl(1)
    enddo
    if(left .eq.MPI_PROC_NULL) then ! x - bottom part
      call set_bc(n(1),n(2),n(3),cbc(0,1),0,1,.true.,bc(0,1),qmin,nh_p,dr,p)
    endif
    if(right.eq.MPI_PROC_NULL) then ! x - upper part
      call set_bc(n(1),n(2),n(3),cbc(1,1),1,1,.true.,bc(1,1),qmin,nh_p,dr,p)
    endif
    !
    ! along y
    !
    do q = 0,qmin
      dr(q) = dl(2)
    enddo
    if(front.eq.MPI_PROC_NULL) then ! y - bottom part
      call set_bc(n(1),n(2),n(3),cbc(0,2),0,2,.true.,bc(0,2),qmin,nh_p,dr,p)
     endif
    if(back .eq.MPI_PROC_NULL) then ! y - upper part
      call set_bc(n(1),n(2),n(3),cbc(1,2),1,2,.true.,bc(1,2),qmin,nh_p,dr,p)
    endif
    !
    ! along z
    !
    do q=0,qmin 
      dr(q) = dzc(-q)
    enddo
    if(bottom.eq.MPI_PROC_NULL) then ! z - bottom boundary
      call set_bc(n(1),n(2),n(3),cbc(0,3),0,3,.true.,bc(0,3),qmin,nh_p,dr,p)
    endif
    do q=0,qmin 
      dr(q) = dzc(n(3)+q)
    enddo
    if(top.eq.MPI_PROC_NULL) then ! z - top boundary
      call set_bc(n(1),n(2),n(3),cbc(1,3),1,3,.true.,bc(1,3),qmin,nh_p,dr,p)
    endif
    !
    return
  end subroutine boundp



  subroutine boundp_IBM(cbc,n,bc,nh_d,nh_p,halo,dl,dzc,dzf,p,IBM_h)
    !
    use mod_common_mpi, only: MPI_PROC_NULL,left,right,front,back,top,bottom
    !
    ! imposes pressure and scalar boundary conditions
    ! on an arbitrary stencil size
    !
    implicit none
    !
    character(len=1), intent(in   ), dimension(0:1,3)                   :: cbc
    integer         , intent(in   ), dimension(3)                       :: n
    real(rp)        , intent(in   ), dimension(0:1,3)                   :: bc
    integer         , intent(in   )                                     :: nh_d,nh_p,IBM_h
    integer         , intent(in   ), dimension(3)                       :: halo 
    real(rp)        , intent(in   ), dimension(3)                       :: dl
    real(rp)        , intent(in   ), dimension(1-nh_d:)                 :: dzc,dzf
    real(rp)        , intent(inout), dimension(1-nh_p:,1-nh_p:,1-nh_p:) :: p
    !
    real(rp), dimension(0:nh_p-1) :: dr
    real(rp), dimension(0:nh_p-1+IBM_h) :: dr_IBM
    integer :: q,qmin
    integer :: ind1,ind2
    !@cuf attributes(managed) :: p, dzc, dzf
    !
    qmin = abs(1-nh_p)
    !
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
    call updthalo(n(1),n(2),n(3),nh_p,halo(ind1),ind1,p)
    call updthalo(n(1),n(2),n(3),nh_p,halo(ind2),ind2,p)
    !
    ! along x
    !
    do q = 0,qmin
      dr(q) = dl(1)
    enddo
    if(left .eq.MPI_PROC_NULL) then ! x - bottom part
      call set_bc_IBM(n(1),n(2),n(3),cbc(0,1),0,1,.true.,bc(0,1),qmin,nh_p,dr,p,IBM_h)
    endif
    if(right.eq.MPI_PROC_NULL) then ! x - upper part
      call set_bc_IBM(n(1),n(2),n(3),cbc(1,1),1,1,.true.,bc(1,1),qmin,nh_p,dr,p,IBM_h)
    endif
    !
    ! along y
    !
    do q = 0,qmin
      dr(q) = dl(2)
    enddo
    if(front.eq.MPI_PROC_NULL) then ! y - bottom part
      call set_bc_IBM(n(1),n(2),n(3),cbc(0,2),0,2,.true.,bc(0,2),qmin,nh_p,dr,p,IBM_h)
     endif
    if(back .eq.MPI_PROC_NULL) then ! y - upper part
      call set_bc_IBM(n(1),n(2),n(3),cbc(1,2),1,2,.true.,bc(1,2),qmin,nh_p,dr,p,IBM_h)
    endif
    !
    ! along z
    !
    qmin = abs(1-nh_p-IBM_h)
    do q=0,qmin 
      dr_IBM(q) = dzc(IBM_h-q)
    enddo
    if(bottom.eq.MPI_PROC_NULL) then ! z - bottom boundary ---> this is changed for IBM
      call set_bc_IBM(n(1),n(2),n(3),cbc(0,3),0,3,.true.,bc(0,3),qmin,nh_p,dr_IBM,p,IBM_h)
    endif
    qmin = abs(1-nh_p)
    do q=0,qmin 
      dr(q) = dzc(n(3)+q)
    enddo
    if(top.eq.MPI_PROC_NULL) then ! z - top boundary
      call set_bc_IBM(n(1),n(2),n(3),cbc(1,3),1,3,.true.,bc(1,3),qmin,nh_p,dr,p,IBM_h)
    endif
    !
    return
  end subroutine boundp_IBM
  !
  subroutine set_bc_mom(nx,ny,nz,ctype,ibound,idir,centered,rvalue,qq_d,nh_p,dr,p,rho)
    ! face centered Dirichlet BC values are multiplied with density
    ! note: to be generalized with non-homogeneous Neumann
    !
    implicit none
    !
    integer         , intent(in   )                                     :: nx,ny,nz
    character(len=1), intent(in   )                                     :: ctype
    integer         , intent(in   )                                     :: ibound,idir
    logical         , intent(in   )                                     :: centered
    real(rp)        , intent(in   )                                     :: rvalue
    integer         , intent(in   )                                     :: qq_d,nh_p
    real(rp)        , intent(in   ), dimension(0:qq_d)                  :: dr
    real(rp)        , intent(inout), dimension(1-nh_p:,1-nh_p:,1-nh_p:) :: p
    real(rp)        , intent(in   ), dimension(0:,0:,0:)                :: rho
    !
    real(rp), dimension(0:qq_d) :: factor
    real(rp) :: factor_value, sgn, rho_f
    integer  :: i,j,k,q
    !@cuf attributes(managed) :: p
    !
    do q=0,qq_d
      factor(q) = rvalue
    enddo
    sgn = 0._rp
    !
    if(ctype.eq.'D'.and.centered) then
      do q=0,qq_d
        factor(q) = +2._rp*factor(q)
      enddo
      sgn = -1._rp
    endif
    if(ctype.eq.'N') then
      if(    ibound.eq.0) then
        do q=0,qq_d
          factor(q) = -dr(q)*factor(q)
        enddo
      elseif(ibound.eq.1) then
        do q=0,qq_d
          factor(q) =  dr(q)*factor(q)
        enddo
      endif
      sgn = +1._rp
    endif
    !
    select case(ctype)
    case('P')
      select case(idir)
      case(1)
        !$OMP PARALLEL DO DEFAULT(none) &
        !$OMP PRIVATE(j,k,q,factor,sgn) &
        !$OMP SHARED(ny,nz,nh_p,p)
        !$acc parallel loop collapse(3)
        do k=1-nh_p,nz+nh_p
          do j=1-nh_p,ny+nh_p
            do q=0,nh_p-1
              p(0-q   ,j,k) = p(nx-q,j,k)
              p(nx+1+q,j,k) = p(1+q ,j,k)
            enddo
          enddo
        enddo
        !$acc end parallel loop
        !$OMP END PARALLEL DO
      case(2)
        !$OMP PARALLEL DO DEFAULT(none) &
        !$OMP PRIVATE(i,k,q,factor,sgn) &
        !$OMP SHARED(nx,nz,nh_p,p)
        !$acc parallel loop collapse(3)
        do k=1-nh_p,nz+nh_p
          do q=0,nh_p-1
            do i=1-nh_p,nx+nh_p
              p(i,0-q   ,k) = p(i,ny-q,k)
              p(i,ny+1+q,k) = p(i,1+q ,k)
            enddo
          enddo
        enddo
        !$acc end parallel loop
        !$OMP END PARALLEL DO
      case(3)
        !$OMP PARALLEL DO DEFAULT(none) &
        !$OMP PRIVATE(i,j,q,factor,sgn) &
        !$OMP SHARED(nx,ny,nz,nh_p,p)
        !$acc parallel loop collapse(3)
        do q=0,nh_p-1
          do j=1-nh_p,ny+nh_p
            do i=1-nh_p,nx+nh_p
              p(i,j,0-q   ) = p(i,j,nz-q)
              p(i,j,nz+1+q) = p(i,j, 1+q)
            enddo
          enddo
        enddo
        !$acc end parallel loop
        !$OMP END PARALLEL DO
      end select
    case('D','N')
      if(centered) then
        select case(idir)
        case(1)
          if    (ibound.eq.0) then
            !$OMP PARALLEL DO DEFAULT(none) &
            !$OMP PRIVATE(j,k,q,factor,sgn) &
            !$OMP SHARED(ny,nz,nh_p,p)
            do q=0,nh_p-1
              factor_value = factor(q)
              !
              !$acc parallel loop collapse(2)
              do k=1-nh_p,nz+nh_p
                do j=1-nh_p,ny+nh_p
                  p(0-q,j,k) = factor_value+sgn*p(1+q,j,k)
                enddo
              enddo
              !$acc end parallel loop
            enddo
            !$OMP END PARALLEL DO
          elseif(ibound.eq.1) then
            !$OMP PARALLEL DO DEFAULT(none) &
            !$OMP PRIVATE(j,k,q,factor,sgn) &
            !$OMP SHARED(nx,ny,nz,nh_p,p)
            do q=0,nh_p-1
              factor_value = factor(q)
              !
              !$acc parallel loop collapse(2)
              do k=1-nh_p,nz+nh_p
                do j=1-nh_p,ny+nh_p
                  p(nx+1+q,j,k) = factor_value+sgn*p(nx-q,j,k)
                enddo
              enddo
              !$acc end parallel loop
            enddo
            !$OMP END PARALLEL DO
          endif
        case(2)
          if    (ibound.eq.0) then
            !$OMP PARALLEL DO DEFAULT(none) &
            !$OMP PRIVATE(i,k,q,factor,sgn) &
            !$OMP SHARED(nx,nz,nh_p,p)
            do q=0,nh_p-1
              factor_value = factor(q)
              !
              !$acc parallel loop collapse(2)
              do k=1-nh_p,nz+nh_p
                do i=1-nh_p,nx+nh_p
                  p(i,0-q  ,k) = factor_value+sgn*p(i,1+q,k)
                enddo
              enddo
              !$acc end parallel loop
            enddo
            !$OMP END PARALLEL DO
          elseif(ibound.eq.1) then
            !$OMP PARALLEL DO DEFAULT(none) &
            !$OMP PRIVATE(i,k,q,factor,sgn) &
            !$OMP SHARED(nx,ny,nz,nh_p,p)
            do q=0,nh_p-1
              factor_value = factor(q)
              !
              !$acc parallel loop collapse(2)
              do k=1-nh_p,nz+nh_p
                do i=1-nh_p,nx+nh_p
                  p(i,ny+1+q,k) = factor_value+sgn*p(i,ny-q,k)
                enddo
              enddo 
              !$acc end parallel loop
            enddo
            !$OMP END PARALLEL DO
          endif
        case(3)
          if    (ibound.eq.0) then
            !$OMP PARALLEL DO DEFAULT(none) &
            !$OMP PRIVATE(i,j,q,factor,sgn) &
            !$OMP SHARED(nx,ny,nh_p,p)
            do q=0,nh_p-1
              factor_value = factor(q)
              !
              !$acc parallel loop collapse(2)
              do j=1-nh_p,ny+nh_p
                do i=1-nh_p,nx+nh_p
                  p(i,j,0-q ) = factor_value+sgn*p(i,j,1+q)
                enddo 
              enddo
              !$acc end parallel loop
            enddo
            !$OMP END PARALLEL DO
          elseif(ibound.eq.1) then
            !$OMP PARALLEL DO DEFAULT(none) &
            !$OMP PRIVATE(i,j,q,factor,sgn) &
            !$OMP SHARED(nx,ny,nz,nh_p,p)
            do q=0,nh_p-1
              factor_value = factor(q)
              !
              !$acc parallel loop collapse(2)
              do j=1-nh_p,ny+nh_p
                do i=1-nh_p,nx+nh_p
                  p(i,j,nz+1+q) = factor_value+sgn*p(i,j,nz-q)
                enddo
              enddo
              !$acc end parallel loop
            enddo
            !$OMP END PARALLEL DO
          endif
        end select
      elseif(.not.centered.and.ctype.eq.'D') then   !! This is the part that is changed here!
        select case(idir)
        case(1)
          if    (ibound.eq.0) then
            !$OMP PARALLEL DO DEFAULT(none) &
            !$OMP PRIVATE(j,k,q,factor,sgn) &
            !$OMP SHARED(ny,nz,nh_p,p)
            do q=0,nh_p-1
              factor_value = factor(q)
              !
              !$acc parallel loop collapse(2)
              do k=1-nh_p,nz+nh_p
                do j=1-nh_p,ny+nh_p
                  rho_f      = 0.5*(rho(1,j,k)+rho(0,j,k))
                  p(0-q,j,k) = factor_value*rho_f
                enddo
              enddo
              !$acc end parallel loop
            enddo
            !$OMP END PARALLEL DO
          elseif(ibound.eq.1) then
            !$OMP PARALLEL DO DEFAULT(none) &
            !$OMP PRIVATE(j,k,q,factor,sgn) &
            !$OMP SHARED(nx,ny,nz,nh_p,p)
            do q=0,nh_p-1
              factor_value = factor(q)
              !
              !$acc parallel loop collapse(2)
              do k=1-nh_p,nz+nh_p
                do j=1-nh_p,ny+nh_p
                  rho_f         = 0.5*(rho(nx,j,k)+rho(nx+1,j,k))
                  p(nx+q  ,j,k) = factor_value*rho_f
                  p(nx+1+q,j,k) = p(nx-1-q,j,k)
                enddo
              enddo
              !$acc end parallel loop
            enddo
            !$OMP END PARALLEL DO
          endif
        case(2)
          if    (ibound.eq.0) then
            !$OMP PARALLEL DO DEFAULT(none) &
            !$OMP PRIVATE(i,k,q,factor,sgn) &
            !$OMP SHARED(nx,nz,nh_p,p)
            do q=0,nh_p-1
              factor_value = factor(q)
              !
              !$acc parallel loop collapse(2)
              do k=1-nh_p,nz+nh_p
                do i=1-nh_p,nx+nh_p
                  rho_f      = 0.5*(rho(i,0,k)+rho(i,1,k))
                  p(i,0-q,k) = factor_value*rho_f
                enddo
              enddo
              !$acc end parallel loop
            enddo
            !$OMP END PARALLEL DO
          elseif(ibound.eq.1) then
            !$OMP PARALLEL DO DEFAULT(none) &
            !$OMP PRIVATE(i,j,q,factor,sgn) &
            !$OMP SHARED(nx,ny,nz,nh_p,p)
            do q=0,nh_p-1
              factor_value = factor(q)
              !
              !$acc parallel loop collapse(2)
              do k=1-nh_p,nz+nh_p
                do i=1-nh_p,nx+nh_p
                  rho_f         = 0.5*(rho(i,ny,k)+rho(i,ny+1,k))
                  p(i,ny+q  ,k) = factor_value*rho_f
                  p(i,ny+1+q,k) = p(i,ny-1-q,k)
                enddo
              enddo
              !$acc end parallel loop
            enddo
            !$OMP END PARALLEL DO
          endif
        case(3)
          if    (ibound.eq.0) then
            !$OMP PARALLEL DO DEFAULT(none) &
            !$OMP PRIVATE(i,j,q,factor,sgn) &
            !$OMP SHARED(nx,ny,nh_p,p)
            do q=0,nh_p-1
              factor_value = factor(q)
              !
              !$acc parallel loop collapse(2)
              do j=1-nh_p,ny+nh_p
                do i=1-nh_p,nx+nh_p
                  rho_f      = 0.5*(rho(i,j,0)+rho(i,j,1))
                  p(i,j,0-q) = factor_value*rho_f
                enddo
              enddo
              !$acc end parallel loop
            enddo
            !$OMP END PARALLEL DO
          elseif(ibound.eq.1) then
            !$OMP PARALLEL DO DEFAULT(none) &
            !$OMP PRIVATE(i,j,q,factor,sgn) &
            !$OMP SHARED(nx,ny,nz,nh_p,p)
            do q=0,nh_p-1
              factor_value = factor(q)
              !
              !$acc parallel loop collapse(2)
              do j=1-nh_p,ny+nh_p
                do i=1-nh_p,nx+nh_p
                  rho_f         = 0.5*(rho(i,j,nz)+rho(i,j,nz+1))
                  p(i,j,nz+q  ) = factor_value*rho_f
                  p(i,j,nz+1+q) = p(i,j,nz-1-q)
                enddo
              enddo
              !$acc end parallel loop
            enddo
            !$OMP END PARALLEL DO
          endif
        end select
      elseif(.not.centered.and.ctype.eq.'N') then
        select case(idir)
        case(1)
          if    (ibound.eq.0) then
            !$OMP PARALLEL DO DEFAULT(none) &
            !$OMP PRIVATE(j,k,q,factor,sgn) &
            !$OMP SHARED(ny,nz,nh_p,p)
            do q=0,nh_p-1
              factor_value = factor(q)
              !
              !$acc parallel loop collapse(2)
              do k=1-nh_p,nz+nh_p
                do j=1-nh_p,ny+nh_p
                  p(0-q,j,k) = 1.0_rp*factor_value + p(1+q,j,k)
                enddo
              enddo
              !$acc end parallel loop
            enddo
            !$OMP END PARALLEL DO
          elseif(ibound.eq.1) then
            !$OMP PARALLEL DO DEFAULT(none) &
            !$OMP PRIVATE(j,k,q,factor,sgn) &
            !$OMP SHARED(nx,ny,nz,nh_p,p)
            do q=0,nh_p-1
              factor_value = factor(q)
              !
              !$acc parallel loop collapse(2)
              do k=1-nh_p,nz+nh_p
                do j=1-nh_p,ny+nh_p
                  p(nx+q  ,j,k) = 1.0_rp*factor_value + p(nx-1-q,j,k)
                  p(nx+1+q,j,k) = 2.0_rp*factor_value + p(nx-1-q,j,k)
                enddo
              enddo
              !$acc end parallel loop
            enddo
            !$OMP END PARALLEL DO
          endif
        case(2)
          if    (ibound.eq.0) then
            !$OMP PARALLEL DO DEFAULT(none) &
            !$OMP PRIVATE(i,k,q,factor,sgn) &
            !$OMP SHARED(nx,nz,nh_p,p)
            do q=0,nh_p-1
              factor_value = factor(q)
              !
              !$acc parallel loop collapse(2)
              do k=1-nh_p,nz+nh_p
                do i=1-nh_p,nx+nh_p
                  p(i,0-q  ,k) = 1.0_rp*factor_value + p(i,1+q  ,k) 
                enddo
              enddo
              !$acc end parallel loop
            enddo
            !$OMP END PARALLEL DO
          elseif(ibound.eq.1) then
            !$OMP PARALLEL DO DEFAULT(none) &
            !$OMP PRIVATE(i,k,q,factor,sgn) &
            !$OMP SHARED(nx,ny,nz,nh_p,p)
            do q=0,nh_p-1
              factor_value = factor(q)
              !
              !$acc parallel loop collapse(2)
              do k=1-nh_p,nz+nh_p
                do i=1-nh_p,nx+nh_p
                  p(i,ny+q  ,k) = 1.0_rp*factor_value + p(i,ny-1-q,k)
                  p(i,ny+1+q,k) = 2.0_rp*factor_value + p(i,ny-1-q,k)
                enddo
              enddo
              !$acc end parallel loop
            enddo
            !$OMP END PARALLEL DO
          endif
        case(3)
          if    (ibound.eq.0) then
            !$OMP PARALLEL DO DEFAULT(none) &
            !$OMP PRIVATE(i,j,q,factor,sgn) &
            !$OMP SHARED(nx,ny,nh_p,p)
            do q=0,nh_p-1
              factor_value = factor(q)
              !
              !$acc parallel loop collapse(2)
              do j=1-nh_p,ny+nh_p
                do i=1-nh_p,nx+nh_p
                  p(i,j,0-q  ) = 1.0_rp*factor_value + p(i,j,1+q  )
                enddo
              enddo
              !$acc end parallel loop
            enddo
            !$OMP END PARALLEL DO
          elseif(ibound.eq.1) then
            !$OMP PARALLEL DO DEFAULT(none) &
            !$OMP PRIVATE(i,j,q,factor,sgn) &
            !$OMP SHARED(nx,ny,nz,nh_p,p)
            do q=0,nh_p-1
              factor_value = factor(q)
              !
              !$acc parallel loop collapse(2)
              do j=1-nh_p,ny+nh_p
                do i=1-nh_p,nx+nh_p
                  p(i,j,nz+q  ) = 1.0_rp*factor_value + p(i,j,nz-1-q)
                  p(i,j,nz+1+q) = 2.0_rp*factor_value + p(i,j,nz-1-q)
                enddo
              enddo
              !$acc end parallel loop
            enddo
            !$OMP END PARALLEL DO
          endif
        end select
      endif
    end select
    !
    return
  end subroutine set_bc_mom
  !
  subroutine set_bc(nx,ny,nz,ctype,ibound,idir,centered,rvalue,qq_d,nh_p,dr,p)
    !
    ! note: to be generalized with non-homogeneous Neumann
    !
    implicit none
    !
    integer         , intent(in   )                                     :: nx,ny,nz
    character(len=1), intent(in   )                                     :: ctype
    integer         , intent(in   )                                     :: ibound,idir
    logical         , intent(in   )                                     :: centered
    real(rp)        , intent(in   )                                     :: rvalue
    integer         , intent(in   )                                     :: qq_d,nh_p
    real(rp)        , intent(in   ), dimension(0:qq_d)                  :: dr
    real(rp)        , intent(inout), dimension(1-nh_p:,1-nh_p:,1-nh_p:) :: p
    !
    real(rp), dimension(0:qq_d) :: factor
    real(rp) :: factor_value, sgn
    integer  :: i,j,k,q
    !@cuf attributes(managed) :: p
    !
    do q=0,qq_d
      factor(q) = rvalue
    enddo
    sgn = 0._rp
    !
    if(ctype.eq.'D'.and.centered) then
      do q=0,qq_d
        factor(q) = +2._rp*factor(q)
      enddo
      sgn = -1._rp
    endif
    if(ctype.eq.'N') then
      if(    ibound.eq.0) then
        do q=0,qq_d
          factor(q) = -dr(q)*factor(q)
        enddo
      elseif(ibound.eq.1) then
        do q=0,qq_d
          factor(q) =  dr(q)*factor(q)
        enddo
      endif
      sgn = +1._rp
    endif
    !
    select case(ctype)
    case('P')
      select case(idir)
      case(1)
        !$OMP PARALLEL DO DEFAULT(none) &
        !$OMP PRIVATE(j,k,q,factor,sgn) &
        !$OMP SHARED(ny,nz,nh_p,p)
        !$acc parallel loop collapse(3)
        do k=1-nh_p,nz+nh_p
          do j=1-nh_p,ny+nh_p
            do q=0,nh_p-1
              p(0-q   ,j,k) = p(nx-q,j,k)
              p(nx+1+q,j,k) = p(1+q ,j,k)
            enddo
          enddo
        enddo
        !$acc end parallel loop
        !$OMP END PARALLEL DO
      case(2)
        !$OMP PARALLEL DO DEFAULT(none) &
        !$OMP PRIVATE(i,k,q,factor,sgn) &
        !$OMP SHARED(nx,nz,nh_p,p)
        !$acc parallel loop collapse(3)
        do k=1-nh_p,nz+nh_p
          do q=0,nh_p-1
            do i=1-nh_p,nx+nh_p
              p(i,0-q   ,k) = p(i,ny-q,k)
              p(i,ny+1+q,k) = p(i,1+q ,k)
            enddo
          enddo
        enddo
        !$acc end parallel loop
        !$OMP END PARALLEL DO
      case(3)
        !$OMP PARALLEL DO DEFAULT(none) &
        !$OMP PRIVATE(i,j,q,factor,sgn) &
        !$OMP SHARED(nx,ny,nz,nh_p,p)
        !$acc parallel loop collapse(3)
        do q=0,nh_p-1
          do j=1-nh_p,ny+nh_p
            do i=1-nh_p,nx+nh_p
              p(i,j,0-q   ) = p(i,j,nz-q)
              p(i,j,nz+1+q) = p(i,j, 1+q)
            enddo
          enddo
        enddo
        !$acc end parallel loop
        !$OMP END PARALLEL DO
      end select
    case('D','N')
      if(centered) then
        select case(idir)
        case(1)
          if    (ibound.eq.0) then
            !$OMP PARALLEL DO DEFAULT(none) &
            !$OMP PRIVATE(j,k,q,factor,sgn) &
            !$OMP SHARED(ny,nz,nh_p,p)
            do q=0,nh_p-1
              factor_value = factor(q)
              !
              !$acc parallel loop collapse(2)
              do k=1-nh_p,nz+nh_p
                do j=1-nh_p,ny+nh_p
                  p(0-q,j,k) = factor_value+sgn*p(1+q,j,k)
                enddo
              enddo
              !$acc end parallel loop
            enddo
            !$OMP END PARALLEL DO
          elseif(ibound.eq.1) then
            !$OMP PARALLEL DO DEFAULT(none) &
            !$OMP PRIVATE(j,k,q,factor,sgn) &
            !$OMP SHARED(nx,ny,nz,nh_p,p)
            do q=0,nh_p-1
              factor_value = factor(q)
              !
              !$acc parallel loop collapse(2)
              do k=1-nh_p,nz+nh_p
                do j=1-nh_p,ny+nh_p
                  p(nx+1+q,j,k) = factor_value+sgn*p(nx-q,j,k)
                enddo
              enddo
              !$acc end parallel loop
            enddo
            !$OMP END PARALLEL DO
          endif
        case(2)
          if    (ibound.eq.0) then
            !$OMP PARALLEL DO DEFAULT(none) &
            !$OMP PRIVATE(i,k,q,factor,sgn) &
            !$OMP SHARED(nx,nz,nh_p,p)
            do q=0,nh_p-1
              factor_value = factor(q)
              !
              !$acc parallel loop collapse(2)
              do k=1-nh_p,nz+nh_p
                do i=1-nh_p,nx+nh_p
                  p(i,0-q  ,k) = factor_value+sgn*p(i,1+q,k)
                enddo
              enddo
              !$acc end parallel loop
            enddo
            !$OMP END PARALLEL DO
          elseif(ibound.eq.1) then
            !$OMP PARALLEL DO DEFAULT(none) &
            !$OMP PRIVATE(i,k,q,factor,sgn) &
            !$OMP SHARED(nx,ny,nz,nh_p,p)
            do q=0,nh_p-1
              factor_value = factor(q)
              !
              !$acc parallel loop collapse(2)
              do k=1-nh_p,nz+nh_p
                do i=1-nh_p,nx+nh_p
                  p(i,ny+1+q,k) = factor_value+sgn*p(i,ny-q,k)
                enddo
              enddo 
              !$acc end parallel loop
            enddo
            !$OMP END PARALLEL DO
          endif
        case(3)
          if    (ibound.eq.0) then
            !$OMP PARALLEL DO DEFAULT(none) &
            !$OMP PRIVATE(i,j,q,factor,sgn) &
            !$OMP SHARED(nx,ny,nh_p,p)
            do q=0,nh_p-1
              factor_value = factor(q)
              !
              !$acc parallel loop collapse(2)
              do j=1-nh_p,ny+nh_p
                do i=1-nh_p,nx+nh_p
                  p(i,j,0-q ) = factor_value+sgn*p(i,j,1+q)
                enddo 
              enddo
              !$acc end parallel loop
            enddo
            !$OMP END PARALLEL DO
          elseif(ibound.eq.1) then
            !$OMP PARALLEL DO DEFAULT(none) &
            !$OMP PRIVATE(i,j,q,factor,sgn) &
            !$OMP SHARED(nx,ny,nz,nh_p,p)
            do q=0,nh_p-1
              factor_value = factor(q)
              !
              !$acc parallel loop collapse(2)
              do j=1-nh_p,ny+nh_p
                do i=1-nh_p,nx+nh_p
                  p(i,j,nz+1+q) = factor_value+sgn*p(i,j,nz-q)
                enddo
              enddo
              !$acc end parallel loop
            enddo
            !$OMP END PARALLEL DO
          endif
        end select
      elseif(.not.centered.and.ctype.eq.'D') then
        select case(idir)
        case(1)
          if    (ibound.eq.0) then
            !$OMP PARALLEL DO DEFAULT(none) &
            !$OMP PRIVATE(j,k,q,factor,sgn) &
            !$OMP SHARED(ny,nz,nh_p,p)
            do q=0,nh_p-1
              factor_value = factor(q)
              !
              !$acc parallel loop collapse(2)
              do k=1-nh_p,nz+nh_p
                do j=1-nh_p,ny+nh_p
                  p(0-q,j,k) = factor_value
                enddo
              enddo
              !$acc end parallel loop
            enddo
            !$OMP END PARALLEL DO
          elseif(ibound.eq.1) then
            !$OMP PARALLEL DO DEFAULT(none) &
            !$OMP PRIVATE(j,k,q,factor,sgn) &
            !$OMP SHARED(nx,ny,nz,nh_p,p)
            do q=0,nh_p-1
              factor_value = factor(q)
              !
              !$acc parallel loop collapse(2)
              do k=1-nh_p,nz+nh_p
                do j=1-nh_p,ny+nh_p
                  p(nx+q  ,j,k) = factor_value
                  p(nx+1+q,j,k) = p(nx-1-q,j,k)
                enddo
              enddo
              !$acc end parallel loop
            enddo
            !$OMP END PARALLEL DO
          endif
        case(2)
          if    (ibound.eq.0) then
            !$OMP PARALLEL DO DEFAULT(none) &
            !$OMP PRIVATE(i,k,q,factor,sgn) &
            !$OMP SHARED(nx,nz,nh_p,p)
            do q=0,nh_p-1
              factor_value = factor(q)
              !
              !$acc parallel loop collapse(2)
              do k=1-nh_p,nz+nh_p
                do i=1-nh_p,nx+nh_p
                  p(i,0-q,k) = factor_value
                enddo
              enddo
              !$acc end parallel loop
            enddo
            !$OMP END PARALLEL DO
          elseif(ibound.eq.1) then
            !$OMP PARALLEL DO DEFAULT(none) &
            !$OMP PRIVATE(i,j,q,factor,sgn) &
            !$OMP SHARED(nx,ny,nz,nh_p,p)
            do q=0,nh_p-1
              factor_value = factor(q)
              !
              !$acc parallel loop collapse(2)
              do k=1-nh_p,nz+nh_p
                do i=1-nh_p,nx+nh_p
                  p(i,ny+q  ,k) = factor_value
                  p(i,ny+1+q,k) = p(i,ny-1-q,k)
                enddo
              enddo
              !$acc end parallel loop
            enddo
            !$OMP END PARALLEL DO
          endif
        case(3)
          if    (ibound.eq.0) then
            !$OMP PARALLEL DO DEFAULT(none) &
            !$OMP PRIVATE(i,j,q,factor,sgn) &
            !$OMP SHARED(nx,ny,nh_p,p)
            do q=0,nh_p-1
              factor_value = factor(q)
              !
              !$acc parallel loop collapse(2)
              do j=1-nh_p,ny+nh_p
                do i=1-nh_p,nx+nh_p
                  p(i,j,0-q) = factor_value 
                enddo
              enddo
              !$acc end parallel loop
            enddo
            !$OMP END PARALLEL DO
          elseif(ibound.eq.1) then
            !$OMP PARALLEL DO DEFAULT(none) &
            !$OMP PRIVATE(i,j,q,factor,sgn) &
            !$OMP SHARED(nx,ny,nz,nh_p,p)
            do q=0,nh_p-1
              factor_value = factor(q)
              !
              !$acc parallel loop collapse(2)
              do j=1-nh_p,ny+nh_p
                do i=1-nh_p,nx+nh_p
                  p(i,j,nz+q  ) = factor_value
                  p(i,j,nz+1+q) = p(i,j,nz-1-q)
                enddo
              enddo
              !$acc end parallel loop
            enddo
            !$OMP END PARALLEL DO
          endif
        end select
      elseif(.not.centered.and.ctype.eq.'N') then
        select case(idir)
        case(1)
          if    (ibound.eq.0) then
            !$OMP PARALLEL DO DEFAULT(none) &
            !$OMP PRIVATE(j,k,q,factor,sgn) &
            !$OMP SHARED(ny,nz,nh_p,p)
            do q=0,nh_p-1
              factor_value = factor(q)
              !
              !$acc parallel loop collapse(2)
              do k=1-nh_p,nz+nh_p
                do j=1-nh_p,ny+nh_p
                  p(0-q,j,k) = 1.0_rp*factor_value + p(1+q,j,k)
                enddo
              enddo
              !$acc end parallel loop
            enddo
            !$OMP END PARALLEL DO
          elseif(ibound.eq.1) then
            !$OMP PARALLEL DO DEFAULT(none) &
            !$OMP PRIVATE(j,k,q,factor,sgn) &
            !$OMP SHARED(nx,ny,nz,nh_p,p)
            do q=0,nh_p-1
              factor_value = factor(q)
              !
              !$acc parallel loop collapse(2)
              do k=1-nh_p,nz+nh_p
                do j=1-nh_p,ny+nh_p
                  p(nx+q  ,j,k) = 1.0_rp*factor_value + p(nx-1-q,j,k)
                  p(nx+1+q,j,k) = 2.0_rp*factor_value + p(nx-1-q,j,k)
                enddo
              enddo
              !$acc end parallel loop
            enddo
            !$OMP END PARALLEL DO
          endif
        case(2)
          if    (ibound.eq.0) then
            !$OMP PARALLEL DO DEFAULT(none) &
            !$OMP PRIVATE(i,k,q,factor,sgn) &
            !$OMP SHARED(nx,nz,nh_p,p)
            do q=0,nh_p-1
              factor_value = factor(q)
              !
              !$acc parallel loop collapse(2)
              do k=1-nh_p,nz+nh_p
                do i=1-nh_p,nx+nh_p
                  p(i,0-q  ,k) = 1.0_rp*factor_value + p(i,1+q  ,k) 
                enddo
              enddo
              !$acc end parallel loop
            enddo
            !$OMP END PARALLEL DO
          elseif(ibound.eq.1) then
            !$OMP PARALLEL DO DEFAULT(none) &
            !$OMP PRIVATE(i,k,q,factor,sgn) &
            !$OMP SHARED(nx,ny,nz,nh_p,p)
            do q=0,nh_p-1
              factor_value = factor(q)
              !
              !$acc parallel loop collapse(2)
              do k=1-nh_p,nz+nh_p
                do i=1-nh_p,nx+nh_p
                  p(i,ny+q  ,k) = 1.0_rp*factor_value + p(i,ny-1-q,k)
                  p(i,ny+1+q,k) = 2.0_rp*factor_value + p(i,ny-1-q,k)
                enddo
              enddo
              !$acc end parallel loop
            enddo
            !$OMP END PARALLEL DO
          endif
        case(3)
          if    (ibound.eq.0) then
            !$OMP PARALLEL DO DEFAULT(none) &
            !$OMP PRIVATE(i,j,q,factor,sgn) &
            !$OMP SHARED(nx,ny,nh_p,p)
            do q=0,nh_p-1
              factor_value = factor(q)
              !
              !$acc parallel loop collapse(2)
              do j=1-nh_p,ny+nh_p
                do i=1-nh_p,nx+nh_p
                  p(i,j,0-q  ) = 1.0_rp*factor_value + p(i,j,1+q  )
                enddo
              enddo
              !$acc end parallel loop
            enddo
            !$OMP END PARALLEL DO
          elseif(ibound.eq.1) then
            !$OMP PARALLEL DO DEFAULT(none) &
            !$OMP PRIVATE(i,j,q,factor,sgn) &
            !$OMP SHARED(nx,ny,nz,nh_p,p)
            do q=0,nh_p-1
              factor_value = factor(q)
              !
              !$acc parallel loop collapse(2)
              do j=1-nh_p,ny+nh_p
                do i=1-nh_p,nx+nh_p
                  p(i,j,nz+q  ) = 1.0_rp*factor_value + p(i,j,nz-1-q)
                  p(i,j,nz+1+q) = 2.0_rp*factor_value + p(i,j,nz-1-q)
                enddo
              enddo
              !$acc end parallel loop
            enddo
            !$OMP END PARALLEL DO
          endif
        end select
      endif
    end select
    !
    return
  end subroutine set_bc

  subroutine set_bc_IBM(nx,ny,nz,ctype,ibound,idir,centered,rvalue,qq_d,nh_p,dr,p,IBM_h)
    !
    ! note: to be generalized with non-homogeneous Neumann
    !
    implicit none
    !
    integer         , intent(in   )                                     :: nx,ny,nz
    character(len=1), intent(in   )                                     :: ctype
    integer         , intent(in   )                                     :: ibound,idir,IBM_h
    logical         , intent(in   )                                     :: centered
    real(rp)        , intent(in   )                                     :: rvalue
    integer         , intent(in   )                                     :: qq_d,nh_p
    real(rp)        , intent(in   ), dimension(0:qq_d)                  :: dr
    real(rp)        , intent(inout), dimension(1-nh_p:,1-nh_p:,1-nh_p:) :: p
    !
    real(rp), dimension(0:qq_d) :: factor
    real(rp) :: factor_value, sgn
    integer  :: i,j,k,q
    !@cuf attributes(managed) :: p
    !
    do q=0,qq_d    !! q_dd is qmin...so it goes from 0,...,abs(1-nh_p-IBM_h)
      factor(q) = rvalue
    enddo
    sgn = 0._rp
    !
    if(ctype.eq.'D'.and.centered) then
      do q=0,qq_d
        factor(q) = +2._rp*factor(q)
      enddo
      sgn = -1._rp
    endif
    if(ctype.eq.'N') then
      if(    ibound.eq.0) then
        do q=0,qq_d
          factor(q) = -dr(q)*factor(q)
        enddo
      elseif(ibound.eq.1) then
        do q=0,qq_d
          factor(q) =  dr(q)*factor(q)
        enddo
      endif
      sgn = +1._rp
    endif
    !
    select case(ctype)
    case('P')
      select case(idir)
      case(1)
        !$OMP PARALLEL DO DEFAULT(none) &
        !$OMP PRIVATE(j,k,q,factor,sgn) &
        !$OMP SHARED(ny,nz,nh_p,p)
        !$acc parallel loop collapse(3)
        do k=1-nh_p,nz+nh_p
          do j=1-nh_p,ny+nh_p
            do q=0,nh_p-1
              p(0-q   ,j,k) = p(nx-q,j,k)
              p(nx+1+q,j,k) = p(1+q ,j,k)
            enddo
          enddo
        enddo
        !$acc end parallel loop
        !$OMP END PARALLEL DO
      case(2)
        !$OMP PARALLEL DO DEFAULT(none) &
        !$OMP PRIVATE(i,k,q,factor,sgn) &
        !$OMP SHARED(nx,nz,nh_p,p)
        !$acc parallel loop collapse(3)
        do k=1-nh_p,nz+nh_p
          do q=0,nh_p-1
            do i=1-nh_p,nx+nh_p
              p(i,0-q   ,k) = p(i,ny-q,k)
              p(i,ny+1+q,k) = p(i,1+q ,k)
            enddo
          enddo
        enddo
        !$acc end parallel loop
        !$OMP END PARALLEL DO
      case(3)
        !$OMP PARALLEL DO DEFAULT(none) &
        !$OMP PRIVATE(i,j,q,factor,sgn) &
        !$OMP SHARED(nx,ny,nz,nh_p,p)
        !$acc parallel loop collapse(3)
        do q=0,nh_p-1
          do j=1-nh_p,ny+nh_p
            do i=1-nh_p,nx+nh_p
              p(i,j,0-q   ) = p(i,j,nz-q)
              p(i,j,nz+1+q) = p(i,j, 1+q)
            enddo
          enddo
        enddo
        !$acc end parallel loop
        !$OMP END PARALLEL DO
      end select
    case('D','N')
      if(centered) then
        select case(idir)
        case(1)
          if    (ibound.eq.0) then
            !$OMP PARALLEL DO DEFAULT(none) &
            !$OMP PRIVATE(j,k,q,factor,sgn) &
            !$OMP SHARED(ny,nz,nh_p,p)
            do q=0,nh_p-1
              factor_value = factor(q)
              !
              !$acc parallel loop collapse(2)
              do k=1-nh_p,nz+nh_p
                do j=1-nh_p,ny+nh_p
                  p(0-q,j,k) = factor_value+sgn*p(1+q,j,k)
                enddo
              enddo
              !$acc end parallel loop
            enddo
            !$OMP END PARALLEL DO
          elseif(ibound.eq.1) then
            !$OMP PARALLEL DO DEFAULT(none) &
            !$OMP PRIVATE(j,k,q,factor,sgn) &
            !$OMP SHARED(nx,ny,nz,nh_p,p)
            do q=0,nh_p-1
              factor_value = factor(q)
              !
              !$acc parallel loop collapse(2)
              do k=1-nh_p,nz+nh_p
                do j=1-nh_p,ny+nh_p
                  p(nx+1+q,j,k) = factor_value+sgn*p(nx-q,j,k)
                enddo
              enddo
              !$acc end parallel loop
            enddo
            !$OMP END PARALLEL DO
          endif
        case(2)
          if    (ibound.eq.0) then
            !$OMP PARALLEL DO DEFAULT(none) &
            !$OMP PRIVATE(i,k,q,factor,sgn) &
            !$OMP SHARED(nx,nz,nh_p,p)
            do q=0,nh_p-1
              factor_value = factor(q)
              !
              !$acc parallel loop collapse(2)
              do k=1-nh_p,nz+nh_p
                do i=1-nh_p,nx+nh_p
                  p(i,0-q  ,k) = factor_value+sgn*p(i,1+q,k)
                enddo
              enddo
              !$acc end parallel loop
            enddo
            !$OMP END PARALLEL DO
          elseif(ibound.eq.1) then
            !$OMP PARALLEL DO DEFAULT(none) &
            !$OMP PRIVATE(i,k,q,factor,sgn) &
            !$OMP SHARED(nx,ny,nz,nh_p,p)
            do q=0,nh_p-1
              factor_value = factor(q)
              !
              !$acc parallel loop collapse(2)
              do k=1-nh_p,nz+nh_p
                do i=1-nh_p,nx+nh_p
                  p(i,ny+1+q,k) = factor_value+sgn*p(i,ny-q,k)
                enddo
              enddo 
              !$acc end parallel loop
            enddo
            !$OMP END PARALLEL DO
          endif
        case(3)
          if    (ibound.eq.0) then
            !$OMP PARALLEL DO DEFAULT(none) &
            !$OMP PRIVATE(i,j,q,factor,sgn) &
            !$OMP SHARED(nx,ny,nh_p,p)
            do q=0,nh_p-1+IBM_h
              factor_value = factor(q)
              ! only this part is changed for IBM from here....!!
              !$acc parallel loop collapse(2)
              do j=1-nh_p,ny+nh_p
                do i=1-nh_p,nx+nh_p
                  p(i,j,0-q + IBM_h ) = factor_value+sgn*p(i,j,1+q + IBM_h)
                enddo 
              enddo
              !$acc end parallel loop
              ! ......to here
            enddo
            !$OMP END PARALLEL DO
          elseif(ibound.eq.1) then
            !$OMP PARALLEL DO DEFAULT(none) &
            !$OMP PRIVATE(i,j,q,factor,sgn) &
            !$OMP SHARED(nx,ny,nz,nh_p,p)
            do q=0,nh_p-1
              factor_value = factor(q)
              !
              !$acc parallel loop collapse(2)
              do j=1-nh_p,ny+nh_p
                do i=1-nh_p,nx+nh_p
                  p(i,j,nz+1+q) = factor_value+sgn*p(i,j,nz-q)
                enddo
              enddo
              !$acc end parallel loop
            enddo
            !$OMP END PARALLEL DO
          endif
        end select
      elseif(.not.centered.and.ctype.eq.'D') then
        select case(idir)
        case(1)
          if    (ibound.eq.0) then
            !$OMP PARALLEL DO DEFAULT(none) &
            !$OMP PRIVATE(j,k,q,factor,sgn) &
            !$OMP SHARED(ny,nz,nh_p,p)
            do q=0,nh_p-1
              factor_value = factor(q)
              !
              !$acc parallel loop collapse(2)
              do k=1-nh_p,nz+nh_p
                do j=1-nh_p,ny+nh_p
                  p(0-q,j,k) = factor_value
                enddo
              enddo
              !$acc end parallel loop
            enddo
            !$OMP END PARALLEL DO
          elseif(ibound.eq.1) then
            !$OMP PARALLEL DO DEFAULT(none) &
            !$OMP PRIVATE(j,k,q,factor,sgn) &
            !$OMP SHARED(nx,ny,nz,nh_p,p)
            do q=0,nh_p-1
              factor_value = factor(q)
              !
              !$acc parallel loop collapse(2)
              do k=1-nh_p,nz+nh_p
                do j=1-nh_p,ny+nh_p
                  p(nx+q  ,j,k) = factor_value
                  p(nx+1+q,j,k) = p(nx-1-q,j,k)
                enddo
              enddo
              !$acc end parallel loop
            enddo
            !$OMP END PARALLEL DO
          endif
        case(2)
          if    (ibound.eq.0) then
            !$OMP PARALLEL DO DEFAULT(none) &
            !$OMP PRIVATE(i,k,q,factor,sgn) &
            !$OMP SHARED(nx,nz,nh_p,p)
            do q=0,nh_p-1
              factor_value = factor(q)
              !
              !$acc parallel loop collapse(2)
              do k=1-nh_p,nz+nh_p
                do i=1-nh_p,nx+nh_p
                  p(i,0-q,k) = factor_value
                enddo
              enddo
              !$acc end parallel loop
            enddo
            !$OMP END PARALLEL DO
          elseif(ibound.eq.1) then
            !$OMP PARALLEL DO DEFAULT(none) &
            !$OMP PRIVATE(i,j,q,factor,sgn) &
            !$OMP SHARED(nx,ny,nz,nh_p,p)
            do q=0,nh_p-1
              factor_value = factor(q)
              !
              !$acc parallel loop collapse(2)
              do k=1-nh_p,nz+nh_p
                do i=1-nh_p,nx+nh_p
                  p(i,ny+q  ,k) = factor_value
                  p(i,ny+1+q,k) = p(i,ny-1-q,k)
                enddo
              enddo
              !$acc end parallel loop
            enddo
            !$OMP END PARALLEL DO
          endif
        case(3)
          if    (ibound.eq.0) then
            !$OMP PARALLEL DO DEFAULT(none) &
            !$OMP PRIVATE(i,j,q,factor,sgn) &
            !$OMP SHARED(nx,ny,nh_p,p)
            do q=0,nh_p-1
              factor_value = factor(q)
              !
              !$acc parallel loop collapse(2)
              do j=1-nh_p,ny+nh_p
                do i=1-nh_p,nx+nh_p
                  p(i,j,0-q) = factor_value 
                enddo
              enddo
              !$acc end parallel loop
            enddo
            !$OMP END PARALLEL DO
          elseif(ibound.eq.1) then
            !$OMP PARALLEL DO DEFAULT(none) &
            !$OMP PRIVATE(i,j,q,factor,sgn) &
            !$OMP SHARED(nx,ny,nz,nh_p,p)
            do q=0,nh_p-1
              factor_value = factor(q)
              !
              !$acc parallel loop collapse(2)
              do j=1-nh_p,ny+nh_p
                do i=1-nh_p,nx+nh_p
                  p(i,j,nz+q  ) = factor_value
                  p(i,j,nz+1+q) = p(i,j,nz-1-q)
                enddo
              enddo
              !$acc end parallel loop
            enddo
            !$OMP END PARALLEL DO
          endif
        end select
      elseif(.not.centered.and.ctype.eq.'N') then
        select case(idir)
        case(1)
          if    (ibound.eq.0) then
            !$OMP PARALLEL DO DEFAULT(none) &
            !$OMP PRIVATE(j,k,q,factor,sgn) &
            !$OMP SHARED(ny,nz,nh_p,p)
            do q=0,nh_p-1
              factor_value = factor(q)
              !
              !$acc parallel loop collapse(2)
              do k=1-nh_p,nz+nh_p
                do j=1-nh_p,ny+nh_p
                  p(0-q,j,k) = 1.0_rp*factor_value + p(1+q,j,k)
                enddo
              enddo
              !$acc end parallel loop
            enddo
            !$OMP END PARALLEL DO
          elseif(ibound.eq.1) then
            !$OMP PARALLEL DO DEFAULT(none) &
            !$OMP PRIVATE(j,k,q,factor,sgn) &
            !$OMP SHARED(nx,ny,nz,nh_p,p)
            do q=0,nh_p-1
              factor_value = factor(q)
              !
              !$acc parallel loop collapse(2)
              do k=1-nh_p,nz+nh_p
                do j=1-nh_p,ny+nh_p
                  p(nx+q  ,j,k) = 1.0_rp*factor_value + p(nx-1-q,j,k)
                  p(nx+1+q,j,k) = 2.0_rp*factor_value + p(nx-1-q,j,k)
                enddo
              enddo
              !$acc end parallel loop
            enddo
            !$OMP END PARALLEL DO
          endif
        case(2)
          if    (ibound.eq.0) then
            !$OMP PARALLEL DO DEFAULT(none) &
            !$OMP PRIVATE(i,k,q,factor,sgn) &
            !$OMP SHARED(nx,nz,nh_p,p)
            do q=0,nh_p-1
              factor_value = factor(q)
              !
              !$acc parallel loop collapse(2)
              do k=1-nh_p,nz+nh_p
                do i=1-nh_p,nx+nh_p
                  p(i,0-q  ,k) = 1.0_rp*factor_value + p(i,1+q  ,k) 
                enddo
              enddo
              !$acc end parallel loop
            enddo
            !$OMP END PARALLEL DO
          elseif(ibound.eq.1) then
            !$OMP PARALLEL DO DEFAULT(none) &
            !$OMP PRIVATE(i,k,q,factor,sgn) &
            !$OMP SHARED(nx,ny,nz,nh_p,p)
            do q=0,nh_p-1
              factor_value = factor(q)
              !
              !$acc parallel loop collapse(2)
              do k=1-nh_p,nz+nh_p
                do i=1-nh_p,nx+nh_p
                  p(i,ny+q  ,k) = 1.0_rp*factor_value + p(i,ny-1-q,k)
                  p(i,ny+1+q,k) = 2.0_rp*factor_value + p(i,ny-1-q,k)
                enddo
              enddo
              !$acc end parallel loop
            enddo
            !$OMP END PARALLEL DO
          endif
        case(3)
          if    (ibound.eq.0) then
            !$OMP PARALLEL DO DEFAULT(none) &
            !$OMP PRIVATE(i,j,q,factor,sgn) &
            !$OMP SHARED(nx,ny,nh_p,p)
            do q=0,nh_p-1
              factor_value = factor(q)
              !
              !$acc parallel loop collapse(2)
              do j=1-nh_p,ny+nh_p
                do i=1-nh_p,nx+nh_p
                  p(i,j,0-q  ) = 1.0_rp*factor_value + p(i,j,1+q  )
                enddo
              enddo
              !$acc end parallel loop
            enddo
            !$OMP END PARALLEL DO
          elseif(ibound.eq.1) then
            !$OMP PARALLEL DO DEFAULT(none) &
            !$OMP PRIVATE(i,j,q,factor,sgn) &
            !$OMP SHARED(nx,ny,nz,nh_p,p)
            do q=0,nh_p-1
              factor_value = factor(q)
              !
              !$acc parallel loop collapse(2)
              do j=1-nh_p,ny+nh_p
                do i=1-nh_p,nx+nh_p
                  p(i,j,nz+q  ) = 1.0_rp*factor_value + p(i,j,nz-1-q)
                  p(i,j,nz+1+q) = 2.0_rp*factor_value + p(i,j,nz-1-q)
                enddo
              enddo
              !$acc end parallel loop
            enddo
            !$OMP END PARALLEL DO
          endif
        end select
      endif
    end select
    !
    return
  end subroutine set_bc_IBM

  subroutine outflow_mom(nx,ny,nz,idir,nh_d,nh_u,dx,dy,dz,dzf,u,v,w,rho,mfx,mfy,mfz,nse_switch,reg_switch)
    !
    use mod_common_mpi, only: MPI_PROC_NULL,left,right,front,back,top,bottom,myid
    use mod_param,      only: rho1,rho2
    !
    implicit none
    !
    integer , intent(in   )                                     :: nx,ny,nz
    integer , intent(in   )                                     :: idir
    integer , intent(in   )                                     :: nh_d,nh_u
    real(rp), intent(in   )                                     :: dx,dy,dz
    real(rp), intent(in   ), dimension(1-nh_d:)                 :: dzf
    real(rp), intent(inout), dimension(1-nh_u:,1-nh_u:,1-nh_u:) :: u,v,w
    real(rp), intent(in   ), dimension(1-nh_u:,1-nh_u:,1-nh_u:) :: mfx,mfy,mfz
    real(rp), intent(in   ), dimension(1-nh_u:,1-nh_u:,1-nh_u:) :: rho
    real(rp), intent(in   )                                     :: nse_switch,reg_switch
    !
    real(rp), dimension(1-nh_d:nz+nh_d) :: dzfi
    real(rp) :: dxi,dyi,rhoxp,rhoxm,rhoyp,rhoym,rhozp,rhozm,div
    integer  :: i,j,k,q,qmin,corr
    !
    qmin = abs(1-nh_u)
    !
    dxi  = dx**(-1)
    dyi  = dy**(-1)
    dzfi = dzf**(-1)
    !
    ! determine face velocity from zero divergence
    ! Note: here we assume that in outflow boundaries, 
    !       there is not phase change or any other volume sources.
    !
    select case(idir)
    case(1) ! x direction, right
      if(right.eq.MPI_PROC_NULL) then
        !i = nx + 1
        i = nx + 0
        !$OMP PARALLEL DO DEFAULT(none) &
        !$OMP PRIVATE(j,k,q) &
        !$OMP SHARED(nz,ny,i,qmin,u,v,w,dx,dyi,dzfi)
        do k=1,nz
          do j=1,ny
            do q=0,qmin+1 ! note that the value at u(i+qmin+1,j,k) is not used
              corr = 0
              IF(q.EQ.1) corr = 1
              rhoxm = 0.5*(rho(i+q-1     ,j  ,k  )+rho(i+q,j,k))
              rhoxp = 0.5*(rho(i+q+1-corr,j  ,k  )+rho(i+q,j,k)) ! assume zero gradient for rho at outflow
              rhoym = 0.5*(rho(i+q       ,j-1,k  )+rho(i+q,j,k))
              rhoyp = 0.5*(rho(i+q       ,j+1,k  )+rho(i+q,j,k))
              rhozm = 0.5*(rho(i+q       ,j  ,k-1)+rho(i+q,j,k))
              rhozp = 0.5*(rho(i+q       ,j  ,k+1)+rho(i+q,j,k))
              div   =-(rho1-rho2)*((mfz(i+q,j,k)/rhozp-mfz(i+q  ,j  ,k-1)/rhozm)*dzfi(k) + &
                                   (mfy(i+q,j,k)/rhoyp-mfy(i+q  ,j-1,k  )/rhoym)*dyi      + &
                                   (mfx(i+q,j,k)/rhoxp-mfx(i+q-1,j  ,k  )/rhoxm)*dxi      )*nse_switch*reg_switch
              u(i+q,j,k) = rhoxp*(u(i-1-q,j,k)/rhoxm + dx*(div-(v(i+q,j,k)/rhoyp-v(i+q,j-1,k)/rhoym)*dyi &
                                                              -(w(i+q,j,k)/rhozp-w(i+q,j,k-1)/rhozm)*dzfi(k)))
            enddo
          enddo
        enddo
        !$OMP END PARALLEL DO
      endif
    case(2) ! y direction, back
      if(back.eq.MPI_PROC_NULL) then
        !j = ny + 1
        j = ny + 0
        !$OMP PARALLEL DO DEFAULT(none) 
        !$OMP PRIVATE(i,k,q) &
        !$OMP SHARED(nx,nz,j,qmin,u,v,w,dy,dxi,dzfi)
        do k=1,nz
          do q=0,qmin+1 ! note that the value at v(i,j+qmin+1,k) is not used
            corr = 0
            IF(q.EQ.1) corr = 1
            do i=1,nx
              rhoxm = 0.5*(rho(i-1,j+q       ,k  )+rho(i,j+q,k))
              rhoxp = 0.5*(rho(i+1,j+q       ,k  )+rho(i,j+q,k))
              rhoym = 0.5*(rho(i  ,j+q-1     ,k  )+rho(i,j+q,k))
              rhoyp = 0.5*(rho(i  ,j+q+1-corr,k  )+rho(i,j+q,k)) ! assume zero gradient for rho at outflow
              rhozm = 0.5*(rho(i  ,j+q       ,k-1)+rho(i,j+q,k))
              rhozp = 0.5*(rho(i  ,j+q       ,k+1)+rho(i,j+q,k))
              !rhoyp = rhoym  ! -> required because otherwise it would collide with the zero gradient velocity BC (ny) = u(ny-1)  -> bullshit
              div   =-(rho1-rho2)*((mfz(i,j+q,k)/rhozp-mfz(i  ,j+q  ,k-1)/rhozm)*dzfi(k) + &
                                   (mfy(i,j+q,k)/rhoyp-mfy(i  ,j+q-1,k  )/rhoym)*dyi      + &
                                   (mfx(i,j+q,k)/rhoxp-mfx(i-1,j+q  ,k  )/rhoxm)*dxi      )*nse_switch*reg_switch
              v(i,j+q,k) = rhoyp*(v(i,j-1-q,k)/rhoym + dy*(div-(u(i,j+q,k)/rhoxp-u(i-1,j+q,k)/rhoxm)*dxi &
                                                              -(w(i,j+q,k)/rhozp-w(i,j+q,k-1)/rhozm)*dzfi(k)))
            enddo
          enddo
        enddo 
        !$OMP END PARALLEL DO
      endif
    case(3) ! z direction, top
      if(top.eq.MPI_PROC_NULL) then
        !k = n(3) + 1
        k = nz + 0
        !$OMP PARALLEL DO DEFAULT(none) &
        !$OMP PRIVATE(i,j,q) &
        !$OMP SHARED(nx,ny,k,qmin,u,v,w,dzf,dxi,dyi)
        do q=0,qmin+1 ! note that the value at w(i,j,k+qmin+1) is not used
          corr = 0
          IF(q.EQ.1) corr = 1
          do j=1,ny
            do i=1,nx
              rhoxm = 0.5*(rho(i-1,j  ,k+q       )+rho(i,j,k+q))
              rhoxp = 0.5*(rho(i+1,j  ,k+q       )+rho(i,j,k+q))
              rhoym = 0.5*(rho(i  ,j-1,k+q       )+rho(i,j,k+q))
              rhoyp = 0.5*(rho(i  ,j+1,k+q       )+rho(i,j,k+q))
              rhozm = 0.5*(rho(i  ,j  ,k+q-1     )+rho(i,j,k+q))
              rhozp = 0.5*(rho(i  ,j  ,k+q+1-corr)+rho(i,j,k+q))  ! assume zero gradient for rho at outflow
              div   =-(rho1-rho2)*((mfz(i,j,k+q)/rhozp-mfz(i  ,j  ,k+q-1)/rhozm)*dzfi(k) + &
                                   (mfy(i,j,k+q)/rhoyp-mfy(i  ,j-1,k+q  )/rhoym)*dyi      + &
                                   (mfx(i,j,k+q)/rhoxp-mfx(i-1,j  ,k+q  )/rhoxm)*dxi      )*nse_switch*reg_switch
              w(i,j,k+q) = rhozp*(w(i,j,k-1-q)/rhozm + dzf(k+q)*(div-(u(i,j,k+q)/rhoxp-u(i-1,j,k+q)/rhoxm)*dxi &
                                                                    -(v(i,j,k+q)/rhoyp-v(i,j-1,k+q)/rhoym)*dyi))
            enddo
          enddo
        enddo 
        !$OMP END PARALLEL DO
      endif
    case(-1) ! x direction, left
      if(left.eq.MPI_PROC_NULL) then
        if(myid.eq.0) print*,'ERROR in outflow_mom(...): This outflow direction is not implemented yet!'
        i = 0
        !$OMP PARALLEL DO DEFAULT(none) &
        !$OMP PRIVATE(j,k,q) &
        !$OMP SHARED(ny,nz,i,qmin,u,v,w,dx,dyi,dzfi)
        do k=1,nz
          do j=1,ny
            do q=0,qmin
              u(i-q,j,k) = u(i+1+q,j,k) + dx*((v(i+1+q,j,k)-v(i+1+q,j-1,k))*dyi+(w(i+1+q,j,k)-w(i+1+q,j,k-1))*dzfi(k))
            enddo
          enddo
        enddo 
        !$OMP END PARALLEL DO
      endif
    case(-2) ! y direction, front
      if(front.eq.MPI_PROC_NULL) then
        if(myid.eq.0) print*,'ERROR in outflow_mom(...): This outflow direction is not implemented yet!'
        j = 0
        !$OMP PARALLEL DO DEFAULT(none) &
        !$OMP PRIVATE(i,k,q) &
        !$OMP SHARED(nx,nz,j,qmin,u,v,w,dy,dxi,dzfi)
        do k=1,nz
          do q=0,qmin
            do i=1,nx
              v(i,j-q,k) = v(i,j+1+q,k) + dy*((u(i,j+1+q,k)-u(i-1,j+1+q,k))*dxi+(w(i,j+1+q,k)-w(i,j+1+q,k-1))*dzfi(k))
            enddo
          enddo
        enddo 
        !$OMP END PARALLEL DO
      endif
    case(-3) ! z direction, bottom
      if(bottom.eq.MPI_PROC_NULL) then
        k = 0
        !$OMP PARALLEL DO DEFAULT(none) &
        !$OMP PRIVATE(i,j,q) &
        !$OMP SHARED(nx,ny,k,qmin,u,v,w,dzf,dxi,dyi)
        do q=0,qmin
          do j=1,ny
            do i=1,nx
              rhoxm = 0.5*(rho(i-1,j  ,k+q+1)+rho(i,j,k+q+1))
              rhoxp = 0.5*(rho(i+1,j  ,k+q+1)+rho(i,j,k+q+1))
              rhoym = 0.5*(rho(i  ,j-1,k+q+1)+rho(i,j,k+q+1))
              rhoyp = 0.5*(rho(i  ,j+1,k+q+1)+rho(i,j,k+q+1))
              rhozm = 0.5*(rho(i  ,j  ,k+q  )+rho(i,j,k+q+1))
              rhozp = 0.5*(rho(i  ,j  ,k+q+2)+rho(i,j,k+q+1))
              div   =-(rho1-rho2)*((mfz(i,j,k+q+1)/rhozp-mfz(i  ,j  ,k+q  )/rhozm)*dzfi(k+1) + &
                                   (mfy(i,j,k+q+1)/rhoyp-mfy(i  ,j-1,k+q+1)/rhoym)*dyi      + &
                                   (mfx(i,j,k+q+1)/rhoxp-mfx(i-1,j  ,k+q+1)/rhoxm)*dxi      )*nse_switch*reg_switch
              w(i,j,k-q) = rhozm*(w(i,j,k+1+q)/rhozp + dzf(k-q+1)*(-div+(u(i,j,k+1+q)/rhoxp-u(i-1,j  ,k+1+q)/rhoxm)*dxi &
                                                                       +(v(i,j,k+1+q)/rhoyp-v(i  ,j-1,k+1+q)/rhoym)*dyi))
            enddo
          enddo
        enddo 
        !$OMP END PARALLEL DO
      endif
    end select
    !
    return
  end subroutine outflow_mom
  
  ! TODO: Add GPU here? 
  subroutine outflow(nx,ny,nz,idir,nh_d,nh_u,dx,dy,dz,dzf,u,v,w,rho,mfx,mfy,mfz,nse_switch,reg_switch)
    !
    use mod_common_mpi, only: MPI_PROC_NULL,left,right,front,back,top,bottom,myid
    use mod_param,      only: rho1,rho2
    !
    implicit none
    !
    integer , intent(in   )                                     :: nx,ny,nz
    integer , intent(in   )                                     :: idir
    integer , intent(in   )                                     :: nh_d,nh_u
    real(rp), intent(in   )                                     :: dx,dy,dz
    real(rp), intent(in   ), dimension(1-nh_d:)                 :: dzf
    real(rp), intent(inout), dimension(1-nh_u:,1-nh_u:,1-nh_u:) :: u,v,w
    real(rp), intent(in   ), dimension(1-nh_u:,1-nh_u:,1-nh_u:) :: mfx,mfy,mfz
    real(rp), intent(in   ), dimension(1-nh_u:,1-nh_u:,1-nh_u:) :: rho
    real(rp), intent(in   )                                     :: nse_switch,reg_switch
    !
    real(rp), dimension(1-nh_d:nz+nh_d) :: dzfi
    real(rp) :: dxi,dyi,rhoxp,rhoxm,rhoyp,rhoym,rhozp,rhozm,div
    integer  :: i,j,k,q,qmin,corr
    !
    qmin = abs(1-nh_u)
    !
    dxi  = dx**(-1)
    dyi  = dy**(-1)
    dzfi = dzf**(-1)
    !
    ! determine face velocity from zero divergence
    ! Note: here we assume that in outflow boundaries, 
    !       there is not phase change or any other volume sources.
    !
    select case(idir)
    case(1) ! x direction, right
      if(right.eq.MPI_PROC_NULL) then
        !i = nx + 1
        i = nx + 0
        !$OMP PARALLEL DO DEFAULT(none) &
        !$OMP PRIVATE(j,k,q) &
        !$OMP SHARED(nz,ny,i,qmin,u,v,w,dx,dyi,dzfi)
        do k=1,nz
          do j=1,ny
            do q=0,qmin ! note that the value at u(i+qmin+1,j,k) is not used
              corr = 0
              IF(q.EQ.1) corr = 1
              rhoxm = 0.5*(rho(i+q-1     ,j  ,k  )+rho(i+q,j,k))
              rhoxp = 0.5*(rho(i+q+1-corr,j  ,k  )+rho(i+q,j,k)) ! assume zero gradient for rho at outflow
              rhoym = 0.5*(rho(i+q       ,j-1,k  )+rho(i+q,j,k))
              rhoyp = 0.5*(rho(i+q       ,j+1,k  )+rho(i+q,j,k))
              rhozm = 0.5*(rho(i+q       ,j  ,k-1)+rho(i+q,j,k))
              rhozp = 0.5*(rho(i+q       ,j  ,k+1)+rho(i+q,j,k))
              div   =-(rho1-rho2)*((mfz(i+q,j,k)/rhozp-mfz(i+q  ,j  ,k-1)/rhozm)*dzfi(k) + &
                                   (mfy(i+q,j,k)/rhoyp-mfy(i+q  ,j-1,k  )/rhoym)*dyi      + &
                                   (mfx(i+q,j,k)/rhoxp-mfx(i+q-1,j  ,k  )/rhoxm)*dxi      )*nse_switch*reg_switch
              u(i+q,j,k) = u(i-1-q,j,k) + dx*(div-(v(i+q,j,k)-v(i+q,j-1,k))*dyi-(w(i+q,j,k)-w(i+q,j,k-1))*dzfi(k))
            enddo
          enddo
        enddo
        !$OMP END PARALLEL DO
      endif
    case(2) ! y direction, back
      if(back.eq.MPI_PROC_NULL) then
        !j = ny + 1
        j = ny + 0
        !$OMP PARALLEL DO DEFAULT(none) &
        !$OMP PRIVATE(i,k,q) &
        !$OMP SHARED(nx,nz,j,qmin,u,v,w,dy,dxi,dzfi)
        do k=1,nz
          do q=0,qmin ! note that the value at v(i,j+qmin+1,k) is not used
            corr = 0
            IF(q.EQ.1) corr = 1
            do i=1,nx
              rhoxm = 0.5*(rho(i-1,j+q       ,k  )+rho(i,j+q,k))
              rhoxp = 0.5*(rho(i+1,j+q       ,k  )+rho(i,j+q,k))
              rhoym = 0.5*(rho(i  ,j+q-1     ,k  )+rho(i,j+q,k))
              rhoyp = 0.5*(rho(i  ,j+q+1-corr,k  )+rho(i,j+q,k)) ! assume zero gradient for rho at outflow
              rhozm = 0.5*(rho(i  ,j+q       ,k-1)+rho(i,j+q,k))
              rhozp = 0.5*(rho(i  ,j+q       ,k+1)+rho(i,j+q,k))
              !rhoyp = rhoym ! -> required because otherwise it would collide with the zero gradient velocity BC (ny) = u(ny-1)  -> bullshit
              div   =-(rho1-rho2)*((mfz(i,j+q,k)/rhozp-mfz(i  ,j+q  ,k-1)/rhozm)*dzfi(k) + &
                                   (mfy(i,j+q,k)/rhoyp-mfy(i  ,j+q-1,k  )/rhoym)*dyi      + &
                                   (mfx(i,j+q,k)/rhoxp-mfx(i-1,j+q  ,k  )/rhoxm)*dxi      )*nse_switch*reg_switch
              v(i,j+q,k) = v(i,j-1-q,k) + dy*(div-(u(i,j+q,k)-u(i-1,j+q,k))*dxi-(w(i,j+q,k)-w(i,j+q,k-1))*dzfi(k))
            enddo
          enddo
        enddo 
        !$OMP END PARALLEL DO
      endif
    case(3) ! z direction, top
      if(top.eq.MPI_PROC_NULL) then
        !k = n(3) + 1
        k = nz + 0
        !$OMP PARALLEL DO DEFAULT(none) &
        !$OMP PRIVATE(i,j,q) &
        !$OMP SHARED(nx,ny,k,qmin,u,v,w,dzf,dxi,dyi)
        do q=0,qmin ! note that the value at w(i,j,k+qmin+1) is not used
          corr = 0
          IF(q.EQ.1) corr = 1
          do j=1,ny
            do i=1,nx
              rhoxm = 0.5*(rho(i-1,j  ,k+q       )+rho(i,j,k+q))
              rhoxp = 0.5*(rho(i+1,j  ,k+q       )+rho(i,j,k+q))
              rhoym = 0.5*(rho(i  ,j-1,k+q       )+rho(i,j,k+q))
              rhoyp = 0.5*(rho(i  ,j+1,k+q       )+rho(i,j,k+q))
              rhozm = 0.5*(rho(i  ,j  ,k+q-1     )+rho(i,j,k+q))
              rhozp = 0.5*(rho(i  ,j  ,k+q+1-corr)+rho(i,j,k+q))  ! assume zero gradient for rho at outflow
              div   =-(rho1-rho2)*((mfz(i,j,k+q)/rhozp-mfz(i  ,j  ,k+q-1)/rhozm)*dzfi(k) + &
                                   (mfy(i,j,k+q)/rhoyp-mfy(i  ,j-1,k+q  )/rhoym)*dyi      + &
                                   (mfx(i,j,k+q)/rhoxp-mfx(i-1,j  ,k+q  )/rhoxm)*dxi      )*nse_switch*reg_switch
              w(i,j,k+q) = w(i,j,k-1-q) + dzf(k+q)*(div-(u(i,j,k+q)-u(i-1,j,k+q))*dxi-(v(i,j,k+q)-v(i,j-1,k+q))*dyi)
            enddo
          enddo
        enddo 
        !$OMP END PARALLEL DO
      endif
    case(-1) ! x direction, left
      if(left.eq.MPI_PROC_NULL) then
        if(myid.eq.0) print*,'ERROR in outflow(...): This outflow direction is not implemented yet!'
        i = 0
        !$OMP PARALLEL DO DEFAULT(none) &
        !$OMP PRIVATE(j,k,q) &
        !$OMP SHARED(ny,nz,i,qmin,u,v,w,dx,dyi,dzfi)
        do k=1,nz
          do j=1,ny
            do q=0,qmin
              u(i-q,j,k) = u(i+1+q,j,k) + dx*((v(i+1+q,j,k)-v(i+1+q,j-1,k))*dyi+(w(i+1+q,j,k)-w(i+1+q,j,k-1))*dzfi(k))
            enddo
          enddo
        enddo 
        !$OMP END PARALLEL DO
      endif
    case(-2) ! y direction, front
      if(front.eq.MPI_PROC_NULL) then
        if(myid.eq.0) print*,'ERROR in outflow(...): This outflow direction is not implemented yet!'
        j = 0
        !$OMP PARALLEL DO DEFAULT(none) &
        !$OMP PRIVATE(i,k,q) &
        !$OMP SHARED(nx,nz,j,qmin,u,v,w,dy,dxi,dzfi)
        do k=1,nz
          do q=0,qmin
            do i=1,nx
              v(i,j-q,k) = v(i,j+1+q,k) + dy*((u(i,j+1+q,k)-u(i-1,j+1+q,k))*dxi+(w(i,j+1+q,k)-w(i,j+1+q,k-1))*dzfi(k))
            enddo
          enddo
        enddo 
        !$OMP END PARALLEL DO
      endif
    case(-3) ! z direction, bottom
      if(bottom.eq.MPI_PROC_NULL) then
        k = 0
        !$OMP PARALLEL DO DEFAULT(none) &
        !$OMP PRIVATE(i,j,q) &
        !$OMP SHARED(nx,ny,k,qmin,u,v,w,dzf,dxi,dyi)
        do q=0,qmin
          do j=1,ny
            do i=1,nx
              rhoxm = 0.5*(rho(i-1,j  ,k+q+1)+rho(i,j,k+q+1))
              rhoxp = 0.5*(rho(i+1,j  ,k+q+1)+rho(i,j,k+q+1))
              rhoym = 0.5*(rho(i  ,j-1,k+q+1)+rho(i,j,k+q+1))
              rhoyp = 0.5*(rho(i  ,j+1,k+q+1)+rho(i,j,k+q+1))
              rhozm = 0.5*(rho(i  ,j  ,k+q  )+rho(i,j,k+q+1))
              rhozp = 0.5*(rho(i  ,j  ,k+q+2)+rho(i,j,k+q+1))
              div   =-(rho1-rho2)*((mfz(i,j,k+q+1)/rhozp-mfz(i  ,j  ,k+q  )/rhozm)*dzfi(k+1+q) + &
                                   (mfy(i,j,k+q+1)/rhoyp-mfy(i  ,j-1,k+q+1)/rhoym)*dyi      + &
                                   (mfx(i,j,k+q+1)/rhoxp-mfx(i-1,j  ,k+q+1)/rhoxm)*dxi      )*nse_switch*reg_switch
              w(i,j,k-q) = w(i,j,k+1+q) + dzf(k+1-q)*(-div+(u(i,j,k+1+q)-u(i-1,j,k+1+q))*dxi+(v(i,j,k+1+q)-v(i,j-1,k+1+q))*dyi)
            enddo
          enddo
        enddo 
        !$OMP END PARALLEL DO
      endif
    end select
    !
    return
  end subroutine outflow
  !
  ! TODO: add GPU here?
  subroutine inflow(nx,ny,nz,idir,nh_u,vel2d,u,v,w)
    !
    use mod_common_mpi, only: MPI_PROC_NULL,left,right,front,back,top,bottom
    !
    implicit none
    !
    integer , intent(in   )                                     :: nx,ny,nz
    integer , intent(in   )                                     :: idir
    integer , intent(in   )                                     :: nh_u
    real(rp), intent(in   ), dimension(     0:,    0:)          :: vel2d
    real(rp), intent(inout), dimension(1-nh_u:,1-nh_u:,1-nh_u:) :: u,v,w
    !
    integer :: i,j,k,q
    !
    select case(idir)
    case(1) ! x direction
      i = 0
      if(left.eq.MPI_PROC_NULL) then
        do k=1,nz
          do j=1,ny
            do q=0,nh_u-1
              u(i-q,j,k) = vel2d(j,k)
            enddo
          enddo
        enddo 
      endif
    case(2) ! y direction
      j = 0
      if(front.eq.MPI_PROC_NULL) then
        do k=1,nz
          do q=0,nh_u-1
            do i=1,nx
              v(i,j-q,k) = vel2d(i,k)
            enddo
          enddo
        enddo 
      endif
    case(3) ! z direction
      k = 0
      if(bottom.eq.MPI_PROC_NULL) then
        do q=0,nh_u-1
          do j=1,ny
            do i=1,nx
              w(i,j,k-q) = vel2d(i,j)
            enddo
          enddo
        enddo
      endif 
    end select
    !
    return
  end subroutine inflow
  !









  subroutine boundsb_array(n,s_bound,nh_p,halo,s)
    !   
    ! imposes scalar boundary conditions (defined on a larger stencil)
    !   
    implicit none
    !   
    integer , intent(in   ), dimension(3)           :: n 
    real(rp), intent(in   ), dimension(1:,1:)       :: s_bound
    integer , intent(in   )                                     :: nh_p
    real(rp), intent(inout), dimension(1-nh_p:,1-nh_p:,1-nh_p:) :: s
    integer , intent(in   ), dimension(3)                       :: halo
    !   
    integer :: q,i,j
    integer :: ind1,ind2
    !@cuf attributes(managed) :: s
    !  
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


    call updthalo(n(1),n(2),n(3),nh_p,halo(ind1),ind1,s)
    call updthalo(n(1),n(2),n(3),nh_p,halo(ind2),ind2,s)


    do q=0,nh_p-1
      do j=1,n(2)
        do i=1,n(1)
          s(i,j,0-q) = s_bound(i,j)
        enddo
      enddo
    enddo  
    !   
    return
  end subroutine boundsb_array

  subroutine boundsb_array_ibm(n,s_bound,nh_p,halo,s,IBM_h)
    !   
    ! imposes scalar boundary conditions (defined on a larger stencil)
    !   
    implicit none
    !   
    integer , intent(in   ), dimension(3)           :: n 
    real(rp), intent(in   ), dimension(1:,1:)       :: s_bound
    integer , intent(in   )                                     :: nh_p,IBM_h
    real(rp), intent(inout), dimension(1-nh_p:,1-nh_p:,1-nh_p:) :: s
    integer , intent(in   ), dimension(3)                       :: halo
    integer :: q,i,j
    integer :: ind1,ind2
    !@cuf attributes(managed) :: s
    !  
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



    do q=0,nh_p-1+IBM_h
      do j=1,n(2)
        do i=1,n(1)
          s(i,j,0-q+ IBM_h ) = s_bound(i,j)
        enddo
      enddo
    enddo
    call updthalo(n(1),n(2),n(3),nh_p,halo(ind1),ind1,s)
    call updthalo(n(1),n(2),n(3),nh_p,halo(ind2),ind2,s)
    !   
    return
  end subroutine boundsb_array_ibm











  subroutine updt_rhs_b(nx,ny,nz,c_or_f,cbc,nh_p,rhsbx,rhsby,rhsbz,p)
    !
    use mod_param     , only: IBM_h
    use mod_common_mpi, only: MPI_PROC_NULL,left,right,front,back,top,bottom
    !
    implicit none
    !
    integer         , intent(in   )                                     :: nx,ny,nz
    character       , intent(in   ), dimension(3)                       :: c_or_f
    character(len=1), intent(in   ), dimension(0:1,3)                   :: cbc
    integer         , intent(in   )                                     :: nh_p
    real(rp)        , intent(in   ), dimension(      :,      :,     0:) :: rhsbx,rhsby,rhsbz
    real(rp)        , intent(inout), dimension(1-nh_p:,1-nh_p:,1-nh_p:) :: p
    !
    integer :: q1, q2, q3
    integer :: i,j,k,idir
    !@cuf attributes(managed) :: p, rhsbx, rhsby, rhsbz
    !
    q1 = 0
    q2 = 0
    q3 = 0
    if(c_or_f(1).eq.'f'.and.cbc(1,1).eq.'D') q1 = 1
    if(c_or_f(2).eq.'f'.and.cbc(1,2).eq.'D') q2 = 1
    if(c_or_f(3).eq.'f'.and.cbc(1,3).eq.'D') q3 = 1
    !
    ! along x
    !
    if(left.eq.MPI_PROC_NULL) then
      !$OMP PARALLEL DO DEFAULT(none) &
      !$OMP PRIVATE(j,k) &
      !$OMP SHARED(ny,nz,p,rhsbx)
      !$acc kernels
      do k=1,nz
        do j=1,ny
          p(1,j,k) = p(1,j,k) + rhsbx(j,k,0)
        enddo
      enddo
      !$acc end kernels
      !$OMP END PARALLEL DO
    endif  
    !
    if(right.eq.MPI_PROC_NULL) then
      !$OMP PARALLEL DO DEFAULT(none) &
      !$OMP PRIVATE(j,k,q) &
      !$OMP SHARED(ny,nx,nz,p,rhsbx)
      !$acc kernels
      do k=1,nz
        do j=1,ny
          p(nx-q1,j,k) = p(nx-q1,j,k) + rhsbx(j,k,1)
        enddo
      enddo
      !$acc end kernels
      !$OMP END PARALLEL DO
    endif
    !
    ! along y
    !
    if(front.eq.MPI_PROC_NULL) then
      !$OMP PARALLEL DO DEFAULT(none) &
      !$OMP PRIVATE(i,k) &
      !$OMP SHARED(nx,nz,p,rhsby)
      !$acc kernels
      do k=1,nz
        do i=1,nx
          p(i,1,k) = p(i,1,k) + rhsby(i,k,0)
        enddo
      enddo
      !$acc end kernels
      !$OMP END PARALLEL DO
    endif
    !
    if(back.eq.MPI_PROC_NULL) then
      !$OMP PARALLEL DO DEFAULT(none) &
      !$OMP PRIVATE(i,k,q) &
      !$OMP SHARED(nx,ny,nz,p,rhsby)
      !$acc kernels
      do k=1,nz
        do i=1,nx
          p(i,ny-q2,k) = p(i,ny-q2,k) + rhsby(i,k,1)
        enddo
      enddo
      !$acc end kernels
      !$OMP END PARALLEL DO
    endif
    !
    ! along z
    !
    if(bottom.eq.MPI_PROC_NULL) then
      !$OMP PARALLEL DO DEFAULT(none) &
      !$OMP PRIVATE(i,k,q) &
      !$OMP SHARED(nx,ny,p,rhsbz)
      !$acc kernels
      do j=1,ny
        do i=1,nx
          p(i,j,1 + IBM_h     ) = p(i,j,1  + IBM_h    ) + rhsbz(i,j,0)
          !p(i,j,1      ) = p(i,j,1      ) + rhsbz(i,j,0)
        enddo
      enddo
      !$acc end kernels
      !$OMP END PARALLEL DO
    endif
    !
    if(top.eq.MPI_PROC_NULL) then
      !$OMP PARALLEL DO DEFAULT(none) &
      !$OMP PRIVATE(i,k,q) &
      !$OMP SHARED(nx,ny,nz,p,rhsbz)
      !$acc kernels
      do j=1,ny
        do i=1,nx
          p(i,j,nz-q3) = p(i,j,nz-q3) + rhsbz(i,j,1)
        enddo
      enddo
      !$acc end kernels
      !$OMP END PARALLEL DO
    endif
    !
    return
  end subroutine updt_rhs_b
  !
  subroutine updthalo(nx,ny,nz,nh,halo,idir,p)
    !
#if defined(_OPENACC)
    use mod_common_mpi, only: xsl_buf, xrl_buf, xsr_buf, xrr_buf, &
                              ysr_buf, yrr_buf, ysl_buf, yrl_buf, &
                              zsr_buf, zrr_buf, zsl_buf, zrl_buf
#endif
    use mod_common_mpi, only: left,right,front,back,top,bottom, &
                              comm_cart,status,ierr
    !
    implicit none
    !
    integer , intent(in   )                               :: nx,ny,nz
    integer , intent(in   )                               :: nh,halo,idir
    real(rp), intent(inout), dimension(1-nh:,1-nh:,1-nh:) :: p
    !@cuf attributes(managed) :: p
    !
    integer :: i,j,k, q, lb1, lb2
    !integer :: requests(4), statuses(MPI_STATUS_SIZE,4)
    !
    !  this subroutine updates the halos that store info
    !  from the neighboring computational sub-domain
    !
    select case(idir)
    case(1) ! x direction
#if defined(_OPENACC)
      !$acc kernels present(xsl_buf,xrl_buf,xsr_buf,xrr_buf)
      xsl_buf(:,:) = - huge(1._rp)
      xrl_buf(:,:) = - huge(1._rp)
      xsr_buf(:,:) = - huge(1._rp)
      xrr_buf(:,:) = - huge(1._rp)
      !$acc end kernels
      !
      do q=0,nh-1
        !$acc kernels present(xsl_buf,xsr_buf)
        do k=1-nh,nz+nh
          do j=1-nh,ny+nh
            xsl_buf(j,k) = p(1+q ,j,k)
            xsr_buf(j,k) = p(nx-q,j,k)
          enddo
        enddo
        !$acc end kernels
        !
        lb1 = lbound(xsl_buf,1) 
        lb2 = lbound(xsl_buf,2) 
        call MPI_SENDRECV(xsl_buf(lb1,lb2), size( xsl_buf ),MPI_REAL_RP,left ,0, &
                          xrr_buf(lb1,lb2), size( xrr_buf ),MPI_REAL_RP,right,0, &
                          comm_cart,status,ierr)
        call MPI_SENDRECV(xsr_buf(lb1,lb2), size( xsr_buf ),MPI_REAL_RP,right,0, &
                          xrl_buf(lb1,lb2), size( xrl_buf ),MPI_REAL_RP,left ,0, &
                          comm_cart,status,ierr)
        !
        ! unpack only from existing neighbours: with MPI_PROC_NULL nothing is received
        ! and the halo must keep the physical boundary value (as in the CPU path below)
        if(right.ne.MPI_PROC_NULL) then
          !$acc kernels present(xrr_buf)
          do k=1-nh,nz+nh
            do j=1-nh,ny+nh
              p(nx+1+q  ,j,k) = xrr_buf(j,k)
            enddo
          enddo
          !$acc end kernels
        endif
        if(left.ne.MPI_PROC_NULL) then
          !$acc kernels present(xrl_buf)
          do k=1-nh,nz+nh
            do j=1-nh,ny+nh
              p(0-q     ,j,k) = xrl_buf(j,k)
            enddo
          enddo
          !$acc end kernels
        endif
      enddo
#else
      call MPI_SENDRECV(p(1      ,1-nh,1-nh),1,halo,left ,0, &
                        p(nx+1   ,1-nh,1-nh),1,halo,right,0, &
                        comm_cart,MPI_STATUS_IGNORE,ierr)
      call MPI_SENDRECV(p(nx+1-nh,1-nh,1-nh),1,halo,right,0, &
                        p(1-nh   ,1-nh,1-nh),1,halo,left ,0, &
                        comm_cart,MPI_STATUS_IGNORE,ierr)
#endif
    case(2) ! y direction
#if defined(_OPENACC)
      !$acc kernels present(ysl_buf,yrl_buf,ysr_buf,yrr_buf)
      ysl_buf(:,:) = - huge(1._rp)
      yrl_buf(:,:) = - huge(1._rp)
      ysr_buf(:,:) = - huge(1._rp)
      yrr_buf(:,:) = - huge(1._rp)
      !$acc end kernels
      !
      do q=0,nh-1
        !$acc kernels present(ysl_buf,ysr_buf)
        do k=1-nh,nz+nh 
          do i=1-nh,nx+nh
            ysl_buf(i,k) = p(i, 1+q,k)
            ysr_buf(i,k) = p(i,ny-q,k)
          enddo
        enddo
        !$acc end kernels
        !
        lb1 = lbound(ysl_buf,1) 
        lb2 = lbound(ysl_buf,2) 
        call MPI_SENDRECV(ysl_buf(lb1,lb2), size( ysl_buf ),MPI_REAL_RP,front,0, &
                          yrr_buf(lb1,lb2), size( yrr_buf ),MPI_REAL_RP,back ,0, &
                          comm_cart,status,ierr)
        call MPI_SENDRECV(ysr_buf(lb1,lb2), size( ysr_buf ),MPI_REAL_RP,back ,0, &
                          yrl_buf(lb1,lb2), size( yrl_buf ),MPI_REAL_RP,front,0, &
                          comm_cart,status,ierr)
        !
        ! unpack only from existing neighbours (see x direction)
        if(back.ne.MPI_PROC_NULL) then
          !$acc kernels present(yrr_buf)
          do k=1-nh,nz+nh 
            do i=1-nh,nx+nh
              p(i,ny+1+q,k) = yrr_buf(i,k)
            enddo
          enddo
          !$acc end kernels
        endif
        if(front.ne.MPI_PROC_NULL) then
          !$acc kernels present(yrl_buf)
          do k=1-nh,nz+nh 
            do i=1-nh,nx+nh
              p(i,   0-q,k) = yrl_buf(i,k)
            enddo
          enddo
          !$acc end kernels
        endif
      enddo
#else
      call MPI_SENDRECV(p(1-nh,1      ,1-nh),1,halo,front,0, &
                        p(1-nh,ny+1   ,1-nh),1,halo,back ,0, &
                        comm_cart,MPI_STATUS_IGNORE,ierr)
      call MPI_SENDRECV(p(1-nh,ny+1-nh,1-nh),1,halo,back ,0, &
                        p(1-nh,1-nh   ,1-nh),1,halo,front,0, &
                        comm_cart,MPI_STATUS_IGNORE,ierr)
#endif
    case(3) ! z direction
#if defined(_OPENACC)
      !$acc kernels present(zsl_buf,zrl_buf,zsr_buf,zrr_buf)
      zsl_buf(:,:) = - huge(1._rp)
      zrl_buf(:,:) = - huge(1._rp)
      zsr_buf(:,:) = - huge(1._rp)
      zrr_buf(:,:) = - huge(1._rp)
      !$acc end kernels
      !
      do q=0,nh-1
        !$acc kernels present(zsl_buf,zsr_buf)
        do j=1-nh,ny+nh 
          do i=1-nh,nx+nh
            zsl_buf(i,j) = p(i,j,  1+q)
            zsr_buf(i,j) = p(i,j, nz-q)
          enddo
        enddo
        !$acc end kernels
        !
        lb1 = lbound(zsl_buf,1) 
        lb2 = lbound(zsl_buf,2) 
        call MPI_SENDRECV(zsl_buf(lb1,lb2), size( zsl_buf ),MPI_REAL_RP,bottom,0, &
                          zrr_buf(lb1,lb2), size( zrr_buf ),MPI_REAL_RP,top   ,0, &
                          comm_cart,status,ierr)
        call MPI_SENDRECV(zsr_buf(lb1,lb2), size( zsr_buf ),MPI_REAL_RP,top   ,0, &
                          zrl_buf(lb1,lb2), size( zrl_buf ),MPI_REAL_RP,bottom,0, &
                          comm_cart,status,ierr)
        !
        ! unpack only from existing neighbours (see x direction)
        if(top.ne.MPI_PROC_NULL) then
          !$acc kernels present(zrr_buf)
          do j=1-nh,ny+nh 
            do i=1-nh,nx+nh
              p(i,j,nz+1+q) = zrr_buf(i,j)
            enddo
          enddo
          !$acc end kernels
        endif
        if(bottom.ne.MPI_PROC_NULL) then
          !$acc kernels present(zrl_buf)
          do j=1-nh,ny+nh 
            do i=1-nh,nx+nh
              p(i,j,   0-q) = zrl_buf(i,j)
            enddo
          enddo
          !$acc end kernels
        endif
      enddo
#else
     call MPI_SENDRECV(p(1-nh,1-nh   ,1   ),1,halo,bottom,0, &
                       p(1-nh,1-nh   ,nz+1),1,halo,top   ,0, &
                       comm_cart,MPI_STATUS_IGNORE,ierr)
     call MPI_SENDRECV(p(1-nh,1-nh,nz+1-nh),1,halo,top   ,0, &
                       p(1-nh,1-nh   ,1-nh),1,halo,bottom,0, &
                       comm_cart,MPI_STATUS_IGNORE,ierr)
#endif
    !
    end select
    !
    return
    !
  end subroutine updthalo


!!!!!!!!!!!Changing it as it can take input for phi and constant parameters!!!!!!!!!!!

subroutine boundp_di(cbc,n,bc,nh_d,nh_p,halo,dl,dzc,dzf,p,cp,h0,tmp)
    !
    use mod_common_mpi, only: MPI_PROC_NULL,left,right,front,back,top,bottom
    !
    ! imposes pressure and scalar boundary conditions
    ! on an arbitrary stencil size
    !
    implicit none
    !
    character(len=1), intent(in   ), dimension(0:1,3)                   :: cbc
    integer         , intent(in   ), dimension(3)                       :: n
    real(rp)        , intent(in   ), dimension(0:1,3)                   :: bc
    integer         , intent(in   )                                     :: nh_d,nh_p 
    integer         , intent(in   ), dimension(3)                       :: halo 
    real(rp)        , intent(in   ), dimension(3)                       :: dl
    real(rp)        , intent(in   ), dimension(1-nh_d:)                 :: dzc,dzf
    real(rp)        , intent(inout), dimension(1-nh_p:,1-nh_p:,1-nh_p:) :: p
    real(rp)        , intent(in), dimension(1-nh_p:,1-nh_p:,1-nh_p:) :: cp,h0,tmp
    !
    !
    real(rp), dimension(0:nh_p-1) :: dr
    integer :: q,qmin, i, j
    integer :: ind1,ind2
    !@cuf attributes(managed) :: p, dzc, dzf
    !
    qmin = abs(1-nh_p)
    !
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
    call updthalo(n(1),n(2),n(3),nh_p,halo(ind1),ind1,p)
    call updthalo(n(1),n(2),n(3),nh_p,halo(ind2),ind2,p)
    !
    ! along x
    !
    do q = 0,qmin
      dr(q) = dl(1)
    enddo
    if(left .eq.MPI_PROC_NULL) then ! x - bottom part
      call set_bc(n(1),n(2),n(3),cbc(0,1),0,1,.true.,bc(0,1),qmin,nh_p,dr,p)
    endif
    if(right.eq.MPI_PROC_NULL) then ! x - upper part
      call set_bc(n(1),n(2),n(3),cbc(1,1),1,1,.true.,bc(1,1),qmin,nh_p,dr,p)
    endif
    !
    ! along y
    !
    do q = 0,qmin
      dr(q) = dl(2)
    enddo
    if(front.eq.MPI_PROC_NULL) then ! y - bottom part
      call set_bc(n(1),n(2),n(3),cbc(0,2),0,2,.true.,bc(0,2),qmin,nh_p,dr,p)
     endif
    if(back .eq.MPI_PROC_NULL) then ! y - upper part
      call set_bc(n(1),n(2),n(3),cbc(1,2),1,2,.true.,bc(1,2),qmin,nh_p,dr,p)
    endif
    !
    ! along z
    !
    do q=0,qmin 
      dr(q) = dzc(-q)
    enddo
    if(bottom.eq.MPI_PROC_NULL) then ! z - bottom boundary
      !call set_bc(n(1),n(2),n(3),cbc(0,3),0,3,.true.,bc(0,3),qmin,nh_p,dr,p,cp,h0,tmp)

    !!!!!!!!!!!!!!!!!here it comes!!!!!!!!!!!!!!!!!!!!
      do j=1-nh_p,n(2)+nh_p
        do i=1-nh_p,n(1)+nh_p
                  p(i,j,0  ) = cp(i,j,0)*tmp(i,j,0)+h0(i,j,0)
        enddo
      enddo 
      
    endif
    do q=0,qmin 
      dr(q) = dzc(n(3)+q)
    enddo
    if(top.eq.MPI_PROC_NULL) then ! z - top boundary
      call set_bc(n(1),n(2),n(3),cbc(1,3),1,3,.true.,bc(1,3),qmin,nh_p,dr,p)
    endif
    !
    return
  end subroutine boundp_di
  !

  !
end module mod_bound
