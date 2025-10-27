!
! SPDX-License-Identifier: MIT
!
module mod_chkdiv
  !
  use mpi
  use mod_common_mpi, only: myid,ierr
  use mod_types
  !@cuf use cudafor
  !
  implicit none
  !
  private
  public  :: chkdiv
  !
  contains
  !
  subroutine chkdiv(nx,ny,nz,dxi,dyi,dzi,nh_d,nh_u,dzfi,u,v,w,divtot,divmax,rho,rhoold,f_t12,div_,src_ph)
    !
    ! checks the divergence of the velocity field
    use mod_param    , only: rho1,rho2
    !
    implicit none
    !
    integer , intent(in )                                     :: nx,ny,nz
    real(rp), intent(in )                                     :: dxi,dyi,dzi,f_t12
    integer , intent(in )                                     :: nh_d,nh_u
    real(rp), intent(in ), dimension(1-nh_d:)                 :: dzfi
    real(rp), intent(in ), dimension(1-nh_u:,1-nh_u:,1-nh_u:) :: u,v,w,rho,rhoold,src_ph
    real(rp), intent(inout ), dimension(1:,1:,1:)             :: div_
    real(rp), intent(out)                                     :: divtot,divmax
    !
    real(rp) :: div,mass_cons,mass_cons_max,mass_cons_tot,rhoxp,rhoyp,rhozp,rhoxm,rhoym,rhozm,resmax
    !@cuf attributes(managed) :: u,v,w,dzfi
    integer :: i,j,k
    !
    divtot        = 0._rp
    resmax        = 0._rp
    divmax        = 0._rp
    mass_cons_tot = 0._rp
    mass_cons_max = 0._rp
    div_(:,:,:) = 0.0
#if defined(_OPENACC)
    !$acc parallel loop collapse(3) reduction(max:divmax) reduction(+:divtot)
#endif
    do k=1,nz
      do j=1,ny
        do i=1,nx
          !
          div       = (w(i,j,k)-w(i,j,k-1))*dzfi(k) + &
                      (v(i,j,k)-v(i,j-1,k))*dyi     + &
                      (u(i,j,k)-u(i-1,j,k))*dxi

          divmax = max(divmax,abs(div))

          divtot = divtot + div
          div_(i,j,k) = div - src_ph(i,j,k)*(1.d0/rho1 - 1.d0/rho2)
          resmax = max(resmax,abs(div_(i,j,k)))
          !
        enddo
      enddo
    enddo
#if defined(_OPENACC)
    !$acc end parallel loop 
#endif
    call mpi_allreduce(MPI_IN_PLACE,divtot,1,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
    call mpi_allreduce(MPI_IN_PLACE,divmax,1,MPI_REAL_RP,MPI_MAX,MPI_COMM_WORLD,ierr)
    call mpi_allreduce(MPI_IN_PLACE,resmax,1,MPI_REAL_RP,MPI_MAX,MPI_COMM_WORLD,ierr)

    call mpi_allreduce(MPI_IN_PLACE,mass_cons_tot,1,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
    call mpi_allreduce(MPI_IN_PLACE,mass_cons_max,1,MPI_REAL_RP,MPI_MAX,MPI_COMM_WORLD,ierr)
    if(myid.eq.0) print*, 'Total divergence = ', divtot, '| Maximum divergence = ', divmax, ' (max. vel. res.: ', resmax, ')'
    !
    return
  end subroutine chkdiv
  !
end module mod_chkdiv
