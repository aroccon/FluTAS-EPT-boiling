!
! SPDX-License-Identifier: MIT
!
module mod_debug
  !
  use mpi
  use mod_common_mpi, only: myid,ierr
  use mod_types
  use my_is_nan_mod, only: my_is_nan_real64
  !@cuf use cudafor
  !
  implicit none
  !
  private
  public  :: cmpt_mean,chk_fields
  !
  contains
  !
  subroutine cmpt_mean(nx,ny,nz,nh_d,nh_p,dx,dy,dzf,lx,ly,lz,p,mean)
    !
    ! compute the mean value of an observable over the entire domain
    !
    implicit none
    !
    integer , intent(in )                                     :: nx,ny,nz
    integer , intent(in )                                     :: nh_d,nh_p
    real(rp), intent(in )                                     :: dx,dy
    real(rp), intent(in ), dimension(1-nh_d:)                 :: dzf
    real(rp), intent(in )                                     :: lx,ly,lz
    real(rp), intent(in ), dimension(1-nh_p:,1-nh_p:,1-nh_p:) :: p
    real(rp), intent(out)                                     :: mean
    !
    !@cuf attributes(managed) :: p, dzf
    integer :: i,j,k
    !
    mean = 0._rp
    !
#if defined(_OPENACC)
    !$acc parallel loop collapse(3) reduction(+:mean)
#else
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP SHARED(n,p,dzf) &
    !$OMP PRIVATE(i,j,k) &
    !$OMP REDUCTION(+:mean)
#endif
    do k=1,nz
      do j=1,ny
        do i=1,nx
          mean = mean + p(i,j,k)*dx*dy*dzf(k)
        enddo
      enddo
    enddo
#if defined(_OPENACC)
    !$acc end parallel loop 
#else
    !$OMP END PARALLEL DO
#endif
    !
    call mpi_allreduce(MPI_IN_PLACE,mean,1,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
    mean = mean/(1._rp*lx*ly*lz)
    !
    return
  end subroutine cmpt_mean
  !
  !
  subroutine chk_fields(lbl,u,v,w,p,psi,ent,mut)
    !
    ! debug: count NaN/Inf in the main fields and print them with label lbl.
    ! The first time a field becomes bad on a rank, that rank also prints the
    ! first bad cell as (i,j,k) of the allocated array (1-based, i.e. halo layer 0 -> 1).
    ! Active only when the environment variable FLUTAS_NANCHECK is set (any value),
    ! otherwise it returns immediately.
    !
    implicit none
    !
    character(len=*), intent(in)                 :: lbl
    real(rp)        , intent(in), dimension(:,:,:) :: u,v,w,p,psi,ent,mut
    !
    character(len=4), parameter :: names(7) = ['u   ','v   ','w   ','p   ','psi ','ent ','mu_t']
    logical, save :: is_first = .true., is_on = .false., reported(7) = .false.
    integer       :: nbad(7),nloc(7),loc(3,7),stat,m
    !@cuf integer :: istat
    !
    if(is_first) then
      is_first = .false.
      call get_environment_variable('FLUTAS_NANCHECK',status=stat)
      is_on = (stat == 0)
    endif
    if(.not.is_on) return
    !
    !@cuf istat = cudaDeviceSynchronize()
    nloc(1) = count_bad(u  ,loc(:,1))
    nloc(2) = count_bad(v  ,loc(:,2))
    nloc(3) = count_bad(w  ,loc(:,3))
    nloc(4) = count_bad(p  ,loc(:,4))
    nloc(5) = count_bad(psi,loc(:,5))
    nloc(6) = count_bad(ent,loc(:,6))
    nloc(7) = count_bad(mut,loc(:,7))
    do m=1,7
      if(nloc(m) > 0 .and. .not.reported(m)) then
        reported(m) = .true.
        write(*,'(a,i4,a,a,a,a,i8,a,3i5,a,3i5)') ' NANCHECK rank',myid,' ',trim(names(m)),' first bad at ',lbl, &
                                                 nloc(m),' cells, first (i,j,k) =',loc(:,m),' of shape',shape(u)
      endif
    enddo
    nbad(:) = nloc(:)
    call mpi_allreduce(MPI_IN_PLACE,nbad,7,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
    if(myid == 0) write(*,'(a,a30,a,7i9)') ' NANCHECK ',lbl,' u,v,w,p,psi,ent,mu_t:',nbad
    !
  end subroutine chk_fields
  !
  integer function count_bad(a,loc)
    implicit none
    real(rp), intent(in), dimension(:,:,:) :: a
    integer , intent(out)                  :: loc(3)
    integer :: i,j,k
    count_bad = 0
    loc(:) = 0
    do k=1,size(a,3)
      do j=1,size(a,2)
        do i=1,size(a,1)
          if(my_is_nan_real64(a(i,j,k)) .or. abs(a(i,j,k)) > huge(1._rp)) then
            count_bad = count_bad + 1
            if(count_bad == 1) loc(:) = [i,j,k]
          endif
        enddo
      enddo
    enddo
  end function count_bad
  !
end module mod_debug
