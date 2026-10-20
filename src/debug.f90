!
! SPDX-License-Identifier: MIT
!
module mod_debug
  !
  use mpi
  use mod_common_mpi, only: myid,ierr,front,back,comm_cart
  use mod_types
  use my_is_nan_mod, only: my_is_nan_real64
  !@cuf use cudafor
  !
  implicit none
  !
  private
  public  :: cmpt_mean,chk_fields,chk_arr
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
  subroutine chk_arr(lbl,a)
    !
    ! debug: print the interior checksum (sum of |a|, one halo layer excluded) of one array.
    ! Active only when the environment variable FLUTAS_NANCHECK is set.
    !
    implicit none
    character(len=*), intent(in)                   :: lbl
    real(rp)        , intent(in), dimension(:,:,:) :: a
    logical, save :: is_first = .true., is_on = .false.
    integer       :: stat
    real(rp)      :: asum
    !@cuf integer :: istat
    !
    if(is_first) then
      is_first = .false.
      call get_environment_variable('FLUTAS_NANCHECK',status=stat)
      is_on = (stat == 0)
    endif
    if(.not.is_on) return
    !@cuf istat = cudaDeviceSynchronize()
    asum = sum_interior(a)
    call mpi_allreduce(MPI_IN_PLACE,asum,1,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
    if(myid == 0) write(*,'(a,a40,es14.6)') ' ARRSUM   ',lbl,asum
  end subroutine chk_arr
  !
  subroutine chk_fields(lbl,u,v,w,p,psi,ent,mut,rho,mu)
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
    real(rp)        , intent(in), dimension(:,:,:) :: u,v,w,p,psi,ent,mut,rho,mu
    !
    character(len=4), parameter :: names(9) = ['u   ','v   ','w   ','p   ','psi ','ent ','mu_t','rho ','mu  ']
    logical, save :: is_first = .true., is_on = .false., reported(9) = .false.
    integer       :: nbad(9),nloc(9),loc(3,9),stat,m
    real(rp)      :: amax(9),asum(9)
    integer       :: nhal(5)
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
    nloc(8) = count_bad(rho,loc(:,8))
    nloc(9) = count_bad(mu ,loc(:,9))
    do m=1,9
      if(nloc(m) > 0 .and. .not.reported(m)) then
        reported(m) = .true.
        write(*,'(a,i4,a,a,a,a,i8,a,3i5,a,3i5)') ' NANCHECK rank',myid,' ',trim(names(m)),' first bad at ',lbl, &
                                                 nloc(m),' cells, first (i,j,k) =',loc(:,m),' of shape',shape(u)
      endif
    enddo
    nbad(:) = nloc(:)
    call mpi_allreduce(MPI_IN_PLACE,nbad,9,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
    if(myid == 0) write(*,'(a,a30,a,9i8)') ' NANCHECK ',lbl,' u,v,w,p,psi,ent,mu_t,rho,mu:',nbad
    !
    ! maximum absolute (finite) value of each field, to spot blow-ups before they become NaN
    amax(1) = max_finite(u  ); amax(2) = max_finite(v  ); amax(3) = max_finite(w  )
    amax(4) = max_finite(p  ); amax(5) = max_finite(psi); amax(6) = max_finite(ent)
    amax(7) = max_finite(mut); amax(8) = max_finite(rho); amax(9) = max_finite(mu )
    call mpi_allreduce(MPI_IN_PLACE,amax,9,MPI_REAL_RP,MPI_MAX,MPI_COMM_WORLD,ierr)
    if(myid == 0) write(*,'(a,a30,a,9es10.2)') ' NANMAX   ',lbl,' u,v,w,p,psi,ent,mu_t,rho,mu:',amax
    !
    ! checksum: sum of |a| over the interior cells (halos excluded), to compare runs
    asum(1) = sum_interior(u  ); asum(2) = sum_interior(v  ); asum(3) = sum_interior(w  )
    asum(4) = sum_interior(p  ); asum(5) = sum_interior(psi); asum(6) = sum_interior(ent)
    asum(7) = sum_interior(mut); asum(8) = sum_interior(rho); asum(9) = sum_interior(mu )
    call mpi_allreduce(MPI_IN_PLACE,asum,9,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
    if(myid == 0) write(*,'(a,a30,a,9es14.6)') ' NANSUM   ',lbl,' u,v,w,p,psi,ent,mu_t,rho,mu:',asum
    !
    ! y-halo check: compare the y-halos with the neighbours' first/last interior rows
    ! (exchanged here on the CPU with plain MPI); counts cells that differ
    nhal(1) = count_yhalo_bad(u  ); nhal(2) = count_yhalo_bad(v  ); nhal(3) = count_yhalo_bad(w  )
    nhal(4) = count_yhalo_bad(p  ); nhal(5) = count_yhalo_bad(psi)
    call mpi_allreduce(MPI_IN_PLACE,nhal,5,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
    if(myid == 0) write(*,'(a,a30,a,5i8)') ' YHALO    ',lbl,' u,v,w,p,psi:',nhal
    !
  end subroutine chk_fields
  !
  integer function count_yhalo_bad(a)
    !
    ! a has one halo layer: 1-based index j=1 is the halo below, j=ny+2 the halo above.
    ! Send the first interior row (j=2) to front and receive the back neighbour's one: it
    ! must equal our upper halo (j=ny+2). Same for the lower halo with the last interior row.
    !
    implicit none
    real(rp), intent(in), dimension(:,:,:) :: a
    real(rp), allocatable, dimension(:,:)  :: sbuf,rbuf
    integer :: nx2,ny2,nz2,i,k,cnt
    !
    nx2 = size(a,1); ny2 = size(a,2); nz2 = size(a,3)
    allocate(sbuf(nx2,nz2),rbuf(nx2,nz2))
    cnt = 0
    !
    sbuf(:,:) = a(:,2,:)
    rbuf(:,:) = 0._rp
    call MPI_SENDRECV(sbuf,nx2*nz2,MPI_REAL_RP,front,11,rbuf,nx2*nz2,MPI_REAL_RP,back ,11,comm_cart,MPI_STATUS_IGNORE,ierr)
    if(back /= MPI_PROC_NULL) then
      do k=2,nz2-1
        do i=2,nx2-1
          if(abs(a(i,ny2,k)-rbuf(i,k)) > 1.e-12_rp*max(1._rp,abs(rbuf(i,k)))) cnt = cnt + 1
        enddo
      enddo
    endif
    !
    sbuf(:,:) = a(:,ny2-1,:)
    rbuf(:,:) = 0._rp
    call MPI_SENDRECV(sbuf,nx2*nz2,MPI_REAL_RP,back ,12,rbuf,nx2*nz2,MPI_REAL_RP,front,12,comm_cart,MPI_STATUS_IGNORE,ierr)
    if(front /= MPI_PROC_NULL) then
      do k=2,nz2-1
        do i=2,nx2-1
          if(abs(a(i,1,k)-rbuf(i,k)) > 1.e-12_rp*max(1._rp,abs(rbuf(i,k)))) cnt = cnt + 1
        enddo
      enddo
    endif
    !
    deallocate(sbuf,rbuf)
    count_yhalo_bad = cnt
  end function count_yhalo_bad
  !
  real(rp) function sum_interior(a)
    ! sum of |a| over the interior (one halo layer excluded), finite entries only
    implicit none
    real(rp), intent(in), dimension(:,:,:) :: a
    integer :: i,j,k
    sum_interior = 0._rp
    do k=2,size(a,3)-1
      do j=2,size(a,2)-1
        do i=2,size(a,1)-1
          if(.not.(my_is_nan_real64(a(i,j,k)) .or. abs(a(i,j,k)) >= huge(1._rp))) &
            sum_interior = sum_interior + abs(a(i,j,k))
        enddo
      enddo
    enddo
  end function sum_interior
  !
  real(rp) function max_finite(a)
    ! max |a| over the finite entries of a
    implicit none
    real(rp), intent(in), dimension(:,:,:) :: a
    integer :: i,j,k
    max_finite = 0._rp
    do k=1,size(a,3)
      do j=1,size(a,2)
        do i=1,size(a,1)
          if(.not.(my_is_nan_real64(a(i,j,k)) .or. abs(a(i,j,k)) >= huge(1._rp))) &
            max_finite = max(max_finite,abs(a(i,j,k)))
        enddo
      enddo
    enddo
  end function max_finite
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
          if(my_is_nan_real64(a(i,j,k)) .or. abs(a(i,j,k)) >= huge(1._rp)) then
            count_bad = count_bad + 1
            if(count_bad == 1) loc(:) = [i,j,k]
          endif
        enddo
      enddo
    enddo
  end function count_bad
  !
end module mod_debug
