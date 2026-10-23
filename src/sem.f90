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
!       Generalization of Synthetic Eddy Method by Franco Magagnato
!       based on the code from J. Canton and L. Hufnagel (KTH)
!       
!       Turbulence at inlet is "generated" WITH divergence-free 
!       synthetic eddy method (DFSyEMiso)
!       Poletto, Ruggero: Divergence free development of the synthetic eddy method 
!       in order to improve synthetic turbulence for embedded les simulations,
!       PhD Thesis, The University of Manchester, 2015
!       
!--------------------------------------------------------------------------------
!       
!       FOR Poletto (SEM)
!       User needs to provide:
!       
!       -A few input parameters in file called sem_input_1.txt
!       (e.g. # of eddies, MAX(sigma), ..
!       containing 
!       -normal distance to the wall
!       -umean(pos)
!       -k(pos)
!--------------------------------------------------------------------------------
!
  MODULE mod_sem
    !
    use mod_types
    use mpi
    use mod_common_mpi, only: myid,ierr,ijk_start,comm_cart,status,n_x,n_y,n_z
    use iso_c_binding , only: C_PTR
    use mod_param     , only: itot,jtot,ktot
    use profiler
    !
    IMPLICIT NONE
    SAVE

    INTEGER, PARAMETER :: sem_mod = 1
    INTEGER, PARAMETER :: max_fam = 1

    INTEGER,ALLOCATABLE :: nlines(:) ! number of elements at sem_input families

    integer ,ALLOCATABLE   :: dir_xyz(:)


    !     Deviation of bulk velocity of generated signal

    real(rp) ,ALLOCATABLE   :: inlet_area(:)

    CHARACTER(80) infile,restart_file

    !     extent, Volume of the virtual SyEM domain
    integer,parameter       :: neddy = 2000   ! must match the eddy count in sem_restart_1.txt (single_nuc_site)
    real(rp) ,ALLOCATABLE   :: bulk_vel_diff(:,:) ! bulk velocity
    real(rp) ,ALLOCATABLE   :: z0(:),u0(:,:) ! bulk velocity
    real(rp) ,ALLOCATABLE   :: sigma_max(:),yplus_cutoff(:)
    real(rp) ,ALLOCATABLE   :: xbmin(:),xbmax(:),ybmin(:),ybmax(:)
    real(rp) ,ALLOCATABLE   :: Vb(:),zbmin(:),zbmax(:),norm(:,:)
    real(rp) ,ALLOCATABLE   :: x_inlet(:) ! x mean value of inlet plane
    real(rp) ,ALLOCATABLE   :: y_inlet(:) ! y mean value of inlet plane
    real(rp) ,ALLOCATABLE   :: z_inlet(:) ! z mean value of inlet plane
    real(rp) ,ALLOCATABLE   :: t_inlet(:) ! z mean value of inlet plane
    real(rp), ALLOCATABLE   :: ex(:),ey(:),ez(:),eps(:,:)
    real(rp), DIMENSION(1:neddy,1:6)   :: eddy
    integer                 :: start_edd,stop_edd

    !     box in rotated frame (SyEM)
    real(rp)               :: smin,smax,rmin,rmax,tmin,tmax
    real(rp)               :: s_box,r_box,t_box

    real(rp) , ALLOCATABLE :: umean_inlet(:,:,:)
    real(rp) , ALLOCATABLE :: sigma_inlet(:,:,:)
    real(rp) , ALLOCATABLE :: u_bc(:,:),v_bc(:,:),w_bc(:,:)
    real(rp) , ALLOCATABLE :: sig_x(:,:,:,:,:)
    !     eddy size at inlet face
    real(rp) , ALLOCATABLE  :: intensity_inlet(:,:,:)
    real(rp) , ALLOCATABLE  :: int_x(:,:,:,:,:)

    !     Coordinates and Area at inlet planes 
    real(rp) ,ALLOCATABLE   :: ydist(:,:,:)
    real(rp) ,ALLOCATABLE   :: xyz_wall(:,:,:,:)
    real(rp) ,ALLOCATABLE   :: xwall(:,:,:,:)

    !INTEGER, ALLOCATABLE :: max_faces(:),max_edges(:)
    INTEGER, ALLOCATABLE :: fam_gl(:,:)
   integer :: iseed=-5
   real, dimension(1:32) :: iv=0.0
   integer :: iy=0
    INTEGER       :: read_sem,max_rank,fgl(100),floc(100)
    real(rp),ALLOCATABLE :: pos(:,:),umean(:,:),tke(:,:),dissip(:,:)
    real(rp),ALLOCATABLE :: rs(:,:,:),lt(:,:,:)
    real(rp),ALLOCATABLE ::  max_range(:)
    CHARACTER(3),ALLOCATABLE ::  cbc_gl(:,:)
    !
    integer             ::  receiver_rank
    integer,ALLOCATABLE ::  sender_ranks(:)
    
    contains 

    SUBROUTINE SEM_setup(n)
      !
      !use ziggurat, only: init_random_seed
      !
      IMPLICIT NONE
      !
      integer                                :: i
      real(rp)                               :: rn
      integer , intent(in ), dimension(3)    :: n
      !
      !
      ALLOCATE(yplus_cutoff(max_fam))
      ALLOCATE(sigma_max(max_fam))
      ALLOCATE(nlines(max_fam))
      ALLOCATE(dir_xyz(max_fam))
      !
      infile      ='sem_input_1.txt'
      restart_file='sem_restart_1.txt'
      !
      if(myid.eq.0) print*, 'Read in SEM data file'
      !
      call SEM_param_in()      
      !
      If(dir_xyz(max_fam).NE.0) then   !  -> if dir_xyz(fam) is zero = no SEM use!
        !
        if(myid.eq.0) print*, 'Setting up SEM module'
        !
        ALLOCATE(inlet_area(max_fam))
        ALLOCATE(xbmin(max_fam))
        ALLOCATE(xbmax(max_fam))
        ALLOCATE(ybmin(max_fam))
        ALLOCATE(ybmax(max_fam))
        ALLOCATE(zbmin(max_fam))
        ALLOCATE(zbmax(max_fam))
        ALLOCATE(vb(max_fam))
        ALLOCATE(u0(3,max_fam))
        ALLOCATE(z0(max_fam))
        ALLOCATE(bulk_vel_diff(3,max_fam))
        ALLOCATE(max_range(max_fam))
        !
        z0 = 0._rp
        u0 = 0._rp
        !
        if(myid.eq.0) print*, 'Determine SEM parameters'
        !
        call SEMinit(n)
        !
        if(myid.eq.0) print*, 'Generate Eddy distribution'
        !
        !
        call random_seed()
        !
        !
        CALL gen_eddy_all()
        !
        CALL SEMrestart('R')   ! if a restart file exists  -> generated eddies are overwritten
        CALL SEMrestart('W')   ! Output the initial eddy distribution!
        !
        if(myid.eq.0) print*, 'SEM module is ready to go...'
        !
      endif
      !
      RETURN
    END SUBROUTINE SEM_setup

!=======================================================================
!     READ parameters SEM
!=======================================================================

    SUBROUTINE SEM_param_in
      !
      IMPLICIT NONE
      !
      LOGICAL           :: semstop
      INTEGER           :: i,e,fid, iface,inlen,iwork
      real(rp)          :: work
      CHARACTER         :: text*132
      INTEGER,PARAMETER :: fam = 1
      !
      dir_xyz   = 0
      semstop   = .TRUE.
      !
      !     CHECK FOR INPUT FILE
      !
      if(myid.eq.0) INQUIRE(file=infile,exist=semstop); semstop=.NOT.semstop
      !
      CALL mpi_bcast(semstop,1,MPI_LOGICAL,0,comm_cart,ierr)
      !
      IF (semstop) THEN
        if(myid.eq.0) print*, 'SEM input file: ', infile , ' does not exist'
        if(myid.eq.0) print*, 'SEM Module is disabled for this simulation'
        dir_xyz = 0
      END IF
      !
      IF(.NOT.semstop) THEN
        !
        IF(myid.EQ.0) then
        OPEN(newunit=fid,file=infile,form='formatted',status='old',action='read')
        !
        READ(fid,*)      ! skip header
        READ(fid,*) nlines(fam),sigma_max(fam),  &
                    yplus_cutoff(fam),dir_xyz(fam)
        !
        READ(fid,*)      ! skip header
        ENDIF
        !
        CALL mpi_bcast(nlines,max_fam,MPI_INTEGER,0,comm_cart,ierr)
        CALL mpi_bcast(sigma_max,max_fam,MPI_REAL_RP,0,comm_cart,ierr)
        CALL mpi_bcast(yplus_cutoff,max_fam,MPI_REAL_RP,0,comm_cart,ierr)
        CALL mpi_bcast(dir_xyz,max_fam,MPI_INTEGER,0,comm_cart,ierr)
        !
        ALLOCATE(pos(nlines(fam),max_fam))
        ALLOCATE(umean(nlines(fam),max_fam))
        ALLOCATE(tke(nlines(fam),max_fam))
        ALLOCATE(dissip(nlines(fam),max_fam))
        !
        pos       = 0._rp
        umean     = 0._rp
        tke       = 0._rp
        dissip    = 0._rp
        !
        IF(myid.EQ.0) then
        !
        DO i=1, nlines(fam)
          READ(fid,*) pos(i,fam), umean(i,fam), tke(i,fam), dissip(i,fam)
        ENDDO
        !
        CLOSE(fid)
        ENDIF
        CALL mpi_bcast(pos,max_fam*nlines(max_fam),MPI_REAL_RP,0,comm_cart,ierr)
        CALL mpi_bcast(umean,max_fam*nlines(max_fam),MPI_REAL_RP,0,comm_cart,ierr)
        CALL mpi_bcast(tke,max_fam*nlines(max_fam),MPI_REAL_RP,0,comm_cart,ierr)
        CALL mpi_bcast(dissip,max_fam*nlines(max_fam),MPI_REAL_RP,0,comm_cart,ierr)
      ENDIF
      !
    RETURN
  END  SUBROUTINE SEM_param_in

!=======================================================================
! Inititialize it
!=======================================================================

  SUBROUTINE SEMinit(n)
    use mod_param, only: dx,IBM_h,lx,lz,dims_in,dy,dz
    use mod_common_mpi, only: MPI_PROC_NULL,left,right,front,back,top,bottom
    IMPLICIT NONE

    real(rp), PARAMETER :: pi = 3.14159265358979323846264338327950288
    integer , intent(in ), dimension(3)    :: n
    INTEGER,PARAMETER :: fam = 1
    INTEGER e, i, j, k, eg,fid, iface,num_inlet,max_fac,test,inx,iny
    INTEGER max_edg,ii,nproc,upper_cnt,lower_cnt,over_full,start_testx,rank,dir_id
    integer start_testy
    real(rp) :: min_sigma

    eddy = 0.0
    rmin = 0.0
    rmax = 0.0
    smin = 0.0
    smax = 0.0
    tmin = 0.0
    tmax = 0.0
    s_box = 0.0
    r_box = 0.0
    t_box = 0.0

    select case(dir_xyz(fam))
    case(2)
          ! Set the box size
          rmin            = 0.0                    !always spanwise
          rmax            = lx                     !always spanwise
          r_box           = lx                     !spanwise box size
          smin            = 0.0 + dz*IBM_h         !always normal
          smax            = lz                     !always normal
          s_box           = (lz-dz*IBM_h)          !wall_normal box size
          sigma_max(fam)  = MIN(smax-smin,rmax-rmin,sigma_max(fam))
          tmin            = 0.0 - sigma_max(fam)   !always streamwise
          tmax            = 0.0 + sigma_max(fam)   !always streamwise
          t_box           = 2*sigma_max(fam)
          inlet_area(fam) = lx*(lz-dz*IBM_h)
          max_range(fam)  = (lz-dz*IBM_h)
          !
          if(myid.EQ.0) print*, 'SEM Box size:'
          if(myid.EQ.0) print*, 'Span-wise (min,max):', rmin,rmax
          if(myid.EQ.0) print*, 'Wall-normal (min,max):', smin,smax
          if(myid.EQ.0) print*, 'Stream-wise (min,max):', tmin,tmax
          !
          min_sigma = 1e6_rp
          !
          !Allocate arrays
          ALLOCATE(sigma_inlet    (n(1),n(3),max_fam))
          ALLOCATE(intensity_inlet(n(1),n(3),max_fam))
          ALLOCATE(umean_inlet    (n(1),n(3),max_fam))
          !
          sigma_inlet     = 0._rp
          intensity_inlet = 0._rp
          umean_inlet     = 0._rp
          CALL SEM_u0_norm(fam,2,n)
          !
          min_sigma = minval(sigma_inlet(:,(1 + IBM_h):n(3),fam))
          !
          if(front.eq.MPI_PROC_NULL) then
            ALLOCATE(u_bc           (n(1),n(3))        )
            ALLOCATE(v_bc           (n(1),n(3))        )
            ALLOCATE(w_bc           (n(1),n(3))        )
            u_bc            = 0.0_rp
            v_bc            = 0.0_rp
            w_bc            = 0.0_rp
          endif
          !
          start_edd = 1
          stop_edd  = neddy
          !
          nproc = dims_in(2)   !dims_in(1)*
          !
          lower_cnt = FLOOR(real(neddy)/real(nproc))
          upper_cnt = lower_cnt + 1
          over_full = MOD(neddy,nproc)
          !
          dir_id = int(ijk_start(2)/n(2))  ! n-th processor in the x2 direction
          !
          ! Distribute the calculation of the eddy contribution over the processors in the x2 direction
          if(dir_id.LT.over_full) then
            start_edd = 1 + dir_id*upper_cnt
            stop_edd  = (dir_id+1)*upper_cnt
          else
            start_edd = over_full*upper_cnt + 1 + (dir_id  -over_full)*lower_cnt
            stop_edd  = over_full*upper_cnt +     (dir_id+1-over_full)*lower_cnt
          endif
          print*, 'Rank: ', myid, ' with directional id: ',dir_id,  ' is calculating the contribution of eddies ', &
                                        start_edd , ' to ', stop_edd
          !
          !
          ! Communicate which ranks are inlet ranks and what the starting address of the global array is for each rank:
          !
          receiver_rank = 0
          allocate(sender_ranks(1:dims_in(2)))
          sender_ranks(:) = 0
          !
          ! determine receiver rank
          do inx = 1,dims_in(1)
            rank       = 0
            start_testx = (inx-1)*n(1)
            if ((ijk_start(1).EQ.start_testx).AND.(ijk_start(2).EQ.0)) then   ! i.e. the rank everyone wants to send their data to
              rank = myid
            endif
            call mpi_allreduce(MPI_IN_PLACE,rank ,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ierr)
            if (ijk_start(1).EQ.start_testx) then   ! all the ranky in that row
              receiver_rank = rank
            endif 
          enddo
          !
          ! determine sender ranks
          do inx = 1,dims_in(1)
            do iny = 1,dims_in(2)
              rank       = 0
              start_testx = (inx-1)*n(1)
              start_testy = (iny-1)*n(2)
              if ((ijk_start(1).EQ.start_testx).AND.(ijk_start(2).EQ.start_testy)) then   ! i.e. go through all ranks
                rank = myid
              endif
              call mpi_allreduce(MPI_IN_PLACE,rank ,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ierr)
              if (ijk_start(1).EQ.start_testx) then   ! i.e. the rank everyone wants to send their data to
                sender_ranks(iny) = rank
              endif
            enddo
          enddo
          if(myid.eq.receiver_rank) print*,'Receiver:' ,receiver_rank, ' senders: ', sender_ranks
          !
          !
          call mpi_allreduce(MPI_IN_PLACE,u0(2,fam) ,1,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
          Vb(fam) = inlet_area(fam)*2.*sigma_max(fam)
          u0(2,fam) = u0(2,fam)/inlet_area(fam)
          z0(fam)   = u0(2,fam)
          !
          if((myid.EQ.0).AND.(min_sigma.NE.0.0_rp)) print*, 'The recommended eddy count is', ((Vb(fam))/(min_sigma**3))
    case(0)
      ! Do nothing
    case default
      if(myid.eq.0) print*, 'ERROR: This SEM direction is not implemented ... yet!'
      call MPI_FINALIZE(ierr)
      call exit
    end select
    

    RETURN
  END SUBROUTINE SEMinit

!=======================================================================
!     SEM Main routine for DFSEM
!=======================================================================

  SUBROUTINE synthetic_eddies(dt,n)

    IMPLICIT NONE
    !
    real(rp) ,ALLOCATABLE                  :: work_array(:,:)
    integer , intent(in ), dimension(3)    :: n
    INTEGER                                :: i
    INTEGER,PARAMETER                      :: fam = 1
    real(rp)                               :: in_area,umaxmax
    real(rp), intent(in )                  :: dt
    !
    !  Update eddies
    If(dir_xyz(max_fam).NE.0) then   !  -> if dir_xyz(fam) is zero = no SEM use!
      call profiler_start("SEM:advect")
      CALL advect_recycle_eddies(neddy,dt)
      call profiler_stop("SEM:advect")
      !  Calculate BC field
      call profiler_start("SEM:iso_dir")
      Call SEM_iso_dir(n)
      call profiler_stop("SEM:iso_dir")
      ! Collect bulk velocities
      !
      call mpi_allreduce(MPI_IN_PLACE,bulk_vel_diff(1,fam) ,1,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      call mpi_allreduce(MPI_IN_PLACE,bulk_vel_diff(2,fam) ,1,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      call mpi_allreduce(MPI_IN_PLACE,bulk_vel_diff(3,fam) ,1,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      !
      in_area = inlet_area(fam)
      !
      bulk_vel_diff(1,fam)=u0(1,fam) - (bulk_vel_diff(1,fam)/in_area)
      bulk_vel_diff(2,fam)=u0(2,fam) - (bulk_vel_diff(2,fam)/in_area)
      bulk_vel_diff(3,fam)=u0(3,fam) - (bulk_vel_diff(3,fam)/in_area)
      ! Correct the inflow BC
      call profiler_start("SEM:bulk_vel")
      CALL SEM_bulk_vel(n)
      call profiler_stop("SEM:bulk_vel")
    endif
    RETURN
  END SUBROUTINE synthetic_eddies

!=======================================================================
!     Convect eddies and recycle by regenerating the location
!=======================================================================

  SUBROUTINE advect_recycle_eddies(n,dt)
    !
    IMPLICIT NONE

    INTEGER, INTENT(in)   :: n
    real(rp), intent(in ) :: dt
    INTEGER,PARAMETER     :: fam = 1
    INTEGER               :: i

       DO i=1,n/2
          eddy(i,3) = eddy(i,3) + z0(fam)*dt
          !     ---- recycle exiting eddies ----
          IF (z0(fam) > 0) THEN
             IF (eddy(i,3).GT.tmax )THEN
                CALL gen_eddy(i,fam)
             ENDIF
          ELSE
             IF (eddy(i,3) .LT. tmin )THEN
                CALL gen_eddy(i,fam)
             ENDIF
          ENDIF
       ENDDO
       !the rest is for the periodicity:
       eddy((  neddy/2+1):neddy,:) = eddy(1:neddy/2,:)
       DO i=n/2+1,n
          IF (eddy(i,1).GT.(0.5_rp*(rmax+rmin))) THEN
            eddy(i,1) = eddy(i,1) - (rmax-rmin)
          ELSE
            eddy(i,1) = eddy(i,1) + (rmax-rmin)
          ENDIF
       ENDDO
    RETURN
  END SUBROUTINE advect_recycle_eddies



!=======================================================================
!     Generate eddy location randomly in bounding box 
!    for SEM (Poletto et al.)  (ONLY on rank=0)
!=======================================================================

  SUBROUTINE gen_eddy(n,fam)
    !use ziggurat, only: rnd_loc
    !
    IMPLICIT NONE


    INTEGER, INTENT(in) :: n,fam
    REAL(rp)            :: rnd, rho, theta
    INTEGER             :: j
    !
    !
    if(myid.EQ.0) THEN
      !
      call random_number(eddy(n,:))
      !
      IF (z0(fam) > 0) THEN
        eddy(n,3) = tmin
      ELSE
        eddy(n,3) = tmax
      ENDIF
      !
      eddy(n,1) =  rmin + eddy(n,1)*(rmax-rmin)
      eddy(n,2) =  smin + eddy(n,2)*(smax-smin)
      !
      DO j=1,3
         IF (eddy(n,3+j).GT.0.5) eddy(n,3+j) =  1
         IF (eddy(n,3+j).LE.0.5) eddy(n,3+j) = -1
      ENDDO
      !
    ENDIF
    !     Broadcast eddy to each processor
    !
    CALL mpi_bcast(eddy(n,:),6,MPI_REAL_RP,0,MPI_COMM_WORLD,ierr)
    !
    RETURN
  END SUBROUTINE gen_eddy


!=======================================================================
!    Generate eddies
!=======================================================================

  SUBROUTINE gen_eddy_all()
    use mod_param, only: dx
    !use mod_random, only: random_value
    !use ziggurat, only: rnd_loc,uni
    IMPLICIT NONE

    INTEGER,PARAMETER :: fam = 1
    REAL(rp)          :: rnd
    INTEGER           :: j,n,ns,root,max_id,process

    !     Generate uniformly distributed random locations
    !     inside BLOCK for SEM of (Poletto et al.) (only on rank-0)
    !
    !
    eddy(:,:) = 0.0
    !  
    call random_number(eddy)  
    !
    root = 0
    !     Broadcast eddy to each processor
    CALL mpi_bcast(eddy(1,1),6*neddy,MPI_REAL_RP,root,MPI_COMM_WORLD,ierr)
    !
    DO n=1,neddy/2
        !
        eddy(n,1) =  rmin + eddy(n,1)*(rmax-rmin)
        eddy(n,2) =  smin + eddy(n,2)*(smax-smin)
        eddy(n,3) =  tmin + eddy(n,3)*(tmax-tmin)
        IF(eddy(n,4).LE.0.5) THEN
          eddy(n,4) = - 1.0
        ELSE
          eddy(n,4) =   1.0
        ENDIF
        
        IF(eddy(n,5).LE.0.5) THEN
          eddy(n,5) = - 1.0
        ELSE
          eddy(n,5) =   1.0
        ENDIF
        
        IF(eddy(n,6).LE.0.5) THEN
          eddy(n,6) = - 1.0
        ELSE
          eddy(n,6) =   1.0
        ENDIF
    ENDDO   
    !the rest is for the periodicity:
    eddy((  neddy/2+1):neddy,:) = eddy(1:neddy/2,:)
    DO n=neddy/2+1,neddy
      IF (eddy(n,1).GT.(0.5_rp*(rmax+rmin))) THEN
        eddy(n,1) = eddy(n,1) - (rmax-rmin)
      ELSE
        eddy(n,1) = eddy(n,1) + (rmax-rmin)
      ENDIF
    ENDDO
    !
    RETURN
    !
  END SUBROUTINE gen_eddy_all
!=======================================================================
!     RETURN velocity/turbulence kinetic energy/dissipation rate
!     at given 1D-position, via linear interpolation.cd fft 
!
!     It is used as normal distance to the nearest wall 
!=======================================================================

  SUBROUTINE SEMinputData(pos_in, vel_out, tke_out, eps_out,fam)

    IMPLICIT NONE

    REAL(rp), INTENT(in ) :: pos_in
    REAL(rp), INTENT(out) :: vel_out, tke_out, eps_out

    INTEGER i,fam
    REAL a1, a2

    DO i=1,nlines(fam)-1
       IF (pos_in.LE.pos(i+1,fam)) THEN
          a1 = (pos(i+1,fam)-pos_in)/(pos(i+1,fam)-pos(i,fam))
          a2 = (pos_in-pos(i,fam))/(pos(i+1,fam)-pos(i,fam))
          vel_out = a1 * umean(i,fam) + a2 * umean(i+1,fam)
          tke_out = a1 * tke(i,fam) + a2 * tke(i+1,fam)
          eps_out = a1 * dissip(i,fam) + a2 * dissip(i+1,fam)
          RETURN
       ENDIF
    ENDDO

    vel_out = umean(nlines(fam),fam)
    tke_out = tke(nlines(fam),fam)
    eps_out = dissip(nlines(fam),fam)

    RETURN
  END SUBROUTINE SEMinputData

!=======================================================================
!    Interpolate inlet data
!=======================================================================

    SUBROUTINE SEM_u0_norm(fam,dir,n)
    use mod_param    ,only: IBM_h,dx,dy,dz
    use mod_common_mpi, only: MPI_PROC_NULL,left,right,front,back,top,bottom
    IMPLICIT NONE
    !
    integer, INTENT(in)                     :: dir
    integer , intent(in ), dimension(3)     :: n
    REAL(rp)                                :: vel_interp, tke_interp
    REAL(rp)                                :: dist,sigmal
    REAL(rp)                                :: dissip_interp
    INTEGER                                 :: ii,e,test,fam,i,j,k
    !
    select case(dir)
    case(2)
      do i=1,n(1)
        do k=1 + IBM_h,n(3)
          dist = (k-IBM_h)*dz-0.5*dz   !distance from wall
          CALL SEMinputData(dist,vel_interp,tke_interp, dissip_interp,fam)
          if(front.eq.MPI_PROC_NULL) u0(2,fam)=u0(2,fam)+vel_interp*(dz*dx)
          ! Implementing formula 8.2.3 of Poletto's PhD here:
          sigmal = (tke_interp**1.5)/dissip_interp
          !
          sigmal = MIN(sigmal,0.41*max_range(fam))
          ! Limit eddy size far away from wall
          ! to Prandtl's mixing length
          ! kappa is .41, max_range is max the distance to the wall
          sigmal = MIN(sigmal,dist)
          ! make eddies no larger than box 
          sigmal = MAX(sigmal, dz)
          ! avoid creating too small eddies
          sigmal = MIN(sigmal,sigma_max(fam))
          !
          sigma_inlet(i,k,fam)     = sigmal
          umean_inlet(i,k,fam)     = vel_interp
          intensity_inlet(i,k,fam) = SQRT(2./3.*tke_interp)
        ENDDO
      enddo
    case(0)
      ! Do nothing
    case default
      if(myid.eq.0) print*, 'ERROR: This SEM direction is not implemented ... yet!'
      call MPI_FINALIZE(ierr)
      call exit
    end select
    !
    RETURN
  END SUBROUTINE SEM_u0_norm
!
!=======================================================================
!    do restart routine
!=======================================================================
!
    SUBROUTINE SEMrestart(action)

    IMPLICIT NONE

    INTEGER                 :: i, j, max_fac, neddy_tmp
    character(1), intent(in):: action
    INTEGER                 :: max_dum,file_id
    CHARACTER(80)           :: readfile
    LOGICAL                 :: file_exist
    INTEGER,PARAMETER       :: fam = 1

    IF (action.EQ.'R') THEN
          ! read in eddy data from files:
          ! ex,ey,ez, eps
          IF (myid.EQ.0) THEN
                INQUIRE (FILE=restart_file, EXIST = file_exist)
                IF (file_exist) THEN
                print*, 'Reading: ', restart_file
                   OPEN(newunit=file_id,form='unformatted',file=restart_file,STATUS='OLD',action='read')
                   READ(file_id) neddy_tmp

                   IF (.NOT.neddy_tmp.EQ.neddy) THEN
                     if(myid.eq.0) print*, 'SEM ERROR: ABORT Trying to restart from restart file with different nEddy'
                     call MPI_FINALIZE(ierr)
                     call exit
                   ENDIF
                   !
                   do j = 1,6
                     do i=1,neddy
                       READ(file_id) eddy(i,j)
                     enddo
                   ENDDO
                   !
                   CLOSE(file_id)
                ENDIF
                !
          ENDIF
          call mpi_barrier(MPI_COMM_WORLD,ierr)
          !     Broadcast eddy to each processor
          CALL mpi_bcast(eddy(1,1),6*neddy,MPI_REAL_RP,0,MPI_COMM_WORLD,ierr)
    ELSEIF ((action.EQ.'W').AND.(dir_xyz(fam).NE.0)) THEN
      !
      ! dump restart files
      ! Data is gathered already
      !
      IF (myid.EQ.0) THEN
          INQUIRE (FILE=restart_file, EXIST = file_exist)
          IF (file_exist) THEN
            print*, 'Writing restart file: ', restart_file
            OPEN(newunit=file_id,form='unformatted',file=restart_file,action='write',status='replace')
          ELSE
            print*, 'Writing restart file: ', restart_file
            OPEN(newunit=file_id,form='unformatted',file=restart_file,action='write',status='new')
          ENDIF
          WRITE(file_id) neddy
          !
          do j = 1,6
            do i=1,neddy
              WRITE(file_id) eddy(i,j)
            enddo
          ENDDO
          !
          CLOSE (file_id)
      ENDIF
    ENDIF

    RETURN

  END SUBROUTINE SEMrestart

!
!=======================================================================
!    Calculate the boundary condition
!=======================================================================
!
  SUBROUTINE SEM_iso_dir(n)
    !
    use mod_common_mpi, only: MPI_PROC_NULL,left,right,front,back,top,bottom
    use mod_param     , only: dx,dy,dz,rho2,cbcvel,IBM_h,itot,jtot,ktot,dims_in
    IMPLICIT NONE

    INTEGER,PARAMETER                                            :: fam = 1
    integer  , intent(in   ), dimension(3)                       :: n 
    INTEGER                                                      :: e,i,j,k,ne,ii,receiver,start,in
    real(rp), PARAMETER                                          :: sqrt32 = SQRT(3./2)
    real(rp)                                                     :: xs,ys,zs,xmin,xmax
    real(rp)                                                     :: intensity,sigma,sigma_i,ff_start,max_sig
    real(rp)                                                     :: ff,fx,fy,fz,rr,rrx,rry,rrz,work,in_area
    real(rp)                                                     :: rho_factor,sqrtn
    real(rp) ,ALLOCATABLE                                        :: work_array(:,:)
    real(rp), PARAMETER                                          :: pi = 3.14159265358979323846264338327950288_rp
    real(rp), parameter                                          :: small_num = 1e-12_rp
    real(rp),allocatable,dimension(:,:,:)                        :: vel_bc_tot,buffer
    integer                                                      :: lb1,lb2,lb3
    !     functions
    ! ---- compute velocity contribution of eddies ------
    !
    bulk_vel_diff(:,fam) = 0
    !
    sqrtn = SQRT(1./REAL(neddy))
    !
    select case(dir_xyz(fam))
    case(2)
      allocate(vel_bc_tot(1:n(1),1:n(3),3))
      allocate(buffer(1:n(1),1:n(3),3))
      vel_bc_tot = 0.0_rp
      buffer = 0.0_rp
      max_sig = maxval(sigma_inlet)
      xmin = (     ijk_start(1))*dx
      xmax = (n(1)+ijk_start(1))*dx
      !
      DO ne=start_edd,stop_edd
        !==================================================
        if(.NOT.(((eddy(ne,1)-xmax).GT.max_sig).OR.((xmin-eddy(ne,1)).GT.max_sig))) then
        !
        do k=1 + IBM_h,n(3)
          zs        = (k-0.5_rp)*dz
          rrz       = abs(zs-eddy(ne,2))     !y is wall normal
          sigma     =  sigma_inlet(1,k,fam)
          IF (rrz.LT.sigma) THEN
            intensity =  intensity_inlet(1,k,fam)
            sigma_i   =  1.0_rp/sigma
            ff_start  = SQRT(16._rp*Vb(fam)/(15._rp*pi*sigma**3.0))
            do i=1,n(1)
              ! first for the x-direction
              ys = -0.5*dy
              xs = (i+ijk_start(1))*dx
              rrx = (xs-eddy(ne,1))     !x is span
              rry = (ys-eddy(ne,3))     !z is streamwise
              rr = SQRT(rrx*rrx + rry*rry + rrz*rrz)
              !
              IF (rr .LT. sigma) THEN
                fx = rry*eddy(ne,6)-rrz*eddy(ne,5)
                !
                fx = fx * intensity*sigma_i
                !
                ff = ff_start * (SIN(pi*rr*sigma_i)**2.0)*rr*sigma_i
                ff = ff/((rr*sigma_i)**3.0)
                ff = ff*sqrtn
                !
                vel_bc_tot(i,k,1)    = vel_bc_tot(i,k,1) + fx*ff
                bulk_vel_diff(1,fam) = bulk_vel_diff(1,fam) + dx*dz*fx*ff
              ENDIF
              ! now for the y-direction
              ys  = 0.0
              xs  = (i+ijk_start(1)-0.5)*dx
              rrx = (xs-eddy(ne,1))     !x is span
              rry = (ys-eddy(ne,3))     !z is streamwise
              rr  = SQRT(rrx*rrx + rry*rry + rrz*rrz)
              !
              IF (rr .LT. sigma) THEN
                fy = rrz*eddy(ne,4)-rrx*eddy(ne,6)
                !
                fy = fy * intensity*sigma_i
                !
                ff = ff_start * (SIN(pi*rr*sigma_i)**2.0)*rr*sigma_i
                ff = ff/((rr*sigma_i)**3.0)
                ff = ff*sqrtn
                !
                vel_bc_tot(i,k,2) = vel_bc_tot(i,k,2) + fy*ff
              ENDIF
            enddo
          endif
        enddo
        !
        do k=1 + IBM_h,n(3)
          zs = k*dz
          ii = 0
          if(k.eq.n(3)) ii = 1
          rrz = abs(zs-eddy(ne,2))     !y is wall normal
          sigma     =  0.5*(sigma_inlet(1,k,fam)+sigma_inlet(1,k+1-ii,fam))
          IF (rrz.LT.sigma) THEN
            intensity =  0.5*(intensity_inlet(1,k,fam)+intensity_inlet(1,k+1-ii,fam))
            sigma_i   =  1.0_rp/sigma
            ff_start  = SQRT(16._rp*Vb(fam)/(15._rp*pi*sigma**3.0))
            do i=1,n(1)
              ! now for the z-direction
              ys = -0.5*dy
              xs  = (i+ijk_start(1)-0.5)*dx
              rrx = (xs-eddy(ne,1))     !x is span
              rry = (ys-eddy(ne,3))     !z is streamwise
              rr = SQRT(rrx*rrx + rry*rry + rrz*rrz)
              !
              IF (rr .LT. sigma) THEN
                fz = rrx*eddy(ne,5)-rry*eddy(ne,4)
                !
                fz = fz * intensity*sigma_i
                !
                ff = ff_start * (SIN(pi*rr*sigma_i)**2.0)*rr*sigma_i
                ff = ff/((rr*sigma_i)**3.0)
                ff = ff*sqrtn
                !
                vel_bc_tot(i,k,3) = vel_bc_tot(i,k,3) + fz*ff
                bulk_vel_diff(3,fam)=bulk_vel_diff(3,fam) + dx*dz*fz*ff
              ENDIF
            enddo
          endif
        enddo
        endif
        !=================================~
      ENDDO
      !
      if((front.eq.MPI_PROC_NULL).AND.(cbcvel(0,2,1).EQ.'D') &
                                 .AND.(cbcvel(0,2,2).EQ.'D').AND.(cbcvel(0,2,3).EQ.'D')) then
        !
        u_bc               = 0.0_rp
        v_bc               = 0.0_rp
        w_bc               = 0.0_rp
      endif
      call profiler_start("SEM_COM", tag = .true., tag_color = COLOR_WHITE)
      !call mpi_allreduce(MPI_IN_PLACE,u_bc_tot ,itot*ktot,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      !call mpi_allreduce(MPI_IN_PLACE,v_bc_tot ,itot*ktot,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      !call mpi_allreduce(MPI_IN_PLACE,w_bc_tot ,itot*ktot,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      !do in = 1,dims_in(1)
      !  receiver   = receiver_rank(in)
      !  start      = start_address(in)
      !  call mpi_reduce(u_bc_tot((1+start):(n(1)+start),1:n(3)),u_bc,n(1)*n(3), &
      !                                     MPI_REAL_RP,MPI_SUM,receiver,MPI_COMM_WORLD,ierr)
      !  call mpi_reduce(v_bc_tot((1+start):(n(1)+start),1:n(3)),v_bc,n(1)*n(3), &
      !                                     MPI_REAL_RP,MPI_SUM,receiver,MPI_COMM_WORLD,ierr)
      !  call mpi_reduce(w_bc_tot((1+start):(n(1)+start),1:n(3)),w_bc,n(1)*n(3), &
      !                                     MPI_REAL_RP,MPI_SUM,receiver,MPI_COMM_WORLD,ierr)
      !enddo
      !
      do in = 1,dims_in(2)
        if(sender_ranks(in).eq.receiver_rank) then
          if(myid.eq.receiver_rank) then
            u_bc(:,:) = u_bc(:,:) + vel_bc_tot(:,:,1)
            v_bc(:,:) = v_bc(:,:) + vel_bc_tot(:,:,2)
            w_bc(:,:) = w_bc(:,:) + vel_bc_tot(:,:,3)
          endif
        else
          if(myid.eq.sender_ranks(in)) then
            lb1 = lbound(vel_bc_tot,1) 
            lb2 = lbound(vel_bc_tot,2)
            lb3 = lbound(vel_bc_tot,3)
            call MPI_SEND(vel_bc_tot(lb1,lb2,lb3), size( vel_bc_tot ),MPI_REAL_RP,receiver_rank ,in, &
                          comm_cart,ierr)
          endif
          !
          if(myid.eq.receiver_rank) then
            lb1 = lbound(buffer,1) 
            lb2 = lbound(buffer,2)
            lb3 = lbound(buffer,3)
            call MPI_RECV(buffer(lb1,lb2,lb3), size( buffer ),MPI_REAL_RP,sender_ranks(in) ,in, &
                          comm_cart,status,ierr)
            u_bc(:,:) = u_bc(:,:) + buffer(:,:,1)
            v_bc(:,:) = v_bc(:,:) + buffer(:,:,2)
            w_bc(:,:) = w_bc(:,:) + buffer(:,:,3)
          endif
        endif
      enddo
      !
      call profiler_stop("SEM_COM")
      !
      if((front.eq.MPI_PROC_NULL).AND.(cbcvel(0,2,1).EQ.'D') &
                                 .AND.(cbcvel(0,2,2).EQ.'D').AND.(cbcvel(0,2,3).EQ.'D')) then
       !
        !u_bc               = 0.0_rp
        !v_bc               = 0.0_rp
        !w_bc               = 0.0_rp
        !u_bc(1:n(1),1:n(3)) = u_bc_tot((1+ijk_start(1)):(n(1)+ijk_start(1)),(1+ijk_start(3)):(n(3)+ijk_start(3)))
        !v_bc(1:n(1),1:n(3)) = v_bc_tot((1+ijk_start(1)):(n(1)+ijk_start(1)),(1+ijk_start(3)):(n(3)+ijk_start(3)))
        !w_bc(1:n(1),1:n(3)) = w_bc_tot((1+ijk_start(1)):(n(1)+ijk_start(1)),(1+ijk_start(3)):(n(3)+ijk_start(3)))
        !!!!!!
        !add mean velocity for y
        do k=1 + IBM_h,n(3)
          do i=1,n(1)
            v_bc(i,k) = v_bc(i,k) + umean_inlet(i,k,fam)  !v_bc_tot(i+ijk_start(1),k+ijk_start(3))
            !u_bc(i,k) = u_bc_tot(i+ijk_start(1),k+ijk_start(3))
            !w_bc(i,k) = w_bc_tot(i+ijk_start(1),k+ijk_start(3))
            bulk_vel_diff(2,fam)=bulk_vel_diff(2,fam) + dx*dz*v_bc(i,k)
          enddo
        enddo
      endif
      deallocate(buffer,vel_bc_tot)
    case(0)
      ! Do nothing
    case default
      if(myid.eq.0) print*, 'ERROR: This SEM direction is not implemented ... yet!'
      call MPI_FINALIZE(ierr)
      call exit
    end select

    RETURN
  END SUBROUTINE SEM_iso_dir
!
!=======================================================================
!    Correct the bulk velocity
!=======================================================================
!
  SUBROUTINE SEM_bulk_vel(n)
    !
    use mod_common_mpi, only: MPI_PROC_NULL,left,right,front,back,top,bottom
    use mod_param     , only: cbcvel,IBM_h
    IMPLICIT NONE

    integer  , intent(in   ), dimension(3)                       :: n 
    INTEGER,PARAMETER                                            :: fam = 1
    INTEGER                                                      :: e,i,j,k,ne,ii

    select case(dir_xyz(fam))
    case(2)
      if((front.eq.MPI_PROC_NULL).AND.(cbcvel(0,2,1).EQ.'D') &
                                 .AND.(cbcvel(0,2,2).EQ.'D').AND.(cbcvel(0,2,3).EQ.'D')) then
        do k=1 + IBM_h,n(3)
          do i=1,n(1)
            u_bc(i,k) =  u_bc(i,k)  + bulk_vel_diff(1,fam)
            v_bc(i,k)  = v_bc(i,k)  + bulk_vel_diff(2,fam)
            w_bc(i,k)  = w_bc(i,k)  + bulk_vel_diff(3,fam)
          enddo
        enddo
      endif
    case(0)
      ! Do nothing
    case default
      if(myid.eq.0) print*, 'ERROR: This SEM direction is not implemented ... yet!'
      call MPI_FINALIZE(ierr)
      call exit
    end select

    RETURN
  END SUBROUTINE SEM_bulk_vel
!
!=======================================================================
!    Apply the boundary condition
!=======================================================================
!
  SUBROUTINE SEM_apply(n,nh_u,u,v,w,rho_switch)
    !
    use mod_common_mpi, only: MPI_PROC_NULL,left,right,front,back,top,bottom
    use mod_param     , only: dx,dy,dz,rho2,cbcvel,IBM_h
    IMPLICIT NONE

    INTEGER,PARAMETER                                            :: fam = 1
    integer  , intent(in   ), dimension(3)                       :: n 
    integer  , intent(in   )                                     :: nh_u
    real(rp) , intent(inout), dimension(1-nh_u:,1-nh_u:,1-nh_u:) :: u,v,w
    logical  , intent(in   )                                     :: rho_switch
    INTEGER                                                      :: e,i,j,k,ne
    real(rp)                                                     :: rho_factor
    !
    If(rho_switch) then
        rho_factor = rho2   !multiply py denity of the liquid phase 
    else 
        rho_factor = 1.0_rp    ! velocity stays velocity
    endif
    !
    !
    select case(dir_xyz(fam))
    case(2)
      if((front.eq.MPI_PROC_NULL).AND.(cbcvel(0,2,1).EQ.'D') &
                                 .AND.(cbcvel(0,2,2).EQ.'D').AND.(cbcvel(0,2,3).EQ.'D')) then
        ! first for the x-direction
        do k=1 + IBM_h,n(3)
          do i=1,n(1)
            u(i,0,k) = u_bc(i,k)*rho_factor
          enddo
        enddo
        ! now for the y-direction
        do k=1 + IBM_h,n(3)
          do i=1,n(1)
            v(i,0,k) = v_bc(i,k)*rho_factor
          enddo
        enddo
        ! now for the z-direction
        do k=1 + IBM_h,n(3)
          do i=1,n(1)
            w(i,0,k) = w_bc(i,k)*rho_factor
          enddo
        enddo
      endif
    case(0)
      ! Do nothing
    case default
      if(myid.eq.0) print*, 'ERROR: This SEM direction is not implemented ... yet!'
      call MPI_FINALIZE(ierr)
      call exit
    end select

    RETURN
  END SUBROUTINE SEM_apply

  END MODULE mod_sem 
