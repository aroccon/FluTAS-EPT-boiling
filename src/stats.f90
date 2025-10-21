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
!
  MODULE mod_stats
    !
    use mod_types
    use mpi
    use mod_common_mpi, only: myid,ierr,ijk_start,comm_cart,status,n_x,n_y,n_z
    use iso_c_binding , only: C_PTR
    use mod_param     , only: ng
    !
    IMPLICIT NONE
    SAVE
    ! Spanwise averaged quantities
    real(rp), allocatable, dimension(:,:)   :: u1_2,u2_2,v1_2,v2_2,w1_2,w2_2
    real(rp), allocatable, dimension(:,:)   :: uu1_2,uu2_2,vv1_2,vv2_2,ww1_2,ww2_2
    real(rp), allocatable, dimension(:,:)   :: temp_2,pc_2,phi_2,pr1,pr2
    ! wall quantities
    real(rp), allocatable, dimension(:,:)   :: phi_inst,phi_avg
    real(rp), allocatable, dimension(:,:)   :: temp_inst,temp_avg,tmptmp_avg
    
    real(rp)                                :: start_time,timespan
    integer                                 :: receiver_rank
    integer, allocatable, dimension(:)      :: sender_ranks
    
    contains 
!=======================================================================
! setup subroutine
!=======================================================================
    SUBROUTINE stats_setup(n,time)
      !
      use mod_param, only: dims_in
      !
      IMPLICIT NONE
      !
      integer                                :: inx,iny,start_testy,start_testx,rank
      real(rp), intent(in )                  :: time
      integer , intent(in ), dimension(3)    :: n
      !
      ! set start time for averaging
      start_time = time
      timespan   = 0.0      !
      receiver_rank = 0
      allocate(sender_ranks(1:dims_in(1)))
      sender_ranks(:) = 0
      !
      ! determine receiver rank
      do iny = 1,dims_in(2)
        rank       = 0
        start_testy = (iny-1)*n(2)
        if ((ijk_start(2).EQ.start_testy).AND.(ijk_start(1).EQ.0)) then   ! i.e. the rank everyone wants to send their data to
          rank = myid
        endif
        call mpi_allreduce(MPI_IN_PLACE,rank ,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ierr)
        if (ijk_start(2).EQ.start_testy) then   ! all the ranks in that row
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
          if (ijk_start(2).EQ.start_testy) then   ! all the ranks in that row
            sender_ranks(inx) = rank
          endif
        enddo
      enddo
      if(myid.eq.receiver_rank) print*,'2D stats receiver:' ,receiver_rank, ' senders: ', sender_ranks
      !
      ! allocate spane-wise avg. arrays
      allocate(u1_2(1:n(2),1:n(3)), &
               u2_2(1:n(2),1:n(3)), &
               v1_2(1:n(2),1:n(3)), &
               v2_2(1:n(2),1:n(3)), &
               w1_2(1:n(2),1:n(3)), &
               w2_2(1:n(2),1:n(3)), &
              uu1_2(1:n(2),1:n(3)), &
              uu2_2(1:n(2),1:n(3)), &
              vv1_2(1:n(2),1:n(3)), &
              vv2_2(1:n(2),1:n(3)), &
              ww1_2(1:n(2),1:n(3)), &
              ww2_2(1:n(2),1:n(3)), &
             temp_2(1:n(2),1:n(3)), &
               pc_2(1:n(2),1:n(3)), &
              phi_2(1:n(2),1:n(3)), &
                pr1(1:n(2),1:n(3)), &
                pr2(1:n(2),1:n(3)))
      !
      ! allocate wall arrays
      allocate(phi_inst(1:n(1),1:n(2))   &
               ,phi_avg(1:n(1),1:n(2))    &
             ,temp_inst(1:n(1),1:n(2))  &
              ,temp_avg(1:n(1),1:n(2))  &
            ,tmptmp_avg(1:n(1),1:n(2)))
      !
      ! set then to zero
      u1_2(:,:)       = 0.0_rp
      u2_2(:,:)       = 0.0_rp
      v1_2(:,:)       = 0.0_rp
      v2_2(:,:)       = 0.0_rp
      w1_2(:,:)       = 0.0_rp
      w2_2(:,:)       = 0.0_rp
      uu1_2(:,:)      = 0.0_rp
      uu2_2(:,:)      = 0.0_rp
      vv1_2(:,:)      = 0.0_rp
      vv2_2(:,:)      = 0.0_rp
      ww1_2(:,:)      = 0.0_rp
      ww2_2(:,:)      = 0.0_rp
      temp_2(:,:)     = 0.0_rp
      pc_2(:,:)       = 0.0_rp
      phi_2(:,:)      = 0.0_rp
      phi_inst(:,:)   = 0.0_rp
      phi_avg(:,:)    = 0.0_rp
      temp_inst(:,:)  = 0.0_rp
      temp_avg(:,:)   = 0.0_rp
      tmptmp_avg(:,:) = 0.0_rp
      pr1(:,:)        = 0.0_rp
      pr2(:,:)        = 0.0_rp
      !
      RETURN
    END SUBROUTINE stats_setup

!=======================================================================
! Output x2x3-field
!=======================================================================

    SUBROUTINE output23field(n,field,name,datadir,time,istep)
      !
      use mod_param, only: dims_in
      use mod_output,only: write_visu_2d
      !
      IMPLICIT NONE
      !
      integer , intent(in ), dimension(3)             :: n
      real(rp), intent(in ), dimension(1:n(2),1:n(3)) :: field
      integer                                         :: lb1,lb2,in,j,k,i,rank
      real(rp), allocatable, dimension(:,:)           :: buffer,span_avg
      real(rp), allocatable, dimension(:,:,:)         :: output_field
      character(len=4),intent(in )                    :: name
      character(len=100),intent(in )                  :: datadir
      real(rp), intent(in )                           :: time
      integer , intent(in )                           :: istep
      character(len=9)                                :: fldnum
      !
      write(fldnum,'(i9.9)') istep
      if(myid.eq.receiver_rank) then
        allocate(buffer(1:n(2),1:n(3)),span_avg(1:n(2),1:n(3)))
        buffer(:,:)   = 0.0_rp
        span_avg(:,:) = 0.0_rp
      endif
      !
      do in = 1,dims_in(1)
        if(sender_ranks(in).eq.receiver_rank) then
          if(myid.eq.receiver_rank) then
            span_avg(:,:) = span_avg(:,:) + field(:,:)
          endif
        else
          if(myid.eq.sender_ranks(in)) then
            lb1 = lbound(field,1) 
            lb2 = lbound(field,2)
            call MPI_SEND(field(lb1,lb2), size( field ),MPI_REAL_RP,receiver_rank ,in, &
                          comm_cart,ierr)
          endif
          !
          if(myid.eq.receiver_rank) then
            lb1 = lbound(buffer,1) 
            lb2 = lbound(buffer,2)
            call MPI_RECV(buffer(lb1,lb2), size( buffer ),MPI_REAL_RP,sender_ranks(in) ,in, &
                          comm_cart,status,ierr)
            span_avg(:,:) = span_avg(:,:) + buffer(:,:)
          endif
        endif
      enddo
      !
      allocate(output_field(1:n(1),1:n(2),1:n(3)))
      output_field = 0.0_rp
      !
      if(myid.eq.receiver_rank) then
        deallocate(buffer)
        span_avg(:,:) = span_avg(:,:)/real(dims_in(1))
        i = 1
        do j = 1,n(2)
          do k = 1,n(3)
            output_field(i,j,k) = span_avg(j,k)
          enddo
        enddo
        deallocate(span_avg)
      endif
      call write_visu_2d(trim(datadir),name//'_'//fldnum//'.bin','log_visu_2d_slice_x2x3.out',name, &
                     1,1,ng,time,istep,output_field(1:n(1),1:n(2),1:n(3)))
      deallocate(output_field)
      RETURN
    END SUBROUTINE output23field
    !
!=======================================================================
! Output x1x2-field on the wall hight
!=======================================================================

    SUBROUTINE outputwallfield(n,field,name,datadir,time,istep)
      !
      use mod_output,only: write_visu_2d
      !
      IMPLICIT NONE
      !
      integer , intent(in ), dimension(3)             :: n
      real(rp), intent(in ), dimension(1:n(1),1:n(2)) :: field
      integer                                         :: j,k,i
      real(rp), allocatable, dimension(:,:,:)         :: output_field
      character(len=4),intent(in )                    :: name
      character(len=100),intent(in )                  :: datadir
      real(rp), intent(in )                           :: time
      integer , intent(in )                           :: istep
      character(len=9)                                :: fldnum
      !
      write(fldnum,'(i9.9)') istep
      !
      allocate(output_field(1:n(1),1:n(2),1:n(3)))
      output_field = 0.0_rp
      !
      k = 1
      do j = 1,n(2)
        do i = 1,n(1)
          output_field(i,j,k) = field(i,j)
        enddo
      enddo
      !
      call write_visu_2d(trim(datadir),name//'_'//fldnum//'.bin','log_visu_2d_slice_wall.out',name, &
                     3,1,ng,time,istep,output_field(1:n(1),1:n(2),1:n(3)))
      !
      deallocate(output_field)
      RETURN
    END SUBROUTINE outputwallfield
    !
!=======================================================================
! Collect the statistics
!=======================================================================

    SUBROUTINE collect_statistics(n,nh_p,nh_u,time,istep,phi,tmp,cond,src_ph,p,u,v,w)
      !
      use mod_param,   only: IBM_h,iout2d,heater_power,dl
      use mod_boiling, only: heater_ind,rex,rey,rez
      !
      IMPLICIT NONE
      !
      integer , intent(in ), dimension(3)                       :: n
      integer                                                   :: j,k,i,ip,jp,kp
      real(rp), intent(in )                                     :: time
      integer , intent(in )                                     :: istep,nh_p,nh_u
      real(rp), intent(in ), dimension(1-nh_p:,1-nh_p:,1-nh_p:) :: phi,tmp,cond,src_ph,p
      real(rp), intent(in ), dimension(1-nh_u:,1-nh_u:,1-nh_u:) :: u,v,w
      real(rp)                                                  :: alpha,beta,dtime
      real(rp)                                                  :: T_wall,phix,phiy,phiz
      real(rp)                                                  :: u1,v1,w1,u2,v2,w2
      real(rp)                                                  :: threshold

      !
      threshold = 1e-4   ! Phasic velocity threshold
      !
      if((mod(istep,25).EQ.0).or.(mod(istep,iout2d).eq.0)) then      
        dtime = time - timespan - start_time
        timespan = time - start_time
        !
        ! Time average is compuated as:
        ! Accumulated=alpha*Accumulated+beta*New
        ! Calculate alpha and beta
        beta = dtime/timespan
        alpha = 1.0_rp - beta
        !
        if(myid.EQ.0) print*,'Perform time-averaging...'
        !===============================================
        ! spanwise-averaged phase field value
        !===============================================
        phi_2(:,:) = phi_2(:,:)*alpha
        do k = 1+IBM_h,n(3)
          do j = 1,n(2)
            do i = 1,n(1)
              phi_2(j,k) = phi_2(j,k) + beta*(1.0_rp/real(n(1),rp))*phi(i,j,k)
            enddo
          enddo
        enddo
        !===============================================
        ! 1. time averaged wall temperature
        ! 2. instantanous wall temperature
        ! 3. time avergaed wall temperature squared
        !===============================================
        temp_avg(:,:)   = temp_avg(:,:)*alpha
        tmptmp_avg(:,:) = tmptmp_avg(:,:)*alpha
        k=IBM_h+1
        do j = 1,n(2)
          do i = 1,n(1)
            if(IBM_h.NE.0) then
              T_wall = (cond(i,j,k)*tmp(i,j,k)+cond(i,j,k-1)*tmp(i,j,k-1)) &
                      /(cond(i,j,k)         +cond(i,j,k-1)          ) &
                      +(heater_power*dl(3)*heater_ind(i,j))          &
                      /(2.0_rp*(cond(i,j,k)+cond(i,j,k-1)))
            else
              T_wall = 0.5_rp*(tmp(i,j,k)+tmp(i,j,k-1))
            endif
            temp_avg(i,j)  = temp_avg(i,j) + beta*T_wall
            temp_inst(i,j) = T_wall
            tmptmp_avg(i,j)  = tmptmp_avg(i,j) + beta*T_wall*T_wall
          enddo
        enddo
        !===============================================
        ! time averaged wall phase field
        !===============================================
        phi_avg(:,:) = phi_avg(:,:)*alpha
        k=IBM_h+1
        do j = 1,n(2)
          do i = 1,n(1)
            phi_avg(i,j) = phi_avg(i,j) + beta*0.5_rp*(phi(i,j,k)+phi(i,j,k-1))
          enddo
        enddo
        !===============================================
        ! instantanous wall phase field
        !===============================================
        k=IBM_h+1
        do j = 1,n(2)
          do i = 1,n(1)
            phi_inst(i,j) = 0.5_rp*(phi(i,j,k)+phi(i,j,k-1))
          enddo
        enddo
        !===============================================
        ! spanwise-averaged mass transfer density
        !===============================================
        pc_2(:,:) = pc_2(:,:)*alpha
        do k = 1+IBM_h,n(3)
          do j = 1,n(2)
            do i = 1,n(1)
              pc_2(j,k) = pc_2(j,k) + beta*(1.0_rp/real(n(1),rp))*src_ph(i,j,k)
            enddo
          enddo
        enddo
        !===============================================
        ! spanwise-averaged temperature field
        !===============================================
        temp_2(:,:) = temp_2(:,:)*alpha
        do k = 1,n(3)
          do j = 1,n(2)
            do i = 1,n(1)
              temp_2(j,k) = temp_2(j,k) + beta*(1.0_rp/real(n(1),rp))*tmp(i,j,k)
            enddo
          enddo
        enddo
        !===============================================
        ! spanwise-averaged pressure field of phase 1
        !===============================================
        pr1(:,:) = pr1(:,:)*alpha
        do k = 1+IBM_h,n(3)
          do j = 1,n(2)
            do i = 1,n(1)
              pr1(j,k) = pr1(j,k) + beta*(1.0_rp/real(n(1),rp))*p(i,j,k)*phi(i,j,k)
            enddo
          enddo
        enddo
        !===============================================
        ! spanwise-averaged temperature field
        !===============================================
        pr2(:,:) = pr2(:,:)*alpha
        do k = 1+IBM_h,n(3)
          do j = 1,n(2)
            do i = 1,n(1)
              pr2(j,k) = pr2(j,k) + beta*(1.0_rp/real(n(1),rp))*p(i,j,k)*(1.0_rp-phi(i,j,k))
            enddo
          enddo
        enddo
        !===============================================
        ! spanwise-averaged velocities (u,v,w) of phase 1
        !===============================================
        u1_2(:,:) = u1_2(:,:)*alpha
        v1_2(:,:) = v1_2(:,:)*alpha
        w1_2(:,:) = w1_2(:,:)*alpha
        do k = 1+IBM_h,n(3)
          do j = 1,n(2)
            do i = 1,n(1)
              ip    = i+1   
              jp    = j+1   
              kp    = k+1  
              phix  = 0.5_rp*(phi(ip,j ,k )+phi(i,j,k))
              phiy  = 0.5_rp*(phi(i ,jp,k )+phi(i,j,k))
              phiz  = 0.5_rp*(phi(i ,j ,kp)+phi(i,j,k))
              u1_2(j,k) = u1_2(j,k) + beta*(1.0_rp/real(n(1),rp))*(u(i,j,k)*phix-rex(i,j,k))
              v1_2(j,k) = v1_2(j,k) + beta*(1.0_rp/real(n(1),rp))*(v(i,j,k)*phiy-rey(i,j,k))
              w1_2(j,k) = w1_2(j,k) + beta*(1.0_rp/real(n(1),rp))*(w(i,j,k)*phiz-rez(i,j,k))
            enddo
          enddo
        enddo
        !===============================================
        ! spanwise-averaged velocities (u,v,w) of phase 2
        !===============================================
        u2_2(:,:) = u2_2(:,:)*alpha
        v2_2(:,:) = v2_2(:,:)*alpha
        w2_2(:,:) = w2_2(:,:)*alpha
        do k = 1+IBM_h,n(3)
          do j = 1,n(2)
            do i = 1,n(1)
              ip    = i+1   
              jp    = j+1   
              kp    = k+1  
              phix  = 0.5_rp*(phi(ip,j ,k )+phi(i,j,k))
              phiy  = 0.5_rp*(phi(i ,jp,k )+phi(i,j,k))
              phiz  = 0.5_rp*(phi(i ,j ,kp)+phi(i,j,k))
              u2_2(j,k) = u2_2(j,k) + beta*(1.0_rp/real(n(1),rp))*(u(i,j,k)*(1.0_rp-phix)+rex(i,j,k))
              v2_2(j,k) = v2_2(j,k) + beta*(1.0_rp/real(n(1),rp))*(v(i,j,k)*(1.0_rp-phiy)+rey(i,j,k))
              w2_2(j,k) = w2_2(j,k) + beta*(1.0_rp/real(n(1),rp))*(w(i,j,k)*(1.0_rp-phiz)+rez(i,j,k))
            enddo
          enddo
        enddo
        !===============================================
        ! spanwise-averaged squared velocities (uu,vv,ww) of phase 1
        !===============================================
        uu1_2(:,:) = uu1_2(:,:)*alpha
        vv1_2(:,:) = vv1_2(:,:)*alpha
        ww1_2(:,:) = ww1_2(:,:)*alpha
        do k = 1+IBM_h,n(3)
          do j = 1,n(2)
            do i = 1,n(1)
              ip    = i+1   
              jp    = j+1   
              kp    = k+1  
              phix  = 0.5_rp*(phi(ip,j ,k )+phi(i,j,k))
              phiy  = 0.5_rp*(phi(i ,jp,k )+phi(i,j,k))
              phiz  = 0.5_rp*(phi(i ,j ,kp)+phi(i,j,k))
              u1    = 0.0_rp
              v1    = 0.0_rp
              w1    = 0.0_rp
              if(phix.GT.threshold) u1 = u(i,j,k)-rex(i,j,k)/phix
              if(phiy.GT.threshold) v1 = v(i,j,k)-rey(i,j,k)/phiy
              if(phiz.GT.threshold) w1 = w(i,j,k)-rez(i,j,k)/phiz
              uu1_2(j,k) = uu1_2(j,k) + beta*(1.0_rp/real(n(1),rp))*(u(i,j,k)*phix-rex(i,j,k))*u1
              vv1_2(j,k) = vv1_2(j,k) + beta*(1.0_rp/real(n(1),rp))*(v(i,j,k)*phiy-rey(i,j,k))*v1
              ww1_2(j,k) = ww1_2(j,k) + beta*(1.0_rp/real(n(1),rp))*(w(i,j,k)*phiz-rez(i,j,k))*w1
            enddo
          enddo
        enddo
        !===============================================
        ! spanwise-averaged squared velocities (uu,vv,ww) of phase 2
        !===============================================
        uu2_2(:,:) = uu2_2(:,:)*alpha
        vv2_2(:,:) = vv2_2(:,:)*alpha
        ww2_2(:,:) = ww2_2(:,:)*alpha
        do k = 1+IBM_h,n(3)
          do j = 1,n(2)
            do i = 1,n(1)
              ip    = i+1   
              jp    = j+1   
              kp    = k+1  
              phix  = 0.5_rp*(phi(ip,j ,k )+phi(i,j,k))
              phiy  = 0.5_rp*(phi(i ,jp,k )+phi(i,j,k))
              phiz  = 0.5_rp*(phi(i ,j ,kp)+phi(i,j,k))
              u2    = 0.0_rp
              v2    = 0.0_rp
              w2    = 0.0_rp
              if(phix.LT.(1.0_rp-threshold)) u2 = u(i,j,k)+rex(i,j,k)/(1.0_rp-phix)
              if(phiy.LT.(1.0_rp-threshold)) v2 = v(i,j,k)+rey(i,j,k)/(1.0_rp-phiy)
              if(phiz.LT.(1.0_rp-threshold)) w2 = w(i,j,k)+rez(i,j,k)/(1.0_rp-phiz)
              uu2_2(j,k) = uu2_2(j,k) + beta*(1.0_rp/real(n(1),rp))*(u(i,j,k)*(1.0_rp-phix)+rex(i,j,k))*u2
              vv2_2(j,k) = vv2_2(j,k) + beta*(1.0_rp/real(n(1),rp))*(v(i,j,k)*(1.0_rp-phiy)+rey(i,j,k))*v2
              ww2_2(j,k) = ww2_2(j,k) + beta*(1.0_rp/real(n(1),rp))*(w(i,j,k)*(1.0_rp-phiz)+rez(i,j,k))*w2
            enddo
          enddo
        enddo
      endif
      !
      RETURN
    END SUBROUTINE collect_statistics

!=======================================================================
! Collect the statistics
!=======================================================================

    SUBROUTINE output_statistics(n,time,istep,datadir)
      !
      IMPLICIT NONE
      !
      integer , intent(in ), dimension(3)                       :: n
      real(rp), intent(in )                                     :: time
      integer , intent(in )                                     :: istep
      character(len=100),intent(in )                            :: datadir
      character(len=9)                                          :: fldnum
      !
      write(fldnum,'(i9.9)') istep
      !
      if(myid.EQ.0) print*,'Output time-averaged files for identifier ', fldnum, &
                         ' corresponding to an averaging timespan of: ', timespan
      !
      call   output23field(n,phi_2     ,'pfx1',datadir,time,istep)
      call   output23field(n,pc_2      ,'pcx1',datadir,time,istep)
      call   output23field(n,temp_2    ,'tmx1',datadir,time,istep)
      call   output23field(n,pr1       ,'p1x1',datadir,time,istep)
      call   output23field(n,pr2       ,'p2x1',datadir,time,istep)
      call   output23field(n,u1_2      ,'u1x1',datadir,time,istep)
      call   output23field(n,v1_2      ,'v1x1',datadir,time,istep)
      call   output23field(n,w1_2      ,'w1x1',datadir,time,istep)
      call   output23field(n,u2_2      ,'u2x1',datadir,time,istep)
      call   output23field(n,v2_2      ,'v2x1',datadir,time,istep)
      call   output23field(n,w2_2      ,'w2x1',datadir,time,istep)
      call   output23field(n,uu1_2     ,'uu11',datadir,time,istep)
      call   output23field(n,vv1_2     ,'vv11',datadir,time,istep)
      call   output23field(n,ww1_2     ,'ww11',datadir,time,istep)
      call   output23field(n,uu2_2     ,'uu21',datadir,time,istep)
      call   output23field(n,vv2_2     ,'vv21',datadir,time,istep)
      call   output23field(n,ww2_2     ,'ww21',datadir,time,istep)
      call outputwallfield(n,temp_avg  ,'twav',datadir,time,istep)
      call outputwallfield(n,temp_inst ,'twin',datadir,time,istep)
      call outputwallfield(n,phi_avg   ,'pwav',datadir,time,istep)
      call outputwallfield(n,phi_inst  ,'pwin',datadir,time,istep)
      call outputwallfield(n,tmptmp_avg,'ttav',datadir,time,istep)
      !
      RETURN
    END SUBROUTINE output_statistics
    !
  END MODULE mod_stats 
