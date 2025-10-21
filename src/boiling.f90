module mod_boiling
use mod_types
use mpi
use mod_common_mpi, only: myid,ierr,ijk_start,comm_cart,status,n_x,n_y,n_z
use iso_c_binding , only: C_PTR
use mod_param     , only: itot,jtot,ktot
use mod_sanity    , only: flutas_error
use mod_bound     , only: updthalo

implicit none

public

real(rp), dimension(:,:)  , allocatable         :: nuc_sites
integer , dimension(:)    , allocatable         :: nuc_sites_chk,nuc_sites_consider,site_activation
real(rp) , dimension(:)    , allocatable        :: nuc_sites_consider_dir
integer                                         :: num_sites,out_dir,in_dir,out_rk_safety
real(rp)                                        :: nuc_rad,out_len,const_out_time
real(rp)                                        :: u_gas,u_liquid,u_cap
logical                                         :: out_forced,in_active,use_mod_nse,use_mod_solver
logical                                         :: in_receive,in_sent,is_const_out_time
logical                                         :: do_phasechnage
real(rp)                                        :: in_force,in_len,max_kappa,pc_factor
real(rp)                                        :: pc_factor1,pc_factor2,pc_factor3
integer                                         :: comm_rank,k_max_nuc
real(rp), dimension(:,:,:), allocatable         :: in_vel,rex,rey,rez
real(rp), dimension(:,:)  , allocatable         :: psi_bound_o,ent_bound_o,mu_bound_o,heater_ind
real(rp), dimension(1:5)                        :: p_lim_array
real(rp)                                        :: bulk_vel,int_start_bulk
real(rp)                                        :: min_gamma, min_eps_prime
real(rp)                                        :: alpha_c, beta_c
character(len=100)                              :: nuc_file
character(len=4)                                :: sgstype
logical, dimension(0:1,3)                       :: wall
real(rp) , parameter                            :: c_smag  = 0.11_rp
real(rp) , parameter                            :: pr_turb = 0.9_rp
character(len=1),            dimension(0:1,3,3) :: cbcreg = reshape((/'D','D','D','D','D','D',&
                                                              'D','D','D','D','D','D',&
                                                              'D','D','D','D','D','D'/),&
                                                               shape(cbcreg))
character(len=1),            dimension(0:1,3) :: cbcpsi_n = reshape((/'N','N','N','N','N','N'/),&
                                                               shape(cbcpsi_n))  ! replacement for all quantities requiering phasefield like bc but do not allow for dirichlet bc's!
real(rp)   , dimension(0:1,3),parameter :: bcpsi_n = reshape((/0._rp,0._rp,0._rp,0._rp,0._rp,0._rp/),&
                                                             shape(bcpsi_n))
real(rp)   , dimension(0:1,3,3),parameter :: bcreg = reshape((/0._rp,0._rp,0._rp,0._rp,0._rp,0._rp,&
                                                             0._rp,0._rp,0._rp,0._rp,0._rp,0._rp,&
                                                             0._rp,0._rp,0._rp,0._rp,0._rp,0._rp/),&
                                                             shape(bcreg))
!
contains
! 
!_______________________________________________________________________________
!
!#################################################
!###     Start of:  IBM subrouties             ###
!#################################################
subroutine init_boiling(n,dl,trigger_out_3D)
  use mod_param, only: constant_dt,heater_switch, heater_x0,heater_x1,heater_y0,heater_y1,heater_power, &
                       cbcpsi
  !
  implicit none
  integer, intent(in   ), dimension(3)            :: n
  real(rp),intent(in   ), dimension(3)            :: dl
  logical, intent(inout)                          :: trigger_out_3D
  integer                                         :: iunit,i,j,q,idir
  real(rp)                                        :: x,y
  !
  trigger_out_3D = .TRUE.
  if(myid.eq.0) print*, 'Reading file: ', 'boiling.in' , ' ...'
  open(newunit=iunit,file='boiling.in',status='old',action='read',iostat=ierr)
  if( ierr.eq.0 ) then
    read(iunit,*) use_mod_nse                 !Use modified set of equations?
    read(iunit,*) use_mod_solver              !
    read(iunit,*) min_gamma, min_eps_prime    !
    read(iunit,*) do_phasechnage, pc_factor1,pc_factor2,pc_factor3  !
    read(iunit,*) is_const_out_time, const_out_time
    read(iunit,*) nuc_file                    ! file for bubble nucleation
    read(iunit,*) cbcreg(0,1,1),cbcreg(1,1,1),cbcreg(0,2,1),cbcreg(1,2,1),cbcreg(0,3,1),cbcreg(1,3,1)
    read(iunit,*) cbcreg(0,1,2),cbcreg(1,1,2),cbcreg(0,2,2),cbcreg(1,2,2),cbcreg(0,3,2),cbcreg(1,3,2)
    read(iunit,*) cbcreg(0,1,3),cbcreg(1,1,3),cbcreg(0,2,3),cbcreg(1,2,3),cbcreg(0,3,3),cbcreg(1,3,3)
    read(iunit,*) out_dir,out_len,out_forced
    read(iunit,*) in_dir,in_force,in_len,int_start_bulk,alpha_c,beta_c
    read(iunit,*) sgstype,wall(0,1),wall(1,1),wall(0,2),wall(1,2),wall(0,3),wall(1,3)
  else
    if(myid.eq.0) print*, 'Error reading the boiling.in input file' 
    if(myid.eq.0) print*, 'Input file missing or incomplete' 
    if(myid.eq.0) print*, 'Aborting...'
    call MPI_FINALIZE(ierr)
    call exit
  end if
  !
  close(iunit)
  !
  pc_factor = max(pc_factor1,pc_factor2,pc_factor3)
  !
  if(constant_dt.AND.is_const_out_time) then
    call flutas_error('ERROR: Constant timestep does not work in combination with constant outputtime!')
  endif
  !New location of the regularisation terms -> Why? -> reduce comutational cost by removing redundancy:
  allocate(rex(0:n(1)+1,0:n(2)+1,0:n(3)+1))
  allocate(rey(0:n(1)+1,0:n(2)+1,0:n(3)+1))
  allocate(rez(0:n(1)+1,0:n(2)+1,0:n(3)+1))
  rex(:,:,:) = 0.0
  rey(:,:,:) = 0.0
  rez(:,:,:) = 0.0
  !
  ! also set the maximum curvature to the grid spacing:
  max_kappa = 1.0/(min(dl(1),dl(2),dl(3))*1.0)
  k_max_nuc = n(3)
  call init_nucleation(n)
  call init_outflow(n)
  call init_inflow(n,dl)
  !
  !setup the IBM heater
  allocate(heater_ind(0:n(1)+1,0:n(2)+1))
  !
  heater_ind = 0.0_rp
  !
  if(heater_switch) then
    do j=0,n(2)+1
      y = (j+ijk_start(2)-0.5_rp)*dl(2)
      do i=0,n(1)+1
        x = (i+ijk_start(1)-0.5_rp)*dl(1)
        !
        if((x.GT.heater_x0).AND.(x.LT.heater_x1).AND. &
           (y.GT.heater_y0).AND.(y.LT.heater_y1)) then
          heater_ind(i,j) = 1.0_rp
        endif
        !
      enddo
    enddo
  endif
  !
  p_lim_array(:) = 0.0_rp
  ! correct cbcpsi_n for periodicity
  do q = 1,3
    do idir = 0,1
      if(cbcpsi(idir,q).EQ.'P') then
        cbcpsi_n(idir,q) = 'P'
      endif
    enddo
  enddo
  !
  return
end subroutine init_boiling
!_______________________________________________________________________________
!
subroutine set_IBM_vel_to_zero(n,nh_u,u,v,w)
  ! sets all velocities inside the solid to zero
  use mod_param,      only: IBM_h
  !   
  implicit none
  !   
  integer , intent(in   ), dimension(3)                       :: n
  integer , intent(in   )                                     :: nh_u
  real(rp), intent(inout), dimension(1-nh_u:,1-nh_u:,1-nh_u:) :: u,v,w
  integer                                                     :: i,j,k
  !
  !
  if(IBM_h.NE.0) then
    do k=1-nh_u,IBM_h
      do j=1-nh_u,n(2)+nh_u
        do i=1-nh_u,n(1)+nh_u
          u(i,j,k) = 0.0
          v(i,j,k) = 0.0
          w(i,j,k) = 0.0
        enddo
      enddo
    enddo
  endif
  !
end subroutine set_IBM_vel_to_zero
! 
!_______________________________________________________________________________
!
subroutine update_property_IBM(n,nh_p,prop12,psi,prop,IBM_cond)
  !
  use mod_param   , only: IBM_h
  !
  ! psi needs to be fully communicated with MPI, so that fluid properties don't need MPI communication
  implicit none
  !
  integer , intent(in ), dimension(3)                       :: n
  integer , intent(in )                                     :: nh_p
  real(rp), intent(in ), dimension(2)                       :: prop12
  real(rp), intent(in ), dimension(1-nh_p:,1-nh_p:,1-nh_p:) :: psi
  real(rp), intent(out), dimension(1-nh_p:,1-nh_p:,1-nh_p:) :: prop
  real(rp), intent(in )                                     :: IBM_cond
  !
  real(rp) :: prop12_1, prop12_2
  integer  :: i,j,k,ii
  integer  :: n1, n2, n3
  ! 
  n1 = n(1)
  n2 = n(2)
  n3 = n(3)
  !
  ii = 0
  !
  prop12_1 = prop12(1)
  prop12_2 = prop12(2)
  !
  !$acc parallel loop collapse(3)
  ! first for IBM part
  if(IBM_h.NE.0) then
    do k=1-nh_p,IBM_h
      do j=1-nh_p,n2+nh_p
        do i=1-nh_p,n1+nh_p
          prop(i,j,k) = IBM_cond
        enddo
      enddo
    enddo
  endif
  !
  if(IBM_h.EQ.0) ii = nh_p
  !
  do k=1+IBM_h-ii,n3+nh_p
    do j=1-nh_p,n2+nh_p
      do i=1-nh_p,n1+nh_p
        prop(i,j,k) = psi(i,j,k)*prop12_1+(1.0_rp-psi(i,j,k))*prop12_2
      enddo
    enddo
  enddo
  !$acc end parallel loop 
  !
  return
end subroutine update_property_IBM
!
!#################################################
!###     End of:  IBM subrouties               ###
!#################################################
!
!#################################################
!###     Start of:  Various initialisation     ###
!###                          subrouties       ###
!#################################################
!
subroutine init_ent(nh_p,h1,h2,h3,h03,n,psi,ind_IBM,ent,dl)
  !
  use mod_param     , only: rho1,cp1,tg0,h01,rho2,cp2,tl0,h02,IBM_rho, &
                            IBM_cp,IBM_iniT,inipsi,eps_int,pi,xc,yc,zc,r,IBM_h
  !
  implicit none
  integer  , intent(in   )                                     :: nh_p
  real(rp) , intent(inout)                                     :: h1,h2,h3,h03
  integer  , intent(in   ), dimension(3)                       :: n 
  real(rp) , intent(in   ), dimension(1-nh_p:,1-nh_p:,1-nh_p:) :: psi,ind_IBM
  real(rp) , intent(  out), dimension(1-nh_p:,1-nh_p:,1-nh_p:) :: ent
  real(rp) , intent(in   ), dimension(3)                       :: dl
  integer                                                      :: iunit,i,temp_dir,temp_count
  integer                                                      :: ix,loc,j,k
  character(len=100),parameter                                 :: file = 'temp_init.in'
  real(rp), dimension(:,:)  , allocatable                      :: init_temps
  real(rp)                                                     :: gamma,temp_loc,z,y,x,d
  !
  !
  h1 = rho1*(cp1*tg0+h01)
  h2 = rho2*(cp2*tl0+h02)
  h3 = h2   ! This is an arbitrary choice
  h03 = h3/IBM_rho - (IBM_cp*IBM_iniT)
  temp_count = 0
  temp_dir = 0
  if(myid.eq.0) print*, 'Reading file: ', trim(file) , ' ...'
  open(newunit=iunit,file=file,status='old',action='read',iostat=ierr)
  if( ierr.eq.0 ) then
    read(iunit,*) temp_dir
    read(iunit,*) temp_count
  elseif(inipsi.eq.'flm') then
    do k=1,n(3) !
      z = (k+ijk_start(3)-0.5_rp)*dl(3)
      do j=1,n(2)
        y = (j+ijk_start(2)-0.5_rp)*dl(2)
        do i=1,n(1)
          x = (i+ijk_start(1)-0.5_rp)*dl(1)
          !
          d = 0.07868440885009152  !Lambda_0 in Sun et al 2014
          !d = 0.0023003957362390023  !Lambda_0 in NTNU_test_4
          d = (d/128.0)*(4.0+cos(2.0*pi*(y-0.5*d)/d))           ! in Sun et al 2014
          !d = (d/1.0)*(0.025+0.0025*cos(2.0*pi*(y-0.5*d)/d))    ! in NTNU_test_4
          d = max((d-z)/d,0.0)
          d = (0.0 + d*1.0)       ! for the Sun et al 2014 test case
          !d = (0.0+d*5.0)       ! for the NTNU_test_4 test case
          ent(i,j,k) = (psi(i,j,k))*(d*(rho1*cp1) + h01*rho1) &
                 + (1.0-psi(i,j,k))*(d*(rho2*cp2) + h02*rho2)
          !
        enddo
      enddo
    enddo
  elseif(inipsi.eq.'sfn') then
    do k=1,n(3) !
      z = (k+ijk_start(3)-0.5_rp)*dl(3)
      do j=1,n(2)
        y = (j+ijk_start(2)-0.5_rp)*dl(2)
        do i=1,n(1)
          x = (i+ijk_start(1)-0.5_rp)*dl(1)
          !
          d = (z-zc(1))
          d = 0.5_rp*(1._rp+tanh( ((r(1))-d)/(2._rp*eps_int) )) !+2*dl(3)
          d = tl0 + (tg0-tl0)*d
          ent(i,j,k) = (psi(i,j,k))*(d*(rho1*cp1) + h01*rho1) &
          + (1.0-psi(i,j,k))*(d*(rho2*cp2) + h02*rho2)
          !
        enddo
      enddo
    enddo
  elseif(inipsi.eq.'one') then
    do k=1,n(3) !
      z = (k+ijk_start(3)-0.5_rp)*dl(3)
      do j=1,n(2)
        y = (j+ijk_start(2)-0.5_rp)*dl(2)
        do i=1,n(1)
          x = (i+ijk_start(1)-0.5_rp)*dl(1)
          !
          d = abs(z-zc(1))
          d = 0.5_rp*(1._rp+tanh( ((r(1))-d)/(2._rp*eps_int) )) !+2*dl(3)
          !d = tg0*d + (1.0-d)*tl0
          if((d.LT.0.5).and.(d.GE.0.183)) then
            d =  tg0 + 0.65*(tl0-tg0) ! was 0.65 -> with wall effect  (w/o antitrapping current 0.35) || 1.0 without wall effect
          elseif((d.LT.0.817).and.(d.GE.0.5)) then
            d =  tg0 + 0.0*(tl0-tg0) ! was 0.0 -> with wall effect || 1.0 without wall effect
          elseif (d.LT.0.183) then
            d = tl0
          else
            d = tl0
          endif
          ent(i,j,k) = (psi(i,j,k))*(d*(rho1*cp1) + h01*rho1) &
          + (1.0-psi(i,j,k))*(d*(rho2*cp2) + h02*rho2)
          !
        enddo
      enddo
    enddo
  elseif(inipsi.eq.'nuc') then
    h2 = rho2*(cp2*tl0+h02)
    h3 = h2   ! This is an arbitrary choice
    h03 = h3/IBM_rho - (IBM_cp*IBM_iniT)
    ent(:,:,1:IBM_h) = IBM_rho*(h03+IBM_cp*IBM_iniT)
    do k=1+IBM_h,n(3)
      do j=1,n(2)
        do i=1,n(1)
          d = psi(i,j,k)
          if((d.LT.0.499).and.(d.GE.0.183)) then
            d =  tg0 + 1.0*(tl0-tg0)
          elseif((d.LT.0.817).and.(d.GE.0.499)) then
            d =  tg0 + 0.5*(tl0-tg0)
          elseif (d.LT.0.183) then
            d = tl0
          else
            d = tg0
          endif
          ent(i,j,k) = (psi(i,j,k))*(d*(rho1*cp1) + h01*rho1) &
          + (1.0-psi(i,j,k))*(d*(rho2*cp2) + h02*rho2)
          !
        enddo
      enddo
    enddo
  elseif(inipsi.eq.'bub') then
    h2 = rho2*(cp2*tl0+h02)
    h3 = h2   ! This is an arbitrary choice
    h03 = h3/IBM_rho - (IBM_cp*IBM_iniT)
    ent(:,:,:) = IBM_rho*(h03+IBM_cp*IBM_iniT)
    do k=IBM_h+1,n(3) !
      z = (k+ijk_start(3)-0.5_rp)*dl(3)
      do j=1,n(2)
        y = (j+ijk_start(2)-0.5_rp)*dl(2)
        do i=1,n(1)
          x = (i+ijk_start(1)-0.5_rp)*dl(1)
          !
#if defined(_TWOD)
          d = sqrt(                  (y-yc(1))**2._rp+(z-zc(1))**2._rp )
#else
          d = sqrt( (x-xc(1))**2._rp+(y-yc(1))**2._rp+(z-zc(1))**2._rp )
#endif
          d = 0.5_rp*(1._rp+tanh( ((r(1))-d)/(2._rp*eps_int) )) !+2*dl(3)
          d = tl0 + (tg0-tl0)*d
          ent(i,j,k) = (psi(i,j,k))*(d*(rho1*cp1) + h01*rho1) &
          + (1.0-psi(i,j,k))*(d*(rho2*cp2) + h02*rho2)
          !
        enddo
      enddo
    enddo
  else
    h1 = rho1*(cp1*tg0+h01)
    h2 = rho2*(cp2*tl0+h02)
    h3 = h2   ! This is an arbitrary choice
    h03 = h3/IBM_rho - (IBM_cp*IBM_iniT)
    call update_property_IBM(n,nh_p,(/h1,h2/),psi,ent,h3)
  end if
  !
  if((myid.eq.0).and.(temp_dir.ne.0)) print*, 'Reading in: ', temp_count , ' Temperature points ...'
  if(temp_count.gt.0) then
    allocate(init_temps(1:temp_count,1:2))
    do i=1,temp_count
      !             loc             temperature       
      read(iunit,*) init_temps(i,1),init_temps(i,2)
    end do
    close(iunit)
  end if
  !
  !
  select case(temp_dir) ! temp_dir -> spatial direction of the profile
  case(123) ! spherical profile
    do k=1,n(3) !
      z = (k+ijk_start(3)-0.5_rp)*dl(3)
      do j=1,n(2)
        y = (j+ijk_start(2)-0.5_rp)*dl(2)
        do i=1,n(1)
          x = (i+ijk_start(1)-0.5_rp)*dl(1)
          d = sqrt( (x-nuc_sites(1,1))**2._rp &
                   +(y-nuc_sites(1,2))**2._rp &
                   +(z-nuc_sites(1,3))**2._rp )
          loc = 0
          do ix=1,temp_count
            if(init_temps(ix,1).LT.d) then
              loc = ix
            else
              exit
           endif
          enddo
          !
          gamma = (d-init_temps(loc,1))/(init_temps(loc+1,1)-init_temps(loc,1))
          temp_loc = (1.0-gamma)*init_temps(loc,2)+gamma*init_temps(loc+1,2)
          !
          if(k.GT.IBM_h) then
            ent(i,j,k) = (psi(i,j,k))*(temp_loc*(rho1*cp1) + h01*rho1) &
                   + (1.0-psi(i,j,k))*(temp_loc*(rho2*cp2) + h02*rho2)
          else
            ent(i,j,k) = (temp_loc*(IBM_rho*IBM_cp) + h03*IBM_rho)
          endif
          !
        enddo
      enddo
    enddo
  case(3)
    do k=1,n(3) !
      z = (k+ijk_start(3)-0.5_rp)*dl(3)
      loc = 0
      do ix=1,temp_count
        if(init_temps(ix,1).LT.z) then
          loc = ix
        else
          exit
        endif
      enddo
      !
      gamma = (z-init_temps(loc,1))/(init_temps(loc+1,1)-init_temps(loc,1))
      temp_loc = (1.0-gamma)*init_temps(loc,2)+gamma*init_temps(loc+1,2)
      !
      do j=1,n(2)
        do i=1,n(1)
          !
          if(k.GT.IBM_h) then
            ent(i,j,k) = (psi(i,j,k))*(temp_loc*(rho1*cp1) + h01*rho1) &
                   + (1.0-psi(i,j,k))*(temp_loc*(rho2*cp2) + h02*rho2)
          else
            ent(i,j,k) = (temp_loc*(IBM_rho*IBM_cp) + h03*IBM_rho)
          endif
          !
        enddo
      enddo
    enddo
  case(0)
    ! do nothing
  case default
    call flutas_error('ERROR: This temperature distribution direction is not implemented ... yet!')
  end select
  !
  if(temp_count.gt.0) deallocate(init_temps)
  !
end subroutine init_ent
!_______________________________________________________________________________

subroutine init_nucleation(n)
  !
  use mod_param,    only: dl,l,cbcpsi
  !
  implicit none
  integer  , intent(in   ), dimension(3)          :: n 
  integer                                         :: iunit,i,j,num_chk,num_consider,p
  real(rp)                                        :: dist_consider,xmin,xmax,ymin,ymax,zmin,zmax
  real(rp)                                        :: dist_from_rad
  !
  ! Load file named "nucleation_sites.in"
  if(myid.eq.0) print*, 'Reading file: ', trim(nuc_file) , ' ...'
  open(newunit=iunit,file=nuc_file,status='old',action='read',iostat=ierr)
  if( ierr.eq.0 ) then
    read(iunit,*) num_sites
    read(iunit,*) nuc_rad
  else
    if(myid.eq.0) print*, 'Error reading the ' , trim(nuc_file) , ' input file' 
    if(myid.eq.0) print*, 'Input file missing or incomplete' 
    if(myid.eq.0) print*, 'Aborting...'
    call MPI_FINALIZE(ierr)
    call exit
  end if
  !
  if(myid.eq.0) print*, 'Reading in: ', num_sites , ' nucleation sites ...'
  if(num_sites.gt.0) then
    allocate(nuc_sites(1:num_sites,1:4))
    allocate(site_activation(1:num_sites))
    site_activation = 0
    do i=1,num_sites
      !               x loc            y loc         z loc       activation temp
      read(iunit,*) nuc_sites(i,1),nuc_sites(i,2),nuc_sites(i,3),nuc_sites(i,4)
      !
    end do
  end if
  close(iunit)
  !
  ! Maximum index of k that needs to be considered for nucleation
  ! to safe computational time
  dist_from_rad = 9
  !
  k_max_nuc = ceiling((maxval(nuc_sites(:,3))+nuc_rad)/dl(3))+dist_from_rad
  !
  ! For paralell optimized nucleation
  ! allocate and array of nucleation sites
  ! (1) to check for activation
  ! (2) to loop over for seeding
  !================================================================================
  ! First: check how many nulceation sites are there for checking for this rank
  !
  num_chk = 0
  do p=1,num_sites
    ! First question: is the nuc_site in this processor?
    xmin = (     ijk_start(1))*dl(1)
    xmax = (n(1)+ijk_start(1))*dl(1)
    ymin = (     ijk_start(2))*dl(2)
    ymax = (n(2)+ijk_start(2))*dl(2)
    zmin = (     ijk_start(3))*dl(3)
    zmax = (n(3)+ijk_start(3))*dl(3)
    if((nuc_sites(p,1).GE.xmin).AND. &
       (nuc_sites(p,1).LT.xmax).AND. &
       (nuc_sites(p,2).GE.ymin).AND. &
       (nuc_sites(p,2).LT.ymax).AND. &
       (nuc_sites(p,3).GE.zmin).AND. &
       (nuc_sites(p,3).LT.zmax)) then
      num_chk = num_chk + 1
    endif
  enddo
  ! Allocate the check array
  if(num_chk.GT.0) then
    allocate(nuc_sites_chk(1:num_chk))
    ! Fill the check array
    num_chk = 0
    do p=1,num_sites
      ! First question: is the nuc_site in this processor?
      xmin = (     ijk_start(1))*dl(1)
      xmax = (n(1)+ijk_start(1))*dl(1)
      ymin = (     ijk_start(2))*dl(2)
      ymax = (n(2)+ijk_start(2))*dl(2)
      zmin = (     ijk_start(3))*dl(3)
      zmax = (n(3)+ijk_start(3))*dl(3)
      if((nuc_sites(p,1).GE.xmin).AND. &
         (nuc_sites(p,1).LT.xmax).AND. &
         (nuc_sites(p,2).GE.ymin).AND. &
         (nuc_sites(p,2).LT.ymax).AND. &
         (nuc_sites(p,3).GE.zmin).AND. &
         (nuc_sites(p,3).LT.zmax)) then
        num_chk = num_chk + 1
        nuc_sites_chk(num_chk) = p
      endif
    enddo
  endif
  !================================================================================
  ! Secondly, check which nucleation sites are this rank has to consider for seeding
  dist_consider = nuc_rad + dist_from_rad*min(dl(3),dl(2),dl(1))
  !
  num_consider = 0
  do p=1,num_sites
    ! First question: is the nuc_site in this processor?
    xmin = (     ijk_start(1))*dl(1)-dist_consider
    xmax = (n(1)+ijk_start(1))*dl(1)+dist_consider
    ymin = (     ijk_start(2))*dl(2)-dist_consider
    ymax = (n(2)+ijk_start(2))*dl(2)+dist_consider
    zmin = (     ijk_start(3))*dl(3)-dist_consider
    zmax = (n(3)+ijk_start(3))*dl(3)+dist_consider
    if((nuc_sites(p,1).GE.xmin).AND. &
       (nuc_sites(p,1).LT.xmax).AND. &
       (nuc_sites(p,2).GE.ymin).AND. &
       (nuc_sites(p,2).LT.ymax).AND. &
       (nuc_sites(p,3).GE.zmin).AND. &
       (nuc_sites(p,3).LT.zmax)) then
      num_consider = num_consider + 1
    elseif(((nuc_sites(p,1)+l(1)).GE.xmin).AND. &
           ((nuc_sites(p,1)+l(1)).LT.xmax).AND. &     !periodicity in the x+dirction
            (nuc_sites(p,2).GE.ymin).AND. &
            (nuc_sites(p,2).LT.ymax).AND. &
            (nuc_sites(p,3).GE.zmin).AND. &
            (nuc_sites(p,3).LT.zmax).AND. &
            (cbcpsi(0,1 )  .EQ.'P' ).AND. &
            (cbcpsi(1,1 )  .EQ.'P' )) then
      num_consider = num_consider + 1
    elseif(((nuc_sites(p,1)-l(1)).GE.xmin).AND. &    !periodicity in the x-dirction
           ((nuc_sites(p,1)-l(1)).LT.xmax).AND. &    
            (nuc_sites(p,2).GE.ymin).AND. &
            (nuc_sites(p,2).LT.ymax).AND. &
            (nuc_sites(p,3).GE.zmin).AND. &
            (nuc_sites(p,3).LT.zmax).AND. &
            (cbcpsi(0,1 )  .EQ.'P' ).AND. &
            (cbcpsi(1,1 )  .EQ.'P' )) then
      num_consider = num_consider + 1
    endif
  enddo
  ! Allocate the check array
  if(num_consider.GT.0) then
    allocate(nuc_sites_consider(1:num_consider))
    allocate(nuc_sites_consider_dir(1:num_consider))   ! for periodicity
    nuc_sites_consider_dir(:) = 0.0_rp                      ! 0 -> no accounting for periodicity
    ! Fill the check array
    num_consider = 0
    do p=1,num_sites
      ! First question: is the nuc_site in this processor?
    xmin = (     ijk_start(1))*dl(1)-dist_consider
    xmax = (n(1)+ijk_start(1))*dl(1)+dist_consider
    ymin = (     ijk_start(2))*dl(2)-dist_consider
    ymax = (n(2)+ijk_start(2))*dl(2)+dist_consider
    zmin = (     ijk_start(3))*dl(3)-dist_consider
    zmax = (n(3)+ijk_start(3))*dl(3)+dist_consider
      if((nuc_sites(p,1).GE.xmin).AND. &
         (nuc_sites(p,1).LT.xmax).AND. &
         (nuc_sites(p,2).GE.ymin).AND. &
         (nuc_sites(p,2).LT.ymax).AND. &
         (nuc_sites(p,3).GE.zmin).AND. &
         (nuc_sites(p,3).LT.zmax)) then
        num_consider = num_consider + 1
        nuc_sites_consider(num_consider) = p
        nuc_sites_consider_dir(num_consider) = 0.0_rp 
      elseif(((nuc_sites(p,1)+l(1)).GE.xmin).AND. &
             ((nuc_sites(p,1)+l(1)).LT.xmax).AND. &     !periodicity in the x+dirction
              (nuc_sites(p,2).GE.ymin).AND. &
              (nuc_sites(p,2).LT.ymax).AND. &
              (nuc_sites(p,3).GE.zmin).AND. &
              (nuc_sites(p,3).LT.zmax).AND. &
              (cbcpsi(0,1 )  .EQ.'P' ).AND. &
              (cbcpsi(1,1 )  .EQ.'P' )) then
        num_consider = num_consider + 1
        nuc_sites_consider(num_consider) = p
        nuc_sites_consider_dir(num_consider) = 1.0_rp  
      elseif(((nuc_sites(p,1)-l(1)).GE.xmin).AND. &    !periodicity in the x-dirction
             ((nuc_sites(p,1)-l(1)).LT.xmax).AND. &  
              (nuc_sites(p,2).GE.ymin).AND. &
              (nuc_sites(p,2).LT.ymax).AND. &
              (nuc_sites(p,3).GE.zmin).AND. &
              (nuc_sites(p,3).LT.zmax).AND. &
              (cbcpsi(0,1 )  .EQ.'P' ).AND. &
              (cbcpsi(1,1 )  .EQ.'P' )) then
        num_consider = num_consider + 1
        nuc_sites_consider(num_consider) = p
        nuc_sites_consider_dir(num_consider) = -1.0_rp 
      endif
    enddo
  else
    allocate(nuc_sites_consider(1:1))
    allocate(nuc_sites_consider_dir(1:1))
    nuc_sites_consider(1) = 1           ! a generic site
    nuc_sites_consider_dir(1) =  0.0_rp
  endif
  !
end subroutine init_nucleation
!_______________________________________________________________________________
!
!
!#################################################
!###     End of:  Various initialisation       ###
!###                          subrouties       ###
!#################################################
!
!#################################################
!###     Start of:  Outflow subrouties         ###
!###     Implementation based on:              ###
!### Akash Dhruv, 2024, "A Vortex Damping      ###
!### Outflow Forcing for Multiphase Flows with ###
!### Sharp Interfacial Jumps",2306.10174,arXiv ###
!### https://doi.org/10.48550/arXiv.2306.10174 ###
!#################################################


subroutine init_outflow(n)
  implicit none
  integer , intent(in   ), dimension(3)     :: n 
  !
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
  out_rk_safety = 5
  !
end subroutine init_outflow

!_______________________________________________________________________________

subroutine cmpt_avg_out_vel(n,dl,l,nh_p,nh_u,psi,u,v,w)
  use mod_param     , only: sigma,rho1,rho2,gacc_x,gacc_y,gacc_z,IBM_h
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


  cache_g1 = 0.0
  cache_g2 = 0.0
  cache_l1 = 0.0
  cache_l2 = 0.0
  !
  select case(out_dir)
  !###################### Positive Y-direction ###################################
  case(2)
    ! x,y,z are face centered locations for y velocity!!!
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
    call mpi_allreduce(MPI_IN_PLACE,cache_g1 ,1,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
    call mpi_allreduce(MPI_IN_PLACE,cache_g2 ,1,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
    call mpi_allreduce(MPI_IN_PLACE,cache_l1 ,1,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
    call mpi_allreduce(MPI_IN_PLACE,cache_l2 ,1,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
    u_gas    = cache_g1/cache_g2
    if (cache_g2.EQ.0.0) u_gas = 0.0
    u_liquid = cache_l1/cache_l2
    u_liquid = max(u_cap,u_liquid)
    u_gas = max(u_cap,u_gas)
  !###################### Positive Z-direction ###################################
  case(3)
    ! x,y,z are face centered locations for z velocity!!!
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
    call mpi_allreduce(MPI_IN_PLACE,cache_g1 ,1,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
    call mpi_allreduce(MPI_IN_PLACE,cache_g2 ,1,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
    call mpi_allreduce(MPI_IN_PLACE,cache_l1 ,1,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
    call mpi_allreduce(MPI_IN_PLACE,cache_l2 ,1,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
    u_gas    = cache_g1/cache_g2
    if (cache_g2.EQ.0.0) u_gas = 0.0
    u_liquid = cache_l1/cache_l2
    u_liquid = max(u_cap,u_liquid)
    u_gas = max(u_cap,u_gas)
  !###################### Negative Z-direction ###################################
  case(-3)
    ! x,y,z are face centered locations for z velocity!!!
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
    call mpi_allreduce(MPI_IN_PLACE,cache_g1 ,1,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
    call mpi_allreduce(MPI_IN_PLACE,cache_g2 ,1,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
    call mpi_allreduce(MPI_IN_PLACE,cache_l1 ,1,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
    call mpi_allreduce(MPI_IN_PLACE,cache_l2 ,1,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
    u_gas    = cache_g1/cache_g2
    if (cache_g2.EQ.0.0) u_gas = 0.0
    u_liquid = cache_l1/cache_l2
    u_liquid = max(u_cap,u_liquid)
    u_gas = max(u_cap,u_gas)
  case(0)
    ! Do nothing
  case default
    call flutas_error('ERROR: This outflow direction is not implemented ... yet!')
  end select

end subroutine cmpt_avg_out_vel

!_______________________________________________________________________________

subroutine corr_u_liquid(n,dl,l,nh_p,nh_u,psi,f_t12,src_ph,u,v,w)
  ! Assumption: all liquid flows out through the outflow
  use mod_param     , only: sigma,rho1,rho2,gacc_x,gacc_y,gacc_z,IBM_h,bcvel,is_outflow
  use mod_common_mpi, only: MPI_PROC_NULL,back
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
  !
  !if(.NOT.(out_dir.EQ.0)) then
  !  do k=IBM_h+1,n(3) !
  !    do j=1,n(2)
   !     do i=1,n(1)
   !       vol_change = vol_change + (1.d0/rho1 - 1.d0/rho2) &
   !                                *dl(3)*dl(1)*dl(2)*src_ph(i,j,k)
   !     enddo
   !   enddo
   ! enddo
    
   ! call mpi_allreduce(MPI_IN_PLACE,vol_change,1,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
 ! endif

  !
  select case(out_dir)
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
  if(.NOT.(out_dir.EQ.0)) then
    if(.NOT.out_forced) then
      call mpi_allreduce(MPI_IN_PLACE,integral1 ,1,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      call mpi_allreduce(MPI_IN_PLACE,integral2 ,1,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      !               Safety margin    ?
      u_liquid_corr = 1.02*abs((compute_influxes(n,nh_u,u,v,w,f_t12,dl) + vol_change - integral1)/integral2) !
      u_liquid_corr = max(u_cap,u_liquid_corr,u_liquid)!
      !
      if(myid.eq.0) print*, 'Change liquid outflow vel by:',   100*(u_liquid_corr-u_liquid)/u_liquid , " %"
    else
      call mpi_allreduce(MPI_IN_PLACE,integral3 ,1,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      u_liquid_corr = 1.4*abs((compute_influxes(n,nh_u,u,v,w,f_t12,dl) + vol_change)/integral3)  !1.2*
      u_liquid_corr = max(u_cap,u_liquid_corr,u_liquid)

      if(myid.eq.0) print*, 'Change liquid outflow vel by:',   100*(u_liquid_corr-u_liquid)/u_liquid , " %"
    endif
    ! Now update
    u_liquid = u_liquid_corr
    u_gas = u_liquid_corr
    !
  endif
  !
end subroutine corr_u_liquid

!_______________________________________________________________________________

subroutine add_outflow_force(n,dl,l,nh_p,nh_u,psi,rho,u,v,w,uold,vold,wold,f_t12)
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
  integer                                                             :: i,j,k,upwind
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
          !h_val    = max((l(2)-y-0.5*dl(2)),0.0_rp) ! is is bad!
          h_val    = max((l(2)-y+0.5*dl(2)),0.0_rp)
          h_val    = (2.0/(1.0+exp(4.0*h_val/out_len)))          
          u_e      = phiz*u_gas + (1.0-phiz)*u_liquid
          ! scaling part:
          !dwdt = h_val*(min(1.0,abs(u_e/wold(i,j,k)))-1.0)*wold(i,j,k)/f_t12
          dwdt = - 0.2*h_val*wold(i,j,k)/f_t12
          ! convection 
          dwdt = dwdt - abs(u_e)*h_val*(wold(i,j,k)-wold(i,j-1,k))/dl(2)   ! first order upwind scheme was important for stability!
          !_______________________________________
          ! Scaling in y
          !h_val    = max((l(2)-y-1.0_rp*dl(2)),0.0_rp) ! is is bad!
          h_val    = max((l(2)-y),0.0_rp)
          h_val    = (2.0/(1.0+exp(4.0*h_val/out_len)))
          u_e      = phiy*u_gas + (1.0-phiy)*u_liquid
          ! scaling part:
          dvdt = h_val*(min(1.0,abs(u_e/vold(i,j,k)))-1.0)*vold(i,j,k)/f_t12
          if(vold(i,j,k).le.(0.1*u_e)) dvdt = - h_val*(vold(i,j,k)-(0.1*u_e))/f_t12
          ! convection 
          dvdt = dvdt - abs(u_e)*h_val*(vold(i,j,k)-vold(i,j-1,k))/dl(2) ! first order upwind scheme was important for stability!
          !_______________________________________
          ! Scaling in x
          !h_val    = max((l(2)-y-0.5*dl(2)),0.0_rp) ! is is bad!
          h_val    = max((l(2)-y+0.5*dl(2)),0.0_rp)
          h_val    = (2.0/(1.0+exp(4.0*h_val/out_len))) 
#if defined(_TWOD)
          u_e      = 0.0
#else
          u_e      = phix*u_gas + (1.0-phix)*u_liquid
#endif
          ! scaling part:
          !dudt = h_val*(min(1.0,abs(u_e/uold(i,j,k)))-1.0)*uold(i,j,k)/f_t12
          dudt = - 0.2*h_val*uold(i,j,k)/f_t12
          ! convection 
          dudt = dudt - abs(u_e)*h_val*(uold(i,j,k)-uold(i,j-1,k))/dl(2) ! first order upwind scheme was important for stability!
          ! Even in the Two-D Case, the dudt term is important for stability
          !_______________________________________
          u(i,j,k) = u(i,j,k) + rhox*dudt*f_t12
          v(i,j,k) = v(i,j,k) + rhoy*dvdt*f_t12
          w(i,j,k) = w(i,j,k) + rhoz*dwdt*f_t12
        enddo
      enddo
    enddo
  !###################### Positive Z-direction ###################################
  case(3)
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
          dwdt = h_val*(min(1.0,abs(u_e/wold(i,j,k)))-1.0)*wold(i,j,k)/f_t12
          ! convection 
          dwdt = dwdt - abs(u_e)*h_val*(wold(i,j,k)-wold(i,j,k-1))/dl(3)   ! first order upwind scheme was important for stability!
          !_______________________________________
          ! Scaling in y
          h_val    = (2.0/(1.0+exp(4.0*(l(3)-z+0.5*dl(3))/out_len)))          
          u_e      = phiy*u_gas + (1.0-phiy)*u_liquid
          ! scaling part:
          dvdt = h_val*(min(1.0,abs(u_e/vold(i,j,k)))-1.0)*vold(i,j,k)/f_t12
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
          dudt = h_val*(min(1.0,abs(u_e/uold(i,j,k)))-1.0)*uold(i,j,k)/f_t12
          ! convection 
          dudt = dudt - abs(u_e)*h_val*(uold(i,j,k)-uold(i,j,k-1))/dl(3) ! first order upwind scheme was important for stability!
          ! Even in the Two-D Case, the dudt term is important for stability
          !_______________________________________
          u(i,j,k) = u(i,j,k) + rhox*dudt*f_t12
          v(i,j,k) = v(i,j,k) + rhoy*dvdt*f_t12
          w(i,j,k) = w(i,j,k) + rhoz*dwdt*f_t12
        enddo
      enddo
    enddo
  !###################### Negative Z-direction ###################################
  case(-3)
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
          dwdt = rhoz*h_val*(min(1.0,abs(u_e/wold(i,j,k)))-1.0)*wold(i,j,k)/f_t12
          ! convection 
          upwind = 1
          !if(wold(i,j,k).LT.0.0) upwind = 1
          dwdt = dwdt + rhoz*abs(u_e)*h_val*(wold(i,j,k+1)-wold(i,j,k))/dl(3)   ! first order upwind scheme was important for stability!!
          !_______________________________________
          ! Scaling in y
          h_val    = (2.0/(1.0+exp(4.0*(z-0.5*dl(3))/out_len)))          
          u_e      = phiy*u_gas + (1.0-phiy)*u_liquid
          ! scaling part:
          dvdt = rhoy*h_val*(min(1.0,abs(u_e/vold(i,j,k)))-1.0)*vold(i,j,k)/f_t12
          ! convection 
          !upwind = 0
          !if((0.25*(wold(i,j,k)+wold(i,j,k-1)+wold(i,j+1,k)+wold(i,j+1,k-1))).LT.0.0) upwind = 1
          dvdt = dvdt + rhoy*abs(u_e)*h_val*(vold(i,j,k+1)-vold(i,j,k))/dl(3)
          !_______________________________________
          ! Scaling in x
          h_val    = (2.0/(1.0+exp(4.0*(z-0.5*dl(3))/out_len)))          
#if defined(_TWOD)
          u_e      = 0.0
#else
          u_e      = phix*u_gas + (1.0-phix)*u_liquid
#endif
          ! scaling part:
          dudt = rhox*h_val*(min(1.0,abs(u_e/uold(i,j,k)))-1.0)*uold(i,j,k)/f_t12
          ! convection 
          !upwind = 0
          !if((0.25*(wold(i,j,k)+wold(i,j,k-1)+wold(i+1,j,k)+wold(i+1,j,k-1))).LT.0.0) upwind = 1
          dudt = dudt + rhox*abs(u_e)*h_val*(uold(i,j,k+1)-uold(i,j,k))/dl(3)
          ! Sign is switched here because of negative direction!!!
          ! Even in the Two-D Case, the dudt term is important for stability
          !_______________________________________
          u(i,j,k) = u(i,j,k) + dudt*f_t12
          v(i,j,k) = v(i,j,k) + dvdt*f_t12
          w(i,j,k) = w(i,j,k) + dwdt*f_t12
        enddo
      enddo
    enddo
  case(0)
    ! Do nothing
  case default
    call flutas_error('ERROR: This outflow direction is not implemented ... yet!')
  end select
end subroutine add_outflow_force

!_______________________________________________________________________________

subroutine outflow_psi(n,dl,nh_v,phiold,u,v,w,dt,bound_o,tmp,ent_switch,back_val)
  !
  use mod_common_mpi, only: MPI_PROC_NULL,left,right,front,back,top,bottom
  use mod_param,      only: h02,rho2
  !                             
  ! Calculates the convected boundary condition of the input field "phiold"
  !
  implicit none
  !
  integer , intent(in ), dimension(3)                       :: n
  logical , intent(in )                                     :: ent_switch
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
          ! Assumption du/dn = 0,dw/dn = 0
          u_loc             = 0.5*(u(i ,j ,k )+u(im,j ,k ))
          w_loc             = 0.5*(w(i ,j ,k )+w(i ,j ,km))
          v_loc             =     (v(i ,j ,k )            )
          !
          if((v_loc.LT.0.0)) then
            aux             = back_val!  0.5*(phiold(i,j,k)+back_val)!
          else
            aux             = phiold(i,j,k)!phiold(i ,jp,k )
          endif
          !
          psi_wall          = 0.5*((phiold(i ,jp,k )+phiold(i ,j ,k ))            )
          !
          phi_grady         =     ((aux             -phiold(i ,j ,k ))/dl(3)      )
          phi_gradz         = 0.5*((phiold(i ,j ,kp)-phiold(i ,j ,km))/dl(2)/(2.0) &
                                  +(phiold(i ,jp,kp)-phiold(i ,jp,km))/dl(2)/(2.0))
          phi_gradx         = 0.5*((phiold(ip,j ,k )-phiold(im,j ,k ))/dl(1)/(2.0) &
                                  +(phiold(ip,jp,k )-phiold(im,jp,k ))/dl(1)/(2.0))

          !
          !if(v_loc.LT.0.0) phi_grady = 0.0
          !
          psi_wall          = max(0.0,psi_wall - dt*(u_loc*phi_gradx &
                                                    +v_loc*phi_grady &
                                                    +w_loc*phi_gradz))
          !                                  
          bound_o(i,k)    =  aux!psi_wall ! because of Dirichlet BC for cell centered scheme
          !
        enddo
      enddo
    endif
  case(3)
    !###################### Positive Z-direction ###################################
    ! if the outflow is in positiv z direction, all ranks are outflow ranks...
    !
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
        ! Assumption du/dn = 0,dv/dn = 0
        u_loc             = 0.5*(u(i ,j ,k )+u(im,j ,k ))
        v_loc             = 0.5*(v(i ,j ,k )+v(i ,jm,k ))
        w_loc             =     (w(i ,j ,k )            )   
        !
        if((w_loc.LT.0.0)) then
          aux             = 0.5*(phiold(i,j,k)+back_val)!back_val
        else
          aux             = phiold(i,j,k)!phiold(i ,j ,kp)
        endif
        !
        psi_wall          = 0.5*((phiold(i ,j ,k )+phiold(i ,j ,kp))            )
        !
        phi_gradz         =     ((aux             -phiold(i ,j ,k ))/dl(3)      )
        phi_grady         = 0.5*((phiold(i ,jp,k )-phiold(i ,jm,k ))/dl(2)/(2.0) &
                                +(phiold(i ,jp,kp)-phiold(i ,jm,kp))/dl(2)/(2.0))
        phi_gradx         = 0.5*((phiold(ip,j ,k )-phiold(im,j ,k ))/dl(1)/(2.0) &
                                +(phiold(ip,j ,kp)-phiold(im,j ,kp))/dl(1)/(2.0))
        !
        !if(w_loc.LT.0.0) phi_gradz = 0.0
        !
        psi_wall          = psi_wall - dt*(u_loc*phi_gradx &
                                                  +v_loc*phi_grady &
                                                  +w_loc*phi_gradz)
        !                                  
        bound_o(i,j)    =  aux !psi_wall ! because of Dirichlet BC for cell centered scheme
        !
      enddo
    enddo
  case(-3)
    !###################### Negative Z-direction ###################################
    ! if the outflow is in negative z direction, all ranks are outflow ranks...
    !
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
        ! Assumption du/dn = 0,dv/dn = 0
        u_loc             = 0.5*(u(i ,j ,k )+u(im,j ,k ))
        v_loc             = 0.5*(v(i ,j ,k )+v(i ,jm,k ))
        w_loc             =     (w(i ,j ,k-1)           )   
        !
        if((w_loc.GT.0.0)) then
          aux             = 0.5*(phiold(i,j,k)+back_val)!back_val
        else
          aux             = phiold(i,j,k)!phiold(i ,j ,km)
        endif
        !
        psi_wall          = 0.5*((phiold(i ,j ,k )+phiold(i ,j ,km))            )
        !
        phi_gradz         =     ((phiold(i ,j ,k )-aux             )/dl(3)      )
        phi_grady         = 0.5*((phiold(i ,jp,k )-phiold(i ,jm,k ))/dl(2)/(2.0) &
                                +(phiold(i ,jp,km)-phiold(i ,jm,km))/dl(2)/(2.0))
        phi_gradx         = 0.5*((phiold(ip,j ,k )-phiold(im,j ,k ))/dl(1)/(2.0) &
                                +(phiold(ip,j ,km)-phiold(im,j ,km))/dl(1)/(2.0))
        !
        !if(w_loc.LT.0.0) phi_gradz = 0.0
        !
        psi_wall          = psi_wall - dt*(u_loc*phi_gradx &
                                                  +v_loc*phi_grady &
                                                  +w_loc*phi_gradz)
        !                                  
        bound_o(i,j)    =  aux!psi_wall ! because of Dirichlet BC for cell centered scheme
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

!_______________________________________________________________________________

subroutine boundsb_arrayo(n,nh_p,halo,s,bound_o)
  !
  use mod_common_mpi, only: MPI_PROC_NULL,left,right,front,back,top,bottom
  !   
  ! imposes Dirichlet BC at the outflow that was calculated from outflow_psi
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


!#################################################
!###     End of:  Outflow subrouties           ###
!#################################################
!
!
!
!
!#################################################
!###     Start of:  Inflow subrouties          ###
!###     Problem: Turbulent inflow required    ###
!###     Flow: single phase inflow, two phase- ###
!###      outflow -> no periodicity possible   ###
!###     Idea: recycle velocities              ###
!#################################################


subroutine init_inflow(n,dl)
  use mod_param     , only: dims_in,itot,jtot,ktot
  implicit none
  integer , intent(in   ), dimension(3)        :: n 
  real(rp), intent(in   ), dimension(3)        :: dl
  integer                                      :: iunit,x_rank,x_start
  integer                                      :: sent,receive
  character(len=100),parameter                 :: file = 'inflow.in'
  real(rp)                                     :: xmin,xmax,ymin,ymax,zmin,zmax
  !
  in_active  = .FALSE.
  in_receive = .FALSE.
  comm_rank  = 0
  ! Load file named "inflow.in"
  !
  ! Following variables define the box location of this rank in global space coordinates
  xmin = (     ijk_start(1))*dl(1)
  xmax = (n(1)+ijk_start(1))*dl(1)
  ymin = (     ijk_start(2))*dl(2)
  ymax = (n(2)+ijk_start(2))*dl(2)
  zmin = (     ijk_start(3))*dl(3)
  zmax = (n(3)+ijk_start(3))*dl(3)
  !
  select case(in_dir)
  case(2)
    if(ijk_start(2).EQ.0) then
      in_active = .TRUE.
      in_receive = .TRUE.
    elseif((ymin.LT.in_len).AND.(ymax.GE.in_len)) then
      in_active = .TRUE.
      in_sent = .TRUE.
    endif
    if(in_active) then
      allocate(in_vel(1:n(1),1:n(3),1:8))
      in_vel(:,:,:) = 0.0
    endif
    do x_rank = 1,dims_in(1)
      x_start = (x_rank-1)*(itot/dims_in(1))
      sent    = 0
      receive = 0
      !
      if(in_active.AND.in_sent.AND.(x_start.EQ.ijk_start(1)))then
        sent = myid
      elseif(in_active.AND.in_receive.AND.(x_start.EQ.ijk_start(1)))then
        receive = myid
      endif
      !
      call mpi_allreduce(MPI_IN_PLACE,sent   ,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
      call mpi_allreduce(MPI_IN_PLACE,receive,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
      !
      if(in_active.AND.in_sent.AND.(x_start.EQ.ijk_start(1)))then
        comm_rank = receive
      elseif(in_active.AND.in_receive.AND.(x_start.EQ.ijk_start(1)))then
        comm_rank = sent
      endif
    enddo
  case(0)
    ! Do nothing
  case default
    call flutas_error('ERROR: This infow direction is not implemented ... yet!')
  end select
end subroutine init_inflow

!_______________________________________________________________________________

subroutine collect_vel(n,dl,l,nh_p,nh_u,psi,u,v,w,p)
  use mod_param     , only: itot,jtot,ktot,IBM_h
  implicit none
  integer         , intent(in ), dimension(3)                       :: n
  real(rp)        , intent(in ), dimension(3)                       :: dl,l
  integer         , intent(in )                                     :: nh_p,nh_u
  real(rp)        , intent(in ), dimension(1-nh_p:,1-nh_p:,1-nh_p:) :: psi,p
  real(rp)        , intent(in ), dimension(1-nh_u:,1-nh_u:,1-nh_u:) :: u,v,w
  integer                                                           :: i,j,k
  integer                                                           :: lb1,lb2,lb3
  real(rp)                                                          :: bulk
  !
  !
  select case(in_dir)
  case(2)
    ! Check if rank is in region of the recycling plane
    if(in_active) then
      ! First set velocities to zero:
      in_vel(:,:,:) = 0.0
    endif
    !
    if(in_active.AND.in_sent) then
      !
      do k=1,n(3) !
        j = nint(in_len/dl(2) - ijk_start(2))
        do i=1,n(1)
          in_vel(i,k,1) = u(i,j  ,k)
          in_vel(i,k,2) = v(i,j  ,k)
          in_vel(i,k,3) = w(i,j  ,k)
          in_vel(i,k,4) = p(i,j  ,k)
          in_vel(i,k,5) = u(i,j+1,k)
          in_vel(i,k,6) = v(i,j+1,k)
          in_vel(i,k,7) = w(i,j+1,k)
          in_vel(i,k,8) = p(i,j+1,k)
        enddo
      enddo
    endif
    ! Sent it!
    if(in_active.AND.in_sent) then
      lb1 = lbound(in_vel,1) 
      lb2 = lbound(in_vel,2)
      lb3 = lbound(in_vel,3)
      call MPI_SEND(in_vel(lb1,lb2,lb3), size( in_vel ),MPI_REAL_RP,comm_rank ,0, &
                        comm_cart,ierr)
    endif
    ! Receive it!
    if(in_active.AND.in_receive) then
      lb1 = lbound(in_vel,1) 
      lb2 = lbound(in_vel,2)
      lb3 = lbound(in_vel,3)
      call MPI_RECV(in_vel(lb1,lb2,lb3), size( in_vel ),MPI_REAL_RP,comm_rank ,0, &
                        comm_cart,status,ierr)
    endif
    bulk = 0.0
    if(in_receive) then
      bulk = sum(in_vel(1:n(1),1+IBM_h:n(3),2)) + sum(in_vel(1:n(1),1+IBM_h:n(3),6))
    endif
    call mpi_allreduce(MPI_IN_PLACE,bulk ,1,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
    bulk_vel = bulk/(itot*(ktot-IBM_h)*2)
  case(0)
    ! Do nothing
  case default
    call flutas_error('ERROR: This inflow direction is not implemented ... yet!')
  end select
  !
end subroutine collect_vel

!_______________________________________________________________________________

subroutine apply_bound_vel(n,dl,l,nh_p,nh_u,u,v,w,halo_u,halo_p,p,mom_switch,rho)
  implicit none
  integer       , intent(in   ), dimension(3)                       :: n
  real(rp)      , intent(in   ), dimension(3)                       :: dl,l
  integer       , intent(in   )                                     :: nh_p,nh_u
  real(rp)      , intent(in   ), dimension(1-nh_p:,1-nh_p:,1-nh_p:) :: rho
  real(rp)      , intent(inout), dimension(1-nh_p:,1-nh_p:,1-nh_p:) :: p
  real(rp)      , intent(inout), dimension(1-nh_u:,1-nh_u:,1-nh_u:) :: u,v,w
  integer                                                           :: i,j,k,q
  real(rp)                                                          :: xmin,xmax,ymin,ymax,zmin,zmax
  real(rp)                                                          :: rhox,rhoy,rhoz
  integer                                                           :: ind1,ind2
  integer       , intent(in   ), dimension(3)                       :: halo_u,halo_p
  logical       , intent(in   )                                     :: mom_switch
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
  !
  select case(in_dir)
  case(2)
    ! Check if rank is region of the recycling plane on the inflow plane
    if((ijk_start(2).EQ.0)) then
      if(mom_switch) then
        do k=1,n(3)
          do i=1,n(1)
            ! rho is assumed to have neuman BC's
            rhox = 0.5d0*(rho(i,1,k)+rho(i+1,1,k  ))
            rhoy =        rho(i,1,k)
            rhoz = 0.5d0*(rho(i,1,k)+rho(i  ,1,k+1))
            u(i,0,k) = in_vel(i,k,1)*rhox
            v(i,0,k) = in_vel(i,k,2)*rhoy
            w(i,0,k) = in_vel(i,k,3)*rhoz
            p(i,0,k) = in_vel(i,k,4)
            u(i,1,k) = in_vel(i,k,5)*rhox
            v(i,1,k) = in_vel(i,k,6)*rhoy
            w(i,1,k) = in_vel(i,k,7)*rhoz
            p(i,1,k) = in_vel(i,k,8)
          enddo
        enddo
      else
        do k=1,n(3)
          do i=1,n(1)
            u(i,0,k) = in_vel(i,k,1)
            v(i,0,k) = in_vel(i,k,2)
            w(i,0,k) = in_vel(i,k,3)
            p(i,0,k) = in_vel(i,k,4)
            u(i,1,k) = in_vel(i,k,5)
            v(i,1,k) = in_vel(i,k,6)
            w(i,1,k) = in_vel(i,k,7)
            p(i,1,k) = in_vel(i,k,8)
          enddo
        enddo
      endif
    endif
    !
    call updthalo(n(1),n(2),n(3),nh_u,halo_u(ind1),ind1,u)
    call updthalo(n(1),n(2),n(3),nh_u,halo_u(ind2),ind2,u)
    !
    call updthalo(n(1),n(2),n(3),nh_u,halo_u(ind1),ind1,v)
    call updthalo(n(1),n(2),n(3),nh_u,halo_u(ind2),ind2,v)
    !
    call updthalo(n(1),n(2),n(3),nh_u,halo_u(ind1),ind1,w)
    call updthalo(n(1),n(2),n(3),nh_u,halo_u(ind2),ind2,w)
        !
    call updthalo(n(1),n(2),n(3),nh_u,halo_p(ind1),ind1,p)
    call updthalo(n(1),n(2),n(3),nh_u,halo_p(ind2),ind2,p)
    !
  case(0)
    ! Do nothing
  case default
    call flutas_error('ERROR: This inflow direction is not implemented ... yet!')
  end select

end subroutine apply_bound_vel

!_______________________________________________________________________________

subroutine add_inflow_force(n,dl,l,nh_p,nh_u,psi,rho,u,v,w,f_t12)
  use mod_param     , only: IBM_h
  !
  implicit none
  integer         , intent(in   ), dimension(3)                       :: n
  real(rp)        , intent(in   )                                     :: f_t12
  real(rp)        , intent(in   ), dimension(3)                       :: dl,l
  integer         , intent(in   )                                     :: nh_p,nh_u
  real(rp)        , intent(in   ), dimension(1-nh_p:,1-nh_p:,1-nh_p:) :: psi,rho
  real(rp)        , intent(inout), dimension(1-nh_u:,1-nh_u:,1-nh_u:) :: u,v,w
  real(rp)                                                            :: x,y,z,h_val,u_e
  real(rp)                                                            :: phiz,phiy,phix
  real(rp)                                                            :: rhoz,rhoy,rhox
  real(rp)                                                            :: dudt,dvdt,dwdt
  integer                                                             :: i,j,k
  real(rp)                                                            :: error,force
  !
  ! #####
  ! This works like a PI controller
  ! #####
  !
  if(in_dir.NE.0) then
    if(myid.eq.0) print*, "The mean inflow velocity is: ", bulk_vel
    if(myid.eq.0) print*, "The integral error is: ", int_start_bulk
    !
    error = in_force - bulk_vel
    int_start_bulk = int_start_bulk + f_t12*error
    !
    force = alpha_c*int_start_bulk + beta_c*error
    if(myid.eq.0) print*, "The driving force is: ", force
  endif
  !
  select case(in_dir)
  case(2)
    dvdt = force
    dudt = 0.0
    dwdt = 0.0
    do k=IBM_h + 1,n(3) !
      do j=1,n(2)
        do i=1,n(1)
          !
          rhoy = 0.5d0*(rho(i,j,k)+rho(i  ,j+1,k  ))
          v(i,j,k) = v(i,j,k) + dvdt*rhoy*f_t12
        enddo
      enddo
    enddo
  case(0)
    ! Do nothing
  case default
    call flutas_error('ERROR: This inflow direction is not implemented ... yet!')
  end select
end subroutine add_inflow_force

!#################################################
!###     End of:  Inflow subrouties            ###
!#################################################
!
!#################################################
!###     Start of:  Boiling subrouties         ###
!#################################################
!_______________________________________________________________________________
!
subroutine trigger_ouput(dt,time,trigger)
  !
  implicit none
  logical,intent(inout)                                      :: trigger
  real(rp),intent(in   )                                     :: dt
  real(rp),intent(in)                                        :: time
  real(rp)                                                   :: next_out_time
  !
  if(trigger) then
    trigger = .FALSE.
  else
    next_out_time = const_out_time*ceiling(time/const_out_time)
    !
    if(((time+dt).GT.(next_out_time)).AND.is_const_out_time) then
      ! Reducing the timestep to exactly match the requested output time did destabilize the simulation.
      ! The following two lines were commented out therefore.
      !dt = next_out_time-time
      !if(myid.eq.0) print*, 'To match output time, the timestep is reduced to ', dt
      trigger = .TRUE.
    endif
  endif

end subroutine trigger_ouput
!
!_______________________________________________________________________________
!
subroutine cmpt_nusselt(nh_p,dl,time,tmp,l,n)
  !
  use mod_param     , only: rho1,rho2,cp1,cp2,cond1,cond2,IBM_h,bctmp,sigma,rho1,rho2,gacc_z
  !
  implicit none
  !
  real(rp), intent(in ), dimension(3)                          :: dl,l
  integer , intent(in ), dimension(3)                          :: n
  real(rp), intent(in )                                        :: time
  integer , intent(in )                                        :: nh_p
  real(rp)                                                     :: nuss,delta_T,length,integral
  integer                                                      :: counter,i,j,k
  real(rp), intent(in ), dimension(1-nh_p:,1-nh_p:,1-nh_p:)    :: tmp
  real(rp),parameter                                           :: small_num = 2e-14_rp
  !
  delta_T = (abs(bctmp(0,3)-0.0)+small_num)        ! Watch out: This is problem dependent
  length = sqrt(sigma/(abs(gacc_z)*(rho2-rho1)+small_num))   ! Watch out: This is problem dependent
  !
  k = IBM_h + 1
  counter = 0
  integral = 0.0
  do j=1,n(2)
    do i=1,n(1)
      !
      integral = integral - (tmp(i,j,k)-tmp(i,j,k-1))/dl(3)
      counter = counter + 1
      !
    enddo
  enddo
  !
  call mpi_allreduce(MPI_IN_PLACE,counter   ,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
  call mpi_allreduce(MPI_IN_PLACE,integral  ,1,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
  !
  integral = integral/counter
  !
  nuss = (length/delta_T)*integral
  !
  if(myid.eq.0) print*, "The Nusselt number is: ", nuss, "   at time: ", time
  !
  return
end subroutine cmpt_nusselt

!_______________________________________________________________________________

subroutine cmpt_phase_change(nh_u,nh_p,nh_v,nh_d,dli,dl,n,src,tmp,psi,f_t12,eps_int,dzc,dzf,halo_v,rhocp,mu_t,rho,cond,nor)
  !
  use mod_common_mpi, only: MPI_PROC_NULL,left,right,front,back,top,bottom
  use mod_param     , only: IBM_h,cond1,cond2,h01,h02,cp1,cp2,pi,cbcpsi,bcpsi,rho1,rho2,is_outflow,heater_power
  use mod_bound     , only: boundp
  !
  implicit none
  !
  integer , intent(in   )                                        :: nh_u,nh_p,nh_v,nh_d  
  real(rp), intent(in   ), dimension(1-nh_d:)                    :: dzc,dzf
  integer , intent(in   ), dimension(3)                          :: halo_v 
  real(rp), intent(in   ), dimension(3)                          :: dli,dl
  integer , intent(in   ), dimension(3)                          :: n
  real(rp), intent(in   )                                        :: eps_int,f_t12
  integer                                                        :: i,j,k
  real(rp), intent(inout), dimension(1-nh_p:,1-nh_p:,1-nh_p:)    :: tmp,psi
  real(rp), intent(in   ), dimension(1-nh_p:,1-nh_p:,1-nh_p:)    :: rhocp,mu_t,rho,cond
  real(rp), intent(out  ), dimension(1-nh_p:,1-nh_p:,1-nh_p:)    :: src
  real(rp), intent(in   ), dimension(1-nh_v:,1-nh_v:,1-nh_v:,1:)   :: nor
  real(rp)                                                       :: aux,temp,t_char,deltaT,T_wall,dx_eq,dlmin
  real(rp)                                                       :: eps_eq,eps1,eps2,eps3,pc_factor_eq,theta_lambda
  integer                                                        :: length,x_start,x_stop,y_start
  integer                                                        :: y_stop,z_start,z_stop,min_dist
  real(rp)                                                       :: delta1,delta2,delta3
  real(rp), parameter                                            :: small = real(2.e-14,rp)
  !
  dlmin   = minval(dl(1:3))
  eps1    = min_eps_prime*dlmin/(dl(1))
  eps2    = min_eps_prime*dlmin/(dl(2))
  eps3    = min_eps_prime*dlmin/(dl(3))
  !_________________________________________________
  ! Phasechnage is switched off in outflow region
  !First give default values:
  x_start = 0
  x_stop  = 0
  y_start = 0
  y_stop  = 0
  z_start = 0
  z_stop  = 0
  min_dist = 3 ! Minimum distance of source term from any outflow BC
  ! Length = number of points within outflow region
  length = int((out_len)/dl(2))
  if(out_dir.NE.0) min_dist=max(min_dist,length)
  ! Taking the is_outflow argument into account too!
  if    (is_outflow(0,1).AND.(left .eq.MPI_PROC_NULL)) then
    x_start = min_dist
  elseif(is_outflow(1,1).AND.(right.eq.MPI_PROC_NULL)) then
    x_stop = min_dist
  elseif(is_outflow(0,2).AND.(front.eq.MPI_PROC_NULL)) then
    y_start = min_dist
  elseif(is_outflow(1,2).AND.(back .eq.MPI_PROC_NULL)) then
    y_stop = min_dist
  elseif(is_outflow(0,3)) then
    z_start = min_dist
  elseif(is_outflow(1,3)) then
    z_stop = min_dist
  endif
  !_________________________________________________
  !
  !
  !Compute the source term
  do k=IBM_h+1+z_start,n(3)-z_stop
    do j=1+y_start,n(2)-y_stop
     do i=x_start+1,n(1)-x_stop
       !
       ! Phase change timescale: Assumption that the heat conduction of the liquid dictates the phase change
       t_char = (dl(1)**2)/(cond2/(rho2*cp2)+(mu_t(i,j,k)/(rho(i,j,k)*pr_turb)))
       ! deltaT -> Change of Temperature during dt due to phase change
       !           Note: Instantanous phase change would mean deltaT = (T_saturation - T). 
       !                 But for the phase field method, additional terms are needed
       !
       !              term1           term2                     term3 
       delta1 =  (f_t12/t_char)   *  pc_factor1   *    (1.0/(((rho1*cp1)/(rho2*cp2))*psi(i,j,k)+(1.0_rp-psi(i,j,k)))) &
                *   (0.5_rp+0.5_rp*cos(2.0_rp*pi*abs(psi(i,j,k)-0.5_rp)))         *    1.0/eps1
       !                              term4                                             term6
       !
       ! Phase change timescale: Assumption that the heat conduction of the liquid dictates the phase change
       t_char = (dl(2)**2)/(cond2/(rho2*cp2)+(mu_t(i,j,k)/(rho(i,j,k)*pr_turb)))
       ! deltaT -> Change of Temperature during dt due to phase change
       !           Note: Instantanous phase change would mean deltaT = (T_saturation - T). 
       !                 But for the phase field method, additional terms are needed
       !
       !              term1           term2                     term3 
       delta2 =  (f_t12/t_char)   *  pc_factor2   *    (1.0/(((rho1*cp1)/(rho2*cp2))*psi(i,j,k)+(1.0_rp-psi(i,j,k)))) &
                *   (0.5_rp+0.5_rp*cos(2.0_rp*pi*abs(psi(i,j,k)-0.5_rp)))         *    1.0/eps2
       !                              term4                                             term6
       !
       ! Phase change timescale: Assumption that the heat conduction of the liquid dictates the phase change
       t_char = (dl(3)**2)/(cond2/(rho2*cp2)+(mu_t(i,j,k)/(rho(i,j,k)*pr_turb)))
       ! deltaT -> Change of Temperature during dt due to phase change
       !           Note: Instantanous phase change would mean deltaT = (T_saturation - T). 
       !                 But for the phase field method, additional terms are needed
       !
       !              term1           term2                     term3 
       delta3 =  (f_t12/t_char)   *  pc_factor3   *    (1.0/(((rho1*cp1)/(rho2*cp2))*psi(i,j,k)+(1.0_rp-psi(i,j,k)))) &
                *   (0.5_rp+0.5_rp*cos(2.0_rp*pi*abs(psi(i,j,k)-0.5_rp)))         *    1.0/eps3
       !                              term4                                             term6
       !  term1 ->  relate phase change speed to current timestep size
       !  term2 ->  model constant (problem independent), 0.7
       !  term3 ->  weighting the phase change speed by rhocp to have faster phase change on gas side
       !            than on the liquid side (Assuming rho2*cp2 > rho1*cp1)
       !  term4 ->  limit phase-change to "interface region".  
       !            where term4(psi=0) = 0;  term4(psi=0.5) = 1;  term4(psi=1.0) = 0
       !  term5 ->  local difference to saturation temperature. i.e. (T - T_saturation).
       !            Note: Reference temperature for tmp(:,:,:) is the saturation temperature already!
       !  term6 ->  adjusting the phase chnage speed (i.e. the temperature increment) to the interface thickness
       !            term5
       deltaT = -   tmp(i,j,k)*(delta1*abs(small+nor(i,j,k,1))+delta2*abs(small+nor(i,j,k,2))+delta3*abs(small+nor(i,j,k,3))) &
                              /(abs(small+nor(i,j,k,1))+abs(small+nor(i,j,k,2))+abs(small+nor(i,j,k,3)))
       T_wall = 0.0
       if((k.EQ.(IBM_h+1)).AND.(IBM_h.NE.0)) T_wall = (cond(i,j,k)*tmp(i,j,k)+cond(i,j,k-1)*tmp(i,j,k-1)) &
                                  /(cond(i,j,k)           +cond(i,j,k-1)             ) &
                                  +(heater_power*dl(3)*heater_ind(i,j))                &
                                  /(2.0_rp*(cond(i,j,k)+cond(i,j,k-1)))
       !
       t_char = (dl(3)**2)/(cond2/(rho2*cp2)+(mu_t(i,j,k)/(rho(i,j,k)*pr_turb)))
       !
       theta_lambda = (cond1/cond2)*psi(i,j,k)+(1.0_rp-psi(i,j,k))
       !
       deltaT = deltaT - (f_t12/t_char)   *    (1.0/(((rho1*cp1)/(rho2*cp2))*psi(i,j,k)+(1.0_rp-psi(i,j,k))))*  theta_lambda  &  !
                                          *    (0.5_rp+0.5_rp*cos(2.0_rp*pi*abs(psi(i,j,k)-0.5_rp)))     *    2.0*  T_wall
       !
       ! Step 1: Update temperature field
       tmp(i,j,k)    = tmp(i,j,k) + deltaT
       !
       ! Step 2: Calculate source term that corresponds to deltaT (from energy equation) !
       src(i,j,k)    = -rhocp(i,j,k)*(deltaT/f_t12)/((cp1-cp2)*tmp(i,j,k)+(h01 - h02))
       !
       !src(i,j,k)    = -0.1*((psi(i,j,k+1)-psi(i,j,k-1))/(2.0*dl(3)))/(1.0/rho1-1.0/rho2)
       !tmp(i,j,k)    = tmp(i,j,k) - f_t12*((cp1-cp2)*tmp(i,j,k)+(h01 - h02))*src(i,j,k)/rhocp(i,j,k)
       !
       ! Step 3: Update phase field 
       psi(i,j,k)    = psi(i,j,k) + src(i,j,k)*f_t12/rho1
       !
     enddo
    enddo
   enddo
   !
   ! Do MPI communication and BC (same as psi for now...)
   !
   call boundp(cbcpsi_n,n,bcpsi_n,nh_d,nh_v,halo_v,dl,dzc,dzf,src)
   !
  return
end subroutine cmpt_phase_change

!_______________________________________________________________________________

subroutine surft_src_FEF(nx,ny,nz,nh_d,nh_u,f_t12,dxi,dyi,dzi,dzci,psi,u,v,w,l,dl,dzc,dzf,halo_v,nh_v,n,FE_type,ls)
  ! Surface tension force based on free energy functional (FEF) to void calculating the curvature
  ! for reducing spourious currents
  !
  ! See [1]  Huang et al. (2020), Eq. 42, DOI: https://doi.org/10.1016/j.jcp.2020.109718
  !     [2]  Badillo      (2013), Eq.  2, DOI: https://doi.org/10.1115/IMECE2013-65925 
  !     [3]  Mirjaili     (2023), Eq.  4, DOI: https://doi.org/10.1016/j.jcp.2022.111795
  !
  ! For discretisation details see [1]
  !
  ! Note: the input velocities are actually rho.u and the output as well
  !
  use mod_param,    only: sigma,IBM_h,eps_int
  use mod_bound,    only: boundp_IBM
  !
  implicit none 
  !
  integer , intent(in   )                                     :: nx,ny,nz
  integer , intent(in   ), dimension(3)                       :: n
  integer , intent(in   )                                     :: nh_d,nh_u,nh_v
  real(rp), intent(in   ), dimension(3)                       :: l,dl
  real(rp), intent(in   )                                     :: f_t12
  real(rp), intent(in   )                                     :: dxi,dyi,dzi
  real(rp), intent(in   ), dimension(1-nh_d:)                 :: dzci
  real(rp), intent(in   ), dimension(0:,0:,0:)                :: psi,ls
  real(rp), intent(inout), dimension(1-nh_u:,1-nh_u:,1-nh_u:) :: u,v,w
  real(rp), allocatable  , dimension(:,:,:)                   :: auxfield
  real(rp), intent(in   ), dimension(1-nh_d:)                 :: dzc,dzf
  integer , intent(in   ), dimension(3)                       :: halo_v 
  character(len=3),  intent(in   )                            :: FE_type
  character(len=1), parameter, dimension(0:1,3)               :: cbcaux= reshape((/'N','N','N' &
                                                                                  ,'N','N','N'/),shape(cbcaux))
  real(rp)        , parameter, dimension(0:1,3)               ::  bcaux= reshape((/0,0,0,0,0,0/),shape(bcaux))
  !
  real(rp) :: phix,phiy,phiz,z,h_val,switch,auxx,auxy,auxz,dmudx,dmudy,dmudz,mag_grad_ls_sq
  integer  :: i,j,k,ip,jp,kp
  !@cuf attributes(managed) ::  psi, dzci, u, v, w
  !
  !$acc kernels
  allocate(auxfield(0:nx+1,0:ny+1,0:nz+1))
  ! switch between sonservative and non conservative version (see [1]):
  select case(FE_type)
  case('nor')  !-> the normal type
    ! First: calculate first term on auxiluary field
    !
    do k=IBM_h+1,nz
      do j=1,ny
        do i=1,nx
          auxfield(i,j,k) = 6.0*sigma*(                                                   &
            -      eps_int*( (psi(i+1,j  ,k  )-2.0*psi(i,j,k)+psi(i-1,j  ,k  ))*(dxi**2)  &
                            +(psi(i  ,j+1,k  )-2.0*psi(i,j,k)+psi(i  ,j-1,k  ))*(dyi**2)  &
                            +(psi(i  ,j  ,k+1)-2.0*psi(i,j,k)+psi(i  ,j  ,k-1))*(dzi**2)) &
            +(1.0/eps_int)*(1.0-psi(i,j,k))*psi(i,j,k)*(1.0-2.0*psi(i,j,k)))
        enddo
      enddo
    enddo
    !
    ! Since the face values of auxfield are required we need to apply 
    ! bouandry conditions and do mpi communication (see [1] page 7)
    !
    call boundp_IBM(cbcaux,n,bcaux,nh_d,nh_v,halo_v,dl,dzc,dzf,auxfield,IBM_h)
    !
    switch = 1.0
    if(out_dir.eq.0) switch = 0.0
    do k=IBM_h+1,nz
      do j=1,ny
        do i=1,nx
          !
          ip = i+1   
          jp = j+1   
          kp = k+1  
          !
          auxx   = 0.5_rp*(auxfield(ip,j,k)+auxfield(i,j,k))
          auxy   = 0.5_rp*(auxfield(i,jp,k)+auxfield(i,j,k))
          auxz   = 0.5_rp*(auxfield(i,j,kp)+auxfield(i,j,k))
#if defined(_TWOD)
          auxx = 0.0
#endif
          !
          z = (k+ijk_start(3))*dl(3)
          h_val    = switch*(2.0/(1.0+exp(4.0*(l(3)-z+0.5*dl(3))/out_len)))
          u(i,j,k) = u(i,j,k) + (1.0-h_val)*f_t12*(auxx*((psi(ip,j,k)-psi(i,j,k))*dxi))
          v(i,j,k) = v(i,j,k) + (1.0-h_val)*f_t12*(auxy*((psi(i,jp,k)-psi(i,j,k))*dyi))
          h_val    = switch*(2.0/(1.0+exp(4.0*(l(3)-z)/out_len)))
          w(i,j,k) = w(i,j,k) + (1.0-h_val)*f_t12*(auxz*((psi(i,j,kp)-psi(i,j,k))*dzi))
          !
          !
        enddo
      enddo
    enddo
  case('con')  !-> the conservative type (see [1])
    !
    ! First: calculate first term on auxiluary field for conservative version
    !
    do k=IBM_h+1,nz
      do j=1,ny
        do i=1,nx
          auxfield(i,j,k) = 6.0*sigma*(                                                   &
            -      eps_int*( (psi(i+1,j  ,k  )-2.0*psi(i,j,k)+psi(i-1,j  ,k  ))*(dxi**2)  &
                            +(psi(i  ,j+1,k  )-2.0*psi(i,j,k)+psi(i  ,j-1,k  ))*(dyi**2)  &
                            +(psi(i  ,j  ,k+1)-2.0*psi(i,j,k)+psi(i  ,j  ,k-1))*(dzi**2)))
        enddo
      enddo
    enddo
    !
    ! Since the face values of auxfield are required we need to apply 
    ! bouandry conditions and do mpi communication (see [1] page 7)
    !
    call boundp_IBM(cbcaux,n,bcaux,nh_d,nh_v,halo_v,dl,dzc,dzf,auxfield,IBM_h)
    !
    switch = 1.0
    if(out_dir.eq.0) switch = 0.0
    do k=IBM_h+1,nz
      do j=1,ny
        do i=1,nx
          !
          ip = i+1   
          jp = j+1   
          kp = k+1  
          !
          auxx   = 0.5_rp*(auxfield(ip,j,k)+auxfield(i,j,k))
          auxy   = 0.5_rp*(auxfield(i,jp,k)+auxfield(i,j,k))
          auxz   = 0.5_rp*(auxfield(i,j,kp)+auxfield(i,j,k))
          dmudx  = (6.0*sigma/eps_int)*((0.5*(psi(i+1,j,k)*(1-psi(i+1,j,k)))**2 &
                                     -   0.5*(psi(i  ,j,k)*(1-psi(i  ,j,k)))**2)*dxi)
          dmudy  = (6.0*sigma/eps_int)*((0.5*(psi(i,j+1,k)*(1-psi(i,j+1,k)))**2 &
                                     -   0.5*(psi(i,j  ,k)*(1-psi(i,j  ,k)))**2)*dyi)
          dmudz  = (6.0*sigma/eps_int)*((0.5*(psi(i,j,k+1)*(1-psi(i,j,k+1)))**2 &
                                     -   0.5*(psi(i,j,k  )*(1-psi(i,j,k  )))**2)*dzi)
#if defined(_TWOD)
          auxx = 0.0
#endif
          !
          z = (k+ijk_start(3))*dl(3)
          h_val    = switch*(2.0/(1.0+exp(4.0*(l(3)-z+0.5*dl(3))/out_len)))
          u(i,j,k) = u(i,j,k) + (1.0-h_val)*f_t12*(dmudx + auxx*((psi(ip,j,k)-psi(i,j,k))*dxi))
          v(i,j,k) = v(i,j,k) + (1.0-h_val)*f_t12*(dmudy + auxy*((psi(i,jp,k)-psi(i,j,k))*dyi))
          h_val    = switch*(2.0/(1.0+exp(4.0*(l(3)-z)/out_len)))
          w(i,j,k) = w(i,j,k) + (1.0-h_val)*f_t12*(dmudz + auxz*((psi(i,j,kp)-psi(i,j,k))*dzi))
          !
          !
        enddo
      enddo
    enddo
  case('brn') !   -> the formulation used in Brown (2024): https://arxiv.org/pdf/2412.16843
    ! First: calculate first term on auxiluary field
    !
    do k=IBM_h+1,nz
      do j=1,ny
        do i=1,nx
          mag_grad_ls_sq  = ((ls(i+1,j  ,k  )-ls(i-1,j  ,k  ))*(0.5*dxi))**2  &
                           +((ls(i  ,j+1,k  )-ls(i  ,j-1,k  ))*(0.5*dyi))**2  &
                           +((ls(i  ,j  ,k+1)-ls(i  ,j  ,k-1))*(0.5*dzi))**2
          !
          auxfield(i,j,k) = 6.0*sigma*( &
            - (1.0-psi(i,j,k))*psi(i,j,k)*( (ls(i+1,j  ,k  )-2.0*ls(i,j,k)+ls(i-1,j  ,k  ))*(dxi**2)  &
                                           +(ls(i  ,j+1,k  )-2.0*ls(i,j,k)+ls(i  ,j-1,k  ))*(dyi**2)  &
                                           +(ls(i  ,j  ,k+1)-2.0*ls(i,j,k)+ls(i  ,j  ,k-1))*(dzi**2)) &
            +(1.0/eps_int)*(1.0-psi(i,j,k))*psi(i,j,k)*(1.0-2.0*psi(i,j,k))*(1.0-mag_grad_ls_sq))
        enddo
      enddo
    enddo
    !
    ! Since the face values of auxfield are required we need to apply 
    ! bouandry conditions and do mpi communication (see [1] page 7)
    !
    call boundp_IBM(cbcaux,n,bcaux,nh_d,nh_v,halo_v,dl,dzc,dzf,auxfield,IBM_h)
    !
    switch = 1.0
    if(out_dir.eq.0) switch = 0.0
    do k=IBM_h+1,nz
      do j=1,ny
        do i=1,nx
          !
          ip = i+1   
          jp = j+1   
          kp = k+1  
          !
          auxx   = 0.5_rp*(auxfield(ip,j,k)+auxfield(i,j,k))
          auxy   = 0.5_rp*(auxfield(i,jp,k)+auxfield(i,j,k))
          auxz   = 0.5_rp*(auxfield(i,j,kp)+auxfield(i,j,k))
#if defined(_TWOD)
          auxx = 0.0
#endif
          !
          z = (k+ijk_start(3))*dl(3)
          h_val    = switch*(2.0/(1.0+exp(4.0*(l(3)-z+0.5*dl(3))/out_len)))
          u(i,j,k) = u(i,j,k) + (1.0-h_val)*f_t12*(auxx*((psi(ip,j,k)-psi(i,j,k))*dxi))
          v(i,j,k) = v(i,j,k) + (1.0-h_val)*f_t12*(auxy*((psi(i,jp,k)-psi(i,j,k))*dyi))
          h_val    = switch*(2.0/(1.0+exp(4.0*(l(3)-z)/out_len)))
          w(i,j,k) = w(i,j,k) + (1.0-h_val)*f_t12*(auxz*((psi(i,j,kp)-psi(i,j,k))*dzi))
          !
          !
        enddo
      enddo
    enddo
  case default
    call flutas_error('ERROR in surface tension: This FE method is not implemented ... yet!')
  end select
  !$acc end kernels
  !
  deallocate(auxfield)
  return
end subroutine surft_src_FEF

!_______________________________________________________________________________

subroutine copy_fields(n,nh_p,in_field,out_field) 
  implicit none
  !
  integer , intent(in ), dimension(3)                       :: n
  integer , intent(in )                                     :: nh_p
  real(rp), intent(in ), dimension(1-nh_p:,1-nh_p:,1-nh_p:) :: in_field
  real(rp), intent(out), dimension(1-nh_p:,1-nh_p:,1-nh_p:) :: out_field
  integer                                                   :: i,j,k
  !
  !
  do k=1-nh_p,n(3)+nh_p
    do j=1-nh_p,n(2)+nh_p
      do i=1-nh_p,n(1)+nh_p
        !
        out_field(i,j,k)    =   in_field(i,j,k)
      enddo
    enddo
  enddo
  !
  return
end subroutine copy_fields

!_______________________________________________________________________________

subroutine add_a_to_b(n,nh_p,field_a,field_b) 
  ! add array a (input) to array b (input & output)
  ! MPI communication required after this step!
  implicit none
  !
  integer , intent(in   ), dimension(3)                       :: n
  integer , intent(in   )                                     :: nh_p
  real(rp), intent(in   ), dimension(1-nh_p:,1-nh_p:,1-nh_p:) :: field_a
  real(rp), intent(inout), dimension(1-nh_p:,1-nh_p:,1-nh_p:) :: field_b
  integer                                                     :: i,j,k
  !
  !
  do k=1,n(3)
    do j=1,n(2)
      do i=1,n(1)
        !
        field_b(i,j,k)    =    field_b(i,j,k)  +  field_a(i,j,k)
      enddo
    enddo
  enddo
  !
  return
end subroutine add_a_to_b

!_______________________________________________________________________________

subroutine mult_by_const(n,nh_p,field,const_val)
  ! multiply array by a constant value, do also halo variables so that no MPI
  ! communication is needed after this operation
  implicit none
  !
  integer , intent(in   ), dimension(3)                       :: n
  integer , intent(in   )                                     :: nh_p
  real(rp), intent(in   )                                     :: const_val
  real(rp), intent(inout), dimension(1-nh_p:,1-nh_p:,1-nh_p:) :: field
  !
  integer  :: i,j,k
  !
  !
  do k=1-nh_p,n(3)+nh_p
    do j=1-nh_p,n(2)+nh_p
      do i=1-nh_p,n(1)+nh_p
        !
        field(i,j,k)    =   const_val*field(i,j,k)
      enddo
    enddo
  enddo
  !
  return
end subroutine mult_by_const

!_______________________________________________________________________________

subroutine static_contact_angle_m_ibm(n,dl,nh_v,eps,theta,ap,ap_bound,ls_bound,IBM_h,kappa)
  !                             
  ! To do: we can treat the term multiplying grad_rho implicitly
  ! 
  implicit none
  !
  integer , intent(in ), dimension(3)                       :: n
  real(rp), intent(in ), dimension(3)                       :: dl
  integer , intent(in )                                     :: nh_v,IBM_h
  real(rp), intent(in )                                     :: eps,theta
  real(rp), intent(in ), dimension(1-nh_v:,1-nh_v:,1-nh_v:) :: ap,kappa
  real(rp), intent(out), dimension( 1:, 1:)                 :: ap_bound,ls_bound
  real(rp)                                                  :: f,mv_wall,ap_wall,tan_abs,tan_x,tan_y,clip_psi
  integer                                                   :: i,j,k,kpp,kp,ip,im,jp,jm
  real(rp), parameter                                       :: small = real(2.e-16,rp)
  real(rp)                                                  :: dim_con,theta_prime,ls2,h,r2,tester,phi_wall,ls_wall
  !
#if defined(_TWOD)
  dim_con = 1.0_rp
#else
  dim_con = 2.0_rp
#endif
  !
  theta_prime = theta
  !
  do j=1,n(2)
    do i=1,n(1)
      !
      ! written by Lorenz Weber 19th DEC 2023 (lorenzw@kth.se)
      ! Further reading in the Literature:
      ! - equation 36 and 37 in Zarareh (2021): Improving the staircase approximation for wettability implementation of phase-field model: Part 1 – Static contact angle
      ! - equation 42, 43 and 44 in Liang (2019): Lattice Boltzmann method for contact-line motion of binary fluids with high density ratio
      ! - equation 18 in Ding (2007): Wetting condition in diffuse interface simulations of contact line motion
      k   =   IBM_h
      kp  = k + 1
      kpp = k + 2
      ip  = i + 1
      im  = i - 1
      jp  = j + 1
      jm  = j - 1
      !
      ls2 = eps*log( (ap(i,j,kp)+small)/(1.0_rp-ap(i,j,kp)+small) )
      r2  = dim_con/(kappa(i,j,kp)+small)
      h   = (r2+ls2)*sind(theta-90.0_rp)
      tester = abs(sind(theta-90.0_rp)*(1+ls2*kappa(i,j,kp)/dim_con))
      !
      if((kappa(i,j,kp).LE.0.0_rp).or.(ap(i,j,kp).GE.0.999_rp).or.(ap(i,j,kp).LE.0.001_rp).or.(tester.GE.1.0_rp)) then
        theta_prime = theta
      elseif(r2.LT.h) then
        theta_prime = 179.0_rp
      else
        theta_prime = 90.0_rp + asind(sind(theta-90.0_rp)*(1+ls2*kappa(i,j,kp)/dim_con))
      endif
      !
      !
      ls_wall = 0.0_rp
      phi_wall = 0.0_rp
      !
      if(.NOT.((ap(i,j,kp).LE.0.001_rp).or.(ap(i,j,kp).GE.0.999_rp))) then
        ls_wall = ls2 - 0.5*dl(3)*cosd(180.0_rp -theta_prime)
        phi_wall = 0.5_rp*(1._rp+tanh((ls_wall)/(2._rp*eps) ))
      endif
      !
      ap_wall  = 0.5*(ap(i,j,k)+ap(i,j,kp))
      f        = ap_wall*(1.d0-ap_wall) ! non-linear term treated explicitly
      !
      tan_x         = 1.5_rp*((ap(ip,j ,kp)-ap(im,j ,kp))/dl(1)/(2.0_rp)) - 0.5_rp*((ap(ip,j ,kpp)-ap(im,j ,kpp))/dl(1)/(2.0_rp))
      tan_y         = 1.5_rp*((ap(i ,jp,kp)-ap(i ,jm,kp))/dl(2)/(2.0_rp)) - 0.5_rp*((ap(i ,jp,kpp)-ap(i ,jm,kpp))/dl(2)/(2.0_rp))
      tan_abs       = sqrt(tan_x**2+tan_y**2)
      ap_bound(i,j) =  ap(i,j,kp) +  min(0.0_rp,max(-0.3*ap(i,j,kp),(tan_abs*TAND(90.0_rp - theta_prime)*dl(3)))) !(0.5*(dl(3)/eps)*cosd(theta)*f)  ! 0.99* + 0.01*phi_wall ! 
      clip_psi      = min(max(0.0_rp,ap_bound(i,j)),1.0_rp) 
      ls_bound(i,j) =   eps*log( (clip_psi+small)/(1.0_rp-clip_psi+small) )   !ls_wall ! 
      !
    enddo
  enddo
  return
end subroutine static_contact_angle_m_ibm

!_______________________________________________________________________________

subroutine adjust_bound_vel(nx,ny,nz,dxi,dyi,dzi,nh_d,nh_u,dzci,dzfi,psi,dl,dzc,dzf,halo_u,f_t12,time)  
    !only works for the z velocity at the top boundary and onlz with a single bubble!
    ! this subroutine works like a control lopp with only the proportinal part :)
    use mod_param    ,only: rho1,rho2,IBM_h,no_outflow,cbcvel,bcvel
    !
    implicit none
    !
    integer , intent(in   )                                     :: nx,ny,nz
    integer , intent(in   ), dimension(3)                       :: halo_u
    real(rp), intent(in   ), dimension(3)                       :: dl
    integer , intent(in   )                                     :: nh_d,nh_u    
    real(rp), intent(in   ), dimension(1-nh_d:)                 :: dzc,dzf
    real(rp), intent(in   )                                     :: dxi,dyi,dzi,f_t12,time
    real(rp), intent(in   ), dimension(1-nh_d:                ) :: dzci,dzfi
    real(rp), intent(in   ), dimension(     0:,     0:,     0:) :: psi
    !
    real(rp) :: rhox,rhoy,rhoz,cog,x,y,z,cog_target,p_factor,bubble_vel,weight,tolerance
    integer  :: im,ip,jm,jp,km,kp,i,j,k
    !
    !@cuf attributes(managed) :: u, v, w, mu, rho, dzci, dzfi, dudt, dvdt, dwdt
    !    
    !$acc kernels
    cog_target = nuc_sites(1,3)
    ! p_factor = max rise vel ~0.15m/s and max. allowed dispacement 0.15mm => 1000.0
    p_factor = 500.0
    tolerance = 0.01*1e-3
    if(cbcvel(1,3,3).NE.'D') then
      !call flutas_error('ERROR: Wrong BC type for adaptive input velocity!')
    endif
    cog = 0.0
    weight = 0.0
    do k=1+IBM_h,nz !
      z = (k+ijk_start(3)-0.5_rp)*dl(3)
      do j=1,ny
        y = (j+ijk_start(2)-0.5_rp)*dl(2)
        do i=1,nx
          x = (i+ijk_start(1)-0.5_rp)*dl(1)
          !
          ip = i+1   
          jp = j+1   
          kp = k+1  
          cog = cog + psi(i,j,k)*z
          weight = weight + psi(i,j,k)
          !
        enddo
      enddo
    enddo
    call mpi_allreduce(MPI_IN_PLACE,cog,1,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
    call mpi_allreduce(MPI_IN_PLACE,weight,1,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
    cog = cog/weight
    bubble_vel = -max((cog-cog_target-tolerance),0.0)*p_factor
    bcvel(1,3,3) = bubble_vel
    if(myid.eq.0) print*, "Inflow velocity: ", bubble_vel, time, f_t12
    !
    !$acc end kernels 
    !$OMP END PARALLEL DO
    !
    return
  end subroutine adjust_bound_vel

!_______________________________________________________________________________


subroutine update_gamma_v(nh_u,nh_d,dli,dzci,dzfi,dxi,dyi,dzi,n,u,v,w,u_reg,eps_int,dtmax,gamma_v,param_update)
  !
  ! Update gamma_v if cdi restriction is not the restrictive parameter
  !
  use mod_param     , only: IBM_h,constant_dt,dt_input,cfl
  !
  ! [1] Jain et al      (2020), https://doi.org/10.1016/j.jcp.2020.109606
  ! [2] Jain            (2022), https://doi.org/10.1016/j.jcp.2022.111529
  ! [3] Mirjalili et al (2020), https://doi.org/10.1016/j.jcp.2019.109006
  !
  implicit none
  !
  real(rp), intent(in   ), dimension(3)                        :: dli
  integer , intent(in   ), dimension(3)                        :: n
  integer , intent(in   )                                      :: nh_u,nh_d
  logical , intent(in   )                                      :: param_update
  real(rp), intent(in   ), dimension(1-nh_d:)                  :: dzci,dzfi
  real(rp), intent(in   )                                      :: dxi,dyi,dzi
  real(rp), intent(in   )                                      :: u_reg
  real(rp), intent(inout)                                      :: gamma_v,eps_int
  real(rp), intent(inout)                                      :: dtmax
  integer                                                      :: i,j,k,ip,im,jp,jm,kp,km
  real(rp), intent(in   ), dimension(1-nh_u:,1-nh_u:,1-nh_u:)  :: u,v,w
  real(rp)                                                     :: gradu,gradv,gradw,mingrad
  real(rp)                                                     :: dlmin,dt_cdi,gamma_old,umax
  real(rp)                                                     :: max_change,best_dt_eps,eps_old
  !
  dlmin         = min(1._rp/dxi,1._rp/dyi,1._rp/dzi)
  dlmin         = min(dlmin,minval(1._rp/dzfi(:))) ! minimum of dzf is an estimate on the safe side  
  !
  gamma_old     = gamma_v
  eps_old       = eps_int
  umax          = (u_reg/gamma_v)
  best_dt_eps   = min_eps_prime*dlmin !max(((1.0+1.0/min_gamma)*0.5*dlmin),min_eps_prime*dlmin)
  max_change    = 1.2
  !
  !______________________________________________________________
  ! First the minimum divergence in the flowfield is determined
  !
  mingrad = 0.0
  !
  do k=1+IBM_h,n(3)
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
        gradu          =   0.d0
#else
        gradu          =   ( u(im,j ,k ) - u(i ,j ,k )  ) * dli(1)
#endif
        !
        gradv          =   ( v(i ,jm,k ) - v(i ,j ,k )  ) * dli(2)
        gradw          =   ( w(i ,j ,km) - w(i ,j ,k )  ) * dli(3)
        mingrad        =   min(mingrad,gradu,gradv,gradw)
        !
      enddo
    enddo
  enddo
  call mpi_allreduce(MPI_IN_PLACE,mingrad,1,MPI_REAL_RP,MPI_MIN,MPI_COMM_WORLD,ierr)
  !
  !______________________________________________________________
  !  now calculate maximum timestep with the lowest allowed regularistaion speed
  !  with the criterion from [1]
  dt_cdi = 1.0/max((((min_gamma*umax*6*best_dt_eps)/(dlmin**2))-mingrad),0.0)
  !
  !______________________________________________________________
  !For adaptive timestepping
  !
  if(.NOT.constant_dt) then
    ! If CDI is restrictive
    if(dt_cdi.LE.dtmax) then
      gamma_v = min_gamma
      eps_int = best_dt_eps
      dt_cdi = dt_cdi  !divide by cfl because cfl is applied after the subroutine separately
      dtmax  = min(dtmax,dt_cdi)
      if(myid.eq.0) print*, "time-step restricted by CDI:",(dt_cdi)
    else ! If CDI not restrictive
      !gamma_v = max(min_gamma,((cfl/dtmax)+mingrad)*dlmin**2/(umax*6*eps_int))
      gamma_v = max(min_gamma,(((1.0/dtmax)+mingrad)*dlmin/(umax*3.0))-1.0) ! Theoretical max gamma for given timestep size
      eps_int = min_eps_prime*dlmin! max(min_eps_prime*dlmin,((1.0+1.0/gamma_v)*0.5*dlmin)) ! Minumim eps to corresponding gamma
      gamma_v = max(min_gamma,((1.0/dtmax)+mingrad)*dlmin**2/(umax*6*eps_int)) !adjust gamma again if eps is restricted
      !
    endif
  endif
  !
  !______________________________________________________________
  !For constant timestepping
  !
  if(constant_dt) then
    ! If CDI is restrictive
    if(dt_cdi.LE.dt_input) then
      gamma_v = min_gamma
      eps_int = best_dt_eps
      dt_cdi = dt_cdi  !divide by cfl although timestep size is obiously not updated
      dtmax  = min(dtmax,dt_cdi)
      if(myid.eq.0) print*, "WARNING: time-step restricted by CDI:",(dt_cdi)
    else ! If CDI not restrictive
      !gamma_v = max(min_gamma,((cfl/dt_input)+mingrad)*dlmin**2/(umax*6*eps_int))
      gamma_v = max(min_gamma,(((1.0/dt_input)+mingrad)*dlmin/(umax*3.0))-1.0) ! Theoretical max gamma for given timestep size
      eps_int = max(min_eps_prime*dlmin,((1.0+1.0/gamma_v)*0.5*dlmin)) ! Minumim eps to corresponding gamma
      gamma_v = max(min_gamma,((1.0/dt_input)+mingrad)*dlmin**2/(umax*6*eps_int)) !adjust gamma again if eps is restricted
      !
    endif
  endif
  !______________________________________________________________
  ! Let gamma increase not too fast for stability
  if(gamma_v.GT.max_change*gamma_old) then
    gamma_v = max_change*gamma_old
    eps_int = min_eps_prime*dlmin !max(min_eps_prime*dlmin,((1.0+1.0/gamma_v)*0.5*dlmin))
  endif
  !______________________________________________________________
  ! If update switch is FALSE (e.g. for simulation start) eps_int and gamm_v are not updated
  if(.NOT.param_update) then
    gamma_v     = gamma_old
    eps_int     = eps_old
  endif
  !______________________________________________________________
  ! Logfile Info
  if(myid.eq.0) print*, "Gamma and epsilon are now:" , gamma_v, eps_int/dlmin
  !
  return
end subroutine update_gamma_v

!_______________________________________________________________________________

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

!_______________________________________________________________________________

real(rp) function surface_integral(nx,ny,nz,idir,nh_u,dx,dy,dz,u,v,w)
  ! facecentered!
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
  surface_integral = integral
  !
  return
end function surface_integral
!#################################################
!###     End of:  Boiling subrouties           ###
!#################################################
!
!#################################################
!###     Start of:  Pressure subrouties        ###
!#################################################
!  
subroutine pres_tw_src(nx,ny,nz,dxi,dyi,dzi,nh_d,nh_u,dzci,rho0i,f_t12,f_t12_o,p,dp,rho,u,v,w,stage,istep)
  !
  ! Stage 1 : Add              old pressure       without splitting      to momentum
  ! Stage 2 : Add extrapolated pressure increment with density splitting to velocity
  ! Stage 3 : Add predicted    pressure increment with density splitting to velocity
  use mod_param     , only: IBM_h
  implicit none 
  !
  integer , intent(in   )                                     :: nx,ny,nz,istep
  real(rp), intent(in   )                                     :: dxi,dyi,dzi
  integer , intent(in   )                                     :: nh_d,nh_u
  real(rp), intent(in   ), dimension(1-nh_d:)                 :: dzci
  real(rp), intent(in   )                                     :: rho0i,f_t12,f_t12_o
  real(rp), intent(in   ), dimension(0:,0:,0:)                :: p,dp
  real(rp), intent(in   ), dimension(0:,0:,0:)                :: rho
  real(rp), intent(inout), dimension(1-nh_u:,1-nh_u:,1-nh_u:) :: u,v,w
  integer , intent(in   )                                     :: stage
  !
  real(rp) :: rhox,rhoy,rhoz,rhoxi,rhoyi,rhozi
  real(rp) :: f2
  integer  :: i,j,k,ip,jp,kp
  !@cuf attributes(managed) :: dzci, p, pold, rho, u, v, w
  !
  if(istep.eq.1) then ! we use 1st order Euler in the first time-step (istep.eq.1)
    f2 = 0._rp
  else
    f2 = (f_t12/f_t12_o)
  endif
  !
  !$acc kernels
  select case(stage)
    case(1)
      do k=1+IBM_h,nz
        do j=1,ny
          do i=1,nx
            !
            ip = i+1
            jp = j+1
            kp = k+1
            !
            u(i,j,k) = u(i,j,k) + f_t12*( - ( p(ip,j,k)-p(i,j,k) )*dxi     )
            v(i,j,k) = v(i,j,k) + f_t12*( - ( p(i,jp,k)-p(i,j,k) )*dyi     )
            w(i,j,k) = w(i,j,k) + f_t12*( - ( p(i,j,kp)-p(i,j,k) )*dzci(k) )
            !
          enddo
        enddo
      enddo
    case(2)
      do k=1+IBM_h,nz
        do j=1,ny
          do i=1,nx
            !
            ip = i+1
            jp = j+1
            kp = k+1
            !
            rhox    = 0.5_rp*(rho(ip,j,k)+rho(i,j,k))
            rhoxi   = 1.0_rp / rhox
            rhoy    = 0.5_rp*(rho(i,jp,k)+rho(i,j,k))
            rhoyi   = 1.0_rp / rhoy
            rhoz    = 0.5_rp*(rho(i,j,kp)+rho(i,j,k))
            rhozi   = 1.0_rp / rhoz
            !
            u(i,j,k) = u(i,j,k) + f_t12*(- (rhoxi  - rho0i)*f2*( dp(ip,j,k)-dp(i,j,k) )*dxi)
            v(i,j,k) = v(i,j,k) + f_t12*(- (rhoyi  - rho0i)*f2*( dp(i,jp,k)-dp(i,j,k) )*dyi)
            w(i,j,k) = w(i,j,k) + f_t12*(- (rhozi  - rho0i)*f2*( dp(i,j,kp)-dp(i,j,k) )*dzci(k))
            !
          enddo
        enddo
      enddo
    case(3)
      do k=1+IBM_h,nz
        do j=1,ny
          do i=1,nx
            !
            ip = i+1
            jp = j+1
            kp = k+1
            !
            rhox    = 0.5_rp*(rho(ip,j,k)+rho(i,j,k))
            rhoxi   = 1.0_rp / rhox
            rhoy    = 0.5_rp*(rho(i,jp,k)+rho(i,j,k))
            rhoyi   = 1.0_rp / rhoy
            rhoz    = 0.5_rp*(rho(i,j,kp)+rho(i,j,k))
            rhozi   = 1.0_rp / rhoz
            !
            u(i,j,k) = u(i,j,k) + f_t12*(- (rhoxi  - rho0i)*( dp(ip,j,k)-dp(i,j,k) )*dxi)
            v(i,j,k) = v(i,j,k) + f_t12*(- (rhoyi  - rho0i)*( dp(i,jp,k)-dp(i,j,k) )*dyi)
            w(i,j,k) = w(i,j,k) + f_t12*(- (rhozi  - rho0i)*( dp(i,j,kp)-dp(i,j,k) )*dzci(k))
            !
          enddo
        enddo
      enddo
    case default
      call flutas_error('ERROR in pres_tw_src: Invalid pressure stage specified!')
  end select
  !$acc end kernels
  !
  return
end subroutine pres_tw_src
!
!_______________________________________________________________________________
!
subroutine calculate_new_pressure(nx,ny,nz,dxi,dyi,dzi,nh_d,nh_u,nh_p,nh_v,dzci,dzfi,u,v,w,p,dp,rho &
     ,src_ph,f_t12,dl,dzc,dzf,halo_u,halo_p,f_t12_o,istep,rhsbx,rhsby,rhsbz, &
     arrplanp,normfftp,lambdaxyp,ap,bp,cp,l,trigger_update )  
  ! solves the pressure with two different methods, only CPU version available right now
  use mod_param     ,only: rho1,rho2,IBM_h,no_outflow,bcvel,cbcvel,is_outflow,cbcpre,bcpre,rho0
  use mod_bound     ,only: bounduvw,updt_rhs_b,boundp,bounduvw_mom,boundp_IBM
  use mod_fillps    ,only: fillps
  use mod_solver_cpu,only: solver_cpu
  use mod_correc    ,only: correc
  ! 
  implicit none
  !
  integer    , intent(in   )                                     :: nx,ny,nz,istep
  integer    , intent(in   ), dimension(3)                       :: halo_u,halo_p
  real(rp)   , intent(in   ), dimension(3)                       :: dl,l
  integer    , intent(in   )                                     :: nh_d,nh_u,nh_p,nh_v
  real(rp)   , intent(in   ), dimension(1-nh_d:)                 :: dzc,dzf
  real(rp   ), intent(in   )                                     :: f_t12,f_t12_o
  real(rp)   , intent(in   )                                     :: dxi,dyi,dzi
  real(rp)   , intent(in   ), dimension(1-nh_d:                ) :: dzci,dzfi
  real(rp)   , intent(inout), dimension(1-nh_u:,1-nh_u:,1-nh_u:) :: u,v,w
  real(rp)   , intent(in   ), dimension(1-nh_p:,1-nh_p:,1-nh_p:) :: rho,src_ph
  real(rp)   , intent(inout), dimension(1-nh_p:,1-nh_p:,1-nh_p:) :: p,dp
  real(rp)   , intent(in   ), dimension(      :,      :,     0:) :: rhsbx,rhsby,rhsbz
  type(C_PTR), intent(in   ), dimension(2,2)                     :: arrplanp
  real(rp)   , intent(in   )                                     :: normfftp
  real(rp)   , intent(in   ), dimension(n_z(1),n_z(2))           :: lambdaxyp
  real(rp)   , intent(in   ), dimension(n_z(3)-IBM_h)            :: ap,bp,cp
  logical    , intent(in   )                                     :: trigger_update
  real(rp)   , allocatable  , dimension(:,:,:)                   :: vel_1_x,vel_1_y,vel_1_z
  real(rp)   , allocatable  , dimension(:,:,:)                   :: vel_2_x,vel_2_y,vel_2_z
  real(rp)                                                       :: nse_switch
  real(rp)                                                       :: rhox,rhoy,rhoz,dti,limiter
  integer                                                        :: i,j,k
  real(rp)   , allocatable  , dimension(:,:,:)                   :: dp1,dp2,dp3,u_grad_rho
  real(rp)                                                       :: rho0i,f_ex
  real(rp) :: dtidxi,dtidyi,rxm,rxp,rym,ryp,rzm,rzp !!!
  !
  nse_switch = 0.0
  if(use_mod_nse)  nse_switch = 1.0
  dti = 1._rp/f_t12
  rho0i = 1._rp/rho0
  f_ex = f_t12/f_t12_o
  !
  ! u,v,w are the momenta here
  !
  allocate( vel_1_x(1-nh_u:nx+nh_u,1-nh_u:ny+nh_u,1-nh_u:nz+nh_u) , &
            vel_1_y(1-nh_u:nx+nh_u,1-nh_u:ny+nh_u,1-nh_u:nz+nh_u) , &
            vel_1_z(1-nh_u:nx+nh_u,1-nh_u:ny+nh_u,1-nh_u:nz+nh_u))
  !extract the first velocity
  call mom_to_vel((/nx,ny,nz/),nh_d,nh_v,nh_u,u,v,w,vel_1_x,vel_1_y,vel_1_z,rho)
  call collect_vel((/nx,ny,nz/),dl,l,nh_p,nh_u,rho,vel_1_x,vel_1_y,vel_1_z,p)
  !
  if(use_mod_solver) then
    allocate( vel_2_x(1-nh_u:nx+nh_u,1-nh_u:ny+nh_u,1-nh_u:nz+nh_u) , &
              vel_2_y(1-nh_u:nx+nh_u,1-nh_u:ny+nh_u,1-nh_u:nz+nh_u) , &
              vel_2_z(1-nh_u:nx+nh_u,1-nh_u:ny+nh_u,1-nh_u:nz+nh_u)) 
    call copy_fields((/nx,ny,nz/),nh_u,vel_1_x,vel_2_x) 
    call copy_fields((/nx,ny,nz/),nh_u,vel_1_y,vel_2_y) 
    call copy_fields((/nx,ny,nz/),nh_u,vel_1_z,vel_2_z)
  endif
  !
  !Add extrapolated pressure increment with density splitting to first velocity
  !if(.NOT. trigger_update) then
    call pres_tw_src(nx,ny,nz,dxi,dyi,dzi,nh_d,nh_u,dzci,rho0i,f_t12,f_t12_o,p,dp,rho,&
                                                    vel_1_x,vel_1_y,vel_1_z,2,istep)
  !endif
  !Communicate before constructing the RHS
  call bounduvw(cbcvel,(/nx,ny,nz/),bcvel,nh_d,nh_u,halo_u,no_outflow,dl,dzc,dzf, &
                                     vel_1_x,vel_1_y,vel_1_z,rho,rex,rey,rez,0.0_rp,1.0_rp,.TRUE.)
  call apply_bound_vel((/nx,ny,nz/),dl,l,nh_p,nh_u,vel_1_x,vel_1_y,vel_1_z,halo_u,halo_p,p,.FALSE.,rho)
  ! Allocate first temporaray pressure increment array dp1
  allocate( dp1(1-nh_p:nx+nh_p,1-nh_p:ny+nh_p,1-nh_p:nz+nh_p))
  ! Fill dp1 with the RHS information 
  call fillps(nx,ny,nz,nh_d,nh_u,dxi,dyi,dzi,dzfi,dti,rho0,rho1,rho2,rho,rho, &   !2nd rho is here as dummy
                     src_ph,vel_1_x,vel_1_y,vel_1_z,dp1,.FALSE.,rex,rey,rez,nse_switch,out_rk_safety)
  !
  !-------------- First FFT -------------------
  call updt_rhs_b(nx,ny,nz,(/'c','c','c'/),cbcpre,nh_p,rhsbx,rhsby,rhsbz,dp1)
  call solver_cpu((/nx,ny,nz/),arrplanp,normfftp,lambdaxyp,ap,bp,cp,cbcpre(:,3),(/'c','c','c'/),dp1)
  call boundp(cbcpre,(/nx,ny,nz/),bcpre,nh_d,nh_p,halo_p,dl,dzc,dzf,dp1)
  !--------------------------------------------
  ! Correct the velocities
  call correc(nx,ny,nz,nh_d,nh_u,dxi,dyi,dzi,dzci,f_t12,rho0,dp1,vel_1_x,vel_1_y,vel_1_z,rho)
  ! If neccesay turn velocity back into a volume avergaed mixture velocity
  if(use_mod_nse) call make_vol_avg_vel((/nx,ny,nz/),nh_d,nh_v,nh_u,vel_1_x,vel_1_y,vel_1_z,rho)
  !Do IBM treatment and communicate velocities
  call bounduvw(cbcvel,(/nx,ny,nz/),bcvel,nh_d,nh_u,halo_u,is_outflow,dl,dzc,dzf, &
                                       vel_1_x,vel_1_y,vel_1_z,rho,rex,rey,rez,0.0_rp,0.0_rp,.TRUE.)
  call apply_bound_vel((/nx,ny,nz/),dl,l,nh_p,nh_u,vel_1_x,vel_1_y,vel_1_z,halo_u,halo_p,p,.FALSE.,rho)    

  !From now on two solvers are possible:
  if(.NOT.use_mod_solver) then 
    call set_IBM_vel_to_zero((/nx,ny,nz/),nh_u, vel_1_x,vel_1_y,vel_1_z)
    ! Copy velocity to output array and deallocate auxiluary velocity arrays
    call copy_fields((/nx,ny,nz/),nh_u,vel_1_x,u) 
    call copy_fields((/nx,ny,nz/),nh_u,vel_1_y,v) 
    call copy_fields((/nx,ny,nz/),nh_u,vel_1_z,w)
    deallocate(vel_1_x,vel_1_y,vel_1_z)
    ! update the pressure increment
    call copy_fields((/nx,ny,nz/),nh_p,dp1,dp)
    ! update the pressure
    call add_a_to_b((/nx,ny,nz/),nh_p,dp1,p) 
    deallocate(dp1)
    call boundp(cbcpre,(/nx,ny,nz/),bcpre,nh_d,nh_p,halo_p,dl,dzc,dzf,p)
    !end of first method -> this is the same as in Frantzis et al (2019)
  else
    !This is the new solver
    !calculate u.grad(rho)
    allocate(u_grad_rho(1:nx,1:ny,1:nz))
    call cmpt_u_grad_rho((/nx,ny,nz/),nh_p,nh_u,vel_1_x,vel_1_y,vel_1_z,rho,u_grad_rho,dxi,dyi,dzi)
    !
    ! velocities from first step are not needed anymore
    deallocate(vel_1_x,vel_1_y,vel_1_z)
    !Communicate before constructing the RHS
    call bounduvw_mom(cbcvel,(/nx,ny,nz/),bcvel,nh_d,nh_u,halo_u,is_outflow,dl,dzc,dzf, &
                                                   u,v,w,rho,rex,rey,rez,0.0_rp,1.0_rp,.TRUE.)
    call apply_bound_vel((/nx,ny,nz/),dl,l,nh_p,nh_u,u,v,w,halo_u,halo_p,p,.TRUE.,rho)
    !allocate dp2 and fill it with RHS information
    allocate( dp2(1-nh_p:nx+nh_p,1-nh_p:ny+nh_p,1-nh_p:nz+nh_p))  
    call fillps(nx,ny,nz,nh_d,nh_u,dxi,dyi,dzi,dzfi,dti,rho0,rho1,rho2,rho,u_grad_rho, & 
               src_ph,u,v,w,dp2,.TRUE.,rex,rey,rez,nse_switch,out_rk_safety)
    ! the array of u.grad(rho) is not needed anymore
    deallocate(u_grad_rho)
    !-------------- Second FFT ------------------
    call updt_rhs_b(nx,ny,nz,(/'c','c','c'/),cbcpre,nh_p,rhsbx,rhsby,rhsbz,dp2)
    call solver_cpu((/nx,ny,nz/),arrplanp,normfftp,lambdaxyp,ap,bp,cp,cbcpre(:,3),(/'c','c','c'/),dp2)
    call boundp(cbcpre,(/nx,ny,nz/),bcpre,nh_d,nh_p,halo_p,dl,dzc,dzf,dp2)
    !--------------------------------------------
    ! calculate the limiter based on dp1 and dp2
    limiter = compute_pres_limiter_v1((/nx,ny,nz/),nh_p,dp1,dp2,dl)
    !p_lim_array(mod(istep,size(p_lim_array))+1) = limiter
    !limiter = sum(p_lim_array)/size(p_lim_array)
    if(myid.eq.0) print*, 'limiter = ',limiter
    ! dp1 is not needed anymore
    deallocate(dp1)
    ! adjust dp2 using 'limiter' as the underrelaxation parameter
    call mult_by_const((/nx,ny,nz/),nh_p,dp2,limiter)
    !Add predicted pressure increment (dp2) with density splitting to second velocity
    call pres_tw_src(nx,ny,nz,dxi,dyi,dzi,nh_d,nh_u,dzci,rho0i,f_t12,f_t12_o,p,dp2,rho,&
                                         vel_2_x,vel_2_y,vel_2_z,3,istep)
    ! dp2 is not needed anymore
    deallocate(dp2)
    !Communicate before constructing the RHS
    call bounduvw(cbcvel,(/nx,ny,nz/),bcvel,nh_d,nh_u,halo_u,no_outflow,dl,dzc,dzf, &
                          vel_2_x,vel_2_y,vel_2_z,rho,rex,rey,rez,0.0_rp,1.0_rp,.TRUE.)
    call apply_bound_vel((/nx,ny,nz/),dl,l,nh_p,nh_u,vel_2_x,vel_2_y,vel_2_z,halo_u,halo_p,p,.FALSE.,rho)
    !allocate dp3 and fill it with RHS information
    allocate( dp3(1-nh_p:nx+nh_p,1-nh_p:ny+nh_p,1-nh_p:nz+nh_p))  
    call fillps(nx,ny,nz,nh_d,nh_u,dxi,dyi,dzi,dzfi,dti,rho0,rho1,rho2,rho,rho, &    !2nd rho is here as dummy
                    src_ph,vel_2_x,vel_2_y,vel_2_z,dp3,.FALSE.,rex,rey,rez,nse_switch,out_rk_safety)
    !-------------- Third FFT -------------------
    call updt_rhs_b(nx,ny,nz,(/'c','c','c'/),cbcpre,nh_p,rhsbx,rhsby,rhsbz,dp3)
    call solver_cpu((/nx,ny,nz/),arrplanp,normfftp,lambdaxyp,ap,bp,cp,cbcpre(:,3),(/'c','c','c'/),dp3)
    call boundp(cbcpre,(/nx,ny,nz/),bcpre,nh_d,nh_p,halo_p,dl,dzc,dzf,dp3)
    !--------------------------------------------
    ! Correct the velocities
    call correc(nx,ny,nz,nh_d,nh_u,dxi,dyi,dzi,dzci,f_t12,rho0,dp3,vel_2_x,vel_2_y,vel_2_z,rho)
    ! If neccesay turn velocity back into a volume avergaed mixture velocity
    if(use_mod_nse) call make_vol_avg_vel((/nx,ny,nz/),nh_d,nh_v,nh_u,vel_2_x,vel_2_y,vel_2_z,rho)
    !Do IBM treatment and communicate velocities
    call set_IBM_vel_to_zero((/nx,ny,nz/),nh_u,vel_2_x,vel_2_y,vel_2_z)
    call bounduvw(cbcvel,(/nx,ny,nz/),bcvel,nh_d,nh_u,halo_u,is_outflow,dl,dzc,dzf, &
                                  vel_2_x,vel_2_y,vel_2_z,rho,rex,rey,rez,0.0_rp,0.0_rp,.TRUE.)
    ! Copy velocity to output array and deallocate auxiluary velocity arrays
    call copy_fields((/nx,ny,nz/),nh_u,vel_2_x,u) 
    call copy_fields((/nx,ny,nz/),nh_u,vel_2_y,v) 
    call copy_fields((/nx,ny,nz/),nh_u,vel_2_z,w)
    deallocate(vel_2_x,vel_2_y,vel_2_z)
    ! update the pressure increment
    call copy_fields((/nx,ny,nz/),nh_p,dp3,dp)
    ! update the pressure
    call add_a_to_b((/nx,ny,nz/),nh_p,dp3,p) 
    deallocate(dp3)
    call boundp(cbcpre,(/nx,ny,nz/),bcpre,nh_d,nh_p,halo_p,dl,dzc,dzf,p)
    call apply_bound_vel((/nx,ny,nz/),dl,l,nh_p,nh_u,u,v,w,halo_u,halo_p,p,.FALSE.,rho)
  endif
  !
  return
end subroutine calculate_new_pressure
!
!_______________________________________________________________________________
!
real(rp) function compute_pres_limiter_v1(n,nh_p,p1,p2,dl)
  ! returns the limiter based on the two first prediction pressures
  use mod_param,      only: is_outflow
  !   
  ! imposes Dirichlet BC at the outflow
  !   
  implicit none
  !   
  integer , intent(in   ), dimension(3)                       :: n
  integer , intent(in   )                                     :: nh_p
  real(rp), intent(in   ), dimension(3)                       :: dl
  real(rp), intent(in   ), dimension(1-nh_p:,1-nh_p:,1-nh_p:) :: p1,p2
  integer                                                     :: i,j,k
  real(rp)                                                    :: relax
  real(rp)                                                    :: maxp1,maxp2
  real(rp)                                                    :: minp1,minp2
  ! To avoid to divide by zero
  real(rp),parameter                                          :: small_num = 1e-12
  !
  maxp1 = 0.0_rp
  maxp2 = 0.0_rp
  minp1 = 0.0_rp
  minp2 = 0.0_rp
  !
  relax = 1.0_rp ! hard coded relaxation factor (0.9<relax<1.0)
  !
  maxp1=maxval(p1(1:n(1),1:n(2),1:n(3)))
  call mpi_allreduce(MPI_IN_PLACE,maxp1,1,MPI_REAL_RP,MPI_MAX,MPI_COMM_WORLD,ierr)
  minp1=minval(p1(1:n(1),1:n(2),1:n(3)))
  call mpi_allreduce(MPI_IN_PLACE,minp1,1,MPI_REAL_RP,MPI_MIN,MPI_COMM_WORLD,ierr)
  maxp2=maxval(p2(1:n(1),1:n(2),1:n(3)))
  call mpi_allreduce(MPI_IN_PLACE,maxp2,1,MPI_REAL_RP,MPI_MAX,MPI_COMM_WORLD,ierr)
  minp2=minval(p2(1:n(1),1:n(2),1:n(3)))
  call mpi_allreduce(MPI_IN_PLACE,minp2,1,MPI_REAL_RP,MPI_MIN,MPI_COMM_WORLD,ierr)
  if(maxp2.EQ.0.0) maxp2=small_num
  if(minp2.EQ.0.0) minp2=small_num
  !
  compute_pres_limiter_v1 = relax*min(1.0,abs(maxp1)/abs(maxp2),abs(minp1)/abs(minp2))
  !
end function compute_pres_limiter_v1
!
!_______________________________________________________________________________
!
real(rp) function compute_pres_limiter_v2(n,nh_p,p1,p2,dl)
  ! returns the limiter based on the two first prediction pressures
  use mod_param,      only: is_outflow
  !   
  ! imposes Dirichlet BC at the outflow
  !   
  implicit none
  !   
  integer , intent(in   ), dimension(3)                       :: n
  integer , intent(in   )                                     :: nh_p
  real(rp), intent(in   ), dimension(3)                       :: dl
  real(rp), intent(in   ), dimension(1-nh_p:,1-nh_p:,1-nh_p:) :: p1,p2
  integer                                                     :: i,j,k
  real(rp)                                                    :: int1,int2,relax
  ! To avoid to divide by zero
  real(rp),parameter                                          :: small_num = 1e-12
  !
  int1 = 0.0
  int2 = 0.0
  !
  relax = 1.0 ! hard coded relaxation factor (0.9<relax<1.0)
  !
  do k=1,n(3)
    do j=1,n(2)
      do i=1,n(1)
        int1=int1+min(1.0,max(0.0,p1(i,j,k)/p2(i,j,k)))*abs(p2(i,j,k)+small_num)
        int2=int2+abs(p2(i,j,k)+small_num);
      enddo
    enddo
  enddo
  !
  call mpi_allreduce(MPI_IN_PLACE,int1,1,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
  call mpi_allreduce(MPI_IN_PLACE,int2,1,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
  !
  compute_pres_limiter_v2 = relax*int1/int2
  !
end function compute_pres_limiter_v2
!
!_______________________________________________________________________________
!
subroutine cmpt_u_grad_rho(n,nh_p,nh_u,u,v,w,rho,u_grad_rho,dxi,dyi,dzi)
  ! compute v.grad(rho)
  use mod_param     , only: IBM_h,rho1,rho2
  implicit none
  !
  integer , intent(in   ), dimension(3)                       :: n
  integer , intent(in   )                                     :: nh_p,nh_u
  real(rp), intent(in   ), dimension(1-nh_u:,1-nh_u:,1-nh_u:) :: u,v,w
  real(rp), intent(in   ), dimension(1-nh_p:,1-nh_p:,1-nh_p:) :: rho
  real(rp), intent(out  ), dimension(1:     ,1:     ,1:     ) :: u_grad_rho
  real(rp), intent(in   )                                     :: dxi,dyi,dzi
  integer                                                     :: i,j,k
  !
  u_grad_rho(:,:,:) = 0.0
  !
  do k=1+IBM_h,n(3)
    do j=1,n(2)
      do i=1,n(1)
        u_grad_rho(i,j,k) = 0.5*(u(i  ,j  ,k  )*(rho(i+1,j  ,k  )-rho(i  ,j  ,k  ))*dxi       &
                                +u(i-1,j  ,k  )*(rho(i  ,j  ,k  )-rho(i-1,j  ,k  ))*dxi)      &
                           +0.5*(v(i  ,j  ,k  )*(rho(i  ,j+1,k  )-rho(i  ,j  ,k  ))*dyi       &
                                +v(i  ,j-1,k  )*(rho(i  ,j  ,k  )-rho(i  ,j-1,k  ))*dyi)      &
                           +0.5*(w(i  ,j  ,k  )*(rho(i  ,j  ,k+1)-rho(i  ,j  ,k  ))*dzi       &
                                +w(i  ,j  ,k-1)*(rho(i  ,j  ,k  )-rho(i  ,j  ,k-1))*dzi)
      enddo
    enddo
  enddo
  !
  return
end subroutine cmpt_u_grad_rho

!_______________________________________________________________________________

subroutine mom_to_vel(n,nh_d,nh_v,nh_u,u_mom,v_mom,w_mom,u_vel,v_vel,w_vel,rho)
  ! turn a momentum (input) into a velocity (output)
  use mod_param     , only: rho1,rho2
  implicit none
  !
  integer , intent(in   ), dimension(3)                       :: n
  integer , intent(in   )                                     :: nh_d,nh_v,nh_u
  real(rp), intent(in   ), dimension(1-nh_u:,1-nh_u:,1-nh_u:) :: u_mom,v_mom,w_mom
  real(rp), intent(out  ), dimension(1-nh_u:,1-nh_u:,1-nh_u:) :: u_vel,v_vel,w_vel
  real(rp), intent(in   ), dimension(     0:,     0:,     0:) :: rho
  integer                                                     :: i,j,k
  real(rp)                                                    :: rhox,rhoy,rhoz
  !
  do k=1,n(3)
    do j=1,n(2)
      do i=1,n(1)
        !
        rhox = 0.5d0*(rho(i,j,k)+rho(i+1,j  ,k  ))
        rhoy = 0.5d0*(rho(i,j,k)+rho(i  ,j+1,k  ))
        rhoz = 0.5d0*(rho(i,j,k)+rho(i  ,j  ,k+1))
        !
        u_vel(i,j,k) = u_mom(i,j,k)/rhox                                     
        v_vel(i,j,k) = v_mom(i,j,k)/rhoy
        w_vel(i,j,k) = w_mom(i,j,k)/rhoz
      enddo
    enddo
  enddo
  !
  return
end subroutine mom_to_vel
!
!#################################################
!###     End of:  Pressure subrouties          ###
!#################################################
!
!#################################################
!###   Start of:  regularisation subrouties    ###
!#################################################
!
!_______________________________________________________________________________

subroutine correct_mom(nx,ny,nz,dxi,dyi,dzi,nh_d,nh_u,dzci,dzfi,u,v,w,mu,rho, &
                              dudt,dvdt,dwdt,f_t12,psi,dl,dzc,dzf,halo_u,f_t12_o)  
    !
    use mod_param    ,only: rho1,rho2,IBM_h,no_outflow,is_outflow
    use mod_bound,    only: bounduvw
    use mod_common_mpi, only: MPI_PROC_NULL,left,right,front,back,top,bottom
    !
    implicit none
    !
    integer , intent(in   )                                     :: nx,ny,nz
    integer , intent(in   ), dimension(3)                       :: halo_u
    real(rp), intent(in   ), dimension(3)                       :: dl
    integer , intent(in   )                                     :: nh_d,nh_u    
    real(rp), intent(in   ), dimension(1-nh_d:)                 :: dzc,dzf
    real(rp), intent(in   )                                     :: f_t12,f_t12_o
    real(rp), intent(in   )                                     :: dxi,dyi,dzi
    real(rp), intent(in   ), dimension(1-nh_d:                ) :: dzci,dzfi
    real(rp), intent(in   ), dimension(1-nh_u:,1-nh_u:,1-nh_u:) :: u,v,w
    real(rp), intent(in   ), dimension(     0:,     0:,     0:) :: mu,rho
    real(rp), intent(in   ), dimension(     0:,     0:,     0:) :: psi
    real(rp), intent(inout), dimension(     1:,     1:,     1:) :: dudt,dvdt,dwdt
    integer                                                     :: length,x_start,x_stop,y_start
    integer                                                     :: y_stop,z_start,z_stop,min_dist
    !
    real(rp) :: dudxp,dudxm,dudyp,dudym,dudzp,dudzm, &
                dvdxp,dvdxm,dvdyp,dvdym,dvdzp,dvdzm, &
                dwdxp,dwdxm,dwdyp,dwdym,dwdzp,dwdzm
    real(rp) :: muxm,muxp,muym,muyp,muzm,muzp, &
                rhoxm,rhoxp,rhoyp,rhoym,rhozp,rhozm,phix,phiy,phiz
    real(rp) :: rhox,rhoy,rhoz
    integer  :: im,ip,jm,jp,km,kp,i,j,k
    real(rp) :: uuip,uuim,uvjp,uvjm,uwkp,uwkm, &
                uvip,uvim,vvjp,vvjm,wvkp,wvkm, &
                uwip,uwim,vwjp,vwjm,wwkp,wwkm
    real(rp) :: term1,term2,term3,overall_corr
    real(rp), parameter :: small_num = 1e-10
    real(rp), allocatable  , dimension(:,:,:,:)                   :: auxfield
    x_start = 0
    x_stop  = 0
    y_start = 0
    y_stop  = 0
    z_start = 0
    z_stop  = 0   
    if    (is_outflow(0,1).AND.(left .eq.MPI_PROC_NULL)) then
      x_start = out_rk_safety
    elseif(is_outflow(1,1).AND.(right.eq.MPI_PROC_NULL)) then
      x_stop = out_rk_safety
    elseif(is_outflow(0,2).AND.(front.eq.MPI_PROC_NULL)) then
      y_start = out_rk_safety
    elseif(is_outflow(1,2).AND.(back .eq.MPI_PROC_NULL)) then
      y_stop = out_rk_safety
    elseif(is_outflow(0,3)) then
      z_start = out_rk_safety
    elseif(is_outflow(1,3)) then
      z_stop = out_rk_safety
    endif
    !
    !@cuf attributes(managed) :: u, v, w, mu, rho, dzci, dzfi, dudt, dvdt, dwdt
    !    
    !$acc kernels
    allocate(auxfield(0:nx+1,0:ny+1,0:nz+1,3))
    auxfield(:,:,:,:) = 0.0
    !
    do k=IBM_h+1,nz
      do j=1,ny
        do i=1,nx
          !
          ip = i+1   
          jp = j+1   
          kp = k+1  
          !
          phix   = 0.5_rp*(psi(ip,j,k)+psi(i,j,k))
          rhox   = phix*rho1+(1-phix)*rho2
          phiy   = 0.5_rp*(psi(i,jp,k)+psi(i,j,k))
          rhoy   = phiy*rho1+(1-phiy)*rho2
          phiz   = 0.5_rp*(psi(i,j,kp)+psi(i,j,k))
          rhoz   = phiz*rho1+(1-phiz)*rho2
#if !defined(_TWOD)
          if((phix.GT.small_num).AND.((1.0-phix).GT.small_num))  auxfield(i,j,k,1) = &
                                   (rho1/(phix) + rho2/(1.0-phix))*rex(i,j,k)
#endif
          if((phiy.GT.small_num).AND.((1.0-phiy).GT.small_num))  auxfield(i,j,k,2) = &
                                   (rho1/(phiy) + rho2/(1.0-phiy))*rey(i,j,k)
          if((phiz.GT.small_num).AND.((1.0-phiz).GT.small_num))  auxfield(i,j,k,3) = &
                                   (rho1/(phiz) + rho2/(1.0-phiz))*rez(i,j,k)
        enddo
      enddo
    enddo
    !
    call bounduvw(cbcreg,(/nx,ny,nz/),bcreg,nh_d,nh_u,halo_u,no_outflow,dl,dzc &
               ,dzf,auxfield(:,:,:,1),auxfield(:,:,:,2),auxfield(:,:,:,3), &
               rho,rex,rey,rez,0.0_rp,0.0_rp,.FALSE.)
    !
    !
    do k=IBM_h+1+z_start,nz-z_stop
      do j=1+y_start,ny-y_stop
        do i=x_start+1,nx-x_stop
          !
          ip = i + 1
          im = i - 1
          jp = j + 1
          jm = j - 1
          kp = k + 1
          km = k - 1
          !
          !
          !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          !
          ! 3. r.h.s. of x-momentum
          !
          !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          term2=0.0
          term3=0.0
          overall_corr = 0.0
          !
          uuip  = 0.25_rp*( rex(ip,j,k)+rex(i,j,k) )*( u(i+1,j  ,k  )+u(i,j  ,k  ) )  
          uuim  = 0.25_rp*( rex(im,j,k)+rex(i,j,k) )*( u(i-1,j  ,k  )+u(i,j  ,k  ) )  
          uvjp  = 0.25_rp*( rex(i,jp,k)+rex(i,j,k) )*( v(i+1,j  ,k  )+v(i,j  ,k  ) )  
          uvjm  = 0.25_rp*( rex(i,jm,k)+rex(i,j,k) )*( v(i+1,j-1,k  )+v(i,j-1,k  ) )  
          uwkp  = 0.25_rp*( rex(i,j,kp)+rex(i,j,k) )*( w(i+1,j  ,k  )+w(i,j  ,k  ) )  
          uwkm  = 0.25_rp*( rex(i,j,km)+rex(i,j,k) )*( w(i+1,j  ,k-1)+w(i,j  ,k-1) )
          !
          term2 =  (rho1-rho2)*(dxi*(uuip - uuim ) +  dyi*(uvjp - uvjm ) + dzci(k)*(uwkp - uwkm ))
          !
          uuip  = 0.25_rp*( rex(ip,j,k)+rex(i,j,k) )*( auxfield(i+1,j  ,k  ,1)+auxfield(i,j  ,k  ,1) )  
          uuim  = 0.25_rp*( rex(im,j,k)+rex(i,j,k) )*( auxfield(i-1,j  ,k  ,1)+auxfield(i,j  ,k  ,1) )  
          uvjp  = 0.25_rp*( rex(i,jp,k)+rex(i,j,k) )*( auxfield(i+1,j  ,k  ,2)+auxfield(i,j  ,k  ,2) )  
          uvjm  = 0.25_rp*( rex(i,jm,k)+rex(i,j,k) )*( auxfield(i+1,j-1,k  ,2)+auxfield(i,j-1,k  ,2) )  
          uwkp  = 0.25_rp*( rex(i,j,kp)+rex(i,j,k) )*( auxfield(i+1,j  ,k  ,3)+auxfield(i,j  ,k  ,3) )  
          uwkm  = 0.25_rp*( rex(i,j,km)+rex(i,j,k) )*( auxfield(i+1,j  ,k-1,3)+auxfield(i,j  ,k-1,3) )  
          !
          term3 =  (dxi*(uuip - uuim ) +  dyi*(uvjp - uvjm ) + dzci(k)*(uwkp - uwkm ))
          !
          !
#if defined(_TWOD)
          overall_corr = 0.0
#else
          overall_corr = - term3   +term2
#endif
          !
          dudt(i,j,k) = dudt(i,j,k) + overall_corr
          !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          !
          ! 3. r.h.s. of y-momentum
          !
          !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          term2=0.0
          term3=0.0
          overall_corr = 0.0
          !
          uvip  = 0.25*( u(i  ,j,k  )+u(i  ,j+1,k  ) )*( rey(i,j,k  )+rey(i+1,j  ,k) )
          uvim  = 0.25*( u(i-1,j,k  )+u(i-1,j+1,k  ) )*( rey(i,j,k  )+rey(i-1,j  ,k) )
          vvjp  = 0.25*( v(i  ,j,k  )+v(i  ,j+1,k  ) )*( rey(i,j,k  )+rey(i  ,j+1,k) )
          vvjm  = 0.25*( v(i  ,j,k  )+v(i  ,j-1,k  ) )*( rey(i,j,k  )+rey(i  ,j-1,k) )
          wvkp  = 0.25*( w(i  ,j,k  )+w(i  ,j+1,k  ) )*( rey(i,j,k+1)+rey(i  ,j  ,k) )
          wvkm  = 0.25*( w(i  ,j,k-1)+w(i  ,j+1,k-1) )*( rey(i,j,k-1)+rey(i  ,j  ,k) )
          !
          term2 =  (rho1-rho2)*(dxi*(uvip - uvim ) +  dyi*(vvjp - vvjm ) + dzci(k)*(wvkp - wvkm ))
          !
          uvip  = 0.25*( auxfield(i  ,j,k  ,1)+auxfield(i  ,j+1,k  ,1) )*( rey(i,j,k)+rey(i+1,j,k) )
          uvim  = 0.25*( auxfield(i-1,j,k  ,1)+auxfield(i-1,j+1,k  ,1) )*( rey(i,j,k)+rey(i-1,j,k) )
          vvjp  = 0.25*( auxfield(i  ,j,k  ,2)+auxfield(i  ,j+1,k  ,2) )*( rey(i,j,k)+rey(i,j+1,k) )
          vvjm  = 0.25*( auxfield(i  ,j,k  ,2)+auxfield(i  ,j-1,k  ,2) )*( rey(i,j,k)+rey(i,j-1,k) )
          wvkp  = 0.25*( auxfield(i  ,j,k  ,3)+auxfield(i  ,j+1,k  ,3) )*( rey(i,j,k+1)+rey(i,j,k) )
          wvkm  = 0.25*( auxfield(i  ,j,k-1,3)+auxfield(i  ,j+1,k-1,3) )*( rey(i,j,k-1)+rey(i,j,k) )
          !
          term3 =  (dxi*(uvip - uvim ) +  dyi*(vvjp - vvjm ) + dzci(k)*(wvkp - wvkm ))
          !
          overall_corr =  - term3  + term2
          !
          dvdt(i,j,k) = dvdt(i,j,k) + overall_corr
          !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          !
          ! 3. r.h.s. of z-momentum
          !
          !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          term2=0.0
          term3=0.0
          overall_corr = 0.0
          !
          uwip  = 0.25*( rez(i,j,k)+rez(i+1,j,k) )*( u(i  ,j  ,k)+u(i  ,j  ,k+1) )
          uwim  = 0.25*( rez(i,j,k)+rez(i-1,j,k) )*( u(i-1,j  ,k)+u(i-1,j  ,k+1) )
          vwjp  = 0.25*( rez(i,j,k)+rez(i,j+1,k) )*( v(i  ,j  ,k)+v(i  ,j  ,k+1) )
          vwjm  = 0.25*( rez(i,j,k)+rez(i,j-1,k) )*( v(i  ,j-1,k)+v(i  ,j-1,k+1) )
          wwkp  = 0.25*( rez(i,j,k)+rez(i,j,k+1) )*( w(i  ,j  ,k)+w(i  ,j  ,k+1) )
          wwkm  = 0.25*( rez(i,j,k)+rez(i,j,k-1) )*( w(i  ,j  ,k)+w(i  ,j  ,k-1) )
          !
          term2 =  (rho1-rho2)*(dxi*(uwip - uwim ) +  dyi*(vwjp - vwjm ) + dzci(k)*(wwkp - wwkm ))
          !
          uwip  = 0.25*( rez(i,j,k)+rez(i+1,j,k) )*( auxfield(i  ,j  ,k,1)+auxfield(i  ,j  ,k+1,1) )
          uwim  = 0.25*( rez(i,j,k)+rez(i-1,j,k) )*( auxfield(i-1,j  ,k,1)+auxfield(i-1,j  ,k+1,1) )
          vwjp  = 0.25*( rez(i,j,k)+rez(i,j+1,k) )*( auxfield(i  ,j  ,k,2)+auxfield(i  ,j  ,k+1,2) )
          vwjm  = 0.25*( rez(i,j,k)+rez(i,j-1,k) )*( auxfield(i  ,j-1,k,2)+auxfield(i  ,j-1,k+1,2) )
          wwkp  = 0.25*( rez(i,j,k)+rez(i,j,k+1) )*( auxfield(i  ,j  ,k,3)+auxfield(i  ,j  ,k+1,3) )
          wwkm  = 0.25*( rez(i,j,k)+rez(i,j,k-1) )*( auxfield(i  ,j  ,k,3)+auxfield(i  ,j  ,k-1,3) )
          !
          term3 =  (dxi*(uwip - uwim ) +  dyi*(vwjp - vwjm ) + dzci(k)*(wwkp - wwkm ))
          !
          overall_corr =  - term3  + term2
          !
          dwdt(i,j,k) = dwdt(i,j,k) + overall_corr
          !
        enddo
      enddo
    enddo
    !
    deallocate(auxfield)
    !$acc end kernels 
    !$OMP END PARALLEL DO
    !
    return
end subroutine correct_mom

!_______________________________________________________________________________

subroutine update_reg_terms(nh_u,nh_p,nh_d,halo_u,nh_v,dli,n,psi,ls,nor,u_reg,eps_int,dl,dzc,dzf)
  !update the regularisation terms and do MPI communication
  !this subroutine should be called only once per timestep
  use mod_param  ,     only: no_outflow,is_outflow,IBM_h
  use mod_common_mpi,  only: MPI_PROC_NULL,left,right,front,back,top,bottom
  use mod_bound  ,     only: bounduvw
  implicit none
  !
  real(rp), intent(in ), dimension(3)                          :: dli
  integer , intent(in ), dimension(3)                          :: n
  integer , intent(in )                                        :: nh_u,nh_p,nh_v,nh_d
  real(rp), intent(in ), dimension(1-nh_d:)                    :: dzc,dzf
  real(rp), intent(in ), dimension(3)                          :: dl
  real(rp), intent(in )                                        :: u_reg,eps_int
  real(rp), intent(in ), dimension(1-nh_p:,1-nh_p:,1-nh_p:)    :: psi,ls
  real(rp), intent(in ), dimension(1-nh_v:,1-nh_v:,1-nh_v:,1:) :: nor
  integer , intent(in ), dimension(3)                          :: halo_u
  real(rp)                                                     :: maxreg
  !
  rex(:,:,:) = 0.0
  rey(:,:,:) = 0.0
  rez(:,:,:) = 0.0
  !
#if defined(_USE_ACDI)    
  call cmpt_regularization_acdi(nh_u,nh_p,nh_v,dli,n,psi,ls,nor,u_reg,eps_int)
#else
  call cmpt_regularization(nh_u,nh_p,nh_v,dli,n,psi,nor,u_reg,eps_int)
#endif
  !
  !if(is_outflow(1,2).AND.(back .eq.MPI_PROC_NULL)) then
   ! rey(:,n(2)-1,:) = 0.0   ! this one was good!
   ! rex(:,n(2),:) = 0.0
    !rez(:,:,IBM_h+1) = 0.0
  !endif
  !
  call bounduvw(cbcreg,n,bcreg,nh_d,nh_u,halo_u,no_outflow,dl,dzc,dzf,rex,rey,rez, &
                                                          psi,rex,rey,rez,0.0_rp,0.0_rp,.FALSE.) ! psi is used as dummy
  !
  return
end subroutine update_reg_terms

!_______________________________________________________________________________

subroutine make_vol_avg_vel(n,nh_d,nh_v,nh_u,u,v,w,rho)
  ! turn a density avergaed mixture velocity (input) into a volume averaged mixture velocity (output)
  use mod_param     , only: IBM_h,rho1,rho2,is_outflow
  use mod_common_mpi,  only: MPI_PROC_NULL,left,right,front,back,top,bottom
  implicit none
  !
  integer , intent(in   ), dimension(3)                       :: n
  integer , intent(in   )                                     :: nh_d,nh_v,nh_u
  real(rp), intent(inout), dimension(1-nh_u:,1-nh_u:,1-nh_u:) :: u,v,w
  real(rp), intent(in   ), dimension(     0:,     0:,     0:) :: rho
  integer                                                     :: i,j,k
  real(rp)                                                    :: rhox,rhoy,rhoz
  !
  integer                                                     :: x_start,x_stop
  integer                                                     :: y_start,y_stop
  integer                                                     :: z_start,z_stop
  x_start = 0
  x_stop  = 0
  y_start = 0
  y_stop  = 0
  z_start = 0
  z_stop  = 0   
  if    (is_outflow(0,1).AND.(left .eq.MPI_PROC_NULL)) then
    x_start = out_rk_safety
  elseif(is_outflow(1,1).AND.(right.eq.MPI_PROC_NULL)) then
    x_stop = out_rk_safety
  elseif(is_outflow(0,2).AND.(front.eq.MPI_PROC_NULL)) then
    y_start = out_rk_safety
  elseif(is_outflow(1,2).AND.(back .eq.MPI_PROC_NULL)) then
    y_stop = out_rk_safety
  elseif(is_outflow(0,3)) then
    z_start = out_rk_safety
  elseif(is_outflow(1,3)) then
    z_stop = out_rk_safety
  endif
  !
  do k=IBM_h+1+z_start,n(3)-z_stop
    do j=1+y_start,n(2)-y_stop
      do i=x_start+1,n(1)-x_stop
        !
        rhox = 0.5d0*(rho(i,j,k)+rho(i+1,j  ,k  ))
        rhoy = 0.5d0*(rho(i,j,k)+rho(i  ,j+1,k  ))
        rhoz = 0.5d0*(rho(i,j,k)+rho(i  ,j  ,k+1))
        !
        u(i,j,k) = u(i,j,k) + (rho1-rho2)*rex(i,j,k)/rhox                                     
        v(i,j,k) = v(i,j,k) + (rho1-rho2)*rey(i,j,k)/rhoy
        w(i,j,k) = w(i,j,k) + (rho1-rho2)*rez(i,j,k)/rhoz
      enddo
    enddo
  enddo
  !
  return
end subroutine make_vol_avg_vel

!_______________________________________________________________________________

subroutine regularization_psi_cdi_ibm(n,dli,nh_d,nh_v,nh_u,dsdt)
  !s is volume fraction/phi
  use mod_param     , only: IBM_h
  implicit none
  !
  integer , intent(in ), dimension(3)                       :: n
  real(rp), intent(in ), dimension(3)                       :: dli 
  integer , intent(in )                                     :: nh_d,nh_v,nh_u
  real(rp), intent(out), dimension(1:,1:,1:)                :: dsdt
  !
  integer  :: i,j,k,ip,im,jp,jm,kp,km
  !
  dsdt(:,:,:) = 0.d0
  !
  !
  do k=1+IBM_h,n(3)
    do j=1,n(2)
      do i=1,n(1)
        !
        ip = i + 1
        im = i - 1
        jp = j + 1
        jm = j - 1
        kp = k + 1
        km = k - 1
        !
        dsdt(i,j,k)    =    (rex(i,j,k) - rex(im,j ,k )) * dli(1) + &
                            (rey(i,j,k) - rey(i ,jm,k )) * dli(2) + &
                            (rez(i,j,k) - rez(i ,j ,km)) * dli(3) 
      enddo
    enddo
  enddo
  !
  return
end subroutine regularization_psi_cdi_ibm

!_______________________________________________________________________________

subroutine cmpt_regularization(nh_u,nh_p,nh_v,dli,n,s,nor,u_reg,eps_int)
  !
  use mod_param     , only: IBM_h
  !
  implicit none
  !
  real(rp), intent(in ), dimension(3)                          :: dli
  integer , intent(in ), dimension(3)                          :: n
  integer , intent(in )                                        :: nh_u,nh_p,nh_v
  real(rp), intent(in )                                        :: u_reg,eps_int
  integer                                                      :: i,j,k,ip,im,jp,jm,kp,km
  real(rp), intent(in ), dimension(1-nh_p:,1-nh_p:,1-nh_p:)    :: s
  real(rp), intent(in ), dimension(1-nh_v:,1-nh_v:,1-nh_v:,1:) :: nor
  real(rp)                                                     :: regxp,regyp,regzp
  real(rp)                                                     :: psifxp,psifyp,psifzp
  real(rp)                                                     :: phix,phiy,phiz
  real(rp)                                                     :: gradpxp,gradpyp,gradpzp
  real(rp)                                                     :: norx,nory,norz
  !
  !
  do k=1+IBM_h,n(3)
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
#else
        phix            =   0.5*( s(i ,j,k  )+ s(ip,j,k  ))
        norx           =   0.5*(nor(i ,j,k,1)+nor(ip,j,k,1))
        gradpxp        =   ( s(ip,j,k) - s(i,j,k)  ) * dli(1)       !grad(phi) at i+0.5
        !
        psifxp         =   ( s(i,j,k)*(1.d0-s(i,j,k))*nor(i,j,k,1) + s(ip,j,k)*(1.d0-s(ip,j,k))*nor(ip,j,k,1))*0.5 !phi*(1-phi)*nor at i+0.5 as p.15 Shahab.M
        !psifxp         =   ( phix*(1.d0-phix)*norx)
        !
        regxp          =   u_reg*( eps_int*gradpxp - psifxp)        !regularization at i+0.5
#endif
        !
        phiy           =   0.5*( s(i ,j,k  )+ s(i,jp,k  ))
        nory           =   0.5*(nor(i,j ,k,2)+nor(i,jp,k,2))
        gradpyp        =   ( s(i,jp,k) - s(i,j,k)  ) * dli(2)       !grad(phi) at j+0.5
        !
        psifyp         =   ( s(i,j,k)*(1.d0-s(i,j,k))*nor(i,j,k,2) + s(i,jp,k)*(1.d0-s(i,jp,k))*nor(i,jp,k,2))*0.5 !phi*(1-phi)*nor at j+0.5 as p.15 Shahab.M
        !psifyp         =   ( phiy*(1.d0-phiy)*nory)
        !
        regyp          =   u_reg*( eps_int*gradpyp - psifyp)        !regularization at j+0.5
        !
        phiz           =   0.5*( s(i ,j,k  )+ s(i,j,kp  ))
        norz           =   0.5*(nor(i,j ,k,3)+nor(i,j,kp,3))
        gradpzp        =   ( s(i,j,kp) - s(i,j,k)  ) * dli(3)       !grad(phi) at k+0.5
        !
        psifzp         =   ( s(i,j,k)*(1.d0-s(i,j,k))*nor(i,j,k,3) + s(i,j,kp)*(1.d0-s(i,j,kp))*nor(i,j,kp,3))*0.5 !phi*(1-phi)*nor at k+0.5 as p.15 Shahab.M
        !psifzp         =   ( phiz*(1.d0-phiz)*norz)
        !
        regzp          =   u_reg*( eps_int*gradpzp - psifzp)        !regularization at k+0.5
        !
        rex(i,j,k) = regxp
        rey(i,j,k) = regyp
        rez(i,j,k) = regzp
        !
      enddo
    enddo
  enddo
  !
  return
end subroutine cmpt_regularization

!_______________________________________________________________________________

subroutine cmpt_regularization_acdi(nh_u,nh_p,nh_v,dli,n,s,ls,nor,u_reg,eps_int)
  !
  use mod_param     , only: IBM_h
  !
  implicit none
  !
  real(rp), intent(in ), dimension(3)                          :: dli
  integer , intent(in ), dimension(3)                          :: n
  integer , intent(in )                                        :: nh_u,nh_p,nh_v
  real(rp), intent(in )                                        :: u_reg,eps_int
  integer                                                      :: i,j,k,ip,im,jp,jm,kp,km
  real(rp), intent(in ), dimension(1-nh_p:,1-nh_p:,1-nh_p:)    :: s,ls
  real(rp), intent(in ), dimension(1-nh_v:,1-nh_v:,1-nh_v:,1:) :: nor
  real(rp)                                                     :: regxp,regyp,regzp
  real(rp)                                                     :: psifxp,psifyp,psifzp
  real(rp)                                                     :: gradpxp,gradpyp,gradpzp
  real(rp)                                                     :: norx,nory,norz
  real(rp)                                                     :: lsx,lsy,lsz
  !
  !
  do k=1+IBM_h,n(3)
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
#else
        lsx            =   0.5*( ls(i ,j,k  )+ ls(ip,j,k  ))
        norx           =   0.5*(nor(i ,j,k,1)+nor(ip,j,k,1))
        gradpxp        =   ( s(ip,j,k) - s(i,j,k)  ) * dli(1)       !grad(phi) at i+0.5
        !
        !psifxp         =   ( 0.25_rp*( 1._rp-tanh(ls(i ,j,k)/(2._rp*eps_int))**2 )*nor(i ,j,k,1) &
        !                   + 0.25_rp*( 1._rp-tanh(ls(ip,j,k)/(2._rp*eps_int))**2 )*nor(ip,j,k,1))*0.5 !phi*(1-phi)*nor at i+0.5 as p.15 Shahab.M
        psifxp         =   ( 0.25_rp*( 1._rp-tanh(lsx/(2._rp*eps_int))**2 )*norx)
        !
        regxp          =   u_reg*( eps_int*gradpxp - psifxp)        !regularization at i+0.5
#endif
        !
        lsy            =   0.5*( ls(i,j ,k  )+ ls(i,jp,k  ))
        nory           =   0.5*(nor(i,j ,k,2)+nor(i,jp,k,2))
        gradpyp        =   ( s(i,jp,k) - s(i,j,k)  ) * dli(2)       !grad(phi) at j+0.5
        !
        !psifyp         =   ( 0.25_rp*( 1._rp-tanh(ls(i,j ,k)/(2._rp*eps_int))**2 )*nor(i,j ,k,2) &
        !                   + 0.25_rp*( 1._rp-tanh(ls(i,jp,k)/(2._rp*eps_int))**2 )*nor(i,jp,k,2))*0.5 !phi*(1-phi)*nor at j+0.5 as p.15 Shahab.M
        psifyp         =   ( 0.25_rp*( 1._rp-tanh(lsy/(2._rp*eps_int))**2 )*nory)
        !
        regyp          =   u_reg*( eps_int*gradpyp - psifyp)        !regularization at j+0.5
        !
        lsz            =   0.5*( ls(i,j ,k  )+ ls(i,j,kp  ))
        norz           =   0.5*(nor(i,j ,k,3)+nor(i,j,kp,3))
        gradpzp        =   ( s(i,j,kp) - s(i,j,k)  ) * dli(3)       !grad(phi) at k+0.5
        !
        !psifzp         =   ( 0.25_rp*( 1._rp-tanh(ls(i,j,k )/(2._rp*eps_int))**2 )*nor(i,j,k ,3) &
        !                   + 0.25_rp*( 1._rp-tanh(ls(i,j,kp)/(2._rp*eps_int))**2 )*nor(i,j,kp,3))*0.5 !phi*(1-phi)*nor at k+0.5 as p.15 Shahab.M
        psifzp         =   ( 0.25_rp*( 1._rp-tanh(lsz/(2._rp*eps_int))**2 )*norz)
        !
        regzp          =   u_reg*( eps_int*gradpzp - psifzp)        !regularization at k+0.5
        !
        rex(i,j,k) = regxp
        rey(i,j,k) = regyp
        rez(i,j,k) = regzp
        !
      enddo
    enddo
  enddo
  !
  return
end subroutine cmpt_regularization_acdi
!
!#################################################
!###     End of:  regularisation subrouties    ###
!#################################################
!
!#################################################
!###   Start of:  postprocessing subrouties    ###
!#################################################
!
subroutine calculate_tot_surface(n,dli,dl,nh_p,psi,surface_density,contact_line_density,contact_area_denity,psi_tr)
  !s is volume fraction/phi
  use mod_param     , only: IBM_h
  implicit none
  !
  integer , intent(in ), dimension(3)                          :: n
  real(rp), intent(in ), dimension(3)                          :: dli,dl
  integer , intent(in )                                        :: nh_p
  real(rp), intent(in ), dimension(1-nh_p:,1-nh_p:,1-nh_p:)    :: psi
  real(rp), intent(out), dimension(1:n(1)  ,1:n(2) ,1:n(3))    :: surface_density,contact_area_denity,contact_line_density
  real(rp), intent(out), dimension(1-nh_p:n(1)+nh_p, &
                                   1-nh_p:n(2)+nh_p, &
                                   1-nh_p:n(3)+nh_p)           :: psi_tr
  real(rp)                                                     :: limiter
  integer                                                      :: i,j,k,ip,im,jp,jm,kp,km
  !   
  surface_density(:,:,:)      = 0.d0
  contact_area_denity(:,:,:)  = 0.d0
  contact_line_density(:,:,:) = 0.d0
  psi_tr(:,:,:)               = 0.0
  !
  limiter = 0.6   ! has to be double the tagging threshold!
  !
  do k=0,n(3)+1
    do j=0,n(2)+1
      do i=0,n(1)+1
        psi_tr(i,j,k) = (1.0/(1.0-limiter))*(min(max(psi(i,j,k),0.5*limiter),1.0-0.5*limiter)-0.5*limiter)
      enddo
    enddo
  enddo
  !
  ! Interfacial surface area
  !
  do k=1+IBM_h,n(3)
    do j=1,n(2)
      do i=1,n(1)
        !
        ip = i + 1
        im = i - 1
        jp = j + 1
        jm = j - 1
        kp = k + 1
        km = k - 1
        ! 
        surface_density(i,j,k) = sqrt((0.5*(psi_tr(ip,j ,k ) - psi_tr(im,j ,k )) * dli(1))**2 + &
                                      (0.5*(psi_tr(i ,jp,k ) - psi_tr(i ,jm,k )) * dli(2))**2 + &
                                      (0.5*(psi_tr(i ,j ,kp) - psi_tr(i ,j ,km)) * dli(3))**2)
      enddo
    enddo
  enddo
  !
  ! contact line length
  !
  k=IBM_h
  do j=1,n(2)
    do i=1,n(1)
      !
      ip = i + 1
      im = i - 1
      jp = j + 1
      jm = j - 1
      kp = k + 1
      km = k - 1
      !
      contact_line_density(i,j,kp) = sqrt((0.5*(0.5*(psi_tr(ip,j ,k ) - psi_tr(im,j ,k ))                 &
                                              + 0.5*(psi_tr(ip,j ,kp) - psi_tr(im,j ,kp))) * dli(1))**2 + &
                                          (0.5*(0.5*(psi_tr(i ,jp,k ) - psi_tr(i ,jm,k ))                 &
                                              + 0.5*(psi_tr(i ,jp,kp) - psi_tr(i ,jm,kp))) * dli(2))**2)
    enddo
  enddo
  !
  ! contact area
  !
  k=IBM_h
  do j=1,n(2)
    do i=1,n(1)
      !
      ip = i + 1
      im = i - 1
      jp = j + 1
      jm = j - 1
      kp = k + 1
      km = k - 1
      !
      contact_area_denity(i,j,kp) = (0.5*(psi_tr(i,j,kp) + psi_tr(i,j,k)))
    enddo
  enddo
  !
  !
  return
end subroutine calculate_tot_surface
!#################################################
!###     End of:  postprocessing subrouties    ###
!#################################################
end module mod_boiling

