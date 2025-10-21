module mod_tagging_acdi_ibm
  !
  ! note: not yet on GPU
  ! 
  use mpi
  use mod_bound     , only: boundp
  use mod_common_mpi, only: myid,comm_cart,ierr,ijk_start, &
                            left,right,front,back,top,bottom
  use mod_boiling   , only: cbcpsi_n,bcpsi_n,calculate_tot_surface
  use mod_types     , only: rp
  !
  implicit none 
  !
  real(rp), parameter :: vof_th = 0.3_rp
  !
  private
  public  :: tagging_acdi_ibm
  !
  contains
  !
  subroutine tagging_acdi_ibm(n,dims,datadir_ta,dl,nh_d,nh_v,nh_u,halo_v,dzc,dzf,phi,u,v,w,istep,time,src_ph)
    !
    use mod_param,   only: dli,IBM_h,l,itot,jtot,ktot,ng
    use mod_boiling, only: rex,rey,rez
    use mod_output,  only: write_visu_3d
    !
    implicit none
    !
    integer         , intent(in), dimension(3)                       :: n
    integer         , intent(in), dimension(3)                       :: dims
    character(len=*), intent(in)                                     :: datadir_ta
    real(rp)        , intent(in), dimension(3)                       :: dl
    integer         , intent(in)                                     :: nh_d,nh_v,nh_u
    integer         , intent(in), dimension(3)                       :: halo_v
    real(rp)        , intent(in), dimension(1-nh_d:)                 :: dzc,dzf
    real(rp)        , intent(in), dimension(1-nh_v:,1-nh_v:,1-nh_v:) :: phi,src_ph
    real(rp)        , intent(in), dimension(1-nh_u:,1-nh_u:,1-nh_u:) :: u,v,w
    integer         , intent(in)                                     :: istep
    real(rp)        , intent(in)                                     :: time
    !
    integer(8), dimension(0:n(1)+1,0:n(2)+1,0:n(3)+1) :: vofTag,vofTag_orig,vofTag_loc
    real(rp)  , dimension(0:n(1)+1,0:n(2)+1,0:n(3)+1) :: phi_f,vofTag_r,aux
    real(rp)  , dimension(1:n(1)  ,1:n(2) ,1:n(3))    :: surface_density,contact_area_denity,contact_line_density
    real(rp)  , dimension(1-nh_v:n(1)+nh_v, &
                          1-nh_v:n(2)+nh_v, &
                          1-nh_v:n(3)+nh_v)           :: psi_tr
    real(rp) :: xmin_,xmax_,ymin_,ymax_,zmin_,zmax_,phix,phiy,phiz
    real(rp) :: rdnInit
    integer  :: i,j,k,b, cellcount,ii,jj,kk,loc
    integer(4)  :: nd
    integer, dimension(1:6) :: taglist
    integer(4) :: id, idst,id_loc
    integer(4), allocatable, dimension(:,:) :: procShare, share,idShare, orderField
    integer(4), allocatable, dimension(:)   :: myid_out, drid,glid
    integer(8), allocatable, dimension(:)   :: idloc
    real(8),    allocatable, dimension(:)   :: xd,yd,zd,ud,vd,wd,vold
    real(8),    allocatable, dimension(:)   :: xmin,xmax,ymin,ymax,zmin,zmax ! min and max location of bubble
    real(8),    allocatable, dimension(:)   :: mtot,mmin,mmax,psimax !mass transfer
    real(8),    allocatable, dimension(:)   :: sur,vol,drya,conl !bubble surface, bubble volume, dryspot area, contactline length
    real(8),    allocatable, dimension(:)   :: voldx,voldy,voldz ! for weighting of velocities
    integer(8)                              :: myid_offset
    integer(8)                              :: search_tag

    character(len=3)                        :: proidch
    character(len=9)                        :: fldnum
    !
    write(fldnum,'(i9.9)') istep
    !@cuf attributes(managed) :: dzc, dzf, vofTag_r
    ! 
    !call random_number(rdnInit)
    ! idst = floor(rdnInit*100000+100000*myid) ! avoid possible ID duplication between processors. May fail for n_d >>10000???
    myid_offset = 1000    ! per rank: mo more than myid_offset bubble counts possible!
    idst        = myid_offset*myid   ! 0 ! avoid possible ID duplication between processors. May fail for n_d >>10000???
    id          = idst+1
    id_loc      = 1
    nd          = 0
    !
    ! initialize vofTag and vof_f
    !
    do k=0,n(3)+1
      do j=0,n(2)+1
        do i=0,n(1)+1
          vofTag(i,j,k)     = 0
          vofTag_loc(i,j,k) = 0
          aux(i,j,k)        = 0.0_rp
          vofTag_r(i,j,k)   = 0.0_rp
          phi_f(i,j,k)      = phi(i,j,k)
        enddo
      enddo
    enddo
    !
    ! main loop for tagging droplets 
    !
    do k=1+IBM_h,n(3)
      do j=1,n(2)
        do i=1,n(1)
          if(phi_f(i,j,k).gt.vof_th) then
            phi_f(i,j,k)  = 0._rp
            vofTag(i,j,k) = id
            vofTag_loc(i,j,k) = id_loc
            call recursive_tagging(n,i,j,k,id,id_loc,phi_f,vofTag,vofTag_loc)
            id_loc = id_loc + 1
            id     = id+1
            nd     = nd+1
          endif
        enddo
      enddo
    enddo
    !
    ! conversion to real(rp) to set_bc
    !
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          vofTag_r(i,j,k) = 1._rp*vofTag(i,j,k)
        enddo
      enddo
    enddo
    call boundp(cbcpsi_n,n,bcpsi_n,nh_d,nh_v,halo_v,dl,dzc,dzf,vofTag_r) ! meaningful bc: N or P
    !
    ! re-conversion to integer
    !
    do k=0,n(3)+1
      do j=0,n(2)+1
        do i=0,n(1)+1
          vofTag_orig(i,j,k) = int(vofTag_r(i,j,k),8)
        enddo
      enddo
    enddo
    !
    !
    ! Now add the surrounding to the bubble
    aux(:,:,:) = vofTag_r(:,:,:)
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          if(aux(i,j,k).eq.0.0_rp) then
            taglist(1) = aux(i-1,j  ,k  )
            taglist(2) = aux(i+1,j  ,k  )
            taglist(3) = aux(i  ,j+1,k  )
            taglist(4) = aux(i  ,j-1,k  )
            taglist(5) = aux(i  ,j  ,k-1)
            taglist(6) = aux(i  ,j  ,k+1)
            vofTag_r(i,j,k) = maxval(taglist)
          endif 
        enddo
      enddo
    enddo
    call boundp(cbcpsi_n,n,bcpsi_n,nh_d,nh_v,halo_v,dl,dzc,dzf,vofTag_r) ! meaningful bc: N or P
    ! do it a second time!
    aux(:,:,:) = vofTag_r(:,:,:)
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          if(aux(i,j,k).eq.0.0_rp) then
            taglist(1) = aux(i-1,j  ,k  )
            taglist(2) = aux(i+1,j  ,k  )
            taglist(3) = aux(i  ,j+1,k  )
            taglist(4) = aux(i  ,j-1,k  )
            taglist(5) = aux(i  ,j  ,k-1)
            taglist(6) = aux(i  ,j  ,k+1)
            vofTag_r(i,j,k) = maxval(taglist)
          endif 
        enddo
      enddo
    enddo
    call boundp(cbcpsi_n,n,bcpsi_n,nh_d,nh_v,halo_v,dl,dzc,dzf,vofTag_r) ! meaningful bc: N or P
    !
    ! re-conversion to integer
    !
    aux(:,:,:) = vofTag_r(:,:,:)
    do k=0,n(3)+1
      do j=0,n(2)+1
        do i=0,n(1)+1
          vofTag(i,j,k) = int(vofTag_r(i,j,k),8)
          if(floor(vofTag_r(i,j,k)/myid_offset).EQ.myid) then
            aux(i,j,k)        = 0.0_rp ! remove native tags so that only forign tags are left after this step!
            vofTag_loc(i,j,k) = mod(vofTag(i,j,k),myid_offset)
          endif
        enddo
      enddo
    enddo
    !
    ! give local tags for forign global tags
    !
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          if(aux(i,j,k).GT.0.0_rp) then
            search_tag = int(aux(i,j,k),8)
            aux(i,j,k)  = 0._rp
            vofTag_loc(i,j,k) = id_loc
            call recursive_remove_tag(n,i,j,k,id_loc,aux,vofTag_loc,search_tag)
            id_loc = id_loc + 1
            nd     = nd+1
          endif
        enddo
      enddo
    enddo
    !
    !communicate local tags
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          aux(i,j,k) = 1._rp*vofTag_loc(i,j,k)
        enddo
      enddo
    enddo
    call boundp(cbcpsi_n,n,bcpsi_n,nh_d,nh_v,halo_v,dl,dzc,dzf,aux) ! meaningful bc: N or P
    !
    ! re-conversion to integer
    !
    do k=0,n(3)+1
      do j=0,n(2)+1
        do i=0,n(1)+1
          vofTag_loc(i,j,k) = int(aux(i,j,k),8)
        enddo
      enddo
    enddo
    ! 
    ! write(proidch,'(i3.3)') myid
    ! open(90+myid,file='proc'//proidch)
    ! do i=0,n(1)+1
    !   do j=0,n(2)+1
    !     do k=0,n(3)+1
    !       write(90+myid,'(4I5.4)') i,j,k,vofTag(i,j,k)
    !     enddo
    !   enddo
    ! enddo
    ! close(92)
    !
    ! print*,'max vof', maxval(vof_f), myid
    !
    ! allocate droplet properties
    !
    call calculate_tot_surface(n,dli,dl,nh_v,phi,surface_density,contact_line_density,contact_area_denity,psi_tr)
    !
    allocate(xd(nd),yd(nd),zd(nd),ud(nd),vd(nd),wd(nd),vold(nd), & 
             idloc(nd),share(nd,6),myid_out(nd),procShare(nd,6), &
             idShare(nd,6), orderField(nd,6),drid(nd),glid(nd),xmin(nd), &
             xmax(nd),ymin(nd),ymax(nd),zmin(nd),zmax(nd) &
             ,mtot(nd),mmin(nd),mmax(nd),sur(nd),vol(nd), &
             drya(nd),conl(nd),voldx(nd),voldy(nd),voldz(nd),psimax(nd))
    xd = 0._rp
    yd = 0._rp
    zd = 0._rp
    ud = 0._rp
    vd = 0._rp
    wd = 0._rp
    xmin = 1e7_rp
    xmax = -1e7_rp
    ymin = 1e7_rp
    ymax = -1e7_rp
    zmin = 1e7_rp
    zmax = -1e7_rp
    psimax = 0.0_rp
    mtot = 0._rp
    mmin = 1e7_rp
    mmax = -1e7_rp
    sur = 0._rp
    vol = 0._rp
    drya = 0._rp
    conl = 0._rp
    vold = 0._rp
    voldx = 0._rp
    voldy = 0._rp
    voldz = 0._rp
    id = 0
    cellcount = 0
    idShare = 0
    procShare = -1
    share = 0
    drid = 0
    glid = 0
    !
    ! Loop for properties characterization
    !
    do k=1+IBM_h,n(3)
      do j=1,n(2)
        do i=1,n(1)
          !
          if(vofTag(i,j,k).gt.0) then
            ! Default values
            xmin_ = 1e7_rp
            xmax_ = -1e7_rp
            ymin_ = 1e7_rp
            ymax_ = -1e7_rp
            zmin_ = 1e7_rp
            zmax_ = -1e7_rp
            if((phi(i,j,k).LE.0.5_rp).AND.(phi(i+1,j,k).GT.0.5_rp)) &
                        xmin_ = (i+ijk_start(1)-0.5_rp)*dl(1)+dl(1)*(0.5_rp      -phi(i,j,k)) &
                                                                   /(phi(i+1,j,k)-phi(i,j,k))
            if((phi(i,j,k).GE.0.5_rp).AND.(phi(i+1,j,k).LT.0.5_rp)) &
                        xmax_ = (i+ijk_start(1)-0.5_rp)*dl(1)+dl(1)*(0.5_rp      -phi(i,j,k)) &
                                                                   /(phi(i+1,j,k)-phi(i,j,k))
            if((phi(i,j,k).LE.0.5_rp).AND.(phi(i,j+1,k).GT.0.5_rp)) &
                        ymin_ = (j+ijk_start(2)-0.5_rp)*dl(2)+dl(2)*(0.5_rp      -phi(i,j,k)) &
                                                                   /(phi(i,j+1,k)-phi(i,j,k))
            if((phi(i,j,k).GE.0.5_rp).AND.(phi(i,j+1,k).LT.0.5_rp)) &
                        ymax_ = (j+ijk_start(2)-0.5_rp)*dl(2)+dl(2)*(0.5_rp      -phi(i,j,k)) &
                                                                   /(phi(i,j+1,k)-phi(i,j,k))
            if((phi(i,j,k).LE.0.5_rp).AND.(phi(i,j,k+1).GT.0.5_rp)) &
                        zmin_ = (k+ijk_start(3)-0.5_rp)*dl(3)+dl(3)*(0.5_rp      -phi(i,j,k)) &
                                                                   /(phi(i,j,k+1)-phi(i,j,k))
            if((phi(i,j,k).GE.0.5_rp).AND.(phi(i,j,k+1).LT.0.5_rp)) &
                        zmax_ = (k+ijk_start(3)-0.5_rp)*dl(3)+dl(3)*(0.5_rp      -phi(i,j,k)) &
                                                                   /(phi(i,j,k+1)-phi(i,j,k))
            ! Now the boundaries
            if((phi(i-1,j,k).GE.0.5_rp).AND.(phi(i,j,k).GE.0.5_rp).AND.((i+ijk_start(1)).EQ.1      )) xmin_ = 0.0_rp
            if((phi(i+1,j,k).GE.0.5_rp).AND.(phi(i,j,k).GE.0.5_rp).AND.((i+ijk_start(1)).EQ.itot   )) xmax_ = l(1)
            !
            if((phi(i,j-1,k).GE.0.5_rp).AND.(phi(i,j,k).GE.0.5_rp).AND.((j+ijk_start(2)).EQ.1      )) ymin_ = 0.0_rp
            if((phi(i,j+1,k).GE.0.5_rp).AND.(phi(i,j,k).GE.0.5_rp).AND.((j+ijk_start(2)).EQ.jtot   )) ymax_ = l(2)
            !
            if((phi(i,j,k).GE.0.5_rp).AND.((k+ijk_start(3)).EQ.IBM_h+1)) zmin_ = dl(3)*IBM_h   !(phi(i,j,k-1).GE.0.5_rp).AND.
            if((phi(i,j,k+1).GE.0.5_rp).AND.(phi(i,j,k).GE.0.5_rp).AND.((k+ijk_start(3)).EQ.ktot   )) zmax_ = l(3)
            !
            phix = 0.5d0*(phi(i,j,k)+phi(i+1,j  ,k  ))
            phiy = 0.5d0*(phi(i,j,k)+phi(i  ,j+1,k  ))
            phiz = 0.5d0*(phi(i,j,k)+phi(i  ,j  ,k+1))
            !
            id        = vofTag_loc(i,j,k)
            drid(id)  = id!vofTag(i,j,k)
            glid(id)  = vofTag(i,j,k)
            idloc(id) = vofTag(i,j,k)
            sur(id)   = sur(id)+ surface_density(i,j,k)*dl(1)*dl(2)*dl(3)
            drya(id)  = drya(id)+ contact_area_denity(i,j,k)*dl(1)*dl(2)
            conl(id)  = conl(id)+ contact_line_density(i,j,k)*dl(1)*dl(2)
            vol(id)   = vol(id) + psi_tr(i,j,k)
            mtot(id)  = mtot(id) + src_ph(i,j,k)*dl(1)*dl(2)*dl(3)
            mmin(id)  = min(mmin(id),src_ph(i,j,k)*dl(1)*dl(2)*dl(3))
            mmax(id)  = max(mmax(id),src_ph(i,j,k)*dl(1)*dl(2)*dl(3))
            psimax(id)= max(psimax(id),phi(i,j,k))
            xmin(id)  = min(xmin(id),xmin_)
            xmax(id)  = max(xmax(id),xmax_)
            ymin(id)  = min(ymin(id),ymin_)
            ymax(id)  = max(ymax(id),ymax_)
            zmin(id)  = min(zmin(id),zmin_)
            zmax(id)  = max(zmax(id),zmax_)
            xd(id)    = xd(id) + psi_tr(i,j,k)*(i+ijk_start(1)-0.5_rp)*dl(1)
            yd(id)    = yd(id) + psi_tr(i,j,k)*(j+ijk_start(2)-0.5_rp)*dl(2)
            zd(id)    = zd(id) + psi_tr(i,j,k)*(k+ijk_start(3)-0.5_rp)*dl(3)
            ud(id)    = ud(id) + u(i,j,k)*phix-(rex(i,j,k))
            vd(id)    = vd(id) + v(i,j,k)*phiy-(rey(i,j,k))
            wd(id)    = wd(id) + w(i,j,k)*phiz-(rez(i,j,k))
            voldx(id) = voldx(id) + phix
            voldy(id) = voldy(id) + phiy
            voldz(id) = voldz(id) + phiz
            vold(id)  = vold(id) + phi(i,j,k)
            !
          endif
          !
          ! find neighbors for all native tags
          if(vofTag_orig(i,j,k).gt.0) then
            id        = vofTag_loc(i,j,k)
            ! print*, vof(i,j,k), i,j,k
            cellcount = 1
            if((i.eq.1).and.share(id,1).eq.0.and.vofTag_orig(0,j,k).ne.0) then
              ! print*, 'i0i0i0i0i'
              share(id,1)     = 1
              procShare(id,1) = left
              idShare(id,1)   = mod(vofTag_orig(0,j,k),myid_offset)  !convert to local tag 
              ! if (vofTag(0,j,k).ne.0)      idShare(id,1)   = vofTag(0,j,k)
            endif
            if((i.eq.n(1)).and.share(id,2).eq.0.and.vofTag_orig(n(1)+1,j,k).ne.0) then
              ! print*, 'n1n1n1n1n1'
              share(id,2)     = 1
              procShare(id,2) = right
              idShare(id,2)   = mod(vofTag_orig(n(1)+1,j,k),myid_offset)  !convert to local tag 
            endif
            if((j.eq.1).and.share(id,3).eq.0.and.vofTag_orig(i,0,k).ne.0) then
              share(id,3)     = 1
              procShare(id,3) = front
              idShare(id,3)   = mod(vofTag_orig(i,0,k),myid_offset)  !convert to local tag 
              ! print*, 'j1', idShare(id,:)
            endif
            if((j.eq.n(2)).and.share(id,4).eq.0.and.vofTag_orig(i,n(2)+1,k).ne.0) then
              share(id,4)     = 1
              procShare(id,4) = back
              idShare(id,4)   = mod(vofTag_orig(i,n(2)+1,k),myid_offset)  !convert to local tag 
              ! print*, 'n2', idShare(id,:)
            endif
            if((k.eq.1).and.share(id,5).eq.0.and.vofTag_orig(i,j,0).ne.0) then
              share(id,5)     = 1
              procShare(id,5) = bottom
              idShare(id,5)   = mod(vofTag_orig(i,j,0),myid_offset)  !convert to local tag 
            endif
            if((k.eq.n(3)).and.share(id,6).eq.0.and.vofTag_orig(i,j,n(3)+1).ne.0) then
              share(id,6)     = 1
              procShare(id,6) = top
              idShare(id,6)   = mod(vofTag_orig(i,j,n(3)+1),myid_offset)  !convert to local tag 
            endif
          endif
        enddo
      enddo
    enddo
    ! neighbors for native tags have priority over forign tags
    do k=1+IBM_h,n(3)
      do j=1,n(2)
        do i=1,n(1)
          ! find neighbors for forign tags
          if((vofTag(i,j,k).gt.0)) then    !.AND.(floor(float(vofTag(i,j,k))/float(myid_offset)).NE.myid)
            id        = vofTag_loc(i,j,k)
            ! print*, vof(i,j,k), i,j,k
            cellcount = 1
            if((i.eq.1).and.share(id,1).eq.0.and.(vofTag(0,j,k).eq.vofTag(i,j,k))) then
              ! print*, 'i0i0i0i0i'
              share(id,1)     = 1
              procShare(id,1) = left
              idShare(id,1)   = vofTag_loc(0,j,k)
              ! if (vofTag(0,j,k).ne.0)      idShare(id,1)   = vofTag(0,j,k)
            endif
            if((i.eq.n(1)).and.share(id,2).eq.0.and.(vofTag(i,j,k).eq.vofTag(n(1)+1,j,k))) then
              ! print*, 'n1n1n1n1n1'
              share(id,2)     = 1
              procShare(id,2) = right
              idShare(id,2)   = vofTag_loc(n(1)+1,j,k)
            endif
            if((j.eq.1).and.share(id,3).eq.0.and.(vofTag(i,0,k).eq.vofTag(i,j,k))) then
              share(id,3)     = 1
              procShare(id,3) = front
              idShare(id,3)   = vofTag_loc(i,0,k)
              ! print*, 'j1', idShare(id,:)
            endif
            if((j.eq.n(2)).and.share(id,4).eq.0.and.(vofTag(i,j,k).eq.vofTag(i,n(2)+1,k))) then
              share(id,4)     = 1
              procShare(id,4) = back
              idShare(id,4)   = vofTag_loc(i,n(2)+1,k) 
              ! print*, 'n2', idShare(id,:)
            endif
            if((k.eq.1).and.share(id,5).eq.0.and.(vofTag(i,j,0).eq.vofTag(i,j,k))) then
              share(id,5)     = 1
              procShare(id,5) = bottom
              idShare(id,5)   = vofTag_loc(i,j,0)
            endif
            if((k.eq.n(3)).and.share(id,6).eq.0.and.(vofTag(i,j,n(3)+1).eq.vofTag(i,j,k))) then
              share(id,6)     = 1
              procShare(id,6) = top
              idShare(id,6)   = vofTag_loc(i,j,n(3)+1)
            endif
          endif
        enddo
      enddo
    enddo
    !
    !do k=0,n(3)+1
    !  do j=0,n(2)+1
    !    do i=0,n(1)+1
    !      aux(i,j,k) = 1._rp*vofTag_loc(i,j,k)
    !    enddo
    !  enddo
    ! enddo
    !
    !call write_visu_3d('data/','tag1_fld_'//fldnum//'.bin','log_visu_3d.out','tag_loc', &
    !                 (/1,1,1/),(/ng(1),ng(2),ng(3)/),(/1,1,1/),time,istep, &
    !                 aux(1:n(1),1:n(2),1:n(3)))
    !call write_visu_3d('data/','tag2_fld_'//fldnum//'.bin','log_visu_3d.out','tag_glob', &
    !                 (/1,1,1/),(/ng(1),ng(2),ng(3)/),(/1,1,1/),time,istep, &
    !                 vofTag_r(1:n(1),1:n(2),1:n(3)))
    
    !do share tags for forign tags
    !do id=1,nd
    !  if(floor(float(glid(id))/float(myid_offset)).NE.myid) then
    !    share(id,1)     = 1
    !    idShare(id,1)   = mod(glid(id),int(myid_offset,4))  !convert to local tag 
    !    if(floor(float(glid(id))/float(myid_offset)).EQ.left)then
    !      procShare(id,1) = left
    !    elseif(floor(float(glid(id))/float(myid_offset)).EQ.right)then
    !      procShare(id,2) = right
    !    elseif(floor(float(glid(id))/float(myid_offset)).EQ.front)then
    !      procShare(id,3) = front
    !    elseif(floor(float(glid(id))/float(myid_offset)).EQ.back)then
     !     procShare(id,4) = back
    !    elseif(floor(float(glid(id))/float(myid_offset)).EQ.bottom)then
     !     procShare(id,5) = bottom
    !    elseif(floor(float(glid(id))/float(myid_offset)).EQ.top)then
    !      procShare(id,6) = top
    !    endif
     ! endif
   ! enddo 
    !
    !
    do i=1,nd
      do j=1,6
        orderField(i,j) = myid*100+i*10+J
        ! print*,'proc', procShare(i,:), myid
        ! print*, nd, myid
        ! print*,'id  ', idShare(i,:), myid
      enddo
    enddo
    !
    myid_out = 0
    do id=1,nd
      myid_out(id) = myid
      if(vol(id).GT.0.0) xd(id)   = xd(id)/vol(id)    
      if(vol(id).GT.0.0) yd(id)   = yd(id)/vol(id)    
      if(vol(id).GT.0.0) zd(id)   = zd(id)/vol(id)    
      ud(id)   = ud(id)/voldx(id)    
      vd(id)   = vd(id)/voldy(id)    
      wd(id)   = wd(id)/voldz(id)   
      vold(id) = vold(id)*dl(1)*dl(2)*dl(3)
      voldx(id)= voldx(id)*dl(1)*dl(2)*dl(3)
      voldy(id)= voldy(id)*dl(1)*dl(2)*dl(3)
      voldz(id)= voldz(id)*dl(1)*dl(2)*dl(3)
      vol(id)  = vol(id)*dl(1)*dl(2)*dl(3)
    enddo 
    !
    ! Writing results
    !call MPI_Comm_set_errhandler(MPI_COMM_WORLD, MPI_ERRORS_RETURN, ierr)
    !
    call write_output_r8(istep,datadir_ta,'xpos',xd  , nd)
    call write_output_r8(istep,datadir_ta,'ypos',yd  , nd)
    call write_output_r8(istep,datadir_ta,'zpos',zd  , nd)
    !
    call write_output_r8(istep,datadir_ta,'volx',voldx, nd)
    call write_output_r8(istep,datadir_ta,'voly',voldy, nd)
    call write_output_r8(istep,datadir_ta,'volz',voldz, nd)
    call write_output_r8(istep,datadir_ta,'vols',vol, nd)
    !
    call write_output_r8(istep,datadir_ta,'xmin',xmin, nd)
    call write_output_r8(istep,datadir_ta,'xmax',xmax, nd)
    call write_output_r8(istep,datadir_ta,'ymin',ymin, nd)
    call write_output_r8(istep,datadir_ta,'ymax',ymax, nd)
    call write_output_r8(istep,datadir_ta,'zmin',zmin, nd)
    call write_output_r8(istep,datadir_ta,'zmax',zmax, nd)
    call write_output_r8(istep,datadir_ta,'pfmx',psimax, nd)
    !
    call write_output_r8(istep,datadir_ta,'mtot',mtot, nd)
    call write_output_r8(istep,datadir_ta,'mmin',mmin, nd)
    call write_output_r8(istep,datadir_ta,'mmax',mmax, nd)
    !
    call write_output_r8(istep,datadir_ta,'surf',sur, nd)
    call write_output_r8(istep,datadir_ta,'conl',conl, nd)
    call write_output_r8(istep,datadir_ta,'drya',drya, nd)
    !
    call write_output_r8(istep,datadir_ta,'uvel',ud  , nd)
    call write_output_r8(istep,datadir_ta,'vvel',vd  , nd)
    call write_output_r8(istep,datadir_ta,'wvel',wd  , nd)
    call write_output_r8(istep,datadir_ta,'vold',vold, nd)
    !
    call write_output_i1(istep,datadir_ta,'drid',drid, nd)
    call write_output_i1(istep,datadir_ta,'idmy',myid_out, nd)
    call write_output_i6(istep,datadir_ta,'proc',transpose(procShare), nd)
    call write_output_i6(istep,datadir_ta,'orde',transpose(orderField), nd)
    call write_output_i6(istep,datadir_ta,'idsh',transpose(idShare), nd)
    !
    !call MPI_Comm_set_errhandler(MPI_COMM_WORLD, MPI_ERRORS_ARE_FATAL, ierr)
    !
    ! deallocate the droplet variables
    !
    deallocate(xd,yd,zd,ud,vd,wd,vold,procShare,idShare)
    !
    return
  end subroutine tagging_acdi_ibm
  !
  recursive subroutine recursive_tagging(n,i,j,k,id,id_loc,vof_f,vofTag,vofTag_loc)
    !
    use mod_param, only: IBM_h
    !
    implicit none
    !
    integer   , intent(in   ), dimension(3)        :: n
    integer   , intent(in   )                      :: i,j,k
    integer(4), intent(inout)                      :: id,id_loc
    real(rp)  , intent(inout), dimension(0:,0:,0:) :: vof_f
    integer(8), intent(inout), dimension(0:,0:,0:) :: vofTag,vofTag_loc
    !
    integer :: ii,jj,kk
    ! 
    do kk=k-1,k+1
      do jj=j-1,j+1
        do ii=i-1,i+1
          !
          if((vof_f(ii,jj,kk).gt.vof_th).and.(ii.ge.1   ).and.(jj.ge.1   ).and.(kk.ge.1+IBM_h).and. &
                                             (ii.le.n(1)).and.(jj.le.n(2)).and.(kk.le.n(3)   )) then
            vof_f(ii,jj,kk)      = 0._rp
            vofTag(ii,jj,kk)     = id
            vofTag_loc(ii,jj,kk) = id_loc
            ! print*, ii,jj,kk, vof_f(ii,jj,kk), vofTag(ii,jj,kk)
            call recursive_tagging(n,ii,jj,kk,id,id_loc,vof_f,vofTag,vofTag_loc)
            !
          endif
          !
        enddo
      enddo
    enddo
    !
    return
  end subroutine recursive_tagging
  !
  recursive subroutine recursive_remove_tag(n,i,j,k,id_loc,vof_f,vofTag_loc,search_tag)
    !
    use mod_param, only: IBM_h
    !
    implicit none
    !
    integer   , intent(in   ), dimension(3)        :: n
    integer   , intent(in   )                      :: i,j,k
    integer(4), intent(inout)                      :: id_loc
    real(rp)  , intent(inout), dimension(0:,0:,0:) :: vof_f
    integer(8), intent(inout), dimension(0:,0:,0:) :: vofTag_loc
    integer(8), intent(in)                         :: search_tag
    !
    integer :: ii,jj,kk
    ! 
    do kk=k-1,k+1
      do jj=j-1,j+1
        do ii=i-1,i+1
          !
          if((int(vof_f(ii,jj,kk),8).eq.search_tag).and.(ii.ge.1   ).and.(jj.ge.1   ).and.(kk.ge.1).and. &
                                             (ii.le.n(1)).and.(jj.le.n(2)).and.(kk.le.n(3))) then
            vof_f(ii,jj,kk)      = 0._rp
            vofTag_loc(ii,jj,kk) = id_loc
            ! print*, ii,jj,kk, vof_f(ii,jj,kk), vofTag(ii,jj,kk)
            call recursive_remove_tag(n,ii,jj,kk,id_loc,vof_f,vofTag_loc,search_tag)
            !
          endif
          !
        enddo
      enddo
    enddo
    !
    return
  end subroutine recursive_remove_tag
  !
  !subroutine write_output_r8(istep,datadir_ta,varname,var,nd)
  !  !
  !  implicit none
!    !
!    integer,          intent(in)                :: istep
!    character(len=*), intent(in)                :: datadir_ta
!    character(len=*), intent(in)                :: varname
!    real(rp),         intent(in), dimension(1:) :: var
!    integer(4),                     intent(in)  :: nd
!    !
!    integer(kind=MPI_OFFSET_KIND) :: filesize,disp,nreals,nreals_myid
!    integer :: offset,fh,SIZE
!    integer(4) :: sum
!    character(9) :: fldnum
!    character(len=1024) :: fname
    !
 !   write(fldnum,'(i9.9)') istep
  !  disp = 0_MPI_OFFSET_KIND
  !  fname = trim(datadir_ta)//varname//'fld'//fldnum//'.bin'
    !
 !   sum = nd
    !print*,nd
    !call MPI_TYPE_SIZE(MPI_INTEGER, SIZE, ierr)
    !if(myid.eq.0) print*, 'MPI_INTEGER SIZE = ', SIZE
    !call mpi_allreduce(MPI_IN_PLACE ,sum,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD,ierr)
!    call MPI_EXSCAN(nd ,sum,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD,ierr)
!    disp = int(sum,kind=MPI_OFFSET_KIND)  !-nd
!    !print*,disp,myid,(sum)  !-nd
!    call MPI_FILE_OPEN(MPI_COMM_WORLD,trim(fname),MPI_MODE_CREATE+MPI_MODE_WRONLY,MPI_INFO_NULL,fh,ierr)
  !  filesize = 0_MPI_OFFSET_KIND
!    call MPI_FILE_SET_SIZE(fh,filesize,ierr)  ! guarantee overwriting
!    call MPI_FILE_SET_VIEW(fh,disp*8, MPI_DOUBLE_PRECISION,MPI_DOUBLE_PRECISION, 'native', &
!         MPI_INFO_NULL, ierr)
!    call MPI_FILE_WRITE(fh,var,nd,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,ierr)
  !  call MPI_FILE_CLOSE(fh,ierr)
  !  !
!    return
 ! end subroutine write_output_r8
  !
  subroutine write_output_r8(istep,datadir_ta,varname,var,nd)
    !
    implicit none
    !
    integer,          intent(in)                :: istep
    character(len=*), intent(in)                :: datadir_ta
    character(len=*), intent(in)                :: varname
    real(rp),         intent(in), dimension(1:) :: var
    integer,          intent(in)                :: nd
    !
    integer(kind=MPI_OFFSET_KIND) :: filesize,disp,nreals,nreals_myid
    integer :: offset,fh
    character(9) :: fldnum
    character(len=1024) :: fname
    !
    write(fldnum,'(i9.9)') istep
    disp = 0
    fname = trim(datadir_ta)//varname//'fld'//fldnum//'.bin'
    !
    call MPI_SCAN(nd,disp,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
    disp = disp-nd
    call MPI_FILE_OPEN(MPI_COMM_WORLD,trim(fname),MPI_MODE_CREATE+MPI_MODE_WRONLY,MPI_INFO_NULL,fh,ierr)
    filesize = 0_MPI_OFFSET_KIND
    call MPI_FILE_SET_SIZE(fh,filesize,ierr)  ! guarantee overwriting
    call MPI_FILE_SET_VIEW(fh,disp*8, MPI_DOUBLE_PRECISION,MPI_DOUBLE_PRECISION, 'native', &
         MPI_INFO_NULL, ierr)
    call MPI_FILE_WRITE(fh,var,nd,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,ierr)
    call MPI_FILE_CLOSE(fh,ierr)
    !
    return
  end subroutine write_output_r8
  !
  subroutine write_output_i6(istep,datadir_ta,varname,var,nd)
    !
    implicit none
    !
    integer,          intent(in)                   :: istep
    integer,  intent(in)     :: nd
    character(len=*), intent(in)                   :: datadir_ta
    character(len=4), intent(in)                   :: varname
    integer(4)      , intent(in), dimension(1:,1:) :: var
    !
    integer(kind=MPI_OFFSET_KIND) :: filesize,nreals,nreals_myid, displacement,disp
    integer :: offset, fh, ndtot,nelem,sum
    integer, save :: subarray
    character(9) :: fldnum
    character(len=1024) :: fname
    ! 
    nelem = 6
    write(fldnum,'(i9.9)') istep
    disp = 0
    fname = trim(datadir_ta)//varname//'fld'//fldnum//'.bin'
    call MPI_SCAN(nd,sum,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
    disp = sum-nd
    call MPI_FILE_OPEN(MPI_COMM_WORLD,trim(fname),MPI_MODE_CREATE+MPI_MODE_WRONLY,MPI_INFO_NULL,fh,ierr)
    filesize = 0_MPI_OFFSET_KIND
    call MPI_FILE_SET_SIZE(fh,filesize,ierr)    
    call MPI_FILE_SET_VIEW(fh,disp*nelem*4, MPI_INTEGER, MPI_INTEGER, 'native', &
         MPI_INFO_NULL, ierr)
    call MPI_FILE_WRITE_all(fh,var,nelem*nd,MPI_INTEGER,MPI_STATUS_IGNORE,ierr)
    call MPI_FILE_CLOSE(fh,ierr)
    !
    return
  end subroutine write_output_i6
  !
  subroutine write_output_i1(istep,datadir_ta,varname,var,nd)
    !
    implicit none
    !
    integer,          intent(in)                :: istep
    integer                   ,  intent(in)     :: nd
    character(len=*), intent(in)                :: datadir_ta
    character(len=4), intent(in)                :: varname
    integer(4),       intent(in), dimension(1:) :: var
    !
    integer(kind=MPI_OFFSET_KIND) :: filesize,disp,nreals,nreals_myid
    integer :: offset,fh,nelem,sum
    character(9) :: fldnum
    character(len=1024) :: fname
    ! 
    nelem = 1
    write(fldnum,'(i9.9)') istep
    disp = 0
    fname = trim(datadir_ta)//varname//'fld'//fldnum//'.bin'
    call MPI_SCAN(nd,sum,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
    disp = sum-nd
    call MPI_FILE_OPEN(MPI_COMM_WORLD,trim(fname),MPI_MODE_CREATE+MPI_MODE_WRONLY,MPI_INFO_NULL,fh,ierr)
    filesize = 0_MPI_OFFSET_KIND
    call MPI_FILE_SET_SIZE(fh,filesize,ierr)  ! guarantee overwriting
    call MPI_FILE_SET_VIEW(fh,disp*4*nelem, MPI_INTEGER,MPI_INTEGER, 'native', &
         MPI_INFO_NULL, ierr)
    call MPI_FILE_WRITE(fh,var,nd,MPI_INTEGER,MPI_STATUS_IGNORE,ierr)
    call MPI_FILE_CLOSE(fh,ierr)
    !
    return
  end subroutine write_output_i1
  !
end module
