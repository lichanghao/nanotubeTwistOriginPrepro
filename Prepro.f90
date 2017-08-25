program PREPRO
  USE data_mesh
  USE data_BC
  USE data_tensor22
  USE data_mat
  USE data_vdw
  USE data_vector2
  USE data_vector3
  implicit REAL(8) (a-h,o-z)
  implicit INTEGER*4 (i-n)
  !include 'mpif.h' !MPI
  TYPE(mesh), ALLOCATABLE :: mesh0(:)    ! array of meshes
  TYPE(mesh) :: meshT                    ! Total mesh
  TYPE(material) :: mat1
  TYPE(vdw_data), ALLOCATABLE :: vdw1(:) ! Array of vdw data
  TYPE(vdw_data):: vdwT                  ! Total vdw data
  TYPE(BC_data), ALLOCATABLE :: BCs(:)   ! Array of BCs
  TYPE(BC_data):: BCsT                   ! Total BCs
  TYPE(tensor22), ALLOCATABLE :: F0(:), temp(:)
  REAL(8), ALLOCATABLE :: x0(:), J0(:)
  REAL(8), ALLOCATABLE :: shapef(:,:,:), weight(:)
  REAL(8), ALLOCATABLE :: x0_BC(:)
  REAL(8) :: crit(2), E_out(4)
  REAL(8), ALLOCATABLE :: ylength(:)
  INTEGER(4), ALLOCATABLE :: nrow(:), ncol(:), numel(:), numno(:)
  INTEGER(4), ALLOCATABLE :: numed(:), nneigh(:), ninrange(:)
  INTEGER(4), ALLOCATABLE :: ntab(:)
  integer, dimension(:,:), allocatable :: mtable
  REAL(8), ALLOCATABLE :: x0T(:), J0T(:)
  TYPE(tensor22), ALLOCATABLE :: F0T(:)
  REAL(8), PARAMETER :: PI = 3.1415926


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!********************************************************************************************************
  REAL(8) :: shapeff(2,12), weightt(2),xneigh(12,3),neigh(12),neigh1(12),shapee(1,12)
  REAL(8) :: xx(3),r1,r2,thta,xxlength,xii(3),xjj(3),xkk(3)
  REAL(8) :: x1,x2,x3,x4,y1,y2,y3,y4,t1,t2,t3,xpoint,ypoint
  REAL(8), ALLOCATABLE :: x0inner(:,:),indexelem(:),shapefff(:,:),x00(:)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!*****************************************************************************************************************************************


  !call MPI_INIT() !MPI
  !call MPI_COMM_SIZE(MPI_COMM_WORLD,NPROCS,ierr)
  !call MPI_COMM_RANK(MPI_COMM_WORLD,ID,ierr)

  write(*,*) ' Skip make_near ? (1=yes)'
  read(*,*) nskip


  open(unit=90,file='nano_dims.dat',status='unknown')
  open(unit=91,file='nano_general.dat',status='unknown')
  open(unit=92,file='nano_BCs.dat',status='unknown')
  open(unit=94,file='nano_Mesh.dat',status='unknown')
  open(unit=95,file='nano_zero.dat',status='unknown')
  open(unit=96,file='nano_config.dat',status='unknown')
  open(unit=56,file='nano_tub_loc.dat',status='unknown')

  open(unit=33,file='data.dat',status='old')
  read(33,*)
  read(33,*) ntubes
  rewind(33)
  ALLOCATE(mesh0(ntubes),BCs(ntubes),vdw1(ntubes),nrow(ntubes),ncol(ntubes),  &
       ylength(ntubes), nneigh(ntubes), ninrange(ntubes), &
       numel(ntubes),numno(ntubes),numed(ntubes),STAT=istat)
  if (istat/=0) STOP '**** Not enough memory ****'

  call read_data(ntubes,nrow,ncol,xlength,ylength, & 
       numel,numno,numed,mat1,stretch_ini, &
       angle,angle2,nloadstep,crit,ngauss,nW_hat,nCodeLoad, &
       imperfect,fact_imp,nborder,vdwT)
  close(33)


  !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  ! Dimensions of Total model
  meshT%numele=0
  meshT%numnods=0
  meshT%nedge=0
  meshT%nelem_ghost=0
  meshT%nnode_ghost=0

  BCsT%nnodBC=0
  BCsT%ndofBC=0
  BCsT%ndofOP=0
  BCsT%nCodeLoad=nCodeLoad
  BCsT%nloadstep=nloadstep

  ylengthT=0.d0

  nelT=0
  nodT=0
  nedT=0
  do imesh=1,ntubes
     nelT=nelT+numel(imesh)
     nodT=nodT+(nrow(imesh))*(ncol(imesh)+1)
     nedT=nedT+2*(nrow(imesh))
  enddo

  ALLOCATE(x0T(3*(nodT+nedT)), &
       meshT%connect(nelT),J0T(nelT), &
       !        meshT%ntable(nodT,maxneigh_vert+1), &
  meshT%nghost_tab(nedT,3), &
       F0T(nelT),STAT=istat)
  if (istat/=0) STOP '**** Not enough memory ****'

  !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


  if (vdwT%nvdw.eq.1) then
     ALLOCATE(vdwT%shapef(vdwT%ngauss_vdw,12),vdwT%weight(vdwT%ngauss_vdw),STAT=istat)
     if (istat/=0) STOP '**** Not enough memory ****'
     call gauss_vdw(vdwT%ngauss_vdw,vdwT%shapef,vdwT%weight)
  endif

  numel_max=maxval(numel)
  ALLOCATE(temp(numel_max),STAT=istat)
  if (istat/=0) STOP '**** Not enough memory ****'


  ALLOCATE(mtable(maxval(numno),13),STAT=istat)
  if (istat/=0) STOP '**** Not enough memory ****'

  do imesh=1,ntubes
     !   if ((imesh.gt.1).and.(nCodeLoad.eq.6)) nCodeLoad=666

     mesh0(imesh)%numele=numel(imesh)
     mesh0(imesh)%numnods=numno(imesh)
     mesh0(imesh)%nedge=numed(imesh)
     BCs(imesh)%nCodeLoad=nCodeLoad
     BCs(imesh)%nloadstep=nloadstep
     BCs(imesh)%value=angle

     ALLOCATE(x0(3*(numno(imesh)+numed(imesh))), &
          mesh0(imesh)%connect(numel(imesh)),J0(mesh0(imesh)%numele), &
          mesh0(imesh)%ntable(numno(imesh),maxneigh_vert+1), &
          mesh0(imesh)%nghost_tab(numed(imesh),3), &
          F0(mesh0(imesh)%numele),STAT=istat)
     if (istat/=0) STOP '**** Not enough memory ****'
     call mesh_gen_square(nrow(imesh),ncol(imesh),xlength,ylength(imesh),x0,mesh0(imesh))
     call Def_Grad_Cart_Conv(temp,mesh0(imesh),x0,J0,F0)
     write(*,*) imesh, ' deallocate x0 '
     DEALLOCATE(x0)
     write(*,*) imesh, ' deallocate connect '
     write(*,*) numel(imesh)

     !   DEALLOCATE(mesh0(imesh)%connect)
     write(*,*) imesh, ' deallocate ntable '
     !   DEALLOCATE(mesh0(imesh)%ntable)
     write(*,*) imesh, ' deallocate nghost_tab '
     !   DEALLOCATE(mesh0(imesh)%nghost_tab)
     write(*,*) imesh, ' deallocate temp '
     !   DEALLOCATE(temp)
!!$   DEALLOCATE(x0,mesh0(imesh)%connect,mesh0(imesh)%ntable,mesh0(imesh)%nghost_tab,temp)
     write(*,*) imesh, ' deallocate end '
     !@@@@@@@@@@ END UNDEFORMED MESH AND RELEVANT DATA @@@@@@@@@@@@@



     !@@@@@@@@@@ Pre-Processing of ghost nodes and elements @@@@@@@@@@@@@
     if (nborder.gt.0) then
        mesh0(imesh)%nnode_ghost=2*nborder*nrow(imesh)
        mesh0(imesh)%nelem_ghost=nborder*4*nrow(imesh)
        ALLOCATE(mesh0(imesh)%elem_ghost(mesh0(imesh)%nelem_ghost), &
             mesh0(imesh)%node_ghost(mesh0(imesh)%nnode_ghost),STAT=istat)
        if (istat/=0) STOP '**** Not enough memory ****'
        ii=0
        jj=0
        do irow=1,nrow(imesh)
           !   do icol=1+nborder,ncol(imesh)-nborder
           do icol=1,nborder
              ielem=((irow-1)*ncol(imesh)+icol)*2 - 1
              ii=ii+1
              mesh0(imesh)%elem_ghost(ii)=ielem
              ii=ii+1
              mesh0(imesh)%elem_ghost(ii)=ielem+1
              jj=jj+1
              inod=(irow-1)*(ncol(imesh)+1)+icol
              mesh0(imesh)%node_ghost(jj)=inod
           enddo
           do icol=ncol(imesh)-nborder+1,ncol(imesh)
              ielem=((irow-1)*ncol(imesh)+icol)*2 - 1
              ii=ii+1
              mesh0(imesh)%elem_ghost(ii)=ielem
              ii=ii+1
              mesh0(imesh)%elem_ghost(ii)=ielem+1
              jj=jj+1
              inod=(irow-1)*(ncol(imesh)+1)+icol+1
              mesh0(imesh)%node_ghost(jj)=inod
           enddo
        end do
     else
        mesh0(imesh)%nnode_ghost=0
        mesh0(imesh)%nelem_ghost=0
        ALLOCATE(mesh0(imesh)%elem_ghost(1),mesh0(imesh)%node_ghost(1),STAT=istat)
        if (istat/=0) STOP '**** Not enough memory ****'
        mesh0(imesh)%elem_ghost(1)=-10
        mesh0(imesh)%node_ghost(1)=-10
     end if
     !@@@@@@@@@@ END Pre-Processing of ghost nodes and elements @@@@@@@@@@@@@

     !@@@@@@@@@@ CYLINDER MESH @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
     ! Correct data for cylinder
     numno(imesh)=(nrow(imesh))*(ncol(imesh)+1)
     numed(imesh)=2*(nrow(imesh))
     !dimensions
     mesh0(imesh)%numele=numel(imesh)
     mesh0(imesh)%numnods=numno(imesh)
     mesh0(imesh)%nedge=numed(imesh)
     write(*,*) '**** New Mesh ****'
     write(*,*) 'numno', numno(imesh), '  numel', numel(imesh), '  numed', numed(imesh)
     ! Allocate memory
     write(*,*) '**** Allocate Memory ****'
     ALLOCATE(x0(3*(numno(imesh)+numed(imesh))), &
          mesh0(imesh)%connect(numel(imesh)), &
          mesh0(imesh)%ntable(numno(imesh),maxneigh_vert+1), &
          mesh0(imesh)%nghost_tab(numed(imesh),3),STAT=istat)
     if (istat/=0) STOP '**** Not enough memory ****1'
     ! Generate the basic mesh, coords, connect and tag on boun. edges
     write(*,*) '**** Generate Mesh ****'
     call mesh_gen_cylin(nrow(imesh),ncol(imesh), &
          xlength,ylength(imesh),x0,mesh0(imesh),stretch_ini)

     write(*,*) '**** connect the mesh ****'
     call connect_mesh(mesh0(imesh),mtable)     
     !@@@@@@@@@@ END CYLINDER MESH @@@@@@@@@@@@@@@@@@@@@@@

     ! arrange the tube
!     if(imesh.ge.2)then
!        rtube = ylength(imesh)/PI + 0.34
!        do inode = 1, numno(imesh)+numed(imesh)
!           ath = (imesh-2)*PI/3.d0
!           x0(3*inode-1) = x0(3*inode-1) + rtube*cos(ath)
!           x0(3*inode)   = x0(3*inode)   + rtube*sin(ath)
!        enddo
!     endif
     

     !@@@@@@@@@@@ INITIALIZATIONS @@@@@@@@@@@@@@@@@@@@@@@@
     ! Van der Waals initializations
     if (vdwT%nvdw.eq.1) then
        write(*,*) ' begin van der Waals initializations '
        ALLOCATE(vdw1(imesh)%shapef(vdwT%ngauss_vdw,12), &
             vdw1(imesh)%weight(vdwT%ngauss_vdw),STAT=istat)
        if (istat/=0) STOP '**** Not enough memory ****3'


        vdw1(imesh)%nvdw=vdwT%nvdw
        vdw1(imesh)%ngauss_vdw=vdwT%ngauss_vdw
        vdw1(imesh)%r_cut=vdwT%r_cut
        vdw1(imesh)%r_bond=vdwT%r_bond
        vdw1(imesh)%a=vdwT%a
        vdw1(imesh)%sig=vdwT%sig
        vdw1(imesh)%y0=vdwT%y0
        vdw1(imesh)%meval=vdwT%meval
        vdw1(imesh)%Vcut=vdwT%Vcut

        call gauss_vdw(vdw1(imesh)%ngauss_vdw,vdw1(imesh)%shapef,vdw1(imesh)%weight)

        vdw1(imesh)%ng_tot=vdwT%ngauss_vdw*mesh0(imesh)%numele
        r_nb=vdwT%r_bond ! this is the distance over which it is bonded
        ! estimate the number of bonded gauss points and gaus points of far away
        ! tubes in the nested case.
        !      nneigh(imesh)=2*int(3.1416*r_nb*r_nb/J0(imesh)*vdwT%ngauss_vdw*2) 
!        nkk=int(3.1416*r_nb*r_nb/(J0(1)/2.)*vdwT%ngauss_vdw*1.3)
        nkk = 0
        if (nskip.eq.1) then
           !nneigh(imesh)=14*vdwT%ngauss_vdw
           nneigh(imesh)=-1
        else
           nneigh(imesh)=nkk
        endif
        ! estim the # of in-range vdw gp
        !      ninrange(imesh)=int(3.1416*(vdwT%r_cut)**2/J0(imesh)*2*vdwT%ngauss_vdw*4) 
        !      ninrange(imesh)=int(3.1416*(vdwT%r_cut)**2/(J0(1)/2.)*vdwT%ngauss_vdw * 2 * 1.4)
!        ninrange(imesh)=int(3.1416*(vdwT%r_cut)**2/(J0(1)/2.)* &
!             vdwT%ngauss_vdw * 3.5-nkk/1.3)

        ninrange(imesh)=int(3.1416*(vdwT%r_cut**2-0.33*0.33)/(J0(1)/2.)* &
             vdwT%ngauss_vdw * 3.5 )*3

        if (nskip.eq.1) then
           continue
        else
           ALLOCATE(vdw1(imesh)%near(vdw1(imesh)%ng_tot,0:nneigh(imesh)), &
                vdw1(imesh)%x(vdw1(imesh)%ng_tot,3),STAT=istat)
           if (istat/=0) STOP '**** Not enough memory ****4'
           call make_near(vdw1(imesh),mesh0(imesh),x0,vdw1(imesh)%r_bond,nskip,nneigh(imesh))
        endif
        write(*,*) ' end van der Waals initializations '
     endif
     !@@@@@@@@@@@ END INITIALIZATIONS @@@@@@@@@@@@@@@@@@@@@@@@



     !@@@@@@@@@@@ PREPROCESSING OF LOADING OF NANOTUBE @@@@@@@@@@@@@@@@@@@
     if (BCs(imesh)%nCodeLoad.eq.0) then
        if (imesh.eq.1) then
           BCs(imesh)%ndofBC=nrow(imesh)+3
           BCs(imesh)%nnodBC=nrow(imesh)
        else
           BCs(imesh)%ndofBC=nrow(imesh)
           BCs(imesh)%nnodBC=nrow(imesh)
        endif
     else if ((BCs(imesh)%nCodeLoad.eq.1).or. &
          (BCs(imesh)%nCodeLoad.eq.2).or. &
          (BCs(imesh)%nCodeLoad.eq.3).or.(BCs(imesh)%nCodeLoad.eq.13)) then
        !if(imesh.eq.1) then
           BCs(imesh)%ndofBC=2*3*nrow(imesh)
           BCs(imesh)%nnodBC=2*nrow(imesh)
        !else
        !   BCs(imesh)%ndofBC=0
        !   BCs(imesh)%nnodBC=0
        !endif
     else if (BCs(imesh)%nCodeLoad.eq.10) then
        BCs(imesh)%ndofBC=2*3*nrow(imesh) + 2*(ncol(imesh)-1)
        BCs(imesh)%nnodBC=2*nrow(imesh) + 2*(ncol(imesh)-1)
        !      BCs(imesh)%nnodBC=1
     else if (BCs(imesh)%nCodeLoad.eq.4) then
        BCs(imesh)%ndofBC=2*nrow(imesh)+3
        BCs(imesh)%nnodBC=2*nrow(imesh)
     else if (BCs(imesh)%nCodeLoad.eq.5) then
        BCs(imesh)%ndofBC=3*nrow(imesh)
        BCs(imesh)%nnodBC=2*nrow(imesh)
     else if ((BCs(imesh)%nCodeLoad.eq.6).and.(imesh.eq.1)) then
        BCs(imesh)%ndofBC=4*(ncol(imesh)+1) + nrow(imesh)
        BCs(imesh)%nnodBC=2*(ncol(imesh)+1) + nrow(imesh) - 2

        write(*,*) 'dof', BCs(imesh)%ndofBC,' nod', BCs(imesh)%nnodBC

     else if ((BCs(imesh)%nCodeLoad.eq.6).and.(imesh.ne.1)) then
        BCs(imesh)%ndofBC=2*(ncol(imesh)+1) + nrow(imesh)
        BCs(imesh)%nnodBC=2*(ncol(imesh)+1) + nrow(imesh) - 2
        !BCs(imesh)%ndofBC=0
        !BCs(imesh)%nnodBC=0
     else if (BCs(imesh)%nCodeLoad.eq.666) then
        BCs(imesh)%ndofBC=2*(ncol(imesh)+1)+2
        BCs(imesh)%nnodBC=2*(ncol(imesh)+1)
     else if (BCs(imesh)%nCodeLoad.eq.7) then
        BCs(imesh)%ndofBC=12
        BCs(imesh)%nnodBC=4
     else if (BCs(imesh)%nCodeLoad.eq.8) then
        BCs(imesh)%ndofBC= 2*(ncol(imesh)-1) + 2 + 2
        BCs(imesh)%nnodBC=2*nrow(imesh)   + 2*(ncol(imesh)-1) + 2
     else
        STOP ' Loading option not implemented1'
     end if
     BCs(imesh)%ndofOP=3*mesh0(imesh)%numnods-BCs(imesh)%ndofBC
     ALLOCATE(BCs(imesh)%mdofBC(BCs(imesh)%ndofBC), &
          BCs(imesh)%mdofOP(BCs(imesh)%ndofOP), &
          BCs(imesh)%mnodBC(BCs(imesh)%nnodBC,2),STAT=istat)
     if (istat/=0) STOP '**** Not enough memory ****5'

     write(*,*) 'We enter load pre'
     call load_pre(x0,mesh0(imesh),BCs(imesh),xlength,ylength(imesh),& 
          nrow(imesh),ncol(imesh),nborder,imesh,angle2)
     !@@@@@@@@@@@ END PREPROCESSING OF LOADING OF NANOTUBE @@@@@@@@@@@@@@@@@@@



     !@@@@@@@@@@@@ Increments and assemblies @@@@@@@@@@@@@@@@@@@@ 
     meshT%numele=meshT%numele+numel(imesh)
     meshT%numnods=meshT%numnods+numno(imesh)
     meshT%nedge=meshT%nedge+numed(imesh)
     meshT%nelem_ghost=meshT%nelem_ghost+mesh0(imesh)%nelem_ghost
     meshT%nnode_ghost=meshT%nnode_ghost+mesh0(imesh)%nnode_ghost

     BCsT%nnodBC=BCsT%nnodBC + BCs(imesh)%nnodBC
     BCsT%ndofBC=BCsT%ndofBC + BCs(imesh)%ndofBC
     BCsT%ndofOP=BCsT%ndofOP + BCs(imesh)%ndofOP

     vdwT%ng_tot=vdwT%ng_tot+vdw1(imesh)%ng_tot

     ylengthT=ylengthT+ylength(imesh)


     istart=meshT%numnods-numno(imesh)+1
     iend  =meshT%numnods 
     x0T(3*istart-2:3*iend)=x0(1:3*numno(imesh))
     istart=meshT%numele-numel(imesh)+1
     iend  =meshT%numele
     J0T(istart:iend)=J0(1:numel(imesh))
     F0T(istart:iend)=F0(1:numel(imesh))

     !%%%% The mesh!!!
     do ielem=1,numel(imesh)
        ielT= meshT%numele-numel(imesh) +  ielem
        ! Connect
        meshT%connect(ielT)=mesh0(imesh)%connect(ielem)
        do i=1,maxneigh_vert
           ikk=meshT%connect(ielT)%neigh_vert(i)-numno(imesh)
           if (ikk.gt.0) then
              meshT%connect(ielT)%neigh_vert(i)=nodT+meshT%nedge-numed(imesh)+ikk
           else
              meshT%connect(ielT)%neigh_vert(i)=meshT%connect(ielT)%neigh_vert(i) &
                   +meshT%numnods-numno(imesh)
           endif
        enddo
        meshT%connect(ielT)%vertices=meshT%connect(ielT)%vertices+meshT%numnods-numno(imesh)
        meshT%connect(ielT)%neigh_elem=meshT%connect(ielT)%neigh_elem+meshT%numele-numel(imesh)
     enddo
     ! Ghost_tab
     do ied=1,numed(imesh)
        iedT= meshT%nedge-numed(imesh) + ied
        meshT%nghost_tab(iedT,:)=mesh0(imesh)%nghost_tab(ied,:)  &
             +meshT%numnods-numno(imesh)
     enddo
     !%%%% The mesh!!!

     !@@@@@@@@@@@@ Increments and assemblies @@@@@@@@@@@@@@@@@@@@ 
     write(*,*) 'We exit load pre3'
     DEALLOCATE(x0,F0,J0)
     write(*,*) 'We exit load pre4'
  enddo


  ALLOCATE(ntab(meshT%numnods),STAT=istat)
  if (istat/=0) STOP '**** Not enough memory ****'


  if (nborder.gt.0) then
     ALLOCATE(meshT%elem_ghost(meshT%nelem_ghost), &
          meshT%node_ghost(meshT%nnode_ghost),STAT=istat)
     n1=0
     n2=0
     n3=0
     n4=0
     do imesh=1,ntubes
        meshT%elem_ghost(n1+1:n1+mesh0(imesh)%nelem_ghost)=mesh0(imesh)%elem_ghost(:)+n3
        meshT%node_ghost(n2+1:n2+mesh0(imesh)%nnode_ghost)=mesh0(imesh)%node_ghost(:)+n4
        n1=n1 +mesh0(imesh)%nelem_ghost
        n2=n2 +mesh0(imesh)%nnode_ghost
        n3=n3 +mesh0(imesh)%numele
        n4=n4 +mesh0(imesh)%numnods
     enddo
  else
     meshT%nnode_ghost=0
     meshT%nelem_ghost=0
     ALLOCATE(meshT%elem_ghost(1),meshT%node_ghost(1),STAT=istat)
     if (istat/=0) STOP '**** Not enough memory ****6'
     meshT%elem_ghost(1)=-10
     meshT%node_ghost(1)=-10
  endif

  ALLOCATE(BCsT%mdofBC(BCsT%ndofBC), &
       BCsT%mdofOP(BCsT%ndofOP), &
       BCsT%mnodBC(BCsT%nnodBC,2),STAT=istat)
  if (istat/=0) STOP '**** Not enough memory ****7'
  nneighT=maxval(nneigh(1:ntubes))
  if (nskip.eq.1) then
     ALLOCATE(vdwT%near(1,0:1),STAT=istat)
  else
     ALLOCATE(vdwT%near(vdwT%ng_tot,0:nneighT),STAT=istat)
  endif
  vdwT%near=0
  BCsT%value=angle

  if (BCsT%nCodeLoad.eq.0) then
     BCsT%value=0.
     BCsT%nloadstep=1
  endif
  BCsT%xc=BCs(1)%xc
  BCsT%rotation=BCs(1)%rotation
  n1=0
  n2=0
  n3=0
  n4=0
  n5=0
  do imesh=1,ntubes
     BCsT%mnodBC(n1+1:n1+BCs(imesh)%nnodBC,1)=BCs(imesh)%mnodBC(:,1)+n4
     BCsT%mnodBC(n1+1:n1+BCs(imesh)%nnodBC,2)=BCs(imesh)%mnodBC(:,2)
     BCsT%mdofBC(n2+1:n2+BCs(imesh)%ndofBC)=BCs(imesh)%mdofBC(:)+3*n4
     BCsT%mdofOP(n3+1:n3+BCs(imesh)%ndofOP)=BCs(imesh)%mdofOP(:)+3*n4
     if  (vdwT%nvdw.eq.1) then
        istart=vdwT%ngauss_vdw*n5+1
        iend  =vdwT%ngauss_vdw*(n5+mesh0(imesh)%numele)
!!$ vdwT%near(istart:iend,0)=vdw1(imesh)%near(:,0)
!!$ vdwT%near(istart:iend,1:nneigh)=vdw1(imesh)%near(:,1:nneigh)
        if (nskip.eq.1) then
           continue
        else
           do ii=istart,iend
              vdwT%near(ii,0)=vdw1(imesh)%near(ii-vdwT%ngauss_vdw*n5,0)
              vdwT%near(ii,1:vdwT%near(ii,0))= &
                   vdw1(imesh)%near(ii-vdwT%ngauss_vdw*n5,1:vdwT%near(ii,0))+vdwT%ngauss_vdw*n5
           enddo
        endif
     endif
     n1=n1 + BCs(imesh)%nnodBC
     n2=n2 + BCs(imesh)%ndofBC
     n3=n3 + BCs(imesh)%ndofOP
     n4=n4 +mesh0(imesh)%numnods
     n5=n5 +mesh0(imesh)%numele
  enddo

  ninrangeT=maxval(ninrange(:))

  !%%%%%%%%%%%%%%% HERE ARE OUTPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !%%%%%%%%%%%%%%% HERE ARE OUTPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !%%%%%%%%%%%%%%% HERE ARE OUTPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !%%%%%%%%%%%%%%% HERE ARE OUTPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!interlayer bond!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !***********************************************************************************************************************
  !***********************************************************************************************************************
  !***********************************************************************************************************************
  
!!$  open(unit=31,file='bondadd.dat',status='new')
!!$
!!$
!!$  call gauss_vdw(2,shapeff,weightt)
!!$
!  do i=1,2
!          write(*,*)"hereeeeeeeeeee",i
!        do j=1,12
!              write(*,*)shapeff(i,j)
!         enddo
!  enddo

!!$
!!$
!!$  ngausstotal= meshT%numele*2
!!$  ngaussfind =(meshT%numele-mesh0(ntubes)%numele)*2
!!$  ALLOCATE(x0inner(ngausstotal,3),indexelem(ngaussfind),shapefff(12,ngaussfind),STAT=istat)
!!$  ALLOCATE(x00(3*(nodT+nedT)),STAT=istat)
!!$  x0inner=0.0
!!$  x00=0.0
!!$
!!$ 
!!$!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$  open(unit=30,file='x0t.dat',status='old')
!!$  !write(*,*)nodT+nedT
!!$  do i=1,nodT+nedT
!!$     ii=i*3
!!$     read(30,*) x00(ii-2),x00(ii-1), x00(ii)
!!$  enddo
!!$!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$  !     write(*,*)x0T
!!$  !       stop
!!$
!!$  do ielem=1,meshT%numele
!!$ !    neigh=meshT%connect(ielem)%neigh_vert
!!$
!!$     if(mod(ielem,2).eq.1)then
!!$    neigh=meshT%connect(ielem)%neigh_vert
!!$     else
!!$       neigh1=meshT%connect(ielem)%neigh_vert
!!$       neigh(1)=  neigh1(6) 
!!$       neigh(2)=  neigh1(10)
!!$       neigh(3)=  neigh1(3)
!!$       neigh(4)=  neigh1(7)
!!$       neigh(5)=  neigh1(11)
!!$       neigh(6)=  neigh1(1)
!!$       neigh(7)=  neigh1(4)
!!$       neigh(8)=  neigh1(8)
!!$       neigh(9)=  neigh1(12)
!!$       neigh(10)= neigh1(2)
!!$       neigh(11)= neigh1(5)
!!$       neigh(12)= neigh1(9)
!!$    endif
!!$
!!$
!!$
!!$
!!$     xneigh(:,1)=x00(3*neigh(:)-2)
!!$     xneigh(:,2)=x00(3*neigh(:)-1)
!!$     xneigh(:,3)=x00(3*neigh(:))
!!$
!!$
!!$
!!$
!!$

!     do i=1,12
!        do j=1,3
!           write(*,*) xneigh(i,j)
!        enddo
!     enddo
!     stop
!!$
!!$
!!$     do ig=1,2
!!$        igg=(ielem-1)*2+ig
!!$        do inod=1,12
!!$           do idof=1,3
!!$              x0inner(igg,idof)= x0inner(igg,idof) + shapeff(ig,inod)*xneigh(inod,idof)
!!$           enddo
!!$!                  write(*,*) inod,mod(inod+12,12)
!!$        enddo
!!$        !             aaa= sqrt(x0inner(igg,2)**2+x0inner(igg,3)**2)
!!$        !   write(*,*)x0inner(igg,1),x0inner(igg,2),x0inner(igg,3
!!$
!!$        !      stop
!!$     enddo
!!$  enddo
!!$
!!$  !      write(*,*)"xxxxxxxxxxxxxxxxx",x0inner(1,1)
!!$
!!$
!!$!  write(*,*)x0inner(3,1),x0inner(3,2),x0inner(3,3)
!!$!  write(*,*)x0inner(4,1),x0inner(4,2),x0inner(4,3)
!!$
!!$
!!$  do ielem=1,meshT%numele-mesh0(ntubes)%numele
!!$     nele=0
!!$
!!$     do i=1,ntubes-1
!!$        if ((ielem.gt.nele).and.(ielem.le.(mesh0(i)%numele+nele)))then
!!$           indextube=i
!!$
!!$           nele=nele+mesh0(i)%numele
!!$
!!$           exit
!!$        endif
!!$
!!$        nele=mesh0(i)%numele+nele
!!$     enddo
!!$
!!$     do ig=1,2
!!$        igg=(ielem-1)*2+ig
!!$        r1=sqrt(x0inner(igg,2)**2+x0inner(igg,3)**2)
!!$        r2=r1+0.34
!!$        xx(1)=x0inner(igg,1)
!!$        xx(2)=x0inner(igg,2)
!!$        xx(3)=x0inner(igg,3)
!!$        !        write(*,*)xx(1),xx(2),xx(3)
!!$        thta=acos(xx(2)/r1)
!!$        if(xx(3).lt.0.0)then
!!$           thta=3.1415926*2-thta
!!$        endif
!!$        !         write(*,*)thta
!!$        !        write(*,*)'radiusssssssssssssssss',r1
!!$        !        x4=(r1+0.37)*thta
!!$        x4=ylength(indextube+1)*thta/(3.1415926*2)
!!$
!!$        y4=xx(1)
 !       if(igg.lt.5)then
 !          write(*,*)x0inner(igg,2)
 !          write(*,*)x0inner(igg,3)
 !          !                                write(*,*)xpoint
 !       endif
!!$        !           ccc=0
!!$        do iiele=nele+1,mesh0(indextube+1)%numele+nele
!!$
!!$
!!$           ii= meshT%connect(iiele)%vertices(1)
!!$           jj= meshT%connect(iiele)%vertices(2)
!!$           kk= meshT%connect(iiele)%vertices(3) 
!!$
!!$           xii(1)=x00(3*ii-2)
!!$           xii(2)=x00(3*ii-1)
!!$           xii(3)=x00(3*ii)
!!$           r1=sqrt(xii(2)**2+xii(3)**2)
!!$           thta=acos(xii(2)/r1)
!!$           if(xii(3).lt.0.0)then
!!$              thta=3.1415926*2-thta
!!$           endif
!!$           if(xii(3).eq.0.0.and.iiele.gt.(nele+1+mesh0(indextube+1)%numele+nele)/2.0)then
!!$              thta=3.1415926*2
!!$           endif
!!$           x1=r1*thta
!!$           y1=xii(1)
!!$
!!$           xjj(1)=x00(3*jj-2)
!!$           xjj(2)=x00(3*jj-1)
!!$           xjj(3)=x00(3*jj)
!!$           r1=sqrt(xjj(2)**2+xjj(3)**2)
!!$           thta=acos(xjj(2)/r1)
!!$           if(xjj(3).lt.0.0)then
!!$              thta=3.1415926*2-thta
!!$           endif
!!$           if(xjj(3).eq.0.0.and.iiele.gt.(nele+1+mesh0(indextube+1)%numele+nele)/2.0)then
!!$              thta=3.1415926*2
!!$           endif
!!$           x2=r1*thta
!!$           y2=xjj(1)
!!$
!!$           xkk(1)=x00(3*kk-2)
!!$           xkk(2)=x00(3*kk-1)
!!$           xkk(3)=x00(3*kk)
!!$           r1=sqrt(xkk(2)**2+xkk(3)**2)
!!$           thta=acos(xkk(2)/r1)
!!$           if(xkk(3).lt.0.0)then
!!$              thta=3.1415926*2-thta
!!$           endif
!!$           if(xkk(3).eq.0.0.and.iiele.gt.(nele+1+mesh0(indextube+1)%numele+nele)/2.0)then
!!$              thta=3.1415926*2
!!$           endif
!!$           x3=r1*thta
!!$           y3=xkk(1)
!!$
!!$
!!$
!!$
!!$           t1=(y4-y1)*(x2-x1)-(x4-x1)*(y2-y1)
!!$           t2=(y4-y2)*(x3-x2)-(x4-x2)*(y3-y2)
!!$           t3=(y4-y3)*(x1-x3)-(x4-x3)*(y1-y3)
!!$
!!$
!!$           !           if(iiele.eq.241)then
!!$           !           write(*,*)x1,x2,x3,x4
!!$           !               write(*,*)y1,y2,y3,y4
!!$           !               write(*,*)t1,t2,t3,t4
!!$           !           endif
!!$
!!$           if((t1.lt.0).and.((t2.lt.0).and.(t3.lt.0)))then
!!$              indexelem(igg)=iiele
!!$              !              ccc=ccc+1
!!$              if(x1.eq.x2)then
!!$                 xpoint=abs((x4-x1)/(x3-x1))
!!$              elseif(x1.eq.x3)then
!!$                 xpoint=abs((x4-x1)/(x2-x1))
!!$              elseif(x2.eq.x3)then
!!$                 xpoint=abs((x4-x2)/(x2-x1))
!!$              else
!!$                 write(*,*)"wrong! no two the same"
!!$              endif
!!$
!!$              if(y1.eq.y2)then
!!$                 ypoint=abs((y4-y1)/(y3-y1))
!!$              elseif(y1.eq.y3)then
!!$                 ypoint=abs((y4-y1)/(y2-y1))
!!$              elseif(y2.eq.y3)then
!!$                 ypoint=abs((y4-y2)/(y2-y1))
!!$              else
!!$                 write(*,*)"wrong! no two the same"
!!$              endif
!!$!                           write(*,*)xpoint,ypoint
!!$              !            stop
!!$
!!$!              if(igg.lt.5)then
!!$!                 write(*,*)x1,x2,x3,x4
!!$!                      write(*,*)Y1,Y2,Y3,Y4
!!$!                                           write(*,*)xpoint,ypoint
!!$!                                          write(*,*)ypoint
!!$!              endif
!!$
!!$              if(mod(iiele,2).eq.0)then
!!$                 change=xpoint
!!$                 xpoint=ypoint 
!!$                 ypoint=change
!!$              endif
!!$
!!$              call gauss_point(1,shapee,weighttt,ypoint,xpoint)
!!$
!!$
!!$              
!!$!             WRITE(*,*)shapee(1,1:12)
!!$                  
!!$
!!$              do indexx=1,12
!!$                 shapefff(indexx,igg)=shapee(1,indexx)
!!$              enddo
!!$!                WRITE(*,*)shapefff(1:12,igg)
!!$
!!$  
!!$           else if((t1.gt.0).and.((t2.gt.0).and.(t3.gt.0)))then
!!$              indexelem(igg)=iiele
!!$              !                ccc=ccc+1
!!$              if(x1.eq.x2)then
!!$                 xpoint=abs((x4-x1)/(x3-x1))
!!$              elseif(x1.eq.x3)then
!!$                 xpoint=abs((x4-x1)/(x2-x1))
!!$              elseif(x2.eq.x3)then
!!$                 xpoint=abs((x4-x2)/(x2-x1))
!!$              else
!!$                 write(*,*)"wrong! no two the same"
!!$              endif
!!$
!!$              if(y1.eq.y2)then
!!$                 ypoint=abs((y4-y1)/(y3-y1))
!!$              elseif(y1.eq.y3)then
!!$                 ypoint=abs((y4-y1)/(y2-y1))
!!$              elseif(y2.eq.y3)then
!!$                 ypoint=abs((y4-y2)/(y2-y1))
!!$              else
!!$                 write(*,*)"wrong! no two the same"
!!$              endif
!!$              !              write(*,*)xpoint,ypoint
!!$
!!$
!!$              if(mod(iiele,2).eq.0)then
!!$                 change=xpoint
!!$                 xpoint=ypoint 
!!$                 ypoint=change
!!$              endif
!!$
!!$              call gauss_point(1,shapee,weighttt,ypoint,xpoint)
!!$
!!$              do indexx=1,12
!!$                 shapefff(indexx,igg)=shapee(1,indexx)
!!$              enddo
!!$
!!$           endif
!!$
!!$        enddo
!!$        !               write(*,*)ccc
!!$
!!$     enddo
!!$
!!$
!!$
!!$  end do
!!$
!!$!!!!!!!!!!!!!!!          output          !!!!!!!!!!!!!!!!!!!
!!$  
!!$  write(31,*) 'Bond data'
!!$  write(31,*) ngaussfind
!!$  write(31,*) 0.1
!!$  do i=1,ngaussfind
!!$     write(31,*)indexelem(i)
!!$  enddo
!!$
!!$  do i=1,ngaussfind
!!$     write(31,'(10d30.17)')shapefff(1:12,i)
!!$  enddo
!!$  
!!$  do i=1,12
!!$     write(31,'(10d30.17)') shapeff(1:2,i)
!!$  enddo
!!$
!!$  write(31,'(10d30.17)') weightt(1:2)
!!$
!!$  close(31)
!!$
!!$

!***********************************************************************************************************************
!***********************************************************************************************************************
!***********************************************************************************************************************


  write(*,*) ' Now we write the mesh'

  call gmsh_out(x0T,meshT,'meshini.msh',ntab)  

  write(*,*) ' We just wrote the mesh'

  ! START DIMS FILE
  write(90,*) 'Dims data'
  write(90,*) '---------'
  write(90,*) 'mesh0%numele'
  write(90,*) meshT%numele
  write(90,*) 'mesh0%numnods'
  write(90,*) meshT%numnods
  write(90,*) 'mesh0%nedge'
  write(90,*) meshT%nedge
  write(90,*) 'mesh0%nelem_ghost'
  write(90,*) meshT%nelem_ghost
  write(90,*) 'mesh0%nnode_ghost'
  write(90,*) meshT%nnode_ghost
  write(90,*) 'ngauss'
  write(90,*) ngauss
  write(90,*) 'BCs%nnodBC'
  write(90,*) BCsT%nnodBC
  write(90,*) 'BCs%ndofBC'
  write(90,*) BCsT%ndofBC
  write(90,*) 'BCs%ndofOP'
  write(90,*) BCsT%ndofOP
  write(90,*) 'vdw1%nvdw'
  write(90,*) vdwT%nvdw
  if (vdwT%nvdw.eq.1) then
     write(90,*) 'vdw1%ngauss_vdw'
     write(90,*) vdwT%ngauss_vdw
     write(90,*) 'vdw1%ng_tot'
     write(90,*) vdwT%ng_tot
     write(90,*) 'nneigh'
     write(90,*) nneighT
     write(90,*) 'ninrange'
     write(90,*) ninrangeT
  endif
  ! END DIMS FILE


  ! START GENERAL FILE
  write(91,*) 'General data'
  write(91,*) '------------'
  write(91,*) 'ylength'
  write(91,'(d30.17)') ylengthT
  write(91,*) 'mat1%A0'
  write(91,'(d30.17)') mat1%A0
  write(91,*) 'mat1%nCode_Pot'
  write(91,*) mat1%nCode_Pot
  write(91,*) 'Material Parameters'
  if (mat1%nCode_Pot.eq.1) then
     ! in this case mat1%A0 should be 0.139
!!$  write(91,*) 1.207
!!$  write(91,*) 19.6279
!!$  write(91,*) 0.74
!!$  write(91,*) 0.754
     write(91,'(d30.17)') 0.603105
     write(91,'(d30.17)') 26.25
     write(91,'(d30.17)') 0.9
     write(91,'(d30.17)') 0.754
  else if (mat1%nCode_Pot.eq.2) then
     write(91,'(d30.17)')  0.139
     write(91,'(d30.17)')  0.9613062
     write(91,'(d30.17)')  21.
     write(91,'(d30.17)')  1.22
     write(91,'(d30.17)')  0.00020813
     write(91,'(d30.17)')  108900.
     write(91,'(d30.17)')  12.25
  else if (mat1%nCode_Pot.eq.22) then
     continue
  else if (mat1%nCode_Pot.eq.3) then
     ! in this case mat1%A0 should be 0.142
     write(91,'(d30.17)') 930.
     write(91,'(d30.17)') 0.74
  endif
  write(91,*) 'mat1%E'
  write(91,'(2d30.17)') mat1%E(1,:)
  write(91,'(2d30.17)') mat1%E(2,:)
  write(91,'(2d30.17)') mat1%E(3,:)
  write(91,*) 'mat1%s0'
  write(91,'(d30.17)') mat1%s0
  write(91,*) 'nW_hat'
  write(91,*) nW_hat
  write(91,*) 'crit'
  write(91,'(d30.17)') crit(1)
  write(91,'(d30.17)') crit(2)
  write(91,*) 'imperfect'
  write(91,*) imperfect
  write(91,*) 'fact_imp'
  write(91,'(d30.17)') fact_imp
  ! END GENERAL FILE


  ! START MESH FILE
  write(94,*) 'Mesh data'
  write(94,*) '---------'
  write(94,*) 'Connect'
  do ielem=1,meshT%numele
     write(94,*) 'New element'
     write(94,*) ielem, meshT%connect(ielem)%vertices(1:3)
     write(94,*) meshT%connect(ielem)%num_neigh_elem
     write(94,*) meshT%connect(ielem)%num_neigh_vert
     do jj=1,12
        write(94,*) meshT%connect(ielem)%neigh_elem(jj), meshT%connect(ielem)%neigh_vert(jj)
     enddo
     write(94,*) meshT%connect(ielem)%code_bc(1:3)
  enddo
  !write(94,*) 'ntable'
  !do inod=1,meshT%numnods
  !  write(94,'(13i6)') meshT%ntable(inod,:)
  !enddo
  write(94,*) 'nghost_tab'
  do inod=1,meshT%nedge
     write(94,'(3i6)') meshT%nghost_tab(inod,:)
  enddo
  write(94,*) 'nelem_ghost'
  do inod=1,meshT%nelem_ghost
     write(94,'(3i6)') meshT%elem_ghost(inod)
  enddo
  write(94,*) 'nnode_ghost'
  do inod=1,meshT%nnode_ghost
     write(94,'(3i6)') meshT%node_ghost(inod)
  enddo
  ! END MESH FILE


  ! START ZERO FILE
  write(95,*) 'Zero data'
  write(95,*) '------------'
  do ielem=1,meshT%numele
     write(95,'(d30.17)') J0T(ielem)
     write(95,'(2d30.17)') F0T(ielem)%val(1,1),F0T(ielem)%val(1,2)
     write(95,'(2d30.17)') F0T(ielem)%val(2,1),F0T(ielem)%val(2,2)
  enddo
  ! END ZERO FILE


  ! START CONFIG FILE
  write(96,*) 'Config Data'
  write(96,*) '-----------'
  write(96,*) 'Nodal positions'
  do i=1,meshT%numnods
     write(96,'(3d30.17)') x0T(3*i-2:3*i)
  enddo
  write(96,*) 'Inner displacements'
  do ielem=1,meshT%numele
     do igauss=1,ngauss
        write(96,'(2d30.17)') 0.d0, 0.d0
     enddo
  enddo
  ! END CONFIG FILE


  ! START BCs FILE
  write(92,*) 'BCs data'
  write(92,*) '--------'
  write(92,*) 'BCs%nloadstep'
  write(92,*) BCsT%nloadstep
  write(92,*) 'BCs%nCodeLoad'
  write(92,*) BCsT%nCodeLoad
  write(92,*) 'BCs%mdofBC'
  do i=1,BCsT%ndofBC
     write(92,*) BCsT%mdofBC(i)
  enddo
  write(92,*) 'BCs%mdofOP'
  do i=1,BCsT%ndofOP
     write(92,*) BCsT%mdofOP(i)
  enddo
  write(92,*) 'BCs%mnodBC'
  do i=1,BCsT%nnodBC
     write(92,*) BCsT%mnodBC(i,1:2)
  enddo
  write(92,*) 'BCs%rotation'
  write(92,'(3d30.17)') BCsT%rotation(1,1:3)
  write(92,'(3d30.17)') BCsT%rotation(2,1:3)
  write(92,'(3d30.17)') BCsT%rotation(3,1:3)
  write(92,*) 'BCs%xc'
  write(92,'(3d30.17)') BCsT%xc(1:3)
  write(92,*) 'BCs%value'
  write(92,'(d30.17)') BCsT%value
  ! END BCs FILE


  ! START VDW FILE
  if (vdwT%nvdw.eq.1) then
     open(unit=93,file='nano_vdw.dat',status='unknown')
     write(93,*) 'vdw1%meval'
     write(93,*) vdwT%meval
     write(93,*) 'vdw1%r_cut'
     write(93,'(d30.17)') vdwT%r_cut
     write(93,*) 'vdw1%r_bond'
     write(93,'(d30.17)') vdwT%r_bond
     write(93,*) 'vdw1%sig'
     write(93,'(d30.17)') vdwT%sig
     write(93,*) 'vdw1%a'
     write(93,'(d30.17)') vdwT%a
     write(93,*) 'vdw1%y0'
     write(93,'(d30.17)') vdwT%y0
     write(93,*) 'vdw1%Vcut(2)'
     write(93,'(2d30.17)') vdwT%Vcut(1:2)
     write(93,*) 'vdw1%shapef'
     do i=1,12
        write(93,'(10d30.17)') vdwT%shapef(1:vdwT%ngauss_vdw,i)
     enddo
     write(93,*) 'vdw1%weight'
     write(93,'(10d30.17)') vdwT%weight(1:vdwT%ngauss_vdw)
     write(93,*) 'vdw1%near'
     if (nskip.eq.1) then
        continue
     else
        do i=1,vdwT%ng_tot
           write(93,*) vdwT%near(i,0:nneighT)
        enddo
     endif
     close(93)
  endif
  ! END VDW FILE

  kk9=0
  write(56,*) ntubes
  do ijk=1,ntubes
     kk9=kk9+numel(ijk)*vdwT%ngauss_vdw
     write(56,*) kk9
  enddo


  !@@@@@@@@@@@ FINALIZE SMOOTHLY @@@@@@@@@@@@@@@@@@@@@@@
  write(*,*) '**** Deallocate Memory ****8'
  DEALLOCATE(BCsT%mdofBC,BCsT%mdofOP,BCsT%mnodBC)
  DEALLOCATE(meshT%connect, &
       !       meshT%ntable, &
  meshT%nghost_tab)
  if (vdwT%nvdw.eq.1) then
     DEALLOCATE(vdwT%shapef,vdwT%weight,vdwT%near)
  endif
  DEALLOCATE(meshT%elem_ghost,meshT%node_ghost)

  DEALLOCATE(mesh0,BCs,vdw1,ntab,temp,mtable)

  DEALLOCATE(nrow,ncol,ylength,nneigh,ninrange,numel,numno,numed)
  DEALLOCATE(x0T,J0T,F0T)
!!$!***************************************************************************************************************
!!$  DEALLOCATE(x0inner,indexelem,shapefff)
!!$!****************************************************************************************************************

  close(90)
  close(91)
  close(92)
  close(94)
  close(95)
  close(96)
  close(56)

end program PREPRO

