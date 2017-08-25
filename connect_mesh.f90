SUBROUTINE connect_mesh(meshh,mtable)
USE data_mesh
implicit REAL(8) (a-h,o-z)
implicit INTEGER*4 (i-n)
logical ask, mask, flag
type(tri) :: triang, tri_master, tri_slave
TYPE(mesh) :: meshh
dimension  ivert(3), mask(maxneigh_vert), iaux(3)
!!$integer, dimension(:,:), allocatable :: mtable
integer :: mtable(meshh%numnods,maxneigh_vert+1)
dimension iperm(2,3), neigh(3), neigh2(3)
dimension ilist_el(4), ilist_ve(4)
dimension :: iright(6), ileft(6), ihelp(3), ikk(3,2)
data iperm(1,1:3) /2, 3, 1/
data iperm(2,1:3) /3, 1, 2/


!!$ALLOCATE(mtable(meshh%numnods,maxneigh_vert+1),STAT=istat)
!!$if (istat/=0) STOP '**** Not enough memory ****'

! Find all nodes connected to a particular node
! Find all the elements connected to particular node
meshh%ntable=0
mtable=0
do inode=1,meshh%numnods
  do ielem=1,meshh%numele
     ivert=meshh%connect(ielem)%vertices
     ask = any(ivert.eq.inode)
	 if (ask) then
       mtable(inode,1)=mtable(inode,1)+1
	   mtable(inode,1+mtable(inode,1))=ielem
	 end if
     if (ask) then
       do i=1,3
          iv=ivert(i)
          mask=meshh%ntable(inode,2:meshh%ntable(inode,1)+1).eq.iv
          ask=any(mask)
          if ((.not.ask).and.(iv/=inode)) then
            meshh%ntable(inode,1)=meshh%ntable(inode,1)+1
		    meshh%ntable(iv,1)=meshh%ntable(iv,1)+1
		    meshh%ntable(inode,1+meshh%ntable(inode,1))=iv
		    meshh%ntable(iv,1+meshh%ntable(iv,1))=inode
          end if
        end do
      end if
  end do
end do


! Fill in gaps in connect
ighost=0
do ielem=1,meshh%numele
  triang=meshh%connect(ielem)                         ! data of ielem stored in triang
  call find_elem_adj(ielem,triang,mtable,neigh,meshh%numnods) 
  ivert=triang%vertices
  iaux=meshh%ntable(ivert(1:3),1)
  if (maxval(iaux).gt.6) then
    write(*,*) 'Error, program not prepared for valences > 6'
  end if

  if (neigh(1)/=0) then
!   Search on first edge to the right
    istep=0
    ilist_el=[3,2,1,12]
    ilist_ve=[3,1,2,5]
    ipivot=ivert(1)             ! we are pivoting around this vertex
    master=ielem                ! master element is for the moment the original one 
    tri_master=meshh%connect(master)
    call glob_loc(ipivot,tri_master,i_master)  ! orient master element wrt ipivot
    do
      call find_elem_adj(master,tri_master,mtable,neigh2,meshh%numnods)  ! find neighbors of master element
      iact=i_master                                ! identify active edge of master
      mslave=neigh2(iact)                               ! identify salve element of master
      flag=(mslave.eq.0).or.(mslave.eq.neigh(3))      ! if active edge of master is boundary or ...
      if (flag) exit                   ! ... reach edge 3, finish with edge 1
      istep=istep+1                                     ! we advance one step
      tri_slave=meshh%connect(mslave)
      call glob_loc(ipivot,tri_slave,i_slave)           ! orient slave element wrt ipivot
      triang%neigh_elem(ilist_el(istep))=mslave         ! set neighbor element data
      triang%neigh_vert(ilist_ve(istep))=tri_slave%vertices(iperm(1,i_slave)) ! set neighbor vertex data
      master=mslave
      tri_master=tri_slave
      i_master=i_slave
    end do

!   Search on the first edge to the left
    istep=0
    ilist_el=[4,5,6,0]
    ilist_ve=[6,10,11,0]
    ipivot=ivert(2)             ! we are pivoting around this vertex
    master=triang%neigh_elem(3)  ! master element is for the moment the original one 
    tri_master=meshh%connect(master)
    call glob_loc(ipivot,tri_master,i_master)                ! orient master element wrt ipivot
    do
      call find_elem_adj(master,tri_master,mtable,neigh2,meshh%numnods)  ! find neighbors of master element
      iact=iperm(2,i_master)                                ! identify active edge of master
      mslave=neigh2(iact)                               ! identify salve element of master
      flag=(mslave.eq.0).or.(mslave.eq.neigh(2))        ! if active edge of master is boundary or ...
      if (flag) exit                                    ! ... we reach edge 2 of ielem, finish with edge 1
      istep=istep+1                                     ! we advance one step
      tri_slave=meshh%connect(mslave)
      call glob_loc(ipivot,tri_slave,i_slave)           ! orient slave element wrt ipivot
      triang%neigh_elem(ilist_el(istep))=mslave         ! set neighbor element data
      triang%neigh_vert(ilist_ve(istep))=tri_slave%vertices(iperm(2,i_slave)) ! set neighbor vertex data

      master=mslave
      tri_master=tri_slave
      i_master=i_slave
    end do
  else
! insert ghost node
    ighost=ighost+1
	meshh%ntable(ivert(1),1)=meshh%ntable(ivert(1),1)+1
	meshh%ntable(ivert(2),1)=meshh%ntable(ivert(2),1)+1
    meshh%ntable(ivert(1),1+meshh%ntable(ivert(1),1))=meshh%numnods+ighost
    meshh%ntable(ivert(2),1+meshh%ntable(ivert(2),1))=meshh%numnods+ighost
    triang%neigh_vert(3)=meshh%numnods+ighost
  endif


  if (neigh(2)/=0) then
!   Search on second edge to the right
    istep=0
    ilist_el=[7,6,5,4]
    ilist_ve=[11,10,6,3]
    ipivot=ivert(2)             ! we are pivoting around this vertex
    master=ielem                ! master element is for the moment the original one 
    tri_master=meshh%connect(master)
    call glob_loc(ipivot,tri_master,i_master)  ! orient master element wrt ipivot
    do
      call find_elem_adj(master,tri_master,mtable,neigh2,meshh%numnods)  ! find neighbors of master element
      iact=i_master                                ! identify active edge of master
      mslave=neigh2(iact)                               ! identify salve element of master
      flag=(mslave.eq.0).or. &
  	     (triang%neigh_elem(ilist_el(istep+1))/=0)    ! if active edge of master is boundary or ...
      if (flag) exit                                    ! ... the data already filled, finish with edge 1
      istep=istep+1                                     ! we advance one step
      tri_slave=meshh%connect(mslave)
      call glob_loc(ipivot,tri_slave,i_slave)           ! orient slave element wrt ipivot
      triang%neigh_elem(ilist_el(istep))=mslave         ! set neighbor element data
      if ((triang%neigh_vert(ilist_ve(istep))/=0).and. &
	     (triang%neigh_vert(ilist_ve(istep))/=tri_slave%vertices(iperm(1,i_slave)))) &
		 write(*,*) ' Gross error in second right '
      triang%neigh_vert(ilist_ve(istep))=tri_slave%vertices(iperm(1,i_slave)) ! set neighbor vertex data
      master=mslave
      tri_master=tri_slave
      i_master=i_slave
    end do
 
!   Search on the second edge left
    istep=0
    ilist_el=[8,9,10,0]
    ilist_ve=[12,9,5,0]
    ipivot=ivert(3)             ! we are pivoting around this vertex
    master=triang%neigh_elem(7)         ! master element 
    tri_master=meshh%connect(master)
    call glob_loc(ipivot,tri_master,i_master)                ! orient master element wrt ipivot
    do
      call find_elem_adj(master,tri_master,mtable,neigh2,meshh%numnods)  ! find neighbors of master element
      iact=iperm(2,i_master)                                ! identify active edge of master
      mslave=neigh2(iact)                               ! identify salve element of master
      flag=(mslave.eq.0).or.(mslave.eq.neigh(3))        ! if active edge of master is boundary or ...
      if (flag) exit                                    ! ... we reach edge 2 of ielem, finish with edge 1
      istep=istep+1                                     ! we advance one step
      tri_slave=meshh%connect(mslave)
      call glob_loc(ipivot,tri_slave,i_slave)           ! orient slave element wrt ipivot
      triang%neigh_elem(ilist_el(istep))=mslave         ! set neighbor element data
      if ((triang%neigh_vert(ilist_ve(istep))/=0).and. &
	     (triang%neigh_vert(ilist_ve(istep))/=tri_slave%vertices(iperm(2,i_slave)))) &
		 write(*,*) ' Gross error in second left '
      triang%neigh_vert(ilist_ve(istep))=tri_slave%vertices(iperm(2,i_slave)) ! set neighbor vertex data

      master=mslave
      tri_master=tri_slave
      i_master=i_slave
    end do
  else
! insert ghost node
    ighost=ighost+1
	meshh%ntable(ivert(3),1)=meshh%ntable(ivert(3),1)+1
	meshh%ntable(ivert(2),1)=meshh%ntable(ivert(2),1)+1
    meshh%ntable(ivert(3),1+meshh%ntable(ivert(3),1))=meshh%numnods+ighost
    meshh%ntable(ivert(2),1+meshh%ntable(ivert(2),1))=meshh%numnods+ighost
    triang%neigh_vert(11)=meshh%numnods+ighost
  end if

  if (neigh(3)/=0) then
!   Search on third edge to the right
    istep=0
    ilist_el=[11,10,9,8]
    ilist_ve=[5,9,12,11]
    ipivot=ivert(3)             ! we are pivoting around this vertex
    master=ielem                ! master element is for the moment the original one 
    tri_master=meshh%connect(master)
    call glob_loc(ipivot,tri_master,i_master)  ! orient master element wrt ipivot
    do
      call find_elem_adj(master,tri_master,mtable,neigh2,meshh%numnods)  ! find neighbors of master element
      iact=i_master                                ! identify active edge of master
      mslave=neigh2(iact)                               ! identify salve element of master
      flag=(mslave.eq.0).or. &
  	     (triang%neigh_elem(ilist_el(istep+1))/=0)    ! if active edge of master is boundary or ...
      if (flag) exit                                    ! ... the data already filled, finish with edge 1
      istep=istep+1                                     ! we advance one step
      tri_slave=meshh%connect(mslave)
      call glob_loc(ipivot,tri_slave,i_slave)           ! orient slave element wrt ipivot
      triang%neigh_elem(ilist_el(istep))=mslave         ! set neighbor element data
      if ((triang%neigh_vert(ilist_ve(istep))/=0).and. &
	     (triang%neigh_vert(ilist_ve(istep))/=tri_slave%vertices(iperm(1,i_slave)))) &
		 write(*,*) ' Gross error in third right '
      triang%neigh_vert(ilist_ve(istep))=tri_slave%vertices(iperm(1,i_slave)) ! set neighbor vertex data
      master=mslave
      tri_master=tri_slave
      i_master=i_slave
    end do
 
!   Search on the third edge left
    istep=0
    ilist_el=[12,1,2,0]
    ilist_ve=[2,1,3,0]
    ipivot=ivert(1)             ! we are pivoting around this vertex
    master=triang%neigh_elem(11)         ! master element 
    tri_master=meshh%connect(master)
    call glob_loc(ipivot,tri_master,i_master)                ! orient master element wrt ipivot
    do
      call find_elem_adj(master,tri_master,mtable,neigh2,meshh%numnods)  ! find neighbors of master element
      iact=iperm(2,i_master)                                ! identify active edge of master
      mslave=neigh2(iact)                               ! identify salve element of master
      flag=(mslave.eq.0).or.(mslave.eq.neigh(1))        ! if active edge of master is boundary or ...
      if (flag) exit                                    ! ... we reach edge 1 of ielem, finish with edge 1
      istep=istep+1                                     ! we advance one step
      tri_slave=meshh%connect(mslave)
      call glob_loc(ipivot,tri_slave,i_slave)           ! orient slave element wrt ipivot
      triang%neigh_elem(ilist_el(istep))=mslave         ! set neighbor element data
      if ((triang%neigh_vert(ilist_ve(istep))/=0).and. &
	     (triang%neigh_vert(ilist_ve(istep))/=tri_slave%vertices(iperm(2,i_slave)))) &
		 write(*,*) ' Gross error in third left '
      triang%neigh_vert(ilist_ve(istep))=tri_slave%vertices(iperm(2,i_slave)) ! set neighbor vertex data

      master=mslave
      tri_master=tri_slave
      i_master=i_slave
    end do
  else
! insert ghost node
    ighost=ighost+1
	meshh%ntable(ivert(3),1)=meshh%ntable(ivert(3),1)+1
	meshh%ntable(ivert(1),1)=meshh%ntable(ivert(1),1)+1
    meshh%ntable(ivert(3),1+meshh%ntable(ivert(3),1))=meshh%numnods+ighost
    meshh%ntable(ivert(1),1+meshh%ntable(ivert(1),1))=meshh%numnods+ighost
    triang%neigh_vert(5)=meshh%numnods+ighost
  end if


  triang%neigh_vert(4)=triang%vertices(1)
  triang%neigh_vert(7)=triang%vertices(2)
  triang%neigh_vert(8)=triang%vertices(3)
  triang%num_neigh_elem=count(triang%neigh_elem/=0)
  triang%num_neigh_vert=count(triang%neigh_vert/=0)

  meshh%connect(ielem)=triang
end do


! Xtra data structure, meshh%nghost_tab, for easiness
jghost=0

do ielem=1,meshh%numele
  ned=count(meshh%connect(ielem)%code_bc==1)
  if (ned>0) then
    ivert=meshh%connect(ielem)%vertices
    do ied=1,3
	  if (meshh%connect(ielem)%code_bc(ied)==1) then
         jghost=jghost+1
         i1=ied
		 i2=iperm(1,ied)
		 i3=iperm(2,ied) 
         meshh%nghost_tab(jghost,:)=[ivert(i1),ivert(i2),ivert(i3)]
      end if	
	end do    
  end if
end do

if (jghost/=meshh%nedge) write(*,*) 'jghost different to nedge'

! Enhance presence of Ghost nodes in connect data structure
! Include completely Ghost nodes in the neigh_vert list, not only in positions 3,11,5
! This part is rather ah-hoc and should be revised for more general topologies
iright=[11,12, 3, 4, 7, 8]
ileft= [ 2, 3, 6, 7,10,11]
ihelp=[3,11,5]
ikk(1,:)=[2, 1]
ikk(2,:)=[6,10]
ikk(3,:)=[12,9]


do ielem=1,meshh%numele
  ivert=meshh%connect(ielem)%vertices
  if (meshh%connect(ielem)%num_neigh_vert/=12) then  ! this element has some empty socket
    do iv=1,3

       ! Look at the first possible free vertex conected to node #iv of ielem
       ivert1=ikk(iv,1)
       if (meshh%connect(ielem)%neigh_vert(ivert1)==0) then !need to fill it
   	     nelemr=meshh%connect(ielem)%neigh_elem(iright(2*iv-1))
	     neleml=meshh%connect(ielem)%neigh_elem(ileft(2*iv-1))
         ! try on the right
         if (nelemr>0) then
            call glob_loc(ivert(iv),meshh%connect(nelemr),ilocr)
            itest=meshh%connect(nelemr)%neigh_vert(ihelp(iperm(2,ilocr)))
            if (itest<meshh%numnods) STOP ' Error 999'
   		    meshh%connect(ielem)%num_neigh_vert=meshh%connect(ielem)%num_neigh_vert + 1
   		    meshh%connect(ielem)%neigh_vert(ivert1)=itest
	     else if (neleml>0) then
         ! try on the left
            call glob_loc(ivert(iv),meshh%connect(neleml),ilocl)
            itest=meshh%connect(neleml)%neigh_vert(ihelp(ilocl))
            if (itest<meshh%numnods) STOP ' Error 999'
   		    meshh%connect(ielem)%num_neigh_vert=meshh%connect(ielem)%num_neigh_vert + 1
   		    meshh%connect(ielem)%neigh_vert(ivert1)=itest
         else
		    STOP ' Error 000'
         end if
       end if

       ! Look at the second possible free vertex conected to node #iv of ielem
       ivert1=ikk(iv,2)
       if (meshh%connect(ielem)%neigh_vert(ivert1)==0) then !need to fill it
   	     nelemr=meshh%connect(ielem)%neigh_elem(iright(2*iv))
	     neleml=meshh%connect(ielem)%neigh_elem(ileft(2*iv))
         ! try on the right
         if (nelemr>0) then
            call glob_loc(ivert(iv),meshh%connect(nelemr),ilocr)
            itest=meshh%connect(nelemr)%neigh_vert(ihelp(iperm(2,ilocr)))
            if (itest<meshh%numnods) STOP ' Error 999'
   		    meshh%connect(ielem)%num_neigh_vert=meshh%connect(ielem)%num_neigh_vert + 1
   		    meshh%connect(ielem)%neigh_vert(ivert1)=itest
	     else if (neleml>0) then
         ! try on the left
            call glob_loc(ivert(iv),meshh%connect(neleml),ilocl)
            itest=meshh%connect(neleml)%neigh_vert(ihelp(ilocl))
            if (itest<meshh%numnods) STOP ' Error 999'
   		    meshh%connect(ielem)%num_neigh_vert=meshh%connect(ielem)%num_neigh_vert + 1
   		    meshh%connect(ielem)%neigh_vert(ivert1)=itest
         else
		    STOP ' Error 000'
         end if
       end if

    end do
  end if
  if (meshh%connect(ielem)%num_neigh_vert/=12) STOP ' Error 555'
end do



if (ighost/=meshh%nedge) write(*,*) 'ighost different to nedge'

!!$DEALLOCATE(mtable)

END SUBROUTINE connect_mesh

!*********************************************************************************
!*********************************************************************************
!*********************************************************************************

! Finds elements that share an edge with a particular element
SUBROUTINE find_elem_adj(ielem,triang,mtable,neigh,numnods)
USE data_tri
implicit REAL(8) (a-h,o-z)
implicit INTEGER*4 (i-n)
logical ask
TYPE(tri) :: triang
dimension mtable(numnods,15), iperm(2,3), neigh(3)
data iperm(1,1:3) /2, 3, 1/
data iperm(2,1:3) /3, 1, 2/

neigh=0

do ivertex=1,3
  iv1=triang%vertices(ivertex)
  iv2=triang%vertices(iperm(1,ivertex))
  if (triang%code_bc(ivertex)/=1) then
! not a boundary edge
    do k=1,mtable(iv1,1)
       jelem=mtable(iv1,k+1)
       ask=any(mtable(iv2,2:(mtable(iv2,1)+1)).eq.jelem)
       if (ask.and.(jelem/=ielem)) then
         neigh(ivertex)=jelem
    	 exit
       end if
    end do
    if (neigh(ivertex).eq.0) write(*,*) ' Error in neighbors of',ielem,ivertex
  end if
end do


END SUBROUTINE find_elem_adj

!*********************************************************************************
!*********************************************************************************
!*********************************************************************************

! given a global numbering, find the local numbering in a given element
! 
SUBROUTINE glob_loc(iglob,triang,iloc)
USE data_tri
implicit REAL(8) (a-h,o-z)
implicit INTEGER*4 (i-n)
TYPE(tri) :: triang

iloc=0

if (iglob.eq.(triang%vertices(1))) iloc=1 
if (iglob.eq.(triang%vertices(2))) iloc=2 
if (iglob.eq.(triang%vertices(3))) iloc=3 
if (iloc.eq.0) write(*,*) 'Error in glob_loc'

END SUBROUTINE glob_loc

!*********************************************************************************
!*********************************************************************************
!*********************************************************************************

SUBROUTINE ghost_nodes(meshh,xx)
USE data_mesh
implicit REAL(8) (a-h,o-z)
implicit INTEGER*4 (i-n)
TYPE(mesh) :: meshh
REAL(8) ::  xx(3*(meshh%numnods+meshh%nedge)), xaux(3)
!REAL(8) :: a(3), b(3), c(3), d(3)
!dimension  ivert(3), iperm(2,3)
!data iperm(1,1:3) /2, 3, 1/
!data iperm(2,1:3) /3, 1, 2/


do ijk=1,meshh%nedge
  i1=meshh%nghost_tab(ijk,1)
  i2=meshh%nghost_tab(ijk,2)
  i3=meshh%nghost_tab(ijk,3)
  xaux=xx(3*i1-2:3*i1)+xx(3*i2-2:3*i2)-xx(3*i3-2:3*i3)    ! parallelogram
  xx(3*(meshh%numnods+ijk)-2:3*(meshh%numnods+ijk))=xaux  ! position of ghost
end do


END SUBROUTINE ghost_nodes

!*********************************************************************************
!*********************************************************************************
!*********************************************************************************
