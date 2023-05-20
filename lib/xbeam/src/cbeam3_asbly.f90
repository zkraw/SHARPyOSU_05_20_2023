!-> Copyright by Imperial College London, 2009
!
!-> Module.- CBEAM3_ASBLY Rafa Palacios. 15Jul2009 - Last Update 20/09/2011
!
!-> Language: FORTRAN90, Free format.
!
!-> Description.-
!
!  Assembly beam equations.
!
!-> Subroutines.-
!
!    -cbeam3_asbly_static : Assembly matrices for nonlinear static problem.
!    -cbeam3_asbly_modal  : Assembly matrices for modal analysis.
!    -cbeam3_asbly_dynamic: Assembly matrices for nonlinear dynamic problems.
!
!-> Remarks.-
!  1) HH (20.09.2011) Changes made according to new version of nlabs r3.0,
!     which include derivatives of follower forces and new slave2master trans
!
!  2) HH (01.11.2013) Need to use full integration in assembly of mass matrix.
!
!  3) ADC (19.12.2016) New applied forces description implemented, it still
!                      supports the old format using function overloading
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module cbeam3_asbly
  use xbeam_shared
  use lib_solv, only: solv_set_vec_rows_zero
  use lib_mat
  use debug_utils
  implicit none
 real(8),parameter,private,dimension(3,3):: Unit= &    ! Unit matrix.
&         reshape((/1.d0,0.d0,0.d0,0.d0,1.d0,0.d0,0.d0,0.d0,1.d0/),(/3,3/))

 interface cbeam3_asbly_static
    module procedure :: cbeam3_asbly_static_old
 end interface cbeam3_asbly_static
interface cbeam3_asbly_dynamic
    module procedure :: cbeam3_asbly_dynamic_old_interface,&
                        cbeam3_asbly_dynamic_new_interface
end interface cbeam3_asbly_dynamic

!interface cbeam3_asbly_modal
!   module procedure :: cbeam3_asbly_modal_updated
!end interface cbeam3_asbly_modal
 contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-> Subroutine CBEAM3_ASBLY_STATIC_OLD
!
!-> Description:
!
!   Assembly matrices for a nonlinear static problem.
!
!-> Remarks.-
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine cbeam3_asbly_static_old (numdof, n_elem, n_node, Elem,Node,Coords,Psi0,PosDefor,PsiDefor,Force, &
&                                Kglobal,Fglobal,Qglobal,Options,Mglobal, MRR)
  use lib_rot
  use lib_rotvect
  use lib_fem
  use lib_sparse
  use lib_cbeam3
  use lib_xbeam
  use debug_utils


! I/O variables.
  integer, intent(IN)    :: numdof
  integer, intent(IN)    :: n_elem
  integer, intent(IN)    :: n_node
  type(xbelem), intent(in) :: Elem(n_elem)           ! Element information.
  type(xbnode), intent(in) :: Node(n_node)           ! List of independent nodes.
  real(8),      intent(in) :: Coords    (n_node,3)   ! Initial coordinates of the grid points.
  real(8),      intent(in) :: Psi0      (n_elem,3,3) ! Initial CRV of the nodes in the elements.
  real(8),      intent(in) :: PosDefor  (n_node,3)   ! Current coordinates of the grid points
  real(8),      intent(in) :: PsiDefor  (n_elem,3,3) ! Current CRV of the nodes in the elements.
  real(8),      intent(in) :: Force     (n_node,6)   ! Force vector.
  real(8), intent(inout):: Kglobal   (numdof,numdof)     ! Sparse stiffness matrix.
  real(8), optional, intent(OUT)       :: Mglobal   (numdof+6,numdof+6)     ! Sparse mass matrix.
  real(8), optional, intent(OUT)        :: MRR(6, 6)
  real(8),      intent(out):: Qglobal   (numdof)     ! Discrete force vector.
  real(8), intent(inout):: Fglobal   (numdof,numdof)     ! Influence coefficients matrix for applied forces.
  type(xbopts), intent(in) :: Options           ! Solver parameters.

! Local variables.
  logical:: Flags(MaxElNod)                ! Auxiliary flags.
  integer:: i,j,i1,j1                      ! Counters.
  integer:: nn, mm                         ! Counter for node number in global ordering
  integer:: rr, cc
    ! Counters for rows (rr) and column (cc) - for Kglobal only - in Kglobal, Qglobal and Fglobal
  integer:: iElem                          ! Counter on the finite elements.
  integer:: NumE                           ! Number of elements in the model.
  integer:: NumNE                          ! Number of nodes in an element.
  integer:: NsphBC                         ! number of hinges
  real(8):: Melem (6*MaxElNod,6*MaxElNod)  ! Element mass matrix.
  real(8):: Felem (6*MaxElNod,6*MaxElNod)  ! Element force influence coefficients.
  real(8):: Kelem (6*MaxElNod,6*MaxElNod)  ! Element tangent stiffness matrix.
  real(8):: Qelem (6*MaxElNod)             ! Total generalized forces on the element.
  real(8):: rElem0(MaxElNod,6)             ! Initial Coordinates/CRV of nodes in the element.
  real(8):: rElem (MaxElNod,6)             ! Current Coordinates/CRV of nodes in the element.
  real(8):: ForceElem (MaxElNod,6)         ! Current forces/moments of nodes in the element.
  real(8):: SB2B1 (6*MaxElNod,6*MaxElNod)  ! Transformation from master to global node orientations.
  integer:: NumGaussMass
  logical:: with_MSS
  logical:: with_MRR
  real(8)   :: MRRelem(6, 6)

  integer, allocatable:: row_sphBC(:)      ! row in global matrices/vector associated with weakly enforced hinge BCs
  Qglobal = 0.0d0
  Fglobal = 0.0d0
  Kglobal = 0.0d0

  with_MSS = present(Mglobal)
  if (with_MSS) then
      Mglobal = 0.0d0
  end if
  with_MRR = present(MRR)
  if (with_MRR) then
      MRR = 0.0d0
  end if

    !   call print_matrix('asbly_force',Force)
! Loop in all elements in the model.
  NumE=size(Elem)


  do iElem=1,NumE
    Melem = 0.0d0
    Kelem=0.d0; Felem=0.d0; Qelem=0.d0; SB2B1=0.d0
    rElem = 0.0d0
    rElem0 = 0.0d0
    MRRelem = 0.0d0

    !print*, iElem
  ! Determine if the element nodes are master (Flag=T) or slave.
    Flags=.false.
    do i=1,Elem(iElem)%NumNodes
        ! print*, '  ', i
        ! print*, Elem(iElem)%Conn(i)
      if (Node(Elem(iElem)%Conn(i))%Master(1).eq.iElem) Flags(i)=.true.
    end do
    !call flush()

  ! Extract components of the displacement and rotation vector at the element nodes
  ! and for the reference and current configurations.
  ! (Given the node nn (global numbering), associated to the node ii (local) of
  ! the element iElem, the function allocates in LocCoords(ii,:) the row Coords(nn,:).
  ! rElem0 and rElem are the output of this call.)
    call fem_glob2loc_extract (Elem(iElem)%Conn,Coords, rElem0(:,1:3),NumNE)
    call fem_glob2loc_extract (Elem(iElem)%Conn,PosDefor,rElem(:,1:3),NumNE)

! Use full integration for mass matrix.
    if (NumNE.eq.2) then
        NumGaussMass=NumNE
    elseif (NumNE.eq.3) then
        NumGaussMass=NumNE
    end if

  ! same is done here but there is no need to call fem_glob2loc_extract as the
  ! Psi arrays are already organised by elements.
    rElem0(:,4:6)= Psi0    (iElem,:,:)
    rElem (:,4:6)= PsiDefor(iElem,:,:)

  ! this adds or subtract 2 pi to the CRVs to ensure that they are all in the
  ! same part of the of the discontinuity [-pi,pi]
    call rotvect_boundscheck2(rElem(1,4:6),rElem(2,4:6))
    if (NumNE.eq.3) &
    &   call rotvect_boundscheck2(rElem(3,4:6),rElem(2,4:6))

  ! Extract current applied forces/moments at the element nodes.
    call fem_glob2loc_extract (Elem(iElem)%Conn,Force,ForceElem,NumNE)
    ! print*, iElem
    ! print*, Elem(iElem)%Conn
    ! print*, ForceElem
    ! pause

  !----------------------------------- ADDED SECTION
    if (with_MSS) then
    ! Contributions of the structural mass to the linearized inertia matrices.
        call cbeam3_rbmass (NumNE,rElem0,rElem,Elem(iElem)%RBMass,Melem)
        call cbeam3_mass (NumNE,rElem0,rElem,Elem(iElem)%Mass,Melem,NumGaussMass)
    end if

    if (with_MRR) then
        call xbeam_mrr  (NumNE,rElem0,rElem              ,Elem(iElem)%Mass,MRRelem,NumGaussMass)
        if (any(abs(Elem(iElem)%RBMass)>0.d0)) then
          call xbeam_rbmrr  (NumNE,rElem0,rElem,                                Elem(iElem)%RBMass,MRRelem)
        end if
        MRR = MRR + MRRelem
    end if

  !----------------------------------- END OF ADDED SECTION
  ! Compute the element tangent stiffness matrix and force vectors.
  ! inportant:
  ! 1. rElem and rEleme0 contain both positions and rotations!
  ! 2. cbeam_fstif adds to Qelem (inout) the contribution given by the internal
  !    forces. Note that Qelem at this point is zero.
    call cbeam3_kmat  (NumNE,rElem0,rElem,Elem(iElem)%Stiff,Kelem,Options%NumGauss)
    call cbeam3_kgeom (NumNE,rElem0,rElem,Elem(iElem)%Stiff,Kelem,Options%NumGauss)
    call cbeam3_fstif (NumNE,rElem0,rElem,Elem(iElem)%Stiff,Qelem,Options%NumGauss)
! Compute the influence coefficients multiplying the vector of external forces.
! - rotate if follower forces and filter out slave nodes.
! - The contribution due to follower forces is added to the stiffness matrix
! in cbeam3_dqext
    call cbeam3_fext  (NumNE,rElem,Flags(1:NumNE),Felem,Options%FollowerForce,Options%FollowerForceRig,Unit)
    call cbeam3_dqext (NumNE,rElem,ForceElem,Flags(1:NumNE),Kelem,Options%FollowerForce)
! Project equations to the orientation of the "master" degrees of freedom.
!    call cbeam3_projs2m (NumNE,Elem(iElem)%Master(:,:),Psi0(iElem,:,:),Psi0,SB2B1)
    call cbeam3_slave2master (NumNE,Elem(iElem)%Master(:,:),rElem0(:,4:6),Psi0,rElem(:,4:6),PsiDefor,SB2B1)
    Kelem=matmul(transpose(SB2B1),matmul(Kelem,SB2B1))
    Felem=matmul(transpose(SB2B1),Felem) !sm: ???
    Qelem=matmul(transpose(SB2B1),Qelem)
    if (with_MSS) then
        ! if (iElem == 1) then
        !     call print_matrix('Melem_before', Melem)
        ! end if
        Melem=matmul(transpose(SB2B1),matmul(Melem,SB2B1))
        ! if (iElem == 1) then
        !     call print_matrix('Melem_after', Melem)
        ! end if
        ! stop
    end if

  ! Add to global matrix Mglobal, not remove columns and rows
    !if (options%gravity_on)
        !do i=1, NumNE
            !nn = Elem(iElem)%Conn(i)
            !i1=Node(nn)%Vdof
            !rr = 6*i1
            !do j=1, NumNE
                !mm=Elem(iElem)%Conn(j)
                !j1=Node(mm)%Vdof
                !cc = 6*j1
                !call sparse_addmat (rr, cc,Melem(6*(i-1)+1:6*i,6*(j-1)+1:6*j), ms,Mglobal)
            !end do
        !end do
    !end if

 ! Add to global matrix. Remove columns and rows at clamped points.
    do i=1,NumNE
      nn = Elem(iElem)%Conn(i)
      i1=Node(nn)%Vdof
      if (i1.ne.0) then ! not clamped: free node (internal or end)or hinged
        rr = 6*(i1-1)
          Qglobal(rr+1:rr+6)= Qglobal(rr+1:rr+6)+Qelem(6*(i-1)+1:6*i)
          do j=1,NumNE
            mm=Elem(iElem)%Conn(j)
            j1=Node(mm)%Vdof
            ! determine position of 1st dof related to node Elem(iElem)%Conn(j)
            ! in global matrices/vectors minus 1
            cc = 6*(j1-1)
            if (j1.ne.0) then ! not clamped
              call mat_addmat(rr, cc,Kelem(6*(i-1)+1:6*i,6*(j-1)+1:6*j),Kglobal)
              call mat_addmat(rr, cc,Felem(6*(i-1)+1:6*i,6*(j-1)+1:6*j),Fglobal)
            end if
          end do
      end if
    end do
 ! Add to global matrix. Dont remove columns and rows at clamped points.
 if (with_MSS) then
    do i=1,NumNE
      nn = Elem(iElem)%Conn(i)
    !   i1=Node(nn)%Vdof + 1
      i1=Node(nn)%Vdof
        ! rr = 6*(i1-1)
        rr = 6*i1
          do j=1,NumNE
            mm=Elem(iElem)%Conn(j)
            ! j1=Node(mm)%Vdof + 1
            j1=Node(mm)%Vdof
            ! cc = 6*(j1-1)
            cc = 6*j1
            call mat_addmat (rr, cc,Melem(6*(i-1)+1:6*i,6*(j-1)+1:6*j),Mglobal)
          end do
    end do
end if
end do

  ! ------------------------------------------------------ Hinged BCs correction
  ! a. determine number of rows to delete
  ! b. determine which rows to eleminate
  ! c. delete themn from sparse matrix and resize it
  ! d. add unitary value to Kglobal and resize)
  ! NsphBC = sum(Node%Sflag)
  ! if (NsphBC>0) then
  !     allocate(row_sphBC (3*NsphBC) )
  !
  !     cc=0 ! here cc is a counter for the number of spherical joints
  !     do nn=1,size(Node)
  !       if (Node(nn)%Sflag == 1) then
  !         i1=Node(nn)%Vdof
  !         rr = 6*(i1-1)
  !         row_sphBC( 3*cc+1:3*(cc+1) ) = (/ rr+1, rr+2, rr+3 /)
  !         cc=cc+1
  !       end if
  !     end do
  !
  !     call sparse_set_rows_zero(row_sphBC,ks,Kglobal)
  !     call sparse_set_rows_zero(row_sphBC,fs,Fglobal)
  !     call sparse_set_rows_unit(row_sphBC,ks,Kglobal)
  !
  !     deallocate(row_sphBC)
  !  end if

    ! call print_matrix('MSS', Mglobal)



  return
 end subroutine cbeam3_asbly_static_old

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-> Subroutine CBEAM3_ASBLY_MODAL_OLD
!
!-> Description:
!
!   Assembly matrices for a modal analysis.
!
!-> Remarks.-
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine cbeam3_asbly_modal_old (Elem,Node,Coords,Psi0,PosDefor,PsiDefor,  &
&                               Vrel,ms,Mglobal,cs,Cglobal,ks,Kglobal,Options)
  use lib_rotvect
  use lib_fem
  use lib_sparse
  use lib_cbeam3

! I/O variables.
  type(xbelem),intent(in) :: Elem(:)               ! Element information.
  type(xbnode),intent(in) :: Node(:)               ! List of independent nodes.
  real(8),     intent(in) :: Coords    (:,:)       ! Initial coordinates of the grid points.
  real(8),     intent(in) :: Psi0      (:,:,:)     ! Initial CRV of the nodes in the elements.
  real(8),     intent(in) :: PosDefor  (:,:)       ! Current coordinates of the grid points
  real(8),     intent(in) :: PsiDefor  (:,:,:)     ! Current CRV of the nodes in the elements.
  real(8),     intent(in) :: Vrel(6)               ! Velocity of reference frame .

  integer,     intent(out):: ms                    ! Size of the sparse mass matrix.
  type(sparse),intent(out):: Mglobal   (:)         ! Sparse mass matrix.
  integer,     intent(out):: cs                    ! Size of the sparse damping matrix.
  type(sparse),intent(out):: Cglobal   (:)         ! Sparse damping matrix.
  integer,     intent(out):: ks                    ! Size of the sparse stiffness matrix.
  type(sparse),intent(out):: Kglobal   (:)         ! Sparse stiffness matrix.
  type(xbopts),intent(in) :: Options               ! Solver parameters.

! Local variables.
  integer:: i,j,i1,j1                      ! Counters.
  integer:: iElem                          ! Counter on the finite elements.
  integer:: NumE                           ! Number of elements in the model.
  integer:: NumNE                          ! Number of nodes in an element.
  integer:: NumGaussMass                   ! Number of Gaussian points in the inertia terms.
  real(8):: Celem (6*MaxElNod,6*MaxElNod)  ! Element damping matrix.
  real(8):: Kelem (6*MaxElNod,6*MaxElNod)  ! Element tangent stiffness matrix.
  real(8):: Melem (6*MaxElNod,6*MaxElNod)  ! Element mass matrix.
  real(8):: rElem0(MaxElNod,6)             ! Initial Coordinates/CRV of nodes in the element.
  real(8):: rElem (MaxElNod,6)             ! Current Coordinates/CRV of nodes in the element.
  real(8):: rElemDot (MaxElNod,6)          ! Current Coordinates/CRV of nodes in the element.
  real(8):: rElemDDot (MaxElNod,6)         ! Current Coordinates/CRV of nodes in the element.
  real(8):: VrelDot  (6)                   ! Time derivative of Vrel.
  real(8):: SB2B1 (6*MaxElNod,6*MaxElNod)  ! Transformation from master to global node orientations.

! Initialize.
  NumE=size(Elem)
  rElemDot = 0.d0
  rElemDDot= 0.d0
  VrelDot  = 0.d0

  do iElem=1,NumE
    Melem=0.d0; Celem=0.d0; Kelem=0.d0; SB2B1=0.d0

! Extract coords of elem nodes and determine if they are master (Flag=T) or slave.
    call fem_glob2loc_extract (Elem(iElem)%Conn,Coords,  rElem0(:,1:3),NumNE)
    call fem_glob2loc_extract (Elem(iElem)%Conn,PosDefor,rElem (:,1:3),NumNE)
    rElem0 (:,4:6)= Psi0     (iElem,:,:)
    rElem  (:,4:6)= PsiDefor (iElem,:,:)

! Use full integration for mass matrix.
    if (NumNE.eq.2) then
        NumGaussMass=NumNE
    elseif (NumNE.eq.3) then
        NumGaussMass=NumNE
    end if

! Compute linearized inertia matrices.
    call cbeam3_mass (NumNE,rElem0,rElem,                                Elem(iElem)%Mass,Melem,NumGaussMass)
    call cbeam3_cgyr (NumNE,rElem0,rElem,rElemDot,Vrel,                  Elem(iElem)%Mass,Celem,Options%NumGauss)
    call cbeam3_kgyr (NumNE,rElem0,rElem,rElemDot,rElemDDot,Vrel,VrelDot,Elem(iElem)%Mass,Kelem,Options%NumGauss)

    ! Add contributions of non-structural (lumped) mass.
    if (any(abs(Elem(iElem)%RBMass)>0.d0)) then
      call cbeam3_rbmass (NumNE,rElem0,rElem,                                Elem(iElem)%RBMass,Melem)
      call cbeam3_rbcgyr (NumNE,rElem0,rElem,rElemDot,          Vrel,        Elem(iElem)%RBMass,Celem)
      call cbeam3_rbkgyr (NumNE,rElem0,rElem,rElemDot,rElemDDot,Vrel,VrelDot,Elem(iElem)%RBMass,Kelem)
    end if

! Compute the element tangent stiffness matrix and force vectors.
    call cbeam3_kmat  (NumNE,rElem0,rElem,Elem(iElem)%Stiff,Kelem,Options%NumGauss)
    call cbeam3_kgeom (NumNE,rElem0,rElem,Elem(iElem)%Stiff,Kelem,Options%NumGauss)

! Project slave degrees of freedom to the orientation of the "master" ones.
    call cbeam3_projs2m (NumNE,Elem(iElem)%Master(:,:),Psi0(iElem,:,:),Psi0,SB2B1)
    Melem=matmul(transpose(SB2B1),matmul(Melem,SB2B1))
    Celem=matmul(transpose(SB2B1),matmul(Celem,SB2B1))
    Kelem=matmul(transpose(SB2B1),matmul(Kelem,SB2B1))

! Add to global matrix. Remove columns and rows at clamped points.
    do i=1,NumNE
      i1=Node(Elem(iElem)%Conn(i))%Vdof
      if (i1.ne.0) then
        do j=1,NumNE
          j1=Node(Elem(iElem)%Conn(j))%Vdof
          if (j1.ne.0) then
            call sparse_addmat (6*(i1-1),6*(j1-1),Melem(6*(i-1)+1:6*i,6*(j-1)+1:6*j),ms,Mglobal)
            call sparse_addmat (6*(i1-1),6*(j1-1),Celem(6*(i-1)+1:6*i,6*(j-1)+1:6*j),cs,Cglobal)
            call sparse_addmat (6*(i1-1),6*(j1-1),Kelem(6*(i-1)+1:6*i,6*(j-1)+1:6*j),ks,Kglobal)
          end if
        end do
      end if
    end do
  end do
  return
end subroutine cbeam3_asbly_modal_old

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-> Subroutine CBEAM3_ASBLY_MODAL_UPDATED
!
!-> Description:
!
!   Assembly matrices for a modal analysis.
!
!-> Remarks.-
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine cbeam3_asbly_modal_updated (num_dof,n_elem,n_node,Elem,Node,Coords,Psi0,PosDefor,PsiDefor,  &
&                               Vrel,FullMglobal,FullCglobal,FullKglobal,Options)
  use lib_rotvect
  use lib_fem
  use lib_mat
  use lib_cbeam3
  use lib_lu

! I/O variables.
  integer,      intent(in) :: num_dof               ! Degrees of freedom
  integer,      intent(in) :: n_elem                 ! Number of elements
  integer,      intent(in) :: n_node                  ! Number of nodes
  type(xbelem),intent(in) :: Elem(n_elem)               ! Element information.
  type(xbnode),intent(in) :: Node(n_node)               ! List of independent nodes.
  real(8),     intent(in) :: Coords   (n_node,3)      ! Initial coordinates of the grid points.
  real(8),     intent(in) :: Psi0      (n_elem,3,3)   ! Initial CRV of the nodes in the elements.
  real(8),     intent(in) :: PosDefor  (n_node,3)        ! Current coordinates of the grid points
  real(8),     intent(in) :: PsiDefor  (n_elem,3,3)    ! Current CRV of the nodes in the elements.
  real(8),     intent(in) :: Vrel(6)               ! Velocity of reference frame .

  real(8),     intent(out)    :: FullMglobal   (num_dof, num_dof)
  real(8),     intent(out)    :: FullCglobal   (num_dof, num_dof)
  real(8),     intent(out)    :: FullKglobal   (num_dof, num_dof)
real(8) :: temp   (num_dof, num_dof)
real(8) :: det

  type(xbopts),intent(in) :: Options               ! Solver parameters.

! Local variables.
  logical:: Flags(MaxElNod)                ! Auxiliary flags.
  integer:: i,j,i1,j1                      ! Counters.
  integer:: iElem                          ! Counter on the finite elements.
  integer:: NumE                           ! Number of elements in the model.
  integer:: NumNE                          ! Number of nodes in an element.
  integer:: NumGaussMass                   ! Number of Gaussian points in the inertia terms.
  real(8):: Celem (6*MaxElNod,6*MaxElNod)  ! Element damping matrix.
  real(8):: Kelem (6*MaxElNod,6*MaxElNod)  ! Element tangent stiffness matrix.
  real(8):: Melem (6*MaxElNod,6*MaxElNod)  ! Element mass matrix.
  real(8):: rElem0(MaxElNod,6)             ! Initial Coordinates/CRV of nodes in the element.
  real(8):: rElem (MaxElNod,6)             ! Current Coordinates/CRV of nodes in the element.
  real(8):: rElemDot (MaxElNod,6)          ! Current Coordinates/CRV of nodes in the element.
  real(8):: rElemDDot (MaxElNod,6)         ! Current Coordinates/CRV of nodes in the element.
  real(8):: VrelDot  (6)                   ! Time derivative of Vrel.
  real(8):: SB2B1 (6*MaxElNod,6*MaxElNod)  ! Transformation from master to global node orientations.

  real(8):: Kgeom (6*MaxElNod,6*MaxElNod)
  real(8):: Kmat (6*MaxElNod,6*MaxElNod)

! Initialize.
  NumE=size(Elem)
  FullMglobal = 0.d0
  FullCglobal = 0.d0
  FullKglobal = 0.d0
  rElemDot = 0.d0
  rElemDDot= 0.d0
  VrelDot  = 0.d0

  do iElem=1,NumE

    Melem=0.d0; Celem=0.d0; Kelem=0.d0; SB2B1=0.d0
    Kgeom=0.d0; Kmat=0.d0

    call fem_glob2loc_extract (Elem(iElem)%Conn,Coords,  rElem0(:,1:3),NumNE)
    call fem_glob2loc_extract (Elem(iElem)%Conn,PosDefor,rElem (:,1:3),NumNE)
    rElem0 (:,4:6)= Psi0     (iElem,:,:)
    rElem  (:,4:6)= PsiDefor (iElem,:,:)

! Use full integration for mass matrix.
    if (NumNE.eq.2) then
        NumGaussMass=NumNE
    elseif (NumNE.eq.3) then
        NumGaussMass=NumNE
    end if

! Compute linearized inertia matrices.
    call cbeam3_mass (NumNE,rElem0,rElem,                                Elem(iElem)%Mass,Melem,NumGaussMass)
    call cbeam3_cgyr (NumNE,rElem0,rElem,rElemDot,Vrel,                  Elem(iElem)%Mass,Celem,Options%NumGauss)
    call cbeam3_kgyr (NumNE,rElem0,rElem,rElemDot,rElemDDot,Vrel,VrelDot,Elem(iElem)%Mass,Kelem,Options%NumGauss)

    ! Add contributions of non-structural (lumped) mass.
    if (any(abs(Elem(iElem)%RBMass)>0.d0)) then
      write(*,*) "computing lumped masses"
      call cbeam3_rbmass (NumNE,rElem0,rElem,                                Elem(iElem)%RBMass,Melem)
      call cbeam3_rbcgyr (NumNE,rElem0,rElem,rElemDot,          Vrel,        Elem(iElem)%RBMass,Celem)
      call cbeam3_rbkgyr (NumNE,rElem0,rElem,rElemDot,rElemDDot,Vrel,VrelDot,Elem(iElem)%RBMass,Kelem)
    end if

! Compute the element tangent stiffness matrix and force vectors.
    call cbeam3_kmat  (NumNE,rElem0,rElem,Elem(iElem)%Stiff,Kelem,Options%NumGauss)
    call cbeam3_kgeom (NumNE,rElem0,rElem,Elem(iElem)%Stiff,Kelem,Options%NumGauss)

! Project slave degrees of freedom to the orientation of the "master" ones.
    !call cbeam3_projs2m (NumNE,Elem(iElem)%Master(:,:),Psi0(iElem,:,:),Psi0,SB2B1)
    call cbeam3_slave2master (NumNE,Elem(iElem)%Master(:,:),rElem0(:,4:6),Psi0,rElem(:,4:6),PsiDefor,SB2B1)
    Melem=matmul(transpose(SB2B1),matmul(Melem,SB2B1))
    Celem=matmul(transpose(SB2B1),matmul(Celem,SB2B1))
    Kelem=matmul(transpose(SB2B1),matmul(Kelem,SB2B1))

! Add to global matrix. Remove columns and rows at clamped points.
    do i=1,NumNE
      i1=Node(Elem(iElem)%Conn(i))%Vdof
      if (i1.ne.0) then
        do j=1,NumNE
          j1=Node(Elem(iElem)%Conn(j))%Vdof
          if (j1.ne.0) then
            call mat_addmat (6*(i1-1),6*(j1-1),Melem(6*(i-1)+1:6*i,6*(j-1)+1:6*j),FullMglobal)
            call mat_addmat (6*(i1-1),6*(j1-1),Celem(6*(i-1)+1:6*i,6*(j-1)+1:6*j),FullCglobal)
            call mat_addmat (6*(i1-1),6*(j1-1),Kelem(6*(i-1)+1:6*i,6*(j-1)+1:6*j),FullKglobal)
          end if
        end do
      end if
    end do
  end do

  return
end subroutine cbeam3_asbly_modal_updated


 subroutine cbeam3_asbly_dynamic_old_interface&
      (Elem,Node,Coords,Psi0,PosDefor,PsiDefor,         &
&                                 PosDeforDot,PsiDeforDot,PosDeforDDot,            &
&                                 PsiDeforDDot,Force,Vrel,VrelDot,ms,Mglobal,Mvel, &
&                                 cs,Cglobal,Cvel,ks,Kglobal,fs,Fglobal,           &
&                                 Qglobal,Options,Cao)
  use lib_rotvect
  use lib_fem
  use lib_sparse
  use lib_cbeam3

! I/O variables.
  type(xbelem),intent(in) :: Elem(:)               ! Element information.
  type(xbnode),intent(in) :: Node(:)               ! List of independent nodes.
  real(8),     intent(in) :: Coords    (:,:)       ! Initial coordinates of the grid points.
  real(8),     intent(in) :: Psi0      (:,:,:)     ! Initial CRV of the nodes in the elements.
  real(8),     intent(in) :: PosDefor  (:,:)       ! Current coordinates of the grid points
  real(8),     intent(in) :: PsiDefor  (:,:,:)     ! Current CRV of the nodes in the elements.
  real(8),     intent(in) :: PosDeforDot  (:,:)    ! Current coordinates of the grid points
  real(8),     intent(in) :: PsiDeforDot  (:,:,:)  ! Current CRV of the nodes in the elements.
  real(8),     intent(in) :: PosDeforDDot (:,:)    ! Current coordinates of the grid points
  real(8),     intent(in) :: PsiDeforDDot (:,:,:)  ! Current CRV of the nodes in the elements.
  real(8),     intent(in) :: Force     (:,:)       ! Force vector.
  real(8),     intent(in) :: Vrel(6), VrelDot(6)   ! Velocity of reference frame and derivative.

  integer,     intent(out):: ms                ! Size of the sparse mass matrix.
  type(sparse),intent(inout):: Mglobal   (:)     ! Sparse mass matrix.
  real(8),     intent(out):: Mvel      (:,:)   ! Reference system mass matrix.
  integer,     intent(out):: cs                ! Size of the sparse damping matrix.
  type(sparse),intent(inout):: Cglobal   (:)     ! Sparse damping matrix.
  real(8),     intent(out):: Cvel      (:,:)   ! Reference system damping matrix.
  integer,     intent(out):: ks                ! Size of the sparse stiffness matrix.
  type(sparse),intent(inout):: Kglobal   (:)     ! Sparse stiffness matrix.
  integer,     intent(out):: fs                ! Size of the sparse stiffness matrix.
  type(sparse),intent(inout):: Fglobal   (:)     ! Influence coefficients matrix for applied forces.
  real(8),     intent(out):: Qglobal   (:)     ! Stiffness and gyroscopic force vector.
  type(xbopts),intent(in) :: Options           ! Solver parameters.
  real(8),     intent(in) :: Cao       (:,:)   ! Rotation operator

  integer                   :: numdof
  integer                   :: n_node
  integer                   :: n_elem
  integer                   :: i

  n_node = size(Node)
  n_elem = size(Elem)

  ms = 0
  cs = 0
  ks = 0
  fs = 0


  numdof = 0
  do i=1, n_node
      if (Node(i)%vdof > 0)  then
          numdof = numdof + 6
      end if
  end do


  call cbeam3_asbly_dynamic_new_interface(numdof, n_node, n_elem, Elem,Node,Coords,Psi0,PosDefor,PsiDefor,         &
&                                 PosDeforDot,PsiDeforDot,PosDeforDDot,            &
&                                 PsiDeforDDot,Force,Vrel,VrelDot,Mglobal(1)%a,Mvel, &
&                                 Cglobal(1)%a,Cvel,Kglobal(1)%a,Fglobal(1)%a,           &
&                                 Qglobal,Options,Cao)

  end subroutine cbeam3_asbly_dynamic_old_interface

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-> Subroutine CBEAM3_ASBLY_DYNAMIC
!
!-> Description:
!
!   Assembly matrices for a dynamic problem.
!
!-> Remarks.-
!
!  1) Acceleration is used to compute the linearized stiffness matrix of the
!     inertia terms, although this can often be neglected.
!  2) RPN (20.04.2011) added lumped masses to the model.
!  3) RPN (20.04.2011) The number of Gauss points has been increased in all the
!     inertia matrices so that they are consistent among them.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine cbeam3_asbly_dynamic_new_interface (numdof, n_node, n_elem, Elem,Node,Coords,Psi0,PosDefor,PsiDefor,         &
&                                 PosDeforDot,PsiDeforDot,PosDeforDDot,            &
&                                 PsiDeforDDot,Force,Vrel,VrelDot,Mglobal,Mvel, &
&                                 Cglobal,Cvel,Kglobal,Fglobal,           &
&                                 Qglobal,Options,Cao)
  use lib_rotvect
  use lib_fem
  use lib_mat
  use lib_cbeam3

! I/O variables.
  integer,      intent(IN)      :: numdof
  integer,      intent(IN)      :: n_node
  integer,      intent(IN)      :: n_elem
  type(xbelem), intent(in)      :: Elem(n_elem)               ! Element information.
  type(xbnode), intent(in)      :: Node(n_node)               ! List of independent nodes.
  real(8),      intent(in)      :: Coords(n_node, 3)       ! Initial coordinates of the grid points.
  real(8),      intent(in)      :: Psi0(n_elem, 3, 3)     ! Initial CRV of the nodes in the elements.
  real(8),      intent(in)      :: PosDefor(n_node, 3)       ! Current coordinates of the grid points
  real(8),      intent(in)      :: PsiDefor(n_elem, 3, 3)     ! Current CRV of the nodes in the elements.
  real(8),      intent(in)      :: PosDeforDot(n_node, 3)    ! Current coordinates of the grid points
  real(8),      intent(in)      :: PsiDeforDot(n_elem, 3, 3)  ! Current CRV of the nodes in the elements.
  real(8),      intent(in)      :: PosDeforDDot(n_node, 3)   ! Current coordinates of the grid points
  real(8),      intent(in)      :: PsiDeforDDot(n_elem, 3, 3) ! Current CRV of the nodes in the elements.
  real(8),      intent(in)      :: Force(n_node, 6)       ! Force vector.
  real(8),      intent(in)      :: Vrel(6)
  real(8),      intent(IN)      :: VrelDot(6)   ! Velocity of reference frame and derivative.

  real(8),     intent(inout)    :: Mglobal   (numdof, numdof)     ! Sparse mass matrix.
  real(8),     intent(out)      :: Mvel      (numdof, 6)   ! Reference system mass matrix.
  real(8),     intent(inout)    :: Cglobal   (numdof, numdof)     ! Sparse damping matrix.
  real(8),     intent(out)      :: Cvel      (numdof, 6)   ! Reference system damping matrix.
  real(8),     intent(inout)    :: Kglobal   (numdof, numdof)     ! Sparse stiffness matrix.
  real(8),     intent(inout)    :: Fglobal   (numdof, numdof)     ! Influence coefficients matrix for applied forces.
  real(8),     intent(out)      :: Qglobal   (numdof)     ! Stiffness and gyroscopic force vector.
  type(xbopts),intent(in)       :: Options           ! Solver parameters.
  real(8),     intent(in)       :: Cao       (3, 3)   ! Rotation operator

! Local variables.
  logical:: Flags(MaxElNod)                ! Auxiliary flags.
  integer:: i,j,i1,j1                      ! Counters.
  integer:: iElem                          ! Counter on the finite elements.
  integer:: NumNE                          ! Number of nodes in an element.
  integer:: NumGaussMass                   ! Number of Gaussian points in the inertia terms.
  integer:: tempi1, tempj1, tempi, tempj
  real(8):: Celem (6*MaxElNod,6*MaxElNod)  ! Element damping matrix.
  real(8):: Felem (6*MaxElNod,6*MaxElNod)  ! Element force influence coefficients.
  real(8):: Kelem (6*MaxElNod,6*MaxElNod)  ! Element tangent stiffness matrix.
  real(8):: Melem (6*MaxElNod,6*MaxElNod)  ! Element mass matrix.
  real(8):: Qelem (6*MaxElNod)             ! Total generalized forces on the element.
  real(8):: Mvelelem(6*MaxElNod,6)         ! Element reference-system mass matrix.
  real(8):: Cvelelem(6*MaxElNod,6)         ! Element reference-system damping matrix.
  real(8):: rElem0(MaxElNod,6)             ! Initial Coordinates/CRV of nodes in the element.
  real(8):: rElem (MaxElNod,6)             ! Current Coordinates/CRV of nodes in the element.
  real(8):: rElemDot (MaxElNod,6)          ! Current Coordinates/CRV of nodes in the element.
  real(8):: rElemDDot (MaxElNod,6)         ! Current Coordinates/CRV of nodes in the element.
  real(8):: ForceElem (MaxElNod,6)         ! Current forces/moments of nodes in the element.
  real(8):: SB2B1 (6*MaxElNod,6*MaxElNod)  ! Transformation from master to global node orientations.

! Loop in all elements in the model.
  Mvel = 0.0d0
  Cvel = 0.0d0
  Qglobal = 0.0d0

  do iElem=1,n_elem
    Melem=0.d0; Celem=0.d0; Kelem=0.d0; Felem=0.d0; Qelem=0.d0
    Mvelelem=0.d0; Cvelelem=0.d0; SB2B1=0.d0
    rElem0 = 0.0d0
    rElem = 0.0d0

! Extract coords of elem nodes and determine if they are master (Flag=T) or slave.
    call fem_glob2loc_extract (Elem(iElem)%Conn,Coords,rElem0(:,1:3),NumNE)

    Flags=.false.
    do i=1,Elem(iElem)%NumNodes
      if (Node(Elem(iElem)%Conn(i))%Master(1).eq.iElem) Flags(i)=.true.
    end do

    call fem_glob2loc_extract (Elem(iElem)%Conn,PosDefor,    rElem    (:,1:3),NumNE)
    call fem_glob2loc_extract (Elem(iElem)%Conn,PosDeforDot, rElemDot (:,1:3),NumNE)
    call fem_glob2loc_extract (Elem(iElem)%Conn,PosDeforDDot,rElemDDot(:,1:3),NumNE)

    rElem0   (:,4:6)= Psi0(iElem,:,:)
    rElem    (:,4:6)= PsiDefor(iElem,:,:)
    rElemDot (:,4:6)= PsiDeforDot(iElem,:,:)
    rElemDDot(:,4:6)= PsiDeforDDot(iElem,:,:)

    call rotvect_boundscheck2(rElem(1,4:6),rElem(2,4:6))
    if (NumNE.eq.3) call rotvect_boundscheck2(rElem(3,4:6),rElem(2,4:6))

! Extract current applied forces/moments at the element nodes.
    call fem_glob2loc_extract (Elem(iElem)%Conn,Force,ForceElem,NumNE)
    ! if (any(isnan(ForceElem))) then
    !     print*, 'Force elem Nan!'
    !     print*, Force
    !     print*, '---'
    !     STOP
    ! end if

! Use full integration for mass matrix.
    NumGaussMass=NumNE

! Compute the element contribution to the mass and damping in the motion of the reference frame.
    call cbeam3_mvel (NumNE,rElem0,rElem,Elem(iElem)%Mass,Mvelelem,NumGaussMass)
    call cbeam3_cvel (NumNE,rElem0,rElem,rElemDot,Vrel,Elem(iElem)%Mass,Cvelelem,Options%NumGauss)
! Contributions of the structural mass to the linearized inertia matrices.
    call cbeam3_mass (NumNE,rElem0,rElem,                                Elem(iElem)%Mass,Melem,NumGaussMass)
    call cbeam3_cgyr (NumNE,rElem0,rElem,rElemDot,Vrel,Elem(iElem)%Mass,Celem,Options%NumGauss)
    call cbeam3_kgyr (NumNE,rElem0,rElem,rElemDot,rElemDDot,Vrel,VrelDot,Elem(iElem)%Mass,Kelem,Options%NumGauss)
    ! if (any(isnan(Kelem))) then
    !     print*, 'Cbeam asbly, line 595'
    !     print*, NumNE
    !     print*, rElem0
    !     print*, rElem
    !     print*, rElemDot
    !     print*, rElemDDot
    !     print*, Vrel
    !     print*, VrelDot
    ! end if

! Compute the gyroscopic force vector.
    call cbeam3_fgyr (NumNE,rElem0,rElem,rElemDot,Vrel,Elem(iElem)%Mass,Qelem,Options%NumGauss)
! Add contributions of non-structural (lumped) mass.
    if (any(abs(Elem(iElem)%RBMass)>0.0d0)) then
      call cbeam3_rbmvel (NumNE,rElem0,rElem,                                Elem(iElem)%RBMass,Mvelelem)
      call cbeam3_rbcvel (NumNE,rElem0,rElem,rElemDot,          Vrel,        Elem(iElem)%RBMass,Cvelelem)
      call cbeam3_rbmass (NumNE,rElem0,rElem,                                Elem(iElem)%RBMass,Melem)
      call cbeam3_rbcgyr (NumNE,rElem0,rElem,rElemDot,          Vrel,        Elem(iElem)%RBMass,Celem)
      call cbeam3_rbkgyr (NumNE,rElem0,rElem,rElemDot,rElemDDot,Vrel,VrelDot,Elem(iElem)%RBMass,Kelem)
    ! if (any(isnan(Kelem))) then
    !     print*, 'Cbeam asbly, line 615'
    !     print*, NumNE
    !     print*, rElem0
    !     print*, rElem
    !     print*, rElemDot
    !     print*, rElemDDot
    !     print*, Vrel
    !     print*, VrelDot
    ! end if
      call cbeam3_rbfgyr (NumNE,rElem0,rElem,rElemDot,          Vrel,        Elem(iElem)%RBMass,Qelem)
    end if
! Compute the element tangent stiffness matrix and force vectors.
    call cbeam3_kmat  (NumNE,rElem0,rElem,Elem(iElem)%Stiff,Kelem,Options%NumGauss)
    ! if (any(isnan(Kelem))) then
    !     print*, 'Cbeam asbly, line 629'
    !     print*, NumNE
    !     print*, rElem0
    !     print*, rElem
    ! end if
    call cbeam3_kgeom (NumNE,rElem0,rElem,Elem(iElem)%Stiff,Kelem,Options%NumGauss)
    ! if (any(isnan(Kelem))) then
    !     print*, 'Cbeam asbly, line 636'
    !     print*, NumNE
    !     print*, rElem0
    !     print*, rElem
    ! end if
    call cbeam3_fstif (NumNE,rElem0,rElem,Elem(iElem)%Stiff,Qelem,Options%NumGauss)

! Add external forces to the residual, and their derivatives with respect to the nodal
! degrees of freedom to the stiffness.
    if (any(abs(ForceElem)>0.0d0)) then
      call cbeam3_fext (NumNE,rElem,Flags(1:NumNE),Felem,Options%FollowerForce,Options%FollowerForceRig,Cao)
    !   print*, 'Elem = ', iElem
      call cbeam3_dqext(NumNE,rElem,ForceElem,Flags(1:NumNE),Kelem,Options%FollowerForce)
        ! if (any(isnan(Kelem))) then
        !     print*, 'Cbeam asbly, line 649'
        !     print*, NumNE
        !     print*, rElem
        !     print*, ForceElem
        !     print*, Flags
        ! end if
    end if

! Project slave degrees of freedom to the orientation of the "master" ones.
    call cbeam3_slave2master (NumNE,Elem(iElem)%Master,rElem0(:,4:6),Psi0,rElem(:,4:6),PsiDefor,SB2B1)
    Melem=matmul(transpose(SB2B1),matmul(Melem,SB2B1))
    Celem=matmul(transpose(SB2B1),matmul(Celem,SB2B1))
    Kelem=matmul(transpose(SB2B1),matmul(Kelem,SB2B1))
    Qelem=matmul(transpose(SB2B1),Qelem)
    Felem=matmul(transpose(SB2B1),Felem)
    Mvelelem=matmul(transpose(SB2B1),Mvelelem)
    Cvelelem=matmul(transpose(SB2B1),Cvelelem)

! Add to global matrix. Remove columns and rows at clamped points.
    do i=1,NumNE
      i1=Node(Elem(iElem)%Conn(i))%Vdof
      if (i1.ne.0) then
        tempi1 = 6*(i1-1)+1
        tempi = 6*(i-1)+1
        Qglobal(tempi1:tempi1+5)   = Qglobal(tempi1:tempi1+5)    + Qelem   (tempi:tempi+5)
        Mvel   (tempi1:tempi1+5,:) = Mvel   (tempi1:tempi1+5,:)  + Mvelelem(tempi:tempi+5,:)
        Cvel   (tempi1:tempi1+5,:) = Cvel   (tempi1:tempi1+5,:)  + Cvelelem(tempi:tempi+5,:)
        do j=1,NumNE
          j1=Node(Elem(iElem)%Conn(j))%Vdof
          if (j1.ne.0) then
            tempi1 = 6*(i1-1)
            tempj1 = 6*(j1-1)
            tempi = 6*(i-1)+1
            tempj = 6*(j-1)+1
            call mat_addmat (tempi1, tempj1, Melem(tempi:tempi+5,tempj:tempj+5),Mglobal)
            call mat_addmat (tempi1, tempj1, Celem(tempi:tempi+5,tempj:tempj+5),Cglobal)
            call mat_addmat (tempi1, tempj1, Kelem(tempi:tempi+5,tempj:tempj+5),Kglobal)
            call mat_addmat (tempi1, tempj1, Felem(tempi:tempi+5,tempj:tempj+5),Fglobal)
          end if
        end do
      end if
    end do
  end do
 end subroutine cbeam3_asbly_dynamic_new_interface
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-> Subroutine CBEAM3_ASBLY_FGLOBAL
!
!-> Description:
!
!   Separate assembly of influence coefficients matrix for
!   applied follower and dead loads
!
!-> Remarks.-
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine cbeam3_asbly_Fglobal (Elem,Node,Coords,Psi0,PosDefor,PsiDefor,Force,           &
&                                 ksf,Kglobal_foll,fsf,Fglobal_foll,fsd,Fglobal_dead,CAG)
  use lib_rotvect
  use lib_fem
  use lib_sparse
  use lib_cbeam3

! I/O variables.
  type(xbelem), intent(in) :: Elem(:)           ! Element information.
  type(xbnode), intent(in) :: Node(:)           ! List of independent nodes.
  real(8),      intent(in) :: Coords    (:,:)   ! Initial coordinates of the grid points.
  real(8),      intent(in) :: Psi0      (:,:,:) ! Initial CRV of the nodes in the elements.
  real(8),      intent(in) :: PosDefor  (:,:)   ! Current coordinates of the grid points
  real(8),      intent(in) :: PsiDefor  (:,:,:) ! Current CRV of the nodes in the elements.
  real(8),      intent(in) :: Force     (:,:)   ! Force vector.
  integer,      intent(out):: ksf               ! Size of the sparse stiffness matrix.
  type(sparse), intent(out):: Kglobal_foll (:)  ! Sparse stiffness matrix.
  integer,      intent(out):: fsf               ! Size of the sparse force matrix, Fglobal_foll.
  type(sparse), intent(out):: Fglobal_foll (:)  ! Influence coefficients matrix for follower forces.
  integer,      intent(out):: fsd               ! Size of the sparse force matrix, Fglobal_dead.
  type(sparse), intent(out):: Fglobal_dead (:)  ! Influence coefficients matrix for dead forces.
  real(8),      intent(in) :: CAG       (:,:)   ! Rotation operator

! Local variables.
  logical:: Flags(MaxElNod)                     ! Auxiliary flags.
  integer:: i,j,i1,j1                           ! Counters.
  integer:: iElem                               ! Counter on the finite elements.
  integer:: NumE                                ! Number of elements in the model.
  integer:: NumNE                               ! Number of nodes in an element.
  real(8):: Kelem_foll (6*MaxElNod,6*MaxElNod)  ! Element tangent stiffness matrix.
  real(8):: Felem_foll (6*MaxElNod,6*MaxElNod)  ! Element force influence coefficients.
  real(8):: Felem_dead (6*MaxElNod,6*MaxElNod)  ! Element force influence coefficients.
  real(8):: rElem0(MaxElNod,6)                  ! Initial Coordinates/CRV of nodes in the element.
  real(8):: rElem (MaxElNod,6)                  ! Current Coordinates/CRV of nodes in the element.
  real(8):: ForceElem (MaxElNod,6)              ! Current forces/moments of nodes in the element.
  real(8):: SB2B1 (6*MaxElNod,6*MaxElNod)       ! Transformation from master to global node orientations.

! Initialise
  call sparse_zero(ksf,Kglobal_foll)
  call sparse_zero(fsf,Fglobal_foll)
  call sparse_zero(fsd,Fglobal_dead)

! Loop in all elements in the model.
  NumE=size(Elem)

  do iElem=1,NumE
    Kelem_foll=0.d0; Felem_foll=0.d0; Felem_dead=0.d0;

    ! Determine if the element nodes are master (Flag=T) or slave.
    Flags=.false.
    do i=1,Elem(iElem)%NumNodes
      if (Node(Elem(iElem)%Conn(i))%Master(1).eq.iElem) Flags(i)=.true.
    end do

    ! Extract components of the displacement and rotation vector at the element nodes
    ! and for the reference and current configurations.
    call fem_glob2loc_extract (Elem(iElem)%Conn,Coords, rElem0(:,1:3),NumNE)
    call fem_glob2loc_extract (Elem(iElem)%Conn,PosDefor,rElem(:,1:3),NumNE)

    rElem0(:,4:6)= Psi0    (iElem,:,:)
    rElem (:,4:6)= PsiDefor(iElem,:,:)
    call rotvect_boundscheck2(rElem(1,4:6),rElem(2,4:6))
    if (NumNE.eq.3) &
&   call rotvect_boundscheck2(rElem(3,4:6),rElem(2,4:6))

    ! Compute the influence coefficients multiplying the vector of external forces.
    call cbeam3_fext (NumNE,rElem,Flags(1:NumNE),Felem_foll,.true._c_bool,.true._c_bool,CAG)
    call cbeam3_fext (NumNE,rElem,Flags(1:NumNE),Felem_dead,.false._c_bool,.false._c_bool,CAG)

    ! Compute contribution to Kglobal from follower forces
    call fem_glob2loc_extract (Elem(iElem)%Conn,Force,ForceElem,NumNE)
    if (any(abs(ForceElem)>0.0d0)) call cbeam3_dqext (NumNE,rElem,ForceElem,Flags(1:NumNE),Kelem_foll,.true._c_bool)

    ! Project equations to the orientation of the "master" degrees of freedom.
    call cbeam3_slave2master (NumNE,Elem(iElem)%Master(:,:),rElem0(:,4:6),Psi0,rElem(:,4:6),PsiDefor,SB2B1)
    Kelem_foll=matmul(transpose(SB2B1),matmul(Kelem_foll,SB2B1))
    Felem_foll=matmul(transpose(SB2B1),Felem_foll)
    Felem_dead=matmul(transpose(SB2B1),Felem_dead)

    ! Add to global matrix. Remove columns and rows at clamped points.
    do i=1,NumNE
      i1=Node(Elem(iElem)%Conn(i))%Vdof
      if (i1.ne.0) then
        do j=1,NumNE
          j1=Node(Elem(iElem)%Conn(j))%Vdof
          if (j1.ne.0) then
            call sparse_addmat (6*(i1-1),6*(j1-1),Kelem_foll(6*(i-1)+1:6*i,6*(j-1)+1:6*j),&
&                               ksf,Kglobal_foll)
            call sparse_addmat (6*(i1-1),6*(j1-1),Felem_foll(6*(i-1)+1:6*i,6*(j-1)+1:6*j),&
&                               fsf,Fglobal_foll)
            call sparse_addmat (6*(i1-1),6*(j1-1),Felem_dead(6*(i-1)+1:6*i,6*(j-1)+1:6*j),&
&                               fsd,Fglobal_dead)
          end if
        end do
      end if
    end do

  end do

  return
 end subroutine cbeam3_asbly_Fglobal

 function cbeam3_asbly_gravity_static(NumDof, options, g) result(gravity_vec)
    use lib_rot
    integer, intent(IN)             :: NumDof
    type(xbopts), intent(IN)        :: options
    real(8), optional, intent(IN)   :: g(3)
    real(8)                         :: gravity_vec(numdof)

    ! local variables
    integer                         :: i
    integer                         :: ii
    real(8)                         :: g_2(3)

    if (present(g) .eqv. .FALSE.) then
        g_2 = rot_unit([options%gravity_dir_x,&
                        options%gravity_dir_y,&
                        options%gravity_dir_z])*options%gravity
        ! print*, 'g static', g_2
    else
        g_2 = g
    end if

    gravity_vec = 0.0d0
    do i=1,NumDof/6
        ii = (i - 1)*6
        gravity_vec(ii+1: ii+3) = g_2
    end do
 end function cbeam3_asbly_gravity_static

 function cbeam3_asbly_gravity_dynamic(NumDof, options, Cao) result(gravity_vec)
    use lib_rot
    integer, intent(IN)         :: NumDof
    type(xbopts), intent(IN)    :: options
    real(8), intent(IN)         :: Cao(3,3)
    real(8)                     :: gravity_vec(numdof)

    ! local vars
    real(8)                     :: g(3)

    g = rot_unit([options%gravity_dir_x,&
                  options%gravity_dir_y,&
                  options%gravity_dir_z])*options%gravity
    g = MATMUL(Cao, g)
    ! print*, 'g dynamic ', g
    gravity_vec = cbeam3_asbly_gravity_static(NumDof, options, g)
 end function cbeam3_asbly_gravity_dynamic

end module cbeam3_asbly
