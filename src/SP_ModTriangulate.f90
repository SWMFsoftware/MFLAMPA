!  Copyright (C) 2002 Regents of the University of Michigan,
!  portions used with permission
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module SP_ModTriangulate

  use ModUtilities, ONLY: CON_stop
#ifdef _OPENACC
  use ModUtilities, ONLY: norm2
#endif
  use SP_ModGrid, ONLY: nLine, nLineAll, Used_B, iLineAll0, &
       search_line, MHData_VIB, X_, Y_, Z_, nP, nMu
  use SP_ModProc, ONLY: iProc, nProc, iComm, iError
  use SP_ModSize, ONLY: nDim
  use SP_ModTurbulence, ONLY: DPerp_CB, DPerp_IIBI

  implicit none

  SAVE

  PRIVATE ! Except

  public:: read_param           ! Read module parameters
  public:: reset_intersect_surf ! Reset arrays used in intersect_surf
  public:: intersect_surf       ! Intersect on a sphere with a given radius
  public:: build_trmesh         ! Triangulation on a sphere
  public:: interpolate_trmesh   ! Interpolate at the specified location

  ! In all names below "Tri" means "Triangulation"
  ! If we use poles in triangulation
  logical, public :: UsePoleTri   = .true.
  logical, public :: UsePlanarTri = .true.
  real,    public, parameter :: iSouthPoleTri_ = -1, iNorthPoleTri_ = -2
  ! # of intersection points of the MFLAMPA grid lines with
  ! the spherical surface, nReachR
  ! # of points involed into triangilation (equals to nReachR or nReach+2)
  integer, public, allocatable :: nTriMesh_I(:)
  ! Locations of the lines intersection with the spherical surface. Regardless
  ! of the sphere radius, the coordinates are divided by this radius, hence,
  ! are the points at the unit sphere. The point number may be equal to
  ! nLineAll, the total number of lines, or two extra points at the poles
  ! may be added. If some lines are unused or do not reach the surface, the
  ! last elements of the array are not used
  real, allocatable :: Xyz_DII(:,:,:)
  real, public, allocatable :: XyzOrig_DII(:,:,:)
  ! Index of the line correspoinding to a given point at the radial surface
  ! It may differ from just the array index because of the presence of polar
  ! points or absence of bad lines not reaching the boundary.
  integer, allocatable :: iLineReach_II(:,:)
  ! Distribution function interpolated to the points at the spherical surface
  real, allocatable    :: DistrR_IIBI(:,:,:,:)
  ! The result of triangulation
  integer, allocatable :: iList_I(:), iPointer_I(:), iEnd_I(:)

  ! Triangulation used by build_trmesh (trmesh_sp) and
  ! interpolate_trmesh (interpolate_trmesh_sp): the nodes are
  ! normalized, nodes closer than TolMergeTri are merged (a crossing at a
  ! pole replaces the added pole), a collinear first triple is reordered,
  ! a great circle is closed, trmesh and the mesh are checked, and the
  ! interpolation weights are finite, >= 0 and sum to 1.
  ! Nodes closer than this angle [rad] are merged into one node
  real :: TolMergeTri = 1.0e-6
  ! The mesh of one shell (like iList_I it is built by build_trmesh
  ! for shell iRMeshSp and used by the next interpolate_trmesh calls):
  ! node coordinates, line index of each node (or -1/-2 for the poles,
  ! -3/-4 for auxiliary nodes), STRIPACK structure, states at the nodes.
  ! Xyz_DII, DistrR_IIBI, DPerp_IIBI and iLineReach_II are not changed.
  integer :: iRMeshSp = -1, nMeshSp = 0
  logical :: IsDPerpMeshSp = .false.
  real,    allocatable :: XyzMeshSp_DI(:,:)
  integer, allocatable :: iCodeMeshSp_I(:)
  integer, allocatable :: iListSp_I(:), iPointerSp_I(:), &
       iEndSp_I(:)
  real,    allocatable :: DistrMeshSp_III(:,:,:), &
       DPerpMeshSp_III(:,:,:)
  ! Input node (index in Xyz_DII) kept as mesh node, 0 for auxiliary nodes
  integer, allocatable :: iSourceSp_I(:)
  ! Notes about merged nodes and fallbacks: at most MaxNoteSp lines
  integer :: nNoteSp = 0
  integer, parameter :: MaxNoteSp = 10

  ! Status codes of build_sp
  integer, parameter :: BuildOk_ = 0, BuildBadNode_ = 1, &
       BuildFewNodes_ = 2, BuildTrmesh_ = 3, BuildMesh_ = 4, &
       BuildCapacity_ = 5
  ! Status codes of locate_sp
  integer, parameter :: LocateInside_ = 0, LocateDegenerate_ = 1, &
       LocateOutside_ = 2, LocateBadPoint_ = -1
  ! Codes of the auxiliary nodes added by build_sp (the poles added by
  ! intersect_surf have the codes -1 and -2)
  integer, parameter :: iAux1Sp_ = -3, iAux2Sp_ = -4

  ! A node is renormalized only if ||x|-1| > TolNormSp, so a second
  ! call does not change a node that is already normalized
  real, parameter :: TolNormSp = 1e-14
  ! A triangle with det(V1,V2,V3) <= DetMinSp is degenerate. The
  ! rounding error of det for unit vectors is a few 1e-16.
  real, parameter :: DetMinSp = 1e-14
  ! Smallest threshold used for "nearly on one great circle" [rad]
  real, parameter :: TolHeightMinSp = 1e-12

  ! Statistics of one build_sp call
  type :: BuildInfoType
     integer :: nNodeIn = 0       ! nodes passed in
     integer :: nMesh = 0         ! nodes of the mesh (merged, + auxiliary)
     integer :: nMerged = 0       ! data nodes merged into another node
     integer :: nPoleReplaced = 0 ! added nodes replaced by a data node
     integer :: nNormalized = 0   ! nodes rescaled to unit length
     logical :: IsReordered = .false. ! first three nodes changed
     logical :: IsGreatCircle = .false. ! all nodes on one great circle
     integer :: nAux = 0          ! auxiliary nodes added
     integer :: nTrmesh = 0       ! number of trmesh calls
     integer :: iErrorTrmesh = 0  ! last trmesh error code
     integer :: nTri = 0          ! triangles
     integer :: nBoundary = 0     ! boundary nodes (0: closed mesh)
     integer :: nDegenerate = 0   ! triangles with det <= DetMinSp
     real    :: MaxNormDev = 0.0  ! max ||x|-1| of the input nodes
     real    :: MinDet = 0.0      ! min det(V1,V2,V3) of the triangles
     real    :: MinHeight3 = 0.0  ! smallest height of nodes 1-3 [rad]
  end type BuildInfoType
contains
  !============================================================================
  subroutine read_param(NameCommand)

    use ModReadParam, ONLY: read_var
    character(len=*), intent(in) :: NameCommand ! From PARAM.in
    character(len=*), parameter:: NameSub = 'read_param'
    !--------------------------------------------------------------------------
    select case(NameCommand)
    case('#TRIANGULATION')
       ! get pole triangulartion flag
       call read_var('UsePoleTriangulation', UsePoleTri)
       ! get the triangulation approach flag
       call read_var('UsePlanarTriangles', UsePlanarTri)
       ! nodes closer than this angle [rad] are merged into one node
       call read_var('TolMergeTriangulation', TolMergeTri)
       if(TolMergeTri < 0.0) call CON_stop(NameSub// &
            ': TolMergeTriangulation must not be negative')
    case default
       call CON_stop(NameSub//': Unknown command '//NameCommand)
    end select

  end subroutine read_param
  !============================================================================
  subroutine reset_intersect_surf(nR)
    ! Since there can be multiple calls (`ModDiffusion` or `ModPlot` or
    ! `ModSatellite`) of `intersect_surf`, another subroutine should be
    ! called from ModDiffusion to reset array bounds:
    integer, intent(in) :: nR ! total # of spherical coordinate surfaces
    !  integer, intent(in) :: iRfirst, iRlast ! First and last # of surfaces
    ! assigned to a given Proc
    character(len=*), parameter:: NameSub = 'reset_intersect_surf'
    !--------------------------------------------------------------------------
    if(allocated(Xyz_DII)) deallocate(Xyz_DII, iLineReach_II, &
         DistrR_IIBI, nTriMesh_I, iList_I, iPointer_I, iEnd_I)
    if(allocated(XyzOrig_DII)) deallocate(XyzOrig_DII, DPerp_IIBI)
    allocate(Xyz_DII(X_:Z_, 1:nLineAll+2, nR))
    allocate(XyzOrig_DII(X_:Z_, 1:nLineAll+2, nR))
    allocate(DistrR_IIBI(0:nP+1, 1:nMu, 1:nLineAll+2, nR))
    allocate(DPerp_IIBI(1:nP, 1:nMu, 1:nLineAll+2, nR))
    allocate(iLineReach_II(1:nLineAll+2, nR))
    allocate(nTriMesh_I(nR)); nTriMesh_I = 0
    allocate(iList_I(6*nLineAll))
    allocate(iPointer_I(6*nLineAll))
    allocate(iEnd_I(nLineAll+2))
  end subroutine reset_intersect_surf
  !============================================================================
  subroutine intersect_surf(rSurf, iRootIn, iRIn, IsIncludeDPerpIn)

    use ModMpi
    use SP_ModDistribution, ONLY: Distribution_CB
    ! Input: radius of the spherical surface
    real,    intent(in) :: rSurf
    integer, optional, intent(in) :: iRootIn, iRIn
    logical, optional, intent(in) :: IsIncludeDPerpIn

    ! called from Dperp (=> true) or ModPlot (=> false, default)
    logical :: IsIncludeDPerp
    ! loop variables
    integer :: iLine, iReachR
    integer :: iRoot, iR, nReachR
    ! indexes of corresponding node, latitude and longitude
    integer :: iLineAll
    ! index of particle just above the radius
    integer :: iAboveR
    ! interpolation weight (in radial direction)
    real    :: WeightR
    ! skip a field line not reaching radius of output sphere
    logical :: DoReachR
    character(len=*), parameter:: NameSub = 'intersect_surf'
    !--------------------------------------------------------------------------
    ! Check whether we need to include Dperp
    if(present(IsIncludeDPerpIn)) then
       IsIncludeDPerp = IsIncludeDPerpIn
    else
       IsIncludeDPerp = .false.
    end if

    ! For full implementation of uniform spherical grid for solving
    ! perpendicular diffusion and drift, the lines below should be
    ! shaped as subroutine init
    if(.not.allocated(Xyz_DII)) then
       allocate(Xyz_DII(X_:Z_, 1:nLineAll+2, 1))
       allocate(DistrR_IIBI(0:nP+1, 1:nMu, 1:nLineAll+2, 1))
       if(IsIncludeDPerp) then
          ! In DPerp
          allocate(XyzOrig_DII(X_:Z_, 1:nLineAll+2, 1))
          allocate(DPerp_IIBI(1:nP, 1:nMu, 1:nLineAll+2, 1))
       end if
       allocate(iLineReach_II(1:nLineAll+2, 1))
       allocate(nTriMesh_I(1)); nTriMesh_I = 0
       allocate(iList_I(6*nLineAll))
       allocate(iPointer_I(6*nLineAll))
       allocate(iEnd_I(nLineAll+2))
    end if
    if(present(iRootIn)) then
       iRoot = iRootIn
    else
       iRoot = 0
    end if
    if(present(iRin)) then
       iR = iRin
    else
       iR = 1
    end if
    Xyz_DII(:,:,iR) = 0.0; iLineReach_II(:,iR) = 0
    DistrR_IIBI(:,:,:,iR) = 0.0
    if(IsIncludeDPerp) DPerp_IIBI(:,:,:,iR) = 0.0

    ! go over all lines on the processor and find the point of
    ! intersection with output sphere if present
    LINE:do iLine = 1, nLine
       if(.not.Used_B(iLine)) CYCLE LINE
       iLineAll = iLineAll0 + iLine

       ! Find the particle just above the given radius
       call search_line(iLine, rSurf, iAboveR, DoReachR, WeightR)
       DoReachR = DoReachR .and. iAboveR /= 1
       ! if no intersection found -> proceed to the next line
       if(.not.DoReachR) CYCLE LINE

       ! Found intersection => get POS & log10(f[IO_Unit]) at that location
       ! Find coordinates and log(Distribution) at intersection
       Xyz_DII(:, iLineAll, iR) = ( &
            MHData_VIB(X_:Z_, iAboveR-1, iLine)*(1-WeightR) +  &
            MHData_VIB(X_:Z_, iAboveR,   iLine)*   WeightR )
       DistrR_IIBI(:, :, iLineAll, iR) = ( &
            Distribution_CB(:, :, iAboveR,   iLine)* WeightR + &
            Distribution_CB(:, :, iAboveR-1, iLine)*(1-WeightR))
       if(IsIncludeDPerp) DPerp_IIBI(:, :, iLineAll, iR) = ( &
            DPerp_CB(:, :, iAboveR,   iLine)* WeightR + &
            DPerp_CB(:, :, iAboveR-1, iLine)*(1-WeightR))
    end do LINE !  iLine

    ! Gather interpolated coordinates on the source processor and Broadcast
    if(nProc > 1) then
       ! The entire coordinate array is sent to iRoot, however, the unused
       ! line intersections are marked with all coordinates set to zero.
       call mpi_reduce_real_array(Xyz_DII(:,:,iR), 3*nLineAll, MPI_SUM, &
            iRoot, iComm, iError)
       call mpi_reduce_real_array(DistrR_IIBI(:,:,:,iR), &
            (nP+2)*nMu*nLineAll, MPI_SUM, iRoot, iComm, iError)
       if(IsIncludeDPerp) call mpi_reduce_real_array(DPerp_IIBI(:,:,:,iR), &
            nP*nMu*nLineAll, MPI_SUM, iRoot, iComm, iError)
    end if
    ! Save original intersection points
    if(IsIncludeDPerp) XyzOrig_DII(:,:,iR) = Xyz_DII(:,:,iR)
    if(iProc /= iRoot) RETURN

    ! Sort out unused points
    iReachR = 0
    do iLineAll = 1, nLineAll
       if(all(Xyz_DII(:,iLineAll,iR)==0.0)) CYCLE
       iReachR = iReachR + 1
       ! Store and normalize coordinates to get points on the unit square
       Xyz_DII(:,iReachR,iR) = Xyz_DII(:,iLineAll,iR)/rSurf
       iLineReach_II(iReachR,iR) = iLineAll
       DistrR_IIBI(:,:,iReachR,iR) = DistrR_IIBI(:,:,iLineAll,iR)
       if(IsIncludeDPerp) &
            DPerp_IIBI(:,:,iReachR,iR) = DPerp_IIBI(:,:,iLineAll,iR)
    end do
    if(iReachR == 0) call CON_stop(NameSub//': No valid intersection points!')
    nReachR = iReachR

    ! For poles
    if(UsePoleTri) then
       ! Add two grid nodes at the poles:
       nTriMesh_I(iR) = nReachR + 2
       Xyz_DII(:, 2:nTriMesh_I(iR)-1,iR)    = Xyz_DII(:, 1:nReachR,iR)
       Xyz_DII(:, 1, iR)                    = [0.0, 0.0, -1.0]
       Xyz_DII(:, nTriMesh_I(iR), iR)       = [0.0, 0.0, +1.0]
       iLineReach_II(2:nTriMesh_I(iR)-1,iR) = iLineReach_II(1:nReachR,iR)
       iLineReach_II(1,iR)                  = iSouthPoleTri_
       iLineReach_II(nTriMesh_I(iR),iR)     = iNorthPoleTri_
       DistrR_IIBI(:,:,2:nTriMesh_I(iR)-1,iR) = DistrR_IIBI(:,:,1:nReachR,iR)
       ! Shift DPerp as DistrR, so node i carries the Dperp of its own line;
       ! pole values are set by fix_state in build_trmesh(IsIncludeDPerpIn)
       if(IsIncludeDPerp) DPerp_IIBI(:,:,2:nTriMesh_I(iR)-1,iR) = &
            DPerp_IIBI(:,:,1:nReachR,iR)
    else
       nTriMesh_I(iR)= nReachR
    end if
  end subroutine intersect_surf
  !============================================================================
  subroutine build_trmesh(iRIn, IsIncludeDPerpIn)

    use ModTriangulateSpherical, ONLY: trmesh, fix_state
    integer, optional, intent(in) :: iRIn
    logical, optional, intent(in) :: IsIncludeDPerpIn

    ! called from Dperp (=> true) or ModPlot (=> false, default)
    logical :: IsIncludeDPerp
    integer :: iR, iMu
    character(len=*), parameter:: NameSub = 'build_trmesh'
    !--------------------------------------------------------------------------
    ! Check whether we need to include Dperp
    if(present(IsIncludeDPerpIn)) then
       IsIncludeDPerp = IsIncludeDPerpIn
    else
       IsIncludeDPerp = .false.
    end if
    if(present(iRIn)) then
       iR = iRIn
    else
       iR = 1
    end if
    iList_I = 0; iPointer_I = 0; iEnd_I = 0

    ! Construct the Triangular mesh used for interpolation
    call trmesh_sp(iR, IsIncludeDPerp, nTriMesh_I(iR), &
         Xyz_DII(X_, 1:nTriMesh_I(iR), iR), &
         Xyz_DII(Y_, 1:nTriMesh_I(iR), iR), &
         Xyz_DII(Z_, 1:nTriMesh_I(iR), iR), &
         iList_I(:6*(nTriMesh_I(iR)-2)),    &
         iPointer_I(:6*(nTriMesh_I(iR)-2)), &
         iEnd_I(:nTriMesh_I(iR)), iError)
    ! trmesh_sp has set the states at the nodes of its mesh (the poles
    ! too), so the fix_state calls below for the old mesh are not used
    RETURN
    ! Fix states (for the given function) at the polar nodes
    if(UsePoleTri) then
       do iMu = 1, nMu
          ! North:
          call fix_state( &
               iNodeToFix = nTriMesh_I(iR),                    &
               nNode      = nTriMesh_I(iR),                    &
               iList_I    = iList_I(:6*(nTriMesh_I(iR)-2)),    &
               iPointer_I = iPointer_I(:6*(nTriMesh_I(iR)-2)), &
               iEnd_I     = iEnd_I(:nTriMesh_I(iR)),           &
               Xyz_DI     = Xyz_DII(:,1:nTriMesh_I(iR),iR),    &
               nVar       = nP+2,                              &
               State_VI   = DistrR_IIBI(:,iMu,1:nTriMesh_I(iR),iR))
          ! South:
          call fix_state( &
               iNodeToFix = 1,                                 &
               nNode      = nTriMesh_I(iR),                    &
               iList_I    = iList_I(:6*(nTriMesh_I(iR)-2)),    &
               iPointer_I = iPointer_I(:6*(nTriMesh_I(iR)-2)), &
               iEnd_I     = iEnd_I(:nTriMesh_I(iR)),           &
               Xyz_DI     = Xyz_DII(:,1:nTriMesh_I(iR),iR),    &
               nVar       = nP+2,                              &
               State_VI   = DistrR_IIBI(:, iMu, 1:nTriMesh_I(iR),iR))
       end do
       if(IsIncludeDPerp) then
          do iMu = 1, nMu
             ! North:
             call fix_state( &
                  iNodeToFix = nTriMesh_I(iR),                    &
                  nNode      = nTriMesh_I(iR),                    &
                  iList_I    = iList_I(:6*(nTriMesh_I(iR)-2)),    &
                  iPointer_I = iPointer_I(:6*(nTriMesh_I(iR)-2)), &
                  iEnd_I     = iEnd_I(:nTriMesh_I(iR)),           &
                  Xyz_DI     = Xyz_DII(:,1:nTriMesh_I(iR),iR),    &
                  nVar       = nP,                                &
                  State_VI   = DPerp_IIBI(:,iMu,1:nTriMesh_I(iR),iR))
             ! South:
             call fix_state( &
                  iNodeToFix = 1,                                 &
                  nNode      = nTriMesh_I(iR),                    &
                  iList_I    = iList_I(:6*(nTriMesh_I(iR)-2)),    &
                  iPointer_I = iPointer_I(:6*(nTriMesh_I(iR)-2)), &
                  iEnd_I     = iEnd_I(:nTriMesh_I(iR)),           &
                  Xyz_DI     = Xyz_DII(:,1:nTriMesh_I(iR),iR),    &
                  nVar       = nP,                                &
                  State_VI   = DPerp_IIBI(:,iMu, 1:nTriMesh_I(iR),iR))
          end do
       end if
    end if
  end subroutine build_trmesh
  !============================================================================
  subroutine interpolate_trmesh(XyzInterp_D, iRIn, &
       DistrInterp_II, DPerp_II, iStencilOut_I, WeightOut_I)

    use ModTriangulateSpherical, ONLY: find_triangle_orig, find_triangle_sph
    ! Input: Xyz coordinates of the point for interpolations
    real,    intent(in)   :: XyzInterp_D(nDim)
    integer, optional, intent(in) :: iRIn
    ! Output: Interpolated values at XyzInterp_D
    real,    optional, intent(out):: DistrInterp_II(0:nP+1,1:nMu)
    real,    optional, intent(out):: DPerp_II(1:nP,1:nMu)
    ! variables for interpolation in a triangular mesh
    integer, optional, intent(out):: iStencilOut_I(3)
    real,    optional, intent(out):: WeightOut_I(3)

    ! local variables
    logical :: IsTriangleFound
    integer :: iStencil_I(3)
    real    :: Weight_I(3)
    ! loop variables
    integer :: i, iR
    character(len=*), parameter:: NameSub = 'interpolate_trmesh'
    !--------------------------------------------------------------------------
    if(present(iRIn)) then
       iR = iRIn
    else
       iR = 1
    end if
    ! Location and weights on the mesh of trmesh_sp, so the old
    ! find_triangle_orig/find_triangle_sph lookup below is not used
    call interpolate_trmesh_sp(XyzInterp_D, iR, DistrInterp_II, &
         DPerp_II, iStencilOut_I, WeightOut_I)
    RETURN
    IsTriangleFound = .false.
    if(present(iStencilOut_I)) iStencilOut_I = -1
    if(present(WeightOut_I)) WeightOut_I = 0.0
    if(present(DistrInterp_II)) DistrInterp_II = 0.0
    if(present(DPerp_II)) DPerp_II = 0.0

    ! Find the triangle where the satellite locates
    if(UsePlanarTri) then
       call find_triangle_orig(               &
            XyzInterp_D/norm2(XyzInterp_D),   &
            nTriMesh_I(iR),                   &
            Xyz_DII(:,1:nTriMesh_I(iR),iR),   &
            iList_I(:6*(nTriMesh_I(iR)-2)),   &
            iPointer_I(:6*(nTriMesh_I(iR)-2)),&
            iEnd_I(:nTriMesh_I(iR)),          &
            Weight_I, IsTriangleFound, iStencil_I)
    else
       call find_triangle_sph(                &
            XyzInterp_D/norm2(XyzInterp_D),   &
            nTriMesh_I(iR),                   &
            Xyz_DII(:,1:nTriMesh_I(iR),iR),   &
            iList_I(:6*(nTriMesh_I(iR)-2)),   &
            iPointer_I(:6*(nTriMesh_I(iR)-2)),&
            iEnd_I(:nTriMesh_I(iR)),          &
            Weight_I(1), Weight_I(2),         &
            Weight_I(3), IsTriangleFound,     &
            iStencil_I(1), iStencil_I(2), iStencil_I(3))
    end if
    if(.not.IsTriangleFound) RETURN
    Weight_I = max(Weight_I, 0.0)

    ! For outputs as requested in the interface
    if(present(DistrInterp_II)) then
       ! Interpolate the distribution function at satellite as outputs
       do i = 1, 3
          DistrInterp_II = DistrInterp_II + &
               DistrR_IIBI(:,:,iStencil_I(i),iR)*Weight_I(i)
       end do
    end if
    if(present(DPerp_II)) then
       ! Interpolate the Dperp coefficient as outputs
       do i = 1, 3
          DPerp_II = DPerp_II + &
               DPerp_IIBI(:,:,iStencil_I(i),iR)*Weight_I(i)
       end do
    end if
    if(present(iStencilOut_I)) then
       do i = 1, 3
          ! Recover iLineAll corresponding to iStencil
          iStencilOut_I(i) = iLineReach_II(iStencil_I(i),iR)
       end do
    end if
    if(present(WeightOut_I)) WeightOut_I = Weight_I
  end subroutine interpolate_trmesh
  !============================================================================
  subroutine trmesh_sp(iR, IsIncludeDPerp, n, x_I, y_I, z_I, &
       iList_I, iPointer_I, iEnd_I, iError)

    ! Replaces "call trmesh(n, x, y, z, list, lptr, lend, iError)" in
    ! build_trmesh: the same arguments, after the shell index iR and
    ! IsIncludeDPerp. The n nodes x_I, y_I, z_I of shell iR are normalized
    ! and merged (TolMergeTri), build_sp calls trmesh, checks iErrorTrmesh
    ! and the mesh. The mesh, the line index (code) and the state of every mesh
    ! node go to the *Sp arrays; the input arrays of intersect_surf are
    ! not changed, so repeated calls for the same shell give the same mesh.
    ! States: a merged node gets the average of its data nodes (DistrR and,
    ! with IsIncludeDPerp, DPerp) and keeps the first line index; the added
    ! poles that are kept and the auxiliary nodes get the inverse-distance
    ! average of their neighbors (fix_state_sp, as fix_state); a pole
    ! replaced by a crossing is not in the mesh.
    ! iList_I, iPointer_I, iEnd_I (the mesh of the n input nodes made by
    ! trmesh) are set to 0: the mesh may have another number of
    ! nodes and is stored in iListSp_I, iPointerSp_I, iEndSp_I.
    ! iError = 0 on return, as from trmesh without error. If build_sp
    ! fails, the run stops with the shell, the node count and the reason.
    integer, intent(in) :: iR
    logical, intent(in) :: IsIncludeDPerp
    integer, intent(in) :: n
    real,    intent(in) :: x_I(n), y_I(n), z_I(n)
    integer, intent(out):: iList_I(6*(n-2)), iPointer_I(6*(n-2)), iEnd_I(n)
    integer, intent(out):: iError

    type(BuildInfoType) :: Info
    integer :: nMax, nMesh, iNode, iMesh, iMu, k, nSpecial
    integer :: iSpecial_I(4), iOrder_I(4)
    integer, allocatable :: iSource_I(:), iGroup_I(:), nGroup_I(:)
    logical, allocatable :: IsSpecial_I(:)
    character(len=300) :: StringError
    character(len=500) :: String
    character(len=*), parameter:: NameSub = 'trmesh_sp'
    !--------------------------------------------------------------------------
    iList_I = 0; iPointer_I = 0; iEnd_I = 0
    nMax = nLineAll + 4 ! room for 2 auxiliary nodes
    if(n < 1 .or. n > nLineAll + 2)then
       write(String,'(a,i0,a,i0)') NameSub//': shell iR=', iR, &
            ' has a bad number of nodes ', n
       call CON_stop(trim(String))
    end if
    if(allocated(XyzMeshSp_DI))then
       if(size(XyzMeshSp_DI, 2) /= nMax .or. &
            size(DistrMeshSp_III, 1) /= nP + 2 .or. &
            size(DistrMeshSp_III, 2) /= nMu) &
            deallocate(XyzMeshSp_DI, iCodeMeshSp_I, iListSp_I, &
            iPointerSp_I, iEndSp_I, DistrMeshSp_III, &
            iSourceSp_I)
    end if
    if(.not.allocated(XyzMeshSp_DI)) allocate( &
         XyzMeshSp_DI(3,nMax), iCodeMeshSp_I(nMax), &
         iListSp_I(6*(nMax-2)), iPointerSp_I(6*(nMax-2)), &
         iEndSp_I(nMax), DistrMeshSp_III(0:nP+1,nMu,nMax), &
         iSourceSp_I(nMax))
    if(IsIncludeDPerp)then
       if(allocated(DPerpMeshSp_III))then
          if(size(DPerpMeshSp_III, 3) /= nMax .or. &
               size(DPerpMeshSp_III, 1) /= nP .or. &
               size(DPerpMeshSp_III, 2) /= nMu) &
               deallocate(DPerpMeshSp_III)
       end if
       if(.not.allocated(DPerpMeshSp_III)) &
            allocate(DPerpMeshSp_III(nP,nMu,nMax))
    end if
    allocate(iSource_I(nMax), iGroup_I(n), nGroup_I(nMax), IsSpecial_I(nMax))

    ! Triangulate a copy of the nodes
    iRMeshSp = -1
    XyzMeshSp_DI(1,1:n) = x_I
    XyzMeshSp_DI(2,1:n) = y_I
    XyzMeshSp_DI(3,1:n) = z_I
    iCodeMeshSp_I(1:n) = iLineReach_II(1:n,iR)
    call build_sp(n, nMax, XyzMeshSp_DI, iCodeMeshSp_I, &
         TolMergeTri, nMesh, iSource_I, iGroup_I, iListSp_I, &
         iPointerSp_I, iEndSp_I, Info, iError, StringError)
    if(iError /= BuildOk_)then
       write(String,'(a,i0,a,i0,2a)') NameSub//': shell iR=', iR, ', ', n, &
            ' nodes (poles included): ', trim(StringError)
       call CON_stop(trim(String))
    end if
    nMeshSp = nMesh
    iSourceSp_I(1:nMesh) = iSource_I(1:nMesh)
    if((Info%nMerged > 0 .or. Info%nPoleReplaced > 0 .or. Info%IsReordered &
         .or. Info%nAux > 0) .and. nNoteSp < MaxNoteSp)then
       nNoteSp = nNoteSp + 1
       write(*,'(a,i0,a,i0,a,i0,a,i0,a,l1,a,i0,a,i0)') NameSub// &
            ': shell iR=', iR, ' nodes=', n, ' merged=', Info%nMerged, &
            ' poles replaced=', Info%nPoleReplaced, ' reordered=', &
            Info%IsReordered, ' auxiliary nodes=', Info%nAux, ' mesh nodes=', &
            nMesh
       if(nNoteSp == MaxNoteSp) write(*,*) NameSub// &
            ': further notes of SP_ModTriangulate are suppressed'
    end if

    ! States at the mesh nodes: average over the data nodes of each group
    nGroup_I(1:nMesh) = 0
    DistrMeshSp_III(:,:,1:nMesh) = 0.0
    if(IsIncludeDPerp) DPerpMeshSp_III(:,:,1:nMesh) = 0.0
    do iNode = 1, n
       ! the added poles are not data, their state is set below
       if(iLineReach_II(iNode,iR) < 0) CYCLE
       iMesh = iGroup_I(iNode)
       nGroup_I(iMesh) = nGroup_I(iMesh) + 1
       DistrMeshSp_III(:,:,iMesh) = DistrMeshSp_III(:,:,iMesh) &
            + DistrR_IIBI(:,:,iNode,iR)
       if(IsIncludeDPerp) DPerpMeshSp_III(:,:,iMesh) = &
            DPerpMeshSp_III(:,:,iMesh) + DPerp_IIBI(:,:,iNode,iR)
    end do
    do iMesh = 1, nMesh
       if(nGroup_I(iMesh) < 2) CYCLE
       DistrMeshSp_III(:,:,iMesh) = &
            DistrMeshSp_III(:,:,iMesh)/nGroup_I(iMesh)
       if(IsIncludeDPerp) DPerpMeshSp_III(:,:,iMesh) = &
            DPerpMeshSp_III(:,:,iMesh)/nGroup_I(iMesh)
    end do

    ! Kept poles and auxiliary nodes (north first, as in build_trmesh)
    IsSpecial_I(1:nMesh) = iCodeMeshSp_I(1:nMesh) < 0
    nSpecial = 0
    iOrder_I = [nint(iNorthPoleTri_), nint(iSouthPoleTri_), iAux1Sp_, &
         iAux2Sp_]
    do k = 1, 4
       do iMesh = 1, nMesh
          if(iCodeMeshSp_I(iMesh) == iOrder_I(k))then
             nSpecial = nSpecial + 1
             iSpecial_I(nSpecial) = iMesh
          end if
       end do
    end do
    do iMu = 1, nMu
       do k = 1, nSpecial
          call fix_state_sp(iSpecial_I(k), nMesh, &
               iListSp_I(:6*(nMesh-2)), iPointerSp_I(:6*(nMesh-2)), &
               iEndSp_I(:nMesh), XyzMeshSp_DI(:,1:nMesh), &
               IsSpecial_I(1:nMesh), nP+2, DistrMeshSp_III(:,iMu,1:nMesh))
          if(IsIncludeDPerp) call fix_state_sp(iSpecial_I(k), nMesh, &
               iListSp_I(:6*(nMesh-2)), iPointerSp_I(:6*(nMesh-2)), &
               iEndSp_I(:nMesh), XyzMeshSp_DI(:,1:nMesh), &
               IsSpecial_I(1:nMesh), nP, DPerpMeshSp_III(:,iMu,1:nMesh))
       end do
    end do
    iRMeshSp = iR
    IsDPerpMeshSp = IsIncludeDPerp
    deallocate(iSource_I, iGroup_I, nGroup_I, IsSpecial_I)
  end subroutine trmesh_sp
  !============================================================================
  subroutine interpolate_trmesh_sp(XyzInterp_D, iR, &
       DistrInterp_II, DPerp_II, iStencilOut_I, WeightOut_I)

    ! Replaces the lookup of interpolate_trmesh (same arguments, iR is
    ! required) on the mesh of trmesh_sp. The weights are finite, >= 0
    ! and sum to 1 (locate_sp). Fallbacks, reported on the screen:
    ! - ill-conditioned triangle: the state of its nearest vertex;
    ! - point outside the triangulated region (only for a mesh with a
    !   boundary, i.e. without the pole nodes and with lines that do not
    !   surround the point): the state of the nearest node, with
    !   iStencilOut_I = -1 and WeightOut_I = 0 as for "not found" before,
    !   so ModPlot and ModSatellite keep their treatment of such points,
    !   while ModDiffusionPerp and ModDrift get the nearest state, not 0.
    ! DPerp_II after a build_trmesh without IsIncludeDPerpIn (as
    ! ModDiffusionPerp does): DPerp_IIBI of the kept input node, as
    ! intersect_surf left it (not averaged, poles not fixed), like the old
    ! lookup; an auxiliary node, which has no input node, adds 0.
    real,    intent(in) :: XyzInterp_D(nDim)
    integer, intent(in) :: iR
    real,    optional, intent(out):: DistrInterp_II(0:nP+1,1:nMu)
    real,    optional, intent(out):: DPerp_II(1:nP,1:nMu)
    integer, optional, intent(out):: iStencilOut_I(3)
    real,    optional, intent(out):: WeightOut_I(3)

    integer :: iStencil_I(3), iStatus, i
    real    :: Weight_I(3)
    character(len=300) :: String
    character(len=*), parameter:: NameSub = 'interpolate_trmesh_sp'
    !--------------------------------------------------------------------------
    if(iR /= iRMeshSp)then
       write(String,'(a,i0,a)') NameSub//': no mesh for shell iR=', &
            iR, ', call build_trmesh first'
       call CON_stop(trim(String))
    end if
    call locate_sp(XyzInterp_D, nMeshSp, &
         XyzMeshSp_DI(:,1:nMeshSp), &
         iListSp_I(:6*(nMeshSp-2)), &
         iPointerSp_I(:6*(nMeshSp-2)), iEndSp_I(:nMeshSp), &
         UsePlanarTri, iStencil_I, Weight_I, iStatus)
    select case(iStatus)
    case(LocateInside_)
    case(LocateDegenerate_, LocateOutside_)
       if(nNoteSp < MaxNoteSp)then
          nNoteSp = nNoteSp + 1
          if(iStatus == LocateDegenerate_)then
             write(*,'(a,i0,a,3es13.5)') NameSub//': shell iR=', iR, &
                  ' ill-conditioned triangle, nearest vertex used at', &
                  XyzInterp_D
          else
             write(*,'(a,i0,a,3es13.5)') NameSub//': shell iR=', iR, &
                  ' point outside the triangulated region, nearest node'// &
                  ' used at', XyzInterp_D
          end if
          if(nNoteSp == MaxNoteSp) write(*,*) NameSub// &
               ': further notes of SP_ModTriangulate are suppressed'
       end if
    case default
       write(String,'(a,i0,a,3es13.5)') NameSub//': shell iR=', iR, &
            ' zero or not finite point', XyzInterp_D
       call CON_stop(trim(String))
    end select

    if(present(DistrInterp_II))then
       DistrInterp_II = 0.0
       do i = 1, 3
          DistrInterp_II = DistrInterp_II + &
               DistrMeshSp_III(:,:,iStencil_I(i))*Weight_I(i)
       end do
    end if
    if(present(DPerp_II))then
       DPerp_II = 0.0
       if(IsDPerpMeshSp)then
          do i = 1, 3
             DPerp_II = DPerp_II + &
                  DPerpMeshSp_III(:,:,iStencil_I(i))*Weight_I(i)
          end do
       else
          do i = 1, 3
             if(iSourceSp_I(iStencil_I(i)) > 0) DPerp_II = DPerp_II + &
                  DPerp_IIBI(:,:,iSourceSp_I(iStencil_I(i)),iR)*Weight_I(i)
          end do
       end if
    end if
    if(present(iStencilOut_I))then
       if(iStatus == LocateOutside_)then
          iStencilOut_I = -1
       else
          ! line index of the node (first line of a merged node), or the
          ! pole (-1, -2) and auxiliary node (-3, -4) codes
          iStencilOut_I = iCodeMeshSp_I(iStencil_I)
       end if
    end if
    if(present(WeightOut_I))then
       if(iStatus == LocateOutside_)then
          WeightOut_I = 0.0
       else
          WeightOut_I = Weight_I
       end if
    end if
  end subroutine interpolate_trmesh_sp
  !============================================================================
  function trmesh_error_meaning(iErrorTrmesh) result(String)

    ! Meaning of the error code of trmesh (Renka, TOMS 772, STRIPACK)
    integer, intent(in) :: iErrorTrmesh
    character(len=120) :: String
    !--------------------------------------------------------------------------
    select case(iErrorTrmesh)
    case(0)
       String = 'no error'
    case(-1)
       String = 'N < 3 nodes'
    case(-2)
       String = 'the first three nodes lie on one great circle'
    case(1:)
       write(String,'(a,i0,a)') 'node ', iErrorTrmesh, &
            ' coincides with a later node (mesh covers only the nodes before)'
    case default
       String = 'unknown trmesh error code'
    end select

  end function trmesh_error_meaning
  !============================================================================
  real function det3(a_D, b_D, p_D)

    ! Same expression and argument order as det(x1,..,z0) of STRIPACK:
    ! det3(a,b,p) = p.(a x b) >= 0 iff p is left of the arc a->b
    real, intent(in) :: a_D(3), b_D(3), p_D(3)
    !--------------------------------------------------------------------------
    det3 = p_D(1)*(a_D(2)*b_D(3) - b_D(2)*a_D(3)) &
         - p_D(2)*(a_D(1)*b_D(3) - b_D(1)*a_D(3)) &
         + p_D(3)*(a_D(1)*b_D(2) - b_D(1)*a_D(2))

  end function det3
  !============================================================================
  real function height3(a_D, b_D, c_D)

    ! Smallest height of the triangle a, b, c of unit vectors: |det| divided
    ! by the longest chord. About the angular distance [rad] of the middle
    ! node from the great circle through the other two. 0 if collinear.
    real, intent(in) :: a_D(3), b_D(3), c_D(3)
    real :: Chord
    !--------------------------------------------------------------------------
    Chord = max(sqrt(sum((a_D - b_D)**2)), sqrt(sum((b_D - c_D)**2)), &
         sqrt(sum((c_D - a_D)**2)))
    if(Chord <= 0.0)then
       height3 = 0.0
    else
       height3 = abs(det3(a_D, b_D, c_D))/Chord
    end if

  end function height3
  !============================================================================
  subroutine sort_index(n, Key_I, iSort_I)

    ! Heap sort: Key_I(iSort_I) is ascending; equal keys by index (so the
    ! result does not depend on the compiler)
    integer, intent(in) :: n
    real,    intent(in) :: Key_I(n)
    integer, intent(out):: iSort_I(n)
    integer :: i, iEnd, iTmp
    !--------------------------------------------------------------------------
    do i = 1, n
       iSort_I(i) = i
    end do
    if(n < 2) RETURN
    do i = n/2, 1, -1
       call sift_down(i, n)
    end do
    do iEnd = n, 2, -1
       iTmp = iSort_I(1); iSort_I(1) = iSort_I(iEnd); iSort_I(iEnd) = iTmp
       call sift_down(1, iEnd - 1)
    end do
  contains
    !==========================================================================
    logical function is_less(i, j)
      integer, intent(in) :: i, j
      !------------------------------------------------------------------------
      is_less = Key_I(i) < Key_I(j) .or. (Key_I(i) == Key_I(j) .and. i < j)
    end function is_less
    !==========================================================================
    subroutine sift_down(iStart, iLast)
      integer, intent(in) :: iStart, iLast
      integer :: iRoot, iChild, iTmp
      !------------------------------------------------------------------------
      iRoot = iStart
      do
         iChild = 2*iRoot
         if(iChild > iLast) EXIT
         if(iChild < iLast)then
            if(is_less(iSort_I(iChild), iSort_I(iChild+1))) &
                 iChild = iChild + 1
         end if
         if(.not.is_less(iSort_I(iRoot), iSort_I(iChild))) EXIT
         iTmp = iSort_I(iRoot); iSort_I(iRoot) = iSort_I(iChild)
         iSort_I(iChild) = iTmp
         iRoot = iChild
      end do
    end subroutine sift_down
    !==========================================================================
  end subroutine sort_index
  !============================================================================
  subroutine merge_nodes(nNode, Xyz_DI, IsAdded_I, Tol, iRep_I)

    ! Single-linkage groups of the nodes closer than the angle Tol [rad].
    ! iRep_I(i) = representative of the group of node i: the data node with
    ! the smallest index, or the added node if the group has no data node.
    ! Cost O(N log N): sort along a generic direction d and sweep; a pair
    ! closer than the chord TolChord has |d.(x_i - x_j)| <= TolChord.
    ! The direction is not a lattice direction, so nodes on one latitude
    ! ring or great circle do not pile up in the sweep window.
    integer, intent(in) :: nNode
    real,    intent(in) :: Xyz_DI(3,nNode)
    logical, intent(in) :: IsAdded_I(nNode)
    real,    intent(in) :: Tol
    integer, intent(out):: iRep_I(nNode)

    real :: Dir_D(3), TolChord, TolChord2
    real, allocatable :: Key_I(:)
    integer, allocatable :: iSort_I(:)
    integer :: iA, iB, i, j
    !--------------------------------------------------------------------------
    do i = 1, nNode
       iRep_I(i) = i
    end do
    if(Tol <= 0.0 .or. nNode < 2) RETURN
    Dir_D = [sqrt(2.0/3.0), sqrt(2.0/15.0), sqrt(0.2)]
    TolChord = 2*sin(0.5*min(Tol, 1.0))
    TolChord2 = TolChord**2
    allocate(Key_I(nNode), iSort_I(nNode))
    do i = 1, nNode
       Key_I(i) = sum(Dir_D*Xyz_DI(:,i))
    end do
    call sort_index(nNode, Key_I, iSort_I)
    do iA = 1, nNode - 1
       i = iSort_I(iA)
       do iB = iA + 1, nNode
          j = iSort_I(iB)
          if(Key_I(j) - Key_I(i) > TolChord) EXIT
          if(sum((Xyz_DI(:,i) - Xyz_DI(:,j))**2) <= TolChord2) &
               call join(i, j)
       end do
    end do
    ! Point every node to its root
    do i = 1, nNode
       iRep_I(i) = find_root(i)
    end do
    deallocate(Key_I, iSort_I)
  contains
    !==========================================================================
    integer function find_root(i)
      integer, intent(in) :: i
      integer :: k
      !------------------------------------------------------------------------
      k = i
      do while(iRep_I(k) /= k)
         iRep_I(k) = iRep_I(iRep_I(k)) ! path halving
         k = iRep_I(k)
      end do
      find_root = k
    end function find_root
    !==========================================================================
    subroutine join(i, j)
      ! The root of every group is its preferred node: a data node before
      ! an added node, then the smaller index
      integer, intent(in) :: i, j
      integer :: iRoot, jRoot
      !------------------------------------------------------------------------
      iRoot = find_root(i); jRoot = find_root(j)
      if(iRoot == jRoot) RETURN
      if(is_preferred(iRoot, jRoot))then
         iRep_I(jRoot) = iRoot
      else
         iRep_I(iRoot) = jRoot
      end if
    end subroutine join
    !==========================================================================
    logical function is_preferred(i, j)
      integer, intent(in) :: i, j
      !------------------------------------------------------------------------
      if(IsAdded_I(i) .neqv. IsAdded_I(j))then
         is_preferred = .not.IsAdded_I(i)
      else
         is_preferred = i < j
      end if
    end function is_preferred
    !==========================================================================
  end subroutine merge_nodes
  !============================================================================
  subroutine order_nodes(nNode, Xyz_DI, TolHeight, iPerm_I, IsReordered, &
       IsGreatCircle, Normal_D, HeightFirst)

    ! Check the first three nodes. If their smallest height is >= TolHeight
    ! nothing changes (iPerm_I = identity). Otherwise node 1 is kept, node
    ! i2 = the node farthest from the axis +-x1 and node i3 = the node
    ! farthest from the great circle through x1 and x_i2 are moved to
    ! positions 2 and 3, the other nodes keep their order. If no node is
    ! farther than TolHeight from that great circle, all nodes lie on one
    ! great circle: IsGreatCircle = .true., Normal_D = its unit normal.
    ! New position k holds the old node iPerm_I(k).
    integer, intent(in) :: nNode
    real,    intent(in) :: Xyz_DI(3,nNode)
    real,    intent(in) :: TolHeight
    integer, intent(out):: iPerm_I(nNode)
    logical, intent(out):: IsReordered, IsGreatCircle
    real,    intent(out):: Normal_D(3), HeightFirst

    integer :: i, k, i2, i3
    real :: Cross_D(3), Cross2, Cross2Max, Dist, DistMax
    !--------------------------------------------------------------------------
    do i = 1, nNode
       iPerm_I(i) = i
    end do
    IsReordered = .false.; IsGreatCircle = .false.; Normal_D = 0.0
    HeightFirst = 0.0
    if(nNode < 3) RETURN
    HeightFirst = height3(Xyz_DI(:,1), Xyz_DI(:,2), Xyz_DI(:,3))
    if(HeightFirst >= TolHeight) RETURN

    ! Node farthest from the axis through node 1
    i2 = 0; Cross2Max = -1.0
    do i = 2, nNode
       Cross_D = cross(Xyz_DI(:,1), Xyz_DI(:,i))
       Cross2 = sum(Cross_D**2)
       if(Cross2 > Cross2Max)then
          Cross2Max = Cross2; i2 = i
       end if
    end do
    if(Cross2Max <= 0.0)then
       ! All nodes at +-x1 (cannot happen after merging if nNode > 2)
       IsGreatCircle = .true.
       Normal_D = any_normal(Xyz_DI(:,1))
       RETURN
    end if
    Normal_D = cross(Xyz_DI(:,1), Xyz_DI(:,i2))/sqrt(Cross2Max)

    ! Node farthest from the great circle through nodes 1 and i2
    i3 = 0; DistMax = -1.0
    do i = 2, nNode
       if(i == i2) CYCLE
       Dist = abs(sum(Normal_D*Xyz_DI(:,i)))
       if(Dist > DistMax)then
          DistMax = Dist; i3 = i
       end if
    end do
    if(DistMax < TolHeight)then
       IsGreatCircle = .true.
       RETURN
    end if
    Normal_D = 0.0

    ! Move i2 and i3 to positions 2 and 3, keep the order of the others
    IsReordered = .true.
    iPerm_I(1) = 1; iPerm_I(2) = i2; iPerm_I(3) = i3
    k = 3
    do i = 2, nNode
       if(i == i2 .or. i == i3) CYCLE
       k = k + 1
       iPerm_I(k) = i
    end do
  end subroutine order_nodes
  !============================================================================
  function cross(a_D, b_D) result(c_D)
    real, intent(in) :: a_D(3), b_D(3)
    real :: c_D(3)
    !--------------------------------------------------------------------------
    c_D = [a_D(2)*b_D(3) - a_D(3)*b_D(2), a_D(3)*b_D(1) - a_D(1)*b_D(3), &
         a_D(1)*b_D(2) - a_D(2)*b_D(1)]
  end function cross
  !============================================================================
  function any_normal(a_D) result(n_D)
    ! A unit vector normal to a_D
    real, intent(in) :: a_D(3)
    real :: n_D(3), e_D(3)
    !--------------------------------------------------------------------------
    e_D = 0.0
    e_D(minloc(abs(a_D), DIM=1)) = 1.0
    n_D = cross(a_D, e_D)
    n_D = n_D/sqrt(sum(n_D**2))
  end function any_normal
  !============================================================================
  subroutine list_tri_sp(nNode, iList_I, iPointer_I, iEnd_I, nTriMax, &
       nTri, iTri_II, nBoundary)

    ! Triangles of a STRIPACK triangulation in O(N), counterclockwise seen
    ! from outside, smallest index first (same rows as TRLIST, which is
    ! O(N^2)). If nTriMax < the number of triangles, only nTri is counted.
    integer, intent(in) :: nNode
    integer, intent(in) :: iList_I(6*(nNode-2)), iPointer_I(6*(nNode-2))
    integer, intent(in) :: iEnd_I(nNode)
    integer, intent(in) :: nTriMax
    integer, intent(out):: nTri, iTri_II(3,nTriMax), nBoundary

    integer :: n1, n2, n3, lp, lpLast
    !--------------------------------------------------------------------------
    nTri = 0; nBoundary = 0
    do n1 = 1, nNode
       lpLast = iEnd_I(n1)
       if(lpLast < 1) CYCLE
       if(iList_I(lpLast) < 0) nBoundary = nBoundary + 1
       lp = lpLast
       do
          lp = iPointer_I(lp)
          n2 = iList_I(lp)
          n3 = abs(iList_I(iPointer_I(lp)))
          ! n2 < 0: last neighbor of a boundary node, no triangle follows
          if(n2 > 0 .and. n1 < n2 .and. n1 < n3)then
             nTri = nTri + 1
             if(nTri <= nTriMax) iTri_II(:,nTri) = [n1, n2, n3]
          end if
          if(lp == lpLast) EXIT
       end do
    end do
  end subroutine list_tri_sp
  !============================================================================
  subroutine build_sp(nNode, nNodeMax, Xyz_DI, iCode_I, Tol, nMesh, &
       iSource_I, iGroup_I, iList_I, iPointer_I, iEnd_I, Info, iError, &
       StringError)

    ! Input: nNode nodes Xyz_DI(:,1:nNode) of any length > 0 and their codes
    ! iCode_I(1:nNode): >= 0 for a data node (e.g. the line index), < 0 for
    ! an added node (the poles of MFLAMPA: -1 and -2). Tol = merge angle.
    ! nNodeMax >= nNode is the capacity; auxiliary nodes need nNode+2.
    ! Output:
    !   Xyz_DI(:,1:nMesh), iCode_I(1:nMesh): mesh nodes (unit vectors) and
    !     codes; auxiliary nodes have the codes iAux1Sp_, iAux2Sp_
    !   iSource_I(1:nMesh): input index of the node kept as mesh node m
    !     (0 for an auxiliary node)
    !   iGroup_I(1:nNode): mesh node that input node i belongs to. A replaced
    !     added node points to the data node that replaced it (its state
    !     must not be averaged).
    !   iList_I, iPointer_I, iEnd_I: STRIPACK structure of the mesh
    !   Info: statistics; iError = BuildOk_ or an error code, StringError
    use ModTriangulateSpherical, ONLY: trmesh
    integer, intent(in)    :: nNode, nNodeMax
    real,    intent(inout) :: Xyz_DI(3,nNodeMax)
    integer, intent(inout) :: iCode_I(nNodeMax)
    real,    intent(in)    :: Tol
    integer, intent(out)   :: nMesh
    integer, intent(out)   :: iSource_I(nNodeMax), iGroup_I(nNode)
    integer, intent(out)   :: iList_I(6*(nNodeMax-2))
    integer, intent(out)   :: iPointer_I(6*(nNodeMax-2)), iEnd_I(nNodeMax)
    type(BuildInfoType), intent(out) :: Info
    integer, intent(out)   :: iError
    character(len=*), intent(out) :: StringError

    integer :: i, m, iTri, iDegenerate
    integer, allocatable :: iRep_I(:), iMeshOf_I(:), iTri_II(:,:)
    real :: r, Normal_D(3), TolHeight, Det, Side, Cross_DI(3,3), Cross2_I(3)
    logical :: IsReordered, IsGreatCircle
    !--------------------------------------------------------------------------
    iError = BuildOk_; StringError = ''; nMesh = 0
    Info%nNodeIn = nNode
    iSource_I = 0; iGroup_I = 0
    if(nNode < 1 .or. nNode > nNodeMax)then
       iError = BuildCapacity_
       write(StringError,'(a,i0,a,i0)') 'bad node count nNode=', nNode, &
            ' nNodeMax=', nNodeMax
       RETURN
    end if

    ! 1. Unit length (STRIPACK assumes unit vectors)
    do i = 1, nNode
       r = sqrt(sum(Xyz_DI(:,i)**2))
       if(.not.(r > 0.0 .and. r <= huge(r)))then
          iError = BuildBadNode_
          write(StringError,'(a,i0,a)') 'node ', i, &
               ' has zero, NaN or infinite coordinates'
          RETURN
       end if
       Info%MaxNormDev = max(Info%MaxNormDev, abs(r - 1))
       if(abs(r - 1) > TolNormSp)then
          Xyz_DI(:,i) = Xyz_DI(:,i)/r
          Info%nNormalized = Info%nNormalized + 1
       end if
    end do

    ! 2. Merge nodes closer than Tol; compact in input order
    allocate(iRep_I(nNode), iMeshOf_I(nNode))
    call merge_nodes(nNode, Xyz_DI(:,1:nNode), iCode_I(1:nNode) < 0, Tol, &
         iRep_I)
    do i = 1, nNode
       if(iRep_I(i) == i)then
          nMesh = nMesh + 1
          iMeshOf_I(i) = nMesh
       else
          iMeshOf_I(i) = 0
          if(iCode_I(i) < 0)then
             Info%nPoleReplaced = Info%nPoleReplaced + 1
          else
             Info%nMerged = Info%nMerged + 1
          end if
       end if
    end do
    do i = 1, nNode
       iGroup_I(i) = iMeshOf_I(iRep_I(i))
    end do
    m = 0
    do i = 1, nNode
       if(iRep_I(i) /= i) CYCLE
       m = m + 1
       Xyz_DI(:,m) = Xyz_DI(:,i); iCode_I(m) = iCode_I(i); iSource_I(m) = i
    end do
    deallocate(iRep_I, iMeshOf_I)

    ! 3. First three nodes; all nodes on one great circle?
    TolHeight = max(Tol, TolHeightMinSp)
    IsGreatCircle = .false.
    if(nMesh < 2)then
       iError = BuildFewNodes_
       write(StringError,'(a,i0,a)') 'only ', nMesh, &
            ' distinct node after merging, cannot triangulate'
       RETURN
    elseif(nMesh == 2)then
       Normal_D = cross(Xyz_DI(:,1), Xyz_DI(:,2))
       r = sqrt(sum(Normal_D**2))
       if(r < TolHeight)then
          iError = BuildFewNodes_
          StringError = 'only 2 antipodal nodes after merging, cannot '// &
               'triangulate'
          RETURN
       end if
       IsGreatCircle = .true.; Normal_D = Normal_D/r
    else
       call order_and_apply
       if(iError /= BuildOk_) RETURN
    end if
    if(IsGreatCircle)then
       ! Close the sphere with two nodes at +-normal (like the poles)
       if(nMesh + 2 > nNodeMax)then
          iError = BuildCapacity_
          StringError = 'all nodes lie on one great circle and there is '// &
               'no room for the 2 auxiliary nodes'
          RETURN
       end if
       Info%IsGreatCircle = .true.
       call add_aux(Normal_D)
       call add_aux(-Normal_D)
       call order_and_apply
       if(iError /= BuildOk_) RETURN
       if(IsGreatCircle)then
          iError = BuildMesh_
          StringError = 'nodes still on one great circle after adding '// &
               'the auxiliary nodes'
          RETURN
       end if
    end if

    ! 4. trmesh and mesh check; a flat side is closed by an auxiliary node
    do
       call run_trmesh
       if(iError /= BuildOk_) RETURN
       if(Info%nDegenerate == 0) EXIT
       if(Info%nAux >= 2 .or. nMesh + 1 > nNodeMax)then
          iError = BuildMesh_
          write(StringError,'(i0,a,es10.3,a,i0,a)') Info%nDegenerate, &
               ' degenerate triangles (min det ', Info%MinDet, ') with ', &
               Info%nAux, ' auxiliary nodes'
          RETURN
       end if
       ! Normal of the plane through the origin and the flat triangle
       do i = 1, 3
          Cross_DI(:,i) = cross(Xyz_DI(:,iTri_II(i,iDegenerate)), &
               Xyz_DI(:,iTri_II(mod(i,3)+1,iDegenerate)))
          Cross2_I(i) = sum(Cross_DI(:,i)**2)
       end do
       i = maxloc(Cross2_I, DIM=1)
       Normal_D = Cross_DI(:,i)/sqrt(Cross2_I(i))
       ! Put the auxiliary node on the empty side of that plane
       Side = 0.0
       do m = 1, nMesh
          Side = Side + sum(Normal_D*Xyz_DI(:,m))
       end do
       if(Side > 0.0) Normal_D = -Normal_D
       call add_aux(Normal_D)
    end do
    Info%nMesh = nMesh
    if(allocated(iTri_II)) deallocate(iTri_II)

  contains
    !==========================================================================
    subroutine add_aux(Aux_D)
      real, intent(in) :: Aux_D(3)
      !------------------------------------------------------------------------
      nMesh = nMesh + 1
      Info%nAux = Info%nAux + 1
      Xyz_DI(:,nMesh) = Aux_D
      if(Info%nAux == 1)then
         iCode_I(nMesh) = iAux1Sp_
      else
         iCode_I(nMesh) = iAux2Sp_
      end if
      iSource_I(nMesh) = 0
    end subroutine add_aux
    !==========================================================================
    subroutine order_and_apply

      integer :: k
      integer, allocatable :: iPerm_I(:), iInv_I(:), iTmp_I(:)
      real, allocatable :: XyzTmp_DI(:,:)
      real :: HeightFirst
      !------------------------------------------------------------------------
      allocate(iPerm_I(nMesh), iInv_I(nMesh), iTmp_I(nMesh), &
           XyzTmp_DI(3,nMesh))
      call order_nodes(nMesh, Xyz_DI(:,1:nMesh), TolHeight, iPerm_I, &
           IsReordered, IsGreatCircle, Normal_D, HeightFirst)
      if(Info%nAux == 0) Info%MinHeight3 = HeightFirst
      if(.not.IsReordered) RETURN
      Info%IsReordered = .true.
      XyzTmp_DI = Xyz_DI(:,1:nMesh)
      do k = 1, nMesh
         Xyz_DI(:,k) = XyzTmp_DI(:,iPerm_I(k))
         iInv_I(iPerm_I(k)) = k
      end do
      iTmp_I = iCode_I(1:nMesh)
      iCode_I(1:nMesh) = iTmp_I(iPerm_I)
      iTmp_I = iSource_I(1:nMesh)
      iSource_I(1:nMesh) = iTmp_I(iPerm_I)
      do k = 1, nNode
         if(iGroup_I(k) > 0) iGroup_I(k) = iInv_I(iGroup_I(k))
      end do
      ! The new first triple must be well conditioned
      if(height3(Xyz_DI(:,1), Xyz_DI(:,2), Xyz_DI(:,3)) < TolHeight)then
         iError = BuildMesh_
         StringError = 'no well-conditioned first triangle found'
      end if
    end subroutine order_and_apply
    !==========================================================================
    subroutine run_trmesh

      real, allocatable :: x_I(:), y_I(:), z_I(:)
      integer :: iErrorTrmesh, nTri, nBoundary, nTriExpected
      !------------------------------------------------------------------------
      allocate(x_I(nMesh), y_I(nMesh), z_I(nMesh))
      x_I = Xyz_DI(1,1:nMesh); y_I = Xyz_DI(2,1:nMesh)
      z_I = Xyz_DI(3,1:nMesh)
      iList_I = 0; iPointer_I = 0; iEnd_I = 0
      call trmesh(nMesh, x_I, y_I, z_I, iList_I(1:6*(nMesh-2)), &
           iPointer_I(1:6*(nMesh-2)), iEnd_I(1:nMesh), iErrorTrmesh)
      deallocate(x_I, y_I, z_I)
      Info%nTrmesh = Info%nTrmesh + 1
      Info%iErrorTrmesh = iErrorTrmesh
      if(iErrorTrmesh /= 0)then
         iError = BuildTrmesh_
         write(StringError,'(a,i0,a,i0,3a)') 'trmesh failed for ', nMesh, &
              ' nodes with iErrorTrmesh=', iErrorTrmesh, ' (', &
              trim(trmesh_error_meaning(iErrorTrmesh)), ')'
         RETURN
      end if
      ! Check the mesh
      if(allocated(iTri_II)) deallocate(iTri_II)
      allocate(iTri_II(3,2*nMesh))
      call list_tri_sp(nMesh, iList_I(1:6*(nMesh-2)), &
           iPointer_I(1:6*(nMesh-2)), iEnd_I(1:nMesh), 2*nMesh, nTri, &
           iTri_II, nBoundary)
      Info%nTri = nTri; Info%nBoundary = nBoundary
      ! Euler: closed mesh 2N-4 triangles, mesh with NB boundary nodes 2N-NB-2
      if(nBoundary == 0)then
         nTriExpected = 2*nMesh - 4
      else
         nTriExpected = 2*nMesh - nBoundary - 2
      end if
      if(nTri /= nTriExpected .or. nTri > 2*nMesh)then
         iError = BuildMesh_
         write(StringError,'(a,i0,a,i0,a,i0,a,i0,a)') 'mesh has ', nTri, &
              ' triangles for ', nMesh, ' nodes and ', nBoundary, &
              ' boundary nodes (expected ', nTriExpected, ')'
         RETURN
      end if
      Info%nDegenerate = 0; Info%MinDet = huge(1.0); iDegenerate = 0
      do iTri = 1, nTri
         Det = det3(Xyz_DI(:,iTri_II(1,iTri)), Xyz_DI(:,iTri_II(2,iTri)), &
              Xyz_DI(:,iTri_II(3,iTri)))
         Info%MinDet = min(Info%MinDet, Det)
         if(.not.(Det > DetMinSp))then
            Info%nDegenerate = Info%nDegenerate + 1
            if(iDegenerate == 0) iDegenerate = iTri
         end if
      end do
    end subroutine run_trmesh
    !==========================================================================
  end subroutine build_sp
  !============================================================================
  real function area_signed(a_D, b_D, c_D)

    ! Signed area of the spherical triangle a, b, c (unit vectors), > 0 if
    ! counterclockwise: Van Oosterom and Strackee (1983),
    ! tan(E/2) = det(a,b,c)/(1 + a.b + b.c + c.a). Accurate also for thin
    ! triangles (unlike the angle excess with acos used by AREAS).
    real, intent(in) :: a_D(3), b_D(3), c_D(3)
    !--------------------------------------------------------------------------
    area_signed = 2*atan2(det3(a_D, b_D, c_D), 1 + sum(a_D*b_D) &
         + sum(b_D*c_D) + sum(c_D*a_D))

  end function area_signed
  !============================================================================
  subroutine locate_sp(Xyz_D, nNode, Xyz_DI, iList_I, iPointer_I, &
       iEnd_I, UsePlanar, iStencil_I, Weight_I, iStatus)

    ! Triangle containing the direction of Xyz_D and interpolation weights.
    ! The triangle is found by STRIPACK (find_triangle_orig, TRFIND). Its
    ! weights are not used: they are NaN with IsTriangleFound=.true. for
    ! nearly coincident nodes. The weights are recomputed here from the
    ! vertex coordinates: planar (UsePlanar: barycentric coordinates of the
    ! central projection onto the triangle plane, the same triple products
    ! as TRFIND) or spherical (signed areas of the 3 sub-triangles, stable
    ! formula). They are clipped at 0 and divided by their sum, so they are
    ! finite, >= 0 and sum to 1.
    ! iStatus = LocateInside_     normal result
    !         = LocateDegenerate_ the triangle is found but its weights are
    !           ill-conditioned (sum < half the triangle det or area, or not
    !           finite): weight 1 at the nearest vertex of the triangle
    !         = LocateOutside_    no triangle (point outside a mesh with a
    !           boundary): iStencil_I = nearest node (3 times), weights
    !           (1,0,0)
    !         = LocateBadPoint_   zero, NaN or infinite Xyz_D: stencil 0
    use ModTriangulateSpherical, ONLY: find_triangle_orig
    real,    intent(in) :: Xyz_D(3)
    integer, intent(in) :: nNode
    real,    intent(in) :: Xyz_DI(3,nNode)
    integer, intent(in) :: iList_I(6*(nNode-2)), iPointer_I(6*(nNode-2))
    integer, intent(in) :: iEnd_I(nNode)
    logical, intent(in) :: UsePlanar
    integer, intent(out):: iStencil_I(3)
    real,    intent(out):: Weight_I(3)
    integer, intent(out):: iStatus

    real    :: r, p_D(3), b_I(3), SumB, DetTri, WeightOrig_I(3), Dot, DotMax
    real    :: v1_D(3), v2_D(3), v3_D(3)
    logical :: IsTriangleFound
    integer :: i, k
    !--------------------------------------------------------------------------
    iStencil_I = 0; Weight_I = 0.0
    r = sqrt(sum(Xyz_D**2))
    if(.not.(r > 0.0 .and. r <= huge(r)))then
       iStatus = LocateBadPoint_
       RETURN
    end if
    p_D = Xyz_D/r
    call find_triangle_orig(p_D, nNode, Xyz_DI, iList_I, iPointer_I, &
         iEnd_I, WeightOrig_I, IsTriangleFound, iStencil_I)
    if(IsTriangleFound)then
       v1_D = Xyz_DI(:,iStencil_I(1))
       v2_D = Xyz_DI(:,iStencil_I(2))
       v3_D = Xyz_DI(:,iStencil_I(3))
       if(UsePlanar)then
          ! As TRFIND: b1 = det(V2,V3,P), b2 = det(V3,V1,P), b3 = det(V1,V2,P)
          b_I(1) = det3(v2_D, v3_D, p_D)
          b_I(2) = det3(v3_D, v1_D, p_D)
          b_I(3) = det3(v1_D, v2_D, p_D)
          DetTri = det3(v1_D, v2_D, v3_D)
       else
          b_I(1) = area_signed(p_D, v2_D, v3_D)
          b_I(2) = area_signed(v1_D, p_D, v3_D)
          b_I(3) = area_signed(v1_D, v2_D, p_D)
          DetTri = area_signed(v1_D, v2_D, v3_D)
       end if
       b_I = max(b_I, 0.0)
       SumB = b_I(1) + b_I(2) + b_I(3)
       ! For P in the triangle SumB >= det(V1,V2,V3) (planar) and
       ! SumB = area (spherical); 0.5 leaves room for rounding
       if(DetTri > 0.0 .and. SumB > 0.5*DetTri .and. SumB <= huge(SumB))then
          Weight_I = b_I/SumB
          iStatus = LocateInside_
          RETURN
       end if
       ! Ill-conditioned weights: nearest vertex of the triangle
       k = 1; DotMax = sum(p_D*v1_D)
       if(sum(p_D*v2_D) > DotMax)then
          k = 2; DotMax = sum(p_D*v2_D)
       end if
       if(sum(p_D*v3_D) > DotMax) k = 3
       Weight_I = 0.0; Weight_I(k) = 1.0
       iStatus = LocateDegenerate_
       RETURN
    end if
    ! Outside the triangulated region: nearest node
    k = 1; DotMax = -huge(1.0)
    do i = 1, nNode
       Dot = sum(p_D*Xyz_DI(:,i))
       if(Dot > DotMax)then
          DotMax = Dot; k = i
       end if
    end do
    iStencil_I = k
    Weight_I = [1.0, 0.0, 0.0]
    iStatus = LocateOutside_

  end subroutine locate_sp
  !============================================================================
  subroutine fix_state_sp(iNodeToFix, nNode, iList_I, iPointer_I, &
       iEnd_I, Xyz_DI, IsSkip_I, nVar, State_VI)

    ! State of an added node (pole or auxiliary node) from its neighbors:
    ! the same inverse-distance average as fix_state of
    ! ModTriangulateSpherical (same arithmetic and neighbor order, the
    ! state of a neighbor closer than cTiny is copied), except that the
    ! neighbors with IsSkip_I = .true. (the other added nodes, whose state
    ! is not known) are left out, unless all neighbors are such nodes.
    use ModNumConst, ONLY: cTiny
    integer, intent(in) :: iNodeToFix, nNode
    integer, intent(in) :: iList_I(6*(nNode-2)), iPointer_I(6*(nNode-2))
    integer, intent(in) :: iEnd_I(nNode)
    real,    intent(in) :: Xyz_DI(3,nNode)
    logical, intent(in) :: IsSkip_I(nNode)
    integer, intent(in) :: nVar
    real, intent(inout) :: State_VI(nVar,nNode)

    integer :: iPass, nUsed, iNeighbor, iNeighborIndex, iNeighborIndexLast
    real :: Xyz_D(3), Dist, DistInv, SumDistInv, State_V(nVar)
    !--------------------------------------------------------------------------
    iNeighborIndexLast = iEnd_I(iNodeToFix)
    if(iNeighborIndexLast < 1) RETURN
    Xyz_D = Xyz_DI(:,iNodeToFix)
    do iPass = 1, 2
       State_V = 0.0; SumDistInv = 0.0; nUsed = 0
       iNeighborIndex = iNeighborIndexLast
       do
          iNeighborIndex = iPointer_I(iNeighborIndex)
          iNeighbor = abs(iList_I(iNeighborIndex))
          if(iPass == 2 .or. .not.IsSkip_I(iNeighbor))then
             Dist = sqrt(sum((Xyz_DI(:,iNeighbor) - Xyz_D)**2))
             if(Dist < cTiny)then
                State_VI(:,iNodeToFix) = State_VI(:,iNeighbor)
                RETURN
             end if
             DistInv = 1/Dist
             State_V = State_V + DistInv*State_VI(:,iNeighbor)
             SumDistInv = SumDistInv + DistInv
             nUsed = nUsed + 1
          end if
          if(iNeighborIndex == iNeighborIndexLast) EXIT
       end do
       if(nUsed > 0)then
          State_VI(:,iNodeToFix) = State_V/SumDistInv
          RETURN
       end if
    end do

  end subroutine fix_state_sp
  !============================================================================
end module SP_ModTriangulate
!==============================================================================
