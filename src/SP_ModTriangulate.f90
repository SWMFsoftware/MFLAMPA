!  Copyright (C) 2002 Regents of the University of Michigan,
!  portions used with permission
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module SP_ModTriangulate

  use ModUtilities, ONLY: CON_stop
  use SP_ModGrid,   ONLY: nLine, nLineAll, Used_B, iLineAll0, &
       search_line, MHData_VIB, X_, Y_, Z_, nP, nMu
  use SP_ModProc,   ONLY: iProc, nProc, iComm, iError
  use SP_ModSize,   ONLY: nDim

  implicit none

  SAVE

  PRIVATE ! Except

  public:: read_param         ! Read module parameters
  public:: intersect_surf     ! Intersect on a sphere with a given radius
  public:: build_trmesh       ! Build the triangulated mesh skeleton
  public:: interpolate_trmesh ! Interpolate at the specified location

  ! Test triangulation
  logical, public :: DoTestTri = .false.
  real,    public :: XyzLocTestTri_I(nDim)

  ! If we use poles in triangulation
  logical, public :: UsePoleTri   = .false.
  logical, public :: UsePlanarTri = .true.
  real,    public, parameter :: iSouthPoleTri_ = -1, iNorthPoleTri_ = -2
contains
  !============================================================================
  subroutine read_param(NameCommand)

    use ModReadParam, ONLY: read_var
    character(len=*), intent(in) :: NameCommand ! From PARAM.in
    character(len=*), parameter:: NameSub = 'read_param'
    !--------------------------------------------------------------------------
    select case(NameCommand)
    case("#TESTTRIANGULATE")
       call read_var('TestTriangulate', DoTestTri)
       if(DoTestTri) then
          call read_var('XLocTestTri', XyzLocTestTri_I(1))
          call read_var('YLocTestTri', XyzLocTestTri_I(2))
          call read_var('ZLocTestTri', XyzLocTestTri_I(3))
       end if
    case('#TRIANGULATION')
       ! get pole triangulartion flag
       call read_var('UsePoleTriangulation', UsePoleTri)
       ! get the triangulation approach flag
       call read_var('UsePlanarTriangles', UsePlanarTri)
    case default
       call CON_stop(NameSub//': Unknown command '//NameCommand)
    end select

  end subroutine read_param
  !============================================================================
  subroutine intersect_surf(rSurf, XyzReachRUnit_DI, &
       iLineReach_I, Log10DistrReachR_IIB, nReachR)

    use ModMpi
    use SP_ModDistribution, ONLY: Distribution_CB
    ! Input: radius of the spherical surface
    real,    intent(in) :: rSurf
    ! Output: useful intersection points on a unit sphere
    ! Here the last index is 2:nLineAll+1 for normal nLineAll lines
    real,    intent(out):: XyzReachRUnit_DI(X_:Z_, 1:nLineAll+2)
    integer, intent(out):: iLineReach_I(1:nLineAll+2)
    real,    intent(out):: Log10DistrReachR_IIB(0:nP+1, 1:nMu, 1:nLineAll+2)
    ! Count how many field lines reach the sphere at the given r=rSurf
    integer, intent(out):: nReachR

    ! loop variables
    integer :: iLine, iReachR
    ! indexes of corresponding node, latitude and longitude
    integer :: iLineAll
    ! index of particle just above the radius
    integer :: iAboveR
    ! interpolation weight (in radial direction)
    real    :: WeightR
    ! skip a field line not reaching radius of output sphere
    logical :: DoReachR_I(nLineAll)
    ! xyz coordinates of all intersection point or average direction
    real    :: XyzUnit_DI(X_:Z_, nLineAll)
    real    :: Log10DistrR_IIB(0:nP+1, 1:nMu, nLineAll)

    character(len=*), parameter:: NameSub = 'triangulate_surf'
    !------------------------------------------------------------------------

    XyzUnit_DI = 0.0; Log10DistrR_IIB = 0.0
    DoReachR_I = .false.; iLineReach_I = 0
    XyzReachRUnit_DI = 0.0; Log10DistrReachR_IIB = 0.0

    ! go over all lines on the processor and find the point of
    ! intersection with output sphere if present
    LINE:do iLine = 1, nLine
       if(.not.Used_B(iLine)) CYCLE LINE
       iLineAll = iLineAll0 + iLine

       ! Find the particle just above the given radius
       call search_line(iLine, rSurf, &
            iAboveR, DoReachR_I(iLineAll), WeightR)
       DoReachR_I(iLineAll) = DoReachR_I(iLineAll) .and. iAboveR /= 1
       ! if no intersection found -> proceed to the next line
       if(.not.DoReachR_I(iLineAll)) CYCLE LINE

       ! Found intersection => get POS & log10(f[IO_Unit]) at that location
       ! Find coordinates and log(Distribution) at intersection
       XyzUnit_DI(:, iLineAll) = ( &
            MHData_VIB(X_:Z_, iAboveR-1, iLine)*(1-WeightR) +  &
            MHData_VIB(X_:Z_, iAboveR,   iLine)*   WeightR )/rSurf
       Log10DistrR_IIB(:, :, iLineAll) = ( &
            log10(Distribution_CB(:, :, iAboveR,   iLine))* WeightR + &
            log10(Distribution_CB(:, :, iAboveR-1, iLine))*(1-WeightR))
    end do LINE !  iLine

    ! Gather interpolated coordinates on the source processor and Broadcast
    call gather_proc

  contains
    !==========================================================================
    subroutine gather_proc

      !------------------------------------------------------------------------

      ! Gather interpolated coordinates on the source processor
      if(nProc > 1) then
         if(iProc == 0) then
            call MPI_REDUCE(MPI_IN_PLACE, XyzUnit_DI, &
                 3*nLineAll, MPI_REAL, MPI_SUM,       &
                 0, iComm, iError)
            call MPI_REDUCE(MPI_IN_PLACE,             &
                 Log10DistrR_IIB, (nP+2)*nMu*nLineAll,&
                 MPI_REAL, MPI_SUM, 0, iComm, iError)
            call MPI_REDUCE(MPI_IN_PLACE, DoReachR_I, &
                 nLineAll, MPI_LOGICAL, MPI_LOR,      &
                 0, iComm, iError)
         else
            call MPI_REDUCE(XyzUnit_DI, XyzUnit_DI,   &
                 3*nLineAll, MPI_REAL, MPI_SUM,       &
                 0, iComm, iError)
            call MPI_REDUCE(Log10DistrR_IIB,          &
                 Log10DistrR_IIB, (nP+2)*nMu*nLineAll,&
                 MPI_REAL, MPI_SUM, 0, iComm, iError)
            call MPI_REDUCE(DoReachR_I, DoReachR_I,   &
                 nLineAll, MPI_LOGICAL, MPI_LOR,      &
                 0, iComm, iError)
         end if
      end if

      ! Send useful interpolated coordinates to the source processor
      if(iProc == 0) then
         ! Calculate the useful points on the sphere
         nReachR = count(DoReachR_I)
         if(nReachR == nLineAll) then
            ! For most cases, the satellite used falls in 1.1-240 AU
            XyzReachRUnit_DI(:, 2:nLineAll+1) = XyzUnit_DI
            do iLineAll = 1, nLineAll
               iLineReach_I(iLineAll+1) = iLineAll
            end do
            Log10DistrReachR_IIB(:, :, 2:nLineAll+1) = Log10DistrR_IIB
         else
            ! Otherwise, we will spend some time reorganizing the points
            iReachR = 1
            do iLineAll = 1, nLineAll
               if(DoReachR_I(iLineAll)) then
                  iReachR = iReachR + 1
                  XyzReachRUnit_DI(:, iReachR) = XyzUnit_DI(:, iLineAll)
                  iLineReach_I(iReachR) = iLineAll
                  Log10DistrReachR_IIB(:, :, iReachR) = &
                       Log10DistrR_IIB(:, :, iLineAll)
               end if
            end do
         end if
      end if

      ! Broadcast the coordinates and flags to all processors
      call MPI_BCAST(XyzReachRUnit_DI, 3*(nLineAll+2), &
           MPI_REAL, 0, iComm, iError)
      call MPI_BCAST(iLineReach_I, nLineAll+2, &
           MPI_INTEGER, 0, iComm, iError)
      call MPI_BCAST(Log10DistrReachR_IIB, (nP+2)*nMu*(nLineAll+2), &
           MPI_REAL, 0, iComm, iError)
      call MPI_BCAST(DoReachR_I, nLineAll, MPI_LOGICAL, 0, iComm, iError)
      call MPI_BCAST(nReachR, 1, MPI_INTEGER, 0, iComm, iError)

    end subroutine gather_proc
    !==========================================================================
  end subroutine intersect_surf
  !============================================================================
  subroutine build_trmesh(nReachR, XyzReachRUnit_DI, iLineReach_I, &
       nTriMesh, lidTri, ridTri, iList_I, iPointer_I, iEnd_I)

    use ModTriangulateSpherical, ONLY: trmesh
    ! Input: how many points reach the sphere
    integer, intent(in)  :: nReachR
    ! In/Outputs: Can be changed if UsePoleTri
    real,    intent(inout):: XyzReachRUnit_DI(X_:Z_, 1:nLineAll+2)
    integer, intent(inout):: iLineReach_I(1:nLineAll+2)
    ! Outputs: Count, left and right indices of the triangulated mesh
    integer, intent(out)  :: nTriMesh, lidTri, ridTri
    ! Outputs: Arrays to construct a triangular mesh on a sphere
    integer, intent(out), allocatable :: iList_I(:), iPointer_I(:), iEnd_I(:)
    !--------------------------------------------------------------------------

    ! For poles
    if(UsePoleTri) then
       ! Add two grid nodes at the poles:
       lidTri = 1
       ridTri = nReachR + 2
       XyzReachRUnit_DI(:, lidTri) = [0.0, 0.0, -1.0]
       XyzReachRUnit_DI(:, ridTri) = [0.0, 0.0, +1.0]
       iLineReach_I(lidTri) = iSouthPoleTri_
       iLineReach_I(ridTri) = iNorthPoleTri_
    else
       lidTri = 2
       ridTri = nReachR + 1
    end if

    nTriMesh = ridTri - lidTri + 1
    ! Allocate and initialize arrays for triangulation and
    ! interpolation; if allocated, first deallocate them
    if(allocated(iList_I))    deallocate(iList_I)
    if(allocated(iPointer_I)) deallocate(iPointer_I)
    if(allocated(iEnd_I))     deallocate(iEnd_I)
    allocate(iList_I(6*(nTriMesh-2)), &
         iPointer_I(6*(nTriMesh-2)), iEnd_I(nTriMesh))
    iList_I = 0; iPointer_I = 0; iEnd_I = 0
    ! Construct the Triangular mesh used for interpolation
    call trmesh(nTriMesh,                     &
         XyzReachRUnit_DI(X_, lidTri:ridTri), &
         XyzReachRUnit_DI(Y_, lidTri:ridTri), &
         XyzReachRUnit_DI(Z_, lidTri:ridTri), &
         iList_I, iPointer_I, iEnd_I, iError)

  end subroutine build_trmesh
  !============================================================================
  subroutine interpolate_trmesh(rSurf, XyzInterp_I, nTriMesh, &
       lidTri, ridTri, iList_I, iPointer_I, iEnd_I, &
       XyzReachRUnit_DI, Log10DistrReachR_IIB, Log10DistrInterp_II, &
       IsTriangleFound, iStencil_I, Weight_I)

    use ModTriangulateSpherical, ONLY: fix_state, &
         find_triangle_orig, find_triangle_sph
    ! Input: Radius of the spherical surface
    real,    intent(in)   :: rSurf
    ! Input: Xyz coordinates of the point for interpolations
    real,    intent(in)   :: XyzInterp_I(nDim)
    ! Inputs: Count, left and right indices of the triangulated mesh
    integer, intent(in)   :: nTriMesh, lidTri, ridTri
    ! Inputs: Arrays to construct a triangular mesh on a sphere
    integer, intent(in)   :: iList_I(6*(nTriMesh-2)), &
         iPointer_I(6*(nTriMesh-2)), iEnd_I(nTriMesh)
    ! Input: Xyz coordinates of the intersection points on the sphere
    real,    intent(in)   :: XyzReachRUnit_DI(X_:Z_, 1:nLineAll+2)
    ! Input: Values at the intersection points on the sphere
    real,    intent(inout):: Log10DistrReachR_IIB(0:nP+1, 1:nMu, 1:nLineAll+2)
    ! Output: Interpolated values at XyzInterp_I
    real,    intent(inout):: Log10DistrInterp_II(0:nP+1, 1:nMu)
    ! variables for interpolation in a triangular mesh
    logical, intent(out) :: IsTriangleFound
    integer, intent(out) :: iStencil_I(3)
    real,    intent(out) :: Weight_I(3)

    ! loop variables
    integer :: iSat, iMu, i
    character(len=*), parameter:: NameSub = 'interpolate_trmesh'
    !--------------------------------------------------------------------------
    IsTriangleFound = .false.
    iStencil_I = -1
    Weight_I = 0.0

    ! Find the triangle where the satellite locates
    if(UsePlanarTri) then
       call find_triangle_orig(                 &
            XyzInterp_I/rSurf, nTriMesh,  &
            XyzReachRUnit_DI(:, lidTri:ridTri), &
            iList_I, iPointer_I, iEnd_I,        &
            Weight_I, IsTriangleFound, iStencil_I)
    else
       call find_triangle_sph(                  &
            XyzInterp_I/rSurf, nTriMesh,  &
            XyzReachRUnit_DI(:, lidTri:ridTri), &
            iList_I, iPointer_I, iEnd_I,        &
            Weight_I(1), Weight_I(2),           &
            Weight_I(3), IsTriangleFound,       &
            iStencil_I(1), iStencil_I(2), iStencil_I(3))
    end if
    if(.not.IsTriangleFound) then
       write(*,*) NameSub, ' WARNING: Interpolation fails on the', &
            'triangulated sphere, at the location of x,y,z=', XyzInterp_I
       ! EXIT TRI_INTERPOLATE
       RETURN
    end if

    ! Fix states (log10 distribution) at the polar nodes
    if(UsePoleTri) then
       do iMu = 1, nMu
          ! North:
          call fix_state( &
               iNodeToFix = ridTri,     &
               nNode      = nTriMesh,   &
               iList_I    = iList_I,    &
               iPointer_I = iPointer_I, &
               iEnd_I     = iEnd_I,     &
               Xyz_DI     = XyzReachRUnit_DI(:, lidTri:ridTri), &
               nVar       = nP+2,       &
               State_VI   = Log10DistrReachR_IIB(:, iMu, lidTri:ridTri))
          ! South:
          call fix_state( &
               iNodeToFix = lidTri,     &
               nNode      = nTriMesh,   &
               iList_I    = iList_I,    &
               iPointer_I = iPointer_I, &
               iEnd_I     = iEnd_I,     &
               Xyz_DI     = XyzReachRUnit_DI(:, lidTri:ridTri), &
               nVar       = nP+2,       &
               State_VI   = Log10DistrReachR_IIB(:, iMu, lidTri:ridTri))
       end do
    end if
    ! Interpolate the log10(distribution) at satellite as outputs
    do i = 1, 3
       Log10DistrInterp_II = Log10DistrInterp_II + &
            Log10DistrReachR_IIB(:, :, iStencil_I(i))*Weight_I(i)
    end do

  end subroutine interpolate_trmesh
  !============================================================================
end module SP_ModTriangulate
!==============================================================================
