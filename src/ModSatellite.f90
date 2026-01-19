!  Copyright (C) 2002 Regents of the University of Michigan,
!  portions used with permission
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module SP_ModSatellite

  use ModUtilities, ONLY: open_file, close_file, CON_stop
  use SP_ModProc, ONLY: iProc, iComm, iError
  use SP_ModSize, ONLY: nDim
  use SP_ModTestFunc, ONLY: lVerbose, test_start, test_stop
  use SP_ModGrid, ONLY: TypeCoordSystem

  implicit none

  save

  private ! Except

  public:: read_param                 ! read satellite-related parameters
  public:: init                       ! Initialize for satellite triangulation
  public:: read_satellite_input_files ! read satellite trajectories
  public:: set_satellite_positions    ! set satellite positions

  ! ----- For saving the distribution function along SATELLITES -----
  integer, public    :: nSat   = 0    ! number of satellites for DistrTraj_
  integer, parameter :: MaxSat = 300  ! Max number of satellites

  ! Names and unit numbers for satellite files
  character(len=50), public:: NameFileSat_I(MaxSat)
  ! Names of the satellite
  character(len=50), public:: NameSat_I(MaxSat)

  ! Current positions
  real, public:: XyzSat_DI(nDim, MaxSat)

  ! variables to record tracked and current satellite position indices
  logical, public:: UseSatellite               = .false.
  logical, public:: DoTrackSatellite_I(MaxSat) = .false.
  integer, public:: iPointCurrentSat_I(MaxSat) = 1

  ! Local variables
  integer, public           :: nPointTraj_I(MaxSat)
  real, allocatable         :: XyzSat_DII(:, :, :)
  real, public, allocatable :: TimeSat_II(:, :)

  ! Unlike ModSatellite in BATSRUS, here the saveplot frequencyis controlled
  ! by #NOUTPUT in SP_ModPlot so we do not read StartTime/EndTime/DtTraj.

  ! Variables for interpolation in a triangular mesh
  logical, public, allocatable :: IsTriangleFoundSat_I(:) ! if we find tri.
  integer, public, allocatable :: iStencilOrigSat_II(:,:) ! orig. sat stencils
  real,    public, allocatable :: WeightSat_II(:,:) ! weights for interpolation
contains
  !============================================================================
  subroutine read_param(NameCommand)

    use ModReadParam, ONLY: read_var

    character(len=*), intent(in) :: NameCommand
    logical :: DoTest
    integer :: iSat                         ! loop variable for satellites
    ! VARs for SATELLITE, specified to save the VDF along trajectories
    integer :: l1, l2                       ! indices to get satellite name
    character(len=*), parameter:: NameSub = 'read_param'
    !--------------------------------------------------------------------------
    call test_start(NameSub, DoTest)

    select case(NameCommand)
    case("#SATELLITE")
       ! Read the number of satellites, for saving the spectrum by DistrTraj_
       call read_var('nSatellite', nSat)
       if(nSat <= 0) RETURN
       UseSatellite = .true.
       if(nSat > MaxSat) call CON_stop(NameSub // ': Number of ' // &
            'output files is too large in #SATELLITE: nSat > MaxSat')

       ! Read satellite input file name and set the satellite name
       do iSat = 1, nSat
          call read_var('NameTrajectoryFile', NameFileSat_I(iSat))
          l1 = index(NameFileSat_I(iSat), '/', back=.true.) + 1
          l2 = index(NameFileSat_I(iSat), '.', back=.true.) - 1
          if(l1-1 <= 0) l1 = 1
          if(l2+1 <= 0) l2 = len_trim(NameFileSat_I(iSat))
          NameSat_I(iSat) = NameFileSat_I(iSat)(l1:l2)
       end do
    case default
       call CON_stop(NameSub//' unknown command='//NameCommand)
    end select

    call test_stop(NameSub, DoTest)
  end subroutine read_param
  !============================================================================
  subroutine init
    ! Initialize arrays for satellite triangulation

    use SP_ModSize, ONLY: nDim
    !--------------------------------------------------------------------------
    if(.not.UseSatellite) RETURN

    if(allocated(IsTriangleFoundSat_I)) deallocate(IsTriangleFoundSat_I)
    allocate(IsTriangleFoundSat_I(nSat))
    IsTriangleFoundSat_I = .false.
    if(allocated(iStencilOrigSat_II)) deallocate(iStencilOrigSat_II)
    allocate(iStencilOrigSat_II(nDim, nSat))
    iStencilOrigSat_II = 0
    if(allocated(WeightSat_II)) deallocate(WeightSat_II)
    allocate(WeightSat_II(nDim, nSat))
    WeightSat_II = 0.0

  end subroutine init
  !============================================================================
  subroutine read_satellite_input_files

    use CON_axes, ONLY: transform_matrix
    use ModTimeConvert, ONLY: time_int_to_real
    use ModIoUnit, ONLY: UnitTmp_, STDOUT_
    use ModKind, ONLY: Real8_
    use ModMpi
    use SP_ModTime, ONLY: StartTime

    ! One line of input
    character(len=100) :: StringLine

    ! Local VARs
    integer            :: nPoint                ! count of location points
    integer            :: MaxPoint              ! maximum location points
    integer            :: i, iSat               ! loop variables
    integer            :: iOffset               ! last point in which time<0
    integer            :: iTime_I(7)            ! time: YYYY/MM/DD-hh:mm:ss:ms
    real               :: Xyz_D(nDim)           ! location coordinates
    real(Real8_)       :: DateTime              ! date and time
    real, allocatable  :: Time_I(:), Xyz_DI(:,:)! time and location arrays
    character(len=100) :: NameFile              ! readin filenames
    character(len=3)   :: TypeCoordSat = 'HGI'  ! 'HGI' in TRAJECTORY files

    logical:: DoTest
    character(len=*), parameter:: NameSub = 'read_satellite_input_files'
    !--------------------------------------------------------------------------
    call test_start(NameSub, DoTest)

    ! Count maximum number of points by reading all satellite files
    MaxPoint = 0
    if(iProc == 0) then
       SATELLITES1: do iSat = 1, nSat

          ! Open the trajectory files
          NameFile = NameFileSat_I(iSat)
          call open_file(file=NameFile, status="old")
          nPoint = 0

          ! Count the number of points
          READFILE1: do
             read(UnitTmp_, '(a)', iostat=iError) StringLine
             if (iError /= 0) EXIT READFILE1

             if(index(StringLine, '#START')>0) then
                READPOINTS1: do
                   read(UnitTmp_, *, iostat=iError) iTime_I, Xyz_D
                   if(iError /= 0) EXIT READFILE1
                   ! Add new point
                   nPoint = nPoint + 1
                end do READPOINTS1
             end if
          end do READFILE1

          call close_file
          MaxPoint = max(MaxPoint, nPoint)
       end do SATELLITES1
    end if

    ! Tell all processors the maximum number of points
    call MPI_BCAST(MaxPoint, 1, MPI_INTEGER, 0, iComm, iError)

    ! allocate arrays depending on number of points
    allocate(Time_I(MaxPoint), Xyz_DI(nDim, MaxPoint))
    allocate(XyzSat_DII(nDim, MaxPoint, nSat))
    allocate(TimeSat_II(MaxPoint, nSat))

    ! Read the trajectories
    SATELLITES: do iSat = 1, nSat
       ! Read file on the root processor
       if(iProc == 0) then

          NameFile = NameFileSat_I(iSat)
          if(lVerbose > 0) &
               write(STDOUT_, *) NameSub, " reading: ", trim(NameFile)

          call open_file(file=NameFile, status="old")
          nPoint  = 0
          iOffset = 0
          ! Read the file: read #COOR TypeCoord, #START and points
          READFILE: do
             read(UnitTmp_,'(a)', iostat=iError) StringLine
             if(iError /= 0) EXIT READFILE

             if(index(StringLine, '#COOR')>0) &
                  read(UnitTmp_,'(a)') TypeCoordSat

             if(index(StringLine, '#START')>0) then
                READPOINTS: do
                   read(UnitTmp_,*, iostat=iError) iTime_I, Xyz_D
                   if (iError /= 0) EXIT READFILE
                   ! Add new point
                   nPoint = nPoint + 1
                   ! Store coordinates
                   Xyz_DI(:, nPoint) = Xyz_D
                   ! Convert integer date/time to simulation time
                   call time_int_to_real(iTime_I, DateTime)
                   Time_I(nPoint) = DateTime - StartTime
                   if(Time_I(nPoint) < 0.0) iOffset = nPoint
                enddo READPOINTS
             end if
          end do READFILE

          call close_file
          if(DoTest) write(*,*) NameSub//': nPoint=', nPoint, &
               ' iOffset=', iOffset
          ! Meaningfull part of trajectory ranges from iOffset to nPoint:
          if(iOffset > 1) then
             ! Temporary: new nPoint
             i = nPoint - iOffset + 1
             Xyz_DI(:, 1:i) = Xyz_DI(:, iOffset:nPoint)
             Time_I(1:i) = Time_I(iOffset:nPoint)
             nPoint = i
          end if
          ! Convert the coordinates if necessary (Out /= In)
          if(TypeCoordSystem /= TypeCoordSat) then
             do i = 1, nPoint
                Xyz_DI(:, i) = matmul(transform_matrix(Time_I(i), &
                     TypeCoordSat, TypeCoordSystem), Xyz_DI(:,i))
             end do
          end if

       end if

       ! Tell the number of points to the other processors
       call MPI_BCAST(nPoint, 1, MPI_INTEGER, 0, iComm, iError)
       nPointTraj_I(iSat) = nPoint

       ! Tell the other processors the satellite time
       call MPI_BCAST(Time_I, nPoint, MPI_REAL, 0, iComm, iError)

       ! Tell the other processors the coordinates
       call MPI_BCAST(Xyz_DI, nDim*nPoint, MPI_REAL, 0, iComm, iError)

       ! Store time and positions for satellite iSat on all PE-s
       TimeSat_II(1:nPoint, iSat) = Time_I(1:nPoint)
       XyzSat_DII(:, 1:nPoint, iSat) = Xyz_DI(:, 1:nPoint)
       iPointCurrentSat_I(1:nSat) = 1

       if(DoTest) then
          nPoint = min(10, nPoint)
          write(*,*) NameSub,': tSat=', TimeSat_II(1:nPoint, iSat)
          write(*,*) NameSub,': xSat=', XyzSat_DII(1, 1:nPoint, iSat)
          write(*,*) NameSub,': ySat=', XyzSat_DII(2, 1:nPoint, iSat)
          write(*,*) NameSub,': zSat=', XyzSat_DII(3, 1:nPoint, iSat)
       end if
    end do SATELLITES

    deallocate(Time_I, Xyz_DI)
    call test_stop(NameSub, DoTest)
  end subroutine read_satellite_input_files
  !============================================================================
  subroutine set_satellite_positions(iSat)

    use ModNumConst, ONLY: cTiny
    use SP_ModTime, ONLY: tSimulation => SPTime

    integer, intent(in) :: iSat
    integer :: i, nPoint
    real    :: dTime

    logical :: DoTest
    character(len=*), parameter:: NameSub = 'set_satellite_positions'
    !--------------------------------------------------------------------------
    call test_start(NameSub, DoTest)

    nPoint = nPointTraj_I(iSat)
    if(nPoint > 0) then
       ! Get the index for current time and satellite
       i = iPointCurrentSat_I(iSat)
       if(DoTest) write(*,*) NameSub, ' nPoint, iPoint, TimeSim, TimeSat=', &
            nPoint, i, tSimulation, TimeSat_II(i, iSat)

       do while (i < nPoint .and. TimeSat_II(i, iSat) <= tSimulation)
          i = i + 1
       end do
       iPointCurrentSat_I(iSat) = i
       if(DoTest) write(*,*) NameSub, ' final iPoint=', i

       ! if the index is out of range: do not track this satellite;
       ! otherwise, make a linear interpolation in time for the location
       if( (i == nPoint .and. tSimulation > TimeSat_II(i, iSat))  &
            .or. i == 1 ) then
          DoTrackSatellite_I(iSat) = .false.
          XyzSat_DI(:, iSat) = 0.0
       else
          DoTrackSatellite_I(iSat) = .true.
          dTime = 1.0 - (TimeSat_II(i, iSat) - tSimulation) /     &
               max(TimeSat_II(i, iSat) - TimeSat_II(i-1, iSat), cTiny)
          XyzSat_DI(:, iSat) = dTime * XyzSat_DII(:, i, iSat) +   &
               (1.0 - dTime) * XyzSat_DII(:, i-1, iSat)
       end if

       if(DoTest) then
          write(*,*) NameSub, ' DoTrackSat =', DoTrackSatellite_I(iSat)
          write(*,*) NameSub, ' XyzSat     =', XyzSat_DI(:, iSat)
       end if
    end if

    call test_stop(NameSub, DoTest)
  end subroutine set_satellite_positions
  !============================================================================
end module SP_ModSatellite
!==============================================================================
