!  Copyright (C) 2002 Regents of the University of Michigan,
!  portions used with permission
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module SP_ModSatellite

  use ModUtilities,   ONLY: open_file, close_file, CON_stop
  use SP_ModProc,     ONLY: iProc, iComm
  use SP_ModTime,     ONLY: StartTime
  use SP_ModTestFunc, ONLY: lVerbose, test_start, test_stop

  implicit none
  save
  private ! Except

  public:: read_param                 ! read satellite file input parameters
  public:: read_satellite_input_files ! read satellite trajectories
  public:: set_satellite_positions    ! set satellite positions

  integer, public    :: nSat   = 0    ! number of satellites
  integer, parameter :: MaxSat = 300  ! Max number of satellites

  ! These variables are public for write_logfile only !!! Should be improved
  ! Names and unit numbers for satellite files
  character(len=50), public:: NameFileSat_I(MaxSat)
  ! Names of the satellite
  character(len=50), public:: NameSat_I(MaxSat)

  ! Coordinate system
  character(len=3), public :: TypeSatCoord_I(MaxSat)
  ! Current positions
  real, public:: XyzSat_DI(3, MaxSat)

  ! variables to record tracked and current satellite position indices
  logical, public:: UseSatellite               = .false.
  logical, public:: DoTrackSatellite_I(MaxSat) = .false.
  integer, public:: iPointCurrentSat_I(MaxSat) = 1

  ! Local variables
  logical:: UseSatFile_I(MaxSat) = .true.
  integer, public           :: nPointTraj_I(MaxSat)
  real, allocatable         :: XyzSat_DII(:, :, :)
  real, public, allocatable :: TimeSat_II(:, :)

  ! Unlike ModSatellite in BASTRUS, here the saveplot frequencyis controlled
  ! by #NOUTPUT in SP_ModPlot so we do not read StartTime/EndTime/DtTraj.
contains
  !============================================================================
  subroutine read_param(NameCommand)

    use SP_ModIO,     ONLY: nFile, MaxFile, Satellite_, TypeCoordPlot_I
    use ModUtilities, ONLY: check_dir
    use ModReadParam, ONLY: read_var

    character(len=*), intent(in) :: NameCommand
    integer           :: iSat, iFile
    character(len=100):: StringSatellite
    integer           :: l1, l2

    logical:: DoTest
    character(len=*), parameter:: NameSub = 'read_param'
    !--------------------------------------------------------------------------
    call test_start(NameSub, DoTest)

    select case(NameCommand)
    case("#SATELLITE")
       call read_var('nSatellite', nSat)
       if(nSat <= 0) RETURN

       UseSatellite = .true.
       nFile = max(nFile, Satellite_ + nSat)
       if (nFile > MaxFile .or. nSat > MaxSat)&
            call CON_stop('The number of output files is too large ' &
            //'in #SATELLITE: nFile > MaxFile .or. nSat > MaxSat')

       do iSat = 1, nSat
          iFile = Satellite_ + iSat
          call read_var('StringSatellite', StringSatellite)

          ! Read satellite input file name and set the satellite name
          call read_var('NameTrajectoryFile', NameFileSat_I(iSat))
          l1 = index(NameFileSat_I(iSat), '/', back=.true.) + 1
          l2 = index(NameFileSat_I(iSat), '.', back=.true.) - 1
          if(l1-1 <= 0) l1 = 1
          if(l2+1 <= 0) l2 = len_trim(NameFileSat_I(iSat))
          NameSat_I(iSat) = NameFileSat_I(iSat)(l1:l2)

          ! whether to use the satellite files
          if(index(StringSatellite,'eqn') > 0 &
               .or. index(StringSatellite,'Eqn') > 0 &
               .or. index(StringSatellite,'EQN') > 0 ) then
             UseSatFile_I(iSat) = .false.
          else
             UseSatFile_I(iSat) = .true.
          end if

          ! Recognize coordinate system name if present
          if (index(StringSatellite,'GEO') > 0) TypeCoordPlot_I(iFile) = 'GEO'
          if (index(StringSatellite,'GSE') > 0) TypeCoordPlot_I(iFile) = 'GSE'
          if (index(StringSatellite,'GSM') > 0) TypeCoordPlot_I(iFile) = 'GSM'
          if (index(StringSatellite,'MAG') > 0) TypeCoordPlot_I(iFile) = 'MAG'
          if (index(StringSatellite,'SMG') > 0) TypeCoordPlot_I(iFile) = 'SMG'
          if (index(StringSatellite,'HGR') > 0) TypeCoordPlot_I(iFile) = 'HGR'
          if (index(StringSatellite,'HGI') > 0) TypeCoordPlot_I(iFile) = 'HGI'
          if (index(StringSatellite,'HGC') > 0) TypeCoordPlot_I(iFile) = 'HGC'
       end do

    case default
       call CON_stop(NameSub//' unknown command='//NameCommand)
    end select

    call test_stop(NameSub, DoTest)

  end subroutine read_param
  !============================================================================
  subroutine read_satellite_input_files

    use CON_axes,       ONLY: transform_matrix
    use ModTimeConvert, ONLY: time_int_to_real
    use ModIoUnit,      ONLY: UnitTmp_
    use ModKind,        ONLY: Real8_
    use ModMpi
    use SP_ModIO,       ONLY: iUnitOut, write_prefix
    use SP_ModTime,     ONLY: StartTime

    integer, parameter :: MaxDim = 3
    integer            :: iError, i, iSat, nPoint

    ! One line of input
    character(len=100) :: StringLine

    ! Local VARs
    integer            :: iTime_I(7)
    real               :: Xyz_D(MaxDim)
    real(Real8_)       :: DateTime
    integer            :: MaxPoint
    real, allocatable  :: Time_I(:), Xyz_DI(:,:)
    character(len=100) :: NameFile
    character(len=3)   :: TypeCoordSystem = 'HGR'  ! 'HGR' in SP as default

    logical:: DoTest
    character(len=*), parameter:: NameSub = 'read_satellite_input_files'
    !--------------------------------------------------------------------------
    call test_start(NameSub, DoTest)

    ! Count maximum number of points by reading all satellite files
    MaxPoint = 0
    if(iProc == 0) then
       SATELLITES1: do iSat = 1, nSat
          if(.not.UseSatFile_I(iSat)) CYCLE SATELLITES1
          NameFile = NameFileSat_I(iSat)
          call open_file(file=NameFile, status="old")
          nPoint = 0

          TypeSatCoord_I(iSat) = TypeCoordSystem
          READFILE1: do
             read(UnitTmp_, '(a)', iostat=iError) StringLine
             if (iError /= 0) EXIT READFILE1

             if(index(StringLine,'#START')>0) then
                READPOINTS1: do
                   read(UnitTmp_, *, iostat=iError) iTime_I, Xyz_D
                   if (iError /= 0) EXIT READFILE1
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
    call MPI_Bcast(MaxPoint, 1, MPI_INTEGER, 0, iComm, iError)

    ! allocate arrays depending on number of points
    allocate(Time_I(MaxPoint), Xyz_DI(MaxDim, MaxPoint))
    allocate(XyzSat_DII(3, nSat, MaxPoint))
    allocate(TimeSat_II(nSat, MaxPoint))

    ! Read the trajectories
    SATELLITES: do iSat = 1, nSat

       if(.not.UseSatFile_I(iSat)) CYCLE SATELLITES

       ! Read file on the root processor
       if(iProc == 0) then

          NameFile = NameFileSat_I(iSat)
          if(lVerbose > 0)then
             call write_prefix
             write(iUnitOut, *) NameSub, " reading: ", trim(NameFile)
          end if

          call open_file(file=NameFile, status="old")
          nPoint = 0

          ! Read the file: read #COOR TypeCoord, #START and points
          ! Default coordinate system is the one used by BATSRUS (or GSM?)
          TypeSatCoord_I(iSat) = TypeCoordSystem
          READFILE: do
             read(UnitTmp_,'(a)', iostat=iError) StringLine
             if(iError /= 0) EXIT READFILE

             if(index(StringLine,'#COOR')>0) &
                  read(UnitTmp_,'(a)') TypeSatCoord_I(iSat)

             if(index(StringLine,'#START')>0) then
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
                enddo READPOINTS
             endif
          enddo READFILE

          call close_file
          if(DoTest) write(*,*) NameSub, ': nPoint=', nPoint

          ! Convert the coordinates if necessary
          if(TypeSatCoord_I(iSat) /= TypeCoordSystem)then
             do i = 1, nPoint
                Xyz_DI(:, i) = matmul(transform_matrix(Time_I(i), &
                     TypeSatCoord_I(iSat), TypeCoordSystem), Xyz_DI(:,i))
             end do
          end if

       end if

       ! Tell the number of points to the other processors
       call MPI_Bcast(nPoint, 1, MPI_INTEGER, 0, iComm, iError)
       nPointTraj_I(iSat) = nPoint

       ! Tell the other processors the satellite time
       call MPI_Bcast(Time_I, nPoint, MPI_REAL, 0, iComm, iError)

       ! Tell the other processors the coordinates
       call MPI_Bcast(Xyz_DI, MaxDim*nPoint, MPI_REAL, 0, iComm, iError)

       ! Store time and positions for satellite iSat on all PE-s
       TimeSat_II(iSat, 1:nPoint) = Time_I(1:nPoint)
       do i = 1, nPoint
          XyzSat_DII(:, iSat, i) = Xyz_DI(:, i)
       end do

       if(DoTest) then
          nPoint = min(10, nPoint)
          write(*,*) NameSub,': tSat=', TimeSat_II(iSat,1:nPoint)
          write(*,*) NameSub,': xSat=', XyzSat_DII(1,iSat,1:nPoint)
          write(*,*) NameSub,': ySat=', XyzSat_DII(2,iSat,1:nPoint)
          write(*,*) NameSub,': zSat=', XyzSat_DII(3,iSat,1:nPoint)
       end if

    end do SATELLITES

    deallocate(Time_I, Xyz_DI)
    call test_stop(NameSub, DoTest)

  end subroutine read_satellite_input_files
  !============================================================================
  subroutine set_satellite_positions(iSat)

    use ModNumConst
    use SP_ModTime, ONLY: tSimulation => SPTime

    integer, intent(in) :: iSat
    integer :: i, nPoint
    real    :: dTime

    logical :: DoTest
    character(len=*), parameter:: NameSub = 'set_satellite_positions'
    !--------------------------------------------------------------------------
    call test_start(NameSub, DoTest)

    if(UseSatFile_I(iSat)) then

       nPoint = nPointTraj_I(iSat)
       if(nPoint > 0) then

          i = iPointCurrentSat_I(iSat)
          if(DoTest) write(*,*) NameSub,  &
               ' nPoint, iPoint, TimeSim, TimeSat=',   &
               nPoint, i, tSimulation, TimeSat_II(iSat,i)

          do while (i < nPoint .and. TimeSat_II(iSat, i) <= tSimulation)
             i = i + 1
          end do
          iPointCurrentSat_I(iSat) = i
          if(DoTest) write(*,*) NameSub, ' final iPoint=', i

          if( (i == nPoint .and. tSimulation > TimeSat_II(iSat, i))  &
               .or. i == 1 ) then
             DoTrackSatellite_I(iSat) = .false.
             XyzSat_DI(:, iSat) = 0.0
          else
             DoTrackSatellite_I(iSat) = .true.
             dTime = 1.0 - (TimeSat_II(iSat, i) - tSimulation) /     &
                  max(TimeSat_II(iSat, i) - TimeSat_II(iSat, i-1), cTiny)
             XyzSat_DI(:, iSat) = dTime * XyzSat_DII(:, iSat, i) +   &
                  (1.0 - dTime) * XyzSat_DII(:, iSat, i-1)
          end if

          if(DoTest) then
             write(*,*) NameSub, ' DoTrackSat =', DoTrackSatellite_I(iSat)
             write(*,*) NameSub, ' XyzSat     =', XyzSat_DI(:,iSat)
          end if

       end if
    else
       call satellite_trajectory_formula(iSat)
    end if

    call test_stop(NameSub, DoTest)

  end subroutine set_satellite_positions
  !============================================================================
  subroutine satellite_trajectory_formula(iSat)

    integer, intent(in) :: iSat
    real    :: XyzSat_D(3)

    logical :: DoTest
    character(len=*), parameter:: NameSub = 'satellite_trajectory_formula'
    !--------------------------------------------------------------------------
    call test_start(NameSub, DoTest)

    XyzSat_D = XyzSat_DI(:,iSat)

    ! Case should be for a specific satellite.  The trajectories can depend
    ! on the 'real' time so that the satellite knows where it is at.  For
    ! example, Cassini could be if'd on the date so that the code knows
    ! whether the run is for a time near Earth, Jupiter or Saturn.

    ! This routine should always set the TrackSatellite(iSat) flag. When
    ! the satellite is at a useless position the time should return a
    ! do not track flag (DoTrackSatellite_I(iSat) = .false.).

    select case(NameFileSat_I(iSat))
    case ('earth', 'cassini')
       XyzSat_D = 5.0
       DoTrackSatellite_I(iSat) = .true.
    case default
       XyzSat_D = 1.0
       DoTrackSatellite_I(iSat) = .false.
    end select

    call test_stop(NameSub, DoTest)

  end subroutine satellite_trajectory_formula
  !============================================================================
end module SP_ModSatellite
!==============================================================================
