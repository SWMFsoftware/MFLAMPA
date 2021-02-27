!  Copyright (C) 2002 Regents of the University of Michigan,
!  portions used with permission
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module SP_wrapper
  use CON_coupler,ONLY: SP_
  use SP_ModSize, ONLY: nVertexMax
  use SP_ModMain, ONLY: run, DoRestart, DoReadMhData
  use SP_ModTime, ONLY: DataInputTime
  implicit none
  save
  private ! except
  public:: SP_set_param
  public:: SP_init_session
  public:: SP_run
  public:: SP_save_restart
  public:: SP_finalize

  ! coupling with MHD components
  public:: SP_do_extract_lines
  public:: SP_put_coupling_param
  public:: SP_adjust_lines

  logical :: DoCheck = .true.
contains
  !============================================================================
  subroutine SP_do_extract_lines(DoExtract)
    ! Interface routine to be called from super-structure on all PEs
    use CON_coupler, ONLY: i_proc0, i_comm, is_proc0
    use ModMpi
    logical, intent(out):: DoExtract

    integer :: iError

    ! when restarting, line data is available, i.e. ready to couple with mh;
    ! get value at SP root and broadcast to all SWMF processors
    character(len=*), parameter:: NameSub = 'SP_do_extract_lines'
    !--------------------------------------------------------------------------
    if(is_proc0(SP_)) DoExtract = .not.DoRestart
    call MPI_Bcast(DoExtract, 1, MPI_LOGICAL, i_proc0(SP_), i_comm(), iError)
  end subroutine SP_do_extract_lines
  !============================================================================
  ! Interface routine to be called from super-structure on all PEs
  subroutine SP_set_param(CompInfo,TypeAction)
    use CON_comp_info
    use SP_ModTime,  ONLY: StartTimeJulian, StartTime, SPTime,&
         time_real_to_julian
    use SP_ModMain,  ONLY: check, read_param
    use SP_ModGrid,  ONLY: TypeCoordSystem
    use CON_coupler, ONLY: is_proc
    use SP_ModProc
    use SP_ModUnit,  ONLY: Si2Io_V, Io2Si_V, UnitX_, UnitEnergy_
    use CON_physics, ONLY: get_time
    use CON_bline,   ONLY: BL_set_grid
    type(CompInfoType),intent(inout):: CompInfo
    character(len=*),  intent(in)   :: TypeAction
    real :: UnitX, EnergyCoeff

    character(len=*), parameter:: NameSub = 'SP_set_param'
    !--------------------------------------------------------------------------
    select case(TypeAction)
    case('VERSION')
       call put(CompInfo,&
            Use        =.true., &
            NameVersion='MFLAMPA', &
            Version    =0.90)
    case('MPI')
       call get(CompInfo, iComm=iComm, iProc=iProc, nProc=nProc)
    case('STDOUT')
       ! placeholder
    case('CHECK')
       if(.not.DoCheck)RETURN
       DoCheck = .false.
       call get_time(tSimulationOut = SPTime, tStartOut = StartTime)
       call time_real_to_julian(StartTime, StartTimeJulian)
       DataInputTime = SPTime
       call check
    case('READ')
       call read_param
    case('GRID')
       if(is_proc(SP_))then
          UnitX = Io2Si_V(UnitX_)
          EnergyCoeff = Si2Io_V(UnitEnergy_)
       end if
       call BL_set_grid(TypeCoordSystem, UnitX, EnergyCoeff)
    case default
       call CON_stop('Cannot call SP_set_param for '//trim(TypeAction))
    end select
  end subroutine SP_set_param
  ! Above two routines are be called from superstructure  on all PEs
  !============================================================================
  ! The following routines are called  on  PEs of the SP model only
  subroutine SP_init_session(iSession,TimeSimulation)
    use SP_ModMain, ONLY: initialize, DoRestart, DoReadMhData
    use SP_ModGrid, ONLY: init_grid=>init, nVar, State_VIB, MHData_VIB
    use SP_ModGrid, ONLY: nLon, nLat, nVertex_B, FootPoint_VB
    use CON_bline,  ONLY: BL_init, BL_get_origin_points
    use SP_ModOriginPoints, ONLY: ROrigin, LonMin, LonMax, LatMin, LatMax
    integer,  intent(in) :: iSession         ! session number (starting from 1)
    real,     intent(in) :: TimeSimulation   ! seconds from start time

    logical, save:: IsInitialized = .false.

    !--------------------------------------------------------------------------
    if(IsInitialized)RETURN
    IsInitialized = .true.
    
    call init_grid

    nullify(MHData_VIB) ;  nullify(State_VIB)
    nullify(nVertex_B);  nullify(FootPoint_VB)
    call BL_init(nVertexMax, nLon, nLat,  &
         MHData_VIB, nVertex_B, nVar, State_VIB, FootPoint_VB)
    call initialize

    if(.not.(DoRestart.or.DoReadMhData))&
         call BL_get_origin_points(ROrigin, LonMin, LonMax, LatMin, LatMax)
  end subroutine SP_init_session
  !============================================================================
  subroutine SP_run(TimeSimulation,TimeSimulationLimit)
    real,intent(inout)::TimeSimulation
    real,intent(in)::TimeSimulationLimit

    !--------------------------------------------------------------------------
    call run(TimeSimulationLimit)
    if(DoReadMhData)then
       TimeSimulation = DataInputTime
    else
       TimeSimulation = TimeSimulationLimit
    end if
  end subroutine SP_run
  !============================================================================
  subroutine SP_finalize(TimeSimulation)
    use SP_ModMain, ONLY: finalize
    real,intent(in)::TimeSimulation
    ! if data are read from files, no special finalization is needed
    !--------------------------------------------------------------------------
    if(.not.DoReadMhData)call run(TimeSimulation)
    call finalize
  end subroutine SP_finalize
  !============================================================================
  subroutine SP_save_restart(TimeSimulation)
    use SP_ModRestart, ONLY: save_restart
    real, intent(in) :: TimeSimulation

    ! if data are read from files, no need for additional run
    !--------------------------------------------------------------------------
    if(.not.DoReadMhData)call run(TimeSimulation)
    call save_restart
  end subroutine SP_save_restart
  !============================================================================
  subroutine SP_put_coupling_param(Source_, TimeIn)
    use SP_ModGrid, ONLY: copy_old_state
    use CON_bline,  ONLY: Lower_
    integer,        intent(in) :: Source_
    real,           intent(in) :: TimeIn
    !--------------------------------------------------------------------------
    if(DataInputTime >= TimeIn)RETURN
    ! New coupling time, get it and save old state
    DataInputTime = TimeIn
    if(Source_==Lower_)then
       call copy_old_state
    else
       call CON_stop("Time in IHSP coupling differs from that in SCSP")
    end if
  end subroutine SP_put_coupling_param
  !============================================================================
  ! Called from coupler after the updated grid point lo<cation are
  ! received from the other component (SC, IH). Determines whether some
  ! grid points should be added/deleted
  subroutine SP_adjust_lines(DoInit, Source_)
    use SP_ModDistribution, ONLY: offset
    use CON_bline,          ONLY: &
         iOffset_B, BL_adjust_lines,  Lower_, Upper_, nLine
    logical, intent(in) :: DoInit
    integer, intent(in) :: Source_
    integer:: iLine  ! loop variable
    character(len=*), parameter:: NameSub = 'SP_adjust_lines'
    !--------------------------------------------------------------------------
    call BL_adjust_lines(DoInit, Source_)
    if(Source_ == Lower_)then
       do iLine = 1, nLine
          !Offset the distribution function array, if needed
          call offset(iLine, iOffset=iOffset_B(iLine))
       end do
    end if
    ! Called after the grid points are received from the
    ! component, nullify offset. 
    if(Source_ == Upper_)iOffset_B(1:nLine) = 0
  end subroutine SP_adjust_lines
  !============================================================================
end module SP_wrapper
