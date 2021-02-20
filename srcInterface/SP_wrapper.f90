!  Copyright (C) 2002 Regents of the University of Michigan,
!  portions used with permission
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module SP_wrapper
  use SP_ModSize, ONLY: nParticleMax
  use SP_ModMain, ONLY: run, DoRestart, DoReadMhData, DataInputTime
  implicit none
  save
  private ! except
  public:: SP_init_session
  public:: SP_set_param
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
  ! Interface routines to be called from super-structure only
  subroutine SP_do_extract_lines(DoExtract)
    use CON_coupler, ONLY: i_proc0, i_comm, is_proc0, SP_
    use ModMpi
    logical, intent(out):: DoExtract

    integer :: iError

    ! when restarting, line data is available, i.e. ready to couple with mh;
    ! get value at SP root and broadcast to all SWMF processors
    character(len=*), parameter:: NameSub = 'SP_check_ready_for_mh'
    !--------------------------------------------------------------------------
    if(is_proc0(SP_)) DoExtract = .not.DoRestart
    call MPI_Bcast(DoExtract, 1, MPI_LOGICAL, i_proc0(SP_), i_comm(), iError)
  end subroutine SP_do_extract_lines
  !============================================================================
  ! Above routines may be called from superstructure only.
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
  subroutine SP_init_session(iSession,TimeSimulation)
    use SP_ModMain, ONLY: initialize, DoRestart, DoReadMhData
    use SP_ModGrid, ONLY: init_grid=>init, nVar, State_VIB, MHData_VIB
    use SP_ModGrid, ONLY: nLon, nLat, nParticle_B, FootPoint_VB
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
    nullify(nParticle_B);  nullify(FootPoint_VB)

    call BL_init(nParticleMax, nLon, nLat,  &
         MHData_VIB, nParticle_B, nVar, State_VIB, FootPoint_VB)
    call initialize

    if(.not.(DoRestart.or.DoReadMhData))&
         call BL_get_origin_points(ROrigin, LonMin, LonMax, LatMin, LatMax)
  end subroutine SP_init_session
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
  subroutine SP_set_param(CompInfo,TypeAction)
    use CON_comp_info
    use SP_ModTime,  ONLY: StartTimeJulian, StartTime, SPTime,&
         time_real_to_julian
    use SP_ModMain,  ONLY: check, read_param
    use SP_ModGrid,  ONLY: TypeCoordSystem
    use CON_coupler, ONLY: SP_, is_proc
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
  !============================================================================
  subroutine SP_save_restart(TimeSimulation)
    use SP_ModMain, ONLY: save_restart
    real, intent(in) :: TimeSimulation

    ! if data are read from files, no need for additional run
    !--------------------------------------------------------------------------
    if(.not.DoReadMhData)call run(TimeSimulation)
    call save_restart
  end subroutine SP_save_restart
  !============================================================================
  subroutine SP_put_coupling_param(iModelIn, rMinIn, rMaxIn, TimeIn,&
       rBufferLoIn, rBufferUpIn)
    use SP_ModMain, ONLY: copy_old_state
    use CON_bline,  ONLY: get_bounds, Lower_
    integer,        intent(in) :: iModelIn
    real,           intent(in) :: rMinIn, rMaxIn
    real,           intent(in) :: TimeIn
    real, optional, intent(in) :: rBufferLoIn, rBufferUpIn
    !--------------------------------------------------------------------------
    call get_bounds(iModelIn, rMinIn, rMaxIn, &
       rBufferLoIn, rBufferUpIn)
    if(DataInputTime >= TimeIn)RETURN
    ! New coupling time, get it and save old state
    DataInputTime = TimeIn
    if(iModelIn==Lower_)then
       call copy_old_state
    else
       call CON_stop("Time in IHSP coupling differs from that in SCSP")
    end if
  end subroutine SP_put_coupling_param
  !============================================================================
  ! Called from coupler after the updated grid point lo<cation are
  ! received from the other component (SC, IH). Determines whether some
  ! grid points should be added/deleted
  subroutine SP_adjust_lines(DoInit)
    use SP_ModDistribution, ONLY: offset
    use SP_ModGrid,         ONLY: R_, State_VIB, nBlock, MHData_VIB, &
         nParticle_B, Length_, FootPoint_VB
    use CON_bline,          ONLY: iOffset_B, rBufferLo, rBufferUp, &
         rInterfaceMin, rInterfaceMax, rMin, adjust_line, LagrID_, &
         X_, Z_, SP_set_line_foot_b
    logical, intent(in) :: DoInit
    ! once new geometry of lines has been put, account for some particles
    ! exiting the domain (can happen both at the beginning and the end)
    integer:: iParticle, iBlock, iBegin,  iEnd, iOffset ! loop variables

    logical:: DoAdjustLo, DoAdjustUp

    character(len=100):: StringError

    character(len=*), parameter:: NameSub = 'SP_adjust_lines'
    !--------------------------------------------------------------------------
    DoAdjustLo = RBufferLo == RInterfaceMin
    DoAdjustUp = RBufferUp == RInterfaceMax
    if(DoInit.and.DoAdjustLo)then
       do iBlock = 1, nBlock
          call SP_set_line_foot_b(iBlock)
       end do
    end if
    BLOCK:do iBlock = 1, nBlock
       call adjust_line(iBlock, iBegin, DoAdjustLo, DoAdjustUp)
       if(iBegin/=1)then
          ! Offset particle arrays
          iEnd   = nParticle_B(iBlock)
          iOffset = 1 - iBegin
          iOffset_B(iBlock) = iOffset
          MHData_VIB(LagrID_:Z_, 1:iEnd+iOffset, iBlock) = &
               MHData_VIB(LagrID_:Z_,iBegin:iEnd,iBlock)
          nParticle_B(iBlock) = nParticle_B(iBlock) + iOffset
          ! need to recalculate footpoints
          call SP_set_line_foot_b(iBlock)
          call offset(iBlock, iOffset)
       end if
    end do BLOCK
    ! may need to add particles to the beginning of lines
    if(DoAdjustLo) call append_particles
    ! Called after the grid points are received from the
    ! component, nullify offset. Alternatively, if the points
    ! are received from SC via initial coupling, there is no
    ! need to apply offset in IH, because there are no points
    ! in IH yet and no need to correct their ID
    if(DoAdjustUp.or.DoInit)iOffset_B(1:nBlock) = 0
  contains

    subroutine append_particles
      use SP_ModGrid, ONLY: R_, State_VIB
      ! appends a new particle at the beginning of lines if necessary
      integer:: iBlock
      real:: DistanceToMin
      real, parameter:: cTol = 1E-06

      character(len=*), parameter:: NameSub = 'append_particles'
      !------------------------------------------------------------------------
      BLOCK:do iBlock = 1, nBlock
         ! check current value of offset: if not zero, adjustments have just
         ! been made, no need to append new particles
         if(iOffset_B(iBlock) /= 0 )CYCLE BLOCK
         ! check if the beginning of the line moved far enough from its
         ! footprint on the solar surface
         DistanceToMin = sqrt(sum((&
              MHData_VIB(X_:Z_,1,iBlock) - FootPoint_VB(X_:Z_,iBlock))**2))
         ! skip the line if it's still close to the Sun
         if(DistanceToMin*(1.0 + cTol) < FootPoint_VB(Length_, iBlock))&
              CYCLE BLOCK
         ! append a new particle
         ! check if have enough space
         if(nParticleMax == nParticle_B( iBlock))call CON_Stop(NameSub//&
              ': not enough memory allocated to append a new particle')
         ! Particles ID as handled by other components keep unchanged
         ! while their order numbers in SP are increased by 1. Therefore,
         iOffset_B(iBlock)  = 1
         MHData_VIB(       LagrID_:Z_,2:nParticle_B(iBlock) + 1, iBlock)&
              = MHData_VIB(LagrID_:Z_,1:nParticle_B(iBlock),     iBlock)
         State_VIB(       R_        ,2:nParticle_B(iBlock) + 1, iBlock)&
              = State_VIB(R_        ,1:nParticle_B(iBlock),     iBlock)
         nParticle_B(iBlock) = nParticle_B(iBlock) + 1
         ! put the new particle just above the lower boundary
         MHData_VIB(LagrID_:Z_,  1, iBlock) = &
              FootPoint_VB(LagrID_:Z_, iBlock)*(1.0 + cTol)
         State_VIB(R_,          1, iBlock) = &
              sqrt(sum((MHData_VIB(X_:Z_,  1, iBlock))**2))
         MHData_VIB(LagrID_,1, iBlock) = MHData_VIB(LagrID_, 2, iBlock) - 1.0
         FootPoint_VB(LagrID_,iBlock) = MHData_VIB(LagrID_, 1, iBlock) - 1.0
         call offset(iBlock, iOffset=iOffset_B(iBlock))
      end do BLOCK
    end subroutine append_particles
    !==========================================================================
  end subroutine SP_adjust_lines
  !============================================================================
end module SP_wrapper
