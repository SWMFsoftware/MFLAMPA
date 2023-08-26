!  Copyright (C) 2002 Regents of the University of Michigan
! portions used with permission
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module SP_ModMain

  use SP_ModProc,       ONLY: iProc
  use SP_ModReadMhData, ONLY: DoReadMhData
  use ModUtilities, ONLY: CON_stop

  implicit none

  SAVE
  private ! except

  ! Indicator of stand alone mode
  logical:: IsStandAlone=.false.

  ! Stopping conditions. These variables are only used in stand alone mode.
  real    :: TimeMax     = -1.0, CpuTimeMax = -1.0
  integer :: nIterMax    = -1
  logical :: UseStopFile = .true.
  logical :: IsLastRead  = .false.

  ! Logicals for actions
  !----------------------
  ! run the component
  logical:: DoRun = .true.
  ! restart the run
  logical:: DoRestart = .false.
  ! Methods and variables from ModReadMhData
  public :: DoReadMhData

  ! Methods and variables from this module
  public:: read_param, initialize, finalize, run, check, DoRestart,        &
       IsLastRead, UseStopFile, CpuTimeMax, TimeMax, nIterMax, IsStandAlone
contains
  !============================================================================
  subroutine read_param

    use SP_ModAdvance,       ONLY: DoTraceShock, UseDiffusion
    use SP_ModAdvance,       ONLY: read_param_adv        =>read_param
    use SP_ModAngularSpread, ONLY: read_param_spread     =>read_param
    use SP_ModDistribution,  ONLY: read_param_dist       =>read_param
    use SP_ModGrid,          ONLY: read_param_grid       =>read_param
    use SP_ModOriginPoints,  ONLY: read_param_origin     =>read_param
    use SP_ModPlot,          ONLY: read_param_plot       =>read_param
    use SP_ModReadMHData,    ONLY: read_param_mhdata     =>read_param
    use SP_ModRestart,       ONLY: read_param_restart    =>read_param
    use SP_ModTime,          ONLY: read_param_time       =>read_param
    use SP_ModTiming,        ONLY: read_param_timing     =>read_param
    use SP_ModTurbulence,    ONLY: read_param_turbulence =>read_param
    use SP_ModUnit,          ONLY: read_param_unit       =>read_param

    ! Read input parameters for SP component
    use ModReadParam, ONLY: &
         read_var, read_line, read_command, i_session_read, read_echo_set

    ! aux variables
    integer:: nParticleCheck, nLonCheck, nLatCheck
    logical:: DoEcho

    ! The name of the command
    character (len=100) :: NameCommand
    character(len=*), parameter:: NameSub = 'read_param'
    !--------------------------------------------------------------------------
    ! Read the corresponding section of input file
    do
       if(.not.read_line() ) then
          IsLastRead = .true.
          EXIT
       end if
       if(.not.read_command(NameCommand)) CYCLE
       select case(NameCommand)
       case('#RESTART')
          call read_var('DoRestart',DoRestart)
          ! read parameters for each module
       case('#ORIGIN')
          if(IsStandAlone)CYCLE
          call read_param_origin
       case('#COORDSYSTEM', '#COORDINATESYSTEM', '#TESTPOS', &
            '#CHECKGRIDSIZE','#DOSMOOTH', '#GRIDNODE')
          if(i_session_read() /= 1)CYCLE
          call read_param_grid(NameCommand)
       case('#MOMENTUMGRID','#FLUXINITIAL', '#FLUXCHANNEL')
          if(i_session_read() /= 1)CYCLE
          call read_param_dist(NameCommand)
       case('#INJECTION','#CFL','#USEFIXEDMFPUPSTREAM', '#DIFFUSION', &
            '#TRACESHOCK', '#SCALETURBULENCE')
          call read_param_adv(NameCommand)
       case('#SAVEPLOT','#USEDATETIME','#SAVEINITIAL','#NTAG')
          call read_param_plot(NameCommand)
       case('#READMHDATA','#MHDATA')
          call read_param_mhdata(NameCommand)
       case('#TURBULENTSPECTRUM', '#INITSPECTRUM', '#ADVECTIONWITHALFVENSPEED')
          call read_param_turbulence(NameCommand)
       case('#PARTICLEENERGYUNIT')
          call read_param_unit(NameCommand)
       case('#SPREADGRID', '#SPREADSOLIDANGLE', '#SPREADSIGMA')
          call read_param_spread(NameCommand)
       case('#DORUN')
          call read_var('DoRun',DoRun)
       case('#TIMING')
          call read_param_timing
       case('#TEST')
          call read_var('DoTraceShock', DoTraceShock)
          call read_var('UseDiffusion', UseDiffusion)
       case('#END')
          call check_stand_alone
          IsLastRead=.true.
          EXIT
       case('#RUN')
          call check_stand_alone
          IsLastRead=.false.
          EXIT
       case('#STOP')
          call check_stand_alone
          call read_var('nIterMax',nIterMax)
          call read_var('TimeMax' ,TimeMax)
       case('#CPUTIMEMAX')
          call check_stand_alone
          call read_var('CpuTimeMax',CpuTimeMax)
       case('#CHECKSTOPFILE')
          call check_stand_alone
          call read_var('UseStopFile',UseStopFile)
       case('#ECHO')
          call check_stand_alone
          call read_var('DoEcho', DoEcho)
          if(iProc==0)call read_echo_set(DoEcho)
       case('#NSTEP','#TIMESIMULATION')
          if(i_session_read() /= 1)CYCLE
          call read_param_time(NameCommand)
       case("#STARTTIME", "#SETREALTIME",'#TIMEACCURATE')
          if(i_session_read() /= 1)CYCLE
          call check_stand_alone
          call read_param_time(NameCommand)
       case('#SAVERESTART')
          if(i_session_read() /= 1)CYCLE
          call check_stand_alone
          call read_param_restart
       case default
          call CON_stop(NameSub//': Unknown command '//NameCommand)
       end select
    end do
  contains
    !==========================================================================
    subroutine check_stand_alone

      ! certain options are only available for stand alone mode;
      ! check whether the mode is used and stop the code if it's no the case
      !------------------------------------------------------------------------
      if(IsStandAlone)RETURN
      call CON_stop(NameSub//': command '//trim(NameCommand)//&
           ' is only allowed in stand alone mode, correct PARAM.in')
    end subroutine check_stand_alone
    !==========================================================================
  end subroutine read_param
  !============================================================================
  subroutine initialize

    use SP_ModAngularSpread, ONLY: init_spread     => init
    use SP_ModDistribution,  ONLY: init_dist       => init
    use SP_ModPlot,          ONLY: init_plot       => init
    use SP_ModReadMhData,    ONLY: init_mhdata     => init
    use SP_ModRestart,       ONLY: &
         read_restart
    use SP_ModTime,          ONLY: init_time       => init
    use SP_ModTurbulence,    ONLY: init_turbulence => init
    use SP_ModUnit,          ONLY: init_unit       => init
    ! initialize the model
    character(len=*), parameter:: NameSub = 'initialize'
    !--------------------------------------------------------------------------
    if(iProc==0)then
       write(*,'(a)')'SP: '
       write(*,'(a)')'SP: initialize'
       write(*,'(a)')'SP: '
    end if
    call init_unit
    call init_dist
    call init_plot
    call init_mhdata
    call init_turbulence
    call init_spread
    if(DoRestart) call read_restart
    if(IsStandAlone) call init_time

  end subroutine initialize
  !============================================================================
  subroutine finalize

    use SP_ModPlot,       ONLY: finalize_plot       =>finalize
    use SP_ModReadMhData, ONLY: finalize_read       =>finalize
    use SP_ModTurbulence, ONLY: finalize_turbulence => finalize
    ! finalize the model
    character(len=*), parameter:: NameSub = 'finalize'
    !--------------------------------------------------------------------------
    call finalize_plot
    call finalize_read
    call finalize_turbulence

  end subroutine finalize
  !============================================================================
  subroutine run(TimeLimit)

    use SP_ModAdvance,       ONLY: DoTraceShock, advance
    use SP_ModAngularSpread, ONLY: get_magnetic_flux, IsReadySpreadPoint
    use SP_ModGrid,          ONLY: get_other_state_var, copy_old_state,  &
       Rho_, nLine, nVertex_B,  DLogRho_, RhoOld_, State_VIB, MHData_VIB
    use SP_ModReadMhData,    ONLY: read_mh_data
    use SP_ModRestart,       ONLY: stand_alone_save_restart
    use SP_ModPlot,          ONLY: save_plot_all
    use SP_ModTime,          ONLY: SPTime, DataInputTime, iIter
    ! advance the solution in time
    real, intent(in)   :: TimeLimit
    logical, save:: IsFirstCall = .true.
    real:: Dt ! time increment in the current call

    ! write the initial background state to the output file
    !--------------------------------------------------------------------------
    if(IsFirstCall)then
       ! recompute the derived components of state vector, e.g.
       ! magnitude of magnetic field and velocity etc. Smooth if needed.
       if(.not.DoReadMhData)call get_other_state_var
       ! print the initial state
       call save_plot_all(IsInitialOutputIn = .true.)
       ! compute magnetic fluxes associated with lines if needed
       if(IsReadySpreadPoint) call get_magnetic_flux

       IsFirstCall = .false.
    end if
    ! May need to read background data from files
    if(DoReadMhData)then
       ! copy some variables from the previous time step
       call copy_old_state
       ! Read the background data from file
       call read_mh_data()
       ! Read from file: MHData_VIB(0:nMHData,::) for the time moment
       ! DataInputTime
    end if
    ! recompute the derived components of state vector, e.g.
    ! magnitude of magnetic field and velocity etc. Smooth if needed.
    call get_other_state_var
    ! if no new background data loaded, don't advance in time
    if(DataInputTime <= SPTime) RETURN

    call lagr_time_derivative
    ! run the model
    if(DoRun) call advance(min(DataInputTime,TimeLimit))

    ! update time & iteration counters
    iIter = iIter + 1
    Dt = min(DataInputTime,TimeLimit) - SPTime
    SPTime = SPTime + Dt
    call save_plot_all

    ! save restart in the stand alone mod
    if(IsStandAlone)call stand_alone_save_restart(Dt)

  contains
    !==========================================================================
    subroutine lagr_time_derivative

      use SP_ModGrid, ONLY: Used_B

      integer:: iLine, iVertex
      !------------------------------------------------------------------------
      do iLine = 1, nLine
         if(.not.Used_B(iLine))CYCLE
         do iVertex = 1, nVertex_B(  iLine)
            ! divergence of plasma velocity
            State_VIB(DLogRho_,iVertex,iLine) = log(&
                 MHData_VIB(Rho_,iVertex,iLine)/&
                 State_VIB(RhoOld_,iVertex,iLine))
         end do
         ! location of shock
         call get_shock_location(iLine)
      end do

    end subroutine lagr_time_derivative
    !==========================================================================
    subroutine get_shock_location(iLine)

      use SP_ModAdvance, ONLY: nWidth
      use SP_ModGrid,    ONLY: R_, NoShock_, Shock_, ShockOld_, iShock_IB

      integer, intent(in) :: iLine
      ! find location of a shock wave on a given line (line)
      ! Do not search too close to the Sun
      integer         :: iShockMin
      real, parameter :: RShockMin = 1.20  ! *RSun
      ! Threshold for shock tracing
      real, parameter :: DLogRhoThreshold = 0.01
      ! Do not search too close to the heliosphere boundary
      integer:: iShockMax
      ! Misc
      integer:: iShockCandidate
      !------------------------------------------------------------------------
      if(.not.DoTraceShock)then
         iShock_IB(Shock_, iLine) = NoShock_
         RETURN
      end if
      ! shock front is assumed to be location of max gradient log(Rho1/Rho2);
      ! shock never moves back
      iShockMin = max(iShock_IB(ShockOld_, iLine), 1 + nWidth )
      iShockMax = nVertex_B(iLine) - nWidth - 1
      iShockCandidate = iShockMin - 1 + maxloc(&
           State_VIB(DLogRho_,iShockMin:iShockMax,iLine),&
           1, MASK = &
           State_VIB(R_,      iShockMin:iShockMax,iLine) > RShockMin .and. &
           State_VIB(DLogRho_,iShockMin:iShockMax,iLine) > DLogRhoThreshold)
      if(iShockCandidate >= iShockMin)&
           iShock_IB(Shock_, iLine) = iShockCandidate

    end subroutine get_shock_location
    !==========================================================================
  end subroutine run
  !============================================================================
  subroutine check

    use ModUtilities, ONLY: make_dir
    use SP_ModPlot,   ONLY: NamePlotDir
    use SP_ModTiming, ONLY: check_timing => check
    ! Make output directory
    character(len=*), parameter:: NameSub = 'check'
    !--------------------------------------------------------------------------
    if(iProc==0) then
       call make_dir(NamePlotDir)
       ! Initialize timing
       call check_timing
    end if

  end subroutine check
  !============================================================================
end module SP_ModMain
!==============================================================================
