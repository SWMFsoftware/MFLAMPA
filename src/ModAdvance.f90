!  Copyright (C) 2002 Regents of the University of Michigan,
!  portions used with permission
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module SP_ModAdvance

  ! The module contains methods for advancing the solution in time
  use SP_ModSize,   ONLY: nVertexMax
  use SP_ModGrid,   ONLY: State_VIB, MHData_VIB, iShock_IB, D_,   &
       Used_B, Shock_, ShockOld_, nLine, nVertex_B, nWidth
  use SP_ModUnit,   ONLY: UnitX_, Io2Si_V
  use ModUtilities, ONLY: CON_stop

  implicit none

  SAVE

  PRIVATE ! except

  ! Public members:
  public:: read_param    ! Read parameters
  public:: advance       ! Advance Distribution_CB through the time interval
  public:: iterate_steady_state ! Iterate to get steady-state Distribution_CB

  ! If the shock wave is traced, the advance algorithms are modified
  logical, public :: DoTraceShock = .true.

  ! Local variables
  real    :: Cfl = 0.9   ! Controls the maximum allowed time step
  logical :: UsePoissonBracket = .false.
contains
  !============================================================================
  subroutine read_param(NameCommand)

    use ModReadParam, ONLY: read_var
    character(len=*), intent(in):: NameCommand ! From PARAM.in
    character(len=*), parameter:: NameSub = 'read_param'
    !--------------------------------------------------------------------------
    select case(NameCommand)
    case('#CFL')
       call read_var('Cfl', Cfl)
    case('#POISSONBRACKET')
       call read_var('UsePoissonBracket', UsePoissonBracket)
    case('#TRACESHOCK')
       call read_var('DoTraceShock', DoTraceShock)
    case default
       call CON_stop(NameSub//': Unknown command '//NameCommand)
    end select

  end subroutine read_param
  !============================================================================
  subroutine advance(TimeLimit)

    ! advance the solution of the diffusive kinetic equation:
    !    if IsMuAvg: Omnidirectional VDF (Parker transport equation):
    !      f_t + [(1/3)*(d(ln rho)/dt]*f_{ln p} = B*d/ds[D/B*df/ds]
    !    if not IsMuAvg: VDF with pitch angle (Focused transport equation):
    !      f_t + {f; p^3/3*(\vec{u}*\vec{B}/|B|)}_{x, p^3/3}
    !          + {f; (\mu^2-1)*p/(2|B|)}_{x, mu}
    !          + {f; (1-mu^2)/2*(mu*(p^3/3)*
    !             (3\vec{b}\vec{b}:\nabla\vec{u} - \nabla\cdot\vec{u})
    !             + p^2*m_i*bDu/Dt)} = I^(s) = B*d/ds[D/B*df/ds]
    ! with accounting for diffusion and Fermi acceleration
    ! from SPTime to TimeLimit
    ! Prototype: FLAMPA/src/SP_main, case("RUN"), Roussev&Sokolov2008
    ! Version: Borovikov&Sokolov, Dec.19 2017, distinctions:
    ! (1) no turbulence (2) new shock finder moved to SP_ModMain,
    ! and (3) new steepen_shock

    use SP_ModTime,             ONLY: SPTime
    use SP_ModGrid,             ONLY: Rho_, RhoOld_, B_, BOld_, U_
    use SP_ModAdvanceAdvection, ONLY: advect_via_log
    use SP_ModAdvancePoisson,   ONLY: advect_via_poisson, &
         init_data_states, advect_via_multi_poisson
    use SP_ModDiffusion,        ONLY: UseDiffusion, set_diffusion_coef
    use SP_ModDistribution,     ONLY: IsMuAvg, IsDistNeg

    real, intent(in):: TimeLimit
    ! Loop variable
    integer :: iLine
    ! For a given line: nVertex_B, iShock_IB:
    integer :: iEnd, iShock, iShockOld
    ! Upper limit and variable for the Loop which makes
    ! a time step so short that the shock wave passes
    ! a single grid interval per a progress step:
    integer :: iProgress, nProgress
    ! coefficient to interpolate "old" and "new"
    real    :: Alpha
    ! Full difference between DataInputTime and SPTime
    real    :: DtFull
    ! Time step in the PROGRESS Loop, DtFull/nProgress
    real    :: DtProgress
    ! Current time: for updating states in multi-Poisson-bracket scheme
    real    :: Time
    ! Local arrays to store the state vectors in SI units
    real, dimension(1:nVertexMax):: nOldSi_I, nSi_I, BSi_I
    ! Lagrangian derivatives
    real, dimension(1:nVertexMax):: dLogRho_I

    !--------------------------------------------------------------------------
    DtFull = TimeLimit - SPTime
    ! go line by line and advance the solution

    LINE:do iLine = 1, nLine
       if(.not.Used_B(iLine)) CYCLE line
       ! the active particles on the line
       iEnd = nVertex_B(iLine)

       ! Various data along the line in SI units.
       ! The IO units of the state vectors could be seen in ModUnit, the
       ! summary is: Length is in the unit of Rs, temperature is in the
       ! unit of kinetic energy, all others are in SI units.
       nOldSi_I(1:iEnd) = State_VIB(RhoOld_, 1:iEnd, iLine)

       ! find how far shock has travelled on this line: nProgress
       iShock    = iShock_IB(Shock_,    iLine)
       iShockOld = iShock_IB(ShockOld_, iLine)
       if(DoTraceShock) then
          ! This is how many steps should be done to allow the shock to
          ! the move not more than one mesh size
          nProgress = max(1, iShock-iShockOld)
          iShockOld = min(iShockOld, iShock-1)
       else
          nProgress = 1
          iShockOld = 0
       end if

       ! each particles shock has crossed should be
       ! processed separately => reduce the time step
       DtProgress = DtFull/nProgress

       ! go over each crossed particle
       PROGRESS:do iProgress = 1, nProgress
          ! account for change in the background up to the current moment
          Alpha = real(iProgress)/real(nProgress)

          ! nSi is needed to set up the distribution at the injection.
          ! It is calculated at the end of the iProgress' time step
          nSi_I(1:iEnd) = State_VIB(RhoOld_, 1:iEnd, iLine) + Alpha* &
               (MhData_VIB( Rho_, 1:iEnd, iLine) - &
               State_VIB(RhoOld_, 1:iEnd, iLine))
          BSi_I(1:iEnd) = State_VIB(  BOld_, 1:iEnd, iLine) + Alpha* &
               (State_VIB(  B_, 1:iEnd, iLine) - &
               State_VIB(BOld_, 1:iEnd, iLine))
          dLogRho_I(1:iEnd) = log(nSi_I(1:iEnd)/nOldSi_I(1:iEnd))

          ! trace shock position and steepen the shock
          iShock = iShockOld + iProgress
          if(iShock < iEnd-nWidth .and. iShock > nWidth  &
               .and. DoTraceShock) call steepen_shock(iEnd)

          ! Advection (2 different schemes) and Diffusion
          ! First, set the diffusion coefficient, from the
          ! given formulae or from the turbulent specrtum, if known
          if(UseDiffusion) call set_diffusion_coef(iLine, iEnd,   &
               iShock, BSi_I(1:iEnd))
          if(UsePoissonBracket) then
             ! Poisson bracket scheme: particle-number-conservative
             if(IsMuAvg) then
                ! Single Poisson bracket: Parker transport equation
                call advect_via_poisson(iLine, iEnd, iShock, DtProgress,&
                     Cfl, nOldSi_I(1:iEnd), nSi_I(1:iEnd), BSi_I(1:iEnd))
             else
                ! Multiple Poisson brackets: Focused transport equation
                ! See the descriptions for development in ModAdvancePoisson.f90
                Time = Alpha*DtFull - DtProgress            ! Tstart this step
                call init_data_states(iLine, iEnd, DtFull)  ! Inital states
                call advect_via_multi_poisson(iLine, iEnd, iShock, &
                     Time, DtProgress, Cfl, nSi_I(1:iEnd), BSi_I(1:iEnd))
             end if
          else
             ! No Poisson bracket scheme, use the default algorithm
             call advect_via_log(iLine, iEnd, iShock, DtProgress, Cfl,  &
                  dLogRho_I(1:iEnd), nSi_I(1:iEnd), BSi_I(1:iEnd))
             if(IsDistNeg) CYCLE LINE
          end if
          ! Store the old density at the end of each iProgress
          nOldSi_I(1:iEnd) = nSi_I(1:iEnd)
       end do PROGRESS
    end do LINE

  contains
    !==========================================================================
    subroutine steepen_shock(iEnd)

      ! change the density profile near the shock front so it
      ! becomes steeper for the current line
      use SP_ModGrid, ONLY: dLogRhoThreshold
      integer, intent(in) :: iEnd
      real   :: DsSi_I(1:iEnd-1)
      real   :: dLogRhoExcess_I(iShock-nWidth:iShock+nWidth-1)
      real   :: dLogRhoExcessIntegral
      ! find the excess of dLogRho within the shock compared
      ! to background averaged over length
      !------------------------------------------------------------------------
      DsSi_I = State_VIB(D_,1:iEnd-1,iLine)*Io2Si_V(UnitX_)
      dLogRhoExcess_I = max(0.5*(dLogRho_I(iShock-nWidth:iShock+nWidth-1) + &
           dLogRho_I(iShock-nWidth+1:iShock+nWidth)) - dLogRhoThreshold, 0.0)
      ! A jump (dLogRhoExcess>0) in velocity accross the shock wave * \Delta t
      dLogRhoExcessIntegral = sum(dLogRhoExcess_I*&
           DsSi_I(iShock-nWidth:iShock+nWidth-1))

      ! check for zero excess
      if(dLogRhoExcessIntegral == 0.0) RETURN
      ! nullify excess within the smoothed shock
      dLogRho_I(iShock-nWidth:iShock+nWidth) = min(dLogRhoThreshold, &
           dLogRho_I(iShock-nWidth:iShock+nWidth))
      ! ... and concentrate it at the shock front, applying the whole jump
      ! in the velocity at a single grid point
      dLogRho_I(iShock) = dLogRhoThreshold + &
           dLogRhoExcessIntegral/DsSi_I(iShock)
      ! also, sharpen the magnetic field magnitude
      ! post shock part
      BSi_I(iShock+1-nWidth:iShock+1)=maxval(BSi_I(iShock+1-nWidth:iShock+1))
      ! pre shock part
      BSi_I(iShock+1:iShock+nWidth  )=minval(BSi_I(iShock+1:iShock+nWidth))

    end subroutine steepen_shock
    !==========================================================================
  end subroutine advance
  !============================================================================
  subroutine iterate_steady_state

    ! advance the solution of the diffusive kinetic equation:
    !     f_t+[(1/3)*(d(ln rho)/dt]*f_{ln p}=B*d/ds[D/B*df/ds]
    ! with accounting for diffusion and Fermi acceleration

    use SP_ModGrid,           ONLY: Rho_, B_
    use SP_ModAdvancePoisson, ONLY: iterate_poisson
    use SP_ModDiffusion,      ONLY: UseDiffusion, set_diffusion_coef

    ! Loop variable
    integer :: iLine
    ! For a given line: nVertex_B, iShock_IB:
    integer :: iEnd, iShock
    ! Local arrays to store the state vectors in SI units
    real, dimension(1:nVertexMax):: nSi_I, BSi_I
    ! go line by line and iterate the solution
    !--------------------------------------------------------------------------

    do iLine = 1, nLine
       if(.not.Used_B(iLine)) CYCLE
       ! the active particles on the line
       iEnd = nVertex_B(iLine)

       ! Various data along the line in SI units. Temperature is in the unit
       ! of kinetic energy, all others are in SI units.
       BSi_I(1:iEnd) = State_VIB(   B_, 1:iEnd, iLine)
       ! nSi is needed to set up the distribution at the injection.
       nSI_I(1:iEnd) = MhData_VIB(Rho_, 1:iEnd, iLine)
       ! find how far shock has travelled on this line
       iShock = iShock_IB(Shock_, iLine)

       ! First, set the diffusion coefficient, from the
       ! given formulae or from the turbulent specrtum, if known
       if(UseDiffusion) call set_diffusion_coef(iLine, iEnd,   &
            iShock, BSi_I(1:iEnd))
       ! Poisson bracket scheme: particle-number-conservative
       call iterate_poisson(iLine, iEnd, iShock, Cfl, BSi_I(1:iEnd), &
            nSi_I(1:iEnd))
    end do

  end subroutine iterate_steady_state
  !============================================================================
end module SP_ModAdvance
!==============================================================================
