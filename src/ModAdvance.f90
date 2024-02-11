!  Copyright (C) 2002 Regents of the University of Michigan,
!  portions used with permission
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module SP_ModAdvance

  ! The module contains methods for advancing the solution in time
  use SP_ModSize, ONLY: nVertexMax
  use SP_ModGrid, ONLY: State_VIB, MHData_VIB, iShock_IB, D_,  &
       Used_B, Shock_, ShockOld_, nLine, nVertex_B, nWidth
  use SP_ModUnit, ONLY: UnitX_, Io2Si_V
  use ModUtilities, ONLY: CON_stop

  implicit none

  SAVE

  PRIVATE ! except

  ! Public members:
  public:: read_param  ! read injection parameters
  public:: advance     ! Advance solution Distribution_IIB

  ! If the shock wave is traced, the advance algorithms are modified
  logical, public :: DoTraceShock = .true.

  ! Local parameters
  real:: Cfl=0.9        ! Controls the maximum allowed time step
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
       call read_var('Cfl',Cfl)
    case('#POISSONBRACKET')
       call read_var('UsePoissonBracket',UsePoissonBracket)
    case('#TRACESHOCK')
       call read_var('DoTraceShock', DoTraceShock)
    case default
       call CON_stop(NameSub//': Unknown command '//NameCommand)
    end select
  end subroutine read_param
  !============================================================================
  subroutine advance(TimeLimit)
    ! advance the solution of the diffusive kinetic equation:
    !            f_t+[(1/3)*(d(ln rho)/dt]*f_{ln p}=B*d/ds[D/B*df/ds]
    ! with accounting for diffusion and Fermi acceleration
    ! from SPTime to TimeLimit
    ! Prototype: FLAMPA/src/SP_main, case("RUN"), Roussev&Sokolov2008
    ! Version: Borovikov&Sokolov, Dec.19 2017, distinctions:
    ! (1) no turbulence (2) new shock finder moved to SP_ModMain,
    ! and (3) new steepen_shock
    use SP_ModTime,             ONLY: SPTime
    use SP_ModGrid,             ONLY: Rho_, RhoOld_, B_, BOld_, U_
    use SP_ModAdvanceAdvection, ONLY: advect_via_log
    use SP_ModAdvancePoisson,   ONLY: advect_via_poisson
    use SP_ModDiffusion,        ONLY: UseDiffusion, set_diffusion_coef
    real, intent(in):: TimeLimit
    ! Loop variables
    integer  :: iP, iVertex, iLine
    ! For a given line: nVertex_B, iShock_IB:
    integer  :: iEnd, iShock, iShockOld
    ! Upper limit and variable for the Loop which makes
    ! a time step so short that the shock wave passes
    ! a single grid interval per a progress step:
    integer  :: iProgress, nProgress
    ! coefficient to interpolate "old" and "new"
    real     :: Alpha
    ! Full difference between DataInputTime and SPTime
    real     :: DtFull
    ! Time step in the PROGRESS Loop, DtFull/nProgress
    real     :: DtProgress
    ! Local arrays to store the state vectors in SI units
    real, dimension(1:nVertexMax):: nSi_I, BSi_I, BOldSi_I, nOldSi_I, uSi_I
    ! Lagrangian derivatives
    real, dimension(1:nVertexMax):: DLogRho_I

    ! Check if any Distribution_IIB < 0 in advect_via_log
    logical  :: IsNeg
    character(len=*), parameter:: NameSub = 'advance'
    !--------------------------------------------------------------------------
    DtFull = TimeLimit - SPTime
    ! go line by line and advance the solution

    line:do iLine = 1, nLine
       if(.not.Used_B(iLine))CYCLE line
       ! the active particles on the line
       iEnd   = nVertex_B( iLine)

       ! Various data along the line in SI units.
       ! The IO units of the state vectors could be seen in ModUnit, the
       ! summary is: Length is in the unit of Rs, Rho is in the unit of
       ! amu/m^3, temperature is in the unit of kinetic energy, all others
       ! are in SI units. So to convert the IO units, the three conversion
       ! factors are needed as follows:
       uSi_I(   1:iEnd) = State_VIB(U_,     1:iEnd,iLine)
       BOldSi_I(1:iEnd) = State_VIB(BOld_,  1:iEnd,iLine)
       nOldSi_I(1:iEnd) = State_VIB(RhoOld_,1:iEnd,iLine)

       ! find how far shock has travelled on this line: nProgress
       iShock    = iShock_IB(Shock_,   iLine)
       iShockOld = iShock_IB(ShockOld_,iLine)
       if(DoTraceShock)then
          ! This is how many steps should be done to allow the shock to
          ! the move not more than one mesh size
          nProgress = MAX(1, iShock - iShockOld)
          iShockOld = MIN(iShockOld, iShock-1)
       else
          nProgress = 1
          iShockOld = 0
       end if

       ! each particles shock has crossed should be
       ! processed separately => reduce the time step
       DtProgress = DtFull / nProgress

       ! go over each crossed particle
       PROGRESS:do iProgress = 1, nProgress
          ! account for change in the background up to the current moment
          Alpha = real(iProgress) / real(nProgress)

          ! Recall that MhData_VIB(Rho_) and State_VIB(RhoOld_) are in
          ! the unit of amu/m^3.
          ! nSi is needed to set up the distribution at the injection.
          ! It is calculated at the end of the iProgress' time step
          nSI_I(1:iEnd) = State_VIB(RhoOld_,1:iEnd,iLine) + Alpha* &
               (MhData_VIB(Rho_,  1:iEnd,iLine) - &
               State_VIB(RhoOld_,1:iEnd,iLine))
          BSi_I(1:iEnd) = State_VIB(BOld_,1:iEnd,iLine) + Alpha* &
               (State_VIB(B_,  1:iEnd,iLine) - &
               State_VIB(BOld_,1:iEnd,iLine))
          DLogRho_I(1:iEnd) = log(nSi_I(1:iEnd)/nOldSi_I(1:iEnd))

          ! trace shock position and steepen the shock
          iShock = iShockOld + iProgress
          if(iShock < iEnd-nWidth .and. iShock > nWidth     &
               .and. DoTraceShock) call steepen_shock(iEnd)

          ! Advection (2 different schemes) and Diffusion
          ! First, set the diffusion coefficient, from the
          ! given formulae or from the turbulent specrtum, if known
          if(UseDiffusion) call set_diffusion_coef(iLine, iEnd,   &
               iShock, BSi_I(1:iEnd))
          if(UsePoissonBracket)then
             ! Poisson bracket scheme: particle-number-conservative
             call advect_via_poisson(iLine, iEnd, iShock, DtProgress,   &
                  Cfl, nOldSi_I(1:iEnd), nSi_I(1:iEnd), BSi_I(1:iEnd))
             ! store density and B-field arrays at the end of this time step
             nOldSi_I(1:iEnd) = nSi_I(1:iEnd)
             BOldSi_I(1:iEnd) = BSi_I(1:iEnd)
          else
             ! No Poisson bracket scheme, use the default algorithm
             call advect_via_log(iLine, iEnd, iShock, DtProgress, Cfl,  &
                  DLogRho_I(1:iEnd), nSi_I(1:iEnd), BSi_I(1:iEnd), IsNeg)
             if(IsNeg) CYCLE line
          end if
       end do PROGRESS
    end do line
  contains
    !==========================================================================
    subroutine steepen_shock(iEnd)
      use SP_ModGrid, ONLY: dLogRhoThreshold
      integer, intent(in) :: iEnd ! To limit the range of search
      ! change the density profile near the shock front so it becomes steeper
      ! for the current line
      real   :: DsSi_I(1:iEnd)
      integer:: iVertex ! loop variable
      real   :: DLogRhoExcessIntegral, DLogRhoExcess
      ! find the excess of DLogRho within the shock compared to background
      ! averaged over length
      !------------------------------------------------------------------------
      DLogRhoExcessIntegral = 0.0
      DsSi_I = State_VIB(D_,1:iEnd,iLine)*Io2Si_V(UnitX_)
      do iVertex = iShock - nWidth, iShock + nWidth - 1
         DLogRhoExcess = 0.5*(DLogRho_I(iVertex) + DLogRho_I(iVertex+1)) &
              - DLogRhoThreshold ! D log(rho)/Dt*\Delta t = -\div U*\Delta t
         if(DLogRhoExcess>0) then
            ! This is a jump in velocity accross the shock wave * \Delta t
            DLogRhoExcessIntegral = DLogRhoExcessIntegral + &
                 DLogRhoExcess*DsSi_I(iVertex)
         end if
      end do

      ! check for zero excess
      if(DLogRhoExcessIntegral == 0.0) RETURN
      ! nullify excess  within the smoothed shock
      DLogRho_I(iShock-nWidth:iShock+nWidth) = min(DLogRhoThreshold, &
           DLogRho_I(iShock-nWidth:iShock+nWidth))
      ! ...and concetrate it at the shock front, applying the whole jump
      ! in the velocity at a single grid point
      DLogRho_I(iShock) = DLogRhoThreshold + DLogRhoExcessIntegral / &
           DsSi_I(iVertex)
      ! also, sharpen the magnetic field magnitude
      ! post shock part
      BSi_I(iShock+1-nWidth:iShock+1)=maxval(BSi_I(iShock+1-nWidth:iShock+1))
      ! pre shock part
      BSi_I(iShock+1:iShock+nWidth  )=minval(BSi_I(iShock+1:iShock+nWidth))
    end subroutine steepen_shock
    !==========================================================================
  end subroutine advance
  !============================================================================
end module SP_ModAdvance
!==============================================================================
