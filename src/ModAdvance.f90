!  Copyright (C) 2002 Regents of the University of Michigan,
!  portions used with permission
!  For more information, see http://csem.engin.umich.edu/tools/swmf

module SP_ModAdvance

  ! The module contains methods for advancing the solution in time

  use ModNumConst, ONLY: cPi
  use ModConst,   ONLY: cMu
  use SP_ModSize, ONLY: nVertexMax
  use SP_ModDistribution, ONLY: nP, Distribution_IIB,                &
       MomentumSI_I, MomentumInjSI, DLogP
  use SP_ModGrid, ONLY: State_VIB, MHData_VIB, iShock_IB,            &
       R_, x_, y_, z_, Used_B, Shock_, NoShock_,                     &
       ShockOld_, DLogRho_, nLine, nVertex_B
  !  use SP_ModTurbulence, ONLY: DoInitSpectrum, UseTurbulentSpectrum,  &
    !   set_wave_advection_rates, reduce_advection_rates, init_spectrum
  use SP_ModUnit, ONLY: UnitX_, UnitEnergy_,                         &
       Io2Si_V, kinetic_energy_to_momentum
  use ModUtilities, ONLY: CON_stop

  implicit none

  SAVE

  PRIVATE ! except

  ! Public members:
  public:: read_param  ! read injection parameters
  public:: advance     ! Advance solution Distribution_IIB
  !||||||||||||||Boundary condition at the injection energy!!!!!!
  ! Injection efficiency and assumed spectral index with the energy
  ! range k_BT_i< Energy < EnergyInjection, to be read from PARAM.in
  real:: CoefInj = 0.25, SpectralIndex = 5.0

!!!!!!!!!!!!!!!!!!!!!!!!! Local parameters!!!!!!!!!!!!!!!
  real:: Cfl=0.9        ! Controls the maximum allowed time step
  integer, public, parameter :: nWidth = 50
  logical :: UsePoissonBracket = .false.
  logical, public:: DoTraceShock = .true., UseDiffusion = .true.

contains
  !============================================================================
  subroutine read_param(NameCommand)
    use ModReadParam, ONLY: read_var
    character(len=*), intent(in):: NameCommand ! From PARAM.in
    character(len=8) :: StringScaleTurbulenceType
    character(len=*), parameter:: NameSub = 'read_param'
    !--------------------------------------------------------------------------
    select case(NameCommand)
    case('#INJECTION')
       call read_var('Efficiency',   CoefInj)
       call read_var('SpectralIndex',SpectralIndex)
    case('#CFL')
       call read_var('Cfl',Cfl)
    case('#POISSONBRACKET')
       call read_var('UsePoissonBracket',UsePoissonBracket)
    case('#TRACESHOCK')
       call read_var('DoTraceShock', DoTraceShock)
    case('#DIFFUSION')
       call read_var('UseDiffusion', UseDiffusion)
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
    use ModConst,             ONLY: cProtonMass, Rsun
    use SP_ModTime,           ONLY: SPTime
    use SP_ModGrid,           ONLY: D_, Rho_, RhoOld_,   &
         B_, BOld_, U_, T_
    use SP_ModAdvection,      ONLY: advance_log_advection
    use SP_ModAdvancePoisson, ONLY: advect_via_poisson_bracket
    use SP_ModDiffusion,      ONLY: diffuse_distribution
    real, intent(in):: TimeLimit
    ! Loop variables
    integer  :: iP, iVertex, iLine
    ! For a given line: nVertex_B, iShock_IB:
    integer  :: iEnd, iShock, iShockOld
    ! Upper limit and variable for the Loop which makes
    ! a time step so short that the shock wave passes
    ! a single grid interval per a progress step:
    integer  :: iProgress, nProgress
    ! Upper limit and variable for the loop which makes a
    ! time step so short that the CFL for the (explicit)
    ! advection is less that CFL declared abobe.
    integer  :: iStep, nStep
    ! coefficient to interpolate  "old" and "new"
    real     :: Alpha
    ! Lower limit to floor the spatial diffusion coefficient For a
    ! given spatial and temporal resolution, the value of the
    ! diffusion coefficient should be artificially increased to get
    ! the diffusion length to be larger than the shock front width,
    ! which even for the steepened shock is as wide as a mesh size
    ! of the Largangian grid, State_VIB(D_,:,:)*RSun. In this way,
    ! the Lagrangian grid resolution is sufficient to resolve a
    ! precursor in the upstream distribution of
    ! DiffCoeffMin=0 (default value), we do NOT enhance
    ! the diffusion coefficient! Physically, DiffCoeffMin
    ! should be given by the product of shock wave speed
    ! and local grid spacing.
    real      :: DtFull
    ! Time step in the PROGRESS Loop, DtFull/nProgress
    real      :: DtProgress
    ! Time step in the STEP Loop, DtProgress/nStep
    real      :: Dt
    real      :: MachAlfven
    real      :: DtReduction

    ! Local arrays to store the position and state vectors in SI units
    real :: XyzSI_DI(3, 1:nVertexMax)
    real, dimension(1:nVertexMax):: RadiusSi_I, DsSI_I,    &
         nSI_I, uSI_I, BSI_I, BOldSI_I, nOldSi_I
    real, dimension(1:nVertexMax):: InvRhoOld_I, InvRho_I
    real, dimension(1:nVertexMax) :: DOuterSI_I, CoefDInnerSI_I

    ! Lagrangian derivatives
    real, dimension(1:nVertexMax):: DLogRho_I, FermiFirst_I

    ! the full time interval
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
       ! UnitX_, UnitRho_, UnitEnergy_
       XyzSI_DI(x_:z_,1:iEnd) = MhData_VIB(x_:z_,1:iEnd,iLine)*IO2SI_V(UnitX_)
       DsSI_I(     1:iEnd) = State_VIB(D_,     1:iEnd,iLine)*IO2SI_V(UnitX_)
       RadiusSI_I( 1:iEnd) = State_VIB(R_,     1:iEnd,iLine)*IO2SI_V(UnitX_)
       uSI_I(      1:iEnd) = State_VIB(U_,     1:iEnd,iLine)
       BOldSI_I(   1:iEnd) = State_VIB(BOld_,  1:iEnd,iLine)
       nOldSI_I(   1:iEnd) = State_VIB(RhoOld_,1:iEnd,iLine)

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
          !
          ! nSI is needed to set up the distribution at the injection.
          ! It is calculated at the end of the iProgress' time step
          nSi_I(1:iEnd) = State_VIB(RhoOld_,1:iEnd,iLine) +Alpha *&
               (MhData_VIB(Rho_,  1:iEnd,iLine) - &
               State_VIB(RhoOld_,1:iEnd,iLine))
          ! dimensionless
          DLogRho_I(1:iEnd)=log(nSi_I(1:iEnd)/nOldSI_I(1:iEnd))

          BSI_I(1:iEnd) = State_VIB(BOld_,1:iEnd,iLine) + Alpha*&
               (State_VIB(B_,  1:iEnd,iLine) - &
               State_VIB(BOld_,1:iEnd,iLine))

          iShock = iShockOld + iProgress

          ! find the shock alfven mach number, also steepen the shock
          if(iShock < iEnd - nWidth.and.iShock > nWidth.and.DoTraceShock)then
             MachAlfven = mach_alfven()
             call steepen_shock
          else
             MachAlfven = 1.0
          end if

          ! 1st order Fermi acceleration is responsible for advection
          ! in momentum space
          ! first order Fermi acceleration for the current line, dimensionless
          FermiFirst_I(1:iEnd) = DLogRho_I(1:iEnd) / (3*DLogP)

          ! if(UseTurbulentSpectrum)then
             ! Calculate the Alfven speed
             ! if (DoInitSpectrum) call init_spectrum(iEnd,              &
             !     XyzSI_DI(x_:z_, 1:iEnd),BSI_I(1:iEnd), MomentumSI_I, &
             !     dLogP, iShock, CoefInj, MachAlfven)
             ! call set_wave_advection_rates(iEnd,     &
             ! BSI_I(1:iEnd),    BOldSI_I(1:iEnd),      &
             ! nSi_I(1:iEnd)*cProtonMass,  nOldSI_I(1:iEnd)*cProtonMass, &
             ! XyzSI_DI(x_:z_, 1:iEnd), DsSI_I(1:iEnd), &
             ! DLogP, DtProgress, DtReduction)

             ! nStep = 1+int(max(DtReduction,                &
             !      maxval(abs(FermiFirst_I(1:iEnd))))/CFL)
             ! nStep = 1+int(maxval(abs(FermiFirst_I(2:iEnd)))/CFL)
          ! else
          ! How many steps should be done to the CFL criterion is fulfilled
          nStep = 1+int(maxval(abs(FermiFirst_I(2:iEnd)))/CFL)
          ! end if

          ! Check if the number of time steps is positive:
          if(nStep < 1)then
             ! if(UseTurbulentSpectrum) &
             !     write(*,*) ' DtReduction               =', DtReduction
             write(*,*) ' maxval(abs(FermiFirst_I)) =', &
                  maxval(abs(FermiFirst_I(2:iEnd)))
             call CON_stop(NameSub//': nStep <= 0????')
          end if
          ! Store the value at the end of the previous time step

          ! set the left boundary condition (for diffusion)
          Distribution_IIB(1:nP+1, 1, iLine) = &
               Distribution_IIB(0, 1, iLine) * &
               (MomentumSI_I(0)/MomentumSI_I(1:nP+1))**SpectralIndex

          ! Advection (with 2 different Algorithms) & Diffusion
          if(UsePoissonBracket)then
             ! update bc for advection
             call set_advection_bc
             ! store/update the inverse rho arrays
             InvRhoOld_I(1:iEnd) = 1.0/nOldSI_I(1:iEnd)
             InvRho_I(1:iEnd)    = 1.0/nSi_I(1:iEnd)
             call advect_via_poisson_bracket(nP, iEnd, DtProgress, &
                  Cfl, InvRhoOld_I(1:iEnd), InvRho_I(1:iEnd),        &
                  Distribution_IIB(:,1:iEnd,iLine), iLine, iShock,   &
                  XyzSI_DI, nSI_I, BSI_I, DsSI_I, RadiusSI_I, UseDiffusion)
             nOldSI_I(1:iEnd) = nSI_I(1:iEnd)
             BOldSI_I(1:iEnd) = BSI_I(1:iEnd)
             ! DoInitSpectrum = .true.
          else
             ! No Poisson bracket, use the default algorithm
             Dt = DtProgress/nStep

             FermiFirst_I(1:iEnd) = FermiFirst_I(1:iEnd) / nStep
             ! if(UseTurbulentSpectrum) call reduce_advection_rates(nStep)

             ! compute diffusion along the field line
             ! we calculate: "Outer diffusion"=BSI_I and
             ! "Inner diffusion at the injection energy (iP=0)
             STEP:do iStep = 1, nStep
                ! update bc for advection
                call set_advection_bc
                ! advection in the momentum space
                do iVertex = 2, iEnd
                   if(any(Distribution_IIB(0:nP+1,iVertex,iLine) < 0.0)) then
                      write(*,*) NameSub, ': Distribution_IIB < 0'
                      Used_B(iLine) = .false.
                      nVertex_B(iLine) = 0
                      CYCLE line
                   end if

                   call advance_log_advection(FermiFirst_I(iVertex), nP,   &
                        1, 1, Distribution_IIB(0:nP+1,iVertex,iLine), .false.)
                end do

                if(UseDiffusion) call diffuse_distribution(iLine, iEnd,    &
                     iShock, Dt, Distribution_IIB(0:nP+1, 1:iEnd, iLine),    &
                     XyzSI_DI, nSI_I, BSI_I, DsSI_I, RadiusSI_I)
             end do STEP
             ! DoInitSpectrum = .true.
          end if
       end do PROGRESS
    end do line
  contains
    !==========================================================================
    function mach_alfven() result(MachAlfven)
      ! alfvenic mach number for the current line
      real:: MachAlfven
      real:: SpeedAlfvenUpstream, SpeedUpstream

      ! speed upstream is relative to the shock:
      ! \rho_u * (U_u - V_{shock}) = \rho_d * (U_d - V_{shock})
      !=> V_{shock}(\rho_d - \rho_u) = \rho_d*U_d -\rho_u*U_u
      !=> V_{shock} - U_u=\rho_d*(U_d - U_u)/(\rho_d - \rho_u)
      ! Hence, below there should be norm2(\vect{U}_d - \vect{U}_u)
      !------------------------------------------------------------------------
      SpeedUpstream = nSi_I(iShock+1-nWidth)*&
           (uSI_I(  iShock + 1 - nWidth) - uSI_I(  iShock + nWidth)) /  &
           (nSi_I(iShock + 1 - nWidth) - nSi_I(iShock + nWidth))
      SpeedAlfvenUpstream = BSI_I(iShock + nWidth)/ &
           sqrt(cMu*nSi_I(iShock + nWidth)*cProtonMass)
      MachAlfven = SpeedUpstream / SpeedAlfvenUpstream
    end function mach_alfven
    !==========================================================================
    subroutine steepen_shock
      ! change the density profile near the shock front so it becomes steeper
      ! for the current line
      integer:: iVertex ! loop variable
      real   :: DLogRhoExcessIntegral, DLogRhoExcess
      real, parameter:: DLogRhoBackground = 0.01
      !------------------------------------------------------------------------
      ! find the excess of DLogRho within the shock compared to background
      ! averaged over length
      DLogRhoExcessIntegral = 0.0
      do iVertex = iShock - nWidth, iShock + nWidth - 1
         DLogRhoExcess = 0.5*(DLogRho_I(iVertex) + DLogRho_I(iVertex+1)) &
              - DLogRhoBackground ! D log(rho)/Dt*\Delta t = -\div U*\Delta t
         if(DLogRhoExcess>0) then
            ! This is a jump in velocity accross the shock wave * \Delta t
            DLogRhoExcessIntegral = DLogRhoExcessIntegral + &
                 DLogRhoExcess*DsSI_I(iVertex)
         end if
      end do

      ! check for zero excess
      if(DLogRhoExcessIntegral == 0.0)RETURN
      ! nullify excess  within the smoothed shock
      DLogRho_I(iShock-nWidth:iShock+nWidth) = min(DLogRhoBackground, &
           DLogRho_I(iShock-nWidth:iShock+nWidth))
      ! ...and concetrate it at the shock front, applying the whole jump
      ! in the velocity at a single grid point
      DLogRho_I(iShock) = DLogRhoBackground + DLogRhoExcessIntegral / &
           DsSI_I(iVertex)
      ! also, sharpen the magnetic field magnitude
      ! post shock part
      BSI_I(iShock+1-nWidth:iShock+1)=maxval(BSI_I(iShock+1-nWidth:iShock+1))
      ! pre shock part
      BSI_I(iShock+1:iShock+nWidth  )=minval(BSI_I(iShock+1:iShock+nWidth))
    end subroutine steepen_shock
    !==========================================================================
    subroutine set_advection_bc
      ! set boundary conditions on grid point on the current line
      ! local variables

      integer:: iVertex     ! loop variable
      real   :: MomentumSI    ! Momentum for the thermal energy k_BTi
      real   :: CoefInjLocal, DistributionBc
      !------------------------------------------------------------------------
      do iVertex = 1, iEnd
         ! injection(Ti, Rho), see Sokolov et al., 2004, eq (3)
         ! f = CoefInj/2/pi * N / (2*m*T_p)^(3/2) * ((2*m*T_p)^(3/2)/p_inj)^5
         !   = CoefInj/2/pi * N / p^3 * (p/p_inj)^5
         ! where p = sqrt(2*m*T_p) is the momentum of thermal ion
         CoefInjLocal = 2.50E-11
         MomentumSI   = kinetic_energy_to_momentum(                     &
              MHData_VIB(T_,iVertex,iLine)*Io2Si_V(UnitEnergy_))

         DistributionBc = (SpectralIndex-3)/(4*cPi)                     &
              * nSI_I(iVertex)/MomentumSI**3                            &
              * (MomentumSI/MomentumInjSI)**SpectralIndex

         if (iShock /= NoShock_           .and.                         &
              iVertex <= iShock + nWidth    .and.                         &
              iVertex >= iShock - nWidth) then
            CoefInjLocal = CoefInj
         endif

         Distribution_IIB(0,iVertex,iLine) = DistributionBc * CoefInjLocal

      end do

    end subroutine set_advection_bc
    !==========================================================================
  end subroutine advance
  !============================================================================
end module SP_ModAdvance
!==============================================================================
