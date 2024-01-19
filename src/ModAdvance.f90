!  Copyright (C) 2002 Regents of the University of Michigan,
!  portions used with permission
!  For more information, see http://csem.engin.umich.edu/tools/swmf

module SP_ModAdvance

  ! The module contains methods for advancing the solution in time

  use ModNumConst, ONLY: cPi
  use ModConst,   ONLY: cLightSpeed, cGEV, cAu, cMu
  use SP_ModSize, ONLY: nVertexMax
  use SP_ModDistribution, ONLY: nP, Distribution_IIB, MomentumSI_I,           &
        EnergySI_I, MomentumMaxSI, MomentumInjSI, DLogP 
  use SP_ModGrid, ONLY: State_VIB, MHData_VIB, iShock_IB,  R_, x_, y_, z_,    &
       iLineAll0, Used_B,     &
       Shock_, NoShock_, ShockOld_, DLogRho_, Wave1_, Wave2_, nLine, nVertex_B
  use SP_ModTurbulence, ONLY: DoInitSpectrum, UseTurbulentSpectrum,           &
       set_wave_advection_rates, reduce_advection_rates, init_spectrum,       &
       update_spectrum
  use SP_ModUnit, ONLY: UnitX_, UnitEnergy_, &
       NameParticle, Io2Si_V, kinetic_energy_to_momentum
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
  ! Diffusion as in Li et al. (2003), doi:10.1029/2002JA009666
  logical, public :: UseFixedMFPUpstream = .false.
  real    :: MeanFreePath0InAu = 1.0
  logical, public:: DoTraceShock = .true., UseDiffusion = .true.

  logical :: UsePoissonBracket = .false.

  ! Parameter characterizing cut-off wavenumber of turbulent spectrum:
  ! value of scale turbulence at 1 AU for any type (const or linear)
  real:: ScaleTurbulenceSI = 0.03 * cAu
  integer:: iScaleTurbulenceType
  integer, parameter:: Const_ = 0, Linear_ = 1

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
       call read_var('SpectralIndex',SpectralIndex)
       call read_var('Efficiency',   CoefInj)
    case('#CFL')
       call read_var('Cfl',Cfl)
    case('#USEFIXEDMFPUPSTREAM')
       call read_var('UseFixedMFPUpstream',UseFixedMFPUpstream)
       if(UseFixedMFPUpstream)then
          ! see Li et al. (2003), doi:10.1029/2002JA009666
          call read_var('MeanFreePath0 [AU]', MeanFreePath0InAu)
       end if
    case('#SCALETURBULENCE')
       ! cut-off wavenumber of turbulence spectrum: k0 = 2 cPi / Scale
       call read_var('ScaleTurbulenceType', StringScaleTurbulenceType)
       select case(StringScaleTurbulenceType)
       case('const', 'constant')
          iScaleTurbulenceType = Const_
       case('linear')
          iScaleTurbulenceType = Linear_
       case default
          call CON_stop(NameSub//": unknown scale turbulence type")
       end select
       call read_var('ScaleTurbulence [AU] at 1 AU', ScaleTurbulenceSI)
       ScaleTurbulenceSI = ScaleTurbulenceSI * cAu
    case('#DIFFUSION')
       call read_var('UseDiffusion', UseDiffusion)
    case('#TRACESHOCK')
       call read_var('DoTraceShock', DoTraceShock)
    case('#POISSONBRACKET')
       call read_var('UsePoissonBracket',UsePoissonBracket)
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
    use SP_ModTime,         ONLY: SPTime
    use SP_ModAdvection,    ONLY: advance_log_advection,&
         advect_via_poisson_bracket
    use ModConst,           ONLY: cProtonMass, cGyroradius, Rsun, &
         cElectronCharge
    use SP_ModGrid,         ONLY: D_, Rho_, RhoOld_, B_, BOld_, U_, T_, &
         iPTest, iParticleTest, iNodeTest
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
    real, dimension(1:nVertexMax):: RadiusSi_I, DsSi_I,    &
         nSi_I, uSI_I, BSI_I, BOldSI_I, nOldSi_I
    real, dimension(1:nVertexMax):: InvRhoOld_I, InvRho_I

    ! Lagrangian derivatives
    real, dimension(1:nVertexMax):: DLogRho_I, FermiFirst_I

    ! Coefficients in the diffusion operator
    ! df/dt = DOuter * d(DInner * df/dx)/dx
    ! DOuter =BSI in the cell center
    ! DInner = DiffusionCoefficient/BSI at the face
    real, dimension(1:nVertexMax):: DOuterSI_I, CoefDInnerSI_I

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
       DtProgress = DtFull/nProgress

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

          if(UseTurbulentSpectrum)then
             ! Calculate the Alfven speed
             if (DoInitSpectrum) call init_spectrum(iEnd,              &
                  XyzSI_DI(x_:z_, 1:iEnd),BSI_I(1:iEnd), MomentumSI_I, &
                  dLogP, iShock, CoefInj, MachAlfven)
             ! call set_wave_advection_rates(iEnd,     &
                 ! BSI_I(1:iEnd),    BOldSI_I(1:iEnd),      &
                 ! nSi_I(1:iEnd)*cProtonMass,  nOldSI_I(1:iEnd)*cProtonMass, &
                 ! XyzSI_DI(x_:z_, 1:iEnd), DsSI_I(1:iEnd), &
                 ! DLogP, DtProgress, DtReduction)

             ! nStep = 1+int(max(DtReduction,                &
             !      maxval(abs(FermiFirst_I(1:iEnd))))/CFL)
             nStep = 1+int(maxval(abs(FermiFirst_I(2:iEnd)))/CFL)
          else
             ! How many steps should be done to the CFL criterion is fulfilled
             nStep = 1+int(maxval(abs(FermiFirst_I(2:iEnd)))/CFL)
          end if

          ! Check if the number of time steps is positive:
          if(nStep < 1)then
             if(UseTurbulentSpectrum) &
                  write(*,*) ' DtReduction               =', DtReduction
             write(*,*) ' maxval(abs(FermiFirst_I)) =', &
                  maxval(abs(FermiFirst_I(2:iEnd)))
             call CON_stop(NameSub//': nStep <= 0????')
          end if
          ! Store the value at the end of the previous time step
          call set_coef_diffusion
          if(UsePoissonBracket)then
             ! update bc for advection
             call set_advection_bc   
             ! store/update the inverse rho arrays
             InvRhoOld_I(1:iEnd) = 1.0/nOldSI_I(1:iEnd)
             InvRho_I(1:iEnd)    = 1.0/nSi_I(1:iEnd)
             call advect_via_poisson_bracket(nP, iEnd, & 
                 DtProgress,  Cfl,                     & 
                 InvRhoOld_I(1:iEnd), InvRho_I(1:iEnd),&
                 Distribution_IIB(:,1:iEnd,iLine))
             nOldSI_I(1:iEnd) = nSI_I(1:iEnd)
             BOldSI_I(1:iEnd) = BSI_I(1:iEnd)
             ! Call diffusion as frequent as you want or at the
             ! end of PROGRESS step only:
             ! diffusion along the field line

             if(UseDiffusion) call diffuse_distribution(nP,    &
                 iLine, iEnd, BSI_I, DOuterSI_I,               &
                 CoefDInnerSI_I, DsSI_I, DtProgress,           &
                 XyzSI_DI, nSI_I)
             DoInitSpectrum = .true.
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
                   if(any(Distribution_IIB(0:nP+1,iVertex,iLine) <0.0 )) then
                      write(*,*) NameSub, ': Distribution_IIB < 0'
                      Used_B(iLine) = .false.
                      nVertex_B(iLine) = 0
                      CYCLE line
                   end if
                   
                   call advance_log_advection(&
                        FermiFirst_I(iVertex), nP, 1, 1,       &
                        Distribution_IIB(0:nP+1,iVertex,iLine),&
                        .false.)
                end do

                if(UseDiffusion) call diffuse_distribution(nP, &
                    iLine, iEnd, BSI_I, DOuterSI_I,            &
                    CoefDInnerSI_I, DsSI_I, Dt,                &
                    XyzSI_DI, nSI_I)
             end do STEP
             DoInitSpectrum = .true.
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
           (uSI_I(  iShock + 1 - nWidth) - uSI_I(  iShock + nWidth))/    &
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
      DLogRho_I(iShock-nWidth:iShock+nWidth) = min(&
           DLogRhoBackground, &
           DLogRho_I(iShock-nWidth:iShock+nWidth))
      ! ...and concetrate it at the shock front, applying the whole jump
      ! in the velocity at a single grid point
      DLogRho_I(iShock) = DLogRhoBackground + DLogRhoExcessIntegral/&
           DsSI_I(iVertex)
      ! also, sharpen the magnetic field magnitude
      ! post shock part
      BSI_I(iShock+1-nWidth:iShock+1)=maxval(BSI_I(iShock+1-nWidth:iShock+1))
      ! pre shock part
      BSI_I(iShock+1:iShock+nWidth  )=minval(BSI_I(iShock+1:iShock+nWidth))
    end subroutine steepen_shock
    !==========================================================================
    subroutine set_coef_diffusion
      ! set diffusion coefficient for the current line
      real, dimension(1:nVertexMax):: ScaleSI_I
      real, parameter:: cCoef = 81./7/cPi/(2*cPi)**(2./3)

      !------------------------------------------------------------------------
      DOuterSI_I(1:iEnd) = BSI_I(1:iEnd)

      if(UseTurbulentSpectrum) RETURN

      ! precompute scale of turbulence along the line
      select case(iScaleTurbulenceType)
      case(Const_)
         ScaleSI_I(1:iEnd) = ScaleTurbulenceSI
      case(Linear_)
         ScaleSI_I(1:iEnd) = ScaleTurbulenceSI*RadiusSI_I(1:iEnd)/cAU
      end select
      ! Compute the diffusion coefficient without the contribution of
      ! v (velocity) and p (momentum), as v and p are different for
      ! different iP
      if(UseFixedMFPUpstream) then
         ! diffusion is different up- and down-stream
         ! Sokolov et al. 2004, paragraphs before and after eq (4)
         where(RadiusSI_I(1:iEnd) > 0.9 * RadiusSI_I(iShock))
            ! upstream: reset the diffusion coefficient to
            ! (1/3)*MeanFreePath0InAu[AU]*(R/1AU)*v*(pc/1GeV)^(1/3)
            ! see Li et al. (2003), doi:10.1029/2002JA009666
            ! ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
            ! 1/AU cancels with unit of Lambda0,no special attention needed;
            ! v (velocity) and (p)^(1/3) are calculated in momentum do loop
            CoefDInnerSI_I(1:iEnd) =                     &
                 (1.0/3)*MeanFreePath0InAU * RadiusSI_I(1:iEnd)&
                 *(cLightSpeed/cGEV)**(1.0/3)
         elsewhere
            CoefDInnerSI_I(1:iEnd) =  (cCoef/3)*BSI_I(1:iEnd)**2 /       &
                 (cMu*sum(MHData_VIB(Wave1_:Wave2_,1:iEnd,iLine),1))*    &
                 (ScaleSI_I(1:iEnd)**2*cGyroRadius/BSI_I(1:iEnd))**(1./3)
         end where
      else    ! .not.UseFixedMFPUpstream
         ! Sokolov et al., 2004: eq (4),
         ! note: Momentum = TotalEnergy * Vel / C**2
         ! Gyroradius = cGyroRadius * momentum / |B|
         ! DInner \propto (B/\delta B)**2*Gyroradius*Vel/|B|
         ! ------------------------------------------------------
         ! effective level of turbulence is different for different momenta:
         ! (\delta B)**2 \propto Gyroradius^(1/3)
         CoefDInnerSI_I(1:iEnd) =  (cCoef/3)*BSI_I(1:iEnd)**2 /       &
              (cMu*sum(MHData_VIB(Wave1_:Wave2_,1:iEnd,iLine),1))*    &
              (ScaleSI_I(1:iEnd)**2*cGyroRadius/BSI_I(1:iEnd))**(1./3)
      end if

      ! Add 1/B as the actual diffusion is D/B
      CoefDInnerSI_I(1:iEnd) = CoefDInnerSI_I(1:iEnd)/BSI_I(1:iEnd)

      ! set the left boundary condition for diffusion
      Distribution_IIB(1:nP+1, 1, iLine) = &
           Distribution_IIB(0, 1, iLine) * &
           (MomentumSI_I(0)/MomentumSI_I(1:nP+1))**SpectralIndex

    end subroutine set_coef_diffusion
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
  subroutine diffuse_distribution(nP, iLine, iEnd, BSI_I,      &
               DOuterSI_I, CoefDInnerSI_I, DsSI_I, Dt,         &
               XyzSI_DI, nSI_I)
      ! set up the diffusion coefficients 
      ! diffuse the distribution function 

      use ModConst, ONLY: cProtonMass
      use SP_ModSize, ONLY: nVertexMax
      use SP_ModDistribution, ONLY: SpeedSI_I,                 &
               Distribution_IIB, MomentumSI_I, DLogP 
      use SP_ModTurbulence, ONLY: UseTurbulentSpectrum, set_dxx, Dxx
      use SP_ModDiffusion, ONLY: advance_diffusion

      ! Variables as inputs (mandatory)
      integer, intent(in) :: nP           ! Momentum (P) grid size
      integer, intent(in) :: iLine, iEnd  ! input line/end indices
      real, intent(in), dimension(nVertexMax) :: BSI_I, &
               DOuterSI_I, CoefDInnerSI_I, DsSI_I
      real, intent(in) :: Dt              ! Time step for diffusion
      ! Variables as inputs (Optional)
      real, optional, intent(in) :: XyzSI_DI(3, nVertexMax)
      real, optional, intent(in) :: nSI_I(nVertexMax)

      ! Variables declared in this subroutine
      integer :: iP, iVertex              ! loop variables
      real :: DInnerSI_I(nVertexMax)      ! DInner Coefficients
      ! Full difference between DataInputTime and SPTime
      real, parameter :: DiffCoeffMinSI = 1.0E+04

      !------------------------------------------------------------------------
      ! diffusion along the field line

      if(UseTurbulentSpectrum)then
         call set_dxx(iEnd, nP, BSI_I(1:iEnd))
      end if

      MOMENTUM:do iP = 1, nP
         ! For each momentum account for dependence
         ! of the diffusion coefficient on momentum
         ! D\propto r_L*v\propto Momentum**2/TotalEnergy
         if (UseTurbulentSpectrum) then
            do iVertex=1, iEnd
               DInnerSI_I(iVertex) = Dxx(iVertex, iP,       &
                  MomentumSI_I(iP), SpeedSI_I(iP),          &
                  BSI_I(iVertex)) / BSI_I(iVertex)
            end do
         else
            ! Add v (= p*c^2/E_total in the relativistic case)
            ! and (p)^(1/3)
            DInnerSI_I(1:iEnd) = CoefDInnerSI_I(1:iEnd)     &
               *SpeedSI_I(iP)*(MomentumSI_I(iP))**(1.0/3)
            
            DInnerSI_I(1:iEnd) = max(DInnerSI_I(1:iEnd),    &
               DiffCoeffMinSI/DOuterSI_I(1:iEnd))
         end if
         
         call advance_diffusion(Dt, iEnd, DsSI_I(1:iEnd),   &
            Distribution_IIB(iP,1:iEnd,iLine),              &
            DOuterSI_I(1:iEnd), DInnerSI_I(1:iEnd))
      end do MOMENTUM
      
      ! if (UseTurbulentSpectrum) then
      !    call update_spectrum(iEnd,nP,MomentumSI_I,DLogP,   &
      !       XyzSI_DI(:,1:iEnd), DsSI_I(1:iEnd),             &
      !       Distribution_IIB(:,1:iEnd,iLine),BSI_I(1:iEnd), &
      !       nSi_I(1:iEnd)*cProtonMass,Dt)
      ! end if

    end subroutine diffuse_distribution
  !============================================================================
end module SP_ModAdvance
!==============================================================================
