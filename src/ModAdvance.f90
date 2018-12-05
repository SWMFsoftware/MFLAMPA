!  Copyright (C) 2002 Regents of the University of Michigan, 
!  portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!==============================================================================
module SP_ModAdvance
  ! The module contains methods for advancing the solution in time
  use ModNumConst,ONLY: cPi
  use ModConst,   ONLY: cLightSpeed, cGEV
  use SP_ModSize, ONLY: nParticleMax
  use SP_ModDistribution, ONLY: nP, Distribution_IIB, MomentumSI_I,           &
        EnergySI_I, MomentumMaxSI, MomentumInjSI, DLogP, SpeedSI_I
  use SP_ModGrid, ONLY: State_VIB, iShock_IB,   R_, x_, y_, z_, iNode_B,      &
       Shock_, ShockOld_, DLogRho_, Wave1_, Wave2_, nBlock, nParticle_B 
  use SP_ModTurbulence, ONLY: DoInitSpectrum, UseTurbulentSpectrum, set_dxx,  &
       set_wave_advection_rates, reduce_advection_rates, dxx, init_spectrum,  &
       update_spectrum
  use SP_ModUnit, ONLY: UnitX_, UnitRho_, UnitEnergy_,                        &
       NameParticle, IO2SI_I, kinetic_energy_to_momentum

  implicit none
  SAVE
  PRIVATE ! except
  !Public members:
  public:: read_param  !read injection parameters
  public:: advance     !Advance solution Distribution_IIB
  !\
  !||||||||||||||Boundary condition at the injection energy!!!!!!
  ! Injection efficiency and assumed spectral index with the energy
  ! range k_BT_i< Energy < EnergyInjection, to be read from PARAM.in 
  real:: CoefInj = 1.0, SpectralIndex = 5.0

  !/
  !\
  !!!!!!!!!!!!!!!!!!!!!!!!!Local parameters!!!!!!!!!!!!!!!
  real:: CFL=0.9        !Controls the maximum allowed time step
  !\
  integer, public, parameter :: nWidth = 50
  !/
  !\
  ! Diffusion as in Li et al. (2003), doi:10.1029/2002JA009666
  logical, public :: UseLiDiffusion = .false.
  real    :: Lambda0InAu = 1.0
  !/
  !\
  logical, public:: DoTraceShock = .true., UseDiffusion = .true.
  !/

  logical, public :: DoTestDiffusion = .false.

contains
  !============================================================================
  subroutine read_param(NameCommand)
    use ModReadParam, ONLY: read_var
    character(len=*), intent(in):: NameCommand ! From PARAM.in  
    character(len=*), parameter :: NameSub='SP:read_param_adv'
    !---------------------------------------------
    select case(NameCommand)
    case('#INJECTION')
       call read_var('SpectralIndex',SpectralIndex)
       call read_var('Efficiency',   CoefInj)
    case('#CFL')
       call read_var('Cfl',CFL)
    case('#USELIDIFFUSION')
       call read_var('UseLiDiffusion',UseLiDiffusion)
       if(UseLiDiffusion)then
          ! see Li et al. (2003), doi:10.1029/2002JA009666
          call read_var('Lambda0InAu [AU]', Lambda0InAu)
       end if
    case('#TESTDIFFUSION')
       call read_var('DoTestDiffusion', DoTestDiffusion)
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
    use SP_ModDiffusion,    ONLY: advance_diffusion
    use SP_ModLogAdvection, ONLY: advance_log_advection
    use ModConst,           ONLY: cMu, cProtonMass, cGyroradius, Rsun, &
         cElectronCharge
    use SP_ModGrid,         ONLY: D_, Rho_, RhoOld_, B_, BOld_, U_, T_, &
         iPTest, iParticleTest, iNodeTest
    real, intent(in):: TimeLimit
    !\
    ! Loop variables
    integer  :: iP, iParticle, iBlock
    !/
    !\
    !For a given block: nParticle_B, iShock_IB:
    integer  :: iEnd, iShock, iShockOld
    !/
    !\
    ! Upper limit and variable for the Loop which makes
    ! a time step so short that the shock wave passes
    ! a single grid interval per a progress step:
    integer  :: iProgress, nProgress
    ! Upper limit and variable for the loop which makes a 
    ! time step so short that the CFL for the (explicit) 
    ! advection is less that CFL declared abobe. 
    integer  :: iStep, nStep
    !/
    !\
    ! coefficient to interpolate  "old" and "new"
    real     :: Alpha
    !/
    !\
    ! Lower imit to floor the spatial diffusion coefficient For a
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
    real, parameter::  DiffCoeffMinSI =1.0E+04
    ! Full difference between DataInputTime and SPTime
    real      :: DtFull
    ! Time step in the PROGRESS Loop, DtFull/nProgress
    real      :: DtProgress 
    ! Time step in the STEP Loop, DtProgress/nStep
    real      :: Dt
    real      :: MachAlfven
    real      :: DtReduction

    ! Local arrays to store the position and state vectors in SI units
    real :: XyzSI_DI(3, 1:nParticleMax)
    real, dimension(1:nParticleMax):: RadiusSI_I, DsSI_I,    &
         RhoSI_I, RhoOldSI_I, nSI_I, uSI_I, BSI_I, BOldSI_I, &
         n_I ! n_I is in the IO unit

    ! Lagrangian derivatives 
    real, dimension(1:nParticleMax):: DLogRho_I, FermiFirst_I

    !\
    ! Coefficients in the diffusion operator
    ! df/dt = DOuter * d(DInner * df/dx)/dx
    real, dimension(1:nParticleMax):: &
         DOuterSI_I, DInnerSI_I, CoefDInnerSI_I
    !/

    character(len=*), parameter:: NameSub = 'SP_advance'
    !--------------------------------------------------------------------------
    ! the full time interval
    DtFull = TimeLimit - SPTime
    ! go line by line and advance the solution

    BLOCK:do iBlock = 1, nBlock
       ! the active particles on the line
       iEnd   = nParticle_B( iBlock)

       ! Various data along the line in SI units.
       ! The IO units of the state vectors could be seen in ModUnit, the 
       ! summary is: Length is in the unit of Rs, Rho is in the unit of
       ! amu/m^3, temperature is in the unit of kinetic energy, all others
       ! are in SI units. So to convert the IO units, only need three
       ! conversion factors are needed:
       ! UnitX_, UnitRho_, UnitEnergy_
       XyzSI_DI(x_,1:iEnd) = State_VIB(x_,     1:iEnd,iBlock)*IO2SI_I(UnitX_)
       XyzSI_DI(y_,1:iEnd) = State_VIB(y_,     1:iEnd,iBlock)*IO2SI_I(UnitX_)
       XyzSI_DI(z_,1:iEnd) = State_VIB(z_,     1:iEnd,iBlock)*IO2SI_I(UnitX_)
       DsSI_I(     1:iEnd) = State_VIB(D_,     1:iEnd,iBlock)*IO2SI_I(UnitX_)
       RadiusSI_I( 1:iEnd) = State_VIB(R_,     1:iEnd,iBlock)*IO2SI_I(UnitX_)
       uSI_I(      1:iEnd) = State_VIB(U_,     1:iEnd,iBlock)
       BOldSI_I(   1:iEnd) = State_VIB(BOld_,  1:iEnd,iBlock)
       RhoOldSI_I( 1:iEnd) = State_VIB(RhoOld_,1:iEnd,iBlock)*IO2SI_I(UnitRho_)

       ! find how far shock has travelled on this line: nProgress
       iShock    = iShock_IB(Shock_,   iBlock)
       iShockOld = iShock_IB(ShockOld_,iBlock)
       if(DoTraceShock)then
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

          ! Recall that State_VIB(Rho_/RhoOld_) is in the unit of amu/m^3.
          ! And we only consider protons, so in principle, 
          ! State_VIB(Rho_/RhoOld_) gives number density in the SI unit.
          ! nSI is needed to set up the distribution at the injection.
          n_I(1:iEnd) = State_VIB(RhoOld_,1:iEnd,iBlock) +Alpha *&
               (State_VIB(Rho_,  1:iEnd,iBlock) - &
               State_VIB(RhoOld_,1:iEnd,iBlock))

          nSI_I(1:iEnd)   = n_I(1:iEnd)
          RhoSI_I(1:iEnd) = n_I(1:iEnd)*IO2SI_I(UnitRho_)

          ! dimensionless
          DLogRho_I(1:iEnd)=log(RhoSI_I(1:iEnd)/RhoOldSI_I(1:iEnd))

          BSI_I(1:iEnd) = State_VIB(BOld_,1:iEnd,iBlock) + Alpha*&
               (State_VIB(B_,  1:iEnd,iBlock) - &
               State_VIB(BOld_,1:iEnd,iBlock))

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
             !Calculate the Alfven speed
             if (DoInitSpectrum) call init_spectrum(iEnd,              &
                  XyzSI_DI(x_:z_, 1:iEnd),BSI_I(1:iEnd), MomentumSI_I, &
                  dLogP, iShock, CoefInj, MachAlfven)
             ! call set_wave_advection_rates(iEnd,     &
             !      BSI_I(1:iEnd),    BOldSI_I(1:iEnd),      &
             !      RhoSI_I(1:iEnd),  RhoOldSI_I(1:iEnd),    &
             !      XyzSI_DI(x_:z_, 1:iEnd), DsSI_I(1:iEnd), &
             !      DLogP, DtProgress, DtReduction)
               
             ! nStep = 1+int(max(DtReduction,                &
             !      maxval(abs(FermiFirst_I(1:iEnd))))/CFL)
             nStep = 1+int(maxval(abs(FermiFirst_I(2:iEnd)))/CFL)
          else
             nStep = 1+int(maxval(abs(FermiFirst_I(2:iEnd)))/CFL)
          end if

          !\
          ! Check if the number of time steps is positive:
          !/
          if(nStep < 1) then
             if(UseTurbulentSpectrum) &
                  write(*,*) ' DtReduction               =', DtReduction
             write(*,*) ' maxval(abs(FermiFirst_I)) =', &
                  maxval(abs(FermiFirst_I(2:iEnd)))
             call CON_stop(NameSub//': nStep <= 0????')
          end if

          RhoOldSI_I(1:iEnd) = RhoSI_I(1:iEnd)
          BOldSI_I(  1:iEnd) = BSI_I(  1:iEnd)

          Dt = DtProgress/nStep

          FermiFirst_I(1:iEnd) = FermiFirst_I(1:iEnd) / nStep
          ! if(UseTurbulentSpectrum) call reduce_advection_rates(nStep)

          ! compute diffusion along the field line
          ! we calculate: "Outer diffusion"=BSI_I and
          ! "Inner diffusion at the injection energy (iP=0)
          call set_coef_diffusion

          STEP:do iStep = 1, nStep

             ! update bc for advection
             call set_advection_bc

             ! advection in the momentum space
             do iParticle = 2, iEnd
                if(any(Distribution_IIB(0:nP+1,iParticle,iBlock) <0.0 )) then
                   write(*,*) NameSub, ': Distribution_IIB < 0'
                   write(*,*) ' Distribution_IIB ='
                   write(*,*) Distribution_IIB(0:nP+1,iParticle,iBlock)
                   call CON_stop(NameSub)
                end if
                
                call advance_log_advection(&
                     FermiFirst_I(iParticle), nP, 1, 1,        &
                     Distribution_IIB(0:nP+1,iParticle,iBlock),&
                     .false.)
             end do
             if(.not.UseDiffusion) CYCLE STEP

             ! diffusion along the field line

             if(UseTurbulentSpectrum)then
                call set_dxx(iEnd, nP, BSI_I(1:iEnd))
             end if

             MOMENTUM:do iP = 1, nP
                ! For each momentum account for dependence
                ! of the diffusion coefficient on momentum
                ! D\propto r_L*v\propto Momentum**2/TotalEnergy
                if(.not. UseTurbulentSpectrum)then
                   ! Add v (= p*c^2/E_total in the relativistic case) 
                   ! and (p)^(1/3)
                   DInnerSI_I(1:iEnd) = CoefDInnerSI_I(1:iEnd) &
                        *SpeedSI_I(iP)*(MomentumSI_I(iP))**(1.0/3)
                   
                   ! Add 1/B as the actual diffusion is D/B
                   DInnerSI_I(1:iEnd) = DInnerSI_I(1:iEnd)/BSI_I(1:iEnd)

                   DInnerSI_I(1:iEnd) = max(DInnerSI_I(1:iEnd), &
                        DiffCoeffMinSI/DOuterSI_I(1:iEnd))
                else
                   do iParticle=1,iEnd
                      DInnerSI_I(iParticle) =                          &
                           Dxx(iParticle, iP, MomentumSI_I(iP),        &
                           SpeedSI_I(iP), BSI_I(iParticle))          / &
                           BSI_I(iParticle)
                   end do
                end if

                if (DoTestDiffusion .and. iNode_B(iBlock) == iNodeTest) then
                   if (iP == iPTest) then
                      write(*,'(a,i3,a,i3,a,es15.7)') NameSub//': iStep =', &
                           iStep, ' iProgress =',  iProgress,               &
                           ' DInnerSI =', DInnerSI_I(iParticleTest)
                   end if
                end if

                call advance_diffusion(Dt, iEnd,              &
                     DsSI_I(1:iEnd),                          &
                     Distribution_IIB(iP,1:iEnd,iBlock),      &
                     DOuterSI_I(1:iEnd), DInnerSI_I(1:iEnd))
             end do MOMENTUM

             ! if(UseTurbulentSpectrum .and. .true.)then
             !    call update_spectrum(iEnd,nP,MomentumSI_I,DLogP,      &
             !         XyzSI_DI(:,1:iEnd), DsSI_I(1:iEnd),              &
             !         Distribution_IIB(:,1:iEnd,iBlock),BSI_I(1:iEnd), &
             !         RhoSI_I(1:iEnd),Dt)
             ! end if

          end do STEP
          DoInitSpectrum = .true.
       end do PROGRESS
    end do BLOCK
  contains
    function mach_alfven() result(MachAlfven)
      ! alfvenic mach number for the current line
      real:: MachAlfven

      real:: SpeedAlfvenUpstream, SpeedUpstream
      !------------------------------------------------------------------------
      ! speed upstream is relative to the shock:
      ! \rho_u * (U_u - V_{shock}) = \rho_d * (U_d - V_{shock})
      SpeedUpstream = RhoSI_I(iShock+1-nWidth)*&
           (uSI_I(  iShock + 1 - nWidth) - uSI_I(  iShock + nWidth))/    &
           (RhoSI_I(iShock + 1 - nWidth) - RhoSI_I(iShock + nWidth))
      SpeedAlfvenUpstream = BSI_I(iShock + nWidth)/ &
           sqrt(cMu*RhoSI_I(iShock + nWidth))
      MachAlfven = SpeedUpstream / SpeedAlfvenUpstream
    end function mach_alfven
    !==========================================================================
    subroutine steepen_shock
      ! change the density profile near the shock front so it becomes steeper
      ! for the current line
      integer:: iParticle ! loop variable
      real   :: DLogRhoExcessIntegral, DLogRhoExcess
      real, parameter:: DLogRhoBackground = 0.01
      !------------------------------------------------------------------------
      ! find the excess of DLogRho within the shock compared to background
      ! averaged over length
      DLogRhoExcessIntegral = 0.0
      do iParticle = iShock - nWidth, iShock + nWidth - 1
         DLogRhoExcess = 0.5*(DLogRho_I(iParticle) + DLogRho_I(iParticle+1)) &
              - DLogRhoBackground !D log(rho)/Dt*\Delta t = -\div U*\Delta t
         if(DLogRhoExcess>0) then
            !This is a jump in velocity accross the shock wave * \Delta t
            DLogRhoExcessIntegral = DLogRhoExcessIntegral + &
                 DLogRhoExcess*DsSI_I(iParticle)
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
           DsSI_I(iParticle)
      ! also, sharpen the magnetic field magnitude
      ! post shock part
      BSI_I(iShock+1-nWidth:iShock+1)=maxval(BSI_I(iShock+1-nWidth:iShock+1))
      ! pre shock part
      BSI_I(iShock+1:iShock+nWidth  )=minval(BSI_I(iShock+1:iShock+nWidth))
    end subroutine steepen_shock
    !==========================================================================
    subroutine set_coef_diffusion
      use ModConst,   ONLY: cMu, cAu
      ! set diffusion coefficient for the current line
      !------------------------------------------------------------------------
      DOuterSI_I(1:iEnd) = BSI_I(1:iEnd)
      !DInner = DiffusionCoeff/B_)

      !\
      ! Compute the diffusion coefficient without the contribution of 
      ! v (velocity) and p (momentum), as v and p are different for 
      ! different iP
      !/
      if (.not. UseTurbulentSpectrum) then
         if(.not.UseLiDiffusion)then
            !\
            ! Sokolov et al., 2004: eq (4), 
            ! note: Momentum = TotalEnergy * Vel / C**2
            ! Gyroradius = cGyroRadius * momentum / |B|
            ! DInner = (B/\delta B)**2*Gyroradius*Vel/|B| 
            !/
            !\
            ! effective level of turbulence is different for different momenta:
            ! (\delta B)**2 \propto Gyroradius^(1/3)
            ! we scale it with reference point:  
            ! mean free path  ~ RSun @ 1AU for 1GeV proton
            ! (B/\delta B)**2 ~ 1    @ 1AU
            ! thus (e1,e2 are turbulent energies): 
            ! DInner=B^2/(2*cMu*(e1+e2))*10*(Gyroradius*RSun**2)**(1/3)*Vel/|B|
            ! where Gyroradius = cGyroRadius * p / |B|
            ! v (velocity) and (p)^(1/3) are calculated in the momentum do loop
            !/

            ! No idea how the Rsun**2 shows up in SI unit...
            CoefDInnerSI_I(1:iEnd) =  BSI_I(1:iEnd)**2 /                  &
                 (2*cMu*sum(State_VIB(Wave1_:Wave2_,1:iEnd,iBlock),1))    &
                 *10*(cGyroRadius/BSI_I(1:iEnd)*Rsun**2)**(1./3)! /Rsun**2
         else
            ! diffusion is different up- and down-stream
            ! Sokolov et al. 2004, paragraphs before and after eq (4)
            where(RadiusSI_I(1:iEnd) > 0.9 * RadiusSI_I(iShock))
               ! upstream: reset the diffusion coefficient to 
               ! (1/3)*Lambda0InAu[AU]*(R/1AU)*v*(pc/1GeV)^(1/3) 
               ! see Li et al. (2003), doi:10.1029/2002JA009666

               ! 1/AU cancels with unit of Lambda0,no special attention needed;
               ! v (velocity) and (p)^(1/3) are calculated in the 
               ! momentum do loop

               CoefDInnerSI_I(1:iEnd) =                    &
                    (1.0/3)*Lambda0InAU * RadiusSI_I(1:iEnd)       &
                    *(cLightSpeed/cGEV)**(1.0/3)
            elsewhere
               CoefDInnerSI_I(1:iEnd) =  BSI_I(1:iEnd)**2 /                  &
                    (2*cMu*sum(State_VIB(Wave1_:Wave2_,1:iEnd,iBlock),1))    &
                    *10*(cGyroRadius/BSI_I(1:iEnd)*Rsun**2)**(1./3)

               !\
               ! LEGACY MODEL PREVIOUSLY USED DOWNSTREAM
               !-------------------------------------------------
               ! downstream we use the estimate for (\delta B/B) 
               ! (see detail in Sokolov, 2004):
               ! (\delta B/B)**2 = (\delta B/B)**2_SW*(R/R_SW)
               ! where SW means "shock wave". For the SW we use an
               ! estimate from  from Lee (1983): 
               ! (\delta B/B)**2_SW = const(=10 below)*CoefInj*MachAlfven

               !!CoefDInnerSI_I(1:iEnd)=&
               !!     cGyroRadius*(MomentumInjSI*cLightSpeed)**2/RSun**2/&
               !!     (BSI_I(1:iEnd)**2*TotalEnergyInjSI)/&
               !!     (10.0*CoefInj*MachAlfven)/&
               !!     min(1.0, 1.0/0.9 * RadiusSI_I(1:iEnd)/RadiusSI_I(iShock))
            end where
         end if
      end if

      ! set the left boundary condition for diffusion
      Distribution_IIB(1:nP+1, 1, iBlock) = &
           Distribution_IIB(0, 1, iBlock) * &
           (MomentumSI_I(0)/MomentumSI_I(1:nP+1))**SpectralIndex

    end subroutine set_coef_diffusion
    !==========================================================================
    subroutine set_advection_bc
      ! set boundary conditions on grid point on the current line
      ! LOCAL VARIABLES:
      integer:: iParticle     ! loop variable
      real   :: MomentumSI    ! Momentum for the thermal energy k_BTi
      !------------------------------------------------------------------------
      do iParticle = 1, iEnd
         ! injection(Ti, Rho), see Sokolov et al., 2004, eq (3)
         ! f = CoefInj/8/pi * N / (2*m*T_p)^(3/2) * ((2*m*T_p)^(3/2)/p_inj)^5
         !   = CoefInj/8/pi * N / p^3 * (p/p_inj)^5
         ! where p = sqrt(2*m*T) is the momentum and T_p is the kinetic 
         ! temperature.
         MomentumSI = kinetic_energy_to_momentum(                  &
              State_VIB(T_,iParticle,iBlock)*IO2SI_I(UnitEnergy_))
         Distribution_IIB(0,iParticle,iBlock) =                    &
              CoefInj*1.0/(4*(SpectralIndex-3)*cPi)                &
              * nSI_I(iParticle)/MomentumSI**3                     &
              * (MomentumSI/MomentumInjSI)**SpectralIndex
      end do
    end subroutine set_advection_bc
    !==========================================================================
  end subroutine advance
  !============================================================================
end module SP_ModAdvance
