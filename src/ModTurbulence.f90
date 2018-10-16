!  Copyright (C) 2002 Regents of the University of Michigan, portions used 
!  with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
Module SP_ModTurbulence
  use ModConst
  use SP_ModDistribution, ONLY: nP
  implicit none
  SAVE

  private

  public :: init, finalize, DoInitSpectrum, UseTurbulentSpectrum, set_dxx, &
       read_param, set_wave_advection_rates, reduce_advection_rates, dxx,  &
       init_spectrum, update_spectrum

  !\
  ! Logicals, all .false. by default
  !/
  logical:: DoInitSpectrum              = .false.
  logical:: UseTurbulentSpectrum        = .false.
  logical:: UseAdvectionWithAlfvenSpeed = .false.

  ! Whether to save gamma as output
  ! logical:: DoOutputGamma  = .false.
  ! integer:: iXOutputStart  = 350
  ! integer:: iXOutputStride = 50
  ! integer:: iXOutputLast   = 1000
  ! integer:: nKOutput       = 30

  integer, parameter :: nK = nP
  real    :: DeltaK

  real, allocatable :: Gamma_I(:,:)
  real, allocatable :: IPlusSI_IX(:,:),IMinusSI_IX(:,:),ICSI_X(:)
  real, allocatable :: vAlfvenSI_I(:)
  real, allocatable :: kOverBSI_I(:)

  ! Rate of advection in k space, neglecting 
  ! a spacial advection with the Alfven speed
  real, allocatable:: DispersionA_I(:)

  !Rate of advection in k space, for I_+/I_- wave 
  real, allocatable:: DispersionPlus_I(:)
  real, allocatable:: DispersionMinus_I(:)

  ! This is the ratio of densities powered 3/2 
  real, allocatable:: RhoCompression_I(:)

  !------------------------------------------------------------------------!
  !          Grid in the momentum space                                    !
  !iP     0     1                         nP   nP+1                        !
  !       |     |    ....                 |     |                          !
  !P      P_inj P_inj*exp(\Delta (Ln P))  P_Max P_Max*exp(\Delta (Ln P))   !
  !             |    Grid in k-space      |     |                          !
  !K/B         KMax                      KMin                              !
  !ik    nP+1   nP                        1     0                          !
  !------------------------------------------------------------------------!

  real,allocatable,private:: AK_II(:,:)
  real,allocatable,private:: BK_II(:,:)

  real,allocatable,private:: CFL_I(:)

  integer,allocatable::      CorrectionMode_X(:)

  !the intensity of the back travelling wave in the initial condition
  real,parameter:: Alpha   = 1.0/10
  real,parameter:: Lambda0 = 4.0/10.0  ![AU]
contains
  !============================================================================
  subroutine read_param(NameCommand)
    use ModReadParam, ONLY: read_var
    character(len=*), intent(in):: NameCommand

    character(len=*), parameter :: NameSub = 'SP:read_param_Turbulence'
    !-------------------------------------------------------------------------

    select case(NameCommand)
    case('#INITSPECTRUM')
       call read_var('DoInitSpectrum', DoInitSpectrum)
    case('#TURBULENTSPECTRUM')
       call read_var('UseTurbulentSpectrum', UseTurbulentSpectrum)
       if (UseTurbulentSpectrum) DoInitSpectrum = .true.
    case('#ADVECTIONWITHALFVENSPEED')
       call read_var('UseAdvectionWithAlfvenSpeed', &
            UseAdvectionWithAlfvenSpeed)
    case default
       call CON_stop(NameSub//' Unknown command '//NameCommand)
    end select
  end subroutine read_param
  !============================================================================
  subroutine init
    ! Init all the allocatable vars
    use SP_ModSize, ONLY: nParticleMax
    use SP_ModProc, ONLY: iProc
    character(len=*), parameter :: NameSub = 'SP_init_turbulent'
    !-------------------------------------------------------------------------

    allocate(IPlusSI_IX(0:nP+1,1:nParticleMax), &
         IMinusSI_IX(0:nP+1,1:nParticleMax))
    allocate(kOverBSI_I(0:nP+1))
    allocate(ICSI_X(1:nParticleMax),CorrectionMode_X(1:nParticleMax))
    allocate(vAlfvenSI_I(1:nParticleMax))
    allocate(DispersionA_I(1:nParticleMax))

    allocate(RhoCompression_I(1:nParticleMax))

    ! if(DoOutputGamma) &
    !     allocate(Gamma_I(nP,(iXOutputLast-iXOutputStart)/iXOutputStride+1))

    allocate(DispersionPlus_I(1:nParticleMax))
    allocate(DispersionMinus_I(1:nParticleMax))
    allocate(CFL_I(1:nParticleMax))

    if (UseTurbulentSpectrum .and. .not. DoInitSpectrum) then
       DoInitSpectrum = .true. 
       if (iProc == 0) then
          write(*,*) NameSub, ': UseTurbulentSpectrum: '
          write(*,*) NameSub, ': DoInitSpectrum is switched to T'
       end if
    end if
  end subroutine init
  !============================================================================
  subroutine finalize
    character(len=*), parameter :: NameSub = 'SP_init_turbulent'

    !-------------------------------------------------------------------------
    deallocate(IPlusSI_IX,IMinusSI_IX)
    deallocate(kOverBSI_I)
    deallocate(ICSI_X,CorrectionMode_X)
    deallocate(vAlfvenSI_I)
    deallocate(DispersionA_I)

    deallocate(RhoCompression_I)

    ! if(DoOutputGamma) &
    !     allocate(Gamma_I(nP,(iXOutputLast-iXOutputStart)/iXOutputStride+1))

    deallocate(DispersionPlus_I)
    deallocate(DispersionMinus_I)
    deallocate(CFL_I)

  end subroutine finalize
  !============================================================================
  subroutine init_spectrum(iEnd, XyzSI_DI, BSI_I, MomentumSI_I, iShock, &
       CoefInj, AlfvenMach)
    !==============Initial spectrum of turbulence=============================!
    ! We recover the initial spectrum of turbulence from the spatial 
    ! distribution of the diffusion coefficient and its dependence on the 
    ! particle energy.

    ! the number of active particles on the line
    integer, intent(in)::  iEnd

    ! Coordinates of Lagrangian Meshes in SI unit [m]
    real,dimension(1:3,1:iEnd),intent(in) :: XyzSI_DI

    ! Magnetic field intensity in SI unit [T]
    real,dimension(1:iEnd),intent(in) :: BSI_I

    ! momentum in SI unit
    real,intent(in)     :: MomentumSI_I(0:nP+1)

    ! coef of injection 
    real,intent(in)     :: CoefInj

    ! Alfven March number
    real,intent(in)     :: AlfvenMach

    ! shock index
    integer, intent(in) :: iShock

    integer :: iParticle,iK
    real    :: ICOldSI, kSI
    real    :: rSI , rShockSI

    !-------------------------------------------------------------------------!

    IPlusSI_IX    = 0.0; IMinusSI_IX     = 0.0
    vAlfvenSI_I = 0.0; DispersionA_I = 0.0

    RhoCompression_I = 0.0

    if(UseAdvectionWithAlfvenSpeed)then
       DispersionPlus_I  = 0.0
       DispersionMinus_I = 0.0
       CFL_I = 1.0
    end if

    !-------------------------------------------------------------------------!
    !          Grid in the momentum space                                     !
    !iP     0     1                         nP   nP+1                         !
    !       |     |    ....                 |     |                           !
    !P      P_inj P_inj*exp(\Delta (Ln P))  P_Max P_Max*exp(\Delta (Ln P))    !
    !             |    Grid in k-space      |     |                           !
    !K/B         KMax                      KMin                               !
    !ik    nP+1   nP                        1     0                           !
    !-------------------------------------------------------------------------!

    kOverBSI_I = cElectronCharge/MomentumSI_I

    rShockSI = sqrt(sum(XyzSI_DI(:,iShock)**2))

    call assign_kolmogorov_spectrum(1,iEnd,0,nP+1,XyzSI_DI(:, 1:iEnd), &
         BSI_I(1:iEnd))

    CorrectionMode_X=1.0

    do iParticle=1,iEnd
       rSI = sqrt(sum(XyzSI_DI(:,iParticle)**2))

       if (rSI < 1.1*rShockSI) then
          ! In this part of the spectrum another equation governs the diffusion
          ICOldSI    = ICSI_X(iParticle)

          ICSI_X(iParticle) = CoefInj *                                 &
               10.0*BSI_I(iParticle)**2*max(AlfvenMach,2.0)/(cMu*3.0) * &
               min(1.0, rSI/rShockSI/(1.0-Alpha))

          if(ICSI_X(iParticle) > ICOldSI)then
             do iK=0,nK+1
                KSI = kOverBSI_I(iK)*BSI_I(iParticle)

                IPlusSI_IX( iK,iParticle) = &
                     (1.0-Alpha)*ICSI_X(iParticle)/KSI
                IMinusSI_IX(iK,iParticle) = &
                     Alpha/(1.0-Alpha)*IPlusSI_IX(iK,iParticle)
             end do
             CorrectionMode_X(iParticle)=2
          else
             ! Do not change if the newly calculated value is less than the 
             ! old one
             ICSI_X(iParticle) = ICOldSI
          end if
       end if
    end do

    ! Need to update when B changes?
    ! DoInitSpectrum = .false.

  end subroutine init_spectrum
  !============================================================================
  subroutine assign_kolmogorov_spectrum(iParticleFirst, iParticleLast, &
       iKFirst, iKLast, XyzSI_DI, BSI_I)

    integer,intent(in):: iParticleFirst, iParticleLast, iKFirst, iKLast
    real,   intent(in):: XyzSI_DI(1:3,iParticleFirst:iParticleLast)
    real,   intent(in):: BSI_I(iParticleFirst:iParticleLast)

    integer :: iParticle, iK
    real    :: kSI, kR0SI, rSI
    !------------------------------------------------------------------------
    do iParticle=iParticleFirst,iParticleLast
       ! We use the formulae:
       ! I_{+}(k)=(1-\Alpha)*IC/K^{5/3} and I_{-}=\Alpha*I_{+}/(1-\Alpha),
       ! where IC=54*B^2/(7*\pi* \Lambda_0* K0^{1/3}*(r/1 [AU])) and 
       ! \Lambda_0=0.4 [AU], and K0=|e|B/1[GeV], in CGS units

       kR0SI = cElectronCharge*BSI_I(iParticle)*cLightSpeed/energy_in('GeV')

       rSI = sqrt(sum(XyzSI_DI(:,iParticle)**2))

       ICSI_X(iParticle)=(216./7.)*(BSI_I(iParticle)**2)/(Lambda0*cMu)    &
            /kR0SI**(1.0/3)/rSI

       do iK=iKFirst,iKLast
          ! kOverBSI_I = cElectronCharge/MomentumSI_I
          kSI = kOverBSI_I(iK)*BSI_I(iParticle)

          IPlusSI_IX( iK,iParticle) = &
               (1.0-Alpha)*ICSI_X(iParticle)/kSI**(5.0/3)
          IMinusSI_IX(iK,iParticle) = &
               Alpha/(1.0-Alpha)*IPlusSI_IX(iK,iParticle) 
       end do
    end do
  end subroutine assign_kolmogorov_spectrum
  !============================================================================
  subroutine set_wave_advection_rates(iEnd, BSI_I, BOldSI_I, RhoSI_I, &
    RhoOldSI_I, XyzSI_DI, DLogP, Dt, DtReduction)
    !=======================Set advection rate in k-space======================
    integer,intent(in) :: iEnd
    real,intent(in)    :: BSI_I(1:iEnd), BOldSI_I(1:iEnd)
    real,intent(in)    :: RhoSI_I(1:iEnd), RhoOldSI_I(1:iEnd)
    real,intent(in)    :: XyzSI_DI(1:3,1:iEnd)
    real,intent(in)    :: DLogP, Dt
    real,intent(out)   :: DtReduction

    ! local vars
    integer :: iParticle, DsSI, DLogRho, DLogB
    !-------------------------------------------------------------------------

    ! Calculate alfven speed in SI unit
    vAlfvenSI_I = BSI_I/sqrt(cMu*RhoSI_I)

    RhoCompression_I=1.50*log(RhoSI_I/RhoOldSI_I) - log(BSI_I/BOldSI_I)

    ! Contribution to the advection in k_space from the Lagrangian derivatives:
    ! Dispersion = Dt * [ (D ln {\rho})/(Dt) - 2 (D ln { B })/(D t) ] + ...

    DispersionA_I = log(RhoSI_I*BOldSI_I**2/(RhoOldSI_I*BSI_I**2)) 

    if(UseAdvectionWithAlfvenSpeed)then
       DsSI        = sqrt(sum((XyzSI_DI(:,2)-XyzSI_DI(:,1))**2))
       CFL_I(1)  = Dt*vAlfvenSI_I(1 )/DsSI

       DsSI        = sqrt(sum((XyzSI_DI(:,iEnd)-XyzSI_DI(:,iEnd-1))**2))
       CFL_I(iEnd) = Dt*vAlfvenSI_I(iEnd)/DsSI

       do iParticle=2,iEnd-1
          DsSI        = sqrt(sum( &
               (XyzSI_DI(:,iParticle+1)-XyzSI_DI(:,iParticle-1))**2 )) /2.0
          CFL_I(iParticle) = Dt*vAlfvenSI_I(iParticle)/DsSI
       end do

       !Now add the contribution from the spacial derivative:
       ! ... + Dt * [ \partial ((1/2) ln {rho})/(\partial s)-2\partial( ln{B})/(\partial s) ]  
       !for the "plus" wave 

       DLogRho = log(RhoSI_I(2)/ RhoSI_I(1))
       DLogB   = log(BSI_I(2)  / BSI_I(1))
       DispersionPlus_I(1)  =                                       &
            DispersionA_I(1)+CFL_I(1 )*(0.5*DLogRho-2.0*DLogB)

       DLogRho = log(RhoSI_I(iEnd)/ RhoSI_I(iEnd-1))
       DLogB   = log(BSI_I(iEnd)  / BSI_I(iEnd-1))
       DispersionPlus_I(iEnd) =                                     &
            DispersionA_I(iEnd)+CFL_I(iEnd)*(0.5*DLogRho-2.0*DLogB)

       do iParticle=2,iEnd-1
          DLogRho = log(RhoSI_I(iParticle+1)/ RhoSI_I(iParticle-1))/2
          DLogB   = log(BSI_I(iParticle+1)  / BSI_I(iParticle-1))/2

          DispersionPlus_I(iParticle) = DispersionA_I(iParticle)    &
               + CFL_I(iParticle)*(0.5*DLogRho-2.0*DLogB)
       end do

       !Now add the contribution from the spacial derivative:
       ! ... - Dt * [ \partial ((1/2) ln {rho})/(\partial s)-2\partial( ln{B})/(\partial s) ]  
       !for the "minus" wave 

       DispersionMinus_I(1:iEnd) = 2.0*DispersionA_I(1:iEnd) &
            - DispersionPlus_I(1:iEnd)

       DispersionPlus_I(1:iEnd)  = DispersionPlus_I(1:iEnd)  / DLogP
       DispersionMinus_I(1:iEnd) = DispersionMinus_I(1:iEnd) / DLogP

       DtReduction=max(maxval(CFL_I(1:iEnd)), &
            maxval(abs(DispersionPlus_I(1:iEnd))), &
            maxval(abs(DispersionMinus_I(1:iEnd))))
    else
       DispersionA_I(1:iEnd) = DispersionA_I(1:iEnd) / DLogP
       DtReduction   = maxval(abs(DispersionA_I(1:iEnd)))
    end if
  end subroutine set_wave_advection_rates
  !============================================================================
  subroutine reduce_advection_rates(nStep)
    integer, intent(in) :: nStep
    !-------------------------------------------------------------------------
    if(UseAdvectionWithAlfvenSpeed) then
       DispersionPlus_I  = DispersionPlus_I  /real(nStep)
       DispersionMinus_I = DispersionMinus_I /real(nStep)

       CFL_I=CFL_I/real(nStep)
    else
       DispersionA_I=DispersionA_I/real(nStep)
    end if
    RhoCompression_I=RhoCompression_I/real(nStep)
  end subroutine reduce_advection_rates
  !============================================================================
  subroutine set_dxx(iEnd, nP, BSI_I)
    ! Calculate the partial sums as described above: 
    !     \sum_{k_{res}}^\infty d(log K)/(k^2*(I_+ + I_-))  (This is AK_I) 
    ! and \sum_{k_{res}}^\infty d(log K)/(k^4*(I_+ + I_-))  (This is BK_I)

    integer,intent(in) :: iEnd,nP
    real,intent(in)    :: BSI_I(iEnd)
    !-------------------------------------------------------------------------
    integer:: iParticle,iK
    real :: F01,F02,F11,F12
    real :: K0,K1,expdLogK
    real :: SpectralIndexAtKMax

    if (.not.allocated(AK_II)) then
       allocate(AK_II(nP,iEnd),BK_II(nP,iEnd)) 
    end if

    !-------------------------------------------------------------------------!
    !          Grid in the momentum space                                     !
    !iP     0     1                         nP   nP+1                         !
    !       |     |    ....                 |     |                           !
    !P      P_inj P_inj*exp(\Delta (Ln P))  P_Max P_Max*exp(\Delta (Ln P))    !
    !             |    Grid in k-space      |     |                           !
    !K/B         KMax                      KMin                               !
    !ik    nP+1   nP                        1     0                           !
    !-------------------------------------------------------------------------!

    do iParticle=1,iEnd
       select case(CorrectionMode_X(iParticle))
       case(1)
          SpectralIndexAtKMax = 1.0+2.0*(1.0/3)
       case(2)
          SpectralIndexAtKMax = 1.0
       end select

       K0=BSI_I(iParticle)*maxval(kOverBSI_I)

       !The integrand for AK_I
       F01=1.0/(K0**2)/(IPlusSI_IX(nP,iParticle)+IMinusSI_IX(nP,iParticle))

       !The integrand for BK_I 
       F02=1.0/(K0**4)/(IPlusSI_IX(nP,iParticle)+IMinusSI_IX(nP,iParticle))

       ExpDLogK= 1 !exp(DeltaLnK)

       !As the starting values for AK_I and BK_I at the maximal momentum,
       !solve the integrals from K_{max} up to \infty, assuming the
       !power law spectrum of turbulence at K>K_{max}
       AK_II(nP,iParticle)=F01/(2.0-SpectralIndexAtKMax)
       BK_II(nP,iParticle)=F02/(4.0-SpectralIndexAtKMax)

       do iK=nP-1,1,-1
          !We calculate the partial sums for a set of the wave number values.
          !The integral is taken from KRes up to infinity, so we start from 
          !the maximal wave number and add the contributions from each of 
          !the wave number intervals.

          !The lower value of the wave number
          K1=K0/expdLogK 

          !The integrands at the lower value of the wave number
          F11=1.0/(K1**2)/(IPlusSI_IX(iK,iParticle)+IMinusSI_IX(iK,iParticle)) 
          F12=1.0/(K1**4)/(IPlusSI_IX(iK,iParticle)+IMinusSI_IX(iK,iParticle))

          !Calculate the new partial sums

          AK_II(iK,iParticle)=AK_II(iK+1,iParticle)+cHalf*(F01+F11)! *DeltaLnk
          BK_II(iK,iParticle)=BK_II(iK+1,iParticle)+cHalf*(F02+F12)! *DeltaLnK 

          ! current values saved as the initial values for the next step in 
          ! the loop
          K0=K1; F01=F11; F02=F12
       end do
    end do
  end subroutine set_dxx
  !===========================================================================
  subroutine update_spectrum(iEnd,nP,PInjectionSI,DLogP,&
       XyzSI_DI,F_II,BSI_I,RhoSI_I,SP_Dt)

    !This is the subroutine which advances the wave spectrum through
    !a time step, by solving the equations
    !dI_+/dt=\gamma_+ I_+
    !dI_-/dt=\gamma_- I_-

    use SP_ModLogAdvection, ONLY: advance_log_advection
    use ModLinearAdvection
    !-----------------------------------------------------------------
    !The number of points and the number of the energy intervals
    integer,intent(in)::iEnd,nP

    real,intent(in) :: PInjectionSI,DLogP,SP_Dt  
    !Coordinates of the Lagrangian points
    real,dimension(3,iEnd),intent(in)::XyzSI_DI
    !The distribution function
    real,dimension(0:nP+1,1:iEnd),intent(in)::F_II
    !The magnetic field intensity and the particle number density, [T],[m^{-3}]
    real,dimension(iEnd),intent(in)::BSI_I,RhoSI_I

    integer::iParticle,iK,iP

    !Resonant value of the particle momentum, for a given K
    real:: PRes

    ! The increment for the I_+ wave, multiplied by (I_+  +  I_-). Below 
    ! the dynamical equations for I_+ and I_- are reformulated in terms of 
    ! \gamma and integrated exactly assuming the constant value of \gamma 
    ! (NOT the increment) through the time step 
    real:: Gamma

    !Wave number
    real:: K

    !Forward and backward values for Ds 
    real:: DsPlusSI,DsMinusSI

    !Spatial derivatives of the distribution function at given values 
    ! of the momentum
    real:: DfDs0,P0,DfDs1,P1

    !Partial sums in the integral for Gamma
    real:: A(nP+1),B(nP+1)

    !Variables used in the spectral model for turbulence
    real:: ExpRhoCompression

    !Miscellaneous 
    real::C1,C2,C3

    do iParticle=1,iEnd

       !Advection in k space:
       if(UseAdvectionWithAlfvenSpeed)then
          call advance_log_advection( &
               DispersionPlus_I( iParticle), nP, 1, 1,         &
               IPlusSI_IX( :,iParticle),IsConservative=.true., &
               DeltaLnP=DLogP)
          call advance_log_advection(DispersionMinus_I(iParticle),nP,1,1,&
               IMinusSI_IX(:,iParticle),IsConservative=.true.,&
               DeltaLnP=DLogP)
       else

          call advance_log_advection(DispersionA_I(iParticle),nP,0,0,&
               IPlusSI_IX( 1:nP,iParticle),IsConservative=.true.,&
               DeltaLnP=DLogP)
          call advance_log_advection(DispersionA_I(iParticle),nP,0,0,&
               IMinusSI_IX(1:nP,iParticle),IsConservative=.true.,&
               DeltaLnP=DLogP)
       end if

       !----------------------------------------------------------------------!
       !          Grid in the momentum space                                  !
       !iP     0     1                         nP   nP+1                      !
       !       |     |    ....                 |     |                        !
       !P      P_inj P_inj*exp(\Delta (Ln P))  P_Max P_Max*exp(\Delta (Ln P)) !
       !             |    Grid in k-space      |     |                        !
       !K/B         KMax                      KMin                            !
       !ik    nP+1   nP                        1     0                        !
       !----------------------------------------------------------------------!

       !Calculate the partial sums in the integral for \Gamma
       !Nullify the sums

       A(nP+1)=0.0; B(nP+1)=0.0

       ! Calculate the spatial derivatives of the distribution function 
       ! at the maximal energy

       if (iParticle==1) then

          DsPlusSI = max(cTiny, sqrt( &
               sum((XyzSI_DI(:,iParticle+1)-XyzSI_DI(:,iParticle))**2)))

          !Use the forward spatial derivative
          DfDs0=(F_II(nP,iParticle+1)-F_II(nP,iParticle))/DsPlusSI

       else if (iParticle==iEnd) then

          DsMinusSI = max(cTiny, sqrt( &
               sum( (XyzSI_DI(:,iParticle)-XyzSI_DI(:,iParticle-1))**2)))

          !Use the backward spatial derivative
          DfDs0=(F_II(nP,iParticle)-F_II(nP,iParticle-1))/DsMinusSI

       else

          DsPlusSI = max(cTiny, sqrt( &
               sum((XyzSI_DI(:,iParticle+1)-XyzSI_DI(:,iParticle))**2)))
          DsMinusSI = max(cTiny, sqrt( &
               sum( (XyzSI_DI(:,iParticle)-XyzSI_DI(:,iParticle-1))**2)))

          !Use the average between the forward and backward spatial derivatives
          DfDs0=cHalf*((F_II(nP,iParticle+1)-F_II(nP,iParticle  ))/DsPlusSI+&
               (F_II(nP,iParticle  )-F_II(nP,iParticle-1))/DsMinusSI)

       end if

       !The momentum at the maximal energy

       P0=PInjectionSI*exp(real(nP+1)*DLogP)  

       ! We calculate the partial sums for a set of the momentum values.
       ! The integral is taken from pRes up to infinity, so we start from the
       ! maximal energy and add the contributions from each of the energy 
       ! intervals.

       do iP=nP,1,-1

          !Calculate the momentum at the lower energy
          P1=PInjectionSI*exp(real(iP)*DLogP)  

          !Calculate the distribution function gradient at the lower energy

          if (iParticle == 1) then
             DfDs1 = (F_II(iP,iParticle+1)-F_II(iP,iParticle  ))/DsPlusSI
          else if (iParticle == iEnd) then
             DfDs1 = (F_II(iP,iParticle  )-F_II(iP,iParticle-1))/DsMinusSI
          else
             DfDs1 = cHalf* &
                  ( (F_II(iP,iParticle+1)-F_II(iP,iParticle  ))/DsPlusSI  +&
                  (  F_II(iP,iParticle  )-F_II(iP,iParticle-1))/DsMinusSI )
          end if

          ! here are the parts for \gamma being integrated 
          ! For I_+ the increment equals \frac{\gamma}{(I_{+}+I_{-})}
          ! where \gamma =-4.0*(cPi**2)*Va/k *Integral/cProtonMass 
          ! Integral=\int_{p_{res}(k)}^{\infty}{\frac{dp}{p}*p^4
          ! (p_{res}(k)-\frac{p_{res}^3}{p^2})\frac{\partial{f}}{\partial{s}}
          !
          ! A=p_{res}^{k}=\int{d(ln{p})*(p^{4}*{\frac{\partial{f}}{\partial{s}}}}
          ! B=p_{res}^{3}\int{d(ln{p})*p^2*{\frac{\partial{f}}{\partial{s}}

          A(iP)=A(iP+1)+cHalf*(DfDs0*(P0**4)+DfDs1*(P1**4))*DLogP
          B(iP)=B(iP+1)+cHalf*(DfDs0*(P0**2)+DfDs1*(P1**2))*DLogP

          ! Save the values for the lower energy to re-use them as the values 
          ! for the higher energy end of the next interval (in CYCLING the 
          ! momentummo DECREASES)

          P0=P1; DfDs0=DfDs1

       end do

       !Calculate the wave increment and update the wave spectra
       ExpRhoCompression=exp(RhoCompression_I(iParticle))

       call assign_kolmogorov_spectrum( &
            iParticle, iParticle, 0, 0, XyzSI_DI(:,iParticle:iParticle),&
            BSI_I(iParticle:iParticle))
       !IPlusSI_IX(    0,iParticle) = 
       !       IPlusSI_IX(  0,iParticle)*ExpRhoCompression
       !IMinusSI_IX(   0,iParticle) = 
       !       IMinusSI_IX( 0,iParticle)*ExpRhoCompression

       !----------------------------------------------------------------------!
       !          Grid in the momentum space                                  !
       !iP     0     1                         nP   nP+1                      !
       !       |     |    ....                 |     |                        !
       !P      P_inj P_inj*exp(\Delta (Ln P))  P_Max P_Max*exp(\Delta (Ln P)) !
       !             |    Grid in k-space      |     |                        !
       !K/B         KMax                      KMin                            !
       !ik    nP+1   nP                        1     0                        !
       !----------------------------------------------------------------------!

       do iK=1,nK
          !The wave number
          K=BSI_I(iParticle)*kOverBSI_I(iK)

          !The resonant momentum, for a given K
          PRes=cElectronCharge*BSI_I(iParticle)/K

          iP=nP+1-iK

          ! here is the dynamic \gamma evaluated

          Gamma=-4.0*2.0*(cPi**2)*vAlfvenSI_I(iParticle)/K*&
               (PRes*A(iP)-(PRes**3)*B(iP))/    &
               cProtonMass 

          ! if(i_output(iParticle)/=0)Gamma_I(iK,i_output(iParticle))=Gamma

          ! We need to integrate the two coupled equations:
          !
          ! DI_+/Dt =(3/2)* (D ln rho/Dt) * I_+  +  \gamma/(I_+ + I_-) * I_+
          ! and
          ! DI_-/Dt =(3/2)* (D ln rho/Dt)) *I_-  -  \gamma/(I_+ + I_-) * I_-
          ! see Eq.(\ref{eq:Lagrangian})
          !
          ! Make a substitution 
          ! I_+ = I_+(new) * rho^{3/2},    I_- = I_-(new) * rho^{3/2}
          !
          ! The modified equations do not involve the density derivative on 
          ! the right hand side and are given below. The solution will be 
          ! expressed in terms of I_+ * I_- and I_+  - I_-, in which we 
          ! multiply the "old" values of intinsity by the compression ratio 
          ! to account for the contribution from the density derivative.

          ! Compression  factor comes twice
          C1=(IPlusSI_IX(iK,iParticle)*IMinusSI_IX(iK,iParticle)) &
               *ExpRhoCompression**2
          C2=(IPlusSI_IX(iK,iParticle)-IMinusSI_IX(iK,iParticle)) &
               *ExpRhoCompression

          ! The solution of the equations
          ! dI_+/dt=(\gamma/(I_+ + I_-)*I_+
          ! and
          ! dI_-/dt=-(\gamma/(I_+  +  I_-)*I_-
          ! reads: (I_+   -   I_-)^{n+1}=(I_+   -   I_-)^n + \gamma*dt,
          ! so calculate the right hand side:

          C3=gamma*SP_Dt+C2

          ! Having in mind the conservation of the c1= I_-*I_+ product, 
          ! solve the quadratic equation for I_+ :   I_+^2-c3*I_+ - c1=0.
          ! To avoid the loss in accuracy, occuring at C3>>C1, calculate that
          ! root of the quadratic equation which may be obataibed as a total 
          ! of two POSITIVE numbers

          if(C3>0)then

             IPlusSI_IX(iK,iParticle)=(C3+sqrt(C3**2+4.0*C1))/2.0 

             ! Now we take again into account the conservation of the 
             ! product I_+*I_-
             ! which is the consequence of the relationship \gamma_-=-\gamma_+

             IMinusSI_IX(iK,iParticle)=C1/IPlusSI_IX(iK,iParticle)
          else

             IMinusSI_IX(iK,iParticle)=(-C3+sqrt(C3**2+4.0*C1))/2.0 

             ! Now we take again into account the conservation of the product 
             ! I_+*I_-
             ! which is the consequence of the relationship \gamma_-=-\gamma_+

             IPlusSI_IX(iK,iParticle)=C1/IMinusSI_IX(iK,iParticle)
          end if

          if(IPlusSI_IX(iK,iParticle)<=0.0)then
             write(*,*)'IPlus(iK,iParticle) < 0, IPlus,iK,iParticle=',&
                  IPlusSI_IX(iK,iParticle),iK,iParticle,c1,c2,c3,gamma,SP_Dt
             stop
          end if

          if(IMinusSI_IX(iK,iParticle)<=0.0)then
             write(*,*)'IMinus(iK,iParticle) < 0, IMinus,iK,iParticle=',&
                  IMinusSI_IX(iK,iParticle),iK,iParticle,c1,c2,c3,gamma,SP_Dt
             stop
          end if

       end do    !cycling iK
       IPlusSI_IX( nP+1,iParticle) = IPlusSI_IX( nP,iParticle)
       IMinusSI_IX(nP+1,iParticle) = IMinusSI_IX(nP,iParticle)
    end do       !cycling iParticle

    if(UseAdvectionWithAlfvenSpeed)then
       do iK=0,nP+1
          IPlusSI_IX(iK,1:iEnd) = &
               vAlfvenSI_I(1:iEnd)*IPlusSI_IX(iK,1:iEnd)/BSI_I( 1:iEnd)
          call advance_lin_advection_plus(CFL_I( 1:iEnd),&
               iEnd,0,0,IPlusSI_IX(iK,1:iEnd))

          IPlusSI_IX(iK,1:iEnd) = &
               BSI_I( 1:iEnd)*IPlusSI_IX(iK,1:iEnd)/vAlfvenSI_I(1:iEnd)

          IMinusSI_IX(iK,1:iEnd)= &
               vAlfvenSI_I(1:iEnd)*IMinusSI_IX(iK,1:iEnd)/BSI_I( 1:iEnd)

          call advance_lin_advection_minus(CFL_I(1:iEnd),&
               iEnd,0,0,IMinusSI_IX(iK,1:iEnd))

          IMinusSI_IX(iK,1:iEnd) = &
               BSI_I( 1:iEnd)*IMinusSI_IX(iK,1:iEnd)/vAlfvenSI_I(1:iEnd)
       end do
    end if
  end subroutine update_spectrum
  !============================================================================
  ! integer function i_output(iParticle)
  !   integer,intent(in)::iX
  !   !---------------------!
  !   i_output=0
  !   if( .not. DoOutputGamma .or. &
  !        iX>iXOutputLast.or.iX<iXOutputStart) return
  !   if(( (iX-iXOutputStart)/iXOutputStride )*iXOutputStride /= &
  !        (iX-iXOutputStart))return
  !   i_output=(iX-iXOutputStart)/iXOutputStride+1
  ! end function i_output
  !============================================================================
  real function Dxx(iX, iP, MomentumSI, BSI, NameParticle)
    !=====The coefficient of a spatial diffusion along the magnetic field===!
    !Given by a formula:

    !D_{xx}=4*v*B^2/(cMu)*(AK_I(k_{res})-k_{res}^2*BK_I(k_{res})),

    !where AK_I=\int_{k_{res}}^\infty{d(\log k)/(k^2*(I_{+}(k) + I_{-}(k) )) }
    !and   BK_I=\int_{k_{res}}^\infty{d(\log k)/(k^4*(I_{+}(k) + I_{-}(k) )) }
    !--------------------------------------------------------------------
    integer,intent(in) :: iX, iP
    real,   intent(in) :: MomentumSI, BSI
    character(LEN=*),intent(in) ::NameParticle
    !--------------------------------------------------------------------
    real              :: kRSI         !The resonant wave number 
    integer           :: iK

    !-------------------------------------------------------------------------!
    !          Grid in the momentum space                                     !
    !iP     0     1                         nP   nP+1                         !
    !       |     |    ....                 |     |                           !
    !P      P_inj P_inj*exp(\Delta (Ln P))  P_Max P_Max*exp(\Delta (Ln P))    !
    !             |    Grid in k-space      |     |                           !
    !K/B         KMax                      KMin                               !
    !ik    nP+1   nP                        1     0                           !
    !-------------------------------------------------------------------------!

    !Calculate the resonant wave number
    kRSI = cElectronCharge*BSI/MomentumSI

    if (abs(kRSI - kOverBSI_I(iP)*BSI) > 1e-9) then
       write(*,*) 'iP =', iP
       write(*,*) 'kOverBSI_I(iP) =', kOverBSI_I(iP)
       write(*,*) 'k/B here       =', cElectronCharge/MomentumSI
       write(*,*) 'MomentumSI     =', MomentumSI
       write(*,*) 'diff           =', abs(kRSI - kOverBSI_I(iP)*BSI)
       write(*,*) 'Shit happens'
       STOP
    end if

    ! The integer representative for the wave number to take the partial sums 
    ! iK=1+nint(log(KR/(B*KminSI)) / DeltaLnK)             

    !Remember that  Kmin=cElectronCharge/(PInjectionSI*exp(real(nP)*DLogP)) 

    !Calculate D_{xx}: KRes-dependent part
    Dxx = (4.0*BSI**2*MomentumSI*cLightSpeed**2)           / &
         (momentum_to_energy(MomentumSI,NameParticle)*cMu) * &
         (AK_II(iP,iX)-BK_II(iP,iX)*(kRSI**2))

  end function Dxx
  !============================================================================
end Module SP_ModTurbulence
