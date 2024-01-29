!  Copyright (C) 2002 Regents of the University of Michigan,
!  portions used
!  with permission
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module SP_ModTurbulence
  use ModConst
  use SP_ModDistribution, ONLY: nP
  use SP_ModGrid,         ONLY: iPTest, iParticleTest
  implicit none
  SAVE

  private

  public :: init, finalize, DoInitSpectrum, UseTurbulentSpectrum, set_dxx, &
       read_param, set_wave_advection_rates, reduce_advection_rates, dxx,  &
       init_spectrum, update_spectrum

  ! Logicals, all .false. by default
  logical:: DoInitSpectrum              = .false.
  logical:: UseTurbulentSpectrum        = .false.
  logical:: UseAdvectionWithAlfvenSpeed = .false.

  integer, parameter :: nK = nP
  real    :: dLogK

  real, allocatable :: Gamma_I(:,:)
  real, allocatable :: IPlusSI_IX(:,:),IMinusSI_IX(:,:),ICSI_X(:)
  real, allocatable :: vAlfvenSI_I(:)
  real, allocatable :: kOverBSI_I(:)
  real, allocatable :: kSI_I(:)

  ! Rate of advection in k space, neglecting
  ! a spacial advection with the Alfven speed
  real, allocatable:: DispersionA_I(:)

  ! Rate of advection in k space, for I_+/I_- wave
  real, allocatable:: DispersionPlus_I(:)
  real, allocatable:: DispersionMinus_I(:)

  ! This is the ratio of densities powered 3/2
  real, allocatable:: RhoCompression_I(:)

  !------------------------------------------------------------------------!
  !          Grid in the momentum space                                    !
  ! iP     0     1                         nP   nP+1                        !
  !       |     |    ....                 |     |                          !
  ! P      P_inj P_inj*exp(\Delta (Ln P))  P_Max P_Max*exp(\Delta (Ln P))   !
  !             |    Grid in k-space      |     |                          !
  ! K/B         KMax                      KMin                              !
  ! ik     0     1                         nP   nP+1                         !
  !------------------------------------------------------------------------!

  real,allocatable,private:: AK_II(:,:)
  real,allocatable,private:: BK_II(:,:)

  real,allocatable,private:: CFL_I(:)

  integer,allocatable::      CorrectionMode_X(:)

  ! the intensity of the back travelling wave in the initial condition
  real :: Alpha       = 1.0/10
  real :: Lambda0InAu = 4.0/10.0  ![AU]
contains
  !============================================================================
  subroutine read_param(NameCommand)
    use ModReadParam, ONLY: read_var
    character(len=*), intent(in):: NameCommand

    character(len=*), parameter:: NameSub = 'read_param'
    !--------------------------------------------------------------------------

    select case(NameCommand)
    case('#INITSPECTRUM')
       call read_var('DoInitSpectrum', DoInitSpectrum)
    case('#TURBULENTSPECTRUM')
       call read_var('UseTurbulentSpectrum', UseTurbulentSpectrum)
       if (UseTurbulentSpectrum) then
          DoInitSpectrum = .true.
          call read_var('Alpha',       Alpha)
          call read_var('Lambda0InAu', Lambda0InAu)
       end if
    case('#ADVECTIONWITHALFVENSPEED')
       call read_var('UseAdvectionWithAlfvenSpeed', &
            UseAdvectionWithAlfvenSpeed)
    case default
       call CON_stop(NameSub//' Unknown command '//NameCommand)
    end select

    if (UseTurbulentSpectrum) UseAdvectionWithAlfvenSpeed = .true.
  end subroutine read_param
  !============================================================================
  subroutine init
    ! Init all the allocatable vars
    use SP_ModSize, ONLY: nVertexMax
    use SP_ModProc, ONLY: iProc
    character(len=*), parameter:: NameSub = 'init'
    !--------------------------------------------------------------------------

    allocate(IPlusSI_IX(0:nP+1,1:nVertexMax), &
         IMinusSI_IX(0:nP+1,1:nVertexMax))
    allocate(kOverBSI_I(0:nP+1), kSI_I(0:nP+1))
    allocate(ICSI_X(1:nVertexMax),CorrectionMode_X(1:nVertexMax))
    allocate(vAlfvenSI_I(1:nVertexMax))
    allocate(DispersionA_I(1:nVertexMax))

    allocate(RhoCompression_I(1:nVertexMax))

    ! if(DoOutputGamma) &
    !     allocate(Gamma_I(nP,(iXOutputLast-iXOutputStart)/iXOutputStride+1))

    allocate(DispersionPlus_I(1:nVertexMax))
    allocate(DispersionMinus_I(1:nVertexMax))
    allocate(CFL_I(1:nVertexMax))

    allocate(AK_II(nP,nVertexMax),BK_II(nP,nVertexMax))

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

    character(len=*), parameter:: NameSub = 'finalize'
    !--------------------------------------------------------------------------
    deallocate(IPlusSI_IX,IMinusSI_IX)
    deallocate(kOverBSI_I, kSI_I)
    deallocate(ICSI_X,CorrectionMode_X)
    deallocate(vAlfvenSI_I)
    deallocate(DispersionA_I)

    deallocate(RhoCompression_I)

    ! if(DoOutputGamma) &
    !     allocate(Gamma_I(nP,(iXOutputLast-iXOutputStart)/iXOutputStride+1))

    deallocate(DispersionPlus_I)
    deallocate(DispersionMinus_I)
    deallocate(CFL_I)

    deallocate(AK_II,BK_II)

  end subroutine finalize
  !============================================================================
  subroutine init_spectrum(iEnd, XyzSI_DI, BSI_I, MomentumSI_I, dLogP, iShock,&
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

    ! delta log p in SI unit
    real,intent(in)     :: dLogP

    ! coef of injection
    real,intent(in)     :: CoefInj

    ! Alfven March number
    real,intent(in)     :: AlfvenMach

    ! shock index
    integer, intent(in) :: iShock

    integer :: iVertex,iK
    real    :: ICOldSI, kSI
    real    :: rSI , rShockSI

    !--------------------------------------------------------------------------

    IPlusSI_IX  = 0.0; IMinusSI_IX   = 0.0
    vAlfvenSI_I = 0.0; DispersionA_I = 0.0

    RhoCompression_I = 0.0

    if(UseAdvectionWithAlfvenSpeed)then
       DispersionPlus_I  = 0.0
       DispersionMinus_I = 0.0
       CFL_I = 1.0
    end if

    !-------------------------------------------------------------------------!
    !          Grid in the momentum space                                     !
    ! iP     0     1                         nP   nP+1                         !
    !       |     |    ....                 |     |                           !
    ! P      P_inj P_inj*exp(\Delta (Ln P))  P_Max P_Max*exp(\Delta (Ln P))    !
    !             |    Grid in k-space      |     |                           !
    ! K/B         KMax                      KMin                               !
    ! ik     0     1                         nP   nP+1                         !
    !-------------------------------------------------------------------------!

    ! k = e*B/p => dlog(k) = - dlog(p)
    dLogK = dLogP

    kOverBSI_I = cElectronCharge/MomentumSI_I

    rShockSI = sqrt(sum(XyzSI_DI(:,iShock)**2))

    call assign_kolmogorov_spectrum(1,iEnd, 0, nP+1, &
         XyzSI_DI(:, 1:iEnd), BSI_I(1:iEnd))

    CorrectionMode_X=1.0

    do iVertex=1,iEnd
       rSI   = sqrt(sum(XyzSI_DI(:,iVertex)**2))
       kSI_I = kOverBSI_I*BSI_I(iVertex)

       if (rSI < 1.1*rShockSI .and. .false.) then
          ! In this part of the spectrum another equation governs the diffusion
          ICOldSI    = ICSI_X(iVertex)

          !!!!
          ICSI_X(iVertex) = CoefInj *                                 &
               10.0*BSI_I(iVertex)**2*max(AlfvenMach,2.0)/(cMu*3.0) * &
               min(1.0, rSI/rShockSI/(1.0-Alpha))

          if(ICSI_X(iVertex) > ICOldSI)then
             IPlusSI_IX(:,  iVertex) = &
                  (1.0 - Alpha)*ICSI_X(iVertex)/kSI_I
             IMinusSI_IX(:, iVertex) = &
                  Alpha/(1.0 - Alpha)*IPlusSI_IX(:,iVertex)

             CorrectionMode_X(iVertex)=2
          else
             ! Do not change if the newly calculated value is less than the
             ! old one
             ICSI_X(iVertex) = ICOldSI
          end if
       end if
    end do

    DoInitSpectrum = .false.

  end subroutine init_spectrum
  !============================================================================
  subroutine assign_kolmogorov_spectrum(iParticleFirst, iParticleLast, &
       iKFirst, iKLast, XyzSI_DI, BSI_I)

    integer, intent(in):: iParticleFirst, iParticleLast, iKFirst, iKLast
    real,    intent(in):: XyzSI_DI(1:3,iParticleFirst:iParticleLast)
    real,    intent(in):: BSI_I(iParticleFirst:iParticleLast)

    integer :: iVertex
    real    :: kSI_I(iKFirst:iKLast), kr0SI, rSI
    !--------------------------------------------------------------------------
    do iVertex=iParticleFirst,iParticleLast
       ! kr0 = c*e*B/(1 GeV) in SI unit.
       kr0SI = cElectronCharge*BSI_I(iVertex)*cLightSpeed/cGeV

       rSI = sqrt(sum(XyzSI_DI(:,iVertex)**2))

       ICSI_X(iVertex)=54*BSI_I(iVertex)**2/(7.0*cPi*cMu*Lambda0InAu*rSi) &
            /kr0SI**(1.0/3)

       ! kOverBSI_I = cElectronCharge/MomentumSI_I
       kSI_I = kOverBSI_I(iKFirst:iKLast)*BSI_I(iVertex)

       ! I_{+} + I_{-} = IC/k^(5/3)
       ! I_{+}(k)=(1-\Alpha)*IC/K^{5/3} and I_{-}=\Alpha*I_{+}/(1-\Alpha),
       ! where IC=54*B^2/(7*pi*mu0*Lambda_0*kr0^{1/3}*(r/1 [AU])) and
       ! \Lambda_0=0.4 [AU], and kr0 = c*e*B/(1 GeV), all in SI unit
       IPlusSI_IX(iKFirst:iKLast,iVertex)  =  &
            (1.0-Alpha)*ICSI_X(iVertex)/kSI_I**(5.0/3)
       IMinusSI_IX(iKFirst:iKLast,iVertex) =  &
            Alpha/(1.0-Alpha)*IPlusSI_IX(iKFirst:iKLast,iVertex)
    end do
  end subroutine assign_kolmogorov_spectrum
  !============================================================================
  subroutine set_wave_advection_rates(iEnd, BSI_I, BOldSI_I, RhoSI_I, &
    RhoOldSI_I, XyzSI_DI, DsSI_I, DLogP, Dt, DtReduction)
    !=======================Set advection rate in k-space======================
    integer, intent(in) :: iEnd
    real, intent(in)    :: BSI_I(1:iEnd), BOldSI_I(1:iEnd)
    real, intent(in)    :: RhoSI_I(1:iEnd), RhoOldSI_I(1:iEnd)
    real, intent(in)    :: XyzSI_DI(1:3,1:iEnd)
    real, intent(in)    :: DsSI_I(1:iEnd)
    real, intent(in)    :: DLogP, Dt
    real, intent(out)   :: DtReduction

    ! local vars
    integer :: iVertex, DsSI, DLogRho, DLogB
    !--------------------------------------------------------------------------

    ! Calculate alfven speed in SI unit
    vAlfvenSI_I = BSI_I/sqrt(cMu*RhoSI_I)

    RhoCompression_I=1.50*log(RhoSI_I/RhoOldSI_I) - log(BSI_I/BOldSI_I)

    ! Contribution to the advection in k_space from the Lagrangian derivatives:
    ! Dispersion = Dt * [ (D ln {\rho})/(Dt) - 2 (D ln { B })/(D t) ] + ...

    DispersionA_I = log(RhoSI_I*BOldSI_I**2/(RhoOldSI_I*BSI_I**2))

    if(UseAdvectionWithAlfvenSpeed)then
       do iVertex =1,iEnd
          ! In this case, only first order accuracy between 2 - iEnd-1
          ! based on L323 in ModGrid.
          if (iVertex /= iEnd) then
             DsSI = DsSI_I(iVertex)
          else
             ! Seems DsSI_I(iEnd) is not defined.
             DsSI = DsSI_I(iEnd-1)
          end if

          CFL_I(iVertex) = Dt*vAlfvenSI_I(iVertex)/DsSI
       end do

       ! Now add the contribution from the spacial derivative:
       ! ... + Dt * [ \partial ((1/2) ln {rho})/(\partial s)-2\partial( ln{B})/(\partial s) ]
       ! for the "plus" wave

       DLogRho = log(RhoSI_I(2)/ RhoSI_I(1))
       DLogB   = log(BSI_I(2)  / BSI_I(1))
       DispersionPlus_I(1)  =                                       &
            DispersionA_I(1)+CFL_I(1 )*(0.5*DLogRho-2.0*DLogB)

       DLogRho = log(RhoSI_I(iEnd)/ RhoSI_I(iEnd-1))
       DLogB   = log(BSI_I(iEnd)  / BSI_I(iEnd-1))
       DispersionPlus_I(iEnd) =                                     &
            DispersionA_I(iEnd)+CFL_I(iEnd)*(0.5*DLogRho-2.0*DLogB)

       do iVertex=2,iEnd-1
          DLogRho = log(RhoSI_I(iVertex+1)/ RhoSI_I(iVertex-1))/2
          DLogB   = log(BSI_I(iVertex+1)  / BSI_I(iVertex-1))/2

          DispersionPlus_I(iVertex) = DispersionA_I(iVertex)    &
               + CFL_I(iVertex)*(0.5*DLogRho-2.0*DLogB)
       end do

       ! Now add the contribution from the spacial derivative:
       ! ... - Dt * [ \partial ((1/2) ln {rho})/(\partial s)-2\partial( ln{B})/(\partial s) ]
       ! for the "minus" wave

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
    !--------------------------------------------------------------------------
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
    integer,intent(in) :: iEnd,nP
    real,intent(in)    :: BSI_I(iEnd)

    integer:: iVertex,iK
    real :: F01,F02,F11,F12
    real :: k0SI,k1SI
    real :: SpectralIndexAtKMax, ISumSI

    logical :: DoTestMe = .false.

    ! The coefficient of a spatial diffusion along the magnetic field is
    ! given by a formula in the SI unit:
    ! D_{xx}=v*B^2/(cPi*cMu)*(AK_I(k_{res})-k_{res}^2*BK_I(k_{res})),

    ! where AK_I=\int_{k_{res}}^\infty{d(\log k)/(k^2*(I_{+}(k) + I_{-}(k) )) }
    ! and   BK_I=\int_{k_{res}}^\infty{d(\log k)/(k^4*(I_{+}(k) + I_{-}(k) )) }

    !          Grid in the momentum space                                     !
    ! iP     0     1                         nP   nP+1                         !
    !       |     |    ....                 |     |                           !
    ! P      P_inj P_inj*exp(\Delta (Ln P))  P_Max P_Max*exp(\Delta (Ln P))    !
    !             |    Grid in k-space      |     |                           !
    ! K/B         kMax(k0)                  kMin                               !
    ! ik     0     1                         nP   nP+1                         !
    !--------------------------------------------------------------------------

    do iVertex=1,iEnd
       select case(CorrectionMode_X(iVertex))
       case(1)
          SpectralIndexAtKMax = 5.0/3
       case(2)
          SpectralIndexAtKMax = 1.0
       end select

       kSI_I = BSI_I(iVertex)*kOverBSI_I

       ! Initially from kMax
       k0SI = kSI_I(1)

       ! The sum of I_{plus}+I_{minus} at P_max
       ISumSI = IPlusSI_IX(1,iVertex)+IMinusSI_IX(1,iVertex)

       ! The integrand for AK_I, BK_I
       F01 = 1.0/(k0SI**2)/ISumSI
       F02 = 1.0/(k0SI**4)/ISumSI

       ! As the starting values for AK_I and BK_I at the minimum momentum,
       ! solve the integrals from K_{max} up to \infty, assuming the
       ! power law spectrum of turbulence at K>K_{max}
       AK_II(1,iVertex)=F01/(2.0-SpectralIndexAtKMax)
       BK_II(1,iVertex)=F02/(4.0-SpectralIndexAtKMax)

       if (iVertex == iParticleTest .and. DoTestMe) then
          write(*,*) 'SpectralIndex   =', SpectralIndexAtKMax
          write(*,*) 'kOverBSI_I(1)   =', kOverBSI_I(1)
          write(*,*) 'BSI_I(iVertex)=', BSI_I(iVertex)
          write(*,*) 'k0SI            =', k0SI
          write(*,*) 'ISumSI          =', ISumSI
          write(*,*) 'expected ISum   =', 54*BSI_I(iVertex)**2 &
               /(7*cPi*cMu*Lambda0InAu*9.3286125374064124E+08) &
               *(cGev/cLightSpeed/cElectronCharge/BSI_I(iVertex))**(1./3)
          write(*,*) 'AK_II(1,iVertex), BK_II(1,iVertex) =', &
               AK_II(1,iVertex), BK_II(1,iVertex)
       end if

       do iK=2,nP
          ! We calculate the partial sums for a set of the wave number values.
          ! The integral is taken from KRes up to infinity, so we start from
          ! the maximal wave number and add the contributions from each of
          ! the wave number intervals.

          ! The current value of the wave number
          k1SI = kSI_I(iK)

          ! The sum of I_{plus}+I_{minus} at P
          ISumSI = IPlusSI_IX(iK,iVertex)+IMinusSI_IX(iK,iVertex)

          ! The integrands at the lower value of the wave number
          F11 = 1.0/(k1SI**2)/ ISumSI
          F12 = 1.0/(k1SI**4)/ ISumSI

          ! Calculate the new partial sums

          AK_II(iK,iVertex)=AK_II(iK-1,iVertex)+0.5*(F01+F11)*dLogK
          BK_II(iK,iVertex)=BK_II(iK-1,iVertex)+0.5*(F02+F12)*dLogK

          ! current values saved as the initial values for the next step in
          ! the loop
          k0SI = k1SI; F01=F11; F02=F12
       end do
    end do
  end subroutine set_dxx
  !============================================================================
  subroutine update_spectrum(iEnd, nP, MomentumSI_I, DLogP, &
       XyzSI_DI, DsSI_I, F_II, BSI_I, RhoSI_I, SP_Dt)

    ! This is the subroutine which advances the wave spectrum through
    ! a time step, by solving the equations
    ! dI_+/dt=\gamma_+ I_+
    ! dI_-/dt=\gamma_- I_-

    use SP_ModAdvection, ONLY: advance_log_advection
    use ModLinearAdvection

    ! The number of points and the number of the energy intervals
    integer,intent(in)::iEnd,nP

    real, intent(in) :: MomentumSI_I(0:nP+1), DLogP, SP_Dt
    ! Coordinates of the Lagrangian points
    real, intent(in) :: XyzSI_DI(3, 1:iEnd)
    real, intent(in) :: DsSI_I(1:iEnd)
    ! The distribution function
    real, intent(in) :: F_II(0:nP+1,1:iEnd)
    ! The magnetic field intensity and the particle number density in SI
    real, intent(in) :: BSI_I(1:iEnd), RhoSI_I(1:iEnd)

    integer :: iVertex,iK,iP

    ! Resonant value of the particle momentum, for a given K
    real:: PRes

    ! The increment for the I_+ wave, multiplied by (I_+  +  I_-). Below
    ! the dynamical equations for I_+ and I_- are reformulated in terms of
    ! \gamma and integrated exactly assuming the constant value of \gamma
    ! (NOT the increment) through the time step
    real:: Gamma

    ! Wave number
    real:: K

    ! Forward and backward values for Ds
    real:: DsPlusSI,DsMinusSI

    ! Spatial derivatives of the distribution function at given values
    ! of the momentum
    real:: DfDs0,P0,DfDs1,P1

    ! Partial sums in the integral for Gamma
    real:: A_I(nP+1),B_I(nP+1)

    ! Variables used in the spectral model for turbulence
    real:: ExpRhoCompression

    ! Miscellaneous
    real:: C1,C2,C3

    !--------------------------------------------------------------------------
    do iVertex=1,iEnd
       ! Advection in k space:
       if(UseAdvectionWithAlfvenSpeed)then
          call advance_log_advection(DispersionPlus_I(iVertex),  nP, 1, 1, &
               IPlusSI_IX( :,iVertex),IsConservative=.true., &
               DeltaLnP=DLogP)
          call advance_log_advection(DispersionMinus_I(iVertex), nP, 1, 1, &
               IMinusSI_IX(:,iVertex),IsConservative=.true., &
               DeltaLnP=DLogP)
       else
          call advance_log_advection(DispersionA_I(iVertex), nP, 0, 0, &
               IPlusSI_IX( 1:nP,iVertex),IsConservative=.true., &
               DeltaLnP=DLogP)
          call advance_log_advection(DispersionA_I(iVertex), nP, 0, 0, &
               IMinusSI_IX(1:nP,iVertex),IsConservative=.true., &
               DeltaLnP=DLogP)
       end if

       !----------------------------------------------------------------------!
       !          Grid in the momentum space                                  !
       ! iP     0     1                         nP   nP+1                      !
       !       |     |    ....                 |     |                        !
       ! P      P_inj P_inj*exp(\Delta (Ln P))  P_Max P_Max*exp(\Delta (Ln P)) !
       !             |    Grid in k-space      |     |                        !
       ! K/B         KMax                      KMin                            !
       ! ik     0     1                         nP   nP+1                      !
       !----------------------------------------------------------------------!

       ! Calculate the partial sums in the integral for \Gamma
       ! Nullify the sums

       A_I(nP+1)=0.0; B_I(nP+1)=0.0

       ! Calculate the spatial derivatives of the distribution function
       ! at the maximal energy

       if (iVertex==1) then
          DsPlusSI = max(cTiny, DsSI_I(1))

          ! Use the forward spatial derivative
          DfDs0=(F_II(nP,iVertex+1)-F_II(nP,iVertex))/DsPlusSI

       else if (iVertex==iEnd) then
          ! again, DsSI_I(iEnd) is not defined...
          DsMinusSI = max(cTiny, DsSI_I(iEnd-1))

          ! Use the backward spatial derivative
          DfDs0=(F_II(nP,iVertex)-F_II(nP,iVertex-1))/DsMinusSI

       else
          DsPlusSI  = max(cTiny, DsSI_I(iVertex))
          DsMinusSI = max(cTiny, DsSI_I(iVertex-1))

          ! Use the average between the forward and backward spatial derivatives
          DfDs0=0.5*((F_II(nP,iVertex+1)-F_II(nP,iVertex  ))/DsPlusSI+&
               (F_II(nP,iVertex  )-F_II(nP,iVertex-1))/DsMinusSI)

       end if

       ! The momentum at the maximal energy

       P0 = MomentumSI_I(nP+1)

       ! We calculate the partial sums for a set of the momentum values.
       ! The integral is taken from pRes up to infinity, so we start from the
       ! maximal energy and add the contributions from each of the energy
       ! intervals.

       do iP=nP,1,-1

          ! Calculate the momentum at the lower energy
          P1 = MomentumSI_I(iP)

          ! Calculate the distribution function gradient at the lower energy

          if (iVertex == 1) then
             DfDs1 = (F_II(iP,iVertex+1)-F_II(iP,iVertex  ))/DsPlusSI
          else if (iVertex == iEnd) then
             DfDs1 = (F_II(iP,iVertex  )-F_II(iP,iVertex-1))/DsMinusSI
          else
             DfDs1 = 0.5* &
                  ( (F_II(iP,iVertex+1)-F_II(iP,iVertex  ))/DsPlusSI  +&
                  (  F_II(iP,iVertex  )-F_II(iP,iVertex-1))/DsMinusSI )
          end if

          ! here are the parts for \gamma being integrated
          ! For I_+ the increment equals \frac{\gamma}{(I_{+}+I_{-})}
          ! where \gamma =-4.0*(cPi**2)*Va/k *Integral/cProtonMass
          ! Integral=\int_{p_{res}(k)}^{\infty}{\frac{dp}{p}*p^4
          ! (p_{res}(k)-\frac{p_{res}^3}{p^2})\frac{\partial{f}}{\partial{s}}
          !
          ! A=p_{res}^{k}=\int{d(ln{p})*(p^{4}*{\frac{\partial{f}}{\partial{s}}}}
          ! B=p_{res}^{3}\int{d(ln{p})*p^2*{\frac{\partial{f}}{\partial{s}}

          A_I(iP)=A_I(iP+1)+0.5*(DfDs0*(P0**4)+DfDs1*(P1**4))*DLogP
          B_I(iP)=B_I(iP+1)+0.5*(DfDs0*(P0**2)+DfDs1*(P1**2))*DLogP

          ! Save the values for the lower energy to re-use them as the values
          ! for the higher energy end of the next interval (in CYCLING the
          ! momentummo DECREASES)

          P0=P1; DfDs0=DfDs1
       end do

       ! Calculate the wave increment and update the wave spectra
       ExpRhoCompression=exp(RhoCompression_I(iVertex))

       call assign_kolmogorov_spectrum( &
            iVertex, iVertex, 0, 0, XyzSI_DI(:,iVertex:iVertex),&
            BSI_I(iVertex:iVertex))
       ! IPlusSI_IX(    0,iVertex) =
       !       IPlusSI_IX(  0,iVertex)*ExpRhoCompression
       ! IMinusSI_IX(   0,iVertex) =
       !       IMinusSI_IX( 0,iVertex)*ExpRhoCompression

       !----------------------------------------------------------------------!
       !          Grid in the momentum space                                  !
       ! iP     0     1                         nP   nP+1                      !
       !       |     |    ....                 |     |                        !
       ! P      P_inj P_inj*exp(\Delta (Ln P))  P_Max P_Max*exp(\Delta (Ln P)) !
       !             |    Grid in k-space      |     |                        !
       ! K/B         KMax                      KMin                            !
       ! ik     0     1                         nP   nP+1                      !
       !----------------------------------------------------------------------!

       do iK = nK,1,-1
          ! The wave number
          K=BSI_I(iVertex)*kOverBSI_I(iK)

          ! The resonant momentum, for a given K
          PRes = cElectronCharge*BSI_I(iVertex)/K

          iP = iK

          ! here is the dynamic \gamma evaluated
          Gamma=-4.0*2.0*(cPi**2)*vAlfvenSI_I(iVertex)/K*&
               (PRes*A_I(iP)-(PRes**3)*B_I(iP))/    &
               cProtonMass

          ! if(i_output(iVertex)/=0)Gamma_I(iK,i_output(iVertex))=Gamma

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
          C1=(IPlusSI_IX(iK,iVertex)*IMinusSI_IX(iK,iVertex)) &
               *ExpRhoCompression**2
          C2=(IPlusSI_IX(iK,iVertex)-IMinusSI_IX(iK,iVertex)) &
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

             IPlusSI_IX(iK,iVertex)=(C3+sqrt(C3**2+4.0*C1))/2.0

             ! Now we take again into account the conservation of the
             ! product I_+*I_-
             ! which is the consequence of the relationship \gamma_-=-\gamma_+

             IMinusSI_IX(iK,iVertex)=C1/IPlusSI_IX(iK,iVertex)
          else

             IMinusSI_IX(iK,iVertex)=(-C3+sqrt(C3**2+4.0*C1))/2.0

             ! Now we take again into account the conservation of the product
             ! I_+*I_-
             ! which is the consequence of the relationship \gamma_-=-\gamma_+

             IPlusSI_IX(iK,iVertex)=C1/IMinusSI_IX(iK,iVertex)
          end if

          if(IPlusSI_IX(iK,iVertex)<=0.0)then
             write(*,*)'IPlus(iK,iVertex) < 0, IPlus,iK,iVertex=',&
                  IPlusSI_IX(iK,iVertex),iK,iVertex,c1,c2,c3,gamma,SP_Dt
             stop
          end if

          if(IMinusSI_IX(iK,iVertex)<=0.0)then
             write(*,*)'IMinus(iK,iVertex) < 0, IMinus,iK,iVertex=',&
                  IMinusSI_IX(iK,iVertex),iK,iVertex,c1,c2,c3,gamma,SP_Dt
             stop
          end if

       end do    ! cycling iK
       IPlusSI_IX( nP+1,iVertex) = IPlusSI_IX( nP,iVertex)
       IMinusSI_IX(nP+1,iVertex) = IMinusSI_IX(nP,iVertex)
    end do       ! cycling iVertex

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
  real function Dxx(iX, iP, MomentumSI, SpeedSI, BSI)
    integer,intent(in) :: iX, iP
    real,   intent(in) :: MomentumSI, SpeedSI, BSI

    real    :: kRSI

    logical :: DoTestMe =.false.

    !          Grid in the momentum space                                     !
    ! iP     0     1                         nP   nP+1                         !
    !       |     |    ....                 |     |                           !
    ! P      P_inj P_inj*exp(\Delta (Ln P))  P_Max P_Max*exp(\Delta (Ln P))    !
    !             |    Grid in k-space      |     |                           !
    ! K/B         KMax                      KMin                               !
    ! ik     0     1                         nP   nP+1                         !
    !--------------------------------------------------------------------------

    ! The coefficient of a spatial diffusion along the magnetic field is
    ! given by a formula in the SI unit:
    ! D_{xx}=v*B^2/(cPi*cMu)*(AK_I(k_{res})-k_{res}^2*BK_I(k_{res})),

    ! where AK_I=\int_{k_{res}}^\infty{d(\log k)/(k^2*(I_{+}(k) + I_{-}(k) )) }
    ! and   BK_I=\int_{k_{res}}^\infty{d(\log k)/(k^4*(I_{+}(k) + I_{-}(k) )) }

    ! The resonant wave number, kr = e*B/p in the SI unit
    kRSI = cElectronCharge*BSI/MomentumSI

    ! Calculate D_{xx}: KRes-dependent part
    Dxx = BSI**2*SpeedSI/(cMu*cPi) * (AK_II(iP,iX)-BK_II(iP,iX)*kRSI**2)

    if (iX == iParticleTest .and. iP == iPTest .and. DoTestMe) then
       write(*,*) 'AK_II(iP,iX), BK_II(iP,iX) =', &
            AK_II(iP,iX), BK_II(iP,iX)
       write(*,*) 'Dxx =', Dxx
       write(*,*) 'D from Li =', 1./3*Lambda0InAu*9.3286125374064124E+08 &
            *(MomentumSI*cLightSpeed/cGeV)**(1./3)*SpeedSI
    end if

  end function Dxx
  !============================================================================
end module SP_ModTurbulence
!==============================================================================
