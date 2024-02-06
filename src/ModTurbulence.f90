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
  real, allocatable :: IPlusSi_IX(:,:),IMinusSi_IX(:,:),ICSi_X(:)
  real, allocatable :: vAlfvenSi_I(:)
  real, allocatable :: kOverBSi_I(:)
  real, allocatable :: kSi_I(:)

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

    allocate(IPlusSi_IX(0:nP+1,1:nVertexMax), &
         IMinusSi_IX(0:nP+1,1:nVertexMax))
    allocate(kOverBSi_I(0:nP+1), kSi_I(0:nP+1))
    allocate(ICSi_X(1:nVertexMax),CorrectionMode_X(1:nVertexMax))
    allocate(vAlfvenSi_I(1:nVertexMax))
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
    deallocate(IPlusSi_IX,IMinusSi_IX)
    deallocate(kOverBSi_I, kSi_I)
    deallocate(ICSi_X,CorrectionMode_X)
    deallocate(vAlfvenSi_I)
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
  subroutine init_spectrum(iEnd, XyzSi_DI, BSi_I, MomentumSi_I, dLogP, iShock,&
       CoefInj, AlfvenMach)
    !==============Initial spectrum of turbulence=============================!
    ! We recover the initial spectrum of turbulence from the spatial
    ! distribution of the diffusion coefficient and its dependence on the
    ! particle energy.

    ! the number of active particles on the line
    integer, intent(in)::  iEnd

    ! Coordinates of Lagrangian Meshes in SI unit [m]
    real,dimension(1:3,1:iEnd),intent(in) :: XyzSi_DI

    ! Magnetic field intensity in SI unit [T]
    real,dimension(1:iEnd),intent(in) :: BSi_I

    ! momentum in SI unit
    real,intent(in)     :: MomentumSi_I(0:nP+1)

    ! delta log p in SI unit
    real,intent(in)     :: dLogP

    ! coef of injection
    real,intent(in)     :: CoefInj

    ! Alfven March number
    real,intent(in)     :: AlfvenMach

    ! shock index
    integer, intent(in) :: iShock

    integer :: iVertex,iK
    real    :: ICOldSi, kSi
    real    :: rSi , rShockSi

    !--------------------------------------------------------------------------

    IPlusSi_IX  = 0.0; IMinusSi_IX   = 0.0
    vAlfvenSi_I = 0.0; DispersionA_I = 0.0

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

    kOverBSi_I = cElectronCharge/MomentumSi_I

    rShockSi = sqrt(sum(XyzSi_DI(:,iShock)**2))

    call assign_kolmogorov_spectrum(1,iEnd, 0, nP+1, &
         XyzSi_DI(:, 1:iEnd), BSi_I(1:iEnd))

    CorrectionMode_X=1.0

    do iVertex=1,iEnd
       rSi   = sqrt(sum(XyzSi_DI(:,iVertex)**2))
       kSi_I = kOverBSi_I*BSi_I(iVertex)

       if (rSi < 1.1*rShockSi .and. .false.) then
          ! In this part of the spectrum another equation governs the diffusion
          ICOldSi    = ICSi_X(iVertex)

!!!!
          ICSi_X(iVertex) = CoefInj *                                 &
               10.0*BSi_I(iVertex)**2*max(AlfvenMach,2.0)/(cMu*3.0) * &
               min(1.0, rSi/rShockSi/(1.0-Alpha))

          if(ICSi_X(iVertex) > ICOldSi)then
             IPlusSi_IX(:,  iVertex) = &
                  (1.0 - Alpha)*ICSi_X(iVertex)/kSi_I
             IMinusSi_IX(:, iVertex) = &
                  Alpha/(1.0 - Alpha)*IPlusSi_IX(:,iVertex)

             CorrectionMode_X(iVertex)=2
          else
             ! Do not change if the newly calculated value is less than the
             ! old one
             ICSi_X(iVertex) = ICOldSi
          end if
       end if
    end do

    DoInitSpectrum = .false.

  end subroutine init_spectrum
  !============================================================================
  subroutine assign_kolmogorov_spectrum(iParticleFirst, iParticleLast, &
       iKFirst, iKLast, XyzSi_DI, BSi_I)

    integer, intent(in):: iParticleFirst, iParticleLast, iKFirst, iKLast
    real,    intent(in):: XyzSi_DI(1:3,iParticleFirst:iParticleLast)
    real,    intent(in):: BSi_I(iParticleFirst:iParticleLast)

    integer :: iVertex
    real    :: kSi_I(iKFirst:iKLast), kr0Si, rSi
    !--------------------------------------------------------------------------
    do iVertex=iParticleFirst,iParticleLast
       ! kr0 = c*e*B/(1 GeV) in SI unit.
       kr0Si = cElectronCharge*BSi_I(iVertex)*cLightSpeed/cGeV

       rSi = sqrt(sum(XyzSi_DI(:,iVertex)**2))

       ICSi_X(iVertex)=54*BSi_I(iVertex)**2/(7.0*cPi*cMu*Lambda0InAu*rSi) &
            /kr0Si**(1.0/3)

       ! kOverBSi_I = cElectronCharge/MomentumSi_I
       kSi_I = kOverBSi_I(iKFirst:iKLast)*BSi_I(iVertex)

       ! I_{+} + I_{-} = IC/k^(5/3)
       ! I_{+}(k)=(1-\Alpha)*IC/K^{5/3} and I_{-}=\Alpha*I_{+}/(1-\Alpha),
       ! where IC=54*B^2/(7*pi*mu0*Lambda_0*kr0^{1/3}*(r/1 [AU])) and
       ! \Lambda_0=0.4 [AU], and kr0 = c*e*B/(1 GeV), all in SI unit
       IPlusSi_IX(iKFirst:iKLast,iVertex)  =  &
            (1.0-Alpha)*ICSi_X(iVertex)/kSi_I**(5.0/3)
       IMinusSi_IX(iKFirst:iKLast,iVertex) =  &
            Alpha/(1.0-Alpha)*IPlusSi_IX(iKFirst:iKLast,iVertex)
    end do
  end subroutine assign_kolmogorov_spectrum
  !============================================================================
  subroutine set_wave_advection_rates(iEnd, BSi_I, BOldSi_I, RhoSi_I, &
       RhoOldSi_I, XyzSi_DI, DsSi_I, DLogP, Dt, DtReduction)
    !=======================Set advection rate in k-space======================
    integer, intent(in) :: iEnd
    real, intent(in)    :: BSi_I(1:iEnd), BOldSi_I(1:iEnd)
    real, intent(in)    :: RhoSi_I(1:iEnd), RhoOldSi_I(1:iEnd)
    real, intent(in)    :: XyzSi_DI(1:3,1:iEnd)
    real, intent(in)    :: DsSi_I(1:iEnd)
    real, intent(in)    :: DLogP, Dt
    real, intent(out)   :: DtReduction

    ! local vars
    integer :: iVertex, DsSi, DLogRho, DLogB
    !--------------------------------------------------------------------------

    ! Calculate alfven speed in SI unit
    vAlfvenSi_I = BSi_I/sqrt(cMu*RhoSi_I)

    RhoCompression_I=1.50*log(RhoSi_I/RhoOldSi_I) - log(BSi_I/BOldSi_I)

    ! Contribution to the advection in k_space from the Lagrangian derivatives:
    ! Dispersion = Dt * [ (D ln {\rho})/(Dt) - 2 (D ln { B })/(D t) ] + ...

    DispersionA_I = log(RhoSi_I*BOldSi_I**2/(RhoOldSi_I*BSi_I**2))

    if(UseAdvectionWithAlfvenSpeed)then
       do iVertex =1,iEnd
          ! In this case, only first order accuracy between 2 - iEnd-1
          ! based on L323 in ModGrid.
          if (iVertex /= iEnd) then
             DsSi = DsSi_I(iVertex)
          else
             ! Seems DsSi_I(iEnd) is not defined.
             DsSi = DsSi_I(iEnd-1)
          end if

          CFL_I(iVertex) = Dt*vAlfvenSi_I(iVertex)/DsSi
       end do

       ! Now add the contribution from the spacial derivative:
       ! ... + Dt * [ \partial ((1/2) ln {rho})/(\partial s)-2\partial( ln{B})/(\partial s) ]
       ! for the "plus" wave

       DLogRho = log(RhoSi_I(2)/ RhoSi_I(1))
       DLogB   = log(BSi_I(2)  / BSi_I(1))
       DispersionPlus_I(1)  =                                       &
            DispersionA_I(1)+CFL_I(1 )*(0.5*DLogRho-2.0*DLogB)

       DLogRho = log(RhoSi_I(iEnd)/ RhoSi_I(iEnd-1))
       DLogB   = log(BSi_I(iEnd)  / BSi_I(iEnd-1))
       DispersionPlus_I(iEnd) =                                     &
            DispersionA_I(iEnd)+CFL_I(iEnd)*(0.5*DLogRho-2.0*DLogB)

       do iVertex=2,iEnd-1
          DLogRho = log(RhoSi_I(iVertex+1)/ RhoSi_I(iVertex-1))/2
          DLogB   = log(BSi_I(iVertex+1)  / BSi_I(iVertex-1))/2

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
  subroutine set_dxx(iEnd, nP, BSi_I)
    integer,intent(in) :: iEnd,nP
    real,intent(in)    :: BSi_I(iEnd)

    integer:: iVertex,iK
    real :: F01,F02,F11,F12
    real :: k0Si,k1Si
    real :: SpectralIndexAtKMax, ISumSi

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

       kSi_I = BSi_I(iVertex)*kOverBSi_I

       ! Initially from kMax
       k0Si = kSi_I(1)

       ! The sum of I_{plus}+I_{minus} at P_max
       ISumSi = IPlusSi_IX(1,iVertex)+IMinusSi_IX(1,iVertex)

       ! The integrand for AK_I, BK_I
       F01 = 1.0/(k0Si**2)/ISumSi
       F02 = 1.0/(k0Si**4)/ISumSi

       ! As the starting values for AK_I and BK_I at the minimum momentum,
       ! solve the integrals from K_{max} up to \infty, assuming the
       ! power law spectrum of turbulence at K>K_{max}
       AK_II(1,iVertex)=F01/(2.0-SpectralIndexAtKMax)
       BK_II(1,iVertex)=F02/(4.0-SpectralIndexAtKMax)

       if (iVertex == iParticleTest .and. DoTestMe) then
          write(*,*) 'SpectralIndex   =', SpectralIndexAtKMax
          write(*,*) 'kOverBSi_I(1)   =', kOverBSi_I(1)
          write(*,*) 'BSi_I(iVertex)=', BSi_I(iVertex)
          write(*,*) 'k0Si            =', k0Si
          write(*,*) 'ISumSi          =', ISumSi
          write(*,*) 'expected ISum   =', 54*BSi_I(iVertex)**2 &
               /(7*cPi*cMu*Lambda0InAu*9.3286125374064124E+08) &
               *(cGev/cLightSpeed/cElectronCharge/BSi_I(iVertex))**(1./3)
          write(*,*) 'AK_II(1,iVertex), BK_II(1,iVertex) =', &
               AK_II(1,iVertex), BK_II(1,iVertex)
       end if

       do iK=2,nP
          ! We calculate the partial sums for a set of the wave number values.
          ! The integral is taken from KRes up to infinity, so we start from
          ! the maximal wave number and add the contributions from each of
          ! the wave number intervals.

          ! The current value of the wave number
          k1Si = kSi_I(iK)

          ! The sum of I_{plus}+I_{minus} at P
          ISumSi = IPlusSi_IX(iK,iVertex)+IMinusSi_IX(iK,iVertex)

          ! The integrands at the lower value of the wave number
          F11 = 1.0/(k1Si**2)/ ISumSi
          F12 = 1.0/(k1Si**4)/ ISumSi

          ! Calculate the new partial sums

          AK_II(iK,iVertex)=AK_II(iK-1,iVertex)+0.5*(F01+F11)*dLogK
          BK_II(iK,iVertex)=BK_II(iK-1,iVertex)+0.5*(F02+F12)*dLogK

          ! current values saved as the initial values for the next step in
          ! the loop
          k0Si = k1Si; F01=F11; F02=F12
       end do
    end do
  end subroutine set_dxx
  !============================================================================
  subroutine update_spectrum(iEnd, nP, MomentumSi_I, DLogP, &
       XyzSi_DI, DsSi_I, F_II, BSi_I, RhoSi_I, SP_Dt)

    ! This is the subroutine which advances the wave spectrum through
    ! a time step, by solving the equations
    ! dI_+/dt=\gamma_+ I_+
    ! dI_-/dt=\gamma_- I_-

    use SP_ModAdvanceAdvection, ONLY: advance_log_advection
    use ModLinearAdvection

    ! The number of points and the number of the energy intervals
    integer,intent(in)::iEnd,nP

    real, intent(in) :: MomentumSi_I(0:nP+1), DLogP, SP_Dt
    ! Coordinates of the Lagrangian points
    real, intent(in) :: XyzSi_DI(3, 1:iEnd)
    real, intent(in) :: DsSi_I(1:iEnd)
    ! The distribution function
    real, intent(in) :: F_II(0:nP+1,1:iEnd)
    ! The magnetic field intensity and the particle number density in SI
    real, intent(in) :: BSi_I(1:iEnd), RhoSi_I(1:iEnd)

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
    real:: DsPlusSi,DsMinusSi

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
          call advance_log_advection(DispersionPlus_I(iVertex), 1, 1, &
               IPlusSi_IX( :,iVertex),IsConservative=.true., &
               DeltaLnP=DLogP)
          call advance_log_advection(DispersionMinus_I(iVertex), 1, 1, &
               IMinusSi_IX(:,iVertex),IsConservative=.true., &
               DeltaLnP=DLogP)
       else
          call advance_log_advection(DispersionA_I(iVertex), 0, 0, &
               IPlusSi_IX( 1:nP,iVertex),IsConservative=.true., &
               DeltaLnP=DLogP)
          call advance_log_advection(DispersionA_I(iVertex), 0, 0, &
               IMinusSi_IX(1:nP,iVertex),IsConservative=.true., &
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
          DsPlusSi = max(cTiny, DsSi_I(1))

          ! Use the forward spatial derivative
          DfDs0=(F_II(nP,iVertex+1)-F_II(nP,iVertex))/DsPlusSi

       else if (iVertex==iEnd) then
          ! again, DsSi_I(iEnd) is not defined...
          DsMinusSi = max(cTiny, DsSi_I(iEnd-1))

          ! Use the backward spatial derivative
          DfDs0=(F_II(nP,iVertex)-F_II(nP,iVertex-1))/DsMinusSi

       else
          DsPlusSi  = max(cTiny, DsSi_I(iVertex))
          DsMinusSi = max(cTiny, DsSi_I(iVertex-1))

          ! Use the average between the forward and backward spatial derivatives
          DfDs0=0.5*((F_II(nP,iVertex+1)-F_II(nP,iVertex  ))/DsPlusSi+&
               (F_II(nP,iVertex  )-F_II(nP,iVertex-1))/DsMinusSi)

       end if

       ! The momentum at the maximal energy

       P0 = MomentumSi_I(nP+1)

       ! We calculate the partial sums for a set of the momentum values.
       ! The integral is taken from pRes up to infinity, so we start from the
       ! maximal energy and add the contributions from each of the energy
       ! intervals.

       do iP=nP,1,-1

          ! Calculate the momentum at the lower energy
          P1 = MomentumSi_I(iP)

          ! Calculate the distribution function gradient at the lower energy

          if (iVertex == 1) then
             DfDs1 = (F_II(iP,iVertex+1)-F_II(iP,iVertex  ))/DsPlusSi
          else if (iVertex == iEnd) then
             DfDs1 = (F_II(iP,iVertex  )-F_II(iP,iVertex-1))/DsMinusSi
          else
             DfDs1 = 0.5* &
                  ( (F_II(iP,iVertex+1)-F_II(iP,iVertex  ))/DsPlusSi  +&
                  (  F_II(iP,iVertex  )-F_II(iP,iVertex-1))/DsMinusSi )
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
            iVertex, iVertex, 0, 0, XyzSi_DI(:,iVertex:iVertex),&
            BSi_I(iVertex:iVertex))
       ! IPlusSi_IX(    0,iVertex) =
       !       IPlusSi_IX(  0,iVertex)*ExpRhoCompression
       ! IMinusSi_IX(   0,iVertex) =
       !       IMinusSi_IX( 0,iVertex)*ExpRhoCompression

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
          K=BSi_I(iVertex)*kOverBSi_I(iK)

          ! The resonant momentum, for a given K
          PRes = cElectronCharge*BSi_I(iVertex)/K

          iP = iK

          ! here is the dynamic \gamma evaluated
          Gamma=-4.0*2.0*(cPi**2)*vAlfvenSi_I(iVertex)/K*&
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
          C1=(IPlusSi_IX(iK,iVertex)*IMinusSi_IX(iK,iVertex)) &
               *ExpRhoCompression**2
          C2=(IPlusSi_IX(iK,iVertex)-IMinusSi_IX(iK,iVertex)) &
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

             IPlusSi_IX(iK,iVertex)=(C3+sqrt(C3**2+4.0*C1))/2.0

             ! Now we take again into account the conservation of the
             ! product I_+*I_-
             ! which is the consequence of the relationship \gamma_-=-\gamma_+

             IMinusSi_IX(iK,iVertex)=C1/IPlusSi_IX(iK,iVertex)
          else

             IMinusSi_IX(iK,iVertex)=(-C3+sqrt(C3**2+4.0*C1))/2.0

             ! Now we take again into account the conservation of the product
             ! I_+*I_-
             ! which is the consequence of the relationship \gamma_-=-\gamma_+

             IPlusSi_IX(iK,iVertex)=C1/IMinusSi_IX(iK,iVertex)
          end if

          if(IPlusSi_IX(iK,iVertex)<=0.0)then
             write(*,*)'IPlus(iK,iVertex) < 0, IPlus,iK,iVertex=',&
                  IPlusSi_IX(iK,iVertex),iK,iVertex,c1,c2,c3,gamma,SP_Dt
             stop
          end if

          if(IMinusSi_IX(iK,iVertex)<=0.0)then
             write(*,*)'IMinus(iK,iVertex) < 0, IMinus,iK,iVertex=',&
                  IMinusSi_IX(iK,iVertex),iK,iVertex,c1,c2,c3,gamma,SP_Dt
             stop
          end if

       end do    ! cycling iK
       IPlusSi_IX( nP+1,iVertex) = IPlusSi_IX( nP,iVertex)
       IMinusSi_IX(nP+1,iVertex) = IMinusSi_IX(nP,iVertex)
    end do       ! cycling iVertex

    if(UseAdvectionWithAlfvenSpeed)then
       do iK=0,nP+1
          IPlusSi_IX(iK,1:iEnd) = &
               vAlfvenSi_I(1:iEnd)*IPlusSi_IX(iK,1:iEnd)/BSi_I( 1:iEnd)
          call advance_lin_advection_plus(CFL_I( 1:iEnd),&
               iEnd,0,0,IPlusSi_IX(iK,1:iEnd))

          IPlusSi_IX(iK,1:iEnd) = &
               BSi_I( 1:iEnd)*IPlusSi_IX(iK,1:iEnd)/vAlfvenSi_I(1:iEnd)

          IMinusSi_IX(iK,1:iEnd)= &
               vAlfvenSi_I(1:iEnd)*IMinusSi_IX(iK,1:iEnd)/BSi_I( 1:iEnd)

          call advance_lin_advection_minus(CFL_I(1:iEnd),&
               iEnd,0,0,IMinusSi_IX(iK,1:iEnd))

          IMinusSi_IX(iK,1:iEnd) = &
               BSi_I( 1:iEnd)*IMinusSi_IX(iK,1:iEnd)/vAlfvenSi_I(1:iEnd)
       end do
    end if
  end subroutine update_spectrum
  !============================================================================
  real function Dxx(iX, iP, MomentumSi, SpeedSi, BSi)
    integer,intent(in) :: iX, iP
    real,   intent(in) :: MomentumSi, SpeedSi, BSi

    real    :: kRSi

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
    kRSi = cElectronCharge*BSi/MomentumSi

    ! Calculate D_{xx}: KRes-dependent part
    Dxx = BSi**2*SpeedSi/(cMu*cPi) * (AK_II(iP,iX)-BK_II(iP,iX)*kRSi**2)

    if (iX == iParticleTest .and. iP == iPTest .and. DoTestMe) then
       write(*,*) 'AK_II(iP,iX), BK_II(iP,iX) =', &
            AK_II(iP,iX), BK_II(iP,iX)
       write(*,*) 'Dxx =', Dxx
       write(*,*) 'D from Li =', 1./3*Lambda0InAu*9.3286125374064124E+08 &
            *(MomentumSi*cLightSpeed/cGeV)**(1./3)*SpeedSi
    end if

  end function Dxx
  !============================================================================
end module SP_ModTurbulence
!==============================================================================
