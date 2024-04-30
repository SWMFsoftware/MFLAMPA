!  Copyright (C) 2002 Regents of the University of Michigan,
!  portions used with permission
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module SP_ModDistribution

  ! The module contains the velocity/momentum distribution function
  ! and methods for initializing it as well as the offset routine.

#ifdef OPENACC
  use ModUtilities, ONLY: norm2
#endif
  use ModUtilities, ONLY: CON_stop
  use ModNumConst, ONLY: cTiny
  use ModConst,    ONLY: cLightSpeed, energy_in, cMeV
  use SP_ModSize,  ONLY: nVertexMax, nP => nMomentum, &
       nMu => nPitchAngle, IsMuAvg => IsPitchAngleAverage
  use SP_ModUnit,  ONLY: NameFluxUnit, NameEnergyFluxUnit,&
       Io2Si_V, Si2Io_V, NameFluxUnit_I, UnitEnergy_, UnitFlux_, &
       kinetic_energy_to_momentum, momentum_to_energy
  use SP_ModGrid,  ONLY: nLine, nVertex_B

  implicit none

  SAVE

  private ! except

  ! Public members:
  public:: init              ! Initialize Distribution_IIB
  public:: read_param        ! Read momentum grid parameters
  public:: offset            ! Sync. index in State_VIB and Dist_IIB
  public:: get_integral_flux ! Calculate Flux_VIB
  public:: nP                ! Number of points in the momentum grid
  public:: nMu               ! Number of points over pitch-angle (\mu)
  public:: IsMuAvg           ! If .true., dist. function is omnidirectional

  ! Injection and maximal energy, in kev (or, in Io energy unit)
  ! To be read from the PARAM.in file
  real, public :: EnergyInjIo = 10.0, EnergyMaxIo = 1.0E+07

  ! Injection momentum, mimimum momentum value in Si
  real, public :: MomentumInjSi

  ! Size of a  log-momentum mesh. For momentum we use both the
  ! denotaion, P, and a word, momentum - whichever is more covenient
  real, public :: dLogP      ! log(MomentumMaxSi/MomentumInjSi)/nP

  ! uniform distribution of pitch angle
  real, public :: DeltaMu = 2.0/nMu
  real, public :: MuFace_I(0:nMu), Mu_I(nMu)

  ! speed, momentum, kinetic energy at the momentum grid points
  real, public, dimension(0:nP+1) :: SpeedSi_I, Momentum_I, &
       KinEnergy_I, VolumeP_I, Background_I
  real, public :: Momentum3_I(-1:nP+1) ! P**3/3, normalized per MomentumInjSi

  ! Total integral (simulated) particle flux
  integer, parameter, public :: Flux0_ = 0
  integer, parameter, public :: FluxFirst_ =1  ! The first channel
  integer, public :: nFluxChannel = 6          ! GOES by default, 6 channels
  integer, public :: FluxLast_ = 6             ! The last channel
  integer, public :: EFlux_    = 7             ! Total integral energy flux
  integer, public :: FluxMax_  = 7
  real, allocatable :: EChannelIo_I(:) ! energy limits of the instrument
  real, public, allocatable  :: Flux_VIB(:,:,:)
  character(len=10), public, allocatable :: NameFluxChannel_I(:)

  !-----------------Grid in the momentum space---------------------------------
  ! iP     0     1                         nP   nP+1
  !        |     |    ....                 |     |
  ! P      P_inj P_inj*exp(\Delta(log P))  P_Max P_Max*exp(DLogP)
  !-----------------Control volumes and faces----------------------------------
  ! Vol  | 0  |  1  |.........          |  nP | nP+1 |
  ! Fcs -1    0     1 ....            nP-1    nP    nP+1
  !----------------------------------------------------------------------------
  ! This is because we put two boundary conditions: the background
  ! value at the right one and the physical condition at the left
  ! one, for the velocity distribution function

  ! Velosity Distribution Function (VDF)
  ! Number of points along the momentum axis is set in ModSize
  ! 1st index - log(momentum)
  ! 2nd index - particle index along the field line
  ! 3rd index - local line number
  real, public, allocatable :: Distribution_CB(:,:,:,:)

  ! distribution is initialized to have integral flux:
  real:: FluxInitIo = 0.01 ! [PFU]
  ! initial values of fluxes in energy channels
  real, public, allocatable:: FluxChannelInit_V(:)

  logical :: DoInit = .true.
contains
  !============================================================================
  subroutine init

    use SP_ModUnit,   ONLY: momentum_to_kinetic_energy
    use ModUtilities, ONLY: check_allocate
    ! loop variables
    integer:: iLine, iVertex, iP, iError, iMu
    ! maximal and current momenta
    real :: MomentumMaxSi, MomentumSi
    ! set the initial distribution on all lines
    character(len=*), parameter:: NameSub = 'init'
    !--------------------------------------------------------------------------

    if(.not.DoInit)RETURN
    DoInit = .false.

    ! convert energies to momenta
    MomentumInjSi= kinetic_energy_to_momentum(EnergyInjIo*Io2Si_V(UnitEnergy_))
    MomentumMaxSi= kinetic_energy_to_momentum(EnergyMaxIo*Io2Si_V(UnitEnergy_))

    ! grid size in the log momentum space
    dLogP = log(MomentumMaxSi/MomentumInjSi)/nP

    ! Functions to convert the grid index to momentum and energy
    Momentum3_I(-1) = exp(-3*(0.50*dLogP))/3 ! P^3/3 at -0.5*dLogP from PInj
    do iP = 0, nP+1
       Momentum_I(iP)    = MomentumInjSi*exp(iP*dLogP)
       SpeedSi_I(iP)     = Momentum_I(iP)*cLightSpeed**2/ &
            momentum_to_energy(Momentum_I(iP))
       Momentum3_I(iP)   = Momentum3_I(iP-1)*exp(3*dLogP)
       VolumeP_I(iP)     = Momentum3_I(iP) - Momentum3_I(iP-1)
       ! Normalize kinetic energy per Unit of energy in SI unit:
       KinEnergy_I(iP)   = momentum_to_kinetic_energy(Momentum_I(iP)) &
            *Si2Io_V(UnitEnergy_)
       ! Normalize momentum per MomentumInjSi
       Momentum_I(iP)    = Momentum_I(iP)/MomentumInjSi
       Background_I(iP)  = FluxInitIo*Io2Si_V(UnitFlux_)/ & ! Integral flux SI
            (EnergyMaxIo-EnergyInjIo) &  ! Energy range
            /Momentum_I(iP)**2           ! Convert from diff flux to VDF
    end do

    ! Calculate all the mu values at cell center and faces
    do iMu = 0, nMu
       MuFace_I(iMu) = -1.0 + real(iMu)*DeltaMu
    end do
    Mu_I = (MuFace_I(0:nMu-1) + MuFace_I(1:nMu))*0.5

    ! Distribution function
    allocate(Distribution_CB(0:nP+1,nMu,nVertexMax,nLine), stat=iError)
    call check_allocate(iError, 'Distribution_IIB')

    ! initialization depends on momentum, however, this corresponds
    ! to a constant differential flux (intensity), thus ensuring
    ! uniform backgound while visualizing this quantity
    do iLine = 1, nLine
       do iVertex = 1, nVertexMax
          do iMu = 1, nMu
             ! Overall density of the fast particles is of the order
             ! of 10^-6 m^-3. Integral flux is less than 100 per
             ! (m^2 ster s). Differential background flux is constant.
             Distribution_CB(:,iMu,iVertex,iLine) = Background_I
          end do
       end do
    end do

    ! GOES by default
    if (.not. allocated(NameFluxChannel_I)) then
       nFluxChannel = 6
       allocate(NameFluxChannel_I(0:nFluxChannel+1))
       NameFluxChannel_I = ['flux_total', 'flux_00005', 'flux_00010', &
            'flux_00030', 'flux_00050', 'flux_00060', 'flux_00100',   &
            'eflux     ']
       allocate (EChannelIo_I(nFluxChannel))
       EChannelIo_I = [5,10,30,50,60,100]  ! in MeV!
    end if
    EChannelIo_I  = EChannelIo_I & ! In MeV Now!!
         *cMeV                   & ! in SI
         *Si2Io_V(UnitEnergy_)     ! in NameUnitEnergy
    if(allocated(NameFluxUnit_I))&
         deallocate(NameFluxUnit_I)
    allocate(NameFluxUnit_I(0:nFluxChannel+1))
    NameFluxUnit_I(0:nFluxChannel) = NameFluxUnit
    NameFluxUnit_I(1+nFluxChannel) = NameEnergyFluxUnit

    if (.not. allocated(Flux_VIB)) then
       allocate(Flux_VIB(Flux0_:FluxMax_,1:nVertexMax,nLine), &
            stat=iError); call check_allocate(iError, 'Flux_VIB')
       Flux_VIB = -1.0
    else
       call CON_stop(NameSub//' Flux_VIB already allocated')
    end if

    ! fill initial values of flux in energy channels
    allocate(FluxChannelInit_V(0:nFluxChannel+1))
    ! for the assumed initial distribution (~1/p^2)
    FluxChannelInit_V(0) = FluxInitIo
    FluxChannelInit_V(1:nFluxChannel) = FluxInitIo * &
         (EnergyMaxIo - EChannelIo_I(:)) / (EnergyMaxIo - EnergyInjIo)
    FluxChannelInit_V(1+nFluxChannel) = FluxInitIo * Io2Si_V(UnitFlux_) * &
         0.5 * (EnergyMaxIo + EnergyInjIo) * Si2Io_V(UnitFlux_)
  end subroutine init
  !============================================================================
  subroutine read_param(NameCommand)

    use ModReadParam, ONLY: read_var
    use SP_ModProc,   ONLY: iProc
    character(len=*), intent(in):: NameCommand ! From PARAM.in
    integer:: nPCheck = nP, nMuCheck = nMu, iFluxChannel
    character(len=5) :: NameFluxChannel
    character(len=*), parameter:: NameSub = 'read_param'
    !--------------------------------------------------------------------------
    select case(NameCommand)
    case('#MOMENTUMGRID')
       ! Read unit to be used for particle energy: eV, keV, GeV
       call read_var('EnergyMin',      EnergyInjIo)
       call read_var('EnergyMax',      EnergyMaxIo)
       call read_var('nP',             nPCheck    )
       if(nP/=nPCheck)then
          if(iProc==0) write(*,'(a,i6,a,i6)') NameSub//' '//     &
               'Code is configured with nMomentum=', nP,         &
               ' while value read from PARAM.in is nP=', nPCheck
          call CON_stop('Modify PARAM.in or reconfigure SP/MFLAMPA')
       end if
       call read_var('nMu', nMuCheck)
       if(nMu/=nMuCheck)then
          if(iProc==0) write(*,'(a,i6,a,i6)') NameSub//' '//     &
               'Code is configured with nMu=', nMu,              &
               ' while value read from PARAM.in is nMu=', nMuCheck
          call CON_stop('Modify PARAM.in or reconfigure SP/MFLAMPA')
       end if
    case('#FLUXINITIAL')
       call read_var('FluxInit', FluxInitIo)
       ! check correctness
       if(FluxInitIo<=0)call CON_stop(NameSub//': flux value must be positive')
    case('#FLUXCHANNEL')
       call read_var('nFluxChannel', nFluxChannel)
       FluxLast_ = nFluxChannel
       EFlux_    = FluxLast_ + 1
       FluxMax_  = EFlux_

       if (allocated(EChannelIo_I)) deallocate(EChannelIo_I)
       allocate(EChannelIo_I(nFluxChannel))
       if (allocated(NameFluxChannel_I)) deallocate(NameFluxChannel_I)
       allocate(NameFluxChannel_I(0:FluxMax_))
       if(allocated(NameFluxUnit_I)) deallocate(NameFluxUnit_I)
       allocate(NameFluxUnit_I(0:FluxMax_))

       NameFluxChannel_I(0)              = 'flux_total'
       NameFluxChannel_I(nFluxChannel+1) = 'eflux'

       do iFluxChannel=1,nFluxChannel
          call read_var('EChannelIo_I [MeV]', EChannelIo_I(iFluxChannel))
          write(NameFluxChannel,'(I5.5)') int(EChannelIo_I(iFluxChannel))
          NameFluxChannel_I(iFluxChannel) = 'flux_'//NameFluxChannel
       end do
       ! Bug: here, Si2No_V is not yet calculated
       ! EChannelIo_I = EChannelIo_I & ! in MeV
       !     *cMeV                  & ! in SI
       !     *Si2Io_V(UnitEnergy_)
       ! NameFluxUnit_I(0:nFluxChannel) = NameFluxUnit
       ! NameFluxUnit_I(EFlux_)         = NameEnergyFluxUnit

    case default
       call CON_stop(NameSub//' Unknown command '//NameCommand)
    end select

  end subroutine read_param
  !============================================================================
  subroutine offset(iLine, iOffset)

    use SP_ModGrid, ONLY: NoShock_, BOld_, RhoOld_, ShockOld_, &
         iShock_IB,  State_VIB, MHData_VIB, X_, Z_, FootPoint_VB
    ! shift in the data arrays is required if the grid point(s) is
    ! appended or removed at the foot point of the magnetic field
    ! line. SHIFTED ARE: State_VIB(/RhoOld_,BOld_),Distribution_IIB
    ! as well as ShockOld_
    integer, intent(in)        :: iLine
    integer, intent(in)        :: iOffset
    real :: Alpha, Distance2ToMin, Distance3To2
    character(len=*), parameter:: NameSub = 'offset'
    !--------------------------------------------------------------------------

    if(iOffset==0)RETURN
    if(iOffset==1)then
       State_VIB([RhoOld_,BOld_],2:nVertex_B(iLine),iLine) &
            = State_VIB([RhoOld_,BOld_],1:nVertex_B(iLine)-1,iLine)
       Distribution_CB(:,:,2:nVertex_B(iLine), iLine)&
            = Distribution_CB(:,:,1:nVertex_B(iLine)-1, iLine)
       ! Extrapolate state vector components and VDF at iVertex=1
       Distance2ToMin = norm2(MHData_VIB(X_:Z_,2,iLine) - &
            FootPoint_VB(X_:Z_,iLine))
       Distance3To2   = norm2(MHData_VIB(X_:Z_,3,iLine) - &
            MHData_VIB(X_:Z_,2,iLine))
       Alpha = Distance2ToMin/(Distance2ToMin + Distance3To2)
       State_VIB([RhoOld_, BOld_], 1, iLine) = &
            (Alpha + 1)*State_VIB([RhoOld_, BOld_], 2, iLine) &
            -Alpha     * State_VIB([RhoOld_, BOld_], 3, iLine)
       Distribution_CB(:,:,1,iLine) = Distribution_CB(:,:,2,iLine) + &
            Alpha*(Distribution_CB(:,:,2,iLine) - &
            Distribution_CB(:,:,3,iLine))
       ! extrapolation may introduced negative values
       ! for strictly positive quantities; such occurences need fixing
       where(State_VIB([RhoOld_,BOld_],1,iLine) <= 0.0)
          State_VIB([RhoOld_,BOld_],1,iLine) = &
               0.01 * State_VIB([RhoOld_,BOld_],2,iLine)
       end where
       where(Distribution_CB(:,:,1,iLine) <= 0.0)
          Distribution_CB(:,:,1,iLine) = &
               0.01 * Distribution_CB(:,:,2,iLine)
       end where
    elseif(iOffset < 0)then
       State_VIB([RhoOld_,BOld_],1:nVertex_B(iLine),iLine) &
            =  State_VIB([RhoOld_,BOld_],1-iOffset:nVertex_B(iLine)&
            - iOffset, iLine)
       Distribution_CB(:,:,1:nVertex_B(iLine), iLine)&
            = Distribution_CB(:,:,1-iOffset:nVertex_B(iLine)-iOffset, &
            iLine)
    else
       call CON_stop('No algorithm for iOffset > 1 in '//NameSub)
    end if
    if(iShock_IB(ShockOld_, iLine)/=NoShock_)&
         iShock_IB(ShockOld_, iLine) = &
         max(iShock_IB(ShockOld_, iLine) + iOffset, 1)
  end subroutine offset
  !============================================================================
  subroutine get_integral_flux

    use ModConst, ONLY: energy_in
    use SP_ModGrid,  ONLY: Used_B
    ! compute the total (simulated) integral flux of particles as well as
    ! particle flux in the 6 GOES channels; also compute total energy flux

    integer:: iLine, iVertex, iP, iFlux ! loop variables
    real   :: EFlux ! the value of energy flux
    real   :: dFlux, dFlux1 ! increments
    real   :: Flux          ! the value of particle flux
    real   :: Flux_I(nFluxChannel)
    !--------------------------------------------------------------------------
    ! energy limits of GOES channels

    do iLine = 1, nLine
       if(.not.Used_B(iLine))CYCLE
       do iVertex = 1, nVertex_B( iLine)
          ! Integration loop with midpoint rule
          ! reset values
          EFlux = 0.0
          Flux_I= 0.0
          Flux  = 0.0
          do iP = 1, nP - 1
             ! the flux increment from iP
             dFlux = 0.5 * &
                  (KinEnergy_I(iP+1) - KinEnergy_I(iP)) * (&
                  Distribution_CB(iP,nMu,  iVertex,iLine)*&
                  Momentum_I(iP)**2 &
                  +&
                  Distribution_CB(iP+1,nMu,iVertex,iLine)*&
                  Momentum_I(iP+1)**2)

             ! increase the total flux
             Flux = Flux + dFlux

             ! increase GOES channels' fluxes
             do iFlux = 1, nFluxChannel
                ! check whether reached the channel's cut-off level
                if(KinEnergy_I(iP+1) < EChannelIo_I(iFlux))&
                     CYCLE

                if(KinEnergy_I(iP) >= EChannelIo_I(iFlux))then
                   Flux_I(iFlux) = Flux_I(iFlux) + dFlux
                else
                   ! channel cutoff level is often in the middle of a bin;
                   ! compute partial flux increment
                   dFlux1 =&
                        ((-0.50*(KinEnergy_I(iP) + EChannelIo_I(iFlux)) + &
                        KinEnergy_I(iP+1) )*&
                        Distribution_CB(iP,nMu,iVertex,iLine)*&
                        Momentum_I(iP)**2  &
                        -0.50*(KinEnergy_I(iP)-EChannelIo_I(iFlux))*&
                        Distribution_CB(iP+1,nMu,iVertex,iLine)*&
                        Momentum_I(iP+1)**2)*&
                        (KinEnergy_I(iP)-EChannelIo_I(iFlux))/&
                        (KinEnergy_I(iP+1)-KinEnergy_I(iP))
                   Flux_I(iFlux) = Flux_I(iFlux) + dFlux1
                end if
             end do

             ! increase total energy flux
             EFlux = EFlux + 0.5 * &
                  (KinEnergy_I(iP+1) - KinEnergy_I(iP)) * (&
                  Distribution_CB(iP,nMu,iVertex,iLine)*&
                  KinEnergy_I(iP) * &
                  Momentum_I(iP)**2 &
                  +&
                  Distribution_CB(iP+1,nMu,iVertex,iLine)*&
                  KinEnergy_I(iP+1) * &
                  Momentum_I(iP+1)**2)
          end do

          ! store the results
          Flux_VIB(Flux0_,              iVertex, iLine) = &
               Flux   * Si2Io_V(UnitFlux_)
          Flux_VIB(FluxFirst_:FluxLast_,iVertex, iLine) = &
               Flux_I * Si2Io_V(UnitFlux_)
          Flux_VIB(EFlux_,              iVertex, iLine) = &
               EFlux  * Si2Io_V(UnitFlux_)
       end do
    end do
  end subroutine get_integral_flux
  !============================================================================
end module SP_ModDistribution
!==============================================================================
