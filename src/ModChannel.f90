!  Copyright (C) 2002 Regents of the University of Michigan,
!  portions used with permission
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module SP_ModChannel

  ! Setting energy channels and get the integral/differential intensity

  use ModUtilities, ONLY: CON_stop
  use SP_ModGrid,   ONLY: nLine, nVertex_B, Used_B
  use SP_ModUnit,   ONLY: Io2Si_V, Si2Io_V, UnitEnergy_, UnitFlux_, &
       NameFluxUnit, NameDiffFluxUnit, NameEnergyFluxUnit

  implicit none

  save

  private ! Except

  public:: init               ! initialize the flux channels
  public:: read_param         ! read satellite-related parameters
  public:: get_integral_flux  ! calculate Flux_VIB

  ! ------------ For saving the intensity in energy channels ------------
  integer, public :: nFluxChannelSat = 1  ! GOES satellite as default
  integer, public :: nFluxChannel    = 6  ! GOES as default, 6 channels
  integer, parameter :: LenNameSat   = 12 ! satellite name length
  type FluxChannelSat
     ! Full set of information, for the energy channels of each satellite
     !
     ! ------ General information ------
     ! Satellite name
     character(len=LenNameSat):: NameSat
     ! Flux channel counts
     integer:: nFluxChannel
     ! Satellite index
     integer:: iKindSat
     !
     ! ------ Saved energy channels and flux type ------
     ! Satellite flux type
     integer:: iKindFlux
     ! Energy channels: Lo & Hi bounds, first in the unit of MeV
     real, allocatable, dimension(:):: EChannelLoIo_I, EChannelHiIo_I
  end type FluxChannelSat

  ! All plot files
  type(FluxChannelSat), public, allocatable :: FluxChannelSat_I(:)
  integer, parameter :: &
       FluxChannelSELF_ = 0, & ! Self-defined
       FluxChannelGOES_ = 1, & ! GOES
       FluxChannelERNE_ = 2, & ! SOHO/ERNE
       FluxChannelEPAM_ = 3    ! ACE/EPAM
  integer, parameter :: &
       FluxInt_  = 1, &         ! Integral intensity
       FluxDiff_ = 2            ! Differential intensity

  ! Total integral (simulated) particle flux
  integer, parameter, public :: Flux0_     = 0 ! Total particle flux
  integer, parameter, public :: FluxFirst_ = 1 ! The first channel
  integer,            public :: FluxLast_  = 6 ! The last channel
  integer,            public :: EFlux_     = 7 ! Total integral energy flux
  integer,            public :: FluxMax_   = 7 ! Maximum flux channel index

  ! Flux saved for this satellite
  ! 1st index: Energy channel
  ! 2nd index: iVertex along the field line
  ! 3rd index: Field line number
  real, public, allocatable :: Flux_VIB(:,:,:)
  ! Initial values of fluxes in energy channels
  real, public, allocatable :: FluxChannelInit_V(:)

  ! Energy channels of the instrument, or defined by users
  integer, public, allocatable, dimension(:) :: EChannelType_I
  real,    public, allocatable, dimension(:) :: EChannelLoIo_I, EChannelHiIo_I
  ! Energy channel bin, for get_integral_flux
  real,            allocatable, dimension(:) :: InvEChannelBinIo_I
  ! Energy channel names and units
  character(len=:),  public, allocatable, dimension(:) :: NameFluxChannel_I
  character(len=12), public, allocatable, dimension(:) :: NameChannelSource_I
  character(len=20), public, allocatable, dimension(:) :: NameFluxUnit_I

  ! Logical variable: whether this is the initial call
  logical :: DoInit = .true.
contains
  !============================================================================
  subroutine read_param(NameCommand)

    use ModReadParam, ONLY: read_var
    use ModUtilities, ONLY: split_string, lower_case, upper_case
    character(len=*), intent(in) :: NameCommand
    ! general loop variables
    integer           :: iFlux
    ! input characters from EChannelIo_I in #ECHANNEL
    integer, parameter:: nStringEChannel = 3
    character(len=30) :: StringEChannel
    character(len=10) :: StringEChannel_I(nStringEChannel)
    ! input characters from EChannelIo_I in #ECHANNELSAT
    ! loop variable for satellites
    integer           :: iSat
    ! local variable of satellite names
    character(len=LenNameSat) :: NameSat
    character(len=*), parameter:: NameSub = 'read_param'
    !--------------------------------------------------------------------------

    select case(NameCommand)
    case('#ECHANNEL')
       call read_var('nFluxChannel', nFluxChannel)
       if(nFluxChannel < 0) then
          call CON_stop(NameSub // ': nFluxChannel must be non-negative')
       else if(nFluxChannel == 0) then
          RETURN
       else
          ! Do nothing
       end if
       FluxLast_ = nFluxChannel

       ! Allocate arrays for flux type, lower and upper bounds
       if(allocated(EChannelType_I)) deallocate(EChannelType_I)
       if(allocated(EChannelLoIo_I)) deallocate(EChannelLoIo_I)
       if(allocated(EChannelHiIo_I)) deallocate(EChannelHiIo_I)
       if(allocated(NameChannelSource_I)) deallocate(NameChannelSource_I)
       allocate(EChannelType_I(FluxFirst_:FluxLast_), &
            EChannelLoIo_I(FluxFirst_:FluxLast_), &
            EChannelHiIo_I(FluxFirst_:FluxLast_), &
            NameChannelSource_I(FluxFirst_:FluxLast_))
       ! Get the input EChannelIo_I in the unit of MeV
       do iFlux = FluxFirst_, FluxLast_
          call read_var('EChannelIo_I', StringEChannel)
          ! in the order of {int/diff} {low bound} {high bound}, in MeV!!
          call split_string(StringEChannel, StringEChannel_I)

          ! 1) Energy channel type: integral or differential intensity
          call lower_case(StringEChannel_I(1))
          select case(StringEChannel_I(1))
          case('int', 'inte', 'integral', 'integrate', 'integrated')
             EChannelType_I(iFlux) = FluxInt_
          case('diff', 'differential', 'spec', 'spectrum', 'spectra')
             EChannelType_I(iFlux) = FluxDiff_
          case default
             call CON_stop(NameSub // ": Unknown EChannelType, " // &
                  StringEChannel_I(1) // " at EChannel ", iFlux)
          end select

          ! 2) Low energy bound: do not know KinEnergyIo_I(1) (min value) here
          read(StringEChannel_I(2), *) EChannelLoIo_I(iFlux)
          if(EChannelLoIo_I(iFlux) < 0.0) call CON_stop(NameSub // &
               ": Energy channel lower bound should not be negative.")

          ! 3) High energy bound: do not know EnergyMaxIo (max value) here
          read(StringEChannel_I(3), *) EChannelHiIo_I(iFlux)
          if(EChannelHiIo_I(iFlux) < EChannelLoIo_I(iFlux)) &
               call CON_stop(NameSub // ": Energy channel upper bound" // &
               " should not be less than lower bound.")

          ! last but not least, name the flux channel source
          NameChannelSource_I(iFlux) = "Self-Defined"
       end do
    case('#ECHANNELSAT')
       ! Read the number of satellites, specified with the energy channels
       call read_var('nFluxChannelSat', nFluxChannelSat)
       ! This means we do not use any energy channel from any satellite
       if(nFluxChannelSat < 0) then
          call CON_stop(NameSub // ': nFluxChannelSat must be non-negative')
       else if(nFluxChannelSat == 0) then
          RETURN
       else
          ! Do nothing
       end if

       if(allocated(FluxChannelSat_I)) deallocate(FluxChannelSat_I)
       allocate(FluxChannelSat_I(1:nFluxChannelSat))

       ! Set KindSat and nFluxChannel for each specified satellite
       do iSat = 1, nFluxChannelSat
          call read_var('NameSatellite', NameSat)
          call upper_case(NameSat)
          FluxChannelSat_I(iSat) % NameSat = trim(NameSat)

          select case(trim(FluxChannelSat_I(iSat) % NameSat))
          case('GOES')
             FluxChannelSat_I(iSat) % iKindSat = FluxChannelGOES_
             FluxChannelSat_I(iSat) % nFluxChannel = 6
          case('ERNE', 'SOHO/ERNE')
             FluxChannelSat_I(iSat) % iKindSat = FluxChannelERNE_
             FluxChannelSat_I(iSat) % nFluxChannel = 20
          case('EPAM', 'ACE/EPAM')
             FluxChannelSat_I(iSat) % iKindSat = FluxChannelEPAM_
             FluxChannelSat_I(iSat) % nFluxChannel = 8
          case default
             FluxChannelSat_I(iSat) % iKindSat = -1
             FluxChannelSat_I(iSat) % nFluxChannel = 0
             call CON_stop(NameSub // ": satellite name (for " // &
                  "energy channels) isn't properly set in PARAM.in")
          end select
       end do
    end select

  end subroutine read_param
  !============================================================================
  subroutine init
    ! To initialize the energy channels of the satellite(s)

    use ModConst,           ONLY: cMeV
    use ModUtilities,       ONLY: check_allocate
    use SP_ModDistribution, ONLY: KinEnergyIo_I, &
         EnergyInjIo, EnergyMaxIo, FluxInitIo
    use SP_ModProc,         ONLY: iError
    use SP_ModSize,         ONLY: nVertexMax
    ! loop variables
    integer :: iSat, iFlux
    ! index VAR for convenience
    integer :: FluxLastSat_
    ! Tmp VARs
    integer :: nFluxChannelTmp = 0
    integer, allocatable, dimension(:) :: EChannelTypeTmp_I
    real,    allocatable, dimension(:) :: EChannelLoIoTmp_I, EChannelHiIoTmp_I
    character(len=20), allocatable, dimension(:) :: NameChannelSourceTmp_I
    ! local NameFluxChannel, saved the channel index
    character(len=2) :: NameFluxChannel

    character(len=*), parameter:: NameSub = 'init'
    !--------------------------------------------------------------------------
    if(.not.DoInit) RETURN
    DoInit = .false.

    ! Not specify satellite energy channels: GOES by default
    ! Not specify self-defined energy channel: 0
    ! 1) To ONLY use FluxChannelSat_I: not set extra nFluxChannel
    ! 2) To include self-defined channels: set nFluxChannel
    ! 3) To ONLY use self-defined channels: set nFluxChannel <= 0
    ! 4) To include both FluxChannelSat_I and self-defined channels: set both
    ! Note that here we allow for nFluxChannel==0 and nFluxChannelSat==0, and
    ! in this case, we only save Flux0_ and EFlux_ in the MH_data files.
    if(allocated(EChannelType_I)) then
       ! Self-defined channels: will drop satellite channels if not specified
       if(.not.allocated(FluxChannelSat_I)) nFluxChannelSat = 0
       ! else: we specify nFluxChannelSat > 0, summed later

       ! Convert units for self-defined energy channels, and the lower
       ! and upper bounds should fall between KinEnergyIo_I(1) and EnergyMaxIo
       EChannelLoIo_I = MIN(EnergyMaxIo, MAX(KinEnergyIo_I(1), &
            EChannelLoIo_I*cMeV*Si2Io_V(UnitEnergy_)))
       EChannelHiIo_I = MIN(EnergyMaxIo, MAX(KinEnergyIo_I(1), &
            EChannelHiIo_I*cMeV*Si2Io_V(UnitEnergy_)))
    else
       ! No self-defined channels
       nFluxChannel = 0
       ! In this case, nFluxChannelSat should be >= 0, either specified by
       ! users, or set to be default. For the former, FluxChannelSat_I is
       ! allocated as long as nFluxChannelSat>0, but when nFluxChannelSat==0,
       ! FluxChannelSat_I is not allocated, same as the latter situation.
       ! However, nFluxChannelSat is > 0, set as default for the latter case,
       ! so we can differentiate them by nFluxChannelSat.

       ! if no specify satellite energy channels: GOES by default
       ! To turn this off: set nFluxChannelSat = 0
       if(.not.allocated(FluxChannelSat_I) .and. nFluxChannelSat>0) then
          allocate(FluxChannelSat_I(1:nFluxChannelSat))
          FluxChannelSat_I(1) % NameSat = 'GOES'
          FluxChannelSat_I(1) % iKindSat = FluxChannelGOES_
          FluxChannelSat_I(1) % nFluxChannel = 6
       end if
    end if

    ! Save self-defined channels if there are plus nFluxChannelSat>0
    if(allocated(EChannelType_I) .and. nFluxChannelSat>0) then
       ! Save the self-defined flux channels
       nFluxChannelTmp = nFluxChannel
       allocate(NameChannelSourceTmp_I(FluxFirst_:nFluxChannelTmp), &
            EChannelTypeTmp_I(FluxFirst_:nFluxChannelTmp), &
            EChannelLoIoTmp_I(FluxFirst_:nFluxChannelTmp), &
            EChannelHiIoTmp_I(FluxFirst_:nFluxChannelTmp))
       NameChannelSourceTmp_I = NameChannelSource_I
       EChannelTypeTmp_I = EChannelType_I
       EChannelLoIoTmp_I = EChannelLoIo_I
       EChannelHiIoTmp_I = EChannelHiIo_I
       deallocate(NameChannelSource_I, EChannelType_I, &
            EChannelLoIo_I, EChannelHiIo_I)
    end if

    ! Next we initialize the details for nFluxChannelSat
    do iSat = 1, nFluxChannelSat
       ! Exclude unreasonable satellites
       if(FluxChannelSat_I(iSat)%iKindSat <= 0) call CON_stop(NameSub // &
            ": Incorrect satellite, " // FluxChannelSat_I(iSat)%NameSat // &
            ", set for energy channels.")

       ! Set the last channel index of satellite, according to its nFluxChannel
       FluxLastSat_ = FluxChannelSat_I(iSat)%nFluxChannel

       ! Allocate the energy channels, names, and units
       ! Lower bound: from 0 => physical KinEnergyIo_I(1) in code
       ! Upper bound: to +\infty => physical EnergyMaxIo in code
       if(allocated(FluxChannelSat_I(iSat)%EChannelLoIo_I)) &
            deallocate(FluxChannelSat_I(iSat)%EChannelLoIo_I)
       allocate(FluxChannelSat_I(iSat)%EChannelLoIo_I(FluxFirst_:FluxLastSat_))
       if(allocated(FluxChannelSat_I(iSat)%EChannelHiIo_I)) &
            deallocate(FluxChannelSat_I(iSat)%EChannelHiIo_I)
       allocate(FluxChannelSat_I(iSat)%EChannelHiIo_I(FluxFirst_:FluxLastSat_))

       ! Get the energy channels, names, and units (still in MeV now)
       select case(FluxChannelSat_I(iSat) % iKindSat)
       case(FluxChannelGOES_)
          FluxChannelSat_I(iSat)%iKindFlux = FluxInt_
          FluxChannelSat_I(iSat)%EChannelLoIo_I = [5, 10, 30, 50, 60, 100]
          FluxChannelSat_I(iSat)%EChannelHiIo_I = &
               EnergyMaxIo*Io2Si_V(UnitEnergy_)/cMeV
       case(FluxChannelERNE_)
          FluxChannelSat_I(iSat)%iKindFlux = FluxDiff_
          FluxChannelSat_I(iSat)%EChannelLoIo_I = [1.3, 1.6, 2.0, 2.5, &
               3.2, 4.0, 5.0, 6.4, 8.0, 10.0, 13.0, 16.0, 20.0, 25.0, &
               32.0, 40.0, 50.0, 64.0, 80.0, 100.0]
          FluxChannelSat_I(iSat)%EChannelHiIo_I(FluxFirst_:FluxLastSat_-1) = &
               FluxChannelSat_I(iSat)%EChannelLoIo_I(FluxFirst_+1:FluxLastSat_)
          FluxChannelSat_I(iSat)%EChannelHiIo_I(FluxLastSat_) = 130.0
       case(FluxChannelEPAM_)
          FluxChannelSat_I(iSat)%iKindFlux = FluxDiff_
          FluxChannelSat_I(iSat)%EChannelLoIo_I = &
               [0.046, 0.067, 0.115, 0.193, 0.315, 0.580, 1.06, 1.88]
          FluxChannelSat_I(iSat)%EChannelHiIo_I(FluxFirst_:FluxLastSat_-1) = &
               FluxChannelSat_I(iSat)%EChannelLoIo_I(FluxFirst_+1:FluxLastSat_)
          FluxChannelSat_I(iSat)%EChannelHiIo_I(FluxLastSat_) = 4.70
       end select

       ! Convert energy channel units: MeV to SI to Io, and the lower
       ! and upper bounds should fall between KinEnergyIo_I(1) and EnergyMaxIo
       FluxChannelSat_I(iSat)%EChannelLoIo_I = &
            MIN(EnergyMaxIo, MAX(KinEnergyIo_I(1), &
            FluxChannelSat_I(iSat)%EChannelLoIo_I*cMeV*Si2Io_V(UnitEnergy_)))
       FluxChannelSat_I(iSat)%EChannelHiIo_I = &
            MIN(EnergyMaxIo, MAX(KinEnergyIo_I(1), &
            FluxChannelSat_I(iSat)%EChannelHiIo_I*cMeV*Si2Io_V(UnitEnergy_)))

       ! Finally sum for nFluxChannel
       nFluxChannel = nFluxChannel + FluxChannelSat_I(iSat)%nFluxChannel
    end do

    ! With summed nFluxChannel, we get indices for flux channels
    FluxLast_ = nFluxChannel
    EFlux_    = FluxLast_ + 1
    FluxMax_  = EFlux_

    ! Re-assign energy channel info to the whole array when nFluxChannelSat>0
    if(nFluxChannelSat > 0) then
       allocate(NameChannelSource_I(FluxFirst_:FluxLast_), &
            EChannelType_I(FluxFirst_:FluxLast_), &
            EChannelLoIo_I(FluxFirst_:FluxLast_), &
            EChannelHiIo_I(FluxFirst_:FluxLast_))

       ! For the first a few channels: self-defined energy channels
       if(nFluxChannelTmp > 0) then
          iFlux = nFluxChannelTmp
          NameChannelSource_I(FluxFirst_:iFlux) = NameChannelSourceTmp_I
          EChannelType_I(FluxFirst_:iFlux) = EChannelTypeTmp_I
          EChannelLoIo_I(FluxFirst_:iFlux) = EChannelLoIoTmp_I
          EChannelHiIo_I(FluxFirst_:iFlux) = EChannelHiIoTmp_I
          ! Ready for the satellite energy channel assignments
          iFlux = iFlux + 1

          deallocate(NameChannelSourceTmp_I, EChannelTypeTmp_I, &
               EChannelLoIoTmp_I, EChannelHiIoTmp_I)
       else
          iFlux = FluxFirst_
       end if

       ! For the rest channels: satellite energy channels
       do iSat = 1, nFluxChannelSat
          NameChannelSource_I(iFlux:iFlux-1+ &
               FluxChannelSat_I(iSat)%nFluxChannel) = &
               trim(FluxChannelSat_I(iSat)%NameSat)
          EChannelType_I(iFlux:iFlux-1+ &
               FluxChannelSat_I(iSat)%nFluxChannel) = &
               FluxChannelSat_I(iSat)%iKindFlux
          EChannelLoIo_I(iFlux:iFlux-1+ &
               FluxChannelSat_I(iSat)%nFluxChannel) = &
               FluxChannelSat_I(iSat)%EChannelLoIo_I
          EChannelHiIo_I(iFlux:iFlux-1+ &
               FluxChannelSat_I(iSat)%nFluxChannel) = &
               FluxChannelSat_I(iSat)%EChannelHiIo_I
          iFlux = iFlux + FluxChannelSat_I(iSat)%nFluxChannel
       end do
    end if

    ! Calculate the energy channel bin width (inverse). Note that we allow
    ! nFluxChannel == 0 to only save Flux0_ and EFlux_, so we should be careful
    if(nFluxChannel > 0) then
       allocate(InvEChannelBinIo_I(FluxFirst_:FluxLast_))
       InvEChannelBinIo_I = 1.0/(EChannelHiIo_I - EChannelLoIo_I)
    end if

    ! Write the header of the energy channels for particle flux
    allocate(character(len=14) :: NameFluxChannel_I(Flux0_:FluxMax_))
    do iFlux = FluxFirst_, FluxLast_
       write(NameFluxChannel,'(I2.2)') iFlux
       NameFluxChannel_I(iFlux) = 'flux_Channel' // NameFluxChannel
    end do
    ! Write the header of total particle and energy fluxes
    NameFluxChannel_I(Flux0_) = 'flux_total    '
    NameFluxChannel_I(EFlux_) = 'eflux         '

    ! Initialize NameFluxUnit_I for different energy channels
    if(allocated(NameFluxUnit_I)) deallocate(NameFluxUnit_I)
    allocate(NameFluxUnit_I(Flux0_:FluxMax_))
    NameFluxUnit_I(Flux0_:FluxLast_) = NameFluxUnit
    ! In case that nFluxChannel == 0
    if(allocated(EChannelType_I)) then
       where(EChannelType_I(FluxFirst_:FluxLast_) == FluxDiff_)
          NameFluxUnit_I(FluxFirst_:FluxLast_) = NameDiffFluxUnit
       end where
    end if
    NameFluxUnit_I(EFlux_) = NameEnergyFluxUnit

    ! Finally allocate Flux_VIB, saved for outputs
    if(.not.allocated(Flux_VIB)) then
       allocate(Flux_VIB(Flux0_:FluxMax_, 1:nVertexMax, nLine), stat=iError)
       call check_allocate(iError, 'Flux_VIB')
       Flux_VIB = -1.0
    else
       call CON_stop(NameSub // ' Flux_VIB already allocated')
    end if

    ! Fill initial values of flux in energy channels
    allocate(FluxChannelInit_V(Flux0_:FluxMax_))
    ! FluxChannelInit_V, for the assumed initial distribution (~1/p^2)
    FluxChannelInit_V(Flux0_) = FluxInitIo
    FluxChannelInit_V(FluxFirst_:FluxLast_) = FluxInitIo * &
         (EChannelHiIo_I - EChannelLoIo_I) / (EnergyMaxIo - EnergyInjIo)
    FluxChannelInit_V(EFlux_) = FluxInitIo*0.5*(EnergyMaxIo + EnergyInjIo)

  end subroutine init
  !============================================================================
  subroutine get_integral_flux
    ! compute the total (simulated) integral flux of particles and particle
    ! fluxes in specified energy channels; also compute total energy flux

    use ModNumConst,        ONLY: cTiny
    use SP_ModDistribution, ONLY: nP, nMu, DeltaMu, Momentum_I, &
         KinEnergyIo_I, Distribution_CB, search_kinetic_energy
    ! loop variables
    integer:: iLine, iVertex, iP, iFlux, iPLo, iPHi
    ! VDF*Momentum_I**2, VDF*Momentum_I**2*KinEnergy
    real   :: DistTimesP2_I(1:nP), DistTimesP2E_I(1:nP)
    ! increments of each bin for particle and energy fluxes
    real   :: dFlux_I(1:nP-1), dEFlux_I(1:nP-1)
    ! increments of the bin where the channel falls in
    real   :: dFluxChannelLo, dFluxChannelHi
    ! logical variables, whether the iKinEnergyOut index is found
    logical:: IsFoundChannelIo, IsFoundChannelHi
    ! particle flux of energy channels
    real   :: FluxLo_I(FluxFirst_:FluxLast_), FluxHi_I(FluxFirst_:FluxLast_)

    !--------------------------------------------------------------------------
    do iLine = 1, nLine
       if(.not.Used_B(iLine)) CYCLE
       do iVertex = 1, nVertex_B(iLine)
          ! Calculate intermediate variables
          DistTimesP2_I = sum(Distribution_CB(1:nP, 1:nMu, iVertex, iLine), &
               DIM=2)*DeltaMu*0.5 * Momentum_I(1:nP)**2
          DistTimesP2E_I = DistTimesP2_I*KinEnergyIo_I(1:nP)
          ! Increment of particle and energy fluxes
          ! Integration loop with midpoint rule
          dFlux_I  = 0.5*(KinEnergyIo_I(2:nP) - KinEnergyIo_I(1:nP-1))* &
               (DistTimesP2_I(1:nP-1) + DistTimesP2_I(2:nP))
          dEFlux_I = 0.5*(KinEnergyIo_I(2:nP) - KinEnergyIo_I(1:nP-1))* &
               (DistTimesP2E_I(1:nP-1) + DistTimesP2E_I(2:nP))

          ! Calculate and store the total particle and energy fluxes
          Flux_VIB(Flux0_, iVertex, iLine) = sum(dFlux_I )*Si2Io_V(UnitFlux_)
          Flux_VIB(EFlux_, iVertex, iLine) = sum(dEFlux_I)*Si2Io_V(UnitFlux_)

          ! Reset values
          FluxLo_I = 0.0; FluxHi_I = 0.0
          ! Calculate the particle fluxes for all energy channels
          do iFlux = FluxFirst_, FluxLast_
             ! 1) For the flux from EnergyLo to EnergyMax
             ! Get the partial flux increments
             call search_kinetic_energy(EChannelLoIo_I(iFlux), iPLo, &
                  IsFoundChannelIo, DistTimesP2_I, dFluxChannelLo)
             ! iPLo: from 1 to nP-1, due to the range of EChannelLoIo_I
             if(iPLo == nP-1) then
                ! Only contributed from nP-1 to nP, when iPLo == nP-1
                FluxLo_I(iFlux) = dFluxChannelLo
             else
                ! For the rest bins: make a summation, when iPLo < nP-1
                FluxLo_I(iFlux) = dFluxChannelLo + sum(dFlux_I(iPLo+1:nP-1))
             end if

             ! 2) For the flux from EnergyHi to EnergyMax
             if(abs(EChannelHiIo_I(iFlux)-KinEnergyIo_I(nP)) < cTiny) then
                FluxHi_I(iFlux) = 0.0
             else
                ! Get the partial flux increments
                call search_kinetic_energy(EChannelHiIo_I(iFlux), iPHi, &
                     IsFoundChannelHi, DistTimesP2_I, dFluxChannelHi)
                ! iPHi: from 1 to nP-1, due to the range of EChannelHiIo_I
                if(iPHi == nP-1) then
                   ! Only contributed from nP-1 to nP, when iPHi == nP-1
                   FluxHi_I(iFlux) = dFluxChannelHi
                else
                   ! For the rest bins: make a summation, when iPHi < nP-1
                   FluxHi_I(iFlux) = dFluxChannelHi + sum(dFlux_I(iPHi+1:nP-1))
                end if
             end if
          end do

          ! Store the results for specified energy channels
          Flux_VIB(FluxFirst_:FluxLast_, iVertex, iLine) = &
               (FluxLo_I - FluxHi_I)*Si2Io_V(UnitFlux_)
          ! In case that nFluxChannel == 0
          if(allocated(EChannelType_I)) then
             where(EChannelType_I(FluxFirst_:FluxLast_) == FluxDiff_)
                Flux_VIB(FluxFirst_:FluxLast_, iVertex, iLine) = &
                     Flux_VIB(FluxFirst_:FluxLast_, iVertex, iLine) &
                     *InvEChannelBinIo_I
             end where
          end if
       end do
    end do
  end subroutine get_integral_flux
  !============================================================================
end module SP_ModChannel
!==============================================================================
