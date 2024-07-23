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
  use SP_ModSize,   ONLY: nVertexMax, nP => nMomentum, &
       nMu => nPitchAngle, IsMuAvg => IsPitchAngleAverage
  use SP_ModGrid,   ONLY: nLine, nVertex_B, Used_B

  implicit none

  SAVE

  private ! except

  ! Public members:
  public:: init              ! Initialize Distribution_CB
  public:: read_param        ! Read momentum grid parameters
  public:: offset            ! Sync. index in State_VIB and Distribution_CB
  public:: check_dist_neg    ! Check any of Distribution_CB is negative
  public:: nP                ! Number of points in the momentum grid
  public:: nMu               ! Number of points over pitch-angle (\mu)
  public:: IsMuAvg           ! If .true., dist. function is omnidirectional

  ! Injection and maximal energy, in kev (or, in Io energy unit)
  ! To be read from the PARAM.in file
  real, public :: EnergyInjIo = 10.0, EnergyMaxIo = 1.0E+07

  ! Injection momentum, mimimum momentum value in Si
  real, public :: MomentumInjSi

  ! Size of the log-momentum mesh. For momentum we use both the
  ! denotaion, P, and a word, momentum - whichever is more convenient
  real, public :: dLogP      ! log(MomentumMaxSi/MomentumInjSi)/nP

  ! uniform distribution of pitch angle
  real, public :: DeltaMu = 2.0/nMu
  real, public :: MuFace_I(0:nMu), Mu_I(nMu)

  ! speed, momentum, kinetic energy at the momentum grid points
  real, public, dimension(0:nP+1) :: SpeedSi_I, Momentum_I, &
       KinEnergyIo_I, VolumeP_I, Background_I
  real, public :: Momentum3_I(-1:nP+1) ! P**3/3, normalized by MomentumInjSi

  !-----------------Grid in the momentum space---------------------------------
  ! iP     0     1                         nP   nP+1
  !        |     |    ....                 |     |
  ! P     P_inj P_inj*exp(dLogP)          P_Max P_Max*exp(dLogP)
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
  ! 2rd index - \mu value = cosine of the pitch angle
  ! 3rd index - particle index along the field line
  ! 4th index - local field line number
  real, public, allocatable :: Distribution_CB(:,:,:,:)

  ! Check if any Distribution_CB < 0 (.true. for such case)
  logical, public :: IsDistNeg = .false.

  ! distribution is initialized to have integral flux:
  real, public :: FluxInitIo = 0.01 ! [PFU]

  ! Logical variable: whether this is the initial call
  logical :: DoInit = .true.
contains
  !============================================================================
  subroutine read_param(NameCommand)

    use ModReadParam, ONLY: read_var
    use SP_ModProc,   ONLY: iProc
    character(len=*), intent(in) :: NameCommand ! From PARAM.in
    integer :: nPCheck = nP, nMuCheck = nMu
    character(len=*), parameter:: NameSub = 'read_param'
    !--------------------------------------------------------------------------
    select case(NameCommand)
    case('#MOMENTUMGRID')
       ! Read energy range for particles, in the unit of NameEnergyUnit
       call read_var('EnergyMin', EnergyInjIo)
       call read_var('EnergyMax', EnergyMaxIo)
       call read_var('nP',        nPCheck    )
       if(nP/=nPCheck) then
          if(iProc==0) write(*,'(a,i6,a,i6)') NameSub//' '//     &
               'Code is configured with nMomentum=', nP,         &
               ' while value read from PARAM.in is nP=', nPCheck
          call CON_stop('Modify PARAM.in or reconfigure SP/MFLAMPA')
       end if
    case('#PITCHANGLEGRID')
       call read_var('nMu', nMuCheck)
       if(nMu/=nMuCheck) then
          if(iProc==0) write(*,'(a,i6,a,i6)') NameSub//' '//     &
               'Code is configured with nMu=', nMu,              &
               ' while value read from PARAM.in is nMu=', nMuCheck
          call CON_stop('Modify PARAM.in or reconfigure SP/MFLAMPA')
       end if
    case('#FLUXINITIAL')
       call read_var('FluxInitIo', FluxInitIo)
       ! check correctness
       if(FluxInitIo <= 0.0) call CON_stop(NameSub // &
            ': flux value must be positive')
    case default
       call CON_stop(NameSub // ' Unknown command ' // NameCommand)
    end select

  end subroutine read_param
  !============================================================================
  subroutine init

    use ModConst,     ONLY: cLightSpeed
    use ModUtilities, ONLY: check_allocate
    use SP_ModProc,   ONLY: iError
    use SP_ModUnit,   ONLY: kinetic_energy_to_momentum, momentum_to_energy, &
         momentum_to_kinetic_energy, Io2Si_V, Si2Io_V, UnitEnergy_, UnitFlux_
    ! loop variables
    integer :: iLine, iVertex, iP, iMu
    ! maximal momentum
    real :: MomentumMaxSi
    ! local FluxChannel, for converting to NameFluxChannel_I
    real :: FluxChannel
    ! local NameFluxChannel and NameUnitChannel, written in headers
    character(len=3) :: NameFluxChannel, NameUnitChannel
    character(len=*), parameter:: NameSub = 'init'
    !--------------------------------------------------------------------------
    if(.not.DoInit) RETURN
    DoInit = .false.

    ! convert energies to momenta
    MomentumInjSi= kinetic_energy_to_momentum(EnergyInjIo*Io2Si_V(UnitEnergy_))
    MomentumMaxSi= kinetic_energy_to_momentum(EnergyMaxIo*Io2Si_V(UnitEnergy_))
    ! grid size in the log momentum space
    dLogP = log(MomentumMaxSi/MomentumInjSi)/nP

    ! Functions to convert the grid index to momentum and energy
    Momentum3_I(-1) = exp(-3*(0.50*dLogP))/3.0 ! P^3/3 at -0.5*dLogP from PInj
    do iP = 0, nP+1
       Momentum_I(iP)    = MomentumInjSi*exp(iP*dLogP)
       SpeedSi_I(iP)     = Momentum_I(iP)*cLightSpeed**2/ &
            momentum_to_energy(Momentum_I(iP))
       Momentum3_I(iP)   = Momentum3_I(iP-1)*exp(3*dLogP)
       VolumeP_I(iP)     = Momentum3_I(iP) - Momentum3_I(iP-1)
       ! Normalize kinetic energy per Unit of energy in SI unit:
       KinEnergyIo_I(iP) = momentum_to_kinetic_energy(Momentum_I(iP)) &
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
    Mu_I = 0.5*(MuFace_I(0:nMu-1) + MuFace_I(1:nMu))

    ! Distribution function
    allocate(Distribution_CB(0:nP+1, nMu, nVertexMax, nLine), stat=iError)
    call check_allocate(iError, 'Distribution_CB')

    ! set the initial distribution on all lines
    ! initialization depends on momentum, however, this corresponds
    ! to a constant differential flux (intensity), thus ensuring
    ! uniform backgound while visualizing this quantity
    Distribution_CB = reshape(spread(Background_I, DIM=2, &
         NCOPIES=nMu*nVertexMax*nLine), [nP+2, nMu, nVertexMax, nLine])
    ! (0:nP+1, nMu*nVertexMax*nLine) -> (0:nP+1, 1:nMu, 1:nVertexMax, 1:nLine)
    ! Overall density of the fast particles is of the order
    ! of 10^-6 m^-3. Integral flux is less than 100 per
    ! (m^2 ster s). Differential background flux is constant.
  end subroutine init
  !============================================================================
  subroutine offset(iLine, iOffset)

    use SP_ModGrid, ONLY: NoShock_, BOld_, RhoOld_, ShockOld_, &
         iShock_IB, State_VIB, MHData_VIB, X_, Z_, FootPoint_VB
    ! shift in the data arrays is required if the grid point(s) is appended
    ! or removed at the foot point of the magnetic field line. SHIFTED ARE:
    ! State_VIB(/RhoOld_,BOld_), Distribution_CB, as well as ShockOld_
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
       Distance2ToMin = norm2(MHData_VIB(X_:Z_, 2, iLine) - &
            FootPoint_VB(X_:Z_, iLine))
       Distance3To2   = norm2(MHData_VIB(X_:Z_, 3, iLine) - &
            MHData_VIB(X_:Z_, 2, iLine))
       Alpha = Distance2ToMin/(Distance2ToMin + Distance3To2)
       State_VIB([RhoOld_, BOld_], 1, iLine) = &
            (Alpha + 1)*State_VIB([RhoOld_, BOld_], 2, iLine) &
            -Alpha     *State_VIB([RhoOld_, BOld_], 3, iLine)
       Distribution_CB(:,:,1,iLine) = Distribution_CB(:,:,2,iLine) + Alpha* &
            (Distribution_CB(:,:,2,iLine) - Distribution_CB(:,:,3,iLine))
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
            = State_VIB([RhoOld_,BOld_],1-iOffset:nVertex_B(iLine)&
            - iOffset, iLine)
       Distribution_CB(:,:,1:nVertex_B(iLine), iLine) &
            = Distribution_CB(:,:,1-iOffset:nVertex_B(iLine)-iOffset, iLine)
    else
       call CON_stop('No algorithm for iOffset > 1 in '//NameSub)
    end if
    if(iShock_IB(ShockOld_, iLine)/=NoShock_)&
         iShock_IB(ShockOld_, iLine) = &
         max(iShock_IB(ShockOld_, iLine) + iOffset, 1)
  end subroutine offset
  !============================================================================
  subroutine check_dist_neg(NameSub, lVertex, rVertex, iLine)

    ! check if any of Distribution_CB(:, :, lVertex:rVertex, iLine) < 0.0
    ! if so, IsDistNeg = .true.; otherwise, IsDistNeg = .false. (normal case)
    character(len=*), intent(in):: NameSub  ! String for the module/subroutine
    integer, intent(in) :: lVertex, rVertex ! Start and end indices of iVertex
    integer, intent(in) :: iLine            ! index of line
    !--------------------------------------------------------------------------
    if(any(Distribution_CB(:, :, lVertex:rVertex, iLine)<0.0)) then
       write(*,*) NameSub, ': Distribution_CB < 0'
       Used_B(iLine) = .false.
       nVertex_B(iLine) = 0
       IsDistNeg = .true.
       ! With negative values in the Poisson bracket scheme: should stop here
       if( index(NameSub, 'poisson')>0 .or. index(NameSub, 'Poisson')>0 &
            .or. index(NameSub, 'POISSON')>0 ) &
            call CON_stop(NameSub//': Distribution_CB < 0')
    end if
  end subroutine check_dist_neg
  !============================================================================
end module SP_ModDistribution
!==============================================================================
