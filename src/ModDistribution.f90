!  Copyright (C) 2002 Regents of the University of Michigan, 
!  portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!==================================================================
module SP_ModDistribution
  ! The module contains the velocity/momentum distribution function
  ! and methods for initializing it as well as the offset routine.
  use ModNumConst,ONLY: cTiny
  use ModConst,   ONLY: cLightSpeed, energy_in
  use SP_ModSize, ONLY: nParticleMax, nP=>nMomentum
  use SP_ModUnit, ONLY: IO2SI_KinEnergy, &
       UnitFluxSI => UnitParticleFluxSI,&
       kinetic_energy_to_momentum, momentum_to_energy
  use SP_ModGrid, ONLY: nBlock, nParticle_B
  implicit none
  SAVE
  PRIVATE

  !Public members:
  public:: init              ! Initialize Distribution_IIB
  public:: read_param        ! Read momentum grid parameters  
  public:: offset            ! Sync. index in State_VIB and Dist_IIB 
  public:: get_integral_flux ! Calculate Flux_VIB
  public:: nP                ! Number of points in the momentum grid 
  public:: MomentumInjSI     ! Mimimum momentum value in SI
  public:: MomentumMaxSI     ! Maximum momentum value in SI
  public:: EnergyInjIo       ! Energy in kev (IO unit)
  public:: EnergyMaxIo       ! Energy in keV (IO unit)
  public:: DLogP             ! Mesh size for log(momentum) grid

  ! Injection and maximal energy in the simulation
  ! To be read from the PARAM.in file: KINETIC energies
  real:: EnergyInjIo=10.0, EnergyMaxIo=1.0E+07

  ! Injection and max momentum in the simulation
  real:: MomentumInjSI, MomentumMaxSI

  ! Size of a  log-momentum mesh. For momentum we use both the 
  ! denotaion, P, and a word, momentum - whichever is more covenient
  real:: DLogP        !log(MomentumMaxSI/MomentumInjSI)/nP

  ! Momentum, kinetic energy and total energy (including the rest mass energy) 
  ! at the momentum grid points
  real, public, dimension(0:nP+1) :: MomentumSI_I, KinEnergySI_I, EnergySI_I

  ! Total integral (simulated) particle flux
  integer, parameter, public :: Flux0_ = 0
  integer, parameter, public :: FluxFirst_ =1  ! The first channel
  integer, public :: nFluxChannel = 6          ! GOES by default, 6 channels
  integer, public :: FluxLast_ = 6             ! The last channel
  integer, public :: EFlux_    = 7             ! Total integral energy flux
  integer, public :: FluxMax_  = 7
  real, allocatable :: EChannel_I(:) ! energy limits of the instrument
  real, public, allocatable  :: Flux_VIB( :,:,:)
  character(len=10), public, allocatable  :: NameFluxChannel_I(:)


  !-----------------Grid in the momentum space---------------------------------
  ! iP     0     1                         nP   nP+1
  !        |     |    ....                 |     | 
  ! P      P_inj P_inj*exp(\Delta(log P))  P_Max P_Max*exp(DLogP)
  !----------------------------------------------------------------------------

  ! This is because we put two boundary conditions: the background
  ! value at the right one and the physical condition at the left  
  ! one, for the velocity distribution function

  ! Velosity Distribution Function (VDF) 
  ! Number of points along the momentum axis is set in ModSize
  ! 1st index - log(momentum)
  ! 2nd index - particle index along the field line
  ! 3rd index - local block number
  real, public, allocatable :: Distribution_IIB(:,:,:)

  ! distribution is initialized to have integral flux:
  real:: FluxInitIo = 0.01 ! [PFU]

  logical :: DoInit = .true.
contains
  !================================================================
  subroutine init
    use SP_ModUnit,   ONLY: momentum_to_kinetic_energy
    use ModUtilities, ONLY: check_allocate
    ! set the initial distribution on all lines
    integer:: iBlock, iParticle, iP, iError
    character(len=*), parameter :: NameSub = 'SP_init'
    !----------------------------------------------------------
    if(.not.DoInit)RETURN
    DoInit = .false.

    ! convert energies to momenta
    MomentumInjSI  = kinetic_energy_to_momentum(EnergyInjIo*IO2SI_KinEnergy)
    MomentumMaxSI  = kinetic_energy_to_momentum(EnergyMaxIo*IO2SI_KinEnergy)  

    ! grid size in the log momentum space
    DLogP = log(MomentumMaxSI/MomentumInjSI)/nP

    !\
    ! Functions to convert the grid index to momentum and energy
    !/
    do iP = 0, nP +1
       MomentumSI_I(iP)  = MomentumInjSI*exp(iP*DLogP)
       KinEnergySI_I(iP) = momentum_to_kinetic_energy(MomentumSI_I(iP))
       EnergySI_I(iP)    = momentum_to_energy(MomentumSI_I(iP))
    end do

    !\
    ! Distribution function
    !/    
    allocate(Distribution_IIB(0:nP+1,1:nParticleMax,nBlock), stat=iError)
    call check_allocate(iError, 'Distribution_IIB')

    ! initialization depends on momentum, however, this corresponds
    ! to a constant differential flux (intensity), thus ensuring
    ! uniform backgound while visualizing this quantity
    do iBlock = 1, nBlock
       do iParticle = 1, nParticleMax
          ! Overall density of the fast particles is of the order 
          ! of 10^-6 m^-3. Integral flux is less than 100 per
          ! (m^2 ster s). Differential background flux is constant.
          do iP = 0, nP +1
             Distribution_IIB(iP,iParticle,iBlock) =                         &
                  FluxInitIo/(EnergyMaxIo-EnergyInjIo)/MomentumSI_I(iP)**2   &
                  *UnitFluxSI/IO2SI_KinEnergy

          end do
       end do
    end do

    if (.not. allocated(EChannel_I)) then
       if (nFluxChannel == 6) then
          allocate (EChannel_I(nFluxChannel))
          EChannel_I = (/5,10,30,50,60,100/)
       else
          call CON_stop(NameSub//' check nFluxChannel ')
       end if
    end if

    if (.not. allocated(Flux_VIB)) then
       allocate(Flux_VIB(Flux0_:FluxMax_,1:nParticleMax,nBlock), &
            stat=iError); call check_allocate(iError, 'Flux_VIB')
       Flux_VIB = -1.0
    else
       call CON_stop(NameSub//' Flux_VIB already allocated')
    end if

  end subroutine init
  !================================================================
  subroutine read_param(NameCommand)
    use ModReadParam, ONLY: read_var
    use ModUtilities, ONLY: lower_case
    use SP_ModProc,   ONLY: iProc
    use SP_ModGrid  , ONLY: T_
    use SP_ModUnit,   ONLY: NameVarUnit_V, NameEnergyUnit
    character(len=*), intent(in):: NameCommand ! From PARAM.in  
    character(len=*), parameter :: NameSub='SP:read_param_dist'
    integer:: nPCheck = nP, iFluxChannel
    character(len=5) :: NameFluxChannel
    !--------------------------------------------------------------
    select case(NameCommand)
    case('#MOMENTUMGRID')
       !Read unit to be used for particle energy: eV, keV, GeV
       call read_var('EnergyMin',      EnergyInjIo)
       call read_var('EnergyMax',      EnergyMaxIo)
       call read_var('nP',             nPCheck    )
       call read_var('NameEnergyUnit', NameEnergyUnit)

       call lower_case(NameEnergyUnit)

       ! Unfortunately this is needed. The reason is that the code does not
       ! allow that the unit of the kinetic energy of the plasma temperature 
       ! to be different from the unit of the kinetic energy of particles.
       ! The default NameVarUnit_V(T_) ='kev', so if NameEnergyUnit is updated,
       ! NameVarUnit_V(T_) should also be updated. In stand along mode, the 
       ! code is happy, but should be tested in the coupling run.
       NameVarUnit_V(T_) = NameEnergyUnit

       if(nP/=nPCheck)then
          if(iProc==0)write(*,'(a,i6,a,i6)')NameSub//' '//         &
               'Code is configured with nMomentum=', nP ,          &
               ' while value read from PARAM.in is nP=',nPCheck 
          call CON_stop('Modify PARAM.in or reconfigure SP/MFLAMPA')
       end if
    case('#FLUXINITIAL')
       call read_var('FluxInit [p.f.u.]',FluxInitIo)
       ! check correctness
       if(FluxInitIo<=0)call CON_stop(NameSub//': flux value must be positive')
    case('#FLUXCHANNEL')
       call read_var('nFluxChannel', nFluxChannel)
       FluxLast_ = nFluxChannel
       EFlux_    = FluxLast_ + 1
       FluxMax_  = EFlux_

       if (allocated(EChannel_I)) deallocate(EChannel_I)
       allocate(EChannel_I(nFluxChannel))
       if (allocated(NameFluxChannel_I)) deallocate(NameFluxChannel_I)
       allocate(NameFluxChannel_I(0:nFluxChannel+1))

       NameFluxChannel_I(0)              = 'flux_total'
       NameFluxChannel_I(nFluxChannel+1) = 'eflux'

       do iFluxChannel=1,nFluxChannel
          call read_var('EChannel_I', EChannel_I(iFluxChannel))
          write(NameFluxChannel,'(I5.5)') int(EChannel_I(iFluxChannel))
          NameFluxChannel_I(iFluxChannel) = 'flux_'//NameFluxChannel
       end do
    case default
       call CON_stop(NameSub//' Unknown command '//NameCommand)
    end select
  end subroutine read_param
  !================================================================
  subroutine offset(iBlock, iOffset)
    use SP_ModGrid, ONLY: NoShock_, BOld_, RhoOld_, ShockOld_, &
         iShock_IB,  State_VIB, X_, Z_, FootPoint_VB
    ! shift in the data arrays is required if the grid point(s) is  
    ! appended or removed at the foot point of the magnetic field 
    ! line. SHIFTED ARE: State_VIB(/RhoOld_,BOld_),Distribution_IIB
    ! as well as ShockOld_
    integer, intent(in)        :: iBlock
    integer, intent(in)        :: iOffset
    real :: Alpha, Distance2ToMin, Distance3To2
    character(len=*), parameter :: NameSub = "SP: offset"
    !------------
    if(iOffset==0)RETURN
    if(iOffset==1)then
       State_VIB((/RhoOld_,BOld_/),2:nParticle_B(iBlock),iBlock) &
            = State_VIB((/RhoOld_,BOld_/),1:nParticle_B(iBlock)-1,iBlock)
       Distribution_IIB(:,2:nParticle_B(iBlock), iBlock)&
            = Distribution_IIB(:,1:nParticle_B(iBlock)-1, iBlock)
       !\
       ! Extrapolate state vector components and VDF at iparticle=1
       Distance2ToMin = sqrt(sum((State_VIB(X_:Z_,2,iBlock) - &
               FootPoint_VB(X_:Z_,iBlock))**2))
          Distance3To2   = sqrt(sum((State_VIB(X_:Z_,3,iBlock) - &
                State_VIB(X_:Z_,2,iBlock))**2))
          Alpha = Distance2ToMin/(Distance2ToMin + Distance3To2)
       State_VIB((/RhoOld_, BOld_/), 1, iBlock) = &
            (Alpha + 1)*State_VIB((/RhoOld_, BOld_/), 2, iBlock) &
            -Alpha     * State_VIB((/RhoOld_, BOld_/), 3, iBlock)
       Distribution_IIB(:,1,iBlock) = Distribution_IIB(:,2,iBlock) + &
            Alpha*(Distribution_IIB(:,2,iBlock) - &
            Distribution_IIB(:,3,iBlock))       
       ! extrapolation may introduced negative values 
       ! for strictly positive quantities; such occurences need fixing
       where(State_VIB((/RhoOld_,BOld_/),1,iBlock) <= 0.0)
          State_VIB((/RhoOld_,BOld_/),1,iBlock) = &
               0.01 * State_VIB((/RhoOld_,BOld_/),2,iBlock)
       end where
       where(Distribution_IIB(:,1,iBlock) <= 0.0)
          Distribution_IIB(:,1,iBlock) = &
               0.01 * Distribution_IIB(:,2,iBlock)
       end where
    elseif(iOffset < 0)then
       State_VIB((/RhoOld_,BOld_/),1:nParticle_B(iBlock),iBlock) &
            =  State_VIB((/RhoOld_,BOld_/),1-iOffset:nParticle_B(iBlock)&
       - iOffset, iBlock)
       Distribution_IIB(:,1:nParticle_B(iBlock), iBlock)&
            = Distribution_IIB(:,1-iOffset:nParticle_B(iBlock)-iOffset, &
            iBlock)
    else
       call CON_stop('No algorithm for iOffset >1 in '//NameSub)
    end if
    if(iShock_IB(ShockOld_, iBlock)/=NoShock_)&
         iShock_IB(ShockOld_, iBlock) = &
         max(iShock_IB(ShockOld_, iBlock) + iOffset, 1)
  end subroutine offset
  !===========================================================================
  subroutine get_integral_flux
    use ModConst, ONLY: energy_in
    ! compute the total (simulated) integral flux of particles as well as
    ! particle flux in the 6 GOES channels; also compute total energy flux
    !------------------------------------------------------------------------
    integer:: iBlock, iParticle, iP, iFlux ! loop variables
    real   :: EFlux ! the value of energy flux
    real   :: dFlux, dFlux1 ! increments
    real   :: Flux          ! the value of particle flux
    real   :: Norm            ! normalization factor
    real, allocatable :: Flux_I(:), EChannelSI_I(:)
    !-------------------------------------------------------------------------
    ! energy limits of GOES channels


    if (.not.allocated(Flux_I)) allocate(Flux_I(nFluxChannel))
    if (.not.allocated(EChannelSI_I)) allocate(EChannelSI_I(nFluxChannel))

    EChannelSI_I = EChannel_I * energy_in('MeV')

    do iBlock = 1, nBlock
       do iParticle = 1, nParticle_B( iBlock)
          !\
          ! Integration loop with midpoint rule
          !/
          ! reset values
          EFlux = 0.0
          Flux_I= 0.0
          Norm  = 0.0
          Flux  = 0.0
          do iP = 1, nP - 1
             ! the flux increment from iP
             dFlux = 0.5 * &
                  (KinEnergySI_I(iP+1) - KinEnergySI_I(iP)) * (&
                  Distribution_IIB(iP,  iParticle,iBlock)*&
                  MomentumSI_I(iP)**2 &
                  +&
                  Distribution_IIB(iP+1,iParticle,iBlock)*&
                  MomentumSI_I(iP+1)**2)

             ! increase the total flux
             Flux = Flux + dFlux

             ! increase FOES channels' fluxes
             do iFlux = 1, nFluxChannel
                ! check whether reached the channel's cut-off level
                if(KinEnergySI_I(iP+1) < EChannelSI_I(iFlux))&
                     CYCLE

                if(KinEnergySI_I(iP) >= EChannelSI_I(iFlux))then
                   Flux_I(iFlux) = Flux_I(iFlux) + dFlux
                else
                   ! channel cutoff level is often in the middle of a bin;
                   ! compute partial flux increment
                   dFlux1 =&
                        ((-0.50*(KinEnergySI_I(iP) + EChannelSI_I(iFlux)) + &
                        KinEnergySI_I(iP+1) )*&
                        Distribution_IIB(iP,iParticle,iBlock)*&
                        MomentumSI_I(iP)**2  &
                        -0.50*(KinEnergySI_I(iP)-EChannelSI_I(iFlux))*&
                        Distribution_IIB(iP+1,iParticle,iBlock)*&
                        MomentumSI_I(iP+1)**2)*&
                        (KinEnergySI_I(iP)-EChannelSI_I(iFlux))/&
                        (KinEnergySI_I(iP+1)-KinEnergySI_I(iP))
                   Flux_I(iFlux) = Flux_I(iFlux) + dFlux1
                end if
             end do

             ! increase total energy flux
             EFlux = EFlux + 0.5 * &
                  (KinEnergySI_I(iP+1) - KinEnergySI_I(iP)) * (&
                  Distribution_IIB(iP,  iParticle,iBlock)*&
                  KinEnergySI_I(iP) * &
                  MomentumSI_I(iP)**2 &
                  +&
                  Distribution_IIB(iP+1,iParticle,iBlock)*&
                  KinEnergySI_I(iP+1) * &
                  MomentumSI_I(iP+1)**2)

             ! normalization factor
             Norm = Norm + 0.5 * &
                  (MomentumSI_I(iP+1) - MomentumSI_I(iP)) * (&
                  Distribution_IIB(iP,  iParticle, iBlock) * &
                  MomentumSI_I(iP)**2 &
                  + &
                  Distribution_IIB(iP+1,iParticle, iBlock) * &
                  MomentumSI_I(iP+1)**2)
          end do

          ! store the results
          Flux_VIB(Flux0_,       iParticle, iBlock) = Flux
          Flux_VIB(FluxFirst_:FluxLast_,iParticle, iBlock) = Flux_I
          Flux_VIB(EFlux_,       iParticle, iBlock) = EFlux
       end do
    end do

    deallocate(Flux_I, EChannelSI_I)

  end subroutine get_integral_flux
end module SP_ModDistribution
