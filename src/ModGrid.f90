module SP_ModGrid

  use SP_ModSize, ONLY: &
       nDim, nLat, nLon, nNode, nMomentumBin, nPitchAngleBin, &
       nParticleMax, Particle_, OriginLat_, OriginLon_

  implicit none

  SAVE

  private ! except

  public:: set_grid_param, init_grid, get_node_indexes, distance_to_next
  public:: append_particles
  public:: iComm, iProc, nProc, nBlock, Proc_, Block_, nBlockIndexes
  public:: LatMin, LatMax, LonMin, LonMax
  public:: RMin, RBufferMin, RBufferMax, RMax, ROrigin
  public:: iGridGlobal_IA, iGridLocal_IB, FootPoint_VB, iNode_II, iNode_B
  public:: State_VIB, Flux_VIB, Distribution_IIB
  public:: MomentumScale_I, LogMomentumScale_I, EnergyScale_I, LogEnergyScale_I
  public:: DMomentumOverDEnergy_I
  public:: nParticle_B, Shock_, ShockOld_, Length_
  public:: nVar, nVarRead,  X_, Y_, Z_, D_, S_, LagrID_, Offset_
  public:: Rho_,T_, Ux_,Uy_,Uz_,U_,DLogRho_, Bx_,By_,Bz_,B_, RhoOld_,BOld_
  public:: EFlux_, Flux0_, Flux1_, Flux2_, Flux3_, Flux4_, Flux5_, Flux6_
  public:: Wave1_, Wave2_,FluxMax_
  public:: NameVar_V
  public:: TypeCoordSystem

  !\
  ! MPI information
  !----------------------------------------------------------------------------
  integer:: iComm = -1
  integer:: iProc = -1
  integer:: nProc = -1
  !/
  !\
  ! Grid info
  ! Containers for coordinates and data
  !----------------------------------------------------------------------------
  ! Starting position of field lines in Rs
  real:: ROrigin = 2.5
  ! Size of angular grid, latitude and longitude, at origin surface R=ROrigin
  real:: LatMin, LatMax, DLat
  real:: LonMin, LonMax, DLon
  ! Lower boundary of the domain in Rs
  real:: RMin=-1.
  ! Upper boundary of the domain in Rs
  real:: RMax=-1.
  ! Boundaries of the buffer layer between SC and IH Rs
  real:: RBufferMin=-1.
  real:: RBufferMax=-1.
  ! Mark that grid or lines' origin have been set
  logical:: IsSetGrid   = .false.
  logical:: IsSetOrigin = .false.
  !----------------------------------------------------------------------------
  !----------------------------------------------------------------------------
  ! Node number based on the field line identified by 2 angular grid indices,
  ! latitude and longitude;
  ! 1st index - latitude index
  ! 2nd index - longitude index
  integer, allocatable:: iNode_II(:,:)
  !----------------------------------------------------------------------------
  ! Number of blocks on this processor
  integer:: nBlock
  !----------------------------------------------------------------------------
  ! Node number based on the local block number
  ! 1st index - block number
  integer, allocatable:: iNode_B(:)
  !----------------------------------------------------------------------------
  ! Various house-keeping information about the node/line;
  ! 1st index - identification of info field
  ! 2nd index - node number / block number
  integer, allocatable:: iGridGlobal_IA(:,:)
  integer, allocatable:: iGridLocal_IB(:,:), nParticle_B(:)
  real,    allocatable:: FootPoint_VB(:,:)
  !----------------------------------------------------------------------------
  ! Number of info fields per node/block and their identifications
  integer, parameter:: nNodeIndexes = 2
  integer, parameter:: &
       Proc_  = 1, & ! Processor that has this line/node
       Block_ = 2    ! Block that has this line/node
  integer, parameter:: nBlockIndexes = 3
  integer, parameter:: &
       Shock_   = 1, & ! Current location of a shock wave
       ShockOld_= 2, & ! Old location of a shock wave
       Offset_  = 3    ! To account for the dymaical grid distinction 
                       ! from that updated in the other components
  integer, parameter:: & 
       Length_ = 4    ! init length of segment 1-2: control for new particles
                      ! being appended to the beginnings of lines
  !----------------------------------------------------------------------------
  ! State vector;
  ! 1st index - identification of variable
  ! 2nd index - particle index along the field line
  ! 3rd index - local block number
  real, allocatable:: State_VIB(:,:,:)
  real, allocatable:: Flux_VIB( :,:,:)
  !----------------------------------------------------------------------------
  ! Number of variables in the state vector and their identifications
  integer, parameter:: nVar     = 28, FluxMax_ = 28
  integer, parameter:: nVarRead = 13
  integer, parameter:: &
       !\
       !-- The following variables MUST be in CONTIGUOUS  order --------------
       !-- used in subroutines read_mh_data, write_restart, read_restart -----
       !-- DO NOT CHANGE WITHOUT CAREFULL CONSIDERATION !!! ------------------
       LagrID_ = 0, & ! Lagrangian id
       X_      = 1, & ! 
       Y_      = 2, & ! Cartesian coordinates
       Z_      = 3, & ! 
       Rho_    = 4, & ! Background plasma density
       T_      = 5, & ! Background temperature
       Ux_     = 6, & !
       Uy_     = 7, & ! Background plasma bulk velocity
       Uz_     = 8, & !
       Bx_     = 9, & !
       By_     =10, & ! Background magnetic field
       Bz_     =11, & !
       Wave1_  =12, & !\
       Wave2_  =13, & ! Alfven wave turbulence
       !-----------------------------------------------------------------------
       D_      =14, & ! Distance to the next particle
       S_      =15, & ! Distance from the beginning of the line
       U_      =16, & ! Magnitude of plasma bulk velocity
       B_      =17, & ! Magnitude of magnetic field
       DLogRho_=18, & ! Dln(Rho), i.e. -div(U) * Dt
       RhoOld_ =19, & ! Background plasma density
       BOld_   =20, & ! Magnitude of magnetic field
       Flux0_  =21, & ! Total integral (simulated) particle flux
       Flux1_  =22, & ! Integral particle flux >  5 MeV (GOES Channel 1)
       Flux2_  =23, & ! Integral particle flux > 10 MeV (GOES Channel 2)
       Flux3_  =24, & ! Integral particle flux > 30 MeV (GOES Channel 3)
       Flux4_  =25, & ! Integral particle flux > 50 MeV (GOES Channel 4)
       Flux5_  =26, & ! Integral particle flux > 60 MeV (GOES Channel 5)
       Flux6_  =27, & ! Integral particle flux >100 MeV (GOES Channel 6)
       EFlux_  =28    ! Total integral energy flux

  ! variable names
  character(len=10), parameter:: NameVar_V(LagrID_:EFlux_) = (/&
       'LagrID    ', &
       'X         ', &
       'Y         ', &
       'Z         ', &
       'Rho       ', &
       'T         ', &
       'Ux        ', &
       'Uy        ', &
       'Uz        ', &
       'Bx        ', &
       'By        ', &
       'Bz        ', &
       'Wave1     ', &
       'Wave2     ', &
       'D         ', &
       'S         ', &
       'U         ', &
       'B         ', &
       'DLogRho   ', &
       'RhoOld    ', &
       'BOld      ', &
       'Flux_Total', &
       'Flux_GOES1', &
       'Flux_GOES2', &
       'Flux_GOES3', &
       'Flux_GOES4', &
       'Flux_GOES5', &
       'Flux_GOES6', &
       'EFlux     '  /)
  !----------------------------------------------------------------------------
  ! Distribution vector;
  ! Number of bins in the distribution is set in ModSize
  ! 1st index - log(momentum) bin
  ! 2nd index - particle index along the field line
  ! 4th index - local block number
  real, allocatable:: Distribution_IIB(:,:,:)
  ! scale with respect to Momentum and log(Momentum)
  real:: MomentumScale_I(nMomentumBin)
  real:: LogMomentumScale_I(nMomentumBin)
  real:: EnergyScale_I(nMomentumBin)
  real:: LogEnergyScale_I(nMomentumBin)
  real:: DMomentumOverDEnergy_I(nMomentumBin)
  !----------------------------------------------------------------------------
  ! Coordinate system and geometry
  character(len=3) :: TypeCoordSystem = 'HGR'
  !/

contains
  
  subroutine set_grid_param(TypeAction)
    use ModReadParam, ONLY: read_var
    use ModNumConst, ONLY: cPi
    character (len=*), intent(in):: TypeAction ! What to do  
    character(len=*), parameter:: NameSub = 'SP:set_grid_param'
    !--------------------------------------------------------------------------
    select case(TypeAction)
    case('#ORIGIN')
       call read_var('ROrigin', ROrigin)
       call read_var('LonMin', LonMin)
       call read_var('LatMin', LatMin)
       call read_var('LonMax', LonMax)
       call read_var('LatMax', LatMax)

       ! check consistency
       if(LonMax <= LonMin .or. LatMax <= LatMin)&
            call CON_stop(NameSub//': Origin surface grid is inconsistent')
       if(ROrigin < 0.0)&
            call CON_stop(NameSub//&
            ': ROrigin, if set, must have a positive values')
       if(any((/RBufferMin, RBufferMax, RMax/) < ROrigin) .and. IsSetGrid)&
            call CON_stop(NameSub//&
            ': inconsistent values of ROrigin, RBufferMin, RBufferMax, RMax')

       ! convert angels from degrees to radians
       LonMax = LonMax * cPi / 180
       LonMin = LonMin * cPi / 180
       ! angular grid's step
       DLon = (LonMax - LonMin) / nLon

       ! convert angels from degrees to radians
       LatMax = LatMax * cPi / 180
       LatMin = LatMin * cPi / 180
       ! angular grid's step
       DLat = (LatMax - LatMin) / nLat

       IsSetOrigin = .true.
    case('#GRID')
       call read_var('RMin',RMin)
       call read_var('RBufferMin', RBufferMin)
       call read_var('RBufferMax', RBufferMax)
       call read_var('RMax',RMax)

       ! check consistency
       if(RBufferMin < 0.0 .or.RBufferMax < 0.0 .or.RMax < 0.0)&
            call CON_stop(NameSub//&
            ': RBufferMin, RBufferMax, RMax must be set to positive values')
       if(any((/RMax, RBufferMax/) < RBufferMin) .or. RMax < RBufferMax .or. &
            any((/RBufferMin, RBufferMax, RMax/) < ROrigin) .and. IsSetOrigin)&
            call CON_stop(NameSub//&
            ': inconsistent values of ROrigin, RBufferMin, RBufferMax, RMax')

       IsSetGrid = .true.
    end select
  end subroutine set_grid_param

  !============================================================================

  subroutine init_grid(DoReadInput)
    ! allocate the grid used in this model
    use ModUtilities, ONLY: check_allocate
    use ModCoordTransform, ONLY: rlonlat_to_xyz

    logical, intent(in):: DoReadInput
    integer:: iError
    integer:: iLat, iLon, iNode, iBlock, iProcNode, iParticle
    character(LEN=*),parameter:: NameSub='SP:init_grid'
    !--------------------------------------------------------------------------
    !\
    ! Check if everything's ready for initialization
    !/
    if(.not.IsSetGrid)&
         call CON_stop(NameSub//': grid is not set in PARAM.in file')
    if(.not.IsSetOrigin .and. .not. DoReadInput)&
         call CON_stop(NameSub//": neither lines' origin is set, "//&
         "nor input files are provided; change PARAM.in file!!!")
    !\
    ! distribute nodes between processors
    !/
    if(nNode < nProc)&
         call CON_stop(NameSub//': There are more processors than field lines')
    nBlock = ((iProc+1)*nNode) / nProc - (iProc*nNode) / nProc
    !\
    ! check consistency
    !/
    if(nLat <= 0 .or. nLon <= 0)&
         call CON_stop(NameSub//': Origin surface grid is invalid')
    !\
    ! allocate data and grid containers
    !/
    allocate(iNode_II(nLon, nLat), stat=iError)
    call check_allocate(iError, NameSub//'iNode_II')
    allocate(iNode_B(nBlock), stat=iError)
    call check_allocate(iError, NameSub//'iNode_B')
    allocate(iGridGlobal_IA(nNodeIndexes, nNode), stat=iError)
    call check_allocate(iError, NameSub//'iGridGlobal_IA')
    allocate(nParticle_B(nBlock), stat=iError)
    call check_allocate(iError, NameSub//'nParticle_B')
    allocate(iGridLocal_IB(nBlockIndexes, nBlock), stat=iError)
    call check_allocate(iError, NameSub//'iGridLocal_IB')
    allocate(FootPoint_VB(LagrID_:Length_, nBlock), stat=iError)
    call check_allocate(iError, NameSub//'FootPoint_VB')
    allocate(State_VIB(LagrID_:nVar,1:nParticleMax,nBlock), &
         stat=iError)
    call check_allocate(iError, NameSub//'State_VIB')
    allocate(Flux_VIB(Flux0_:FluxMax_,1:nParticleMax,nBlock), &
         stat=iError); call check_allocate(iError, 'Flux_VIB')
    allocate(Distribution_IIB(&
         nMomentumBin,1:nParticleMax,nBlock), &
         stat=iError)
    call check_allocate(iError, NameSub//'Distribution_IIB')
    !\
    ! fill grid containers
    !/
    iBlock = 1
    do iLat = 1, nLat
       do iLon = 1, nLon
          iNode = iLon + nLon * (iLat-1)
          iNode_II(iLon, iLat) = iNode
          iProcNode = ceiling(real(iNode*nProc)/nNode) - 1
          if(iProcNode==iProc)then
             iNode_B(iBlock) = iNode
             nParticle_B(     iBlock) = 1
             iGridLocal_IB(:, iBlock) = 0
          end if
          iGridGlobal_IA(Proc_,   iNode)  = iProcNode
          iGridGlobal_IA(Block_,  iNode)  = iBlock
          if(iNode == ((iProcNode+1)*nNode)/nProc)then
             iBlock = 1
          else
             iBlock = iBlock + 1
          end if
       end do
    end do
    !\
    ! reset and fill data containers
    !/
    Distribution_IIB = tiny(1.0)
    State_VIB = -1; Flux_VIB = -1; FootPoint_VB = -1
    
    !\
    ! reset lagrangian ids
    !/
    do iParticle = 1, nParticleMax
       State_VIB(LagrID_, iParticle, 1:nBlock) = real(iParticle)
    end do

    if(DoReadInput)then
       if(IsSetOrigin)&
            write(*,*)NameSub//": input files are provided, "//&
            "but lines' origin is set in PARAM.in. "//&
            "The latter will be IGNORED!!!"
       RETURN
    end if
    
    do iLat = 1, nLat
       do iLon = 1, nLon
          iNode = iNode_II(iLon, iLat)
          iBlock = iGridGlobal_IA(Block_, iNode)
          if(iProc == iGridGlobal_IA(Proc_, iNode))then
             call rlonlat_to_xyz(&
                  (/ROrigin, LonMin+(iLon-0.5)*DLon, LatMin+(iLat-0.5)*DLat/),&
                  State_VIB(X_:Z_,1,iBlock))
          end if
       end do
    end do
  end subroutine init_grid
  !============================================================================


  subroutine get_node_indexes(iNodeIn, iLonOut, iLatOut)
    ! return angular grid's indexes corresponding to this node
    integer, intent(in) :: iNodeIn
    integer, intent(out):: iLonOut
    integer, intent(out):: iLatOut
    !---------------------------------------------------------
    iLatOut = 1 + (iNodeIn-1) / nLon
    iLonOut = iNodeIn - nLon * (iLatOut-1)
  end subroutine get_node_indexes

  !============================================================================

  function distance_to_next(iParticle, iBlock) result(Distance)
    ! the function returns distance to the next particle measured in Rs;
    ! formula for distance between 2 points in rlonlat system:
    !  Distance**2 = R1**2 + R2**2 - 
    !     2*R1*R2*(Cos(Lat1)*Cos(Lat2) * Cos(Lon1-Lon2) + Sin(Lat1)*Sin(Lat2))
    ! NOTE: function doesn't check whether iParticle is last on the field line
    integer, intent(in):: iParticle
    integer, intent(in):: iBlock
    real               :: Distance
    !--------------------------------------------------------------------
    Distance = sqrt(sum((&
         State_VIB(X_:Z_, iParticle,   iBlock) - &
         State_VIB(X_:Z_, iParticle+1, iBlock))**2))
  end function distance_to_next

  !============================================================================

  subroutine append_particles
    !appends a new particle at the beginning of lines if necessary
    integer:: iBlock
    real:: DistanceToMin, Alpha
    real, parameter:: cTol = 1E-06

    character(len=*), parameter:: NameSub = 'append_particles'
    !--------------------------------------------------------------------
    do iBlock = 1, nBlock
       ! check current value of offset: if not zero, adjustments have just
       ! been made, no need to append new particles
       if(iGridLocal_IB(Offset_, iBlock) /= 0 )&
            CYCLE
       ! check if the beginning of the line moved far enough from its 
       ! footprint on the solar surface
       DistanceToMin = sqrt(sum((&
            State_VIB(X_:Z_,1,iBlock) - FootPoint_VB(X_:Z_,iBlock))**2))
       ! skip the line if it's still close to the Sun
       if(DistanceToMin * (1.0 + cTol) < FootPoint_VB(Length_, iBlock)) CYCLE
       ! append a new particle
       !-----------------------
       ! check if have enough space
       if(nParticleMax == nParticle_B( iBlock))&
            call CON_Stop(NameSub//&
            ': not enough memory allocated to append a new particle')
       ! shift the grid:
       State_VIB(     :,2:nParticle_B( iBlock)+1, iBlock) = &
            State_VIB(:,1:nParticle_B( iBlock),   iBlock)
       Distribution_IIB(     :,2:nParticle_B( iBlock)+1, iBlock) = &
            Distribution_IIB(:,1:nParticle_B( iBlock),   iBlock)
       nParticle_B( iBlock)  = nParticle_B( iBlock) + 1
       !Particles ID as handled by other components keep unchanged
       !while their order numbers in SP are increased by 1. Therefore,
       iGridLocal_IB(Offset_, iBlock)  = iGridLocal_IB(Offset_, iBlock) + 1
       ! put the new particle just above the lower boundary
       State_VIB(X_:Z_,  1, iBlock) = &
            FootPoint_VB(X_:Z_, iBlock) * (1.0 + cTol)
       State_VIB(LagrID_,1, iBlock) = State_VIB(LagrID_, 2, iBlock) - 1.0
       FootPoint_VB(LagrID_,iBlock) = State_VIB(LagrID_, 1, iBlock) - 1.0
       ! for old values of background parameters use extrapolation
       Alpha = DistanceToMin / (DistanceToMin + State_VIB(D_, 2, iBlock))
       State_VIB((/RhoOld_, BOld_/), 1, iBlock) = &
            (Alpha + 1)*State_VIB((/RhoOld_, BOld_/), 2, iBlock) &
            -Alpha     * State_VIB((/RhoOld_, BOld_/), 3, iBlock)
       Distribution_IIB(:,1,iBlock) = Distribution_IIB(:,2,iBlock) + &
            Alpha*(Distribution_IIB(:,2,iBlock) - Distribution_IIB(:,3,iBlock))
    end do
  end subroutine append_particles

end module SP_ModGrid
