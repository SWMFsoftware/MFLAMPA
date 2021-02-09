!  Copyright (C) 2002 Regents of the University of Michigan,
!  portions used with permission
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module SP_wrapper

  use SP_ModMain, ONLY: &
       run, save_restart, &
       DoRestart, DoReadMhData, &
       nDim, nLat, nLon, nBlock, nParticleMax, &
       RMin=>RScMin, RBufferMin=>RIhMin, &
       RBufferMax=>RScMax, RMax=>RIhMax, &
       MHData_VIB, iNode_B, FootPoint_VB, DataInputTime, &
       nParticle_B, Length_,&
       LagrID_,X_, Y_, Z_, Rho_, Bx_, Bz_,  Ux_, Uz_, T_, &
       Wave1_, Wave2_,  SI2IO_I, UnitEnergy_
  use CON_comp_info
  use CON_router, ONLY: IndexPtrType, WeightPtrType
  use CON_coupler, ONLY: &
       set_coord_system, SP_, is_proc0, i_comm, i_proc0, &
       init_decomposition, get_root_decomposition, bcast_decomposition, &
       iVar_V, DoCoupleVar_V, &
       Density_, RhoCouple_, Pressure_, PCouple_, &
       Momentum_, RhoUxCouple_, RhoUzCouple_, &
       BField_, BxCouple_, BzCouple_, &
       Wave_, WaveFirstCouple_, WaveLastCouple_
  use ModConst, ONLY: rSun, cProtonMass
  use ModMpi
  use CON_world, ONLY: is_proc0, is_proc, n_proc
  implicit none
  save
  private ! except
  public:: SP_set_param
  public:: SP_init_session
  public:: SP_run
  public:: SP_save_restart
  public:: SP_finalize

  ! coupling with MHD components
  public:: SP_put_from_mh
  public:: SP_interface_point_coords
  public:: SP_put_line
  public:: SP_adjust_lines
  public:: SP_get_bounds_comp
  public:: SP_n_particle
  public:: SP_check_ready_for_mh
  public:: SP_put_coupling_param
  ! variables requested via coupling: coordinates,
  ! field line and particles indexes
  character(len=*), parameter:: NameVarCouple =&
       'rho p mx my mz bx by bz i01 i02 pe'
  integer :: Model_ = -1
  integer, parameter:: Lower_=0, Upper_=1
  ! coupling parameters:
  ! domain boundaries
  real :: rInterfaceMin, rInterfaceMax
  ! buffer boundaries located near lower (Lo) or upper (Up) boudanry of domain
  real :: rBufferLo, rBufferUp
  integer, allocatable :: iOffset_B(:)
  logical :: DoCheck = .true.
contains
  !============================================================================
  ! Interface routines to be called from super-structure only
  subroutine SP_check_ready_for_mh(IsReady)
    use ModMpi
    logical, intent(out):: IsReady

    integer :: iError
    character(len=*), parameter:: NameSub = 'SP_check_ready_for_mh'
    !--------------------------------------------------------------------------
    ! when restarting, line data is available, i.e. ready to couple with mh;
    ! get value at SP root and broadcast to all SWMF processors
    if(is_proc0(SP_)) IsReady = DoRestart
    call MPI_Bcast(IsReady, 1, MPI_LOGICAL, i_proc0(SP_), i_comm(), iError)
  end subroutine SP_check_ready_for_mh
  !============================================================================
  subroutine SP_get_bounds_comp(ThisModel_, RMinOut, RMaxOut)
    use ModMpi
    ! return the MHD boundaries as set in SP component
    integer, intent(in )  :: ThisModel_
    real,    intent(out)  :: RMinOut, RMaxOut
    integer :: iError
    real    :: rAux_I(2)
    character(len=*), parameter:: NameSub = 'SP_get_bounds_comp'
    !--------------------------------------------------------------------------
    if(is_proc0(SP_))then
       select case(ThisModel_)
       case(Lower_)
          rAux_I(1) = RMin
          rAux_I(2) = RBufferMax
       case(Upper_)
          rAux_I(1) = RBufferMin
          rAux_I(2) = RMax
       case default
          call CON_stop('Incorrect model ID in '//NameSub)
       end select
    end if
    call MPI_Bcast(rAux_I(1), 2, MPI_REAL, i_proc0(SP_), i_comm(), iError)
       RMinOut = rAux_I(1); RMaxOut = rAux_I(2)
  end subroutine SP_get_bounds_comp
  !============================================================================
  ! Above routines may be called from superstructure only.
  integer function SP_n_particle(iBlockLocal)
    integer, intent(in) :: iBlockLocal
    !--------------------------------------------------------------------------
    SP_n_particle = nParticle_B(  iBlockLocal)
  end function SP_n_particle
  !============================================================================
  subroutine SP_run(TimeSimulation,TimeSimulationLimit)
    real,intent(inout)::TimeSimulation
    real,intent(in)::TimeSimulationLimit

    !--------------------------------------------------------------------------
    call run(TimeSimulationLimit)
    if(DoReadMhData)then
       TimeSimulation = DataInputTime
    else
       TimeSimulation = TimeSimulationLimit
    end if
  end subroutine SP_run
  !============================================================================
  subroutine SP_init_session(iSession,TimeSimulation)
    use SP_ModMain, ONLY: initialize
    use SP_ModGrid, ONLY: init_grid=>init
    use CON_mflampa, ONLY: set_state_pointer
    integer,  intent(in) :: iSession         ! session number (starting from 1)
    real,     intent(in) :: TimeSimulation   ! seconds from start time

    logical, save:: IsInitialized = .false.
    !--------------------------------------------------------------------------
    if(IsInitialized)&
         RETURN
    IsInitialized = .true.
    call init_grid
    nullify(MHData_VIB)
    call set_state_pointer(MHData_VIB, nBlock, nParticleMax)
    call initialize
    allocate(iOffset_B(nBlock)); iOffset_B = 0
  end subroutine SP_init_session
  !============================================================================
  subroutine SP_finalize(TimeSimulation)
    use SP_ModMain, ONLY: finalize
    real,intent(in)::TimeSimulation

    ! if data are read from files, no special finalization is needed
    !--------------------------------------------------------------------------
    if(.not.DoReadMhData)call run(TimeSimulation)
    call finalize
  end subroutine SP_finalize
  !============================================================================
  subroutine SP_set_param(CompInfo,TypeAction)
    use SP_ModTime,  ONLY: StartTimeJulian, StartTime, SPTime,&
         time_real_to_julian
    use SP_ModMain,  ONLY: check, read_param
    use SP_ModProc
    use CON_physics, ONLY: get_time
    type(CompInfoType),intent(inout):: CompInfo
    character(len=*),  intent(in)   :: TypeAction

    character(len=*), parameter:: NameSub = 'SP_set_param'
    !--------------------------------------------------------------------------
    select case(TypeAction)
    case('VERSION')
       call put(CompInfo,&
            Use        =.true., &
            NameVersion='MFLAMPA', &
            Version    =0.90)
    case('MPI')
       call get(CompInfo, iComm=iComm, iProc=iProc, nProc=nProc)
    case('STDOUT')
       ! placeholder
    case('CHECK')
       if(.not.DoCheck)RETURN
       DoCheck = .false.
       call get_time(tSimulationOut = SPTime, tStartOut = StartTime)
       call time_real_to_julian(StartTime, StartTimeJulian)
       DataInputTime = SPTime
       call check
    case('READ')
       call read_param
    case('GRID')
       call SP_set_grid
    case default
       call CON_stop('Cannot call SP_set_param for '//trim(TypeAction))
    end select
  end subroutine SP_set_param
  !============================================================================
  subroutine SP_save_restart(TimeSimulation)
    real, intent(in) :: TimeSimulation
    !--------------------------------------------------------------------------
    ! if data are read from files, no need for additional run
    if(.not.DoReadMhData)call run(TimeSimulation)
    call save_restart
  end subroutine SP_save_restart
  !============================================================================

  subroutine SP_put_from_mh(nPartial,iPutStart,Put,W,DoAdd,Buff_I,nVar)
    integer,intent(in)::nPartial,iPutStart,nVar
    type(IndexPtrType),intent(in)::Put
    type(WeightPtrType),intent(in)::W
    logical,intent(in)::DoAdd
    real,dimension(nVar),intent(in)::Buff_I
    integer:: iRho, iP, iMx, iMz, iBx, iBz, iWave1, iWave2
    integer:: i, j, k, iBlock
    integer:: iPartial
    real:: Weight
    real:: R, Aux

    character(len=100):: StringError
    character(len=*), parameter:: NameSub = 'SP_put_from_mh'
    !--------------------------------------------------------------------------
    ! check consistency of DoCoupleVar_V
    if(.not. DoCoupleVar_V(Density_) .and. &
         (DoCoupleVar_V(Pressure_) .or. DoCoupleVar_V(Momentum_)))&
         call CON_Stop(NameSub//': pressure or momentum is coupled,'&
         //' but density is not')
    ! indices of variables in the buffer
    iRho  = iVar_V(RhoCouple_)
    iP    = iVar_V(PCouple_)
    iMx   = iVar_V(RhoUxCouple_)
    iMz   = iVar_V(RhoUzCouple_)
    iBx   = iVar_V(BxCouple_)
    iBz   = iVar_V(BzCouple_)
    iWave1= iVar_V(WaveFirstCouple_)
    iWave2= iVar_V(WaveLastCouple_)
    Aux = 0
    if(DoAdd)Aux = 1.0
    do iPartial = 0, nPartial-1
       ! cell and block indices
       i      = Put%iCB_II(1, iPutStart + iPartial)
       j      = Put%iCB_II(2, iPutStart + iPartial)
       k      = Put%iCB_II(3, iPutStart + iPartial)
       iBlock = Put%iCB_II(4, iPutStart + iPartial)
       ! interpolation weight
       Weight = W%Weight_I(   iPutStart + iPartial)
       if(is_in_buffer_xyz(Lower_,MHData_VIB(X_:Z_,i,iBlock)))then
          R = sqrt(sum(MHData_VIB(X_:Z_,i,iBlock)**2))
          Aux = 1.0
          Weight = Weight * (0.50 + 0.50*tanh(2*(2*R - &
               RBufferLo - RInterfaceMin)/(RBufferLo - RInterfaceMin)))
       end if
       if(is_in_buffer_xyz(Upper_,MHData_VIB(X_:Z_,i,iBlock)))then
          R = sqrt(sum(MHData_VIB(X_:Z_,i,iBlock)**2))
          Weight = Weight * (0.50 - 0.50*tanh(2*(2*R - &
               RInterfaceMax - RBufferUp)/(RInterfaceMax - RBufferUp)))
       end if
       ! put the data
       if(DoCoupleVar_V(Density_))&
            MHData_VIB(Rho_,i,iBlock) = Aux*MHData_VIB(Rho_,i,iBlock) &
            + Buff_I(iRho)/cProtonMass*Weight
       if(DoCoupleVar_V(Pressure_))&
            MHData_VIB(T_,i,iBlock) = Aux*MHData_VIB(T_,i,iBlock) + &
            Buff_I(iP)/Buff_I(iRho)*cProtonMass*SI2IO_I(UnitEnergy_)*Weight
       if(DoCoupleVar_V(Momentum_))&
            MHData_VIB(Ux_:Uz_,i,iBlock) = Aux*MHData_VIB(Ux_:Uz_,i,iBlock) + &
            Buff_I(iMx:iMz) / Buff_I(iRho) * Weight
       if(DoCoupleVar_V(BField_))&
            MHData_VIB(Bx_:Bz_,i,iBlock) = Aux*MHData_VIB(Bx_:Bz_,i,iBlock) + &
            Buff_I(iBx:iBz) * Weight
       if(DoCoupleVar_V(Wave_))&
            MHData_VIB(Wave1_:Wave2_,i,iBlock) = &
            Aux*MHData_VIB(Wave1_:Wave2_,i,iBlock) + &
            Buff_I(iWave1:iWave2)*Weight
    end do
  end subroutine SP_put_from_mh
  !============================================================================
  subroutine SP_set_grid
    use SP_ModGrid, ONLY: iGridGlobal_IA, Block_, Proc_, &
          TypeCoordSystem
    logical, save:: IsInitialized = .false.
    !--------------------------------------------------------------------------
    if(IsInitialized)RETURN
    IsInitialized = .true.
    ! Initialize 3D grid with NON-TREE structure
    call init_decomposition(&
         GridID_ = SP_,&
         CompID_ = SP_,&
         nDim    = nDim)
    ! Construct decomposition
    if(is_proc0(SP_))&
         call get_root_decomposition(&
         GridID_       = SP_,&
         iRootMapDim_D = [1, nLon, nLat],&
         CoordMin_D    = [0.50, 0.50, 0.50],&
         CoordMax_D    = [nParticleMax, nLon, nLat] + 0.50,&
         nCells_D      = [nParticleMax, 1, 1],&
         PE_I          = iGridGlobal_IA(Proc_,:),&
         iBlock_I      = iGridGlobal_IA(Block_,:))
    call bcast_decomposition(SP_)
    ! Coordinate system is Heliographic Inertial Coordinate System (HGI)
    ! with length measured in solar radii
    call set_coord_system(&
         GridID_      = SP_, &
         TypeCoord    = TypeCoordSystem, &
         TypeGeometry = 'cartesian', &
         NameVar      = NameVarCouple, &
         UnitX        = rSun)
  end subroutine SP_set_grid
  !============================================================================
  subroutine SP_put_coupling_param(iModelIn, rMinIn, rMaxIn, TimeIn,&
       rBufferLoIn, rBufferUpIn)
    use SP_ModMain, ONLY: copy_old_state
    integer,        intent(in) :: iModelIn
    real,           intent(in) :: rMinIn, rMaxIn
    real,           intent(in) :: TimeIn
    real, optional, intent(in) :: rBufferLoIn, rBufferUpIn
    !--------------------------------------------------------------------------
    ! set domain boundaries
    rInterfaceMin = rMinIn; rInterfaceMax = rMaxIn
    ! set buffer boundaries
    if(present(rBufferLoIn))then
       rBufferLo = rBufferLoIn
    else
       rBufferLo = rMinIn
    end if
    if(present(rBufferUpIn))then
       rBufferUp = rBufferUpIn
    else
       rBufferUp = rMaxIn
    end if
    Model_ = iModelIn
    if(DataInputTime >= TimeIn)RETURN
    ! New coupling time, get it and save old state
    DataInputTime = TimeIn
    if(Model_==Lower_)then
       call copy_old_state
    else
       call CON_stop("Time in IHSP coupling differs from that in SCSP")
    end if
  end subroutine SP_put_coupling_param
  !============================================================================
  subroutine SP_interface_point_coords(nDim, Xyz_D, &
       nIndex, iIndex_I, IsInterfacePoint)
    ! interface points (request), which needed to be communicated
    ! to other components to perform field line extraction and
    ! perform further coupling with SP:
    ! the framework tries to determine Xyz_D of such points,
    ! SP changes them to the correct values
    integer,intent(in)   :: nDim
    real,   intent(inout):: Xyz_D(nDim)
    integer,intent(in)   :: nIndex
    integer,intent(inout):: iIndex_I(nIndex)
    logical,intent(out)  :: IsInterfacePoint
    integer:: iParticle, iBlock
    real:: R2
    character(len=*), parameter:: NameSub = 'SP_interface_point_coords'
    !--------------------------------------------------------------------------
    iParticle = iIndex_I(1); iBlock    = iIndex_I(4)
    ! Check whether the particle is within interface bounds
    R2 = sum(MHData_VIB(X_:Z_,iParticle,iBlock)**2)
    IsInterfacePoint = &
         R2 >= rInterfaceMin**2 .and. R2 < rInterfaceMax**2
    ! Fix coordinates to be used in mapping
    if(IsInterfacePoint)&
         Xyz_D = MHData_VIB(X_:Z_, iParticle, iBlock)
  end subroutine SP_interface_point_coords
  !============================================================================
  subroutine SP_put_line(nPartial, iPutStart, Put,&
       Weight, DoAdd, Coord_D, nVar)
    integer, intent(in) :: nPartial, iPutStart, nVar
    type(IndexPtrType), intent(in) :: Put
    type(WeightPtrType),intent(in) :: Weight
    logical,            intent(in) :: DoAdd
    real,               intent(in) :: Coord_D(nVar) ! nVar=nDim

    ! indices of the particle
    integer:: iBlock, iParticle
    ! Misc
    real :: R2
    character(len=*), parameter:: NameSub = 'SP_put_line'
    !--------------------------------------------------------------------------
    R2 = sum(Coord_D(1:nDim)**2)
    ! Sort out particles out of the SP domain
    if(R2<RMin**2.or.R2>=RMax**2)RETURN
    ! store passed particles
    iBlock    = Put%iCB_II(4,iPutStart)
    iParticle = Put%iCB_II(1,iPutStart) + iOffset_B(iBlock)
    ! put coordinates
    MHData_VIB(X_:Z_,iParticle, iBlock) = Coord_D(1:nDim)
    nParticle_B(iBlock) = MAX(nParticle_B(iBlock), iParticle)
  end subroutine SP_put_line
  !============================================================================
  ! Called from coupler after the updated grid point location are
  ! received from the other component (SC, IH). Determines whether some
  ! grid points should be added/deleted
  subroutine SP_adjust_lines(DoInit)
    use SP_ModDistribution, ONLY: offset
    use SP_ModGrid,         ONLY: R_, State_VIB
    logical, intent(in) :: DoInit
    ! once new geometry of lines has been put, account for some particles
    ! exiting the domain (can happen both at the beginning and the end)
    integer:: iParticle, iBlock, iBegin,  iEnd, iOffset ! loop variables
    integer:: iParticle_I(2), iLoop
    logical:: DoAdjustLo, DoAdjustUp
    logical:: IsMissing

    integer, parameter:: Lo_ = 1, Up_ = 2
    integer, parameter:: iIncrement_II(2,2) =reshape([1,0,0,-1],[2,2])

    character(len=100):: StringError

    character(len=*), parameter:: NameSub = 'SP_adjust_lines'
    !--------------------------------------------------------------------------
    DoAdjustLo = RBufferLo == RInterfaceMin
    DoAdjustUp = RBufferUp == RInterfaceMax
    if(DoInit.and.DoAdjustLo)then
       do iBlock = 1, nBlock
          call SP_set_line_foot_b(iBlock)
       end do
    end if
    BLOCK:do iBlock = 1, nBlock
       ! Called after the grid points are received from the
       ! component, nullify offset
       if(DoAdjustLo)iOffset_B(iBlock) = 0
       iBegin = 1
       iEnd   = nParticle_B(  iBlock)
       iParticle_I(:) = [iBegin, iEnd]
       if(DoAdjustUp) then
          iLoop = Up_
       else
          iLoop = Lo_
       end if
       PARTICLE: do while(iParticle_I(1) < iParticle_I(2))
          iParticle = iParticle_I(iLoop)
          ! account for all missing partiles along the line;
          ! --------------------------------------------------------
          ! particle import MUST be performed in order from lower to upper
          ! components (w/respect to radius), e.g. SC, then IH, then OH
          ! adjustment are made after each import;
          ! --------------------------------------------------------
          ! lines are assumed to start in lowest model
          ! lowest model should NOT have the Lo buffer,
          ! while the upper most should NOT have Up buffer,
          ! buffer are REQUIRED between models
          ! --------------------------------------------------------
          ! adjust as follows:
          ! LOWEST model may loose particles at the start of lines
          ! and (if line ends in the model) in the tail;
          ! loop over particles Lo-2-Up and mark losses (increase iBegin)
          ! until particles reach UP buffer;
          ! if line ends in the model, loop over particles Up-2-Lo
          ! and mark losses (decrease nParticle_B) until UP buffer is reached;
          ! after that, if line reenters the model,
          ! particles within the model may not be lost,
          !
          ! MIDDLE models are not allowed to loose particles
          !
          ! HIGHEST model may only loose particles in the tail;
          ! loop Up-2-Lo and mark losses  (decrease nParticle_B)
          ! until the bottom of Lo buffer is reaced
          ! --------------------------------------------------------
          ! whenever a particle is lost in lower models -> ERROR

          ! when looping Up-2-Lo and particle is in other model ->
          ! adjustments are no longer allowed
          if(iLoop == Up_ .and. (&
               State_VIB(R_,iParticle,iBlock) <  RInterfaceMin&
               .or.&
               State_VIB(R_,iParticle,iBlock) >= RInterfaceMax)&
               )then
             DoAdjustLo = .false.
             DoAdjustUp = .false.
          end if

          ! determine whether particle is missing
          IsMissing = all(MHData_VIB(X_:Z_,iParticle,iBlock)==0.0)
          if(.not.IsMissing) then
             iParticle_I = iParticle_I + iIncrement_II(:,iLoop)
             CYCLE PARTICLE
          end if

          ! missing point in the lower part of the domain -> ERROR
          if(State_VIB(R_,iParticle,iBlock) < RInterfaceMin)&
               call CON_stop(NameSub//": particle has been lost")

          ! missing point in the upper part of the domain -> IGNORE;
          ! if needed to adjust beginning, then it is done,
          ! switch left -> right end of range and start adjusting
          ! tail of the line, if it has reentered current part of the domain
          if(State_VIB(R_,iParticle,iBlock) >= RInterfaceMax)then
             if(iLoop == Lo_)&
                  iLoop = Up_
             iParticle_I = iParticle_I + iIncrement_II(:,iLoop)
             if(DoAdjustLo)then
                DoAdjustLo = .false.
                DoAdjustUp = .true.
             end if
             CYCLE PARTICLE
          end if

          ! if point used to be in a upper buffer -> IGNORE
          if(is_in_buffer_r(Upper_, State_VIB(R_, iParticle, iBlock)))then
             iParticle_I = iParticle_I + iIncrement_II(:,iLoop)
             CYCLE PARTICLE
          end if

          ! if need to adjust lower, but not upper boundary -> ADJUST
          if(DoAdjustLo .and. .not.DoAdjustUp)then
             ! push iBegin in front of current particle;
             ! it will be pushed until it finds a non-missing particle
             iBegin = iParticle + 1
             iParticle_I = iParticle_I + iIncrement_II(:,iLoop)
             CYCLE PARTICLE
          end if

          ! if need to adjust upper, but not lower boundary -> ADJUST
          if(DoAdjustUp .and. .not.DoAdjustLo)then
             ! push nParticle_B() below current particle;
             ! it will be pushed until it finds a non-missing particle
             nParticle_B(iBlock) = iParticle - 1
             iParticle_I = iParticle_I + iIncrement_II(:,iLoop)
             CYCLE PARTICLE
          end if

          ! remaining case:
          ! need to adjust both boudnaries -> ADJUST,but keep longest range
          if(iParticle - iBegin > nParticle_B(iBlock) - iParticle)then
             nParticle_B(iBlock) = iParticle - 1
             EXIT PARTICLE
          else
             iBegin = iParticle + 1
          end if
          iParticle_I = iParticle_I + iIncrement_II(:,iLoop)
       end do PARTICLE

       DoAdjustLo = RBufferLo == RInterfaceMin
       DoAdjustUp = RBufferUp == RInterfaceMax

       if(iBegin/=1)then
          ! Offset particle arrays
          iEnd   = nParticle_B(iBlock)
          iOffset = 1 - iBegin
          iOffset_B(iBlock) = iOffset
          MHData_VIB(LagrID_:Z_, 1:iEnd+iOffset, iBlock) = &
               MHData_VIB(LagrID_:Z_,iBegin:iEnd,iBlock)
          nParticle_B(iBlock) = nParticle_B(iBlock) + iOffset
          ! need to recalculate footpoints
          call SP_set_line_foot_b(iBlock)
          call offset(iBlock, iOffset)
       end if
    end do BLOCK
    ! may need to add particles to the beginning of lines
    if(DoAdjustLo) call append_particles
    ! Called after the grid points are received from the
    ! component, nullify offset. Alternatively, if the points
    ! are received from SC via initial coupling, there is no
    ! need to apply offset in IH, because there are no points
    ! in IH yet and no need to correct their ID
    if(DoAdjustUp.or.DoInit)iOffset_B(1:nBlock) = 0
  contains
    !==========================================================================
    subroutine SP_set_line_foot_b(iBlock)
      integer, intent(in) :: iBlock

      ! existing particle with lowest index along line
      real:: Xyz1_D(nDim)
      ! direction of the field at Xyz1_D and segment vectors between particles
      real, dimension(nDim):: Dir0_D, Dist1_D, Dist2_D
      ! dot product Xyz1 and Dir1 and its sign
      real:: Dot, S
      ! distances between particles
      real:: Dist1, Dist2
      ! variable to compute coords of the footprints
      real:: Alpha
      !------------------------------------------------------------------------
      ! get the coordinates of lower particle
      Xyz1_D = MHData_VIB(X_:Z_, 1, iBlock)

      ! generally, field direction isn't known
      ! approximate it using directions of first 2 segments of the line
      Dist1_D = MHData_VIB(X_:Z_, 1, iBlock) - &
           MHData_VIB(X_:Z_, 2, iBlock)
      Dist1 = sqrt(sum(Dist1_D**2))
      Dist2_D = MHData_VIB(X_:Z_, 2, iBlock) - &
           MHData_VIB(X_:Z_, 3, iBlock)
      Dist2 = sqrt(sum(Dist2_D**2))
      Dir0_D = ((2*Dist1 + Dist2)*Dist1_D - Dist1*Dist2_D)/(Dist1 + Dist2)

      Dir0_D = Dir0_D/sqrt(sum(Dir0_D**2))

      ! dot product and sign: used in computation below
      Dot = sum(Dir0_D*Xyz1_D)
      S   = sign(1.0, Dot)

      ! there are 2 possible failures of the algorithm:
      ! Failure (1):
      ! no intersection of smoothly extended line with the sphere R = RMin
      if(Dot**2 - sum(Xyz1_D**2) + RMin**2 < 0)then
         ! project first particle for new footpoint
         FootPoint_VB(X_:Z_,iBlock) = Xyz1_D * RMin / sqrt(sum(Xyz1_D**2))
      else
         ! Xyz0, the footprint, is distance Alpha away from Xyz1:
         ! Xyz0 = Xyz1 + Alpha * Dir0 and R0 = RMin =>
         Alpha = S * sqrt(Dot**2 - sum(Xyz1_D**2) + RMin**2) - Dot
         ! Failure (2):
         ! intersection is too far from the current beginning of the line,
         ! use distance between 2nd and 3rd particles on the line as measure
         if(abs(Alpha) > Dist2)then
            ! project first particle for new footpoint
            FootPoint_VB(X_:Z_,iBlock) = Xyz1_D * RMin / sqrt(sum(Xyz1_D**2))
         else
            ! store newly found footpoint of the line
            FootPoint_VB(X_:Z_,iBlock) = Xyz1_D + Alpha * Dir0_D
         end if
      end if

      ! length is used to decide when need to append new particles:
      ! use distance between 2nd and 3rd particles on the line
      FootPoint_VB(Length_,    iBlock) = Dist2
      FootPoint_VB(LagrID_,    iBlock) = MHData_VIB(LagrID_,1,iBlock) - 1.0
    end subroutine SP_set_line_foot_b
    !==========================================================================
    subroutine append_particles
      use SP_ModGrid, ONLY: R_, State_VIB
      ! appends a new particle at the beginning of lines if necessary
      integer:: iBlock
      real:: DistanceToMin
      real, parameter:: cTol = 1E-06

      character(len=*), parameter:: NameSub = 'append_particles'
      !------------------------------------------------------------------------
      BLOCK:do iBlock = 1, nBlock
         ! check current value of offset: if not zero, adjustments have just
         ! been made, no need to append new particles
         if(iOffset_B(iBlock) /= 0 )CYCLE BLOCK
         ! check if the beginning of the line moved far enough from its
         ! footprint on the solar surface
         DistanceToMin = sqrt(sum((&
              MHData_VIB(X_:Z_,1,iBlock) - FootPoint_VB(X_:Z_,iBlock))**2))
         ! skip the line if it's still close to the Sun
         if(DistanceToMin*(1.0 + cTol) < FootPoint_VB(Length_, iBlock))&
              CYCLE BLOCK
         ! append a new particle
         ! check if have enough space
         if(nParticleMax == nParticle_B( iBlock))call CON_Stop(NameSub//&
              ': not enough memory allocated to append a new particle')
         ! Particles ID as handled by other components keep unchanged
         ! while their order numbers in SP are increased by 1. Therefore,
         iOffset_B(iBlock)  = 1
         MHData_VIB(       LagrID_:Z_,2:nParticle_B(iBlock) + 1, iBlock)&
              = MHData_VIB(LagrID_:Z_,1:nParticle_B(iBlock),     iBlock)
         State_VIB(       R_        ,2:nParticle_B(iBlock) + 1, iBlock)&
              = State_VIB(R_        ,1:nParticle_B(iBlock),     iBlock)
         nParticle_B(iBlock) = nParticle_B(iBlock) + 1
         ! put the new particle just above the lower boundary
         MHData_VIB(LagrID_:Z_,  1, iBlock) = &
              FootPoint_VB(LagrID_:Z_, iBlock)*(1.0 + cTol)
         State_VIB(R_,          1, iBlock) = &
              sqrt(sum((MHData_VIB(X_:Z_,  1, iBlock))**2))
         MHData_VIB(LagrID_,1, iBlock) = MHData_VIB(LagrID_, 2, iBlock) - 1.0
         FootPoint_VB(LagrID_,iBlock) = MHData_VIB(LagrID_, 1, iBlock) - 1.0
         call offset(iBlock, iOffset=iOffset_B(iBlock))
      end do BLOCK
    end subroutine append_particles
    !==========================================================================
  end subroutine SP_adjust_lines
  !============================================================================
  function is_in_buffer_r(iBuffer, R) Result(IsInBuffer)
    integer,intent(in) :: iBuffer
    real,   intent(in) :: R
    logical:: IsInBuffer
    !--------------------------------------------------------------------------
    select case(iBuffer)
    case(Lower_)
       IsInBuffer = R >= rInterfaceMin .and. R < rBufferLo
    case(Upper_)
       IsInBuffer = R >= rBufferUp .and. R < rInterfaceMax
    case default
        call CON_stop("ERROR: incorrect call of SP_wrapper:is_in_buffer")
    end select
  end function is_in_buffer_r
  !============================================================================
  function is_in_buffer_xyz(iBuffer, Xyz_D) Result(IsInBuffer)
    integer,intent(in):: iBuffer
    real,   intent(in) :: Xyz_D(nDim)
    logical:: IsInBuffer
    real:: R
    !--------------------------------------------------------------------------
    R = sqrt(sum(Xyz_D**2))
    IsInBuffer = is_in_buffer_r(iBuffer, R)
  end function is_in_buffer_xyz
  !============================================================================
end module SP_wrapper
!==============================================================================
