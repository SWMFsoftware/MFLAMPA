!  Copyright (C) 2002 Regents of the University of Michigan, 
!  portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!=============================================================!
module SP_wrapper

  use ModNumConst, ONLY: cHalfPi
  use ModConst, ONLY: rSun, cProtonMass
  use ModCoordTransform, ONLY: xyz_to_rlonlat, rlonlat_to_xyz
  use ModMain, ONLY: &
       run, initialize, finalize, check, read_param,&
       get_node_indexes, &
       iComm, iProc, nProc, &
       nDim, nNode, nLat, nLon, nBlock,&
       iParticleMin, iParticleMax, nParticle,&
       RSc, LatMin, LatMax, LonMin, LonMax, &
       iGridGlobal_IA, iGridLocal_IB, State_VIB, iNode_B, TypeCoordSystem,&
       Block_, Proc_, Begin_, End_, &
       R_, Lat_, Lon_, Rho_, Bx_,By_,Bz_,B_, Ux_,Uy_,Uz_, T_, RhoOld_, BOld_
  use ModBufferQueue, ONLY: TypeBufferQueue, &
       init_queue, reset_queue, add_to_queue, reset_peek_queue, peek_queue
  use CON_comp_info
  use CON_router, ONLY: IndexPtrType, WeightPtrType
  use CON_coupler, ONLY: &
       set_coord_system, &
       init_decomposition, get_root_decomposition, bcast_decomposition, &
       iVar_V, DoCoupleVar_V, &
       Density_, RhoCouple_, Pressure_, PCouple_, &
       Momentum_, RhoUxCouple_, RhoUzCouple_, &
       BField_, BxCouple_, BzCouple_
  use CON_world, ONLY: is_proc0
  use CON_comp_param, ONLY: SP_

  implicit none

  save

  private ! except

  public:: SP_set_param
  public:: SP_init_session
  public:: SP_run
  public:: SP_save_restart
  public:: SP_finalize

  ! coupling with MHD components
  public:: SP_put_input_time
  public:: SP_put_from_mh
  public:: SP_get_request
  public:: SP_put_line
  public:: SP_get_grid_descriptor_param
  public:: SP_get_line_all
  public:: SP_get_solar_corona_boundary

  ! variables requested via coupling: coordinates, 
  ! field line and particles indexes
  character(len=*), parameter:: NameVarCouple = 'rho p mx my mz bx by bz'

  ! particles that need to be requested from MH component
  type(TypeBufferQueue):: QueueRequest
  ! indices that define a request: field line and particle indices
  integer, parameter  :: nRequestIndex = 2
  ! estimated number of requests per field line
  integer, parameter  :: nRequestPerLine = nParticle / 100

contains

  subroutine SP_run(TimeSimulation,TimeSimulationLimit)
    real,intent(inout)::TimeSimulation
    real,intent(in)::TimeSimulationLimit
    !--------------------------------------------------------------------------
    call run(TimeSimulationLimit)
    TimeSimulation = TimeSimulationLimit
  end subroutine SP_run

  !========================================================================

  subroutine SP_init_session(iSession,TimeSimulation)


    integer,  intent(in) :: iSession         ! session number (starting from 1)
    real,     intent(in) :: TimeSimulation   ! seconds from start time

    logical, save:: IsInitialized = .false.
    !--------------------------------------------------------------------------
    if(IsInitialized)&
         RETURN
    IsInitialized = .true.
     call initialize(TimeSimulation)
  end subroutine SP_init_session

  !======================================================================

  subroutine SP_finalize(TimeSimulation)


    real,intent(in)::TimeSimulation
    !--------------------------------------------------------------------------
    call finalize
  end subroutine SP_finalize

  !=========================================================

  subroutine SP_set_param(CompInfo,TypeAction)

    type(CompInfoType),intent(inout):: CompInfo
    character(len=*),  intent(in)   :: TypeAction

    character(len=*), parameter :: NameSub='SP_set_param'
    !-------------------------------------------------------------------------
    select case(TypeAction)
    case('VERSION')
       call put(CompInfo,&
            Use        =.true., &
            NameVersion='Empty', &
            Version    =0.0)
    case('MPI')
       call get(CompInfo, iComm=iComm, iProc=iProc, nProc=nProc)
    case('STDOUT')
       ! placeholder
    case('CHECK')
       call check
    case('READ')
       call read_param(TypeAction)
    case('GRID')
       call SP_set_grid
    case default
       call CON_stop('Can not call SP_set_param for '//trim(TypeAction))
    end select
  end subroutine SP_set_param

  !=========================================================

  subroutine SP_save_restart(TimeSimulation) 

    real,     intent(in) :: TimeSimulation 
    call CON_stop('Can not call SP_save restart')
  end subroutine SP_save_restart

  !=========================================================

  subroutine SP_put_input_time(TimeIn)

    real,     intent(in)::TimeIn
    call CON_stop('Can not call SP_get_input_time')
  end subroutine SP_put_input_time

  !===================================================================

  subroutine SP_put_from_mh(nPartial,iPutStart,Put,W,DoAdd,Buff_I,nVar)

    integer,intent(in)::nPartial,iPutStart,nVar
    type(IndexPtrType),intent(in)::Put
    type(WeightPtrType),intent(in)::W
    logical,intent(in)::DoAdd
    real,dimension(nVar),intent(in)::Buff_I
    integer:: iRho, iP, iMx, iMz, iBx, iBz
    integer:: i, j, k, iBlock
    integer:: iPartial
    real:: Weight
    real:: Aux

    real, external:: energy_in

    character(len=*), parameter:: NameSub='SP_put_from_mh'
    !------------------------------------------------------------
    ! check consistency: momentum and pressure are needed together with density
    if(.not. DoCoupleVar_V(Density_) .and. &
         (DoCoupleVar_V(Pressure_) .or. DoCoupleVar_V(Momentum_)))&
         call CON_Stop(NameSub//': pressure or momentum is coupled,'//&
         ' but density is not')

    ! indices of variables in the buffer
    iRho= iVar_V(RhoCouple_)
    iP  = iVar_V(PCouple_)
    iMx = iVar_V(RhoUxCouple_)
    iMz = iVar_V(RhoUzCouple_)
    iBx = iVar_V(BxCouple_)
    iBz = iVar_V(BzCouple_)   
    ! auxilary factor to account for value of DoAdd
    Aux = 0.0
    if(DoAdd) Aux = 1.0

    do iPartial = 0, nPartial-1
       ! cell and block indices
       i      = Put%iCB_II(1, iPutStart + iPartial)
       j      = Put%iCB_II(2, iPutStart + iPartial)
       k      = Put%iCB_II(3, iPutStart + iPartial)
       iBlock = Put%iCB_II(4, iPutStart + iPartial)
       ! interpolation weight
       Weight = W%Weight_I(   iPutStart + iPartial)
       ! put the data
       ! NOTE: State_VIB must be reset to zero before putting coupled data
       if(DoCoupleVar_V(Density_))&
            State_VIB(Rho_,i,iBlock) = Aux * State_VIB(Rho_,i,iBlock) + &
            Buff_I(iRho)/cProtonMass * Weight
       if(DoCoupleVar_V(Pressure_))&
            State_VIB(T_,i,iBlock) = Aux * State_VIB(T_,i,iBlock) + &
            Buff_I(iP)/Buff_I(iRho)*cProtonMass/energy_in('kev') * Weight
       if(DoCoupleVar_V(Momentum_))&
            State_VIB(Ux_:Uz_,i,iBlock) = Aux * State_VIB(Ux_:Uz_,i,iBlock) + &
            Buff_I(iMx:iMz) / Buff_I(iRho) * Weight
       if(DoCoupleVar_V(BField_))&
            State_VIB(Bx_:Bz_,i,iBlock) = Aux * State_VIB(Bx_:Bz_,i,iBlock) + &
            Buff_I(iBx:iBz) * Weight
    end do
  end subroutine SP_put_from_mh

  !===================================================================

  subroutine SP_set_grid

    ! Initialize 3D grid with NON-TREE structure
    call init_decomposition(&
         GridID_ = SP_,&
         CompID_ = SP_,&
         nDim    = nDim)

    ! Construct decomposition
    if(is_proc0(SP_))&
         call get_root_decomposition(&
         GridID_       = SP_,&
         iRootMapDim_D = (/1, nLat, nLon/),&
         XyzMin_D      = (/real(iParticleMin), LatMin, LonMin/),&
         XyzMax_D      = (/real(iParticleMax), LatMax, LonMax/),&
         nCells_D      = (/nParticle , 1, 1/),&
         PE_I          = iGridGlobal_IA(Proc_,:),&
         iBlock_I      = iGridGlobal_IA(Block_,:))
    call bcast_decomposition(SP_)

    ! Coordinate system is Heliographic Inertial Coordinate System (HGI)
    ! with length measured in solar radii
    call set_coord_system(&
         GridID_      = SP_, &
         TypeCoord    = TypeCoordSystem, &
         TypeGeometry = 'spherical', &
         NameVar      = NameVarCouple, &
         UnitX        = rSun)
  end subroutine SP_set_grid

  !===================================================================

  subroutine SP_get_solar_corona_boundary(RScOut)
    ! return the value of the solar corona boundary as set in SP component
    real, intent(out):: RScOut
    !-----------------------------------------------------------------
    RScOut = RSc
  end subroutine SP_get_solar_corona_boundary

  !===================================================================
  subroutine SP_get_request( &
       nRequestOut, nCoordOut, CoordOut_DI, iIndexOut_II, nAux, AuxOut_VI)
    ! request coordinates & indices of field lines' beginning/origin/end
    ! for the current processor
    !---------------------------------------------------------------
    integer,              intent(out):: nRequestOut
    integer,              intent(out):: nCoordOut
    real,    allocatable, intent(out):: CoordOut_DI(:, :)
    integer, allocatable, intent(out):: iIndexOut_II(:,:)
    integer,              intent(out):: nAux
    real,    allocatable, intent(out):: AuxOut_VI(:,:)

    ! loop variables
    integer:: iParticle, iBlock, iNode, iRequest
    integer:: iRequest_I(nRequestIndex)

    logical, save:: IsFirstCall = .true.
    character(len=*), parameter:: NameSub='SP_get_request'
    !----------------------------------------------------------------
    if(IsFirstCall)then
       IsFirstCall = .false.
       ! need to initalize request buffer:
       ! size of request buffer is estimated average number of requests 
       ! per line times number of lines on this proc
       call init_queue(QueueRequest, nRequestPerLine * nBlock, nRequestIndex)
       ! the initial request contains the origin points only
       do iBlock = 1, nBlock
          call add_to_queue(QueueRequest, (/iBlock, 0/))
       end do
    end if

    ! size of the request
    nRequestOut = QueueRequest % nRecordAll
    nCoordOut   = nDim
    nAux        = 2

    ! prepare containers to hold the request
    if(allocated(CoordOut_DI)) deallocate(CoordOut_DI)
    allocate(CoordOut_DI(nDim, nRequestOut))
    if(allocated(iIndexOut_II)) deallocate(iIndexOut_II)
    allocate(iIndexOut_II(nDim+1, nRequestOut))! 3 cell + 1 block index
    if(allocated(AuxOut_VI)) deallocate(AuxOut_VI)
    allocate(AuxOut_VI(nAux, nRequestOut))
    
    ! go over the request queue
    call reset_peek_queue(QueueRequest)
    do iRequest = 1, nRequestOut
       call peek_queue(QueueRequest, iRequest_I)
       iBlock    = iRequest_I(1)
       iParticle = iRequest_I(2)
       iNode     = iNode_B(iBlock)
       CoordOut_DI(:, iRequest) = &
            State_VIB((/R_,Lon_,Lat_/), iParticle, iBlock)
       iIndexOut_II(1, iRequest) = iParticle
       call get_node_indexes(iNode, &
            iIndexOut_II(2, iRequest), iIndexOut_II(3, iRequest))
       iIndexOut_II(4, iRequest) = iBlock
       AuxOut_VI(1, iRequest) = real(iNode)
       AuxOut_VI(2, iRequest) = real(iParticle)
    end do

    ! request is complete => reset queue
    call reset_queue(QueueRequest)

  end subroutine SP_get_request

  !===================================================================

  subroutine SP_put_line(nPut, Coord_DI, iIndex_II)
    use ModMpi
    ! store particle coordinates extracted elsewhere
    !---------------------------------------------------------------
    integer, intent(in):: nPut
    real,    intent(in):: Coord_DI( nDim,   nPut)
    integer, intent(in):: iIndex_II(nDim+1, nPut)

    ! cartesian coordinates
    real:: Xyz_D(nDim)
    ! radius-lon-lat coordinates
    real:: Coord_D(nDim)
    ! loop variables
    integer:: iPut, iBlock, iNode
    ! indices of the particle
    integer:: iLine, iParticle
    integer:: iMin_A(nNode),iMax_A(nNode)
    integer:: iError
    ! which field lines are being extracted
    logical:: WasInSc, IsInSc, WasUndef
    logical, save:: IsFirstCall = .true.
    logical, save,allocatable:: DoneExtractSolarCorona_B(:)
    integer, parameter:: nVarReset  = 8
    integer, parameter:: &
         VarReset_I(nVarReset) = (/Rho_,Bx_,By_,Bz_,T_,Ux_,Uy_,Uz_/)
    character(len=*), parameter:: NameSub='SP_put_line'
    !----------------------------------------------------------------
    if(IsFirstCall)then
       ! initialize information of which 
       IsFirstCall = .false.
       allocate(DoneExtractSolarCorona_B(nBlock))
       DoneExtractSolarCorona_B = .false.
    end if

    ! store passed particles
    do iPut = 1, nPut
       iBlock = iIndex_II(4, iPut)
       iLine  = iNode_B(iBlock)
       iParticle = iIndex_II(1, iPut)
       if(iParticle < iParticleMin)&
            call CON_stop(NameSub//': particle index is below limit')
       if(iParticle > iParticleMax)&
            call CON_stop(NameSub//': particle index is above limit')
       iGridLocal_IB(Begin_,iBlock)=MIN(iGridLocal_IB(Begin_,iBlock),iParticle)
       iGridLocal_IB(End_,  iBlock)=MAX(iGridLocal_IB(End_,  iBlock),iParticle)
       if(iGridGlobal_IA(Proc_, iLine) /= iProc)&
            call CON_stop(NameSub//': Incorrect message pass')
       
       ! check if the particle has crossed the solar corona boundary
       WasUndef = State_VIB(R_, iParticle, iBlock) < 0.0
       WasInSc  = State_VIB(R_, iParticle, iBlock) < RSc
       IsInSc   = Coord_DI( R_, iPut)              < RSc
       if(.not.WasUndef .and. WasInSc .and. .not. IsInSc)then
          ! particle existed and crossed SC boundary
          call add_to_queue(QueueRequest, (/iBlock, iParticle/))
       elseif(.not.DoneExtractSolarCorona_B(iBlock) .and. &
            WasUndef .and. .not. IsInSc)then
          ! particle didn't exist, is the 1st beyond SC
          call add_to_queue(QueueRequest, (/iBlock, iParticle/))
          DoneExtractSolarCorona_B(iBlock) = .true.
       end if

       ! keep some variables as "old" state
       State_VIB((/RhoOld_,BOld_/),  iParticle, iBlock) = &
            State_VIB((/Rho_,B_/),   iParticle, iBlock)
       ! reset others 
       State_VIB(VarReset_I,         iParticle, iBlock) = 0.0
       ! put coordinates
       State_VIB((/R_, Lon_, Lat_/), iParticle, iBlock) = &
            Coord_DI(1:nDim, iPut)
    end do
  end subroutine SP_put_line

  !===================================================================

  subroutine SP_get_grid_descriptor_param(&
       iGridMin_D, iGridMax_D, Displacement_D)
    integer, intent(out):: iGridMin_D(nDim)
    integer, intent(out):: iGridMax_D(nDim)
    real,    intent(out):: Displacement_D(nDim)
    !-----------------------------------------
    iGridMin_D = (/iParticleMin, 1, 1/)
    iGridMax_D = (/iParticleMax, 1, 1/)
    Displacement_D = 0.0
  end subroutine SP_get_grid_descriptor_param

  !===================================================================

  subroutine SP_get_line_all(Xyz_DI)
    use ModMpi
    real, pointer:: Xyz_DI(:, :)

    integer:: iNode, iParticle, iBlock, iError
    ! radius-lon-lat coordinates
    real:: Coord_D(nDim)
    !-----------------------------------------
    Xyz_DI = 0.0
    do iNode = 1, nNode
       if(iGridGlobal_IA(Proc_, iNode) /= iProc)&
            CYCLE
       iBlock = iGridGlobal_IA(Block_, iNode)
       do iParticle = iParticleMin, iParticleMax
          if(  iParticle < iGridLocal_IB(Begin_, iBlock) .or. &
               iParticle > iGridLocal_IB(End_,   iBlock)) &
               CYCLE
          Coord_D = State_VIB((/R_,Lon_,Lat_/), iParticle, iBlock)
          call rlonlat_to_xyz(Coord_D, &
               Xyz_DI(:, (iNode-1)*nParticle+iParticle-iParticleMin+1) )
       end do
    end do
    call MPI_Allreduce(MPI_IN_PLACE, Xyz_DI, nParticle*nNode*nDim, MPI_REAL, &
         MPI_SUM, iComm, iError)

  end subroutine SP_get_line_all

end module SP_wrapper
