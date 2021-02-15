!  Copyright (C) 2002 Regents of the University of Michigan,
!  portions used with permission
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module SP_wrapper

  use SP_ModMain, ONLY: &
       run, save_restart, &
       DoRestart, DoReadMhData, &
       nDim, nBlock, nParticleMax, &
       RMin=>RScMin, RBufferMin=>RIhMin, &
       RBufferMax=>RScMax, RMax=>RIhMax, &
       MHData_VIB, iNode_B, FootPoint_VB, DataInputTime, &
       nParticle_B, Length_, LagrID_,X_, Y_, Z_
  use CON_comp_info
  use CON_coupler, ONLY: &
       set_coord_system, SP_, is_proc0, i_comm, i_proc0, &
       init_decomposition, get_root_decomposition, bcast_decomposition, &
       iVar_V, DoCoupleVar_V, &
       Density_, RhoCouple_, Pressure_, PCouple_, &
       Momentum_, RhoUxCouple_, RhoUzCouple_, &
       BField_, BxCouple_, BzCouple_, &
       Wave_, WaveFirstCouple_, WaveLastCouple_
  use ModConst, ONLY: rSun
  use ModMpi
  use CON_world,  ONLY: is_proc0, is_proc, n_proc
  implicit none
  save
  private ! except
  public:: SP_set_param
  public:: SP_init_session
  public:: SP_run
  public:: SP_save_restart
  public:: SP_finalize

  ! coupling with MHD components
  public:: SP_adjust_lines
  public:: SP_get_bounds_comp
  public:: SP_check_ready_for_mh
  public:: SP_put_coupling_param
  ! variables requested via coupling: coordinates,
  ! field line and particles indexes
  character(len=*), parameter:: NameVarCouple =&
       'rho p mx my mz bx by bz i01 i02 pe'
  integer, parameter:: Lower_=0, Upper_=1
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
    use SP_ModMain, ONLY: initialize, DoRestart, DoReadMhData
    use SP_ModOriginPoints, ONLY: get_origin_points
    use SP_ModGrid, ONLY: init_grid=>init
    use CON_mflampa, ONLY: set_state_pointer
    use SP_ModUnit, ONLY:  NameEnergyUnit
    use ModConst,   ONLY: energy_in
    integer,  intent(in) :: iSession         ! session number (starting from 1)
    real,     intent(in) :: TimeSimulation   ! seconds from start time

    logical, save:: IsInitialized = .false.
    !--------------------------------------------------------------------------
    if(IsInitialized)&
         RETURN
    IsInitialized = .true.
    call init_grid
    nullify(MHData_VIB); nullify(nParticle_B)
    call set_state_pointer(MHData_VIB, nParticle_B, &
         nBlock, nParticleMax,                      &
         rMin, rMax,                                &
         1/energy_in(NameEnergyUnit))
    call initialize
    if(.not.(DoRestart.or.DoReadMhData))&
         call get_origin_points
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
  subroutine SP_set_grid
    use SP_ModGrid, ONLY: iGridGlobal_IA, Block_, Proc_, &
          TypeCoordSystem, nLon, nLat
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
    use CON_mflampa, ONLY: get_bounds
    integer,        intent(in) :: iModelIn
    real,           intent(in) :: rMinIn, rMaxIn
    real,           intent(in) :: TimeIn
    real, optional, intent(in) :: rBufferLoIn, rBufferUpIn
    !--------------------------------------------------------------------------
    call get_bounds(iModelIn, rMinIn, rMaxIn, &
       rBufferLoIn, rBufferUpIn)
    if(DataInputTime >= TimeIn)RETURN
    ! New coupling time, get it and save old state
    DataInputTime = TimeIn
    if(iModelIn==Lower_)then
       call copy_old_state
    else
       call CON_stop("Time in IHSP coupling differs from that in SCSP")
    end if
  end subroutine SP_put_coupling_param
  !============================================================================
  ! Called from coupler after the updated grid point lo<cation are
  ! received from the other component (SC, IH). Determines whether some
  ! grid points should be added/deleted
  subroutine SP_adjust_lines(DoInit)
    use SP_ModDistribution, ONLY: offset
    use SP_ModGrid,         ONLY: R_, State_VIB
    use CON_mflampa, ONLY: iOffset_B, rBufferLo, rBufferUp, &
         rInterfaceMin, rInterfaceMax
    logical, intent(in) :: DoInit
    ! once new geometry of lines has been put, account for some particles
    ! exiting the domain (can happen both at the beginning and the end)
    integer:: iParticle, iBlock, iBegin,  iEnd, iOffset ! loop variables
    integer:: iParticle_I(2), iLoop
    logical:: DoAdjustLo, DoAdjustUp
    logical:: IsMissing

    real              :: R
    
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
          R = State_VIB(R_,iParticle,iBlock)
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
               R <  RInterfaceMin&
               .or.&
               R >= RInterfaceMax)&
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
          if(R < RInterfaceMin)&
               call CON_stop(NameSub//": particle has been lost")

          ! missing point in the upper part of the domain -> IGNORE;
          ! if needed to adjust beginning, then it is done,
          ! switch left -> right end of range and start adjusting
          ! tail of the line, if it has reentered current part of the domain
          if(R >= RInterfaceMax)then
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
          if(R >= rBufferUp .and. R < rInterfaceMax)then
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
end module SP_wrapper
