!  Copyright (C) 2002 Regents of the University of Michigan,
!  portions used with permission
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module SP_ModShock

  ! This module contains subroutines for determining the shock location,
  ! and steepening the density and magnetic field strength at shock front.
  use SP_ModGrid,   ONLY: nLine, Used_B, nVertex_B, State_VIB, &
       MhData_VIB, iShock_IB, NoShock_, Shock_, ShockOld_
  use SP_ModSize,   ONLY: nVertexMax
  use ModUtilities, ONLY: CON_stop

  implicit none

  SAVE

  PRIVATE ! except

  ! Public members:
  public:: read_param           ! Read parameters
  public:: init                ! Initialize arrays on the grid
  public:: get_divU             ! get divergence of velocity \vec{u}
  public:: get_shock_location   ! finds shock location on all lines
  public:: steepen_shock        ! steepen the density profile at the shock

  ! If the shock wave is traced, the advance algorithms are modified
  logical, public :: DoTraceShock = .true.
  ! divergence of velocity \vec{U}: for determining the shock locations
  real, public, allocatable :: divU_II(:,:)

  ! Shock algorithm parameters:
  real,    public, parameter :: dLogRhoThreshold = 0.01
  integer, public, parameter :: nShockWidth = 50
contains
  !============================================================================
  subroutine read_param(NameCommand)

    use ModReadParam, ONLY: read_var
    character(len=*), intent(in):: NameCommand ! From PARAM.in
    character(len=*), parameter:: NameSub = 'read_param'
    !--------------------------------------------------------------------------
    select case(NameCommand)
    case('#TRACESHOCK')
       call read_var('DoTraceShock', DoTraceShock)
    case default
       call CON_stop(NameSub//': Unknown command '//NameCommand)
    end select

  end subroutine read_param
  !============================================================================
  subroutine init

    ! initialize the divU array along s_L axis for all field lines

    use ModUtilities, ONLY: check_allocate
    use SP_ModProc,   ONLY: iError
    character(len=*), parameter:: NameSub = 'init'
    !--------------------------------------------------------------------------

    if(allocated(divU_II)) deallocate(divU_II)
    allocate(divU_II(1:nVertexMax, 1:nLine))
    call check_allocate(iError, 'divU_II')

  end subroutine init
  !============================================================================
  subroutine get_divU

    ! divU = dLogRho_I for time-accurate run
    ! divU = B*d(u/B)/ds in steady-state run
    use SP_ModGrid, ONLY: U_, B_, D_, Rho_, RhoOld_
    use SP_ModTime, ONLY: IsSteadyState
    use SP_ModUnit, ONLY: Io2Si_V, UnitX_

    ! Local VARs
    ! Array of u=\vec{u}*\vec{B}/|B| and u/B at the face center
    real :: uSi_F(nVertexMax), uOverBSi_F(0:nVertexMax)
    ! Array of B and 1/B at the cell- and face-center
    real, dimension(nVertexMax) :: BSi_I, InvBSi_C, InvBSi_F
    ! Distance between adjacent meshes and faces
    real, dimension(nVertexMax) :: DsMeshSi_I, DsFaceSi_I
    ! d(u/B)/ds variable at the cell center
    real, dimension(nVertexMax) :: DuOverBDsSi_C
    ! Loop variables
    integer :: iLine, iEnd
    !--------------------------------------------------------------------------

    if(IsSteadyState) then
       do iLine = 1, nLine
          ! go line by line and get divU if active
          if(.not.Used_B(iLine)) then
             iShock_IB(Shock_,iLine) = NoShock_
             CYCLE
          end if
          ! Number of the active particles on the line
          iEnd = nVertex_B(iLine)

          ! Calculate u/B = \vec{u}*\vec{B} / |\vec{B}|**2 at cell face
          ! uSi_F with the index of "i" is the value at the face between
          ! the mesh "i" and "i+1"
          uSi_F(1:iEnd-1) = State_VIB(U_, 1:iEnd-1, iLine)
          ! Get B at cell center
          BSi_I(1:iEnd)   = State_VIB(B_, 1:iEnd, iLine)
          ! Calculate 1/B at cell- and face-center
          InvBSi_C(1:iEnd)   = 1.0/BSi_I(1:iEnd)
          InvBSi_F(1:iEnd-1) = 0.5*(InvBSi_C(1:iEnd-1) + InvBSi_C(2:iEnd))

          ! u/B with the index of "i" is the value at the face between
          ! the mesh "i" and "i+1"
          ! Average 1/B and multiply by uSi at face centers
          uOverBSi_F(1:iEnd-1) = uSi_F(1:iEnd-1)*InvBSi_F(1:iEnd-1)
          uOverBSi_F(0 )       = uSi_F(1)*InvBSi_C(1)
          uOverBSi_F(iEnd)     = uSi_F(iEnd-1)*InvBSi_C(iEnd)

          ! In M-FLAMPA DsMeshSi_I(i) is the distance between meshes i and i+1
          DsMeshSi_I(1:iEnd-1) = State_VIB(D_, 1:iEnd-1, iLine)*Io2Si_V(UnitX_)
          ! Within the framework of finite volume method, the cell volume
          ! is used, which is proportional to the distance between the faces
          ! bounding the volume with an index, i, which is half of sum of
          ! distance between meshes i-1 and i, i.e. DsMeshSi_I(i-1), and that
          ! between meshes i and i+1, which is DsMeshSi_I(i):
          DsFaceSi_I(2:iEnd-1)= 0.5*(DsMeshSi_I(1:iEnd-2)+DsMeshSi_I(2:iEnd-1))
          DsFaceSi_I(1)       = DsMeshSi_I(1)
          DsFaceSi_I(iEnd)    = DsMeshSi_I(iEnd-1)

          ! In each cell, d(u/B)/ds = [(u/B at the right face) - (u/B at the
          ! left face)]/(distance between the left and right faces)
          DuOverBDsSi_C(1:iEnd) = (uOverBSi_F(1:iEnd) - &
               uOverBSi_F(0:iEnd-1))/DsFaceSi_I(1:iEnd)
          ! divergence of plasma \vec{u} = B * d(u/B)/ds at cell center
          divU_II(1:iEnd, iLine) = BSi_I(1:iEnd)*DuOverBDsSi_C(1:iEnd)
       end do
    else
       do iLine = 1, nLine
          ! go line by line and get divU if active
          if(.not.Used_B(iLine)) then
             iShock_IB(Shock_,iLine) = NoShock_
             CYCLE
          end if
          ! Number of the active particles on the line
          iEnd = nVertex_B(iLine)

          ! divergence of plasma \vec{u} = -d(ln(rho))/dt at cell center
          divU_II(1:iEnd, iLine) = -log(MhData_VIB(Rho_,1:iEnd,iLine)/ &
               State_VIB(RhoOld_,1:iEnd,iLine))
       end do
    end if

  end subroutine get_divU
  !============================================================================
  subroutine get_shock_location

    ! find location of a shock wave on a given line (line)
    ! shock front is assumed to be location of max log(Rho/RhoOld)
    use SP_ModGrid, ONLY: R_

    ! Do not search too close to the Sun
    real, parameter :: RShockMin = 1.20  ! *RSun
    integer         :: iShockMin
    ! Do not search too close to the heliosphere boundary
    integer :: iShockMax
    ! Misc
    integer :: iShockCandidate
    ! Loop variables
    integer :: iLine, iEnd
    !--------------------------------------------------------------------------

    do iLine = 1, nLine
       ! go line by line and get the shock location if active
       if(.not.Used_B(iLine)) CYCLE
       ! Number of the active particles on the line
       iEnd = nVertex_B(iLine)

       ! if it is time accurate:
       ! shock front is assumed to be location of max log(Rho/RhoOld);
       ! divergence of plasma velocity (negatively) propto dLogRho
       ! if it is steady state:
       ! shock front is assumed to be location of min B*d(u/B)/ds;
       ! divergence of plasma velocity (positively) propto dLogRho
       ! note: shock never moves back
       iShockMin = max(iShock_IB(ShockOld_, iLine), nShockWidth+1)
       iShockMax = iEnd - nShockWidth - 1
       iShockCandidate = iShockMin - 1 + minloc( &
            divU_II(iShockMin:iShockMax, iLine), DIM=1, MASK= &
            State_VIB(R_,iShockMin:iShockMax,iLine) > RShockMin .and. &
            divU_II(iShockMin:iShockMax, iLine) < -dLogRhoThreshold)
       if(iShockCandidate >= iShockMin) &
            iShock_IB(Shock_, iLine) = iShockCandidate
    end do

  end subroutine get_shock_location
  !============================================================================
  subroutine steepen_shock(iLine, nX, iShock, BSi_I, dLogRhoIn_I)

    ! change the density profile near the shock front
    ! so it becomes steeper for the current line
    use SP_ModGrid, ONLY: D_
    use SP_ModUnit, ONLY: Io2Si_V, UnitX_

    ! INPUTs
    integer, intent(in) :: iLine, iShock ! Indices of line and shock front
    integer, intent(in) :: nX            ! Number of meshes along s_L axis
    real, intent(inout) :: BSi_I(nX)     ! Magnetic field strength
    real, optional, intent(inout) :: dLogRhoIn_I(nX) ! for time-accurate run
    ! Local VARs
    real :: DsSi_I(1:nX-1), dLogRho_I(nX)
    real :: dLogRhoExcess_I(iShock-nShockWidth:iShock+nShockWidth-1)
    real :: dLogRhoExcessSum
    !--------------------------------------------------------------------------
    DsSi_I(1:nX-1) = State_VIB(D_, 1:nX-1, iLine)*Io2Si_V(UnitX_)
    ! get dLogRho_I if given; otherwise = -divU
    if(present(dLogRhoIn_I)) then
        dLogRho_I = dLogRhoIn_I
    else
        dLogRho_I = -divU_II(1:nX, iLine)
    end if

    ! find the excess of dLogRho within the shock compared
    ! to background averaged over length
    dLogRhoExcess_I = max(0.5*( &
         dLogRho_I(iShock-nShockWidth:iShock+nShockWidth-1) + &
         dLogRho_I(iShock-nShockWidth+1:iShock+nShockWidth)) - &
         dLogRhoThreshold, 0.0)
    ! a jump (dLogRhoExcess>0) in velocity accross the shock wave * \Delta t
    dLogRhoExcessSum = sum(dLogRhoExcess_I* &
         DsSi_I(iShock-nShockWidth:iShock+nShockWidth-1))

    ! check for zero excess
    if(dLogRhoExcessSum == 0.0) RETURN
    ! nullify excess within the smoothed shock
    dLogRho_I(iShock-nShockWidth:iShock+nShockWidth) = min( &
         dLogRhoThreshold, dLogRho_I(iShock-nShockWidth:iShock+nShockWidth))
    ! ... and concentrate it at the shock front, applying the whole jump
    ! in the velocity at a single grid point
    dLogRho_I(iShock) = dLogRhoThreshold + &
         dLogRhoExcessSum/DsSi_I(iShock)
    ! also, sharpen the magnetic field magnitude
    ! post shock part
    BSi_I(iShock-nShockWidth+1:iShock+1) = &
         maxval(BSi_I(iShock+1-nShockWidth:iShock+1))
    ! pre shock part
    BSi_I(iShock+1:iShock+nShockWidth  ) = &
         minval(BSi_I(iShock+1:iShock+nShockWidth))
    
    ! update dLogRhoIn (if given) and divU
    if(present(dLogRhoIn_I)) then
        dLogRhoIn_I = dLogRho_I
        divU_II(1:nX, iLine) = -dLogRhoIn_I
    else
        divU_II(1:nX, iLine) = -dLogRho_I
    end if

  end subroutine steepen_shock
  !============================================================================
end module SP_ModShock
