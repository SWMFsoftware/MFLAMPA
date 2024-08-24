!  Copyright (C) 2002 Regents of the University of Michigan,
!  portions used with permission
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module SP_ModShock

  use SP_ModGrid,   ONLY: nLine, Used_B, nVertex_B, &
       State_VIB, MhData_VIB, NoShock_, Shock_, ShockOld_
  use ModUtilities, ONLY: CON_stop

  implicit none

  SAVE

  PRIVATE ! except

  ! Public members:
  public:: read_param           ! Read parameters
  public:: get_divU             ! get divergence of velocity \vec{u}
  public:: get_shock_location   ! finds shock location on all lines

  ! If the shock wave is traced, the advance algorithms are modified
  logical, public :: DoTraceShock = .true.

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
  subroutine get_divU(iLine, nX, divU_I)

    ! divU = dLogRho_I for time-accurate run
    ! divU = B*d(u/B)/ds in steady-state run
    use SP_ModGrid, ONLY: U_, B_, D_, Rho_, RhoOld_
    use SP_ModTime, ONLY: IsSteadyState
    use SP_ModUnit, ONLY: Io2Si_V, UnitX_

    ! Input: line index
    integer, intent(in) :: iLine
    ! Input: nX = iEnd = nVertex_B(iLine)
    integer, intent(in) :: nX
    ! Output: divU for determining the shock locations
    ! divU = -dLogRho_I for time-accurate run, and = B*d(u/B)/ds in steady-state
    real,   intent(out) :: divU_I(nX)
    ! Local VARs
    ! Physical variables 
    ! Array of u=\vec{u}*\vec{B}/|B| and u/B at the face center
    real :: uSi_F(nX), uOverBSi_F(0:nX)
    ! Array of B and 1/B at the cell- and face-center
    real, dimension(nX) :: BSi_C, InvBSi_C, InvBSi_F
    ! Distance between adjacent meshes and faces
    real, dimension(nX) :: DsMeshSi_I, DsFaceSi_I
    ! d(u/B)/ds variable at the cell center
    real, dimension(nX) :: DuOverBDsSi_C
    !--------------------------------------------------------------------------

    if(IsSteadyState) then
       ! Calculate u/B = \vec{u}*\vec{B} / |\vec{B}|**2 at cell face
       ! uSi_F with the index of "i" is the value at the face between
       ! the mesh "i" and "i+1"
       uSi_F(1:nX-1) = State_VIB(U_, 1:nX-1, iLine)
       ! Get B at cell center
       BSi_C(1:nX)   = State_VIB(B_, 1:nX, iLine)
       ! Calculate 1/B at cell- and face-center
       InvBSi_C(1:nX)   = 1.0/BSi_C(1:nX)
       InvBSi_F(1:nX-1) = 0.5*(InvBSi_C(1:nX-1) + InvBSi_C(2:nX))

       ! u/B with the index of "i" is the value at the face between
       ! the mesh "i" and "i+1"
       ! Average 1/B and multiply by uSi at face centers
       uOverBSi_F(1:nX-1) = uSi_F(1:nX-1)*InvBSi_F(1:nX-1)
       uOverBSi_F(0 )     = uSi_F(1)*InvBSi_C(1)
       uOverBSi_F(nX)     = uSi_F(nX-1)*InvBSi_C(nX)

       ! In M-FLAMPA DsMeshSi_I(i) is the distance between meshes i and i+1
       DsMeshSi_I(1:nX-1) = State_VIB(D_, 1:nX-1, iLine)*Io2Si_V(UnitX_)
       ! Within the framework of finite volume method, the cell volume
       ! is used, which is proportional to the distance between the faces
       ! bounding the volume with an index, i, which is half of sum of
       ! distance between meshes i-1 and i, i.e. DsMeshSi_I(i-1), and that
       ! between meshes i and i+1, which is DsMeshSi_I(i):
       DsFaceSi_I(2:nX-1) = 0.5*(DsMeshSi_I(1:nX-2)+DsMeshSi_I(2:nX-1))
       DsFaceSi_I(1)      = DsMeshSi_I(1)
       DsFaceSi_I(nX)     = DsMeshSi_I(nX-1)

       ! In each cell, d(u/B)/ds = [(u/B at the right face) - (u/B at the
       ! left face)]/(distance between the left and right faces)
       DuOverBDsSi_C(1:nX)  = (uOverBSi_F(1:nX) - &
            uOverBSi_F(0:nX-1))/DsFaceSi_I(1:nX)
       ! divergence of plasma \vec{u} = B * d(u/B)/ds at cell center
       divU_I(1:nX) = BSi_C(1:nX)*DuOverBDsSi_C(1:nX)
    else
       ! divergence of plasma \vec{u} = -d(ln(rho))/dt at cell center
       divU_I(1:nX) = -log(MhData_VIB(Rho_,1:nX,iLine)/ &
            State_VIB(RhoOld_,1:nX,iLine))
    end if

  end subroutine get_divU
  !============================================================================
  subroutine get_shock_location

    ! find location of a shock wave on a given line (line)
    ! shock front is assumed to be location of max log(Rho/RhoOld)
    use SP_ModGrid, ONLY: iShock_IB, R_
    use SP_ModSize, ONLY: nVertexMax

    ! divU_I (at shock) is importent for particle acceleration processes
    real :: divU_I(nVertexMax)
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
       if(.not.Used_B(iLine)) then
          iShock_IB(Shock_,iLine) = NoShock_
          CYCLE
       end if
       ! Number of the active particles on the line
       iEnd = nVertex_B(iLine)

       ! if it is time accurate:
       ! shock front is assumed to be location of max log(Rho/RhoOld);
       ! divergence of plasma velocity (negatively) propto dLogRho
       ! if it is steady state:
       ! shock front is assumed to be location of min B*d(u/B)/ds;
       ! divergence of plasma velocity (positively) propto dLogRho
       call get_divU(iLine, iEnd, divU_I(1:iEnd))
       ! shock never moves back
       iShockMin = max(iShock_IB(ShockOld_, iLine), nShockWidth+1)
       iShockMax = iEnd - nShockWidth - 1
       iShockCandidate = iShockMin - 1 + minloc( &
            divU_I(iShockMin:iShockMax), DIM=1, MASK= &
            State_VIB(R_,iShockMin:iShockMax,iLine) > RShockMin .and. &
            abs(divU_I(iShockMin:iShockMax)) > dLogRhoThreshold)
       if(iShockCandidate >= iShockMin) &
            iShock_IB(Shock_, iLine) = iShockCandidate
    end do

  end subroutine get_shock_location
  !============================================================================

  !============================================================================
end module SP_ModShock
