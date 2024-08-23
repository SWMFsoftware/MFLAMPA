!  Copyright (C) 2002 Regents of the University of Michigan,
!  portions used with permission
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module SP_ModAdvance

  ! The module contains methods for advancing the solution in time
  use SP_ModDiffusion,    ONLY: UseDiffusion, set_diffusion_coef
  use SP_ModDistribution, ONLY: IsDistNeg
  use SP_ModGrid,         ONLY: State_VIB, MHData_VIB, &
       B_, U_, D_, iShock_IB, Shock_, ShockOld_, Rho_, &
       Used_B, nLine, nVertex_B, IsMuAvg, dLogRhoThreshold
  use SP_ModSize,         ONLY: nVertexMax
  use ModUtilities,       ONLY: CON_stop

  implicit none

  SAVE

  PRIVATE ! except

  ! Public members:
  public:: read_param           ! Read parameters
  public:: get_shock_location   ! finds shock location on all lines
  public:: advance              ! Advance time-accurate Distribution_CB
  public:: iterate_steady_state ! Iterate to get steady-state Distribution_CB

  ! If the shock wave is traced, the advance algorithms are modified
  logical, public :: DoTraceShock = .true.

  ! Local variables
  real    :: Cfl = 0.9   ! Controls the maximum allowed time step
  ! Logical variable for the advection scheme
  logical :: UsePoissonBracket = .false.
contains
  !============================================================================
  subroutine read_param(NameCommand)

    use ModReadParam, ONLY: read_var
    character(len=*), intent(in):: NameCommand ! From PARAM.in
    character(len=*), parameter:: NameSub = 'read_param'
    !--------------------------------------------------------------------------
    select case(NameCommand)
    case('#CFL')
       call read_var('Cfl', Cfl)
    case('#ADVECTION')
       call read_var('UsePoissonBracket', UsePoissonBracket)
    case('#TRACESHOCK')
       call read_var('DoTraceShock', DoTraceShock)
    case default
       call CON_stop(NameSub//': Unknown command '//NameCommand)
    end select

  end subroutine read_param
  !============================================================================
  subroutine get_shock_location

    ! find location of a shock wave on a given line (line)
    ! shock front is assumed to be location of max log(Rho/RhoOld)
    use SP_ModGrid, ONLY: NoShock_, RhoOld_, BOld_, R_, nWidth
    use SP_ModTime, ONLY: IsSteadyState
    use SP_ModUnit, ONLY: Io2Si_V, UnitX_

    ! Physical variables for determining the shock locations
    real :: dLogRho_I(nVertexMax), divU_I(nVertexMax)
    ! Array of u=\vec{u}*\vec{B}/|B| and u/B at the face center
    real :: uSi_F(nVertexMax), uOverBSi_F(0:nVertexMax)
    ! Array of B and 1/B at the cell- and face-center
    real :: BSi_C(nVertexMax), InvBSi_C(nVertexMax), InvBSi_F(nVertexMax)
    ! Distance between adjacent meshes and faces
    real :: DsMeshSi_I(nVertexMax), DsFaceSi_I(nVertexMax)
    ! d(u/B)/ds variable at the cell center
    real :: duOverBSi_C(nVertexMax)
    ! Do not search too close to the Sun
    real, parameter :: RShockMin = 1.20  ! *RSun
    integer         :: iShockMin
    ! Do not search too close to the heliosphere boundary
    integer:: iShockMax
    ! Misc
    integer:: iShockCandidate
    ! Loop variables
    integer:: iLine, iEnd
    !--------------------------------------------------------------------------

    if(IsSteadyState) then
       do iLine = 1, nLine
          if(.not.Used_B(iLine))then
             iShock_IB(Shock_,iLine) = NoShock_
             CYCLE
          end if

          ! shock front is assumed to be location of min B*d(u/B)/ds;
          ! divergence of plasma velocity (positively) propto dLogRho
          iEnd = nVertex_B(iLine)

          ! Calculate u/B = \vec{u}*\vec{B} / |\vec{B}|**2 at cell face
          ! uSi_F with the index of "i" is the value at the face between
          ! the mesh "i" and "i+1"
          uSi_F(1:iEnd-1) = State_VIB(U_, 1:iEnd-1, iLine)
          ! Get B at cell center
          BSi_C(1:iEnd)   = State_VIB(B_, 1:iEnd, iLine)
          ! Calculate 1/B at cell- and face-center
          InvBSi_C           = 1.0/BSi_C(1:iEnd)
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
          duOverBSi_C(1:iEnd) = (uOverBSi_F(1:iEnd) - &
               uOverBSi_F(0:iEnd-1))/DsFaceSi_I(1:iEnd)
          ! divergence of plasma \vec{u} = B * d(u/B)/ds at cell center
          divU_I(1:iEnd) = BSi_C(1:iEnd)*duOverBSi_C(1:iEnd)

          ! shock never moves back
          iShockMin = max(iShock_IB(ShockOld_, iLine), 1 + nWidth)
          iShockMax = iEnd - nWidth - 1
          iShockCandidate = iShockMin - 1 + minloc(&
               divU_I(   iShockMin:iShockMax), DIM=1, MASK = &
               State_VIB(R_,iShockMin:iShockMax,iLine) > RShockMin .and. &
               divU_I(   iShockMin:iShockMax)       > dLogRhoThreshold)
          if(iShockCandidate >= iShockMin)&
               iShock_IB(Shock_, iLine) = iShockCandidate
       end do
    else
       do iLine = 1, nLine
          if(.not.Used_B(iLine))then
             iShock_IB(Shock_,iLine) = NoShock_
             CYCLE
          end if

          ! shock front is assumed to be location of max log(Rho/RhoOld);
          ! divergence of plasma velocity (negatively) propto dLogRho
          iEnd = nVertex_B(iLine)
          dLogRho_I(1:iEnd) = log(MhData_VIB(Rho_,1:iEnd,iLine)/&
               State_VIB(RhoOld_,1:iEnd,iLine))

          ! shock never moves back
          iShockMin = max(iShock_IB(ShockOld_, iLine), 1 + nWidth)
          iShockMax = iEnd - nWidth - 1
          iShockCandidate = iShockMin - 1 + maxloc(&
               DLogRho_I(   iShockMin:iShockMax), DIM=1, MASK = &
               State_VIB(R_,iShockMin:iShockMax,iLine) > RShockMin .and. &
               DLogRho_I(   iShockMin:iShockMax)       > dLogRhoThreshold)
          if(iShockCandidate >= iShockMin)&
               iShock_IB(Shock_, iLine) = iShockCandidate
       end do
    end if

  end subroutine get_shock_location
  !============================================================================
  subroutine advance(TimeLimit)
    ! advance the solution of the diffusive kinetic equation:
    !    if IsMuAvg: Omni-directional VDF (Parker transport equation):
    !      f_t + [(1/3)*(d(ln rho)/dt]*f_{ln p} = B*d/ds[Dxx/B*df/ds]
    !    if not IsMuAvg: VDF with pitch angle (Focused transport equation):
    !      f_t + {f; p**3/3*(\vec{u}*\vec{B}/|B|)}_{x, p**3/3}
    !          + {f; (mu**2-1)*v/(2|B|)}_{x, mu}
    !          + {f; (1-mu**2)/2*(mu*(p**3/3)*
    !             (3\vec{b}\vec{b}:\nabla\vec{u} - \nabla\cdot\vec{u})
    !             + p**2*gamma*m_i*bDu/Dt)} = I**(scattering)
    ! with accounting for scattering and first-order Fermi acceleration
    ! from SPTime to TimeLimit
    ! Prototype: FLAMPA/src/SP_main, case("RUN"), Roussev&Sokolov2008
    ! Version: Borovikov&Sokolov, Dec.19 2017, distinctions:
    ! (1) no turbulence (2) new shock finder moved to SP_ModMain,
    ! and (3) new steepen_shock

    use SP_ModAdvanceAdvection, ONLY: advect_via_log
    use SP_ModAdvancePoisson,   ONLY: advect_via_poisson_parker, &
         calc_states_poisson_focused, advect_via_poisson_focused
    use SP_ModGrid,             ONLY: RhoOld_, BOld_, nWidth
    use SP_ModTime,             ONLY: SPTime
    use SP_ModUnit, ONLY: Io2Si_V, UnitX_

    real, intent(in):: TimeLimit
    ! Loop variable
    integer :: iLine
    ! For a given line: nVertex_B, iShock_IB:
    integer :: iEnd, iShock, iShockOld
    ! Upper limit and variable for the loop which makes a time step so short
    ! that the shock wave passes a single grid interval per a progress step:
    integer :: iProgress, nProgress
    real    :: InvnProgress
    ! Coefficient to interpolate "old" and "new"
    real    :: Alpha
    ! Full difference between DataInputTime and SPTime
    real    :: DtFull
    ! Time step in the PROGRESS Loop, DtFull/nProgress
    real    :: DtProgress
    ! Local arrays to store the state vectors in SI units
    real, dimension(1:nVertexMax) :: nOldSi_I, nSi_I, BOldSi_I, BSi_I, Mass_C
    real    :: DsSi_I(1:nVertexMax-1)
    ! Lagrangian derivatives
    real, dimension(1:nVertexMax) :: dLogRho_I

    !--------------------------------------------------------------------------
    DtFull = TimeLimit - SPTime

    ! Go line by line and advance the solution
    LINE:do iLine = 1, nLine
       if(.not.Used_B(iLine)) CYCLE LINE
       ! Number of the active particles on the line
       iEnd = nVertex_B(iLine)
       ! In M-FLAMPA DsSi_I(i) is the distance between meshes i and i+1
       DsSi_I(1:iEnd-1) = State_VIB(D_, 1:iEnd-1, iLine)*Io2Si_V(UnitX_)
       Mass_C(2:iEnd-1) = 0.5*(DsSi_I(2:iEnd-1) + DsSi_I(1:iEnd-2))*&
            MhData_VIB(Rho_, 2:iEnd-1, iLine)/State_VIB(B_, 2:iEnd-1, iLine)
       Mass_C(1)        = DsSi_I(1)*MhData_VIB(Rho_,1,iLine)&
            /State_VIB(B_, 1, iLine)
       Mass_C(iEnd)     = DsSi_I(iEnd-1)*MhData_VIB(Rho_,iEnd,iLine)&
            /State_VIB(B_, iEnd, iLine)
       ! For focused transport equation:
       ! initialize states for the multiple Poisson bracket scheme
       if(UsePoissonBracket .and. .not.IsMuAvg) &
            call calc_states_poisson_focused(iLine, iEnd, DtFull)

       ! Various data along the line in SI units.
       ! The IO units of the state vectors could be seen in ModUnit, the
       ! summary is: Length is in the unit of Rsun, temperature is in the
       ! unit of kinetic energy, all others are in SI units.
       nOldSi_I(1:iEnd) = State_VIB(RhoOld_, 1:iEnd, iLine)
       BOldSi_I(1:iEnd) = State_VIB(BOld_, 1:iEnd, iLine)

       ! Find how far shock has travelled on this line: nProgress
       iShock    = iShock_IB(Shock_,    iLine)
       iShockOld = iShock_IB(ShockOld_, iLine)
       if(DoTraceShock) then
          ! This is how many steps should be done to allow the shock to
          ! the move not more than one mesh size
          nProgress = max(1, iShock-iShockOld)
          iShockOld = min(iShockOld, iShock-1)
       else
          nProgress = 1
          iShockOld = 0
       end if

       ! Each particle shock has crossed should be
       ! processed separately => reduce the time step
       InvnProgress = 1.0/nProgress
       DtProgress = DtFull*InvnProgress

       ! Go over each crossed particle
       PROGRESS:do iProgress = 1, nProgress
          ! Account for change in the background up to the current moment
          Alpha = iProgress*InvnProgress

          ! nSi is needed to set up the distribution at the injection.
          ! It is calculated at the end of the iProgress' time step
          nSi_I(1:iEnd) = State_VIB(RhoOld_, 1:iEnd, iLine) + &
               Alpha*(MhData_VIB(Rho_, 1:iEnd, iLine) - &
               State_VIB(RhoOld_, 1:iEnd, iLine))
          BSi_I(1:iEnd) = State_VIB(BOld_, 1:iEnd, iLine) + &
               Alpha*(State_VIB(B_, 1:iEnd, iLine) - &
               State_VIB(BOld_, 1:iEnd, iLine))
          dLogRho_I(1:iEnd) = log(nSi_I(1:iEnd)/nOldSi_I(1:iEnd))

          ! trace shock position and steepen the shock
          iShock = iShockOld + iProgress
          if(iShock < iEnd-nWidth .and. iShock > nWidth &
               .and. DoTraceShock) call steepen_shock

          ! Advection (2 different schemes) and Diffusion
          ! First, set the diffusion coefficient, from the
          ! given formulae or from the turbulent specrtum, if known
          if(UseDiffusion) call set_diffusion_coef(iLine, iEnd, &
               iShock, BSi_I(1:iEnd))
          if(UsePoissonBracket) then
             ! Poisson bracket scheme: particle-number-conservative
             if(IsMuAvg) then
                ! Single Poisson bracket: Parker transport equation
                call advect_via_poisson_parker(iLine, iEnd, &
                     iShock, DtProgress, Cfl, &
                     nSi_I(1:iEnd)*exp(-dLogRho_I(1:iEnd)), &
                     nSi_I(1:iEnd), BSi_I(1:iEnd), Mass_C(1:iEnd))
             else
                ! Multiple Poisson brackets: Focused transport equation
                ! See descriptions for development in ModAdvancePoisson.f90
                call advect_via_poisson_focused(iLine, iEnd, &
                     iShock, DtProgress, Cfl, &
                     nSi_I(1:iEnd)*exp(-dLogRho_I(1:iEnd)), nSi_I(1:iEnd), &
                     BOldSi_I(1:iEnd), BSi_I(1:iEnd), Mass_C(1:iEnd))
             end if
          else
             ! No Poisson bracket scheme, use the default algorithm
             call advect_via_log(iLine, iEnd, iShock, DtProgress, Cfl, &
                  dLogRho_I(1:iEnd), nSi_I(1:iEnd), BSi_I(1:iEnd))
          end if

          ! For any scheme, check if any VDF along this line is negative
          if(IsDistNeg) CYCLE LINE
          ! Store the old density at the end of each iProgress
          nOldSi_I(1:iEnd) = nSi_I(1:iEnd)
          BOldSi_I(1:iEnd) = BSi_I(1:iEnd)
       end do PROGRESS
    end do LINE

  contains
    !==========================================================================
    subroutine steepen_shock
      ! change the density profile near the shock front so it
      ! becomes steeper for the current line

      real :: dLogRhoExcess_I(iShock-nWidth:iShock+nWidth-1)
      real :: dLogRhoExcessIntegral
      ! find the excess of dLogRho within the shock compared
      ! to background averaged over length
      !------------------------------------------------------------------------
      dLogRhoExcess_I = max(0.5*(dLogRho_I(iShock-nWidth:iShock+nWidth-1) + &
           dLogRho_I(iShock-nWidth+1:iShock+nWidth)) - dLogRhoThreshold, 0.0)
      ! A jump (dLogRhoExcess>0) in velocity accross the shock wave * \Delta t
      dLogRhoExcessIntegral = sum(dLogRhoExcess_I*&
           DsSi_I(iShock-nWidth:iShock+nWidth-1))

      ! check for zero excess
      if(dLogRhoExcessIntegral == 0.0) RETURN
      ! nullify excess within the smoothed shock
      dLogRho_I(iShock-nWidth:iShock+nWidth) = min(dLogRhoThreshold, &
           dLogRho_I(iShock-nWidth:iShock+nWidth))
      ! ... and concentrate it at the shock front, applying the whole jump
      ! in the velocity at a single grid point
      dLogRho_I(iShock) = dLogRhoThreshold + &
           dLogRhoExcessIntegral/DsSi_I(iShock)
      ! also, sharpen the magnetic field magnitude
      ! post shock part
      BSi_I(iShock+1-nWidth:iShock+1) = maxval(BSi_I(iShock+1-nWidth:iShock+1))
      ! pre shock part
      BSi_I(iShock+1:iShock+nWidth  ) = minval(BSi_I(iShock+1:iShock+nWidth))

    end subroutine steepen_shock
    !==========================================================================
  end subroutine advance
  !============================================================================
  subroutine iterate_steady_state
    ! advance sol. of the diffusive kinetic equation (from Parker Eqn.):
    !     f_t+[(1/3)*(d(ln rho)/dt]*f_{ln p}=B*d/ds[D/B*df/ds]
    ! with accounting for diffusion and first-order Fermi acceleration

    use SP_ModAdvancePoisson, ONLY: iterate_poisson_parker, &
         iterate_poisson_focused
    ! Loop variable
    integer :: iLine
    ! For a given line: nVertex_B, iShock_IB:
    integer :: iEnd, iShock
    ! Local arrays to store the state vectors in SI units
    real, dimension(1:nVertexMax):: nSi_I, BSi_I

    !--------------------------------------------------------------------------
    LINE:do iLine = 1, nLine
       ! go line by line and iterate the solution
       if(.not.Used_B(iLine)) CYCLE LINE
       ! the active particles on the line
       iEnd = nVertex_B(iLine)

       ! get data along the line in SI units: Temperature is in the unit
       ! of kinetic energy, all others are in SI units.
       BSi_I(1:iEnd) = State_VIB(   B_, 1:iEnd, iLine)
       ! nSi is needed to set up the distribution at the injection.
       nSi_I(1:iEnd) = MhData_VIB(Rho_, 1:iEnd, iLine)
       ! find how far shock has travelled on this line
       iShock = iShock_IB(Shock_, iLine)

       ! first, set the diffusion coefficient, from the
       ! given formulae or from the turbulent specrtum, if known
       if(UseDiffusion) call set_diffusion_coef(iLine, iEnd, &
            iShock, BSi_I(1:iEnd))
       ! Poisson bracket scheme: particle-number-conservative
       if(IsMuAvg) then
          ! Single Poisson bracket: Parker transport equation
          call iterate_poisson_parker(iLine, iEnd, iShock, Cfl, &
               BSi_I(1:iEnd), nSi_I(1:iEnd))
       else
          call iterate_poisson_focused(iLine, iEnd, iShock, Cfl, &
               BSi_I(1:iEnd), nSi_I(1:iEnd))
       end if

       ! for any scheme, check if any VDF along this line is negative
       if(IsDistNeg) CYCLE LINE
    end do LINE

  end subroutine iterate_steady_state
  !============================================================================
end module SP_ModAdvance
!==============================================================================
