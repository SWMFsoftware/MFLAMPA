!  Copyright (C) 2002 Regents of the University of Michigan,
!  portions used with permission
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module SP_ModAdvancePoisson

  ! High resolution finite volume method for kinetic equations
  ! with Poisson brackets (Sokolov et al., 2023)
  ! See https://doi.org/10.1016/j.jcp.2023.111923

  use SP_ModSize,         ONLY: nVertexMax
  use SP_ModBc,           ONLY: set_momentum_bc, set_VDF, &
       UseUpperEndBc, UseLowerEndBc, iStart
  use SP_ModDistribution, ONLY: nP, VolumeP_I, Momentum3_I, &
       nMu, IsMuAvg, Distribution_CB, IsDistNeg, check_dist_neg
  use SP_ModDiffusion,    ONLY: UseDiffusion, diffuse_distribution
  use ModUtilities,       ONLY: CON_stop
  use ModPoissonBracket,  ONLY: explicit

  implicit none

  PRIVATE ! Except

  SAVE

  ! Public members:
  public:: read_param                  ! Read parameters
  ! For solving the Parker transport equation
  public:: advect_via_poisson_parker   ! Time-accurate advance through given Dt
  public:: iterate_poisson_parker      ! Local time-stepping for steady states
  ! For solving the focused transport equation with pitch angles
  public:: init_states_poisson_focused ! Initialize the states for multi-P.B.
  public:: advect_via_poisson_focused  ! Time-accurate advance through given Dt

  ! Whether to include the betatron acceleration in focused transport equation
  logical, public :: UseBetatron = .false.
  ! Whether to include the inertial force in focused transport equation
  logical, public :: UseInertialForce = .false.
  ! Time-derivatives over DtFull (linear interpolation for time):
  real, public, dimension(nVertexMax) :: bDuDt_C
contains
  !============================================================================
  subroutine read_param(NameCommand)

    use ModReadParam, ONLY: read_var
    use SP_ModProc,   ONLY: iProc
    character(len=*), intent(in):: NameCommand ! From PARAM.in
    logical :: UseBetatronIn, UseInertialForceIn 
    character(len=*), parameter:: NameSub = 'read_param'
    !--------------------------------------------------------------------------
    select case(NameCommand)
    case('#FOCUSEDTRANSPORT')
       call read_var('UseBetatron',      UseBetatronIn)
       call read_var('UseInertialForce', UseInertialForceIn)
       if(.not.IsMuAvg) then
          UseBetatron = UseBetatronIn
          UseInertialForce = UseInertialForceIn
       else
          if(iProc==0) write(*,'(a,i6,a)') NameSub//' Code is configured '// &
               'with nMu=', nMu, ' with no use of focused transport equation.'
          call CON_stop('Modify PARAM.in or reconfigure SP/MFLAMPA')
       end if
    case default
       call CON_stop(NameSub//': Unknown command '//NameCommand)
    end select

  end subroutine read_param
  !============================================================================
  subroutine advect_via_poisson_parker(iLine, nX, iShock, &
       tFinal, CflIn, nOldSi_I, nSi_I, BSi_I)
    ! advect via Possion Bracket scheme: (p**3/3, s_L)
    ! diffuse the distribution function at each time step

    ! INPUTS:
    integer, intent(in) :: iLine, iShock ! Indices of line and shock
    integer, intent(in) :: nX            ! Number of meshes along s_L axis
    real,    intent(in) :: tFinal        ! Time interval to advance through
    real,    intent(in) :: CflIn         ! Input CFL number
    ! Input variables for diffusion, assuming at the time of tFinal
    real,    intent(in) :: nOldSi_I(nX), nSi_I(nX), BSi_I(nX)
    ! Extended arrays for implementation of the Poisson bracket scheme
    ! VolumeXStart_I: geometric volume when the subroutine starts
    ! VolumeXEnd_I: geometric volume when the subroutine ends
    ! dVolumeXDt_I: time derivative of geometric volume
    real, dimension(0:nX+1) :: VolumeXStart_I, VolumeXEnd_I, dVolumeXDt_I
    ! VolumeOld_G: total control volume of last each iteration
    ! Volume_G: total control volume at the end of each iteration
    ! dVolumeDt_G: total control time derivative
    real, dimension(0:nP+1, 0:nX+1) :: VolumeOld_G, Volume_G, dVolumeDt_G
    ! DeltaHamiltonian
    real :: dHamiltonian01_FX(-1:nP+1, 0:nX+1)
    ! Extended array for distribution function
    real :: VDF_G(-1:nP+2, -1:nX+2)
    ! Advection term
    real :: Source_C(nP, nX)
    ! Time, ranging from 0 to tFinal
    real :: Time
    ! Time step
    real :: Dt
    ! Prediction for the next time step:
    real :: DtNext

    ! Now this is the particle-number-conservative advection scheme
    character(len=*), parameter:: NameSub = 'advect_via_single_poisson'
    !--------------------------------------------------------------------------
    IsDistNeg = .false.

    ! Initialize arrays
    ! Geometric volume: for DeltaS/B, with 1 ghost point at the boundary
    ! Start volume
    VolumeXStart_I(1:nX) = 1.0/nOldSi_I(1:nX)
    VolumeXStart_I(0)    = VolumeXStart_I(1)
    VolumeXStart_I(nX+1) = VolumeXStart_I(nX)
    ! End volume
    VolumeXEnd_I(1:nX)   = 1.0/nSi_I(1:nX)
    VolumeXEnd_I(0)      = VolumeXEnd_I(1)
    VolumeXEnd_I(nX+1)   = VolumeXEnd_I(nX)
    ! Time derivative
    dVolumeXDt_I         = (VolumeXEnd_I - VolumeXStart_I)/tFinal
    ! Total control volume: initial and time derivative
    Volume_G             = spread(VolumeP_I, DIM=2, NCOPIES=nX+2)* &
         spread(VolumeXStart_I, DIM=1, NCOPIES=nP+2)
    dVolumeDt_G          = spread(VolumeP_I, DIM=2, NCOPIES=nX+2)* &
         spread(dVolumeXDt_I, DIM=1, NCOPIES=nP+2)
    ! Calculate 1st Hamiltonian function used in the time-dependent
    ! poisson bracket: {f_jk, (p**3/3)*(DeltaS/B)}_{tau, p**3/3}
    dHamiltonian01_FX    = -spread(Momentum3_I, DIM=2, NCOPIES=nX+2)* &
         spread(dVolumeXDt_I, DIM=1, NCOPIES=nP+3)
    ! Time initialization
    Time   = 0.0

    ! Update Bc for VDF at minimal energy, at nP = 0
    call set_momentum_bc(iLine, nX, nSi_I, iShock)
    ! Trial time step: Get DtNext
    call set_VDF(iLine, nX, VDF_G) ! Set the VDF first
    call explicit(nP, nX, VDF_G, Volume_G, Source_C,  &
         dHamiltonian01_FX = dHamiltonian01_FX,       &
         dVolumeDt_G = dVolumeDt_G,                   &
         CFLIn = CflIn, DtOut = DtNext)

    ! Advection by the single-Poisson-bracket scheme
    do
       ! Update Time step
       Dt = min(DtNext, tFinal - Time)
       ! Update Volumes
       VolumeOld_G = Volume_G
       Volume_G    = VolumeOld_G + Dt*dVolumeDt_G

       ! Update Bc for VDF
       call set_VDF(iLine, nX, VDF_G)
       ! Advance by single-Poisson-bracket scheme
       call explicit(nP, nX, VDF_G, Volume_G, Source_C,  &
            dHamiltonian01_FX = dHamiltonian01_FX,       &
            dVolumeDt_G = dVolumeDt_G,                   &
            DtIn = Dt,         & ! Input time step, which may be reduced
            CFLIn = CflIn,     & ! Input CFL to calculate next time step
            DtOut = DtNext)

       ! May need to correct the volume if the time step has been reduced
       Volume_G = VolumeOld_G + Dt*dVolumeDt_G
       ! Update velocity distribution function
       Distribution_CB(1:nP, 1, iStart:nX, iLine) = &
            Distribution_CB(1:nP, 1, iStart:nX, iLine) + &
            Source_C(1:nP, iStart:nX)
       ! Check if the VDF includes negative values
       call check_dist_neg(NameSub, 1, nX, iLine)
       if(IsDistNeg) RETURN

       ! Diffuse the distribution function
       if(UseDiffusion) then
          if(UseUpperEndBc) then
             if(UseLowerEndBc) then
                ! with lower or upper end BCs
                call diffuse_distribution(iLine, nX, iShock, Dt,      &
                     nSi_I, BSi_I, LowerEndSpectrum_I=VDF_G(1:nP, 0), &
                     UpperEndSpectrum_I=VDF_G(1:nP, nX+1))
             else
                ! with upper end BC but no lower end BC
                call diffuse_distribution(iLine, nX, iShock, Dt,      &
                     nSi_I, BSi_I, UpperEndSpectrum_I=VDF_G(1:nP, nX+1))
             end if
          else
             if(UseLowerEndBc) then
                ! with lower end BC but no upper end BC
                call diffuse_distribution(iLine, nX, iShock, Dt,      &
                     nSi_I, BSi_I, LowerEndSpectrum_I=VDF_G(1:nP, 0))
             else
                ! with no lower or upper end BCs
                call diffuse_distribution(iLine, nX, iShock, Dt, nSi_I, BSi_I)
             end if
          end if
          ! Check if the VDF includes negative values after diffusion
          call check_dist_neg(NameSub//' after diffusion', 1, nX, iLine)
          if(IsDistNeg) RETURN
       end if

       ! Update time
       Time = Time + Dt
       if(Time > tFinal - 1.0e-8*DtNext) EXIT
    end do

  end subroutine advect_via_poisson_parker
  !============================================================================
  subroutine iterate_poisson_parker(iLine, nX, iShock, CflIn, BSi_I, nSi_I)
    ! Advect via Possion Bracket scheme to the steady state
    ! Diffuse the distribution function at each time step

    use SP_ModGrid, ONLY: State_VIB, D_, U_
    use SP_ModUnit, ONLY: UnitX_, Io2Si_V
    ! INPUTS:
    integer, intent(in) :: iLine, iShock ! Indices of line and shock
    integer, intent(in) :: nX            ! Number of meshes along s_L axis
    real,    intent(in) :: CflIn         ! Input CFL number
    ! Input variables for diffusion
    real,    intent(in) :: BSi_I(nX), nSi_I(nX)
    real                :: uSi_I(nX-1), DsSi_I(nX-1)
    ! Volume_G: global space volume = product of distance in each dimension
    real :: Volume_G(0:nP+1, 0:nX+1)
    ! VolumeX_I: geometric volume = distance between two geometric faces
    real :: VolumeX_I(0:nX+1)
    ! u/B variable at face
    real :: uOverBSiFace_I(-1:nX+1)
    ! Hamiltonian at cell face
    real :: Hamiltonian_N(-1:nP+1, -1:nX+1)
    ! Extended array for distribution function
    real :: VDF_G(-1:nP+2, -1:nX+2)
    ! Advection term
    real :: Source_C(nP, nX)
    ! Time step
    real :: Dt_C(nP, nX)

    ! Now particle-number-conservative advection scheme for steady-state soln.
    character(len=*), parameter:: NameSub = 'iterate_single_poisson'
    !--------------------------------------------------------------------------

    ! In M-FLAMPA DsSi_I(i) is the distance between meshes i and i+1
    DsSi_I(1:nX-1) = State_VIB(D_, 1:nX-1, iLine)*Io2Si_V(UnitX_)
    ! Initialize arrays
    VolumeX_I(2:nX-1) = 0.5*(DsSi_I(2:nX-1) + DsSi_I(1:nX-2))/BSi_I(2:nX-1)
    VolumeX_I(1)      = DsSi_I(1)/BSi_I(1)
    VolumeX_I(0)      = VolumeX_I(1)
    VolumeX_I(nX)     = DsSi_I(nX-1)/BSi_I(nX)
    VolumeX_I(nX+1)   = VolumeX_I(nX)
    ! Phase volume: initial values
    Volume_G          = spread(VolumeP_I, DIM=2, NCOPIES=nX+2)* &
         spread(VolumeX_I, DIM=1, NCOPIES=nP+2)

    ! Calculate u/B = \vec{u}*\vec{B} / |\vec{B}|**2 at cell face
    ! u/B with the index of "i" is the value at the face between
    ! the mesh "i" and "i+1"
    ! uSi_I with the index of "i" is the value at the face between
    ! the mesh "i" and "i+1"
    uSi_I(1:nX-1) = State_VIB(U_, 1:nX-1, iLine)
    ! Average 1/B and multiply by uSi
    uOverBSiFace_I(1:nX-1) = (0.50/BSi_I(2:nX) + 0.50/BSi_I(1:nX-1))*&
         uSi_I(1:nX-1)
    uOverBSiFace_I(0 )     = uSi_I(1)/BSi_I(1)
    uOverBSiFace_I(-1)     = uOverBSiFace_I(0)
    uOverBSiFace_I(nX)     = uSi_I(nX-1)/BSi_I(nX)
    uOverBSiFace_I(nX+1)   = uOverBSiFace_I(nX)

    ! Hamiltonian = -(u/B)*(p**3/3) at cell face, for {s_L, p**3/3}
    Hamiltonian_N          = -spread(Momentum3_I, DIM=2, NCOPIES=nX+3)* &
         spread(uOverBSiFace_I, DIM=1, NCOPIES=nP+3)

    ! Update bc for at minimal energy, at nP = 0
    call set_momentum_bc(iLine, nX, nSi_I, iShock)
    ! Update Bc for VDF
    call set_VDF(iLine, nX, VDF_G)
    call explicit(nP, nX, VDF_G, Volume_G, Source_C,   &
         Hamiltonian12_N=Hamiltonian_N, CFLIn=CflIn,   &
         IsSteadyState=.true., DtOut_C=Dt_C)

    ! Update velocity distribution function
    Distribution_CB(1:nP, 1, iStart:nX, iLine) = &
         Distribution_CB(1:nP, 1, iStart:nX, iLine) + Source_C(1:nP, iStart:nX)
    ! Check if the VDF includes negative values
    call check_dist_neg(NameSub, 1, nX, iLine)
    if(IsDistNeg) RETURN

    ! Diffuse the distribution function
    if(UseDiffusion) then
       if(UseUpperEndBc) then
          if(UseLowerEndBc) then
             ! with lower or upper end BCs
             call diffuse_distribution(iLine, nX, iShock, Dt_C,    &
                  nSi_I, BSi_I, LowerEndSpectrum_I=VDF_G(1:nP, 0), &
                  UpperEndSpectrum_I=VDF_G(1:nP, nX+1))
          else
             ! with upper end BC but no lower end BC
             call diffuse_distribution(iLine, nX, iShock, Dt_C,    &
                  nSi_I, BSi_I, UpperEndSpectrum_I=VDF_G(1:nP, nX+1))
          end if
       else
          if(UseLowerEndBc) then
             ! with lower end BC but no upper end BC
             call diffuse_distribution(iLine, nX, iShock, Dt_C,    &
                  nSi_I, BSi_I, LowerEndSpectrum_I=VDF_G(1:nP, 0))
          else
             ! with no lower or upper end BCs
             call diffuse_distribution(iLine, nX, iShock, Dt_C, nSi_I, BSi_I)
          end if
       end if
       ! Check if the VDF includes negative values after diffusion
       call check_dist_neg(NameSub//' after diffusion', 1, nX, iLine)
       if(IsDistNeg) RETURN
    end if

  end subroutine iterate_poisson_parker
  !============================================================================
  subroutine init_states_poisson_focused(iLine, nX, DtFull)
    ! Calculate time derivatives from the input files, used for focused
    ! transport equation, if UseBetatron or UseInertialForce is .true.

    use SP_ModGrid, ONLY: UOld_, U_, State_VIB
    ! Line number and number along grid axis
    integer, intent(in) :: iLine, nX
    ! Time difference over the entire loop: DtFull
    real,    intent(in) :: DtFull
    ! Calculate values with time derivatives: Use linear interpolation later
    ! to get the data of each time step from every to consecutive files
    !--------------------------------------------------------------------------
    if(.not.UseInertialForce) RETURN

    ! Now we use the inertial force in focused transport equation
    ! Calculate b*Du/Dt, for Dt in nProgress, i.e., DtFull
    bDuDt_C(1:nX) = (State_VIB(U_, 1:nX, iLine) - &
         State_VIB(UOld_, 1:nX, iLine))/DtFull

  end subroutine init_states_poisson_focused
  !============================================================================
  subroutine advect_via_poisson_focused(iLine, nX, iShock, &
       tFinal, CflIn, nOldSi_I, nSi_I, BOldSi_I, BSi_I)
    ! Advect via multiple Possion Bracket scheme for the focused transport
    ! equation considering pitch angle, may including: adiabatic focusing
    ! (conservation of the magnetic moment), cooling or accleration (with
    ! first-order Fermi accleration, betatron acceleration, by setting
    ! UseBetatron, and inertial forces, by setting UseInertialForce), and
    ! pitch-angle scattering terms, with (p**3/3, mu, s_L).
    ! Here, we diffuse the distribution function at each time step.

    use ModConst,           ONLY: cProtonMass
    use SP_ModDistribution, ONLY: DeltaMu, Mu_I, MuFace_I, SpeedSi_I
    ! INPUTS:
    integer, intent(in) :: iLine, iShock ! Indices of line and shock
    integer, intent(in) :: nX            ! Number of meshes along s_L axis
    real,    intent(in) :: tFinal        ! Time interval to advance through
    real,    intent(in) :: CflIn         ! Input CFL number
    ! Input variables at last and next time steps
    real,    intent(in) :: nOldSi_I(nX), nSi_I(nX), BOldSi_I(nX), BSi_I(nX)
    ! Extended arrays for implementation of the Poisson bracket scheme
    ! Loop variables
    integer :: iP, iMu, iX
    ! For time derivatives and interpolations
    real, dimension(nX) :: InvBOldSi_I, InvBSi_I, dInvBSiDt_C
    ! 1/B at the cell center and face center
    real :: InvBSi_C(nX), InvBFaceSi_C(0:nX)
    ! ------------ Volumes ------------
    ! VolumeXStart_I: geometric volume when the subroutine starts
    ! VolumeXEnd_I: geometric volume when the subroutine ends
    ! dVolumeXDt_I: time derivative of geometric volume
    real, dimension(0:nX+1) :: VolumeXStart_I, VolumeXEnd_I, dVolumeXDt_I
    ! VolumeStart_G: total control volume at the very beginning
    ! Volume_G: total control volume at the end of each iteration
    ! dVolumeDt_G: total control time derivative
    real, dimension(0:nP+1, 0:nMu+1, 0:nX+1) :: &
         VolumeStart_G, Volume_G, dVolumeDt_G
    ! ------------ Hamiltonian functions ------------
    ! Poisson bracket with regard to the first and second vars
    ! considering the case when there are more than one Poisson bracket
    ! and when there is a Poisson bracket with respect to the time
    real :: dHamiltonian01_FX(-1:nP+1,  0:nMu+1, 0:nX+1)
    real :: dHamiltonian02_FY( 0:nP+1, -1:nMu+1, 0:nX+1)
    real :: Hamiltonian2_N( 0:nP+1, -1:nMu+1, -1:nX+1)
    real :: Hamiltonian3_N(-1:nP+1, -1:nMu+1,  0:nX+1)
    ! Extended array for distribution function
    real    :: VDF_G(-1:nP+2, -1:nMu+2, -1:nX+2)
    ! Advection term
    real    :: Source_C(nP, nMu, nX)
    ! Time, ranging from TimeStart to tFinal
    real    :: Time
    ! Inverse of tFinal (DtProgress)
    real    :: InvtFinal
    ! Time step
    real    :: Dt
    ! Prediction of next time step:
    real    :: DtNext
    ! Mark whether this is the last run:
    logical :: IsExit

    ! Now this is the particle-number-conservative advection scheme with mu
    character(len=*), parameter:: NameSub = 'advect_via_triple_poisson'
    !--------------------------------------------------------------------------
    IsDistNeg = .false.
    InvtFinal = 1.0/tFinal

    ! Initialize arrays
    ! Geometric volume: for DeltaS/B, with 1 ghost point at the boundary
    ! Start volume
    VolumeXStart_I(1:nX) = 1.0/nOldSi_I
    VolumeXStart_I(0)    = VolumeXStart_I(1)
    VolumeXStart_I(nX+1) = VolumeXStart_I(nX)
    ! End volume
    VolumeXEnd_I(1:nX)   = 1.0/nSi_I
    VolumeXEnd_I(0)      = VolumeXEnd_I(1)
    VolumeXEnd_I(nX+1)   = VolumeXEnd_I(nX)
    ! Time derivative
    dVolumeXDt_I         = (VolumeXEnd_I - VolumeXStart_I)*InvtFinal

    ! Magnetic field: cell-centered 1/B, D(1/B)/Dt, and face-centered 1/B
    InvBOldSi_I = 1.0/BOldSi_I
    InvBSi_I    = 1.0/BSi_I
    dInvBSiDt_C = (InvBSi_I - InvBOldSi_I)*InvtFinal
    InvBSi_C    = 0.5*(InvBOldSi_I + InvBSi_I)
    InvBFaceSi_C(1:nX-1) = (InvBSi_C(2:nX) + InvBSi_C(1:nX-1))*0.5
    InvBFaceSi_C(0 )     = InvBSi_C(1 ) - (InvBSi_C(2 ) - InvBSi_C(1   ))*0.5
    InvBFaceSi_C(nX)     = InvBSi_C(nX) + (InvBSi_C(nX) - InvBSi_C(nX-1))*0.5

    ! Calculate total control volumes and time-related Hamiltonian functions
    do iX = 0, nX+1
       ! Total control volume at the very beginning
       VolumeStart_G(:, :, iX) = VolumeXStart_I(iX)*DeltaMu* &
            spread(VolumeP_I, DIM=2, NCOPIES=nMu+2)
       ! Time derivative of total control volume
       dVolumeDt_G(:, :, iX) = dVolumeXDt_I(iX)*DeltaMu* &
            spread(VolumeP_I, DIM=2, NCOPIES=nMu+2)

       ! Calculate 1st Hamiltonian function used in the time-dependent Poisson
       ! bracket: {f_jk; (p**3/3)*(DeltaS/B)}_{tau, p**3/3} => \tilde\DeltaH
       dHamiltonian01_FX(:, 1:nMu, iX) = -dVolumeXDt_I(iX)*DeltaMu* &
            3.0*spread(Mu_I**2, DIM=1, NCOPIES=nP+3)* &
            spread(Momentum3_I, DIM=2, NCOPIES=nMu)
       ! Calculate 2nd Hamiltonian function used in the time-dependent Poisson
       ! bracket: {f_jk; mu*(1-mu**2)*(DeltaS/B)}_{tau, mu} => \tilde\DeltaH
       dHamiltonian02_FY(:, 0:nMu, iX) = -dVolumeXDt_I(iX)* &
            spread(MuFace_I*(1.0-MuFace_I**2), DIM=1, NCOPIES=nP+2)* &
            spread(VolumeP_I, DIM=2, NCOPIES=nMu+1)
    end do
    ! Bcs for time-related Hamiltonian functions: {f_jk; H}_{tau, ...}
    dHamiltonian01_FX(:,     0, :) = dHamiltonian01_FX(:,   1, :)
    dHamiltonian01_FX(:, nMu+1, :) = dHamiltonian01_FX(:, nMu, :)
    dHamiltonian02_FY(:,    -1, :) = dHamiltonian02_FY(:,   0, :)
    dHamiltonian02_FY(:, nMu+1, :) = dHamiltonian02_FY(:, nMu, :)

    ! Time initialization
    Time = 0.0

    ! Here we would like to get the first trial of DtNext
    ! Update Bc for VDF at minimal energy, at nP = 0
    call set_momentum_bc(iLine, nX, nSi_I, iShock)
    call set_VDF(iLine, nX, VDF_G) ! Set the VDF first
    call advance_triple_poisson    ! Now we get DtNext

    ! Advect by the triple-Poisson-bracket scheme
    do
       ! Update Time step
       if(DtNext > tFinal - Time) then
          Dt = tFinal - Time
          IsExit = .true.          ! Last step
       else
          Dt = DtNext
          IsExit = .false.         ! Intermediate steps
       end if

       ! Update Bc for at minimal energy, at nP = 0
       call set_momentum_bc(iLine, nX, nSi_I, iShock)
       ! Update Bc for VDF
       call set_VDF(iLine, nX, VDF_G)
       ! Advance by triple-Poisson-bracket scheme
       call advance_triple_poisson(DtIn=Dt)

       if(IsExit) then
          ! This step is the last step
          Distribution_CB(1:nP, 1:nMu, iStart:nX, iLine) =      &
               Distribution_CB(1:nP, 1:nMu, iStart:nX, iLine) + &
               Source_C(1:nP, 1:nMu, iStart:nX)
       else
          ! This step is not the last step
          ! Update VDF_G considering the time-dependent Volume_G

          Source_C = Source_C*Volume_G(1:nP, 1:nMu, 1:nX)
          ! Update DeltaSOverB_C for the calculation of volume
          Volume_G = VolumeStart_G + (Time+Dt)*dVolumeDt_G

          ! Update VDF_G to CURRENT time: no BCs for (1:nQ, 1:nP, iStart:nR)
          Distribution_CB(1:nP, 1:nMu, iStart:nX, iLine) =      &
               Distribution_CB(1:nP, 1:nMu, iStart:nX, iLine) + &
               Source_C(1:nP, 1:nMu, iStart:nX)/ &
               Volume_G(1:nP, 1:nMu, iStart:nX)
       end if
       ! Check if the VDF includes negative values after scattering
       call check_dist_neg(NameSub, 1, nX, iLine)
       if(IsDistNeg) RETURN

       ! ------------ Future Work ------------
       ! This is the first version of draft implementing multi-Poisson-bracket
       ! scheme in SP/MFLAMPA made by Weihao Liu. It needs more development
       ! in the future. Here, we will also include some scattering functions.
       ! (Clearly, UseDiffusion and diffuse_distrition is missing here.)
       ! Moreover, the structure should be optimized, with some bugs fixed.
       ! One can refer to test_multi_poisson in share/Library/test/
       ! One should specify the TypeScatter in the PARAM.in file.
       ! We will work on this further later since it is more advanced.
       ! ------------ Thank you! ------------
       ! Check if the VDF includes negative values
       call check_dist_neg(NameSub//' after scattering', 1, nX, iLine)
       if(IsDistNeg) RETURN

       ! Update time
       Time = Time + Dt
       if(Time > tFinal - 1.0e-8*DtNext) EXIT
    end do

  contains
    !==========================================================================
    subroutine calc_triple_hamiltonian_2
      ! Calculate the 2nd Hamiltonian function at each fixed time:
      ! (mu**2-1)*v/(2B) at face of s_L and mu, regarding to s_L and mu

      ! Considering the law of relativity, v = p*c**2/sqrt(p**2+(m*c**2)**2), 
      ! calculated as a function of p. Note that momentum (or speed) in this
      ! term is non-related to the Lagrangian coordinates (s_L, mu), so it
      ! should be cell-centered, by using SpeedSi_I.

      ! Calculate the 2nd hamiltonian = (mu**2-1)*v/(2B) => \tilde\deltaH
      do iX = 0, nX
         Hamiltonian2_N(:, 0:nMu, iX) = 0.5*InvBFaceSi_C(iX)*( &
              spread(MuFace_I**2, DIM=1, NCOPIES=nP+2) - 1.0)* &
              spread(SpeedSi_I*VolumeP_I, DIM=2, NCOPIES=nMu+1)
      end do

      ! Boundary condition of Hamiltonian function
      Hamiltonian2_N(:,    -1, :) = Hamiltonian2_N(:,     1, :)
      Hamiltonian2_N(:, nMu+1, :) = Hamiltonian2_N(:, nMu-1, :)

    end subroutine calc_triple_hamiltonian_2
    !==========================================================================
    subroutine calc_triple_hamiltonian_3
      ! Calculate the 3rd Hamiltonian function at each fixed time:
      ! (1-mu**2)/2*[ mu*(p**3/3)*(3\vec{b}\vec{b}:\nabla\vec{u} - div\vec{u})
      ! + ProtonMass*p**2*(\vec{b}*d\vec{u}/dt) ], regarding to p**3/3 and mu,
      ! so the coordinates of p**3/3 and mu are at face, and s_L is at center.
      ! Here we list the variables in the analytical function:
      ! (3\vec{b}\vec{b}:\nabla\vec{u} - div\vec{u}) = d(ln(B*ds**2))/dt
      ! ProtonMass = cRmeProtonGeV, in the unit of GeV/c**2
      ! \vec{b}*d\vec{u}/dt = bDuDt_C

      !------------------------------------------------------------------------

     !  do iX = 1, nX
     !     do iMu = 0, nMu
     !        Hamiltonian3_N(:, iMu, iX) = 0.5*(1.0 - MuFace_I(iMu)**2)* &
     !             (MuFace_I(iMu)*Momentum3_I*dLnBDeltaS2Dt_C(iX) +      &
     !             cProtonMass*(Momentum3_I*3.0)**(2.0/3.0)*bDuDt_C(iX))
     !     end do
     !     ! Here, what we use actually is: \tilde\deltaH
     !     Hamiltonian3_N(:, 0:nMu, iX) = &
     !          Hamiltonian3_N(:, 0:nMu, iX)*DeltaSOverB_C(iX)
     !  end do

      ! Boundary condition of Hamiltonian function
      Hamiltonian3_N(:,     :,    0) = Hamiltonian3_N(:,     :,  1)
      Hamiltonian3_N(:,     :, nX+1) = Hamiltonian3_N(:,     :, nX)
      Hamiltonian3_N(:,    -1,    :) = Hamiltonian3_N(:,     1,  :)
      Hamiltonian3_N(:, nMu+1,    :) = Hamiltonian3_N(:, nMu-1,  :)

    end subroutine calc_triple_hamiltonian_3
    !==========================================================================
    subroutine advance_triple_poisson(DtIn)
      ! advance by the triple-Poisson-bracket scheme for each time step

      real, optional, intent(inout) :: DtIn ! Input time step
      ! Get DeltaSOverB_C at current time: for calculating
      ! the total control volume and Hamiltonian functions
      !------------------------------------------------------------------------

      ! Clean 2nd and 3rd Hamiltonian functions
      Hamiltonian2_N = 0.0
      Hamiltonian3_N = 0.0

      ! Calculate 2nd Hamiltonian function used in the time-dependent
      ! poisson bracket: {f_jk, (mu**2-1)*v/(2B)}_{s_L, mu}
      call calc_triple_hamiltonian_2
      ! Calculate 3rd Hamiltonian function used in the time-dependent
      ! poisson bracket: {f_jk, (1-mu**2)/2 * [ mu*(p**3/3)*
      ! d(ln(B*ds**2))/dt + ProtonMass*p**2*bDuDt_C ]}_{p**3/3, mu}
      call calc_triple_hamiltonian_3

      ! Update the total control volume
      Volume_G = VolumeStart_G + Time*dVolumeDt_G

      ! Here we have three Lagrangian coordinates: (p**3/3, mu, s_L)
      if(present(DtIn)) then
         call explicit(nP, nMu, nX, VDF_G, Volume_G, Source_C, &
              dHamiltonian01_FX = dHamiltonian01_FX,           &
              Hamiltonian12_N = Hamiltonian3_N,                &
              Hamiltonian23_N = -Hamiltonian2_N,               &
              dVolumeDt_G = dVolumeDt_G,                       &
              DtIn = DtIn, DtOut = DtNext, CFLIn = CflIn)
      else
         ! Not present DtIn for the first call
         call explicit(nP, nMu, nX, VDF_G, Volume_G, Source_C, &
              dHamiltonian01_FX = dHamiltonian01_FX,           &
              Hamiltonian12_N = Hamiltonian3_N,                &
              Hamiltonian23_N = -Hamiltonian2_N,               &
              dVolumeDt_G = dVolumeDt_G,                       &
              DtOut = DtNext, CFLIn = CflIn)
      end if

    end subroutine advance_triple_poisson
    !==========================================================================
  end subroutine advect_via_poisson_focused
  !============================================================================
end module SP_ModAdvancePoisson
!==============================================================================
