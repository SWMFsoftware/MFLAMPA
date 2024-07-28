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
  use SP_ModDistribution, ONLY: nP, nMu, Distribution_CB, &
       dLogP, VolumeP_I, Momentum3_I, DeltaMu, Mu_I,      &
       IsDistNeg, check_dist_neg
  use SP_ModDiffusion,    ONLY: UseDiffusion, diffuse_distribution
  use ModUtilities,       ONLY: CON_stop
  use ModPoissonBracket,  ONLY: explicit
  implicit none

  PRIVATE ! Except

  SAVE
  public:: advect_via_single_poisson ! Time-accurate advance through given Dt
  public:: iterate_single_poisson    ! Local time-stepping for steady states
  public:: advect_via_double_poisson ! Advance with 2 Poisson brackets
  public:: advect_via_triple_poisson ! Advance with 3 Poisson brackets
  public:: init_states_for_poisson   ! Initialize the states for multi_poisson

  ! \Deltas/b, ln(B\deltas^2) at Old time
  real, public, dimension(nVertexMax) :: DeltaSOverBOld_C, LnBDeltaS2Old_C
  ! \Deltas/b, ln(B\deltas^2) at New time
  real, public, dimension(nVertexMax) :: DeltaSOverBNew_C, LnBDeltaS2New_C
  ! Time-derivatives:
  real, public, dimension(nVertexMax) :: DeltaSOverB_C, &
       dDeltaSOverBDt_C, dLnBdeltaS2Dt_C, bDuDt_C
contains
  !============================================================================
  subroutine advect_via_single_poisson(iLine, nX, iShock, &
       tFinal, CflIn, nOldSi_I, nSi_I, BSi_I)
    ! advect via Possion Bracket scheme
    ! diffuse the distribution function at each time step

    ! INPUTS:
    integer, intent(in) :: iLine, iShock ! Indices of line and shock
    integer, intent(in) :: nX            ! Number of meshes along s_L axis
    real,    intent(in) :: tFinal        ! Time interval to advance through
    real,    intent(in) :: CflIn         ! Input CFL number
    ! Input variables for diffusion
    real,    intent(in) :: nOldSi_I(nX), nSi_I(nX), BSi_I(nX)
    ! Extended arrays for implementation of the Poisson bracket scheme
    ! VolumeXStart_I: geometric volume when the subroutine starts
    ! VolumeXEnd_I: geometric volume when the subroutine ends
    ! dVolumeXDt_I: time derivative of geometric volume
    real, dimension(0:nX+1) :: VolumeXStart_I, VolumeXEnd_I, dVolumeXDt_I
    ! Volume_G: total control volume at the end of each iteration
    ! VolumeOld_G: total control volume at the end of each iteration
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
    ! Geometric volume: use 1 ghost point at each side of the boundary
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
    ! calculate: dHamiltonian/dVolumeSubX
    dHamiltonian01_FX    = -spread(Momentum3_I, DIM=2, NCOPIES=nX+2)* &
         spread(dVolumeXDt_I, DIM=1, NCOPIES=nP+3)
    Time   = 0.0
    ! Trial timestep
    DtNext = CflIn/maxval(abs(dVolumeXDt_I)/ &
         max(VolumeXEnd_I, VolumeXStart_I))/(3.0*dLogP)

    ! Update Bc for VDF at minimal energy, at nP = 0
    call set_momentum_bc(iLine, nX, nSi_I, iShock)
    ! Advection by the single-Poisson-bracket scheme
    do
       ! Update Time step
       Dt = min(DtNext, tFinal - Time)
       ! Update Volumes
       VolumeOld_G = Volume_G
       Volume_G    = VolumeOld_G + Dt*dVolumeDt_G
       ! Update Bc
       call set_VDF(iLine, nX, VDF_G)

       ! Advance by Poisson bracket scheme
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

  end subroutine advect_via_single_poisson
  !============================================================================
  subroutine iterate_single_poisson(iLine, nX, iShock, CflIn, BSi_I, nSi_I)
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
    real :: uOverBNodeSi_I(-1:nX+1)
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

    ! Calculate u/B = \vec{u}*\vec{B} / |\vec{B}|^2 at cell face
    ! u/B with the index of "i" is the value at the face between
    ! the mesh "i" and "i+1"
    ! uSi_I with the index of "i" is the value at the face between
    ! the mesh "i" and "i+1"
    uSi_I(1:nX-1) = State_VIB(U_, 1:nX-1, iLine)
    ! Average 1/B and multiply by uSi
    uOverBNodeSi_I(1:nX-1) = (0.50/BSi_I(2:nX) + 0.50/BSi_I(1:nX-1))*&
         uSi_I(1:nX-1)
    uOverBNodeSi_I(0 )     = uSi_I(1)/BSi_I(1)
    uOverBNodeSi_I(-1)     = uOverBNodeSi_I(0)
    uOverBNodeSi_I(nX)     = uSi_I(nX-1)/BSi_I(nX)
    uOverBNodeSi_I(nX+1)   = uOverBNodeSi_I(nX)

    ! Hamiltonian = -(u/B)*(p**3/3) at cell face, for {s_L, p^3/3}
    Hamiltonian_N          = -spread(Momentum3_I, DIM=2, NCOPIES=nX+3)* &
         spread(uOverBNodeSi_I, DIM=1, NCOPIES=nP+3)

    ! Update bc for at minimal energy, at nP = 0
    call set_momentum_bc(iLine, nX, nSi_I, iShock)
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

  end subroutine iterate_single_poisson
  !============================================================================
  subroutine advect_via_double_poisson(iLine, nX, iShock, &
       tFinal, CflIn, nSi_I, BSi_I)
    ! Advect via multiple Possion Bracket scheme for the focused transport
    ! equation considering pitch angle, including: adiabatic cooling or
    ! accleration, and pitch-angle scattering terms, but without adiabatic
    ! focusing term (partial f/pattial mu) for the pitch angle.
    ! Here, we diffuse the distribution function at each time step.

    ! INPUTS:
    integer, intent(in) :: iLine, iShock ! Indices of line and shock
    integer, intent(in) :: nX            ! Number of meshes along s_L axis
    real,    intent(in) :: tFinal        ! Time interval to advance through
    real,    intent(in) :: CflIn         ! Input CFL number
    ! Input variables for diffusion
    real,    intent(in) :: nSi_I(nX), BSi_I(nX)
    ! Extended arrays for implementation of the Poisson bracket scheme
    ! Loop variables
    integer :: iX, iMu, iP
    ! ------------ Control volume ------------
    ! Volume_G: total control volume at the end of each iteration
    ! dVolumeDt_G: total control time derivative
    real, dimension(0:nP+1, 0:nMu+1, 0:nX+1) :: Volume_G
    real :: dVolumeDt_G(0:nP+1, 0:nMu+1, 0:nX+1)
    ! ------------ Hamiltonian functions ------------
    ! Poisson bracket with regard to the first and second vars
    ! considering the case when there are more than one Poisson bracket
    ! and when there is a Poisson bracket with respect to the time
    real :: dHamiltonian01_FX(0:nP+1, 0:nMu+1, -1:nX+1)
    real :: Hamiltonian2_N(0:nP+1, -1:nMu+1, -1:nX+1)
    ! Extended array for distribution function
    real    :: VDF_G(-1:nP+2, -1:nMu+2, -1:nX+2)
    ! Advection term
    real    :: Source_C(nP, nMu, nX)
    ! Time, ranging from TimeStart to tFinal
    real    :: Time
    ! Time step
    real    :: Dt
    ! Prediction of next time step:
    real    :: DtNext
    ! Mark whether this is the last run:
    logical :: IsExit

    ! Now this is the particle-number-conservative advection scheme with \mu
    character(len=*), parameter:: NameSub = 'advect_via_double_poisson'
    !--------------------------------------------------------------------------
    IsDistNeg = .false.

    ! Calculate time derivative of total control volume
    do iX = 1, nX
       dVolumeDt_G(:, :, iX) = dDeltaSOverBDt_C(iX)*DeltaMu* &
            spread(VolumeP_I, DIM=2, NCOPIES=nMu+2)
    end do
    ! BCs of the time derivative of total control volume
    dVolumeDt_G(:, :,    0) = dVolumeDt_G(:, :,  1)
    dVolumeDt_G(:, :, nX+1) = dVolumeDt_G(:, :, nX)
    ! Time initialization
    Time = 0.0

    ! Here we would like to get the first trial of DtNext
    call set_VDF(iLine, nX, VDF_G) ! Set the VDF first
    call advance_double_poisson    ! Now we get DtNext

    ! Advect by the double-Poisson-bracket scheme
    do
       ! Update Time step
       if(DtNext > tFinal - Time) then
          Dt = tFinal - Time
          IsExit = .true.          ! Last step
       else
          Dt = DtNext
          IsExit = .false.         ! Intermediate steps
       end if

       ! Update Bc for VDF at minimal energy, at nP = 0
       call set_momentum_bc(iLine, nX, nSi_I, iShock)
       ! Update Bc for VDF
       call set_VDF(iLine, nX, VDF_G)
       ! Advance by multi-Poisson-bracket scheme
       call advance_double_poisson(DtIn=Dt)

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
          call update_states_for_poisson(iLine, nX, Time)
          ! Calculate the total control volume
          do iX = 1, nX
             Volume_G(:, :, iX) = DeltaSOverB_C(iX)*DeltaMu* &
                  spread(VolumeP_I, DIM=2, NCOPIES=nMu+2)
          end do

          ! Update VDF_G to CURRENT time: no BCs for (1:nQ, 1:nP, iStart:nR)
          Distribution_CB(1:nP, 1:nMu, iStart:nX, iLine) =      &
               Distribution_CB(1:nP, 1:nMu, iStart:nX, iLine) + &
               Source_C(1:nP, 1:nMu, iStart:nX)/ &
               Volume_G(1:nP, 1:nMu, iStart:nX)
       end if
       ! Check if the VDF includes negative values
       call check_dist_neg(NameSub, 1, nX, iLine)
       if(IsDistNeg) RETURN

       ! !! SCATTERING and DIFFUSION here !! !

       ! Update time
       Time = Time + Dt
       if(Time > tFinal - 1.0e-8*DtNext) EXIT
    end do

  contains
    !==========================================================================
    subroutine advance_double_poisson(DtIn)
      ! Advance by the double-Poisson-bracket scheme for each time step

      real, optional, intent(inout) :: DtIn ! Input time step
      ! Get DeltaSOverB_C at current time: for calculating
      ! the total control volume and Hamiltonian functions
      !------------------------------------------------------------------------

      call update_states_for_poisson(iLine, nX, Time)
      Hamiltonian2_N = 0.0

      ! Calculate 1st Hamiltonian function used in the time-dependent
      ! poisson bracket: {f_jk, (mu^2 * p^3)*(DeltaS/B)}_{tau, p^3/3}
      call calc_double_hamiltonian_1(nX, dHamiltonian01_FX)
      ! Calculate 2nd Hamiltonian function used in the time-dependent
      ! poisson bracket: {f_jk, (mu^2-1)*v/(2B)}_{s_L, mu}
      call calc_triple_hamiltonian_2(nX, BSi_I, Hamiltonian2_N)

      ! Calculate the total control volume
      do iX = 1, nX
         Volume_G(:, :, iX) = DeltaSOverB_C(iX)*DeltaMu* &
              spread(VolumeP_I, DIM=2, NCOPIES=nMu+2)
      end do

      ! Boundary conditions for total control volume
      Volume_G(:, :,    0) = Volume_G(:, :,  1)
      Volume_G(:, :, nX+1) = Volume_G(:, :, nX)

      ! Here we have three Lagrangian coordinates: p^3/3, mu, s_L
      if(present(DtIn)) then
         call explicit(nP, nMu, nX, VDF_G, Volume_G, Source_C, &
              dHamiltonian01_FX = dHamiltonian01_FX,           &
              Hamiltonian23_N = -Hamiltonian2_N,               &
              dVolumeDt_G = dVolumeDt_G,                       &
              DtIn = DtIn, DtOut = DtNext, CFLIn = CflIn)
      else
         ! Not present DtIn for the first call
         call explicit(nP, nMu, nX, VDF_G, Volume_G, Source_C, &
              dHamiltonian01_FX = dHamiltonian01_FX,           &
              Hamiltonian23_N = -Hamiltonian2_N,               &
              dVolumeDt_G = dVolumeDt_G,                       &
              DtOut = DtNext, CFLIn = CflIn)
      end if

    end subroutine advance_double_poisson
    !==========================================================================
  end subroutine advect_via_double_poisson
  !============================================================================
  subroutine calc_double_hamiltonian_1(nX, dHamiltonian01_FX)
    ! Calculate the 1st Hamiltonian function with time:
    ! mu^2 * p^3/3 * \deltas/B at p^3/3 face, regarding to tau and p^3/3

    integer, intent(in) :: nX        ! Number of s_L grid
    real, intent(inout) :: dHamiltonian01_FX(-1:nP+1, 0:nMu+1, 0:nX+1)
    integer             :: iX, iMu   ! Loop variables
    ! Calculate the first Hamiltonian function for the 2-Poisson-bracket scheme
    !--------------------------------------------------------------------------

    do iX = 1, nX
       dHamiltonian01_FX(:, 1:nMu, iX) = -dDeltaSOverBDt_C(iX)* &
            spread(Momentum3_I, DIM=2, NCOPIES=nMu)*3.0* & ! momentum at face
            spread(Mu_I, DIM=1, NCOPIES=nP+3)              ! mu at cell center
    end do

    ! Calculate the Hamiltonian function used actually：\tilde\deltaH
    dHamiltonian01_FX(:, 1:nMu, 1:nX) = &
         dHamiltonian01_FX(:, 1:nMu, 1:nX)*DeltaMu

    ! Boundary condition of Hamiltonian function
    dHamiltonian01_FX(:,    :,    0) = dHamiltonian01_FX(:,  :,  1)
    dHamiltonian01_FX(:,    :, nX+1) = dHamiltonian01_FX(:,  :, nX)
    dHamiltonian01_FX(:,    0,    :) = dHamiltonian01_FX(:,  1,  :)
    dHamiltonian01_FX(:, nX+1,    :) = dHamiltonian01_FX(:, nX,  :)

  end subroutine calc_double_hamiltonian_1
  !============================================================================
  subroutine advect_via_triple_poisson(iLine, nX, iShock, &
       tFinal, CflIn, nSi_I, BSi_I)
    ! Advect via multiple Possion Bracket scheme for the focused transport
    ! equation considering pitch angle, including: adiabatic focusing,
    ! cooling or accleration, and pitch-angle scattering terms.
    ! Here, we diffuse the distribution function at each time step.

    ! INPUTS:
    integer, intent(in) :: iLine, iShock ! Indices of line and shock
    integer, intent(in) :: nX            ! Number of meshes along s_L axis
    real,    intent(in) :: tFinal        ! Time interval to advance through
    real,    intent(in) :: CflIn         ! Input CFL number
    ! Input variables for diffusion
    real,    intent(in) :: nSi_I(nX), BSi_I(nX)
    ! Loop variables
    integer :: iX, iMu, iP
    ! ------------ Control volume ------------
    ! Volume_G: total control volume at the end of each iteration
    ! dVolumeDt_G: total control time derivative
    real, dimension(0:nP+1, 0:nMu+1, 0:nX+1) :: Volume_G
    real :: dVolumeDt_G(0:nP+1, 0:nMu+1, 0:nX+1)
    ! ------------ Hamiltonian functions ------------
    ! Poisson bracket with regard to the first and second vars
    ! considering the case when there are more than one Poisson bracket
    ! and when there is a Poisson bracket with respect to the time
    real :: dHamiltonian01_FX(0:nP+1, 0:nMu+1, -1:nX+1)
    real :: Hamiltonian2_N(0:nP+1, -1:nMu+1, -1:nX+1)
    real :: Hamiltonian3_N(-1:nP+1, -1:nMu+1, 0:nX+1)
    ! Extended array for distribution function
    real    :: VDF_G(-1:nP+2, -1:nMu+2, -1:nX+2)
    ! Advection term
    real    :: Source_C(nP, nMu, nX)
    ! Time, ranging from TimeStart to tFinal
    real    :: Time
    ! Time step
    real    :: Dt
    ! Prediction of next time step:
    real    :: DtNext
    ! Mark whether this is the last run:
    logical :: IsExit

    ! Now this is the particle-number-conservative advection scheme with mu
    character(len=*), parameter:: NameSub = 'advect_via_triple_poisson'
    !--------------------------------------------------------------------------

    ! Calculate time derivative of total control volume
    do iX = 1, nX
       dVolumeDt_G(:, :, iX) = dDeltaSOverBDt_C(iX)*DeltaMu* &
            spread(VolumeP_I, DIM=2, NCOPIES=nMu+2)
    end do
    ! BCs of the time derivative of total control volume
    dVolumeDt_G(:, :,    0) = dVolumeDt_G(:, :,  1)
    dVolumeDt_G(:, :, nX+1) = dVolumeDt_G(:, :, nX)
    ! Time initialization
    Time = 0.0

    ! Here we would like to get the first trial of DtNext
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

       ! Update bc for at minimal energy, at nP = 0
       call set_momentum_bc(iLine, nX, nSi_I, iShock)
       ! Update Bc for VDF
       call set_VDF(iLine, nX, VDF_G)
       ! Advance by multi-Poisson-bracket scheme
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
          call update_states_for_poisson(iLine, nX, Time)
          ! Calculate the total control volume
          do iX = 1, nX
             Volume_G(:, :, iX) = DeltaSOverB_C(iX)*DeltaMu* &
                  spread(VolumeP_I, DIM=2, NCOPIES=nMu+2)
          end do

          ! Update VDF_G to CURRENT time: no BCs for (1:nQ, 1:nP, iStart:nR)
          Distribution_CB(1:nP, 1:nMu, iStart:nX, iLine) =      &
               Distribution_CB(1:nP, 1:nMu, iStart:nX, iLine) + &
               Source_C(1:nP, 1:nMu, iStart:nX)/ &
               Volume_G(1:nP, 1:nMu, iStart:nX)
       end if
       ! Check if the VDF includes negative values
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

       ! Update time
       Time = Time + Dt
       if(Time > tFinal - 1.0e-8*DtNext) EXIT
    end do

  contains
    !==========================================================================
    subroutine advance_triple_poisson(DtIn)
      ! advance by the triple-Poisson-bracket scheme for each time step

      real, optional, intent(inout) :: DtIn ! Input time step
      ! Get DeltaSOverB_C at current time: for calculating
      ! the total control volume and Hamiltonian functions
      !------------------------------------------------------------------------

      call update_states_for_poisson(iLine, nX, Time)
      Hamiltonian2_N = 0.0
      Hamiltonian3_N = 0.0

      ! Calculate 1st Hamiltonian function used in the time-dependent
      ! poisson bracket: {f_jk, (p^3/3)*(DeltaS/B)}_{tau, p^3/3}
      call calc_triple_hamiltonian_1(nX, dHamiltonian01_FX)
      ! Calculate 2nd Hamiltonian function used in the time-dependent
      ! poisson bracket: {f_jk, (mu^2-1)*v/(2B)}_{s_L, mu}
      call calc_triple_hamiltonian_2(nX, BSi_I, Hamiltonian2_N)
      ! Calculate 3rd Hamiltonian function used in the time-dependent
      ! poisson bracket: {f_jk, (1-mu^2)/2 * [ mu*(p^3/3)*
      ! d(ln(B*ds^2))/dt + ProtonMass*p^2*bDuDt_C ]}_{p^3/3, mu}
      call calc_triple_hamiltonian_3(nX, Hamiltonian3_N)

      ! Calculate the total control volume
      do iX = 1, nX
         Volume_G(:, :, iX) = DeltaSOverB_C(iX)*DeltaMu* &
              spread(VolumeP_I, DIM=2, NCOPIES=nMu+2)
      end do

      ! Boundary conditions for total control volume
      Volume_G(:, :,    0) = Volume_G(:, :,  1)
      Volume_G(:, :, nX+1) = Volume_G(:, :, nX)

      ! Here we have three Lagrangian coordinates: p^3/3, mu, s_L
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
  end subroutine advect_via_triple_poisson
  !============================================================================
  subroutine init_states_for_poisson(iLine, nX, DtProgress, &
       XyzOldSi_DI, XyzSi_DI, BOldSi_I, BSi_I)
    ! Calculate data states from the input files

    use SP_ModGrid, ONLY: X_, Z_, UOld_, U_, State_VIB
    ! Line number and number along grid axis
    integer, intent(in) :: iLine, nX
    ! Time difference between Old (iProgress-1) and New (iProgress) States
    real,    intent(in) :: DtProgress
    ! Given Old and New coordinates
    real,    intent(in) :: XyzOldSi_DI(X_:Z_, nX), XyzSi_DI(X_:Z_, nX)
    ! Old and New |B| states, and we include steepen_shock to make updates
    real,    intent(in) :: BOldSi_I(nX), BSi_I(nX)
    ! Inverse of DtProgress, the reduced time step
    real :: InvDtProgress
    ! \delta_s for consecutive points, at old and new states
    real :: DeltaSOld_I(nX), DeltaSNew_I(nX)

    !--------------------------------------------------------------------------
    InvDtProgress = 1.0/DtProgress

    ! Calculate values at OLD time
    DeltaSOld_I = calc_DeltaS_I(XyzOldSi_DI)
    ! Calculate \DeltaS/B at grid center
    DeltaSOverBOld_C(1:nX) = DeltaSOld_I/BOldSi_I(1:nX)
    ! Calculate ln(B*\DeltaS^2) at grid center
    LnBDeltaS2Old_C(1:nX) = log(BOldSi_I(1:nX)) + 2.0*log(DeltaSOld_I)

    ! Calculate values at NEW time
    DeltaSNew_I = calc_DeltaS_I(XyzSi_DI)
    ! Calculate \DeltaS/B at grid center
    DeltaSOverBNew_C(1:nX) = DeltaSNew_I/BSi_I(1:nX)
    ! Calculate ln(B*\DeltaS^2) at grid center
    LnBDeltaS2New_C(1:nX) = log(BSi_I(1:nX)) + 2.0*log(DeltaSNew_I)

    ! Calculate the time-derivative physical quantities
    ! Calculate \DeltaS/B time derivative = D[delta(s_L)/B]/Dt
    dDeltaSOverBDt_C(1:nX) = (DeltaSOverBNew_C(1:nX) - &
         DeltaSOverBOld_C(1:nX))*InvDtProgress
    ! Calculate Dln(B\DeltaS^2)/Dt
    dLnBdeltaS2Dt_C(1:nX)  = (LnBDeltaS2New_C(1:nX) - &
         LnBDeltaS2Old_C(1:nX))*InvDtProgress
    ! Calculate b*Du/Dt
    bDuDt_C(1:nX)          = (State_VIB(U_, 1:nX, iLine) - &
         State_VIB(UOld_, 1:nX, iLine))*InvDtProgress

  contains
    !==========================================================================
    function calc_DeltaS_I(XyzSi_DI)
      ! Calculate DeltaS_I at given time-interpolation coefficient, Alpha

      real, intent(in) :: XyzSi_DI(X_:Z_, nX) ! Current coordinates
      real :: calc_DeltaS_I(nX)               ! Output results: DeltaS_I
      real :: MidPoint_IB(X_:Z_, nX)          ! Midpoints
      ! Get DeltaS_I at current time: for calculating other old and new states
      !------------------------------------------------------------------------

      ! Calculate midpoints
      MidPoint_IB(:, 1:nX-1) = (XyzSi_DI(:, 2:nX) + XyzSi_DI(:, 1:nX-1))*0.5

      ! Calculate DeltaS_I
      calc_DeltaS_I(2:nX-1) = sqrt(sum((MidPoint_IB(:, 2:nX-1) &
           - MidPoint_IB(:, 1:nX-2))**2, dim=1))

      ! Linear interpolate the deltas such that there will be nQ deltas
      calc_DeltaS_I(1 ) = 2.0*sqrt(sum(( &
           MidPoint_IB(:, 1) - XyzSi_DI(:,  1))**2))
      calc_DeltaS_I(nX) = 2.0*sqrt(sum(( &
           MidPoint_IB(:, nX-1) - XyzSi_DI(:, nX))**2))

    end function calc_DeltaS_I
    !==========================================================================
  end subroutine init_states_for_poisson
  !============================================================================
  subroutine update_states_for_poisson(iLine, nX, Time)
    ! Update states according to the current time

    integer, intent(in) :: iLine  ! Index of the current field line
    integer, intent(in) :: nX     ! Number of grid along s_L axis
    real,    intent(in) :: Time   ! Current time
    ! Calculate values for CURRENT time: here we use linear interpolation
    ! to get the data of each time step from every to consecutive files
    !--------------------------------------------------------------------------

    ! Update \deltas/B
    DeltaSOverB_C(1:nX) = DeltaSOverBOld_C(1:nX) + dDeltaSOverBDt_C(1:nX)*Time
  end subroutine update_states_for_poisson
  !============================================================================
  subroutine calc_triple_hamiltonian_1(nX, dHamiltonian01_FX)
    ! Calculate the 1st Hamiltonian function with time:
    ! p^3/3*\deltas/B at cell face of p^3/3, regarding to tau and p^3/3

    integer, intent(in) :: nX        ! Number of s_L grid
    real, intent(inout) :: dHamiltonian01_FX(-1:nP+1, 0:nMu+1, 0:nX+1)
    integer             :: iX, iMu   ! Loop variables
    ! Calculate the 1st Hamiltonian function for the 3-Poisson-bracket scheme
    !--------------------------------------------------------------------------

    do iX = 1, nX
       dHamiltonian01_FX(:, :, iX) = -dDeltaSOverBDt_C(iX)* &
            spread(Momentum3_I, DIM=2, NCOPIES=nMu+2)
    end do

    ! Calculate the Hamiltonian function used actually：\tilde\deltaH
    dHamiltonian01_FX(0:nP, 0:nMu+1, 1:nX) =    &
         dHamiltonian01_FX(0:nP, 0:nMu+1, 1:nX)*DeltaMu

    ! Boundary condition of Hamiltonian function
    dHamiltonian01_FX(:, :,    0) = dHamiltonian01_FX(:, :,  1)
    dHamiltonian01_FX(:, :, nX+1) = dHamiltonian01_FX(:, :, nX)

  end subroutine calc_triple_hamiltonian_1
  !============================================================================
  subroutine calc_triple_hamiltonian_2(nX, BSi_I, Hamiltonian2_N)
    ! Calculate the 2nd Hamiltonian function at each fixed time:
    ! (mu^2-1)*v/(2B) at face of s_L and mu, regarding to s_L and mu

    use ModConst,           ONLY: cProtonMass, cLightSpeed
    use SP_ModDistribution, ONLY: Momentum_I, MuFace_I
    integer, intent(in) :: nX        ! Number of s_L grid
    real, intent(in)    :: BSi_I(nX) ! B-field strength at cell center
    real, intent(inout) :: Hamiltonian2_N(0:nP+1, -1:nMu+1, -1:nX+1)
    real    :: Velocity_I(0:nP+1)    ! Particle velocity array at cell face
    real    :: InvBSi_I(nX), InvBFaceSi_I(0:nX) ! 1/B at center and face
    integer :: iX, iMu               ! Loop variables
    ! Calculate the 2nd hamiltonian function for the 3-Poisson-bracket scheme
    !--------------------------------------------------------------------------

    Velocity_I = 1.0/sqrt(1.0 + (cProtonMass*cLightSpeed/Momentum_I)**2)
    ! Considering the law of relativity, v = 1/sqrt(1+(m*c/p)**2), v can be
    ! calculated as a function of p. Note that light speed is the unit of
    ! speed here, so we do not need to multiply c^2 in the following steps
    ! Note that momentum (or velocity) in this term is non-related to the
    ! Lagrangian coordinates {s_L, mu}, so it should be cell-centered.

    ! Calculate 1/B on the boundary of grid
    InvBSi_I             = 1.0/BSi_I
    InvBFaceSi_I(1:nX-1) = (InvBSi_I(2:nX) + InvBSi_I(1:nX-1))*0.5
    InvBFaceSi_I(0 )     = InvBSi_I(1 ) - (InvBSi_I(2 ) - InvBSi_I(1   ))*0.5
    InvBFaceSi_I(nX)     = InvBSi_I(nX) + (InvBSi_I(nX) - InvBSi_I(nX-1))*0.5

    ! Calculate the second hamiltonian function = (mu^2-1)*v/(2B)
    do iX = 1, nX
       do iMu = 0, nMu
          Hamiltonian2_N(:, iMu, iX) = 0.5*(MuFace_I(iMu)**2 - 1.0)* &
               Velocity_I*InvBFaceSi_I(iX)
          ! Here, what we use actually is: \tilde\deltaH
          Hamiltonian2_N(:, iMu, iX) = Hamiltonian2_N(:, iMu, iX)*VolumeP_I
       end do
    end do

    ! Boundary condition of Hamiltonian function
    Hamiltonian2_N(:,     :,   -1) = Hamiltonian2_N(:,     :,  0)
    Hamiltonian2_N(:,     :, nX+1) = Hamiltonian2_N(:,     :, nX)
    Hamiltonian2_N(:,    -1,    :) = Hamiltonian2_N(:,     1,  :)
    Hamiltonian2_N(:, nMu+1,    :) = Hamiltonian2_N(:, nMu-1,  :)

  end subroutine calc_triple_hamiltonian_2
  !============================================================================
  subroutine calc_triple_hamiltonian_3(nX, Hamiltonian3_N)
    ! Calculate the 3rd Hamiltonian function at each fixed time:
    ! (1-mu^2)/2 * [ mu*(p^3/3)*(3\vec{b}\vec{b}:\nabla\vec{u} - div\vec{u})
    ! + ProtonMass*p^2*(\vec{b}*d\vec{u}/dt) ], regarding to p^3/3 and mu,
    ! so the coordinates of p^3/3 and mu are at face, and s_L is at center.
    ! Here we list the variables in the analytical function:
    ! (3\vec{b}\vec{b}:\nabla\vec{u} - div\vec{u}) = d(ln(B*ds^2))/dt
    ! ProtonMass = cRmeProtonGeV, in the unit of GeV/c^2
    ! \vec{b}*d\vec{u}/dt = bDuDt_C

    use ModConst,           ONLY: cProtonMass
    use SP_ModDistribution, ONLY: MuFace_I
    integer, intent(in) :: nX        ! Number of s_L grid
    real, intent(inout) :: Hamiltonian3_N(-1:nP+1, -1:nMu+1, 0:nX+1)
    integer             :: iX, iMu   ! Loop variables
    ! Calculate the 3rd hamiltonian function for the 3-Poisson-bracket scheme
    !--------------------------------------------------------------------------

    do iX = 1, nX
       do iMu = 0, nMu
          Hamiltonian3_N(:, iMu, iX) = 0.5*(1.0 - MuFace_I(iMu)**2)* &
               (MuFace_I(iMu)*Momentum3_I*dLnBDeltaS2Dt_C(iX) +      &
               cProtonMass*(Momentum3_I*3.0)**(2.0/3.0)*bDuDt_C(iX))
       end do
       ! Here, what we use actually is: \tilde\deltaH
       Hamiltonian3_N(:, 0:nMu, iX) = &
            Hamiltonian3_N(:, 0:nMu, iX)*DeltaSOverB_C(iX)
    end do

    ! Boundary condition of Hamiltonian function
    Hamiltonian3_N(:,     :,    0) = Hamiltonian3_N(:,     :,  1)
    Hamiltonian3_N(:,     :, nX+1) = Hamiltonian3_N(:,     :, nX)
    Hamiltonian3_N(:,    -1,    :) = Hamiltonian3_N(:,     1,  :)
    Hamiltonian3_N(:, nMu+1,    :) = Hamiltonian3_N(:, nMu-1,  :)

  end subroutine calc_triple_hamiltonian_3
  !============================================================================
end module SP_ModAdvancePoisson
!==============================================================================
