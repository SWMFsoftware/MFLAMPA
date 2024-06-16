!  Copyright (C) 2002 Regents of the University of Michigan,
!  portions used with permission
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module SP_ModAdvancePoisson

  ! High resolution finite volume method for kinetic equations
  ! with Poisson brackets (Sokolov et al., 2023)
  ! See https://doi.org/10.1016/j.jcp.2023.111923

  use SP_ModSize,         ONLY: nVertexMax
  use SP_ModGrid,         ONLY: nLine
  use SP_ModDistribution, ONLY: nP, nMu, Distribution_CB, &
       Background_I, IsDistNeg, check_dist_neg
  use ModUtilities,       ONLY: CON_stop
  implicit none

  PRIVATE ! Except

  SAVE
  public:: advect_via_poisson       ! Time-accurate advance through given Dt
  public:: iterate_poisson          ! Local time-stepping for steady states
  public:: advect_via_multi_poisson ! Time-accurate advance with pitch angle
  public:: init_data_states         ! Initialize the states for multi_poisson

  ! \Deltas/b, ln(B\deltas^2) at Old time
  real, public, dimension(nVertexMax) :: DeltaSOverBOld_C, LnBDeltaS2Old_C
  ! \Deltas/b, ln(B\deltas^2) at New time
  real, public, dimension(nVertexMax) :: DeltaSOverBNew_C, LnBDeltaS2New_C
  ! Time-derivatives:
  real, public, dimension(nVertexMax) :: DeltaSOverB_C, &
       dDeltaSOverBDt_C, dLnBdeltaS2Dt_C, bDuDt_C
contains
  !============================================================================
  subroutine advect_via_poisson(iLine, nX, iShock, &
       tFinal, CflIn, nOldSi_I, nSi_I, BSi_I)
    ! advect via Possion Bracket scheme
    ! diffuse the distribution function at each time step

    use ModPoissonBracket,  ONLY: explicit
    use SP_ModDistribution, ONLY: dLogP, VolumeP_I, Momentum3_I
    use SP_ModDiffusion,    ONLY: UseDiffusion, diffuse_distribution
    use SP_ModBc,           ONLY: set_momentum_bc, UseUpperEndBc
    ! INPUTS:
    integer, intent(in):: iLine, iShock ! Indices of line and shock
    integer, intent(in):: nX            ! Number of meshes along s_L axis
    real,    intent(in):: tFinal        ! Time interval to advance through
    real,    intent(in):: CflIn         ! Input CFL number
    ! Input variables for diffusion
    real,    intent(in):: nOldSi_I(nX), nSi_I(nX), BSi_I(nX)
    ! Loop variables
    integer :: iX
    ! Extended arrays for implementation of the Poisson bracket scheme
    ! VolumeXStart_I: geometric volume when the subroutine starts
    ! VolumeXEnd_I: geometric volume when the subroutine ends
    ! dVolumeXDt_I: time derivative of geometric volume
    real, dimension(0:nX+1):: VolumeXStart_I, VolumeXEnd_I, dVolumeXDt_I
    ! Volume_G: total control volume at the end of each iteration
    ! VolumeOld_G: total control volume at the end of each iteration
    ! dVolumeDt_G: total control time derivative
    real, dimension(0:nP+1, 0:nX+1):: VolumeOld_G, Volume_G, dVolumeDt_G
    ! DeltaHamiltonian
    real    :: dHamiltonian01_FX(-1:nP+1, 0:nX+1)
    ! Extended array for distribution function
    real    :: VDF_G(-1:nP+2, -1:nX+2)
    ! Advection term
    real    :: Source_C(nP, nX)
    ! Time, ranging from 0 to tFinal
    real    :: Time
    ! Time step
    real    :: Dt
    ! Prediction for the next time step:
    real    :: DtNext
    ! Now this is the particle-number-conservative advection scheme
    character(len=*), parameter:: NameSub = 'advect_via_poisson'
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
    ! Advection by Poisson bracket scheme
    do
       ! Time Updates
       Dt = min(DtNext, tFinal - Time)
       ! Volume Updates
       VolumeOld_G = Volume_G
       Volume_G    = VolumeOld_G + Dt*dVolumeDt_G
       ! BC Updates
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
       Distribution_CB(1:nP, 1, 1:nX, iLine) = &
            Distribution_CB(1:nP, 1, 1:nX, iLine) + Source_C
       ! Check if the VDF includes negative values
       call check_dist_neg(NameSub, 1, nX, iLine)
       if(IsDistNeg)RETURN

       ! Diffuse the distribution function
       if(UseDiffusion) then
          if(UseUpperEndBc) then
             call diffuse_distribution(iLine, nX, iShock, Dt,         &
                  nSi_I, BSi_I, LowerEndSpectrum_I=VDF_G(1:nP, 0),    &
                  UpperEndSpectrum_I=VDF_G(1:nP, nX+1))
          else
             call diffuse_distribution(iLine, nX, iShock, Dt,         &
                  nSi_I, BSi_I, LowerEndSpectrum_I=VDF_G(1:nP, 0))
          end if
          ! Check if the VDF includes negative values after diffusion
          call check_dist_neg(NameSub//' after diffusion', 1, nX, iLine)
          if(IsDistNeg)RETURN
       end if

       ! Update time
       Time = Time + Dt
       if(Time > tFinal - 1.0e-8*DtNext) EXIT
    end do

  end subroutine advect_via_poisson
  !============================================================================
  subroutine iterate_poisson(iLine, nX, iShock, CflIn, BSi_I, nSi_I)
    ! Advect via Possion Bracket scheme to the steady state
    ! Diffuse the distribution function at each time step

    use ModPoissonBracket,  ONLY: explicit
    use SP_ModDistribution, ONLY: VolumeP_I, Momentum3_I
    use SP_ModGrid,         ONLY: State_VIB, D_, U_
    use SP_ModUnit,         ONLY: UnitX_, Io2Si_V
    use SP_ModDiffusion,    ONLY: UseDiffusion, diffuse_distribution
    use SP_ModBc,           ONLY: set_momentum_bc, UseUpperEndBc

    integer, intent(in):: iLine, iShock ! Indices of line and shock
    integer, intent(in):: nX            ! Number of meshes along s_L axis
    real,    intent(in):: CflIn         ! Input CFL number
    ! Input variables for diffusion
    real,    intent(in):: BSi_I(nX), nSi_I(nX)
    real               :: uSi_I(nX-1), DsSi_I(nX-1)
    ! Loop variable
    integer :: iX
    ! Volume_G: global space volume = product of distance in each dimension
    real    :: Volume_G(0:nP+1, 0:nX+1)
    ! VolumeX_I: geometric volume = distance between two geometric faces
    real    :: VolumeX_I(0:nX+1)
    ! u/B variable at face
    real    :: uOverBNodeSi_I(-1:nX+1)
    ! Hamiltonian at cell face
    real    :: Hamiltonian_N(-1:nP+1, -1:nX+1)
    ! Extended array for distribution function
    real    :: VDF_G(-1:nP+2, -1:nX+2)
    ! Advection term
    real    :: Source_C(nP, nX)
    ! Time step
    real    :: Dt_C(nP, nX)
    ! Now particle-number-conservative advection scheme for steady-state soln.
    character(len=*), parameter:: NameSub = 'iterate_poisson'
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
    Distribution_CB(1:nP, 1, 1:nX, iLine) = &
         Distribution_CB(1:nP, 1, 1:nX, iLine) + Source_C
    ! Check if the VDF includes negative values
    call check_dist_neg(NameSub, 1, nX, iLine)
    if(IsDistNeg)RETURN

    ! Diffuse the distribution function
    if(UseDiffusion) then
       if(UseUpperEndBc) then
          call diffuse_distribution(iLine, nX, iShock, Dt_C,     &
               nSi_I, BSi_I, LowerEndSpectrum_I=VDF_G(1:nP, 0),  &
               UpperEndSpectrum_I=VDF_G(1:nP, nX+1))
       else
          call diffuse_distribution(iLine, nX, iShock, Dt_C,     &
               nSi_I, BSi_I, LowerEndSpectrum_I=VDF_G(1:nP, 0))
          ! Check if the VDF includes negative values after diffusion
          call check_dist_neg(NameSub//' after diffusion', 1, nX, iLine)
          if(IsDistNeg)RETURN
       end if
    end if

  end subroutine iterate_poisson
  !============================================================================
  subroutine advect_via_multi_poisson(iLine, nX, iShock,  &
       TimeStart, DtFinal, CflIn, nSi_I, BSi_I)
    ! advect via multiple Possion Bracket scheme for the
    ! focused transport equation considering pitch angle
    ! diffuse the distribution function at each time step

    use ModPoissonBracket,  ONLY: explicit
    use SP_ModDistribution, ONLY: DeltaMu, VolumeP_I
    use SP_ModDiffusion,    ONLY: UseDiffusion, diffuse_distribution
    use SP_ModBc,           ONLY: set_momentum_bc, UseUpperEndBc

    integer, intent(in):: iLine, iShock ! Indices of line and shock
    integer, intent(in):: nX        ! Number of meshes along s_L axis
    real,    intent(in):: TimeStart ! Start time before advancing
    real,    intent(in):: DtFinal   ! Time interval to advance through
    real,    intent(in):: CflIn     ! Input CFL number
    ! Input variables for diffusion
    real,    intent(in):: nSi_I(nX), BSi_I(nX)
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
    real    :: Time, tFinal
    ! Time step
    real    :: Dt
    ! Prediction of next time step:
    real    :: DtNext
    ! Mark whether this is the last run:
    logical :: IsExit
    ! Now this is the particle-number-conservative advection scheme with \mu
    character(len=*), parameter:: NameSub = 'advect_via_multi_poisson'
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
    tFinal = TimeStart + DtFinal
    Time = TimeStart

    ! Here we would like to get the first trial of DtNext
    call set_VDF(iLine, nX, VDF_G)  ! Set the VDF first
    call advance_multi_poisson      ! Now we get DtNext

    ! Advection by multiple-Poisson-bracket scheme
    do
       ! Time Updates
       if(DtNext > tFinal - Time) then
          Dt = tFinal - Time
          IsExit = .true.          ! Last step
       else
          Dt = DtNext
          IsExit = .false.         ! Intermediate steps
       end if

       ! Update bc for at minimal energy, at nP = 0
       call set_momentum_bc(iLine, nX, nSi_I, iShock)
       call set_VDF(iLine, nX, VDF_G)
       ! Advance by multi-Poisson-bracket scheme
       call advance_multi_poisson(DtIn=Dt)

       if(IsExit) then
          ! This step is the last step
          Distribution_CB(1:nP, 1:nMu, 1:nX, iLine) =       &
               Distribution_CB(1:nP, 1:nMu, 1:nX, iLine) + Source_C
       else
          ! This step is not the last step
          ! Update VDF_G considering the time-dependent Volume_G

          Source_C = Source_C*Volume_G(1:nP, 1:nMu, 1:nX)
          ! Update DeltaSOverB_C for the calculation of volume
          call update_states(iLine, nX, Time)
          ! Calculate the total control volume
          do iX = 1, nX
             Volume_G(:, :, iX) = DeltaSOverB_C(iX)*DeltaMu* &
                  spread(VolumeP_I, DIM=2, NCOPIES=nMu+2)
          end do

          ! Update VDF_G to CURRENT time: no BCs for (1:nQ, 1:nP, 1;nR)
          Distribution_CB(1:nP, 1:nMu, 1:nX, iLine) =       &
               Distribution_CB(1:nP, 1:nMu, 1:nX, iLine) +  &
               Source_C/Volume_G(1:nP, 1:nMu, 1:nX)
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
    subroutine advance_multi_poisson(DtIn)
      ! advance by multiple Poisson brackets for each time step

      real, optional :: DtIn
      ! Get DeltaSOverB_C at current time: for calculating
      ! the total control volume and Hamiltonian functions
      !------------------------------------------------------------------------

      call update_states(iLine, nX, Time)
      Hamiltonian2_N = 0.0
      Hamiltonian3_N = 0.0

      ! Calculate 1st Hamiltonian function used in the time-dependent
      ! poisson bracket: {f_jk, (p^3/3)*(DeltaS/B)}_{tau, p^3/3}
      call calc_hamiltonian_1(nX, dHamiltonian01_FX)
      ! Calculate 2nd Hamiltonian function used in the time-dependent
      ! poisson bracket: {f_jk, (mu^2-1)*v/(2B)}_{s_L, mu}
      call calc_hamiltonian_2(nX, BSi_I, Hamiltonian2_N)
      ! Calculate 3rd Hamiltonian function used in the time-dependent
      ! poisson bracket: {f_jk, (1-mu^2)/2 * [ mu*(p^3/3)*
      ! d(ln(B*ds^2))/dt + ProtonMass*p^2*bDuDt_C ]}_{p^3/3, mu}
      call calc_hamiltonian_3(nX, Hamiltonian3_N)

      ! Calculate the total control volume
      do iX = 1, nX
         Volume_G(:, :, iX) = DeltaSOverB_C(iX)*DeltaMu* &
              spread(VolumeP_I, DIM=2, NCOPIES=nMu+2)
      end do

      ! Boundary conditions for total control volume
      Volume_G(:, :,    0) = Volume_G(:, :,  1)
      Volume_G(:, :, nX+1) = Volume_G(:, :, nX)

      call explicit(nX, nMu, nX, VDF_G, Volume_G, Source_C, &
           Hamiltonian12_N = Hamiltonian2_N,                &
           Hamiltonian23_N = -Hamiltonian3_N,               &
           dHamiltonian03_FZ = dHamiltonian01_FX,           &
           dVolumeDt_G = dVolumeDt_G,                       &
           DtIn = DtIn, DtOut = DtNext, CFLIn=CflIn)

    end subroutine advance_multi_poisson
    !==========================================================================
  end subroutine advect_via_multi_poisson
  !============================================================================
  subroutine set_VDF(iLine, nX, VDF_G)
    ! We need the VDF on the extended grid with two layers of ghost cells,
    ! to solve the second order scheme. Add solution in physical cells and
    ! in a single layer of the ghost cells along the momentum coordinate.

    use SP_ModBc, ONLY: UseUpperEndBc, set_upper_end_bc, &
         UpperEndBc_I, set_lower_end_bc, LowerEndBc_I
    integer, intent(in):: iLine     ! Indices of line and shock
    integer, intent(in):: nX        ! Number of meshes along s_L axis
    real, intent(inout):: VDF_G(-1:nP+2, -1:nX+2)

    !--------------------------------------------------------------------------
    VDF_G(0:nP+1, 1:nX) = Distribution_CB(:, 1, 1:nX, iLine)

    ! Manipulate the LowerEndBc along the line coordinate:
    call set_lower_end_bc(iLine)
    VDF_G(0:nP+1,    0) = max(LowerEndBc_I, Background_I)

    ! Manipulate the UpperEndBc along the line coordinate:
    if(UseUpperEndBc) then
       call set_upper_end_bc(iLine, nX)
       VDF_G(1:nP, nX+1) = max(UpperEndBc_I, Background_I(1:nP))
       VDF_G(0   , nX+1) = max(VDF_G(0, nX), Background_I(0))
       VDF_G(nP+1, nX+1) = max(VDF_G(nP+1,nX), Background_I(nP+1))
    else
       VDF_G(0:nP+1, nX+1) = Background_I
    end if

    ! Add a second layer of the ghost cells along the line coordinate:
    VDF_G(0:nP+1,   -1) = VDF_G(0:nP+1,    0)
    VDF_G(0:nP+1, nX+2) = VDF_G(0:nP+1, nX+1)
    ! Add a second layer of the ghost cells along the momentum coordinate:
    VDF_G(-1  , :) = VDF_G(0   , :)
    VDF_G(nP+2, :) = VDF_G(nP+1, :)
  end subroutine set_VDF
  !============================================================================
  subroutine init_data_states(iLine, nX, DtFull)
    ! Calculate data states from the input files

    use SP_ModGrid, ONLY: State_VIB, MHData_VIB, x_, z_, BOld_, B_, UOld_, U_
    ! Line number and number along grid axis
    integer, intent(in) :: iLine, nX
    ! Time difference between Old (State_VIB) and New (MHData_VIB) States
    real, intent(in)    :: DtFull
    ! Midpoint for to consecutive points, \deltas
    real                :: MidPoint_IB(x_:z_, nX), DeltaS_I(nX)
    real                :: InvDtFull
    !--------------------------------------------------------------------------
    InvDtFull = 1.0/DtFull

    ! Calculate values at OLD time
    ! Calculate midpoints
    MidPoint_IB(x_:z_, 1:nX-1) = (State_VIB(x_:z_, 2:nX, iLine)   &
         + State_VIB(x_:z_, 1:nX-1, iLine))*0.5
    ! Calculate DeltaS
    DeltaS_I(2:nX-1) = sqrt(sum((MidPoint_IB(x_:z_, 2:nX-1) &
         - MidPoint_IB(x_:z_, 1:nX-2))**2, dim=1))
    ! Linear interpolate the deltas such that there will be nQ deltas
    DeltaS_I(1 ) = 2*sqrt(sum((MidPoint_IB(x_:z_, 1)     &
         - State_VIB(x_:z_,  1, iLine))**2))
    DeltaS_I(nX) = 2*sqrt(sum((MidPoint_IB(x_:z_, nX-1)  &
         - State_VIB(x_:z_, nX, iLine))**2))
    ! Calculate \DeltaS/B at grid center
    DeltaSOverBOld_C(1:nX) = DeltaS_I/State_VIB(BOld_, 1:nX, iLine)
    ! Calculate ln(B*\DeltaS^2) at grid center
    LnBDeltaS2Old_C(1:nX) = log(State_VIB(BOld_, 1:nX, iLine)*DeltaS_I**2)

    ! Calculate values at NEW time
    MidPoint_IB(x_:z_, 1:nX-1) = (MHData_VIB(x_:z_, 2:nX, iLine)  &
         + MhData_VIB(x_:z_, 1:nX-1, iLine))*0.5
    ! Calculate DeltaS
    DeltaS_I(2:nX-1) = sqrt(sum((MidPoint_IB(x_:z_, 2:nX-1) &
         - MidPoint_IB(x_:z_, 1:nX-2))**2, dim=1))
    ! Linear interpolate the deltas such that there will be nQ deltas
    DeltaS_I(1 ) = 2*sqrt(sum((MidPoint_IB(x_:z_, 1)     &
         - MHData_VIB(x_:z_,  1, iLine))**2))
    DeltaS_I(nX) = 2*sqrt(sum((MidPoint_IB(x_:z_, nX-1)  &
         - MHData_VIB(x_:z_, nX, iLine))**2))
    ! Calculate \DeltaS/B at grid center
    DeltaSOverBNew_C(1:nX) = DeltaS_I/State_VIB(B_, 1:nX, iLine)
    ! Calculate ln(B*\DeltaS^2) at grid center
    LnBDeltaS2New_C(1:nX) = log(State_VIB(B_, 1:nX, iLine)*DeltaS_I**2)

    ! Calculate the time-derivative physical quantities
    ! Calculate \deltas/B time derivative = D[delta(s_L)/B]/Dt
    dDeltaSOverBDt_C(1:nX) = (DeltaSOverBNew_C(1:nX) -   &
         DeltaSOverBOld_C(1:nX))*InvDtFull
    ! Calculate Dln(B\deltas^2)/Dt
    dLnBdeltaS2Dt_C(1:nX)  = (LnBDeltaS2New_C(1:nX) -    &
         LnBDeltaS2Old_C(1:nX))*InvDtFull
    ! Calculate b*Du/Dt
    bDuDt_C(1:nX)          = (State_VIB(U_, 1:nX, iLine) -  &
         State_VIB(UOld_, 1:nX, iLine))*InvDtFull

  end subroutine init_data_states
  !============================================================================
  subroutine update_states(iLine, nX, Time)
    ! Update states according to the current time

    integer, intent(in) :: iLine  ! Index of the current field line
    integer, intent(in) :: nX     ! Number of grid along s_L axis
    real, intent(in)    :: Time   ! Current time
    ! Calculate values for CURRENT time: here we use linear interpolation
    ! to get the data of each time step from every to consecutive files
    !--------------------------------------------------------------------------

    ! Update \deltas/B
    DeltaSOverB_C(1:nX) = DeltaSOverBOld_C(1:nX) + dDeltaSOverBDt_C(1:nX)*Time
  end subroutine update_states
  !============================================================================
  subroutine calc_hamiltonian_1(nX, dHamiltonian01_FX)
    ! Calculate the 1st Hamiltonian function with time:
    ! p^3/3*\deltas/B at cell face of p^3/3, regarding to tau and p^3/3

    use SP_ModDistribution, ONLY: DeltaMu, Momentum3_I
    integer, intent(in) :: nX                   ! Number of s_L grid
    real, intent(inout) :: dHamiltonian01_FX(-1:nP+1, 0:nMu+1, 0:nX+1)
    integer             :: iX, iMu              ! Loop variables
    ! Calculate the first Hamiltonian function
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
  end subroutine calc_hamiltonian_1
  !============================================================================
  subroutine calc_hamiltonian_2(nX, BSi_I, Hamiltonian2_N)
    ! Calculate the 2nd Hamiltonian function at each fixed time:
    ! (mu^2-1)*v/(2B) at face of s_L and mu, regarding to s_L and mu

    use ModConst, ONLY: cProtonMass, cLightSpeed
    use SP_ModDistribution, ONLY: Momentum_I, VolumeP_I, MuFace_I
    integer, intent(in) :: nX             ! Number of s_L grid
    real, intent(in)    :: BSi_I(nX)      ! B-field strength at cell center
    real, intent(inout) :: Hamiltonian2_N(0:nP+1, -1:nMu+1, -1:nX+1)
    real    :: Velocity_I(0:nP+1)   ! Particle velocity array at cell face
    real    :: InvBSi_I(nX), InvBFaceSi_I(0:nX) ! 1/B at center and face
    integer :: iX, iMu              ! Loop variables

    !--------------------------------------------------------------------------
    Velocity_I = 1.0/sqrt(1.0 + (cProtonMass*cLightSpeed/Momentum_I)**2)
    ! Considering the law of relativity, v=1/sqrt(1+m^2*c^2/p^2), v can be
    ! calculated as a function of p. Note that light speed is the unit of
    ! speed here, so we do not need to multiply c^2 in the following steps

    ! Calculate 1/B on the boundary of grid
    InvBSi_I             = 1.0/BSi_I
    InvBFaceSi_I(1:nX-1) = (InvBSi_I(2:nX) + InvBSi_I(1:nX-1))*0.5
    InvBFaceSi_I(0 )     = InvBSi_I(1 ) - (InvBSi_I(2 ) - InvBSi_I(1   ))*0.5
    InvBFaceSi_I(nX)     = InvBSi_I(nX) + (InvBSi_I(nX) - InvBSi_I(nX-1))*0.5

    ! Calculate the second hamiltonian function = (mu^2-1)*v/(2B)
    do iX = 1, nX
       do iMu = 0, nMu
          Hamiltonian2_N(:, iMu, iX) = (MuFace_I(iMu)**2 - 1.0)*  &
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
  end subroutine calc_hamiltonian_2
  !============================================================================
  subroutine calc_hamiltonian_3(nX, Hamiltonian3_N)
    ! Calculate the 3rd Hamiltonian function at each fixed time:
    ! (1-mu^2)/2 * [ mu*(p^3/3)*(3\vec{b}\vec{b}:\nabla\vec{u} - div\vec{u})
    ! + ProtonMass*p^2*(\vec{b}*d\vec{u}/dt) ], regarding to p^3/3 and mu,
    ! so the coordinates of p^3/3 and mu are at face, and s_L is at center.
    ! Here we list the variables in the analytical function:
    ! (3\vec{b}\vec{b}:\nabla\vec{u} - div\vec{u}) = d(ln(B*ds^2))/dt
    ! ProtonMass = cRmeProtonGeV, in the unit of GeV/c^2
    ! \vec{b}*d\vec{u}/dt = bDuDt_C

    use ModConst, ONLY: cProtonMass
    use SP_ModDistribution, ONLY: Momentum3_I, MuFace_I
    integer, intent(in) :: nX             ! Number of s_L grid
    real, intent(inout) :: Hamiltonian3_N(-1:nP+1, -1:nMu+1, 0:nX+1)
    integer             :: iX, iMu        ! Loop variables
    ! Calculate the third hamiltonian function

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
  end subroutine calc_hamiltonian_3
  !============================================================================
end module SP_ModAdvancePoisson
!==============================================================================
