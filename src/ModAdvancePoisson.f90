!  Copyright (C) 2002 Regents of the University of Michigan,
!  portions used with permission
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module SP_ModAdvancePoisson

  ! High resolution finite volume method for kinetic equations
  ! with Poisson brackets (Sokolov et al., 2023)
  ! See https://doi.org/10.1016/j.jcp.2023.111923
  implicit none

  PRIVATE ! Except

  SAVE
  public :: advect_via_poisson ! Time-accurate advance through given Dt
  public :: iterate_poisson    ! Local time-stepping
contains
  !============================================================================
  subroutine advect_via_poisson(iLine, nX, iShock, &
       tFinal, CflIn, nOldSi_I, nSi_I, BSi_I)
    ! advect via Possion Bracket scheme
    ! diffuse the distribution function at each time step

    use ModPoissonBracket, ONLY: explicit
    use SP_ModDistribution, ONLY: nP, Momentum3_I, VolumeP_I, DLogP, &
         Distribution_IIB, Momentum_I, Background_I
    use SP_ModDiffusion, ONLY: UseDiffusion, diffuse_distribution
    use SP_ModBc,   ONLY: set_momentum_bc, SpectralIndex, &
         UseUpperEndBc, set_upper_end_bc, UpperEndBc_I
    integer, intent(in):: iLine, iShock ! indices of line and shock
    integer, intent(in):: nX        ! # of meshes along lnp-coordinate
    real,    intent(in):: tFinal    ! time interval to advance through
    real,    intent(in):: CflIn     ! input CFL number
    ! Input variables for diffusion
    real,    intent(in):: nOldSi_I(nX), nSi_I(nX), BSi_I(nX)
    ! Loop variables
    integer :: iP
    ! Extended arrays for implementation of the Poisson Bracket Alg.
    ! VolumeXStart_I: geometric volume when the subroutine starts
    ! VolumeXEnd_I: geometric volume when the subroutine ends
    ! dVolumeXDt_I: time derivative of geometric volume
    real, dimension(0:nX+1):: VolumeXStart_I, VolumeXEnd_I, dVolumeXDt_I
    ! Volume_G: phase space volume at the end of each iteration
    ! VolumeOld_G : phase space volume at the end of each iteration
    ! dVolumeDt_G : phase space time derivative
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
    !--------------------------------------------------------------------------
    ! Initialize arrays
    ! Geometric volume: use 1 ghost point at each side of the boundary
    ! Start volume
    VolumeXStart_I(1:nX) = 1/nOldSi_I(1:nX)
    VolumeXStart_I(0)    = VolumeXStart_I(1)
    VolumeXStart_I(nX+1) = VolumeXStart_I(nX)
    ! End volume
    VolumeXEnd_I(1:nX)   = 1/nSi_I(1:nX)
    VolumeXEnd_I(0)      = VolumeXEnd_I(1)
    VolumeXEnd_I(nX+1)   = VolumeXEnd_I(nX)
    ! Time derivative
    dVolumeXDt_I         = (VolumeXEnd_I - VolumeXStart_I)/tFinal
    ! Phase volume: initial and time derivative
    do iP = 0, nP+1
       Volume_G(iP,:)    = VolumeP_I(iP)*VolumeXStart_I
       dVolumeDt_G(iP,:) = VolumeP_I(iP)*dVolumeXDt_I
    end do
    ! calculate: dHamiltonian/dVolumeSubX
    do iP = -1, nP+1
       dHamiltonian01_FX(iP,:) = - Momentum3_I(iP)*dVolumeXDt_I
    end do
    ! Time initialization
    Time   = 0.0
    ! Trial timestep
    DtNext = CflIn/maxval(abs(dVolumeXDt_I)/ &
         max(VolumeXEnd_I, VolumeXStart_I))/(3*DLogP)

    ! Advection by Poisson Bracket Algorithm
    do
       ! Time Updates
       Dt = min(DtNext, tFinal - Time)
       ! Volume Updates
       VolumeOld_G = Volume_G
       Volume_G    = VolumeOld_G + Dt*dVolumeDt_G
       ! update bc for at minimal and maximal energy
       call set_momentum_bc(iLine, nX, nSi_I, iShock)
       call set_VDF
       call explicit(nP, nX, VDF_G, Volume_G, Source_C,  &
            dHamiltonian01_FX=dHamiltonian01_FX,         &
            dVolumeDt_G = dVolumeDt_G,                   &
            DtIn=Dt,           & ! Input time step, which may be reduced
            CFLIn=CflIn,       & ! Input CFL to calculate next time step
            DtOut=DtNext)
       ! May need to correct the volume if the time step has been reduced
       Volume_G = VolumeOld_G + Dt*dVolumeDt_G
       ! Update velocity distribution function
       Distribution_IIB(1:nP, 1:nX, iLine) = &
            Distribution_IIB(1:nP, 1:nX, iLine) + Source_C

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
       end if
       ! Update time
       Time = Time + Dt
       if(Time > tFinal - 1.0e-8*DtNext) EXIT
    end do
  contains
    !==========================================================================
    subroutine set_VDF
      ! We need the VDF on the extended grid with two layers of ghost cells,
      ! to solve the second order scheme. Add solution in physical cells and
      ! in a single layer of the ghost cells along the momentum coordinate:
      !------------------------------------------------------------------------
      VDF_G(0:nP+1, 1:nX) = Distribution_IIB(:, 1:nX, iLine)
      ! Apply bc along the line coordinate:
      VDF_G(0:nP+1,    0) = Distribution_IIB(0, 1, iLine)*  &
           (Momentum_I(0)/Momentum_I(0:nP+1))**SpectralIndex
      if(UseUpperEndBc) then
         call set_upper_end_bc(iLine, nX)
         VDF_G(1:nP, nX+1) = UpperEndBc_I
         VDF_G(0, nX+1) = VDF_G(0, nX); VDF_G(nP+1, nX+1) = VDF_G(nP+1, nX)
      else
         VDF_G(0:nP+1, nX+1) = Background_I
      end if
      ! Add a second layer of the ghost cells along the line coordinate:
      VDF_G(0:nP+1,   -1) = VDF_G(0:nP+1,    0)
      VDF_G(0:nP+1, nX+2) = VDF_G(0:nP+1, nX+1)
      ! Add a second layer of the ghost cells along the momentum coordinate:
      VDF_G(-1,:) = VDF_G(0,:); VDF_G(nP+2,:) = VDF_G(nP+1,:)
    end subroutine set_VDF
    !==========================================================================
  end subroutine advect_via_poisson
  !============================================================================
  subroutine iterate_poisson(iLine, nX, iShock, CflIn, &
       uSi_I, BSi_I, nSi_I, DsSi_I)
    ! Advect via Possion Bracket scheme to the steady state
    ! Diffuse the distribution function at each time step

    use ModNumConst,        ONLY: cTiny
    use ModPoissonBracket,  ONLY: explicit
    use SP_ModDistribution, ONLY: nP, Momentum3_I, VolumeP_I,         &
         Distribution_IIB, Momentum_I, Background_I
    use SP_ModDiffusion,    ONLY: UseDiffusion, diffuse_distribution
    use SP_ModBc,           ONLY: set_momentum_bc, SpectralIndex,     &
         UseUpperEndBc, set_upper_end_bc, UpperEndBc_I
    integer, intent(in):: iLine, iShock ! indices of line and shock
    integer, intent(in):: nX        ! # of meshes along lnP-coordinate
    real,    intent(in):: CflIn     ! input CFL number
    ! Input variables for diffusion
    real,    intent(in):: uSi_I(nX), BSi_I(nX), nSi_I(nX), DsSi_I(nX)
    ! Loop variable
    integer :: iP
    ! Volume_G: phase space volume
    real    :: Volume_G(0:nP+1, 0:nX+1)
    ! VolumeX_I: geometric volume
    real    :: VolumeX_I(0:nX+1)
    ! Node-centered u and B variables
    real    :: uNodeSi_I(-1:nX+1), BNodeSi_I(-1:nX+1)
    ! Hamiltonian
    real    :: Hamiltonian_N(-1:nP+1, -1:nX+1)
    ! Extended array for distribution function
    real    :: VDF_G(-1:nP+2, -1:nX+2)
    ! Advection term
    real    :: Source_C(nP, nX)
    ! Time step
    real    :: Dt_C(nP, nX)
    ! Now particle-number-conservative advection scheme for steady-state soln.
    !--------------------------------------------------------------------------
    ! Initialize arrays
    VolumeX_I(2:nX-1) = max(0.5*(DsSi_I(2:nX-1) +      &
         DsSi_I(1:nX-2)), cTiny)/BSi_I(2:nX-1)
    VolumeX_I(1)      = DsSi_I(1)/BSi_I(1)
    VolumeX_I(0)      = VolumeX_I(1)
    VolumeX_I(nX)     = DsSi_I(nX-1)/BSi_I(nX)
    VolumeX_I(nX+1)   = VolumeX_I(nX)
    ! Phase volume: initial values
    do iP = 0, nP+1
       Volume_G(iP,:) = VolumeP_I(iP)*VolumeX_I(0:nX+1)
    end do
    ! Calculate node-centered u=|vector{u}*vector{B}|/|vector{B}|
    uNodeSi_I(2:nX-1) = 0.5*(uSi_I(2:nX-1) + uSi_I(1:nX-2))
    uNodeSi_I(1)      = uSi_I(1)
    uNodeSi_I(-1:0)   = uNodeSi_I(1)
    uNodeSi_I(nX)     = uSi_I(nX)
    uNodeSi_I(nX+1)   = uNodeSi_I(nX)
    ! Calculate node-centered B=|vector{B}|
    BNodeSi_I(2:nX-1) = 0.5*(BSi_I(2:nX-1) + BSi_I(1:nX-2))
    BNodeSi_I(1)      = BSi_I(1)
    BNodeSi_I(-1:0)   = BNodeSi_I(1)
    BNodeSi_I(nX)     = BSi_I(nX)
    BNodeSi_I(nX+1)   = BNodeSi_I(nX)
    ! Calculate Hamiltonian = (u/B)*(p**3/3), node-centered
    do iP = -1, nP+1
       Hamiltonian_N(iP, :) = -uNodeSi_I/BNodeSi_I*Momentum3_I(iP)
    end do

    ! Update bc for at minimal and maximal energy (left BC)
    call set_momentum_bc(iLine, nX, nSi_I, iShock)
    call set_VDF
    call explicit(nP, nX, VDF_G, Volume_G, Source_C,  &
         Hamiltonian12_N=Hamiltonian_N, CFLIn=CflIn,  &
         IsSteadyState=.true., DtOut_C=Dt_C)

    ! Update velocity distribution function
    Distribution_IIB(1:nP, 1:nX, iLine) = &
         Distribution_IIB(1:nP, 1:nX, iLine) + Source_C

    ! Diffuse the distribution function
    if(UseDiffusion) then
       if(UseUpperEndBc) then
          call diffuse_distribution(iLine, nX, iShock, Dt_C,     &
               nSi_I, BSi_I, LowerEndSpectrum_I=VDF_G(1:nP, 0),  &
               UpperEndSpectrum_I=VDF_G(1:nP, nX+1))
       else
          call diffuse_distribution(iLine, nX, iShock, Dt_C,     &
               nSi_I, BSi_I, LowerEndSpectrum_I=VDF_G(1:nP, 0))
       end if
    end if
  contains
    !==========================================================================
    subroutine set_VDF
      ! We need the VDF on the extended grid with two layers of ghost cells,
      ! to solve the second order scheme. Add solution in physical cells and
      ! in a single layer of the ghost cells along the momentum coordinate:
      !------------------------------------------------------------------------
      VDF_G(0:nP+1, 1:nX) = Distribution_IIB(:, 1:nX, iLine)
      ! Apply bc along the line coordinate:
      VDF_G(0:nP+1,    0) = Distribution_IIB(0, 1, iLine)*  &
           (Momentum_I(0)/Momentum_I(0:nP+1))**SpectralIndex
      if(UseUpperEndBc) then
         call set_upper_end_bc(iLine, nX)
         VDF_G(1:nP, nX+1) = UpperEndBc_I
         VDF_G(0, nX+1) = VDF_G(0, nX); VDF_G(nP+1, nX+1) = VDF_G(nP+1, nX)
      else
         VDF_G(0:nP+1, nX+1) = Background_I
      end if
      ! Add a second layer of the ghost cells along the line coordinate:
      VDF_G(0:nP+1,   -1) = VDF_G(0:nP+1,    0)
      VDF_G(0:nP+1, nX+2) = VDF_G(0:nP+1, nX+1)
      ! Add a second layer of the ghost cells along the momentum coordinate:
      VDF_G(-1,:) = VDF_G(0,:); VDF_G(nP+2,:) = VDF_G(nP+1,:)
    end subroutine set_VDF
    !==========================================================================
  end subroutine iterate_poisson
  !============================================================================
end module SP_ModAdvancePoisson
!==============================================================================
