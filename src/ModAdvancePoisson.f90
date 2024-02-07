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
  public :: advect_via_poisson
contains
  !============================================================================
  subroutine advect_via_poisson(iLine, nX, iShock, &
       tFinal, CflIn, nOldSi_I, nSi_I, BSi_I)
    ! advect via Possion Bracket scheme
    ! diffuse the distribution function at each time step

    use ModPoissonBracket, ONLY: explicit
    use SP_ModSize, ONLY: nVertexMax
    use SP_ModDistribution, ONLY: nP, Momentum3Si_I, VolumeP_I, DLogP, &
         Distribution_IIB, MomentumSi_I, MomentumInjSi, Background_I
    use SP_ModDiffusion, ONLY: UseDiffusion, diffuse_distribution
    use SP_ModBc,   ONLY: set_momentum_bc, SpectralIndex
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
       dHamiltonian01_FX(iP,:) = - Momentum3Si_I(iP)*dVolumeXDt_I
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
       ! set the left boundary condition (for diffusion)
       if(UseDiffusion) call diffuse_distribution(iLine, nX, iShock, &
            Dt, nSi_I, BSi_I, LowerEndSpectrum_I=VDF_G(1:nP, 0))
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
      VDF_G(0:nP+1,    0) = Distribution_IIB(0, 1, iLine) * &
           (MomentumSi_I(0)/MomentumSi_I(0:nP+1))**SpectralIndex
      VDF_G(0:nP+1, nX+1) = Background_I
      ! Add a second layer of the ghost cells along the line coordinate:
      VDF_G(0:nP+1,   -1) = VDF_G(0:nP+1,    0)
      VDF_G(0:nP+1, nX+2) = VDF_G(0:nP+1, nX+1)
      ! Add a second layer of the ghost cells along the momentum coordinate:
      VDF_G(-1,:) = VDF_G(0,:); VDF_G(nP+2,:) = VDF_G(nP+1,:)
    end subroutine set_VDF
    !==========================================================================
  end subroutine advect_via_poisson
  !============================================================================
end module SP_ModAdvancePoisson
!==============================================================================
