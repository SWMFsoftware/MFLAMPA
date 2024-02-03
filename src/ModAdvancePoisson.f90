!  Copyright (C) 2002 Regents of the University of Michigan,
!  portions used with permission
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module SP_ModAdvancePoisson
  !   High resolution finite volume method for kinetic equations
  !   with Poisson brackets (Sokolov et al., 2023)
  !   https://doi.org/10.1016/j.jcp.2023.111923
  implicit none
  public :: advect_via_poisson_bracket
contains
  !============================================================================
  subroutine advect_via_poisson_bracket(nX, tFinal, CflIn,  &
       iLine, iShock, XyzSI_DI, nOldSI_I,      &
       nSI_I, BSI_I, DsSI_I, RadiusSI_I, UseDiffusion)
    ! advect via Possion Bracket + diffusion by encapsulation

    use ModPoissonBracket, ONLY: explicit
    use SP_ModSize, ONLY: nVertexMax
    use SP_ModDistribution, ONLY: nP, Momentum3SI_I,        &
         VolumeP_I, DLogP, Distribution_IIB
    use SP_ModDiffusion, ONLY: diffuse_distribution

    integer,intent(in):: nX         ! # of meshes along lnp-coordinate
    real,   intent(in):: tFinal     ! time interval to advance through
    real,   intent(in):: CflIn      !
    ! Variables for diffusion
    integer, intent(in) :: iLine, iShock
    real, intent(in) :: XyzSI_DI(3, 1:nVertexMax)
    real, intent(in), dimension(1:nVertexMax) :: nOldSI_I, nSI_I,  &
         BSI_I, DsSI_I, RadiusSI_I
    logical, intent(in) :: UseDiffusion   ! diffuse_Distribution or not
    ! Loop variables
    integer :: iP
    ! Extended arrays for the implementation of the Poisson Bracket Alg.
    ! VolumeXStart_I: geometric volume when the subroutine starts
    ! VolumeXEnd_I: geometric volume when the subroutine ends
    ! dVolumeXDt_I: time derivative of geometric volume
    real, dimension(0:nX+1) :: VolumeXStart_I, VolumeXEnd_I, dVolumeXDt_I
    ! Volume_G: phase space volume at the end of each iteration
    ! VolumeOld_G : phase space volume at the end of each iteration
    ! dVolumeDt_G : phase space time derivative
    real, dimension(0:nP+1, 0:nX+1) :: VolumeOld_G, Volume_G, dVolumeDt_G
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
    ! Now this is the conservative form for the particle number

    ! Initialize arrays
    !--------------------------------------------------------------------------
    Source_C   = 0.0
    VDF_G      = 1.0e-8
    VDF_G(0:nP+1, 1:nX) = Distribution_IIB(:, 1:nX, iLine)

    ! Geometric volume: use 1 ghost point at each side of the boundary
    ! Start volume
    VolumeXStart_I(1:nX)  = 1/nOldSI_I(1:nX)
    VolumeXStart_I(0)     = VolumeXStart_I(1)
    VolumeXStart_I(nX+1)  = VolumeXStart_I(nX)
    ! End volume
    VolumeXEnd_I(1:nX)    = 1/nSI_I(1:nX)
    VolumeXEnd_I(0)       = VolumeXEnd_I(1)
    VolumeXEnd_I(nX+1)    = VolumeXEnd_I(nX)
    ! Time derivative
    dVolumeXDt_I       = (VolumeXEnd_I - VolumeXStart_I)/tFinal
    ! Phase volume: initial and time derivative
    do iP = 0, nP+1
       Volume_G(iP,:)    = VolumeP_I(iP)*VolumeXStart_I
       dVolumeDt_G(iP,:) = VolumeP_I(iP)*dVolumeXDt_I
    end do
    ! calculate: dHamiltonian/dVolumeSubX
    do iP = -1, nP+1
       dHamiltonian01_FX(iP,:) = - Momentum3SI_I(iP)*dVolumeXDt_I
    end do
    ! Time initialization
    Time    = 0.0
    ! Trial timestep
    DtNext  = CflIn/maxval(abs(dVolumeXDt_I)/&
         max(VolumeXEnd_I, VolumeXStart_I))/(3*DLogP)
    ! Advection by Poisson Bracket Algorithm
    do
       ! Time Updates
       Dt = min(DtNext, tFinal - Time)
       ! Volume Updates
       VolumeOld_G = Volume_G
       Volume_G = VolumeOld_G + Dt*dVolumeDt_G

       call explicit(nP, nX, VDF_G, Volume_G, Source_C,  &
            dHamiltonian01_FX=dHamiltonian01_FX,         &
            dVolumeDt_G = dVolumeDt_G,                   &
            DtIn=Dt,           & ! Input time step, which may be reduced
            CFLIn=CflIn,       & ! Input CFL to calculate next time step
            DtOut=DtNext)
       ! May need to correct the volume if the time step has been reduced
       Volume_G = VolumeOld_G + Dt*dVolumeDt_G

       ! Update velocity distribution function
       VDF_G(1:nP, 1:nX) = VDF_G(1:nP, 1:nX) + Source_C
       if(UseDiffusion)call diffuse_distribution(iLine, nX, iShock,   &
            Dt, VDF_G(0:nP+1, 1:nX), XyzSI_DI, nSI_I, &
            BSI_I, DsSI_I, RadiusSI_I)
       ! Update time
       Time = Time + Dt
       if(Time > tFinal - 1.0e-8*DtNext) EXIT

       ! BCs
       VDF_G(-1, 1:nX)        = VDF_G(0, 1:nX)
       VDF_G(nP+1:nP+2, 1:nX) = 1.0e-8
       VDF_G(:, 0   )         = VDF_G(:, 1   )
       VDF_G(:, -1  )         = VDF_G(:, 1   )
       VDF_G(:, nX+1:nX+2)    = 1.0e-8
    end do

    Distribution_IIB(:, 1:nX, iLine) = VDF_G(0:nP+1, 1:nX)

  end subroutine advect_via_poisson_bracket
  !============================================================================
end module SP_ModAdvancePoisson
!==============================================================================
