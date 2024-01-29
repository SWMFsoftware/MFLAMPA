!  Copyright (C) 2002 Regents of the University of Michigan,
!  portions used with permission
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module SP_ModAdvancePoisson
  ! DESCRIPTION:
  !   High resolution finite volume method for kinetic equations 
  !   with Poisson brackets (Sokolov et al., 2023) 
  !   https://doi.org/10.1016/j.jcp.2023.111923
  implicit none
  public :: advect_via_poisson_bracket
contains
  !============================================================================
  subroutine advect_via_poisson_bracket(nP, nX,        &
       tFinal, CflIn, InvRhoOld_I, InvRho_I, FInOut_I, &
       iLine, iShock, XyzSI_DI, nSI_I, BSI_I, DsSI_I,  &
       RadiusSI_I, UseDiffusion)
    ! advect via Possion Bracket + diffusion by encapsulation

    use ModPoissonBracket, ONLY: explicit
    use SP_ModSize, ONLY: nVertexMax
    use SP_ModDistribution, ONLY: Momentum3SI_I, VolumeP_I, DLogP
    use SP_ModDiffusion, ONLY: diffuse_distribution

    integer,intent(in):: nP, nX     ! # of meshes along lnp-coordinate
    real,   intent(in):: tFinal     ! time interval to advance through
    real,   intent(in):: CflIn      ! FermiFirsr_I in ModAdvance
    real,   intent(in):: InvRhoOld_I(1:nX) ! Old density (inverse)
    real,   intent(in):: InvRho_I(1:nX)    ! New density (inverse)
    real,intent(inout):: FInOut_I(0:nP+1, 1:nX)  ! Solution
    ! Variables for diffusion
    integer, intent(in) :: iLine, iShock
    real, intent(in) :: XyzSI_DI(3, 1:nVertexMax)
    real, intent(in), dimension(1:nVertexMax) :: nSI_I,  &
         BSI_I, DsSI_I, RadiusSI_I
    logical, intent(in) :: UseDiffusion   ! diffuse_Distribution or not

    ! Loop variables
    integer :: iP
    ! Extended arrays for the implementation of the Poisson Bracket Alg.
    ! VolumeXOld_I: geometric volume when the subroutine starts
    ! VolumeX_I: geometric volume when the subroutine ends
    ! Linearly interpolate the geometric volume from the start to the end
    ! VolumeSubX_I: current geometric volume at each iteration
    real, dimension(0:nX+1) :: VolumeXOld_I, VolumeX_I,  &
         VolumeSubXOld_I, VolumeSubX_I, dVolumeSubXDt_I
    ! VolumeSub_G: current phase space volume at each iteration
    real, dimension(0:nX+1, 0:nP+1) :: VolumeSubOld_G,   &
         VolumeSub_G, dVolumeSubDt_G
    real    :: dHamiltonian02_FY(0:nX+1, -1:nP+1)
    real    :: VDF_G(-1:nX+2, -1:nP+2), Source_C(nX, nP)
    real    :: Time, Dt, DtTrial, DtInv, DtNext
    !--------------------------------------------------------------------------
    ! Now this is the conservative form for the particle number

    ! Initialize the CFL number, (sub-)nStep and arrays 
    Source_C   = 0.0 
    VDF_G      = 1.0e-8 
    VDF_G(1:nX, 0:nP+1) = transpose(FInOut_I)

    ! Volume initialization: use 1 ghost point at each side of the boundary
    VolumeXOld_I(1:nX)  = InvRhoOld_I
    VolumeXOld_I(0)     = VolumeXOld_I(1)
    VolumeXOld_I(nX+1)  = VolumeXOld_I(nX)
    VolumeX_I(1:nX)     = InvRho_I 
    VolumeX_I(0)        = VolumeX_I(1)
    VolumeX_I(nX+1)     = VolumeX_I(nX)
    VolumeSubX_I        = VolumeXOld_I
    dVolumeSubXDt_I     = (VolumeX_I-VolumeXOld_I)/tFinal

    ! Time initialization
    Time    = 0.0
    DtTrial = CflIn/maxval(abs(dVolumeSubXDt_I)/max(VolumeX_I, VolumeXOld_I))/(3*DLogP)
    DtInv   = 1.0/DtTrial
    DtNext  = DtTrial

    ! Advection by Poisson Bracket Algorithm
    do 
       ! Time Updates
       Dt = min(DtNext, tFinal - Time); DtInv = 1.0/Dt

       ! Volume Updates
       VolumeSubXOld_I = VolumeSubX_I
       VolumeSubX_I = VolumeSubXOld_I + Dt*dVolumeSubXDt_I
       do iP = 0, nP+1 
          VolumeSubOld_G(:,iP) = VolumeP_I(iP)*VolumeSubXOld_I
          VolumeSub_G(:,iP)    = VolumeP_I(iP)*VolumeSubX_I
       end do
       dVolumeSubDt_G = DtInv*(VolumeSub_G - VolumeSubOld_G)

       ! calculate: dHamiltonian/dVolumeSubX
       do iP = -1, nP+1
          dHamiltonian02_FY(:,iP) = - Momentum3SI_I(iP)*dVolumeSubXDt_I
       end do

       call explicit(nX, nP, VDF_G, VolumeSub_G, Source_C,  &
            dHamiltonian02_FY=dHamiltonian02_FY,            &
            dVolumeDt_G = dVolumeSubDt_G,                   &
            DtIn=Dt,           & ! Input time step, which may be reduced
            CFLIn=CflIn,       & ! Input CFL to calculate next time step
            DtOut=DtNext)
       VolumeSubX_I = VolumeSubXOld_I + Dt*dVolumeSubXDt_I

       ! Update velocity distribution function 
       VDF_G(1:nX, 1:nP) = VDF_G(1:nX, 1:nP) + Source_C 
       ! Diffuse the distribution function 
       if(UseDiffusion) then
          FInOut_I = transpose(VDF_G(1:nX, 0:nP+1))
          call diffuse_distribution(iLine, nX, iShock,      &
               Dt, FInOut_I, XyzSI_DI, nSI_I,              &
               BSI_I, DsSI_I, RadiusSI_I)
          VDF_G(1:nX, 0:nP+1) = transpose(FInOut_I)
       end if
       ! Update the time 
       Time = Time + Dt
       if(Time > tFinal - 1.0e-8*DtNext) EXIT

       ! BCs
       VDF_G(1:nX, -1)        = VDF_G(1:nX, 0)
       VDF_G(1:nX,nP+1:nP+2)  = 1.0e-8
       VDF_G( 0,      :)      = VDF_G(1,       :)
       VDF_G(-1,      :)      = VDF_G(1,       :)
       VDF_G(nX+1:nX+2,:)     = 1.0e-8
    end do

    FInOut_I = transpose(VDF_G(1:nX, 0:nP+1))

  end subroutine advect_via_poisson_bracket
  !============================================================================
end module SP_ModAdvancePoisson
