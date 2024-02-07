!  Copyright (C) 2002 Regents of the University of Michigan,
!  portions used with permission
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module SP_ModAdvanceAdvection

  use ModUtilities, ONLY: CON_stop
  ! Solves advection in the momentum space (=first order Fermi acceleration)
  ! First way - solve advection over log P coordinate.
  ! In space physics applications one often needs to solve the
  ! "log-advection" equations of the kind of:
  !         f_t+A*f_lnp=0             (NC)
  ! or
  !         f_t+A*(pf)_p=0 .          (C)
  ! The CONSERVATION LAW Eq.(C) can be also written as
  !         f_t+A*f_lnp+ A*f=0,
  ! which is the total of the NONCONSERVATIVE Eq.(NC) plus Af
  ! Both Eqs.(NC,C) propagate the initial values of f along the lines
  !         p=p_0*exp(At),
  ! which results in a uniform expansion (A>0) or shrinking (A<0) of
  ! the distribution, keeping its shape unchanged. For Eq.(C) the solution
  ! also scales accordingly, keeping unchanged the integral \int{f dp}
  !
  ! EXAMPLES: Non-conservative formulation: first-order Fermi acceleration,
  !           Alfven-wave turbulence evolution in the solar wind
  !           Conservative formulation: Kolmogorov-type cascade,
  !           wave-particle interactions
  !
  ! In all these applications, an acceleration rate, A, does not depend on the
  ! phase coordinate, p
  use SP_ModDistribution, ONLY: nP, MomentumSi_I, Distribution_IIB, dLogP
  implicit none
  ! Revision history
  ! Prototype: Sokolov&Roussev, FLAMPA code, 2004
  ! Version: Sokolov& Roussev, Jan,2008, SP/FLAMPA/src/ModLogAdvection.f90

  SAVE

  PRIVATE ! Except

  public :: advect_via_log
contains
  !============================================================================
  subroutine advect_via_log(iLine, nX, iShock, DtProgress, Cfl, &
       dLogRho_I, nSi_I, BSi_I, IsNeg)
    ! The subroutine encapsulates the logarithmic advection scheme,
    ! which is non-conservative, and the diffusion

    use SP_ModDiffusion, ONLY: UseDiffusion, diffuse_distribution
    use SP_ModBc,   ONLY: set_momentum_bc, SpectralIndex
    use SP_ModGrid, ONLY: Used_B, nVertex_B
    use SP_ModTurbulence, ONLY: advance_log_advection
    ! INPUTS:
    ! id of line, particle #, and Shock location
    integer, intent(in):: iLine, nX, iShock
    ! input time step
    real,    intent(in):: DtProgress
    ! CFL number
    real,    intent(in):: Cfl
    ! Ratio of densities at upper and lower level
    real,    intent(in):: dLogRho_I(1:nX)
    ! Density and magnetic field at the upper level
    real,    intent(in):: nSI_I(1:nX), BSI_I(1:nX)
    ! OUTPUT:
    ! Check if any distrubution function value is negative
    logical, intent(out):: IsNeg
    ! local variables, declared in this subroutine
    integer  :: iStep, iVertex      ! loop variables
    ! time step is split for nStep intervals,  so short that the CFL for
    ! (explicit) advection operator is less that CFL declared abobe.
    integer  :: nStep
    ! Time step in the STEP Loop, DtProgress/nStep
    real     :: Dt
    ! Lagrangian derivatives
    real     :: FermiFirst_I(1:nX)

    character(len=*), parameter:: NameSub = 'advect_via_log'
    !--------------------------------------------------------------------------
    ! 1st order Fermi acceleration is responsible for advection
    ! in momentum space
    IsNeg = .false.
    ! first order Fermi acceleration for the current line
    FermiFirst_I = DLogRho_I / (3*DLogP)
    ! if(UseTurbulentSpectrum)then
    ! nStep = 1+int(max(DtReduction,                &
    !      maxval(abs(FermiFirst_I(1:nX))))/CFL)
    ! else
    ! How many steps should be done to the CFL criterion is fulfilled
    nStep = 1+int(maxval(abs(FermiFirst_I(2:nX)))/CFL)
    ! end if

    ! Check if the number of time steps is positive:
    if(nStep < 1)then
       ! if(UseTurbulentSpectrum) &
       !     write(*,*) ' DtReduction =', DtReduction
       write(*,*) ' maxval(abs(FermiFirst_I)) =', &
            maxval(abs(FermiFirst_I(2:nX)))
       call CON_stop(NameSub//': nStep <= 0????')
    end if
    Dt = DtProgress / nStep
    FermiFirst_I = FermiFirst_I / nStep
    ! if(UseTurbulentSpectrum) call reduce_advection_rates(nStep)
    STEP:do iStep = 1, nStep
       ! update bc for advection
       call set_momentum_bc(iLine, nX, nSi_I(1:nX), iShock)
       ! advection in the momentum space
       do iVertex = 1, nX
          if(any(Distribution_IIB(0:nP+1,iVertex,iLine) < 0.0)) then
             write(*,*) NameSub, ': Distribution_IIB < 0'
             Used_B(iLine) = .false.
             nVertex_B(iLine) = 0
             ! CYCLE line
             IsNeg = .true.
             EXIT STEP
          end if

          call advance_log_advection(FermiFirst_I(iVertex), &
               1, 1, Distribution_IIB(0:nP+1,iVertex,iLine), .false.)
       end do
       ! compute diffusion along the field line
       ! set the left boundary condition (for diffusion)
       if(UseDiffusion) call diffuse_distribution(iLine, nX,    &
            iShock, Dt, nSi_I, BSi_I, LowerEndSpectrum_I= &
            Distribution_IIB(0, 1, iLine) * &
            (MomentumSi_I(0)/MomentumSi_I(1:nP))**SpectralIndex)
    end do STEP
  end subroutine advect_via_log
  !============================================================================
end module SP_ModAdvanceAdvection
!==============================================================================
