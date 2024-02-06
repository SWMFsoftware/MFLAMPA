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

  public :: advect_via_log, advance_log_advection
contains
  !============================================================================
  subroutine advect_via_log(iLine, nX, iShock, DtProgress, Cfl, &
       dLogRho_I, nSi_I, BSi_I, IsNeg)
    ! The subroutine encapsulates the logarithmic advection scheme,
    ! which is non-conservative, and the diffusion

    use SP_ModDiffusion, ONLY: UseDiffusion, diffuse_distribution
    use SP_ModBc,   ONLY: set_momentum_bc, SpectralIndex
    use SP_ModGrid, ONLY: Used_B, nVertex_B
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
  ! The procedure integrates the log-advection equation, in
  ! the conservative or non-conservative formulation, at a logarithmic grid,
  ! using a single-stage second order scheme
  subroutine advance_log_advection(CFLIn, nGCLeft, nGCRight,        &
       FInOut_I, IsConservative, DeltaLnP)

    real,   intent(in):: CFLIn        ! Time step * acceleration rate/(Dlnp)
    integer,intent(in):: nGCLeft      ! The solution in the ghost cells is not
    integer,intent(in):: nGCRight     ! advanced in time but used as the BCs
    real,intent(inout):: FInOut_I(1-nGCLeft:nP+nGCRight)  ! Solution
    ! sol. index 0 is to set boundary condition at the injection energy
    logical,intent(in):: IsConservative  ! Solve (C) if .true. (NC) otherwise
    real,optional,intent(in)::DeltaLnP   ! Used only for NonConservative=.true.!
    ! Subcycling to achieve the condition CFL<1, if nStep>1
    integer :: nStep
    ! Loop variables
    integer :: iStep
    ! Extended version of the sulution array to implement BCc
    real    :: F_I(1 - max(nGCLeft,2):nP+max(nGCRight,2))
    real    :: FSemiintUp_I(0:nP+1), FSemiintDown_I(0:nP+1)
    real    :: CFL, HalfADtIfNeeded
    !-------------------------NonConservative---------------------------
    ! This is a one-stage second-order scheme for the advection equation:
    !             f_t+A f_{ln p}=0,
    !
    ! Herewith CFL = A*Delta t/Delta(ln p):
    ! f^(n+1)_i-CFL*(f^n_(i-1/2)-f^n_(i+1/2)=f^n_i, where
    ! For CFL>0:
    ! f^n_(i-1/2)=f^n_(i-1)+cHalf*(cOne-CFL)*df_lim^n_(i-1)
    ! For CFL<0:
    ! f^n_(i-1/2)=f^n_(i  )-cHalf*(cOne+CFL)*df_lim^n_(i  )
    !-------------------------   Conservative---------------------------
    ! This is a one-stage second-order scheme for the advection equation:
    !             f_t+A (f p}_p=0,
    ! We now need to use both
    ! A*Delta t = CFLIn*Delta (LnP)
    ! and
    ! CFL = A*Delta t/(2 tanh(Delta (ln p)/2)
    ! f^(n+1)_i = CFL*(f^n_(i-1/2)-f^n_(i+1/2))-
    !        (1/2)(f^n_(i-1/2)+f^n_(i+1/2))A*Delta t + f^n_i,
    ! where
    ! For CFL>0:
    ! f^n_(i-1/2)=f^n_(i-1)*(1-(1/2)*A*Delta t)+cHalf*(cOne-CFL)*df_lim^n_(i-1)
    ! For CFL<0:
    ! f^n_(i-1/2)=f^n_(i  )*(1-(1/2)*A*Delta t)-cHalf*(cOne+CFL)*df_lim^n_(i  )

    !--------------------------------------------------------------------------
    if(IsConservative)then
       HalfADtIfNeeded = 0.50*CFLIn*DeltaLnP
       CFL = HalfADtIfneeded/tanh(0.50 * DeltaLnP)
       HalfADtIfNeeded = -HalfADtIfNeeded
    else
       HalfADtIfNeeded = 0.0; CFL = CFLIn
    end if

    nStep = 1 + int(abs(CFL)); CFL = CFL/real(nStep)
    HalfADtIfNeeded=HalfADtIfNeeded/real(nStep)
    F_I(1-nGCLeft:nP+nGCRight) = FInOut_I(1 - nGCLeft:nP+nGCRight)

    ! Check for positivity
    if(any(F_I(1-nGCLeft:nP+nGCRight)<=0.0))then
       write(*,*)'Before advection F_I <=0'
       write(*,*)IsConservative
       write(*,*)F_I
       stop
    end if

    ! One stage second order upwind scheme

    if (CFL>0.0) then
       do iStep = 1, nStep
          ! Boundary condition at the left boundary
          if(nGCLeft<2)F_I(            -1:0-nGCLeft) = F_I( 1-nGCLeft )
          ! Boundary condition at the right boundary
          if(nGCRight<2)F_I(nP+1-nGCRight:nP+2     ) = F_I(nP+nGCRight)

          ! f_(i+1/2):
          FSemiintUp_I(0:nP) = F_I(0:nP)*(1.0 + HalfADtIfNeeded)&
               + 0.5*(1.0 - CFL)*df_lim_array(0, nP)
          ! f_(i-1/2):
          FSemiintDown_I(1:nP) = FSemiintUp_I(0:nP-1)

          ! Update the solution from f^(n) to f^(n+1):
          F_I(1:nP) = F_I(1:nP)+CFL*(FSemiintDown_I(1:nP)-FSemiintUp_I(1:nP))+&
               HalfADtIfNeeded*(FSemiintDown_I(1:nP)+FSemiintUp_I(1:nP))
       end do
    else
       do iStep = 1, nStep
          ! Boundary condition at the left boundary
          if(nGCLeft<2)F_I(            -1:0-nGCLeft) = F_I( 1-nGCLeft)
          ! Boundary condition at the right boundary
          if(nGCRight<2)F_I(nP+1-nGCRight:nP+2     ) = F_I(nP+nGCRight)

          ! f_(i-1/2):
          FSemiintDown_I(1:nP+1) = F_I(1:nP+1)*(1.0 + HalfADtIfNeeded)&
               - 0.5*(1.0 + CFL)*df_lim_array(1, nP+1)
          ! f_(i+1/2):
          FSemiintUp_I(1:nP) = FSemiintDown_I(2:nP+1)

          ! Update the solution from f^(n) to f^(n+1):
          F_I(1:nP) = F_I(1:nP)+CFL*(FSemiintDown_I(1:nP)-FSemiintUp_I(1:nP))+&
               HalfADtIfNeeded*(FSemiintDown_I(1:nP)+FSemiintUp_I(1:nP))
       end do
    end if

    FInOut_I(1:nP) = F_I(1:nP)
    if(any(FInOut_I(1-nGCLeft:nP+nGCRight)<=0.0))then
       write(*,*)'After advection F_I <=0, for CFLFermi= ',CFL
       write(*,*)F_I
       call CON_stop('Negative distribution function')
    end if
  contains
    !==========================================================================
    function df_lim_array(iLeft, iRight) result(df_lim_arr)
      integer, intent(in):: iLeft, iRight ! start/left and end/right indices
      real :: df_lim_arr(iLeft:iRight)    ! output results
      real :: dF1(iLeft:iRight), dF2(iLeft:iRight)

      !------------------------------------------------------------------------
      dF1 = F_I(iLeft+1:iRight+1)-F_I(iLeft:iRight)
      dF2 = F_I(iLeft:iRight)-F_I(iLeft-1:iRight-1)

      ! df_lim=0 if dF1*dF2<0, sign(dF1) otherwise:
      df_lim_arr = sign(0.50,dF1)+sign(0.50,dF2)
      dF1 = abs(dF1)
      dF2 = abs(dF2)
      df_lim_arr = df_lim_arr*min(max(dF1,dF2),2.0*dF1,2.0*dF2)
    end function df_lim_array
    !==========================================================================
  end subroutine advance_log_advection
  !============================================================================
end module SP_ModAdvanceAdvection
!==============================================================================
