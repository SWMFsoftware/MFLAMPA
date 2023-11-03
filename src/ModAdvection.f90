!  Copyright (C) 2002 Regents of the University of Michigan,
!  portions used with permission
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module SP_ModAdvection
  !DESCRIPTION:
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
  ! In all these applications, A does not depend on the phase coordinate, p
  implicit none
  ! Revision history
  ! Prototype: Sokolov&Roussev, FLAMPA code, 2004
  ! Version: Sokolov& Roussev, Jan,2008, SP/FLAMPA/src/ModLogAdvection.f90
contains
  !============================================================================
  !===========advance_log_advection=======================================
  ! DESCRIPTION: the procedure integrates the log-advection equation, in
  ! the conservative or non-conservative formulation, at a logarithmic grid,
  ! using a single-stage second order scheme
  subroutine advance_log_advection(CFLIn, nP, nGCLeft, nGCRight,        &
       FInOut_I, IsConservative, DeltaLnP)

    real,   intent(in):: CFLIn        ! Time step * acceleration rate/(Dlnp)
    integer,intent(in):: nP           ! Number of meshes along lnp-coordinate
    integer,intent(in):: nGCLeft      ! The solution in the ghost cells is not
    integer,intent(in):: nGCRight     ! advanced in time but used as the BCs
    real,intent(inout):: FInOut_I(1-nGCLeft:nP+nGCRight)  ! Solution
    ! sol. index 0 is to set boundary condition at the injection energy
    logical,intent(in):: IsConservative  ! Solve (C) if .true. (NC) otherwise
    real,optional,intent(in)::DeltaLnP   ! Used only for NonConservative=.true.!
    ! Loop variables
    integer :: iP, iStep
    ! Subcycling to achive the condition CFL<1, if nStep>1
    integer :: nStep
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

    if(IsConservative)then
       HalfADtIfNeeded = 0.50*CFLIn*DeltaLnP
       CFL = HalfADtIfneeded/tanh(0.50 * DeltaLnP)
       HalfADtIfNeeded = -HalfADtIfNeeded
    else
       HalfADtIfNeeded = 0.0; CFL = CFLIn
    end if

    nStep = 1 + int(CFL); CFL = CFL/real(nStep)
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
       do iStep=1,nStep

          ! Boundary condition at the left boundary
          if(nGCLeft<2)F_I(            -1:0-nGCLeft) = F_I( 1-nGCLeft )
          ! Boundary condition at the right boundary
          if(nGCRight<2)F_I(nP+1-nGCRight:nP+2     ) = F_I(nP+nGCRight)

          do iP=0,nP

             ! f_(i+1/2):
             FSemiintUp_I(iP) = F_I(iP)*(1.00 + HalfADtIfNeeded)+&
                  0.50*(1.00 - CFL)*df_lim(iP)
          end do
          ! f_(i-1/2):
          FSemiintDown_I(1:nP) = FSemiintUp_I(0:nP-1)

          ! Update the solution from f^(n) to f^(n+1):

          F_I(1:nP) = F_I(1:nP)+CFL*(FSemiintDown_I(1:nP)-FSemiintUp_I(1:nP))+&
               HalfADtIfNeeded*(FSemiintDown_I(1:nP)+FSemiintUp_I(1:nP))
       end do
    else
       do iStep=1,nStep
          ! Boundary condition at the left boundary
          if(nGCLeft<2)F_I(            -1:0-nGCLeft) = F_I( 1-nGCLeft)
          ! Boundary condition at the right boundary
          if(nGCRight<2)F_I(nP+1-nGCRight:nP+2     ) = F_I(nP+nGCRight)
          do iP=1,nP+1
             ! f_(i-1/2):
             FSemiintDown_I(iP) = F_I(iP)*(1.00 + HalfADtIfNeeded)&
                  -0.50 * (1.00 + CFL)*df_lim(iP)
          end do
          ! f_(i+1/2):
          FSemiintUp_I(1:nP) = FSemiintDown_I(2:nP+1)

          ! Update the solution from f^(n) to f^(n+1):
          F_I(1:nP) = F_I(1:nP)+CFL*(FSemiintDown_I(1:nP)-FSemiintUp_I(1:nP))+&
               HalfADtIfNeeded*(FSemiintDown_I(1:nP)+FSemiintUp_I(1:nP))
       end do
    end if

    FInOut_I(1:nP)=F_I(1:nP)
    if(any(FInOut_I(1-nGCLeft:nP+nGCRight)<=0.0))then
       write(*,*)'After advection F_I <=0, for CFLFermi= ',CFL
       write(*,*)F_I
       stop
    end if
    !------------------------------------ DONE -------------------------------!
  Contains
    real function df_lim(i)
      integer, intent(in):: i

      real:: dF1,dF2

      !------------------------------------------------------------------------
      dF1 = F_I(i+1)-F_I(i)
      dF2 = F_I(i)-F_I(i-1)

      ! df_lim=0 if dF1*dF2<0, sign(dF1) otherwise:

      df_lim = sign(0.50,dF1)+sign(0.50,dF2)
      dF1 = abs(dF1)
      dF2 = abs(dF2)
      df_lim = df_lim*min(max(dF1,dF2),2.0*dF1,2.0*dF2)
      !---------------------------------- DONE -------------------------------!
    end function df_lim
    !==========================================================================
  end subroutine advance_log_advection
  !============================================================================
  subroutine advect_via_poisson_bracket(nP, nX, nGCLeft, nGCRight, &
       CflIn, tFinal, InvRhoOld_I, InvRho_I, &
       Momentum3SI_I, VolumeP_I, VolumeX_I, VolumeNewX_I, &
       FInOut_I, DeltaLnP)
    use ModPoissonBracket, ONLY: explicit
    integer,intent(in):: nP, nX       ! Number of meshes along lnp-coordinate
    integer,intent(in):: nGCLeft      ! The solution in the ghost cells is not
    integer,intent(in):: nGCRight     ! advanced in time but used as the BCs
    real,   intent(in):: CFlIn        ! FermiFirsr_I in ModAdvance
    real,   intent(in):: tFinal       ! time interval to advance through
    real,   intent(in):: InvRhoOld_I(nX) ! Old density (inverse)
    real,   intent(in):: InvRho_I(nX)    ! New density (inverse)
    real,   intent(in):: Momentum3SI_I(-1:nP+1) ! Momentum**3/3 in SI
    real,   intent(in):: VolumeP_I(0:nP+1) ! Volume of P in SI
    real,   intent(in):: VolumeX_I(0:nX+1) ! Volume of X in SI
    real,   intent(in):: VolumeNewX_I(0:nX+1) ! New Volume of X in SI
    real,intent(inout):: FInOut_I(1-nGCLeft:nP+nGCRight, 1:nX)  ! Solution
    ! sol. index 0 is to set boundary condition at the injection energy
    real,optional,intent(in)::DeltaLnP  
    ! Loop variables
    integer :: iP, iStep, nStep
    ! Extended arrays for the implementation of the Poisson Bracket Alg.
    real    :: VDF_G(-1:nX+2, -1:nP+2), Dt_C(nX,nP)
    real    :: Volume_G(0:nX+1, 0:nP+1), VolumeNew_G(0:nX+1, 0:nP+1) 
    real    :: dVolumeDt_G(0:nX+1, 0:nP+1), dVolumeXDt_G(0:nX+1)
   !  real    :: Hamiltonian_N(-1:nX+1, -1:nP+1)
    real    :: dHamiltonian02_FY(0:nX+1, -1:nP+1)
    real    :: Source_C(nX, nP)
    real    :: CFL, Time, Dt, DtTrial, DtInv, DtNext
    !--------------------------------------------------------------------------
    ! Now this is the the conservative (-- for the particle number) form
    
    ! Initialization
    nStep = 1 + int(CFlIn); CFL = CFlIn/real(nStep)
    Source_C = 0.0
    VDF_G = 1.0e-8 !VDF_G(:,1) = 1/VolumeP_I(1)
    VDF_G(1:nX, 1-nGCLeft:nP+nGCRight) = transpose(FInOut_I)
    
    ! Time 
    Time = 0.0
    DtTrial = tFinal/nStep
    DtInv = 1/DtTrial
    DtNext = DtTrial

    ! Volume
    do iP = 0, nP + 1
       Volume_G(0:nX+1,iP) = VolumeP_I(iP)*VolumeX_I(0:nX+1)
       VolumeNew_G(0:nX+1,iP) = VolumeP_I(iP)*VolumeNewX_I(0:nX+1)
    end do
    
    ! Advection by Poisson Bracket Alg.
    do iStep = 1, nStep
       Dt = min(DtNext, tFinal - Time); DtInv = 1/Dt
       dVolumeDt_G  = DtInv*(VolumeNew_G  - Volume_G)
       dVolumeXDt_G = DtInv*(VolumeNewX_I - VolumeX_I)
       do iP = -1, nP+1
          dHamiltonian02_FY(:,iP) = - Momentum3SI_I(iP)*dVolumeXDt_G
       end do
       call explicit(nX, nP, VDF_G, Volume_G, Source_C, &
            dHamiltonian02_FY=dHamiltonian02_FY,        &
            dVolumeDt_G = dVolumeDt_G,                  &
            DtIn=Dt,           & ! Input time step, which may be reduced
            CFLIn=CFL,         & ! Input CFL to calculate next time step
            DtOut=DtNext)
       ! Update
       VDF_G(1:nX, 1:nP) = VDF_G(1:nX, 1:nP) + Source_C 
       Time = Time + Dt

       ! BCs
       VDF_G(1:nX, -1:0-nGCLeft) = 1.0e-8
       VDF_G(1:nX,nP+nGCRight:nP+2) = 1.0e-8
       VDF_G( 0,      :) = VDF_G(1,       :)
       VDF_G(-1,      :) = VDF_G(1,       :)
       VDF_G(nX+1:nX+2,:) = 1.0e-8
      !  VDF_G(nX+1:nX+2,1) = 1/VolumeP_I(1)
    end do

    FInOut_I = transpose(VDF_G(1:nX, 1-nGCLeft:nP+nGCRight))
    
  end subroutine advect_via_poisson_bracket
  !============================================================================
end module SP_ModAdvection
