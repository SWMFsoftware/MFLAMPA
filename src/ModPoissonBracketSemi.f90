!  Copyright (C) 2002 Regents of the University of Michigan,
!  portions used with permission
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!
! Semi-implicit Poisson bracket schemes for a single bracket in a 2D
! phase plane, with an optional diffusion along either coordinate:
!
! df/dt + \{f;H\}_{x,y} = div(D grad f).
!
! calc_implicit_source solves the two-level theta method in the increment
! form,
!
! V*(f^{n+1} - f^n)/Dt = alpha*R(f^{n+1}) + (1 - alpha)*R(f^n),
!
! R being the net inflow through the faces of the cell. The time
! weight alpha is 0 (forward, explicit Euler), 1/2 (central,
! trapezoidal or Crank-Nicolson) or 1 (backward, implicit Euler). The
! face value in the Poisson bracket flux \delta H*f_face is weighted
! with sigma:
!
! f_face = (1 - sigma)*f_upwind + sigma*f_downwind,
!
! sigma being 0 (backward in space, upwind), 1/2 (central) or 1
! (forward in space, downwind). For sigma = 0 the limited second order
! correction of explicit4 (Sokolov et al. 2023, JCP 476, 111923) can be
! added. Its Courant weighting, max(0, 1 - 2*alpha)*CFL, recovers
! explicit4 at alpha = 0 and implicit2 at alpha = 1, with the same
! limiter; the weighting (1 - alpha)*CFL, the Courant number of the
! explicit part, makes the total variation decrease in every step at
! 0 < alpha < 1 too, at the cost of the second order in time. The face
! values are limited with the Koren limiter, as in the manuscript
! (UseKorenSemi; .false.: the beta-limiter, the default of
! ModPoissonBracket). The correction is taken at the old time level or,
! with Picard iterations, weighted with alpha as the rest of the
! scheme. Optionally, the face correction from a cell at an extremum
! along the trajectory is bounded by (1 - alpha) times the jump, the
! sufficient TVD condition of the semi-implicit scheme.
! The linear system is solved by banded Gaussian elimination with
! partial pivoting, which also works for the non-M-matrices of the
! central and downwind schemes; the factors are reused while the matrix
! does not change. Alternatively, for the M-matrices of the upwind
! scheme, by line Gauss-Seidel iterations: the Thomas algorithm along
! coordinate 1, the sweeps along coordinate 2 alternating in direction,
! or in one direction with the stopping rule of implicit2.
!
! semi_implicit_one, semi_implicit_two and semi_implicit_tvd are the schemes
! of the manuscript on the semi-implicit Poisson bracket scheme (2026),
! with the arguments of explicit2 and the time weight alpha:
!
! (I + alpha*Dt*L1)(f^{n+1} - f^n) = -Dt*L1 f^n,                     (1)
! (I + alpha*Dt*L2)(f^{n+1} - f^n) = -Dt*L1 f^n + Dt*A(f^n),         (2)
! (I + alpha*Dt*L1)(f^{n+1} - f^n) = -Dt*L1 f^n
!                   + Dt*(alpha*A(f^{n+1}) + (1 - alpha)*A(f^n)),  (TVD)
!
! L1 and L2 being the upwind and the central operators, and A the
! second order correction of explicit4 with the Courant number
! max(0, 1 - 2*alpha)*CFL. In (2) it is not limited: the scheme is
! explicit4 without limiter at alpha = 0 and the theta method with the
! central face values at alpha >= 1/2. In (TVD) it is limited and
! bounded at the extrema as in calc_implicit_source, and Picard iterations
! take it at the new time level. The time step is reduced to keep the
! explicit part stable, (1 - alpha)*CFL <= CflMax, and for (2) also the
! Gauss-Seidel iterations convergent, max(sqrt(1 - 2*alpha), alpha)*CFL
! <= CflMax. The linear systems are solved by symmetric Gauss-Seidel
! iterations, to the round-off level or to a relative tolerance of the
! increment (GsTolIn), or, for a 1D problem, by the Thomas algorithm
! along each line.

module ModPoissonBracketSemi

  use ModUtilities, ONLY: CON_stop
  implicit none
  SAVE
  private ! except
  !PUBLIC MEMBER FUNCTION:
  public :: calc_implicit_source
  ! First order, second order and TVD schemes of the manuscript
  public :: semi_implicit_one, semi_implicit_two, semi_implicit_tvd

  ! The limited correction of all schemes uses the Koren limiter, as the
  ! manuscript; .false.: the beta-limiter (the default of ModPoissonBracket)
  logical, public :: UseKorenSemi = .true.

  character(len=*), parameter:: NameMod = 'ModPoissonBracketSemi'
  ! Tolerance of the increments of the Hamiltonian; of the limiters,
  ! relative to the magnitude of f around the cell
  real, parameter :: cTol = 1.0e-14
  ! Beta-limiter with the minmod beta-parameter, as in ModPoissonBracket
  logical, parameter :: UseMinmodBeta = .true.
  ! Relative convergence tolerance of the Picard iterations
  real, parameter :: PicardTol = 1.0e-10
  ! Relative tolerance and maximum number of Gauss-Seidel sweeps
  real,    parameter :: GsTol = 1.0e-12
  integer, parameter :: nSweepMax = 200
  ! Gauss-Seidel of implicit2: at most 10 sweeps in one direction, until
  ! the change of the lagged flux relative to the diagonal term is 1e-5
  integer, parameter :: nSweepImplicit2 = 10
  real,    parameter :: GsTolImplicit2 = 1.0e-5
  ! Schemes of the manuscript: first order, second order and TVD
  integer, parameter :: FirstOrder_ = 1, SecondOrder_ = 2, Tvd_ = 3
  character(len=17), parameter :: NameScheme_I(3) = &
       ['semi_implicit_one', 'semi_implicit_two', 'semi_implicit_tvd']
  ! If DtIn results in a Courant number of the explicit part above
  ! CflMax, the time step is reduced
  real, parameter :: CflMax = 0.990
  ! Maximum number and relative tolerance of the Picard iterations
  ! of semi_implicit_tvd
  integer, parameter :: nPicardTvd = 100
  real,    parameter :: PicardTolTvd = 1.0e-7
  ! Faces of a cell: left, right, bottom, top
  integer, parameter :: iShift_DS(2,4) = reshape(&
       [-1,  0, 1,  0, 0, -1, 0, 1], [2,4])

  ! Matrix of the last factorization, its band LU factors and pivots
  real,    allocatable :: MatSaved_VC(:,:,:), Lu_II(:,:)
  integer, allocatable :: iPivot_I(:)

contains
  !============================================================================
  subroutine calc_implicit_source(nI, nJ, VDF_G, Volume_G, Source_C,    &
       Hamiltonian12_N, AlphaIn, SigmaIn, UseTvdIn, UseSemiBoundIn, &
       nPicardIn, Diffusion_FX, Diffusion_FY, AlphaDiffIn,          &
       IsSteadyState, DtIn, CFLIn, DtOut, DtOut_C, TypeSolverIn,    &
       nPicardOut, nSweepOut, UseTvdWeightIn, nSweepMaxIn, GsTolIn, &
       GsErrorOut)

    integer, intent(in) :: nI    ! # of cells along coordinate 1
    integer, intent(in) :: nJ    ! # of cells along coordinate 2

    ! Distribution function with two layers of ghost cells, the second
    ! used by the limited correction only. The ghost cell values are
    ! the boundary data, fixed during the time step
    real, intent(in) :: VDF_G(-1:nI+2,-1:nJ+2)

    ! Cell volume. One layer of face ghost cells is used
    real, intent(in) :: Volume_G(0:nI+1,0:nJ+1)

    ! Increment over the time step, f(t+dt) - f(t) = Source_C
    real, intent(out) :: Source_C(1:nI,1:nJ)

    ! Hamiltonian function in nodes, one layer of ghost nodes
    real, intent(in) :: Hamiltonian12_N(-1:nI+1,-1:nJ+1)

    ! Time weight: 0 forward, 0.5 central, 1 backward (default)
    real, optional, intent(in) :: AlphaIn

    ! Face weight: 0 backward (upwind, default), 0.5 central,
    ! 1 forward (downwind)
    real, optional, intent(in) :: SigmaIn

    ! Add the limited second order correction (only with sigma = 0)
    logical, optional, intent(in) :: UseTvdIn

    ! Where \delta^-f and f_j - f have opposite signs, bound the
    ! correction from the cell to its neighbor j by (1-alpha)*|f_j - f|.
    ! Default: .false., as in implicit2
    logical, optional, intent(in) :: UseSemiBoundIn

    ! Picard iterations weighting the correction with alpha. The
    ! default, 0, takes the correction at the old time level
    integer, optional, intent(in) :: nPicardIn

    ! Diffusion coefficients multiplied by the face area over the
    ! distance between the cell centers, as in implicit2: the flux
    ! through the face i+1/2 is Diffusion_FX(i,j)*(f(i+1,j) - f(i,j))
    real, optional, intent(in) :: Diffusion_FX(0:nI,1:nJ)
    real, optional, intent(in) :: Diffusion_FY(1:nI,0:nJ)

    ! Time weight of the diffusion, the default being AlphaIn
    real, optional, intent(in) :: AlphaDiffIn

    ! Local time steps, with CFLIn in each cell
    logical, optional, intent(in) :: IsSteadyState

    real, optional, intent(in)  :: DtIn    ! Options to set time step
    real, optional, intent(in)  :: CFLIn   ! Options to set time step
    real, optional, intent(out) :: DtOut   ! Option to report time step
    real, optional, intent(out) :: DtOut_C(nI,nJ)

    ! Linear solver: 'LU' (default, band LU of the whole system) or line
    ! Gauss-Seidel, converging for the M-matrices of sigma = 0. In every
    ! line Gauss-Seidel variant (prefix 'l') each row j is one line: the
    ! Thomas algorithm solves it along i with the neighboring rows at their
    ! latest values, and the sweeps go over the rows. 'lGS': the sweeps
    ! alternate in direction until the largest change is GsTol of the
    ! largest increment. 'lGS1': the iterations of implicit2, in one
    ! direction,
    ! upward (increasing j) if the strongest coupling to the row below is
    ! stronger than that to the row above, else downward, until the
    ! change of the lagged flux into the row before, relative to the
    ! diagonal term f*M, is below 1e-5 in every cell (cells with f = 0
    ! are skipped), but at most 10 sweeps. 'lGS1U' and 'lGS1D': the same
    ! as 'lGS1', always upward and always downward
    character(len=*), optional, intent(in) :: TypeSolverIn

    ! Number of Picard iterations done and of Gauss-Seidel sweeps in
    ! the last solve, negative if not converged
    integer, optional, intent(out) :: nPicardOut, nSweepOut

    ! Maximum number of Gauss-Seidel sweeps and their tolerance, instead
    ! of 200 and GsTol ('lGS') or 10 and 1e-5 (implicit2), which zero and
    ! negative values keep; with the tolerance 0 the solver does
    ! nSweepMaxIn sweeps
    integer, optional, intent(in) :: nSweepMaxIn
    real,    optional, intent(in) :: GsTolIn

    ! Error measure of the last Gauss-Seidel sweep: the largest change
    ! relative to the largest increment ('lGS') or that of implicit2
    real,    optional, intent(out) :: GsErrorOut

    ! Weight the limited correction with the Courant number of the
    ! explicit part, (1 - alpha)*CFL, instead of max(0, 1 - 2*alpha)*CFL:
    ! the total variation then decreases in every step with the trapezoidal
    ! weight too, but the scheme is first order in time for 0 < alpha < 1
    logical, optional, intent(in) :: UseTvdWeightIn

    ! Face-centered variations of the Hamiltonian function, one layer
    ! of ghost faces. For each cell the value of DeltaH for the face in
    ! the positive direction of i or j is DeltaH_DG(:,i,j), for the
    ! negative direction it is the other sign of the neighbor's value
    real :: DeltaH_DG(2,-1:nI+1,-1:nJ+1)

    ! Sums of the outgoing and incoming DeltaH, the latter being
    ! negative
    real, dimension(0:nI+1,0:nJ+1) :: SumDeltaHPlus_G, SumDeltaHMinus_G

    ! Local time step, CFL = Dt*A/V and its value in the correction
    real, dimension(0:nI+1,0:nJ+1) :: Dt_G, CFL_G, CFLW_G

    ! Coefficients at the four faces of each cell: Poisson bracket,
    ! (1-sigma)*min(0,DeltaH) + sigma*max(0,DeltaH), and diffusion
    real, allocatable :: Coef_SC(:,:,:), Diff_SC(:,:,:)

    ! Matrix: main diagonal (0) and couplings to the neighbors (1:4),
    ! each row divided by V/Dt
    real, allocatable :: Mat_VC(:,:,:)

    ! Net inflow of the linear scheme at time level n, of the limited
    ! correction at level n and at the Picard iterate
    real, dimension(nI,nJ) :: Res_C, Flux2Old_C, Flux2_C, Rhs_C, &
         SourceOld_C
    real :: Fld_G(-1:nI+2,-1:nJ+2)

    real    :: Alpha, Sigma, AlphaDiff, Dt, Tau, Res, DeltaH_S(4)
    logical :: UseTvd, UseSemiBound, UseLocalDt, IsImplicit, UseGs, &
         UseTvdWeight
    integer :: nPicard, iPicard, nB, nSweep, i, j, iSide, iD_D(2)
    ! Gauss-Seidel: the direction of the sweeps (0 alternating, 1 upward,
    ! -1 downward, 2 chosen as in implicit2), their maximum number and
    ! tolerance, and the strongest couplings to the rows below and above
    integer :: iDirSweep, nSweepLim
    real    :: GsTolUse, CouplingDown, CouplingUp
    character(len=5) :: TypeSolver
    character(len=*), parameter:: NameSub = 'calc_implicit_source'
    !--------------------------------------------------------------------------
    Alpha = 1.0
    if(present(AlphaIn)) Alpha = AlphaIn
    Sigma = 0.0
    if(present(SigmaIn)) Sigma = SigmaIn
    AlphaDiff = Alpha
    if(present(AlphaDiffIn)) AlphaDiff = AlphaDiffIn
    UseTvd = .false.
    if(present(UseTvdIn)) UseTvd = UseTvdIn
    UseSemiBound = .false.
    if(present(UseSemiBoundIn)) UseSemiBound = UseSemiBoundIn
    UseTvdWeight = .false.
    if(present(UseTvdWeightIn)) UseTvdWeight = UseTvdWeightIn
    nPicard = 0
    if(present(nPicardIn)) nPicard = nPicardIn
    UseLocalDt = .false.
    if(present(IsSteadyState)) UseLocalDt = IsSteadyState
    TypeSolver = 'LU'
    if(present(TypeSolverIn)) TypeSolver = TypeSolverIn
    select case(TypeSolver)
    case('LU', 'lGS')
       ! Band LU, or line Gauss-Seidel sweeps alternating in direction
       iDirSweep = 0
    case('lGS1')
       ! Line Gauss-Seidel of implicit2, direction from the couplings
       iDirSweep = 2
    case('lGS1U')
       ! Line Gauss-Seidel of implicit2, always upward (increasing j)
       iDirSweep = 1
    case('lGS1D')
       ! Line Gauss-Seidel of implicit2, always downward
       iDirSweep = -1
    case default
       call CON_stop(NameSub//': TypeSolverIn should be LU, lGS, lGS1, '// &
            'lGS1U or lGS1D')
    end select
    UseGs = TypeSolver /= 'LU'
    nSweepLim = nSweepMax; GsTolUse = GsTol
    if(UseGs .and. TypeSolver /= 'lGS')then
       nSweepLim = nSweepImplicit2; GsTolUse = GsTolImplicit2
    end if
    if(present(nSweepMaxIn))then
       if(nSweepMaxIn > 0) nSweepLim = nSweepMaxIn
    end if
    if(present(GsTolIn))then
       if(GsTolIn >= 0.0) GsTolUse = GsTolIn
    end if
    if(present(nPicardOut)) nPicardOut = 0
    if(present(nSweepOut)) nSweepOut = 0
    if(present(GsErrorOut)) GsErrorOut = 0.0
    if(min(Alpha, Sigma, AlphaDiff) < 0.0 .or.             &
         max(Alpha, Sigma, AlphaDiff) > 1.0) call CON_stop(&
         NameSub//': AlphaIn, SigmaIn and AlphaDiffIn should be in [0,1]')
    if(UseTvd .and. Sigma > 0.0) call CON_stop(&
         NameSub//': the limited correction needs SigmaIn = 0')
    if(.not.(present(DtIn) .or. present(CFLIn))) call CON_stop(&
         'Either CflIn or DtIn should be provided in '//NameSub)

    ! Bracket {F,H12}_{x,y}, faces as in implicit2:
    ! y
    ! 1---------<---------
    ! |                  |
    ! v       1,1        ^
    ! |                  |
    ! |                  |
    ! 0--------->--------1x
    DeltaH_DG = 0.0
    DeltaH_DG(1,-1:nI+1, 0:nJ+1)          = &
         Hamiltonian12_N(-1:nI+1, 0:nJ+1) - &
         Hamiltonian12_N(-1:nI+1,-1:nJ)
    DeltaH_DG(2, 0:nI+1,-1:nJ+1)          = &
         Hamiltonian12_N(-1:nI  ,-1:nJ+1) - &
         Hamiltonian12_N(0:nI+1 ,-1:nJ+1)
    ! Cleanup
    where(abs(DeltaH_DG) <= cTol) DeltaH_DG = 0.0

    ! Sums of outgoing and incoming DeltaH in physical and ghost cells
    do j = 0, nJ+1; do i = 0, nI+1
       SumDeltaHPlus_G(i,j) = max(0.0, DeltaH_DG(1,i,j)) &
            + max(0.0, -DeltaH_DG(1,i-1,j))              &
            + max(0.0,  DeltaH_DG(2,i,j))                &
            + max(0.0, -DeltaH_DG(2,i,j-1))
       SumDeltaHMinus_G(i,j) = min(0.0, DeltaH_DG(1,i,j)) &
            + min(0.0, -DeltaH_DG(1,i-1,j))               &
            + min(0.0,  DeltaH_DG(2,i,j))                 &
            + min(0.0, -DeltaH_DG(2,i,j-1))
    end do; end do
    ! CFL coefficient, A/V
    CFL_G = -SumDeltaHMinus_G/Volume_G

    ! Time step
    if(present(DtIn))then
       Dt = DtIn
       if(present(DtOut) .and. present(CFLIn)) &
            DtOut = CFLIn/maxval(CFL_G(1:nI,1:nJ))
       Dt_G = Dt
    elseif(UseLocalDt)then
       ! Local time steps; where there is no flow use the longest one
       where(CFL_G > 0.0)
          Dt_G = CFLIn/max(CFL_G, 1e-30)
       elsewhere
          Dt_G = -1.0
       end where
       ! Without any flow the local time steps are not defined
       if(all(Dt_G < 0.0)) call CON_stop(NameSub// &
            ': IsSteadyState needs a nonzero Hamiltonian increment')
       where(Dt_G < 0.0) Dt_G = maxval(Dt_G)
       if(present(DtOut)) DtOut = minval(Dt_G(1:nI,1:nJ))
    else
       Dt = CFLIn/maxval(CFL_G(1:nI,1:nJ))
       if(present(DtOut)) DtOut = Dt
       Dt_G = Dt
    end if
    if(present(DtOut_C)) DtOut_C = Dt_G(1:nI,1:nJ)
    ! Local Courant number and its part in the limited correction
    CFL_G = Dt_G*CFL_G
    if(UseTvdWeight)then
       CFLW_G = (1.0 - Alpha)*CFL_G
    else
       CFLW_G = max(0.0, 1.0 - 2.0*Alpha)*CFL_G
    end if

    ! Coefficients at the faces: left, right, bottom, top
    allocate(Coef_SC(4,nI,nJ), Diff_SC(4,nI,nJ))
    do j = 1, nJ; do i = 1, nI
       DeltaH_S = [-DeltaH_DG(1,i-1,j), DeltaH_DG(1,i,j), &
            -DeltaH_DG(2,i,j-1), DeltaH_DG(2,i,j)]
       ! Coefficient of f_neighbor - f in the net outflow of the cell:
       ! inflow (upwind) and outflow (downwind) parts
       Coef_SC(:,i,j) = (1.0 - Sigma)*min(0.0, DeltaH_S) &
            + Sigma*max(0.0, DeltaH_S)
    end do; end do
    Diff_SC = 0.0
    if(present(Diffusion_FX))then
       Diff_SC(1,:,:) = Diffusion_FX(0:nI-1,1:nJ)
       Diff_SC(2,:,:) = Diffusion_FX(1:nI,  1:nJ)
    end if
    if(present(Diffusion_FY))then
       Diff_SC(3,:,:) = Diffusion_FY(1:nI,0:nJ-1)
       Diff_SC(4,:,:) = Diffusion_FY(1:nI,1:nJ)
    end if

    ! Net inflow at time level n, the ghost cells included
    do j = 1, nJ; do i = 1, nI
       Res = 0.0
       do iSide = 1, 4
          iD_D = [i,j] + iShift_DS(:,iSide)
          Res = Res + (Diff_SC(iSide,i,j) - Coef_SC(iSide,i,j))* &
               (VDF_G(iD_D(1),iD_D(2)) - VDF_G(i,j))
       end do
       Res_C(i,j) = Res
    end do; end do
    Rhs_C = Res_C
    if(UseTvd)then
       call calc_flux2(nI, nJ, VDF_G, DeltaH_DG, SumDeltaHPlus_G,       &
            SumDeltaHMinus_G, CFLW_G, .true., UseSemiBound, 1.0 - Alpha, &
            Flux2Old_C)
       Rhs_C = Rhs_C + Flux2Old_C
    end if

    IsImplicit = Alpha > 0.0 .or. (AlphaDiff > 0.0 .and. &
         (present(Diffusion_FX) .or. present(Diffusion_FY)))
    if(.not.IsImplicit)then
       ! Forward Euler
       Source_C = Dt_G(1:nI,1:nJ)/Volume_G(1:nI,1:nJ)*Rhs_C
       RETURN
    end if

    ! Matrix, each row divided by V/Dt. The couplings to the ghost
    ! cells enter the main diagonal only: their increments are zero
    allocate(Mat_VC(0:4,nI,nJ))
    do j = 1, nJ; do i = 1, nI
       Tau = Dt_G(i,j)/Volume_G(i,j)
       Mat_VC(1:4,i,j) = Tau*(Alpha*Coef_SC(:,i,j) &
            - AlphaDiff*Diff_SC(:,i,j))
       Mat_VC(0,i,j) = 1.0 - sum(Mat_VC(1:4,i,j))
    end do; end do
    ! The direction of the one-way sweeps chosen as in implicit2, whose
    ! couplings are those of Mat_VC times the volume
    if(iDirSweep == 2)then
       CouplingDown = 0.0; CouplingUp = 0.0
       if(nJ > 1)then
          CouplingDown = min(0.0, &
               minval(Mat_VC(3,:,2:nJ)*Volume_G(1:nI,2:nJ)))
          CouplingUp = min(0.0, &
               minval(Mat_VC(4,:,1:nJ-1)*Volume_G(1:nI,1:nJ-1)))
       end if
       iDirSweep = -1
       if(CouplingDown < CouplingUp) iDirSweep = 1
    end if
    ! The shorter coordinate is the fast index, nB is the band width
    nB = min(nI, nJ)
    if(.not.UseGs) call factorize
    Source_C = 0.0
    call solve(Rhs_C, Source_C)

    if(.not.UseTvd .or. Alpha == 0.0) RETURN
    ! The correction at the new time level, weighted with alpha
    do iPicard = 1, nPicard
       Fld_G = VDF_G
       Fld_G(1:nI,1:nJ) = VDF_G(1:nI,1:nJ) + Source_C
       call calc_flux2(nI, nJ, Fld_G, DeltaH_DG, SumDeltaHPlus_G,       &
            SumDeltaHMinus_G, CFLW_G, .true., UseSemiBound, 1.0 - Alpha, &
            Flux2_C)
       Rhs_C = Res_C + (1.0 - Alpha)*Flux2Old_C + Alpha*Flux2_C
       SourceOld_C = Source_C
       call solve(Rhs_C, Source_C)
       if(present(nPicardOut)) nPicardOut = iPicard
       if(maxval(abs(Source_C - SourceOld_C)) <= &
            PicardTol*maxval(abs(Source_C))) EXIT
    end do
  contains
    !==========================================================================
    integer function k_cell(i, j)
      ! Index of the unknown, the shorter coordinate being the fast one
      integer, intent(in) :: i, j
      !------------------------------------------------------------------------
      if(nI >= nJ)then
         k_cell = (i - 1)*nJ + j
      else
         k_cell = (j - 1)*nI + i
      end if
    end function k_cell
    !==========================================================================
    subroutine factorize
      ! Band LU factorization, skipped if the matrix is the same as for
      ! the saved factors
      integer :: i, j, k, iSide, iD_D(2)
      !------------------------------------------------------------------------
      if(allocated(MatSaved_VC))then
         if(size(MatSaved_VC,2) == nI .and. size(MatSaved_VC,3) == nJ)then
            if(all(MatSaved_VC == Mat_VC)) RETURN
         end if
         deallocate(MatSaved_VC, Lu_II, iPivot_I)
      end if
      allocate(MatSaved_VC(0:4,nI,nJ), Lu_II(3*nB+1,nI*nJ), &
           iPivot_I(nI*nJ))
      MatSaved_VC = Mat_VC
      ! Band storage: A(k,k') = Lu_II(2*nB+1+k-k',k')
      Lu_II = 0.0
      do j = 1, nJ; do i = 1, nI
         k = k_cell(i,j)
         Lu_II(2*nB+1,k) = Mat_VC(0,i,j)
         do iSide = 1, 4
            iD_D = [i,j] + iShift_DS(:,iSide)
            ! Skip the ghost cells
            if(any(iD_D < 1) .or. iD_D(1) > nI .or. iD_D(2) > nJ) CYCLE
            Lu_II(2*nB+1+k-k_cell(iD_D(1),iD_D(2)), &
                 k_cell(iD_D(1),iD_D(2))) = Mat_VC(iSide,i,j)
         end do
      end do; end do
      call band_lu(nI*nJ, nB, Lu_II, iPivot_I)
    end subroutine factorize
    !==========================================================================
    subroutine solve(Rhs_C, Sol_C)
      ! Solve the system with the right hand side Rhs_C*Dt/V. Sol_C is
      ! the initial guess for the Gauss-Seidel iterations
      real, intent(in)    :: Rhs_C(nI,nJ)
      real, intent(inout) :: Sol_C(nI,nJ)
      real :: R_C(nI,nJ), Rhs_I(nI), Old_I(nI), dMax, Defect, Den, Error
      integer :: iSweep, i, j, j1, j2, jStep, jPrev, iSidePrev
      !------------------------------------------------------------------------
      R_C = Dt_G(1:nI,1:nJ)/Volume_G(1:nI,1:nJ)*Rhs_C
      if(.not.UseGs)then
         ! The unknowns are ordered with the fast index nB
         if(nI >= nJ)then
            R_C = reshape(transpose(R_C), [nI,nJ])
            call band_solve(nI*nJ, nB, Lu_II, iPivot_I, R_C)
            Sol_C = transpose(reshape(R_C, [nJ,nI]))
         else
            call band_solve(nI*nJ, nB, Lu_II, iPivot_I, R_C)
            Sol_C = R_C
         end if
         RETURN
      end if
      Error = 0.0
      do iSweep = 1, nSweepLim
         ! Alternating sweeps: upward in the odd ones
         jStep = iDirSweep
         if(iDirSweep == 0) jStep = 1 - 2*mod(iSweep - 1, 2)
         if(jStep > 0)then
            j1 = 1;  j2 = nJ
         else
            j1 = nJ; j2 = 1
         end if
         ! The row solved before row j in this sweep, and its coupling to j
         iSidePrev = 4
         if(jStep < 0) iSidePrev = 3
         dMax = 0.0; Defect = 0.0
         do j = j1, j2, jStep
            Rhs_I = R_C(:,j)
            if(j > 1)  Rhs_I = Rhs_I - Mat_VC(3,:,j)*Sol_C(:,j-1)
            if(j < nJ) Rhs_I = Rhs_I - Mat_VC(4,:,j)*Sol_C(:,j+1)
            Old_I = Sol_C(:,j)
            call thomas_algorithm(nI, Mat_VC(1,:,j), Mat_VC(0,:,j), &
                 Mat_VC(2,:,j), Rhs_I, Sol_C(:,j))
            dMax = max(dMax, maxval(abs(Sol_C(:,j) - Old_I)))
            jPrev = j - jStep
            if(iDirSweep == 0 .or. jPrev < 1 .or. jPrev > nJ) CYCLE
            ! implicit2: the change of the flux from row j into the row
            ! before, which used the old values, relative to f*M of row j
            do i = 1, nI
               Den = abs((VDF_G(i,j) + Sol_C(i,j))*Mat_VC(0,i,j) &
                    *Volume_G(i,j))
               if(Den > 0.0) Defect = max(Defect, abs((Old_I(i) &
                    - Sol_C(i,j))*Mat_VC(iSidePrev,i,jPrev)      &
                    *Volume_G(i,jPrev))/Den)
            end do
         end do
         if(iDirSweep == 0)then
            if(maxval(abs(Sol_C)) > 0.0) Error = dMax/maxval(abs(Sol_C))
            if(dMax <= GsTolUse*maxval(abs(Sol_C))) EXIT
         else
            Error = Defect
            if(Defect <= GsTolUse) EXIT
         end if
      end do
      nSweep = iSweep
      if(nSweep > nSweepLim) nSweep = -nSweepLim
      if(present(nSweepOut)) nSweepOut = nSweep
      if(present(GsErrorOut)) GsErrorOut = Error
    end subroutine solve
    !==========================================================================
  end subroutine calc_implicit_source
  !============================================================================
  subroutine semi_implicit_one(nI, nJ, VDF_G, Volume_G, Source_C, Alpha, &
       Hamiltonian12_N, dHamiltonian01_FX, dHamiltonian02_FY,        &
       DtIn, CFLIn, DtOut, DtOut_C, IsPeriodicIn_D, nSweepMaxIn,     &
       GsTolIn, nSweepOut, GsErrorOut)
    ! First order (upwind) scheme, scheme (1) in the header:
    ! alpha = 0 is the explicit, 1/2 the Crank-Nicolson and 1 the
    ! implicit first order scheme

    integer, intent(in) :: nI, nJ
    real,    intent(in) :: VDF_G(-1:nI+2,-1:nJ+2)
    real,    intent(in) :: Volume_G(0:nI+1,0:nJ+1)
    real,    intent(out):: Source_C(nI,nJ)
    real,    intent(in) :: Alpha ! time weight in [0,1]
    real, optional, intent(in) :: Hamiltonian12_N(-1:nI+1,-1:nJ+1)
    real, optional, intent(in) :: dHamiltonian01_FX(-1:nI+1,0:nJ+1)
    real, optional, intent(in) :: dHamiltonian02_FY(0:nI+1,-1:nJ+1)
    real, optional, intent(inout) :: DtIn
    real, optional, intent(in)    :: CFLIn
    real, optional, intent(out)   :: DtOut
    real, optional, intent(out)   :: DtOut_C(nI,nJ)
    logical, optional, intent(in) :: IsPeriodicIn_D(:)
    ! Gauss-Seidel solver of the 2D problems (see solve_gs): maximum
    ! number of sweep pairs, relative tolerance of the increment (default:
    ! the round-off level), the sweep pairs done (negative if not
    ! converged) and the error bound, all of the last solve
    integer, optional, intent(in)  :: nSweepMaxIn
    real,    optional, intent(in)  :: GsTolIn
    integer, optional, intent(out) :: nSweepOut
    real,    optional, intent(out) :: GsErrorOut
    !--------------------------------------------------------------------------
    call semi_implicit(FirstOrder_, nI, nJ, VDF_G, Volume_G, Source_C, &
         Alpha, Hamiltonian12_N, dHamiltonian01_FX, dHamiltonian02_FY,  &
         DtIn, CFLIn, DtOut, DtOut_C, IsPeriodicIn_D, nSweepMaxIn,      &
         GsTolIn, nSweepOut, GsErrorOut)
  end subroutine semi_implicit_one
  !============================================================================
  subroutine semi_implicit_two(nI, nJ, VDF_G, Volume_G, Source_C, Alpha, &
       Hamiltonian12_N, dHamiltonian01_FX, dHamiltonian02_FY,        &
       DtIn, CFLIn, DtOut, DtOut_C, IsPeriodicIn_D, nSweepMaxIn,     &
       GsTolIn, nSweepOut, GsErrorOut)
    ! Second order scheme without limiter, scheme (2) in the header:
    ! alpha = 0 is explicit4 without limiter, 1/2 the Crank-Nicolson
    ! and 1 the implicit scheme with the central face values

    integer, intent(in) :: nI, nJ
    real,    intent(in) :: VDF_G(-1:nI+2,-1:nJ+2)
    real,    intent(in) :: Volume_G(0:nI+1,0:nJ+1)
    real,    intent(out):: Source_C(nI,nJ)
    real,    intent(in) :: Alpha ! time weight in [0,1]
    real, optional, intent(in) :: Hamiltonian12_N(-1:nI+1,-1:nJ+1)
    real, optional, intent(in) :: dHamiltonian01_FX(-1:nI+1,0:nJ+1)
    real, optional, intent(in) :: dHamiltonian02_FY(0:nI+1,-1:nJ+1)
    real, optional, intent(inout) :: DtIn
    real, optional, intent(in)    :: CFLIn
    real, optional, intent(out)   :: DtOut
    real, optional, intent(out)   :: DtOut_C(nI,nJ)
    logical, optional, intent(in) :: IsPeriodicIn_D(:)
    ! Gauss-Seidel solver of the 2D problems (see solve_gs): maximum
    ! number of sweep pairs, relative tolerance of the increment (default:
    ! the round-off level), the sweep pairs done (negative if not
    ! converged) and the error bound, all of the last solve
    integer, optional, intent(in)  :: nSweepMaxIn
    real,    optional, intent(in)  :: GsTolIn
    integer, optional, intent(out) :: nSweepOut
    real,    optional, intent(out) :: GsErrorOut
    !--------------------------------------------------------------------------
    call semi_implicit(SecondOrder_, nI, nJ, VDF_G, Volume_G, Source_C, &
         Alpha, Hamiltonian12_N, dHamiltonian01_FX, dHamiltonian02_FY,   &
         DtIn, CFLIn, DtOut, DtOut_C, IsPeriodicIn_D, nSweepMaxIn,       &
         GsTolIn, nSweepOut, GsErrorOut)
  end subroutine semi_implicit_two
  !============================================================================
  subroutine semi_implicit_tvd(nI, nJ, VDF_G, Volume_G, Source_C, Alpha, &
       Hamiltonian12_N, dHamiltonian01_FX, dHamiltonian02_FY,           &
       DtIn, CFLIn, DtOut, DtOut_C, IsPeriodicIn_D, nSweepMaxIn,     &
       GsTolIn, nSweepOut, GsErrorOut)
    ! TVD scheme, scheme (TVD) in the header: alpha = 0 is explicit4,
    ! 1/2 the Crank-Nicolson and 1 the implicit TVD scheme

    integer, intent(in) :: nI, nJ
    real,    intent(in) :: VDF_G(-1:nI+2,-1:nJ+2)
    real,    intent(in) :: Volume_G(0:nI+1,0:nJ+1)
    real,    intent(out):: Source_C(nI,nJ)
    real,    intent(in) :: Alpha ! time weight in [0,1]
    real, optional, intent(in) :: Hamiltonian12_N(-1:nI+1,-1:nJ+1)
    real, optional, intent(in) :: dHamiltonian01_FX(-1:nI+1,0:nJ+1)
    real, optional, intent(in) :: dHamiltonian02_FY(0:nI+1,-1:nJ+1)
    real, optional, intent(inout) :: DtIn
    real, optional, intent(in)    :: CFLIn
    real, optional, intent(out)   :: DtOut
    real, optional, intent(out)   :: DtOut_C(nI,nJ)
    logical, optional, intent(in) :: IsPeriodicIn_D(:)
    ! Gauss-Seidel solver of the 2D problems (see solve_gs): maximum
    ! number of sweep pairs, relative tolerance of the increment (default:
    ! the round-off level), the sweep pairs done (negative if not
    ! converged) and the error bound, all of the last solve
    integer, optional, intent(in)  :: nSweepMaxIn
    real,    optional, intent(in)  :: GsTolIn
    integer, optional, intent(out) :: nSweepOut
    real,    optional, intent(out) :: GsErrorOut
    !--------------------------------------------------------------------------
    call semi_implicit(Tvd_, nI, nJ, VDF_G, Volume_G, Source_C,        &
         Alpha, Hamiltonian12_N, dHamiltonian01_FX, dHamiltonian02_FY, &
         DtIn, CFLIn, DtOut, DtOut_C, IsPeriodicIn_D, nSweepMaxIn,     &
         GsTolIn, nSweepOut, GsErrorOut)
  end subroutine semi_implicit_tvd
  !============================================================================
  subroutine semi_implicit(iScheme, nI, nJ, VDF_G, Volume_G, Source_C, &
       Alpha, Hamiltonian12_N, dHamiltonian01_FX, dHamiltonian02_FY,    &
       DtIn, CFLIn, DtOut, DtOut_C, IsPeriodicIn_D, nSweepMaxIn,     &
       GsTolIn, nSweepOut, GsErrorOut)
    ! The schemes (1), (2) and (TVD) of the header

    integer, intent(in) :: iScheme ! FirstOrder_, SecondOrder_ or Tvd_
    integer, intent(in) :: nI, nJ
    real,    intent(in) :: VDF_G(-1:nI+2,-1:nJ+2)
    real,    intent(in) :: Volume_G(0:nI+1,0:nJ+1)
    real,    intent(out):: Source_C(nI,nJ)
    real,    intent(in) :: Alpha
    real, optional, intent(in) :: Hamiltonian12_N(-1:nI+1,-1:nJ+1)
    real, optional, intent(in) :: dHamiltonian01_FX(-1:nI+1,0:nJ+1)
    real, optional, intent(in) :: dHamiltonian02_FY(0:nI+1,-1:nJ+1)
    real, optional, intent(inout) :: DtIn
    real, optional, intent(in)    :: CFLIn
    real, optional, intent(out)   :: DtOut
    real, optional, intent(out)   :: DtOut_C(nI,nJ)
    logical, optional, intent(in) :: IsPeriodicIn_D(:)
    integer, optional, intent(in)  :: nSweepMaxIn
    real,    optional, intent(in)  :: GsTolIn
    integer, optional, intent(out) :: nSweepOut
    real,    optional, intent(out) :: GsErrorOut

    ! Face-centered variations of the Hamiltonian function, as in
    ! calc_implicit_source
    real :: DeltaH_DG(2,-1:nI+1,-1:nJ+1)
    ! Inverse volume, sums of the outgoing and incoming DeltaH, upwind
    ! difference \delta^-f^n, Courant number and its weight in A
    real, dimension(0:nI+1,0:nJ+1) :: vInv_G, SumDeltaHPlus_G, &
         SumDeltaHMinus_G, DeltaMinusF_G, CFL_G, CFLW_G
    ! Couplings of the matrix to the left, right, bottom and top cells
    real, dimension(nI,nJ) :: Aax_II, Ccx_II, Aay_II, Ccy_II
    ! First order explicit increment, Dt*A at time level n and at the
    ! Picard iterate, right hand side, previous iterate
    real, dimension(nI,nJ) :: Source1_C, AdOld_C, Ad_C, Rhs_C, &
         SourceOld_C
    real :: Fld_G(-1:nI+2,-1:nJ+2)
    real :: Dt, CflCoefMax, CapCoef, Coef
    ! Tolerances of the linear solves and of the Picard iterations, the
    ! changes of the last three Picard iterates and their ratio
    real :: GsTolSolve, PicardTolUse, dSource, dSourceOld, dSourceOld2, &
         Rate
    logical :: IsPeriodic_D(2), UseRate
    integer :: nDim, jStart, jLast, i, j, iPicard
    character(len=17) :: NameSub
    !--------------------------------------------------------------------------
    NameSub = NameScheme_I(iScheme)
    if(present(nSweepOut)) nSweepOut = 0
    if(present(GsErrorOut)) GsErrorOut = 0.0
    ! With a tolerance GsTolIn of the increment at least that of the
    ! Picard iterations of (TVD), they stop when their estimated error,
    ! the change times q/(1 - q), q being the larger ratio of successive
    ! changes among the last three (at most 0.9), is GsTolIn of the
    ! increment, and their linear solves at a hundredth of it. A
    ! smaller GsTolIn > 0 only tightens the linear solves, to GsTolIn or
    ! a hundredth of the Picard tolerance
    GsTolSolve = -1.0
    if(present(GsTolIn)) GsTolSolve = GsTolIn
    UseRate = iScheme == Tvd_ .and. GsTolSolve >= PicardTolTvd
    PicardTolUse = PicardTolTvd
    if(UseRate)then
       PicardTolUse = GsTolSolve
       GsTolSolve = 0.01*PicardTolUse
    elseif(iScheme == Tvd_ .and. GsTolSolve > 0.0)then
       GsTolSolve = min(GsTolSolve, 0.01*PicardTolTvd)
    end if
    if(Alpha < 0.0 .or. Alpha > 1.0) call CON_stop(trim(NameSub)// &
         ': alpha must be in [0,1]')
    if(.not.(present(DtIn) .or. present(CFLIn))) call CON_stop( &
         trim(NameSub)//': either DtIn or CFLIn must be provided')
    IsPeriodic_D = .false.
    if(present(IsPeriodicIn_D)) IsPeriodic_D = IsPeriodicIn_D(1:2)
    ! Without Hamiltonian12_N and dHamiltonian02_FY the problem is 1D
    if(present(Hamiltonian12_N) .or. present(dHamiltonian02_FY))then
       nDim = 2; jStart = 0; jLast = nJ + 1
    else
       nDim = 1; jStart = 1; jLast = nJ
    end if
    vInv_G = 1.0/Volume_G

    DeltaH_DG = 0.0
    if(present(Hamiltonian12_N))then
       DeltaH_DG(1,-1:nI+1, 0:nJ+1)          = &
            Hamiltonian12_N(-1:nI+1, 0:nJ+1) - &
            Hamiltonian12_N(-1:nI+1,-1:nJ)
       DeltaH_DG(2, 0:nI+1,-1:nJ+1)          = &
            Hamiltonian12_N(-1:nI  ,-1:nJ+1) - &
            Hamiltonian12_N(0:nI+1 ,-1:nJ+1)
    end if
    if(present(dHamiltonian01_FX)) DeltaH_DG(1,-1:nI+1,jStart:jLast) = &
         DeltaH_DG(1,-1:nI+1,jStart:jLast)                            &
         + dHamiltonian01_FX(-1:nI+1,jStart:jLast)
    if(present(dHamiltonian02_FY)) DeltaH_DG(2,0:nI+1,-1:nJ+1) = &
         DeltaH_DG(2,0:nI+1,-1:nJ+1) + dHamiltonian02_FY
    ! Cleanup
    where(abs(DeltaH_DG) <= cTol) DeltaH_DG = 0.0

    SumDeltaHPlus_G = 0.0; SumDeltaHMinus_G = 0.0
    do j = jStart, jLast; do i = 0, nI+1
       SumDeltaHPlus_G(i,j) = max(0.0, DeltaH_DG(1,i,j)) &
            + max(0.0, -DeltaH_DG(1,i-1,j))              &
            + max(0.0,  DeltaH_DG(2,i,j))                &
            + max(0.0, -DeltaH_DG(2,i,j-1))
       SumDeltaHMinus_G(i,j) = min(0.0, DeltaH_DG(1,i,j)) &
            + min(0.0, -DeltaH_DG(1,i-1,j))               &
            + min(0.0,  DeltaH_DG(2,i,j))                 &
            + min(0.0, -DeltaH_DG(2,i,j-1))
    end do; end do
    call get_delta_minus(nI, nJ, VDF_G, DeltaH_DG, SumDeltaHMinus_G, &
         DeltaMinusF_G)
    ! CFL coefficient, A/V
    CFL_G = -SumDeltaHMinus_G*vInv_G

    ! Time step, reduced if CapCoef*CFL > CflMax
    if(iScheme == SecondOrder_)then
       CapCoef = max(sqrt(max(0.0, 1.0 - 2.0*Alpha)), Alpha)
    else
       CapCoef = 1.0 - Alpha
    end if
    CflCoefMax = maxval(CFL_G(1:nI,1:nJ))
    if(present(DtIn))then
       Dt = DtIn
       if(CapCoef > cTol .and. CapCoef*Dt*CflCoefMax > CflMax)then
          Dt = CflMax/(CapCoef*CflCoefMax)
          DtIn = Dt
       end if
    else
       Dt = CFLIn/max(CflCoefMax, cTol)
       if(CapCoef > cTol .and. CapCoef*Dt*CflCoefMax > CflMax) &
            Dt = CflMax/(CapCoef*CflCoefMax)
    end if
    if(present(DtOut))   DtOut   = Dt
    if(present(DtOut_C)) DtOut_C = Dt
    ! Local Courant number and its part in the second order correction
    CFL_G = Dt*CFL_G
    CFLW_G = max(0.0, 1.0 - 2.0*Alpha)*CFL_G

    ! Explicit first order increment, -Dt*L1 f^n
    Source1_C = -CFL_G(1:nI,1:nJ)*DeltaMinusF_G(1:nI,1:nJ)
    ! Second order correction at time level n, Dt*A(f^n)
    AdOld_C = 0.0
    if(iScheme /= FirstOrder_)then
       call calc_flux2(nI, nJ, VDF_G, DeltaH_DG, SumDeltaHPlus_G,      &
            SumDeltaHMinus_G, CFLW_G, iScheme == Tvd_, Alpha > cTol,   &
            1.0 - Alpha, AdOld_C, IsPeriodic_D)
       AdOld_C = Dt*vInv_G(1:nI,1:nJ)*AdOld_C
    end if
    if(Alpha <= cTol)then
       ! Explicit scheme
       Source_C = Source1_C + AdOld_C
       RETURN
    end if

    if(any(IsPeriodic_D)) call CON_stop(trim(NameSub)// &
         ': periodic boundaries are not supported by the implicit solve')
    ! Couplings of I + alpha*Dt*L, L being L2 for scheme (2), else L1.
    ! Those to the ghost cells are not used: their increments are zero
    do j = 1, nJ; do i = 1, nI
       if(iScheme == SecondOrder_)then
          Coef = 0.5*Alpha*Dt*vInv_G(i,j)
          Aax_II(i,j) = -Coef*DeltaH_DG(1,i-1,j)
          Ccx_II(i,j) =  Coef*DeltaH_DG(1,i,j)
          Aay_II(i,j) = -Coef*DeltaH_DG(2,i,j-1)
          Ccy_II(i,j) =  Coef*DeltaH_DG(2,i,j)
       else
          Coef = Alpha*Dt*vInv_G(i,j)
          Aax_II(i,j) = Coef*min(0.0, -DeltaH_DG(1,i-1,j))
          Ccx_II(i,j) = Coef*min(0.0,  DeltaH_DG(1,i,j))
          Aay_II(i,j) = Coef*min(0.0, -DeltaH_DG(2,i,j-1))
          Ccy_II(i,j) = Coef*min(0.0,  DeltaH_DG(2,i,j))
       end if
    end do; end do
    if(iScheme /= Tvd_)then
       Rhs_C = Source1_C + AdOld_C
       call solve(Rhs_C, Source_C)
       RETURN
    end if

    ! Picard iterations, the limited correction weighted with alpha
    Source_C = 0.0
    dSourceOld = 0.0; dSourceOld2 = 0.0
    do iPicard = 1, nPicardTvd
       Fld_G = VDF_G
       Fld_G(1:nI,1:nJ) = VDF_G(1:nI,1:nJ) + Source_C
       SourceOld_C = Source_C
       call calc_flux2(nI, nJ, Fld_G, DeltaH_DG, SumDeltaHPlus_G, &
            SumDeltaHMinus_G, CFLW_G, .true., .true., 1.0 - Alpha, Ad_C)
       Ad_C = Dt*vInv_G(1:nI,1:nJ)*Ad_C
       Rhs_C = Source1_C + Alpha*Ad_C + (1.0 - Alpha)*AdOld_C
       call solve(Rhs_C, Source_C)
       dSource = maxval(abs(Source_C - SourceOld_C))
       if(UseRate)then
          if(iPicard > 2)then
             Rate = min(0.9, max(dSource/max(dSourceOld, tiny(1.0)), &
                  dSourceOld/max(dSourceOld2, tiny(1.0))))
             if(dSource*Rate/(1 - Rate) < &
                  PicardTolUse*(maxval(abs(Source_C)) + cTol)) RETURN
          end if
          dSourceOld2 = dSourceOld; dSourceOld = dSource
       elseif(dSource < PicardTolUse*(maxval(abs(Source_C)) + cTol))then
          RETURN
       end if
    end do
    call CON_stop(trim(NameSub)//': Picard iterations did not converge')
  contains
    !==========================================================================
    subroutine solve(R_C, S_C)
      ! Solve (I + alpha*Dt*L)*S_C = R_C
      real, intent(in)  :: R_C(nI,nJ)
      real, intent(out) :: S_C(nI,nJ)
      integer :: j
      !------------------------------------------------------------------------
      if(nDim == 2)then
         call solve_gs(nI, nJ, Aax_II, Ccx_II, Aay_II, Ccy_II, R_C, &
              S_C, nSweepMaxIn, GsTolSolve, nSweepOut, GsErrorOut)
      else
         do j = 1, nJ
            call thomas_algorithm(nI, Aax_II(:,j),               &
                 1.0 - Aax_II(:,j) - Ccx_II(:,j), Ccx_II(:,j), R_C(:,j), &
                 S_C(:,j))
         end do
      end if
    end subroutine solve
    !==========================================================================
  end subroutine semi_implicit
  !============================================================================
  subroutine get_delta_minus(nI, nJ, Fld_G, DeltaH_DG, SumDeltaHMinus_G, &
       DeltaMinusF_G)
    ! Upwind difference \delta^-f in the physical and first ghost cells:
    ! the mean of f - f_j over the inflow faces, weighted with DeltaH

    integer, intent(in) :: nI, nJ
    real,    intent(in) :: Fld_G(-1:nI+2,-1:nJ+2)
    real,    intent(in) :: DeltaH_DG(2,-1:nI+1,-1:nJ+1)
    real,    intent(in) :: SumDeltaHMinus_G(0:nI+1,0:nJ+1)
    real,    intent(out):: DeltaMinusF_G(0:nI+1,0:nJ+1)

    real    :: SumDF
    integer :: i, j
    !--------------------------------------------------------------------------
    do j = 0, nJ+1; do i = 0, nI+1
       if(SumDeltaHMinus_G(i,j) == 0.0)then
          DeltaMinusF_G(i,j) = 0.0
          CYCLE
       end if
       SumDF = min(0.0, DeltaH_DG(1,i,j))*(Fld_G(i+1,j) - Fld_G(i,j))  &
            + min(0.0,-DeltaH_DG(1,i-1,j))*(Fld_G(i-1,j) - Fld_G(i,j)) &
            + min(0.0, DeltaH_DG(2,i,j))*(Fld_G(i,j+1) - Fld_G(i,j))   &
            + min(0.0,-DeltaH_DG(2,i,j-1))*(Fld_G(i,j-1) - Fld_G(i,j))
       DeltaMinusF_G(i,j) = -SumDF/SumDeltaHMinus_G(i,j)
    end do; end do
  end subroutine get_delta_minus
  !============================================================================
  subroutine calc_flux2(nI, nJ, Fld_G, DeltaH_DG, SumDeltaHPlus_G, &
       SumDeltaHMinus_G, CFLW_G, UseLimiter, UseBound, OneMAlpha,  &
       Flux2_C, IsPeriodicIn_D)
    ! Second order correction of explicit4 for the state Fld_G with the
    ! Courant number CFLW_G; Flux2_C is the net inflow into the physical
    ! cells. With UseLimiter, the face values are limited (Koren or beta,
    ! UseKorenSemi) and, if UseBound, the correction from a cell at an
    ! extremum along the trajectory to its neighbor j is bounded by
    ! OneMAlpha*|f_j - f|. Otherwise they are half differences, as in
    ! explicit4 without limiter.

    integer, intent(in) :: nI, nJ
    real,    intent(in) :: Fld_G(-1:nI+2,-1:nJ+2)
    real,    intent(in) :: DeltaH_DG(2,-1:nI+1,-1:nJ+1)
    real, dimension(0:nI+1,0:nJ+1), intent(in) :: SumDeltaHPlus_G, &
         SumDeltaHMinus_G, CFLW_G
    logical, intent(in) :: UseLimiter, UseBound
    real,    intent(in) :: OneMAlpha
    real,    intent(out):: Flux2_C(nI,nJ)
    logical, optional, intent(in) :: IsPeriodicIn_D(2)

    real :: DeltaMinusF_G(0:nI+1,0:nJ+1), Flux2_G(0:nI+1,0:nJ+1)
    ! DeltaH and limited flux for each outflow face
    integer, parameter :: DeltaH_ = 1, Flux_ = 2
    real    :: Buff_VS(DeltaH_:Flux_,4)
    integer :: iSide_S(4), nFlux, iFlux, iDim, iSide, i, j, iD_D(2)
    real    :: SumFluxPlus, SumMajor, DeltaPlusFLimited, Gamma, Cfl, UnitF
    logical :: IsPeriodic_D(2)
    !--------------------------------------------------------------------------
    IsPeriodic_D = .false.
    if(present(IsPeriodicIn_D)) IsPeriodic_D = IsPeriodicIn_D
    call get_delta_minus(nI, nJ, Fld_G, DeltaH_DG, SumDeltaHMinus_G, &
         DeltaMinusF_G)
    Flux2_G = 0.0
    do j = 1, nJ; do i = 1, nI
       ! Limit and store fluxes across delta plus H faces
       if(SumDeltaHPlus_G(i,j) == 0.0) CYCLE
       ! The power of two above the magnitude of f around the cell: the
       ! limiters, homogeneous of degree one, are applied to f/UnitF, so
       ! that their tolerance cTol is relative and the scaling is exact
       UnitF = scale(1.0, exponent(max(abs(Fld_G(i,j)), abs(Fld_G(i-1,j)), &
            abs(Fld_G(i+1,j)), abs(Fld_G(i,j-1)), abs(Fld_G(i,j+1)),     &
            tiny(1.0))))
       nFlux = 0
       do iDim = 1, 2
          ! Faces i/i+1 and j/j+1
          iSide = 2*iDim
          if(DeltaH_DG(iDim,i,j) > 0.0)then
             nFlux = nFlux + 1; iSide_S(nFlux) = iSide
             Buff_VS(:,nFlux) = DeltaH_DG(iDim,i,j)* &
                  [1.0, face_value(iSide,i,j,UnitF)]
          end if
          ! Faces i-1/i and j-1/j
          iSide = 2*iDim - 1
          iD_D = [i,j]; iD_D(iDim) = iD_D(iDim) - 1
          if(DeltaH_DG(iDim,iD_D(1),iD_D(2)) < 0.0)then
             nFlux = nFlux + 1; iSide_S(nFlux) = iSide
             Buff_VS(:,nFlux) = -DeltaH_DG(iDim,iD_D(1),iD_D(2))* &
                  [1.0, face_value(iSide,i,j,UnitF)]
          end if
       end do
       ! Gamma-limiter in the given cell
       SumFluxPlus = sum(Buff_VS(Flux_,1:nFlux))
       ! The value of limited \delta^+f/2 to be achieved
       DeltaPlusFLimited = limited(SumFluxPlus/SumDeltaHPlus_G(i,j), &
            DeltaMinusF_G(i,j), UnitF)
       ! The major flux, having the same sign as SumFluxPlus
       SumMajor = sum(Buff_VS(Flux_,1:nFlux), &
            MASK=SumFluxPlus*Buff_VS(Flux_,1:nFlux) > 0.0)
       if(abs(SumMajor) > 0.0)then
          Gamma = 1.0 + (DeltaPlusFLimited*SumDeltaHPlus_G(i,j) &
               - SumFluxPlus)/SumMajor
          where(SumFluxPlus*Buff_VS(Flux_,1:nFlux) > 0.0) &
               Buff_VS(Flux_,1:nFlux) = Buff_VS(Flux_,1:nFlux)*Gamma
       end if
       ! Give the fluxes to the neighbors, the rest leaves the cell
       Cfl = CFLW_G(i,j)
       do iFlux = 1, nFlux
          iD_D = [i,j] + iShift_DS(:,iSide_S(iFlux))
          Flux2_G(iD_D(1),iD_D(2)) = Flux2_G(iD_D(1),iD_D(2)) &
               + Buff_VS(Flux_,iFlux)                         &
               - Buff_VS(DeltaH_,iFlux)*Cfl*DeltaPlusFLimited
       end do
       Flux2_G(i,j) = Flux2_G(i,j) &
            - (1.0 - Cfl)*DeltaPlusFLimited*SumDeltaHPlus_G(i,j)
    end do; end do
    ! Second order fluxes across the periodic boundaries, or from the
    ! ghost cells
    if(IsPeriodic_D(1))then
       Flux2_G(1,1:nJ)  = Flux2_G(1,1:nJ)  + Flux2_G(nI+1,1:nJ)
       Flux2_G(nI,1:nJ) = Flux2_G(nI,1:nJ) + Flux2_G(0,1:nJ)
    else
       do j = 1, nJ
          if(DeltaH_DG(1,0,j) > 0.0) Flux2_G(1,j) = Flux2_G(1,j) &
               + DeltaH_DG(1,0,j)*(1.0 - CFLW_G(0,j))*ghost_value(2,0,j)
          if(DeltaH_DG(1,nI,j) < 0.0) Flux2_G(nI,j) = Flux2_G(nI,j) &
               - DeltaH_DG(1,nI,j)*(1.0 - CFLW_G(nI+1,j))           &
               *ghost_value(1,nI+1,j)
       end do
    end if
    if(IsPeriodic_D(2))then
       Flux2_G(1:nI,1)  = Flux2_G(1:nI,1)  + Flux2_G(1:nI,nJ+1)
       Flux2_G(1:nI,nJ) = Flux2_G(1:nI,nJ) + Flux2_G(1:nI,0)
    else
       do i = 1, nI
          if(DeltaH_DG(2,i,0) > 0.0) Flux2_G(i,1) = Flux2_G(i,1) &
               + DeltaH_DG(2,i,0)*(1.0 - CFLW_G(i,0))*ghost_value(4,i,0)
          if(DeltaH_DG(2,i,nJ) < 0.0) Flux2_G(i,nJ) = Flux2_G(i,nJ) &
               - DeltaH_DG(2,i,nJ)*(1.0 - CFLW_G(i,nJ+1))           &
               *ghost_value(3,i,nJ+1)
       end do
    end if
    Flux2_C = Flux2_G(1:nI,1:nJ)
  contains
    !==========================================================================
    real function face_value(iSide, i, j, UnitF)
      ! Half difference toward the neighbor across iSide, limited with
      ! UseLimiter; zero below cTol*UnitF
      integer, intent(in) :: iSide, i, j
      real,    intent(in) :: UnitF
      integer :: iD_D(2), iU_D(2)
      real    :: DeltaF
      !------------------------------------------------------------------------
      iD_D = [i,j] + iShift_DS(:,iSide)
      iU_D = [i,j] - iShift_DS(:,iSide)
      DeltaF = Fld_G(iD_D(1),iD_D(2)) - Fld_G(i,j)
      if(.not.UseLimiter)then
         face_value = 0.0
         if(abs(DeltaF) > cTol*UnitF) face_value = 0.5*DeltaF
         RETURN
      end if
      face_value = UnitF*betalimiter(                                   &
           DeltaF = DeltaF/UnitF,                                       &
           UpwindDeltaF = (Fld_G(i,j) - Fld_G(iU_D(1),iU_D(2)))/UnitF,  &
           Cfl = CFLW_G(i,j),                                           &
           DownwindDeltaMinusF = DeltaMinusF_G(iD_D(1),iD_D(2))/UnitF,  &
           UpwindDeltaMinusF = DeltaMinusF_G(i,j)/UnitF)
      ! Semi-implicit bound at an extremum along the trajectory
      if(UseBound .and. DeltaMinusF_G(i,j)*DeltaF < 0.0 .and. &
           abs(face_value) > OneMAlpha*abs(DeltaF))           &
           face_value = sign(OneMAlpha*abs(DeltaF), face_value)
    end function face_value
    !==========================================================================
    real function limited(Arg1, Arg2, UnitF)
      ! minmod with UseLimiter, otherwise Arg1; zero below cTol*UnitF
      real, intent(in) :: Arg1, Arg2, UnitF
      !------------------------------------------------------------------------
      if(UseLimiter)then
         limited = UnitF*minmod(Arg1/UnitF, Arg2/UnitF)
      elseif(abs(Arg1) <= cTol*UnitF)then
         limited = 0.0
      else
         limited = Arg1
      end if
    end function limited
    !==========================================================================
    real function ghost_value(iSide, i, j)
      ! Limited face value of the ghost cell (i,j) across iSide, with the
      ! magnitude of f along the face direction
      integer, intent(in) :: iSide, i, j
      integer :: iD_D(2), iU_D(2)
      real    :: UnitF
      !------------------------------------------------------------------------
      iD_D = [i,j] + iShift_DS(:,iSide)
      iU_D = [i,j] - iShift_DS(:,iSide)
      UnitF = scale(1.0, exponent(max(abs(Fld_G(i,j)),          &
           abs(Fld_G(iD_D(1),iD_D(2))), abs(Fld_G(iU_D(1),iU_D(2))), &
           tiny(1.0))))
      ghost_value = limited(face_value(iSide,i,j,UnitF), &
           DeltaMinusF_G(i,j), UnitF)
    end function ghost_value
    !==========================================================================
  end subroutine calc_flux2
  !============================================================================
  subroutine solve_gs(nI, nJ, Aax_II, Ccx_II, Aay_II, Ccy_II, Rhs_II, &
       Sol_II, nSweepMaxIn, GsTolIn, nSweepOut, GsErrorOut)
    ! Solve the five-point system
    !   Diag*s(i,j) + Aax*s(i-1,j) + Ccx*s(i+1,j)
    !               + Aay*s(i,j-1) + Ccy*s(i,j+1) = Rhs(i,j),
    !   Diag = 1 - (Aax + Ccx + Aay + Ccy),
    ! with zero increments in the ghost cells, by symmetric Gauss-Seidel
    ! sweeps. The matrices of the schemes (1) and (TVD) at any time step,
    ! and of the scheme (2) at alpha*CFL < 1, are strictly diagonally
    ! dominant, so the sweeps converge. By default they stop at the
    ! round-off level (backward error) of the residual,
    !   max|Rhs - A*s| <= TolExact*(max|Rhs| + |A|*max|s|).
    ! With the tolerance tau = GsTolIn > 0 they stop when the error of
    ! the increment is at most tau times its largest value, by the bound
    !   max|s - s_exact| <= max|Rhs - A*s|/Margin,
    ! Margin = min(Diag - sum|couplings|): 1 for the upwind matrices of
    ! (1) and (TVD), 1 - alpha*CFL for (2). The sweep pairs needed are
    ! about sweeps_per_e*ln(1/tau), and at most 3*sweeps_per_e*
    ! ln(1/(tau*Margin)) + 20 are done.

    integer, intent(in) :: nI, nJ
    real,    intent(in) :: Aax_II(nI,nJ), Ccx_II(nI,nJ), Aay_II(nI,nJ), &
         Ccy_II(nI,nJ), Rhs_II(nI,nJ)
    real,    intent(out):: Sol_II(nI,nJ)

    ! Maximum number of sweep pairs instead of that of the tolerance; at
    ! this number the iterations stop without an error
    integer, optional, intent(in) :: nSweepMaxIn

    ! Relative tolerance tau: negative (default) for the round-off level,
    ! 0 for exactly nSweepMaxIn sweep pairs
    real,    optional, intent(in) :: GsTolIn

    ! Number of sweep pairs done, negative if not converged, and the
    ! bound of the error relative to the largest value of the increment
    integer, optional, intent(out) :: nSweepOut
    real,    optional, intent(out) :: GsErrorOut

    real, parameter :: TolExact = 1.0e-14
    real, allocatable :: S_G(:,:), DiagInv_II(:,:)
    real    :: sNew, ResidMax, ResidSweep, ScaleF, qJMax, NormA, Tol, &
         Margin
    integer :: i, j, iSweep, iLo, iHi, jLo, jHi, iStr, iIter, MaxIter
    logical :: IsRoundOff, DoStop
    character(len=*), parameter:: NameSub = 'solve_gs'
    !--------------------------------------------------------------------------
    allocate(S_G(0:nI+1,0:nJ+1), DiagInv_II(nI,nJ))
    DiagInv_II = 1.0/(1.0 - (Aax_II + Ccx_II + Aay_II + Ccy_II))
    ! Jacobi contraction factor and norm of the matrix
    qJMax = maxval((abs(Aax_II) + abs(Ccx_II) + abs(Aay_II) &
         + abs(Ccy_II))*DiagInv_II)
    NormA = maxval(1.0/DiagInv_II + abs(Aax_II) + abs(Ccx_II) &
         + abs(Aay_II) + abs(Ccy_II))
    ! Margin of the diagonal dominance
    Margin = minval(1.0/DiagInv_II - abs(Aax_II) - abs(Ccx_II) &
         - abs(Aay_II) - abs(Ccy_II))
    Tol = -1.0
    if(present(GsTolIn)) Tol = GsTolIn
    IsRoundOff = Tol < 0.0
    if(IsRoundOff)then
       MaxIter = 400 + nint(64.0/max(1.0 - qJMax, 1.0e-3))
    elseif(Tol > 0.0)then
       if(Margin <= 0.0) call CON_stop(NameMod//':'//NameSub// &
            ': the matrix is not diagonally dominant')
       ! Three times the pairs of the empirical formula for the bound,
       ! the implicit Courant number being the largest sum of the
       ! couplings
       MaxIter = nint(3*sweeps_per_e(maxval(abs(Aax_II) + abs(Ccx_II) &
            + abs(Aay_II) + abs(Ccy_II)))*log(1/(Tol*Margin))) + 20
    else
       MaxIter = 0
    end if
    DoStop = .true.
    if(present(nSweepMaxIn))then
       if(nSweepMaxIn > 0)then
          MaxIter = nSweepMaxIn; DoStop = .false.
       end if
    end if
    if(MaxIter < 1) call CON_stop(NameMod//':'//NameSub// &
         ': GsTolIn = 0 needs nSweepMaxIn > 0')
    ! The right hand side is the initial guess, the ghost cells are 0
    S_G = 0.0
    S_G(1:nI,1:nJ) = Rhs_II
    do iIter = 1, MaxIter
       ! Forward and backward sweep; the largest residual at the visit,
       ! in both sweeps (ResidMax) and in the backward one (ResidSweep)
       ResidMax = 0.0
       do iSweep = 1, 2
          if(iSweep == 1)then
             iLo = 1;  iHi = nI; jLo = 1;  jHi = nJ; iStr = 1
          else
             iLo = nI; iHi = 1;  jLo = nJ; jHi = 1;  iStr = -1
          end if
          ResidSweep = 0.0
          do j = jLo, jHi, iStr; do i = iLo, iHi, iStr
             sNew = (Rhs_II(i,j)                                     &
                  - Aax_II(i,j)*S_G(i-1,j) - Ccx_II(i,j)*S_G(i+1,j)  &
                  - Aay_II(i,j)*S_G(i,j-1) - Ccy_II(i,j)*S_G(i,j+1)) &
                  *DiagInv_II(i,j)
             ResidSweep = max(ResidSweep, &
                  abs(sNew - S_G(i,j))/DiagInv_II(i,j))
             S_G(i,j) = sNew
          end do; end do
          ResidMax = max(ResidMax, ResidSweep)
       end do
       if(IsRoundOff)then
          ScaleF = maxval(abs(Rhs_II)) &
               + NormA*maxval(abs(S_G(1:nI,1:nJ))) + tiny(1.0)
          if(ResidMax > TolExact*ScaleF) CYCLE
          ! The residual of the final iterate: the backward error test
          if(residual() <= TolExact*ScaleF) EXIT
       elseif(Tol > 0.0)then
          ! The bound of the error relative to the largest increment, by
          ! the residual of the final iterate, computed once that at the
          ! visit in the backward sweep is within a factor of 10
          ScaleF = Tol*Margin*maxval(abs(S_G(1:nI,1:nJ)))
          if(ResidSweep > 10*ScaleF) CYCLE
          if(residual() <= ScaleF) EXIT
       end if
    end do
    if(iIter > MaxIter)then
       if(DoStop) call CON_stop(NameMod//':'//NameSub// &
            ': Gauss-Seidel iterations did not converge')
       iIter = -MaxIter
    end if
    Sol_II = S_G(1:nI,1:nJ)
    if(present(nSweepOut)) nSweepOut = iIter
    if(present(GsErrorOut)) GsErrorOut = residual() &
         /(max(Margin, tiny(1.0))*(maxval(abs(Sol_II)) + tiny(1.0)))
  contains
    !==========================================================================
    real function residual()
      ! The largest residual of the iterate
      real :: Resid
      integer :: i, j
      !------------------------------------------------------------------------
      residual = 0.0
      do j = 1, nJ; do i = 1, nI
         Resid = Rhs_II(i,j) - S_G(i,j)/DiagInv_II(i,j)        &
              - Aax_II(i,j)*S_G(i-1,j) - Ccx_II(i,j)*S_G(i+1,j) &
              - Aay_II(i,j)*S_G(i,j-1) - Ccy_II(i,j)*S_G(i,j+1)
         residual = max(residual, abs(Resid))
      end do; end do
    end function residual
    !==========================================================================
  end subroutine solve_gs
  !============================================================================
  real function sweeps_per_e(AlphaCfl)
    ! Symmetric Gauss-Seidel sweep pairs of solve_gs per factor e of
    ! the L1 error of the increment, at the implicit Courant number
    ! AlphaCfl (the largest sum of the couplings): 0.25*AlphaCfl**0.68,
    ! fitted for the first order scheme with backward Euler on the
    ! gyration test, AlphaCfl = 0.25 to 64, within a factor of 1.6
    ! (testSemiImplicit.ipynb, Section 8.6). Crank-Nicolson needs the
    ! same at the same AlphaCfl, the rotation test fewer, and the scheme
    ! (2) up to 2.7 times more near AlphaCfl = 1
    real, intent(in) :: AlphaCfl
    !--------------------------------------------------------------------------
    sweeps_per_e = 0.25*AlphaCfl**0.68
  end function sweeps_per_e
  !============================================================================
  subroutine thomas_algorithm(n, Lower_I, Diag_I, Upper_I, Rhs_I, Sol_I)
    ! Solve the tridiagonal system Lower_I(i)*Sol_I(i-1) + Diag_I(i)*Sol_I(i) +
    ! Upper_I(i)*Sol_I(i+1) = Rhs_I(i); Lower_I(1) and Upper_I(n) are not used.
    ! From: Numerical Recipes, Chapter 2.6, p.40.

    integer, intent(in)  :: n
    real,    intent(in)  :: Lower_I(n), Diag_I(n), Upper_I(n), Rhs_I(n)
    real,    intent(out) :: Sol_I(n)

    integer :: i
    real    :: Aux, Aux_I(2:n)
    !--------------------------------------------------------------------------
    Aux = Diag_I(1)
    Sol_I(1) = Rhs_I(1)/Aux
    do i = 2, n
       Aux_I(i) = Upper_I(i-1)/Aux
       Aux = Diag_I(i) - Lower_I(i)*Aux_I(i)
       Sol_I(i) = (Rhs_I(i) - Lower_I(i)*Sol_I(i-1))/Aux
    end do
    do i = n - 1, 1, -1
       Sol_I(i) = Sol_I(i) - Aux_I(i+1)*Sol_I(i+1)
    end do
  end subroutine thomas_algorithm
  !============================================================================
  subroutine band_lu(n, nB, Ab_II, iPivot_I)
    ! LU factorization with partial pivoting of the n*n band matrix with
    ! nB subdiagonals and nB superdiagonals, stored as in LAPACK dgbtrf:
    ! A(i,j) = Ab_II(2*nB+1+i-j,j). The rows 1:nB of Ab_II, zero on input,
    ! receive the fill-in due to the row interchanges.

    integer, intent(in)    :: n, nB
    real,    intent(inout) :: Ab_II(3*nB+1,n)
    integer, intent(out)   :: iPivot_I(n)

    integer :: j, jCol, jP, jU, m, kV
    real    :: Aux
    character(len=*), parameter:: NameSub = 'band_lu'
    !--------------------------------------------------------------------------
    kV = 2*nB
    jU = 1
    do j = 1, n
       ! Number of the subdiagonal elements in column j
       m = min(nB, n - j)
       ! Pivot: the largest element on and below the diagonal
       jP = maxloc(abs(Ab_II(kV+1:kV+1+m,j)), DIM=1)
       iPivot_I(j) = j + jP - 1
       if(Ab_II(kV+jP,j) == 0.0) call CON_stop(NameMod//':'//NameSub// &
            ': singular matrix')
       jU = max(jU, min(j + nB + jP - 1, n))
       if(jP /= 1)then
          ! Interchange the rows j and j+jP-1 in the columns j:jU
          do jCol = j, jU
             Aux = Ab_II(kV+jP+j-jCol,jCol)
             Ab_II(kV+jP+j-jCol,jCol) = Ab_II(kV+1+j-jCol,jCol)
             Ab_II(kV+1+j-jCol,jCol) = Aux
          end do
       end if
       if(m == 0) CYCLE
       ! Multipliers
       Ab_II(kV+2:kV+1+m,j) = Ab_II(kV+2:kV+1+m,j)/Ab_II(kV+1,j)
       ! Update the next columns within the band
       do jCol = j + 1, jU
          Aux = Ab_II(kV+1+j-jCol,jCol)
          if(Aux /= 0.0) Ab_II(kV+2+j-jCol:kV+1+m+j-jCol,jCol) = &
               Ab_II(kV+2+j-jCol:kV+1+m+j-jCol,jCol) - Aux*Ab_II(kV+2:kV+1+m,j)
       end do
    end do
  end subroutine band_lu
  !============================================================================
  subroutine band_solve(n, nB, Ab_II, iPivot_I, Rhs_I)
    ! Solve A*x = Rhs with the factors from band_lu; x overwrites Rhs_I

    integer, intent(in)    :: n, nB, iPivot_I(n)
    real,    intent(in)    :: Ab_II(3*nB+1,n)
    real,    intent(inout) :: Rhs_I(n)

    integer :: j, m, l, i1, kV
    real    :: Aux
    !--------------------------------------------------------------------------
    kV = 2*nB
    ! Forward substitution with the row interchanges
    do j = 1, n - 1
       m = min(nB, n - j)
       l = iPivot_I(j)
       Aux = Rhs_I(l)
       if(l /= j)then
          Rhs_I(l) = Rhs_I(j); Rhs_I(j) = Aux
       end if
       Rhs_I(j+1:j+m) = Rhs_I(j+1:j+m) - Aux*Ab_II(kV+2:kV+1+m,j)
    end do
    ! Back substitution, the upper factor having 2*nB superdiagonals
    do j = n, 1, -1
       Rhs_I(j) = Rhs_I(j)/Ab_II(kV+1,j)
       i1 = max(1, j - kV)
       Rhs_I(i1:j-1) = Rhs_I(i1:j-1) - Rhs_I(j)*Ab_II(kV+1+i1-j:kV,j)
    end do
  end subroutine band_solve
  !============================================================================
  real function minmod(Arg1, Arg2)
    ! Copy of minmod of ModPoissonBracket. Its .not.UseLimiter branch is
    ! dropped (UseLimiter=.true. there; here the callers handle it)
    real, intent(in) :: Arg1, Arg2
    !--------------------------------------------------------------------------
    if(abs(Arg1)<=cTol)then
       minmod = 0.0
    else
       minmod = (sign(0.5, Arg1) + sign(0.5, Arg2))*min(abs(Arg1), abs(Arg2))
    end if
  end function minmod
  !============================================================================
  real function half_beta(Cfl, DownwindDeltaMinusF, UpwindDeltaMinusF)
    ! Copy of half_beta of ModPoissonBracket
    !  To switch between two choices of beta-parameter:
    !  1. minmod, if UseMinmodBeta=.true.
    !  2. DownwindDeltaMinusF  otherwise
    real, intent(in) :: Cfl, DownwindDeltaMinusF, UpwindDeltaMinusF
    !--------------------------------------------------------------------------
    half_beta = 1.0 - 1.0e-14 - Cfl*(1.0 - minmod( &
         DownwindDeltaMinusF/UpwindDeltaMinusF, 1.0))
  end function half_beta
  !============================================================================
  real function betalimiter(DeltaF, UpwindDeltaF,&
       Cfl, DownwindDeltaMinusF, UpwindDeltaMinusF)
    ! Copy of betalimiter of ModPoissonBracket, with UseKorenSemi for its
    ! UseKoren and without its .not.UseLimiter branch (see minmod)
    !  \delta^-f in the neighboring cells
    real :: HalfBeta
    real, intent(in) :: DeltaF                ! f_j -f
    real, intent(in) :: UpwindDeltaF     ! f - f_j^\prime at the opposite face
    real, intent(in) :: Cfl, DownwindDeltaMinusF, UpwindDeltaMinusF
    real :: SignDeltaF, AbsDeltaF, AbsDeltaFLimited
    !--------------------------------------------------------------------------
    if(abs(DeltaF) <=cTol)then
       betalimiter = 0.0
    else
       SignDeltaF = sign(1.0, DeltaF); AbsDeltaF = abs(DeltaF)
       if(UseKorenSemi)then
          AbsDeltaFLimited = max(2*AbsDeltaF + &
               SignDeltaF*UpwindDeltaF,0.0)/6.0
       else
          AbsDeltaFLimited = 0.5*max(AbsDeltaF, SignDeltaF*UpwindDeltaF)
       end if
       if(UseMinmodBeta.and.CFL > 0.0.and.UpwindDeltaMinusF*DeltaF < 0.0)then
          HalfBeta = max(half_beta(Cfl, (1.0 - Cfl)*DownwindDeltaMinusF, &
               Cfl*UpwindDeltaMinusF), &
               1.0 - 1.0e-14 + CFL*UpwindDeltaMinusF/DeltaF)
       else
          HalfBeta = 1.0 - 1.0e-14
       end if
       betalimiter = SignDeltaF*min(AbsDeltaF*HalfBeta,AbsDeltaFLimited)
    end if
  end function betalimiter
  !============================================================================
end module ModPoissonBracketSemi
!==============================================================================
