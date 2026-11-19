!  Copyright (C) 2002 Regents of the University of Michigan,
!  portions used with permission
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!===========================TESTS==============================================
module ModDiffusion

  ! Solve the diffusion term in the Parker equation
  ! Revision history:
  ! Prototype:SP/FLAMPA/src/SP_ModDiffusion.f90
  ! Adapted for the use in MFLAMPA (Dist_I is an input paramater,
  ! fixed contributions to M_I in the end points)-D.Borovikov, 2017
  ! Updated (identation, comments):  I.Sokolov, Dec.17, 2017
  use ModUtilities, ONLY: CON_stop
  implicit none

  PRIVATE
  public:: advance_diffusion1, tridiag
  interface advance_diffusion1
     module procedure advance_diffusion_s
     module procedure advance_diffusion_arr
  end interface advance_diffusion1

contains
  !============================================================================
  subroutine advance_diffusion_s(Dt, GcnL, nI, GcnR, Dist_I, F_I, &
       DOuter_I, DInner_I)

    ! Solve the diffusion equation:
    !         f_t-D_outer(D_inner*f_x)_x=0,
    ! with zero Neumann boundary condition, if GcnL=GcnR=0.
    ! With given Dirichlet CB at the left bondary, if GcnL=1.
    ! With given Dirichlet CB at the rightt bondary, if GcnR=1.
    ! The solution is advanced in time using fully implicit scheme.

    real,   intent(in   ):: Dt     ! Time step
    integer,intent(in   ):: GcnL, GcnR ! Number of left and right ghost cells
    integer,intent(in   ):: nI     ! Number of physical cells meshes
    ! Distance to the next mesh
    real,   intent(in   ):: Dist_I(1-GcnL:nI-1+GcnR)
    ! In:sol.to be advanced (including ghost cells)
    ! Out:advanced sol (only physical cells)
    real,   intent(inout):: F_I(1-GcnL:nI+GcnR)

    ! Laplace multiplier.
    real,   intent(in   ):: DOuter_I(nI)
    ! Diffusion coefficient
    real,   intent(in   ):: DInner_I(1-GcnL:nI+GcnR)
    real :: Dt_I(nI)
    !--------------------------------------------------------------------------
    Dt_I = Dt
    call advance_diffusion_arr(Dt_I, GcnL, nI, GcnR, Dist_I, F_I, &
         DOuter_I, DInner_I)
  end subroutine advance_diffusion_s
  !============================================================================
  subroutine advance_diffusion_arr(Dt_I, GcnL, nI, GcnR, Dist_I, F_I, &
       DOuter_I, DInner_I)

    ! Solve the diffusion equation:
    !         f_t-D_outer(D_inner*f_x)_x=0,
    ! with zero Neumann boundary condition. The solution is advanced in time
    ! using fully implicit scheme.

    use ModNumConst, ONLY: cTiny

    integer,intent(in   ):: GcnL, GcnR ! Number of left and right ghost cells
    integer,intent(in   ):: nI     ! Number of physical cells meshes
    real,   intent(in   ):: Dt_I(nI)   ! Time step
    ! Distance to the next mesh
    real,   intent(in   ):: Dist_I(1-GcnL:nI-1+GcnR)
    ! In:sol.to be advanced (including ghost cells)
    ! Out:advanced sol (only physical cells)
    real,   intent(inout):: F_I(1-GcnL:nI+GcnR)

    ! Laplace multiplier.
    real,   intent(in   ):: DOuter_I(nI)
    ! Diffusion coefficient
    real,   intent(in   ):: DInner_I(1-GcnL:nI+GcnR)
    ! Mesh spacing and face spacing.
    real                 :: DsFace_I(nI)

    ! Main, upper, and lower diagonals.
    real, dimension(nI)  :: Main_I, Upper_I, Lower_I, R_I
    integer:: i
    real:: Aux1, Aux2
    !--------------------------------------------------------------------------
    Lower_I = 0.0; Main_I  = 1.0; Upper_I = 0.0; R_I = F_I(1:nI)
    ! Within the framework of finite volume method, the cell
    ! volume is used, which is proportional to  the distance between
    ! the faces bounding the volume with an index, i, which is half of
    ! sum of distance between meshes i-1 and i (i.e. D_I(i-1) and that
    ! between meshes i and i+1 (which is D_I(i)):
    do i=2, nI-1
       DsFace_I(i) = max(0.5*(Dist_I(i) + Dist_I(i-1)),cTiny)
       ! In flux coordinates, the control volume associated with the
       ! given cell has a cross-section equal to (Magnetic Flux)/B,
       ! where the flux is a constant along the magnetic field line,
       ! set to one hereafter. Therefore, the diffusion equation has
       ! a following form:
       ! (DsFace_i/B_i)(f^(n+1) - f^n) = Flux_(i-1/2) - Flux_(i+1/2),
       ! where the particle density flux should be multiplied by the
       ! cross-section area too (magnetic flux factor is one!):
       !  Flux_(i-1/2) = (diffusion coefficient)_(i-1/2)/B_(i-1/2)*&
       !                 (f^(n+1)_(i-1) - f^(n+1)_i),
       !  The multiplier, DsFace_i/B_i, is denoted as DsFace_i/DOuter_i
       !  The face-centered combination,
       ! f^(n+1)_i-Dt*DOuter_I/DsFace_I*(&
       !     DInner_(i+1/2)*(f^(n+1)_(i+1)-f^(n+1)_i)/DsMesh_(i+1)-&
       !     DInner_(i-1/2)*(f^(n+1)_i -f^(n+1)_(i-1)/DsMesh_(i ))=f^n_i
       Aux1 = Dt_I(i)*DOuter_I(i)*0.50*(DInner_I(i  ) + DInner_I(i+1))/&
            (Dist_I(i)*DsFace_I(i))
       Aux2 = Dt_I(i)*DOuter_I(i)*0.50*(DInner_I(i-1) + DInner_I(i  ))/&
            (Dist_I(i-1)*DsFace_I(i))
       Main_I(i)  = Main_I(i) + Aux1 + Aux2
       Upper_I(i) = -Aux1
       Lower_I(i) = -Aux2
    end do
    ! For i=1:
    if(GcnL==0) then
       Aux1 = Dt_I(1)*DOuter_I(1)*0.50*(DInner_I(1) + DInner_I(2))/&
            (Dist_I(1)**2)
       Main_I( 1) = Main_I(1) + Aux1
       Upper_I(1) = -Aux1
    elseif(GcnL>0) then
       DsFace_I(1) = max(0.5*(Dist_I(1) + Dist_I(0)),cTiny)
       Aux1 = Dt_I(1)*DOuter_I(1)*0.50*(DInner_I(1) + DInner_I(2))/&
            (Dist_I(1)*DsFace_I(1))
       Aux2 = Dt_I(1)*DOuter_I(1)*0.50*(DInner_I(0) + DInner_I(1))/&
            (Dist_I(0)*DsFace_I(1))
       Main_I(1)  = Main_I(1) + Aux1 + Aux2
       Upper_I(1) = -Aux1
       Lower_I(1) = -Aux2
       R_I(1) = R_I(1) - Lower_I(1)*F_I(0)
    else
       call CON_stop('SP_advance_diffusion: ghost cell number can be 0 or 1!')
    end if
    ! For i=nI
    if(GcnR==0) then
       Aux2 = Dt_I(nI)*DOuter_I(nI)*0.50*(DInner_I(nI-1) + DInner_I(nI))/&
            (Dist_I(nI-1)**2)
       Main_I( nI) = Main_I(nI) + Aux2
       Lower_I(nI) = -Aux2
    elseif(GcnR>0) then
       DsFace_I(nI) = max(0.5*(Dist_I(nI) + Dist_I(nI-1)),cTiny)
       Aux1 = Dt_I(nI)*DOuter_I(nI)*0.50*(DInner_I(nI) + DInner_I(nI+1))/&
            (Dist_I(nI)*DsFace_I(nI))
       Aux2 = Dt_I(nI)*DOuter_I(nI)*0.50*(DInner_I(nI-1) + DInner_I(nI))/&
            (Dist_I(nI-1)*DsFace_I(nI))
       Main_I(nI)  = Main_I(nI) + Aux1 + Aux2
       Upper_I(nI) = -Aux1
       Lower_I(nI) = -Aux2
       R_I(nI) = R_I(nI) - Upper_I(nI)*F_I(nI+1)
    else
       call CON_stop('SP_advance_diffusion: ghost cell number can be 0 or 1!')
    end if
    call tridiag(nI,Lower_I,Main_I,Upper_I,R_I,F_I(1:nI))
  end subroutine advance_diffusion_arr
  !============================================================================
  subroutine tridiag(n, L_I, M_I, U_I, R_I, W_I)

    ! Solve tri-diagonal system of equations:
    !  ||m_1 u_1  0....        || ||w_1|| ||r_1||
    !  ||l_2 m_2 u_2...        || ||w_2|| ||r_2||
    !  || 0  l_3 m_3 u_3       ||.||w_3||=||r_3||
    !  ||...                   || ||...|| ||...||
    !  ||.............0 l_n m_n|| ||w_n|| ||r_n||
    ! From: Numerical Recipes, Chapter 2.6, p.40.

    ! input parameters
    integer,            intent(in):: n
    real, dimension(n), intent(in):: L_I, M_I ,U_I ,R_I
    ! Output parameters
    real, intent(out):: W_I(n)
    ! Misc
    integer:: j
    real:: Aux, Aux_I(2:n)
    !--------------------------------------------------------------------------
    if (M_I(1)==0.0) then
       write(*,*)' Error in tridiag: M_I(1)=0'
       stop
    end if
    Aux = M_I(1)
    W_I(1) = R_I(1)/Aux
    do j=2,n
       Aux_I(j) = U_I(j-1)/Aux
       Aux = M_I(j)-L_I(j)*Aux_I(j)
       if (Aux == 0.0) then
          write(*,*)'M_I(j), L_I(j), Aux_I(j) = ',&
               M_I(j),L_I(j),Aux_I(j)
          write(*,*)'  For j=',j
          write(*,*)'Tridiag failed'
          stop
       end if
       W_I(j) = (R_I(j)-L_I(j)*W_I(j-1))/Aux
    end do
    do j=n-1,1,-1
       W_I(j) = W_I(j)-Aux_I(j+1)*W_I(j+1)
    end do
  end subroutine tridiag
  !============================================================================
end module ModDiffusion
!==============================================================================
module ModTestPoissonBracket
  use ModPoissonBracket, ONLY: explicit, implicit2
  use ModPoissonBracketSemi, ONLY: calc_implicit_source, semi_implicit_one, &
       semi_implicit_two, semi_implicit_tvd, UseKorenSemi
  use ModUtilities, ONLY: CON_stop
  use ModNumConst, ONLY: cTwoPi
  use ModPlotFile, ONLY: save_plot_file
  use ModIoUnit, ONLY: io_unit_new
  use ModConst
  implicit none
  ! Pass the VDF from test_dsa_poisson to restart test_dsa_sa_mhd
  integer, parameter :: nX = 1000  !# of Lagrangian points
  integer :: nP = 40    ! 80  - for Fig5.Right Panel!# momentum bins
  real, allocatable ::  VDFOld_G(:,:)
  ! Momentum grid of the shock acceleration tests: the momenta at the
  ! lower and upper boundaries, the floor of the distribution function
  ! in the initial and boundary values, and the bin of the injection,
  ! whose lower face is at p = 1. The nightly tests use 40 bins from 1 to
  ! 100; set_dsa_grid changes them. The spectral index is fitted for
  ! 0.5 < log10(p) < LogPFitMaxDsa
  real    :: pMinDsa = 1.0, pMaxDsa = 100.0, FloorDsa = 1.0e-8
  real    :: LogPFitMaxDsa = 1.5
  integer :: iPInj = 1
  ! The injected particles enter upstream with the distribution function
  ! vdf_inj in the injection bin [p0, p1]. By default it is a unit density
  ! in the bin, and the analytic spectrum behind a thin shock of
  ! compression 4 is C (p0/p)^4, C = 3(p1^4 - p0^4)/(p1^3 - p0^3) (4.25 at
  ! 20 bins per decade). LevelDsa > 0 scales it to C = LevelDsa
  real    :: LevelDsa = 0.0
  ! Linear solver of test_dsa_semi where it uses the Gauss-Seidel
  ! iterations (lGS; lGS1 for those of implicit2, see calc_implicit_source),
  ! their maximum number (0: that of the solver), and the order of the
  ! coordinates: with IsTransposedDsa the solver lines are along the
  ! momentum, as the coordinates of ModAdvancePoisson, instead of x
  character(len=5) :: TypeSolverDsa = 'lGS'
  integer :: nSweepDsa = 0
  logical :: IsTransposedDsa = .false.
  ! Iterations of the explicit scheme in test_dsa_sa_mhd and the CFL of
  ! its local time steps: 1000 at CFL = 0.98 in the nightly test, while
  ! its steady state needs about 8000 (./test_poisson.exe Comprho4Shock
  ! steps=); the steady state depends on the CFL (Comprho4Shock cfl=)
  integer :: nStepSaMhd = 1000
  real :: CflSaMhd = 0.98
  ! Implicit DSA steps and CFL (Comprho4Shock isteps=, cin=)
  integer :: nStepImpl = 1000
  real :: CflImpl = 98.0
  ! test_dsa_sa_mhd and test_dsa_impl write their convergence histories
  ! (./test_poisson.exe Comprho4Shock) in NameTest_history.txt
  logical :: DoHistoryDsa = .false.
  ! Rotation tests of the semi-implicit schemes: width of the domain
  ! and names of the time weights 0, 1/2 and 1
  real, parameter :: WidthRotation = 6.0
  character(len=14), parameter :: NameAlpha_I(3) = &
       ['explicit a=0  ', 'Crank-Nic a0.5', 'implicit a=1  ']
contains
  !============================================================================
  subroutine test_poisson_bracket(tFinal)

    real, intent(in) :: tFinal

    ! Misc:
    ! Gyration in a uniform magnetic field, nQ numper of points
    ! for the momentum grid, nP is number of point over azimuthal angle
    integer, parameter::  nQ = 300,  nP = 360

    ! Loop variables
    integer           ::  iQ, iStep

    ! Momentum max and min, in units of mc
    real, parameter   :: qMax = 10.0, qMin = 0.01

    ! Mesh size in \phi
    real, parameter   :: DeltaPhi = cTwoPi/nP
    real :: MomentumRatio, MomentumMin, MomentumMax
    real :: VDF_G(-1:nQ+2, -1:nP+2), Volume_G(0:nQ+1, 0:nP+1)
    real :: Hamiltonian_N(-1:nQ+1, -1:nP+1)
    real :: LogMomentum_I(0:nQ+1), Momentum2_I(-1:nQ+1)
    real :: Time, Dt, Source_C(nQ,nP)
    !--------------------------------------------------------------------------
    MomentumRatio = exp(log(qMax/qMin)/nQ)
    MomentumMax =  qMin/MomentumRatio
    Momentum2_I(-1) = MomentumMax**2
    Hamiltonian_N(-1,:) = Hamiltonian(Momentum2_I(-1))
    do iQ = 0, nQ+1
       MomentumMin = MomentumMax
       MomentumMax = MomentumMin*MomentumRatio
       Momentum2_I( iQ) = MomentumMax**2
       Hamiltonian_N(iQ,:) = Hamiltonian(Momentum2_I(iQ))
       Volume_G(iQ,:) = 0.5*DeltaPhi*&
            (Momentum2_I(iQ) - Momentum2_I(iQ-1))
       LogMomentum_I(iQ)   = 0.50*log10(MomentumMin*MomentumMax)
    end do
    ! Initial distribution function
    VDF_G = 0.0; VDF_G(1:nQ, 172:189) = 1.0; Source_C = 0.0
    ! Compiutation
    Time = 0.0; iStep = 0
    do
       call explicit(nQ, nP, VDF_G, Volume_G, Source_C, &
            Hamiltonian_N,   &
            CFLIn=0.99, DtOut = Dt, IsPeriodicIn_D=[.false.,.true.])
       iStep = iStep +1
       if(Time + Dt >= tFinal)then
          Dt = tFinal - Time
          call explicit(nQ, nP, VDF_G, Volume_G, Source_C, &
               Hamiltonian_N,   &
               DtIn = Dt, IsPeriodicIn_D=[.false.,.true.])
          VDF_G(1:nQ, 1:nP) = VDF_G(1:nQ, 1:nP) + Source_C
          EXIT
       else
          Time = Time + Dt
          VDF_G(1:nQ, 1:nP) = VDF_G(1:nQ, 1:nP) + Source_C
          ! Periodic boundary conditions
          VDF_G(1:nQ,-1:0 ) = VDF_G(1:nQ, nP-1:nP)
          VDF_G(1:nQ, nP+1:nP+2 ) = VDF_G(1:nQ, 1:2)
       end if
    end do
    call save_plot_file(NameFile='test_poisson.out', &
         TypeFileIn='ascii', TimeIn=tFinal, nStepIn = iStep, &
         NameVarIn='-p_y  p_x VDF', &
         CoordMinIn_D=[10.0**LogMomentum_I(1 ), 0.50*DeltaPhi], &
         CoordMaxIn_D=[10.0**LogMomentum_I(nQ), cTwoPi - 0.50*DeltaPhi], &
         Coord1In_I = 10.0**LogMomentum_I(1:nQ), &
         VarIn_II = VDF_G(1:nQ,1:nP) )
  end subroutine test_poisson_bracket
  !============================================================================
  real function Hamiltonian(P2)
    real, intent(in) :: P2 ! momentum squared
    !--------------------------------------------------------------------------
    Hamiltonian = sqrt(1.0 + P2)
  end function Hamiltonian
  !============================================================================
  subroutine test_poisson_2d(tFinal)

    real, intent(in) :: tFinal

    ! Misc:
    ! Harmonic oscillators, nQ is the number of meshes over coordinate,
    ! nP is number of meshes overr generalized momentum
    integer, parameter::  nQ = 300,  nP = 300

    ! Loop variables
    integer           ::  iQ, iP, iStep

    ! Mesh size
    real, parameter   :: DeltaQ = 24.0/nQ, DeltaP = 24.0/nP
    real :: VDF_G(-1:nQ+2, -1:nP+2)
    real :: Volume_G(0:nQ+1, 0:nP+1)
    real :: Hamiltonian_N(-1:nQ+1, -1:nP+1)
    real :: Time, Dt, Source_C(nQ,nP)
    real :: NormL2Init, NormL2, EnergyInit, Energy, qNode, pNode, Q, P
    !--------------------------------------------------------------------------
    ! Control volume, for a uniform rectangular grid
    Volume_G = DeltaQ*DeltaP
    ! Hamiltonian at the nodes
    do iP = -1, nP+1
       pNode = DeltaP*(iP - nP/2)
       do iQ = -1, nQ+1
          qNode = DeltaQ*(iQ - nQ/2)
          Hamiltonian_N(iQ,iP) = hamiltonian(qNode**2 + pNode**2)
       end do
    end do
    ! Initial distribution function
    VDF_G = 0.0
    do iP = -1, nP+2
       P = DeltaP*(iP - nP/2 - 0.50)
       do iQ = -1, nQ+2
          Q = DeltaQ*(iQ - nQ/2 - 0.50)
          if(abs(P) < -Q*tan(cPi/20.0).and.P**2+Q**2<=100)VDF_G(iQ,iP) = 1.0
       end do
    end do
    Source_C = 0.0
    ! Compiutation
    Time = 0.0; iStep = 0
    do
       call explicit(nQ, nP, VDF_G, Volume_G, Source_C, &
            Hamiltonian_N,   &
            CFLIn=0.99, DtOut = Dt)
       iStep = iStep +1
       if(Time + Dt >= tFinal)then
          Dt = tFinal - Time
          call explicit(nQ, nP, VDF_G, Volume_G, Source_C, &
               Hamiltonian_N,   &
               DtIn = Dt)
          VDF_G(1:nQ, 1:nP) = VDF_G(1:nQ, 1:nP) + Source_C
          EXIT
       else
          Time = Time + Dt
          VDF_G(1:nQ, 1:nP) = VDF_G(1:nQ, 1:nP) + Source_C
          write(*,*)Time
       end if
    end do
    call save_plot_file(NameFile='test_poisson2d.out', &
         TypeFileIn='ascii', TimeIn=tFinal, nStepIn = iStep, &
         NameVarIn='-Py    Px VDF', &
         CoordMinIn_D=[-12.0 + 0.50*DeltaQ, -12.0 +  0.50*DeltaP], &
         CoordMaxIn_D=[12.0 - 0.50*DeltaQ, 12.0 - 0.50*DeltaP], &
         VarIn_II = VDF_G(1:nQ,1:nP))

  end subroutine test_poisson_2d
  !============================================================================
  subroutine test_poisson_2d_smooth(tFinal)
    real, intent(in) :: tFinal

    ! Misc:
    ! Anharmonic oscillators, nQ is the number of meshes over coordinate,
    ! nP is number of meshes overr generalized momentum
    integer, parameter::  nQ = 300,  nP = 300

    ! Loop variables
    integer           ::  iQ, iP, iStep

    ! Mesh size
    real, parameter   :: DeltaQ = 24.0/nQ, DeltaP = 24.0/nP
    real :: VDF_G(-1:nQ+2, -1:nP+2), VDFFinal_G(-1:nQ+2,-1:nP+2)
    real :: Plot_VC(3,1:nQ,1:nP)
    real :: Volume_G(0:nQ+1, 0:nP+1)
    real :: Hamiltonian_N(-1:nQ+1, -1:nP+1)
    real :: Time, Dt, Source_C(nQ,nP), Error
    real :: NormL2Init, NormL2, EnergyInit, Energy, qNode, pNode, Q, P
    real, parameter :: WidthX = 6.0, WidthY = 6.0
    ! Control volume, for a uniform rectangular grid
    !--------------------------------------------------------------------------
    Volume_G = DeltaQ*DeltaP
    ! Hamiltonian at the nodes
    do iP = -1, nP+1
       pNode = DeltaP*(iP - nP/2)
       do iQ = -1, nQ+1
          qNode = DeltaQ*(iQ - nQ/2)
          Hamiltonian_N(iQ,iP) = hamiltonian(qNode**2 + pNode**2)
       end do
    end do
    ! Initial distribution function
    VDF_G = 0.0
    do iP = -1, nP+2
       P = DeltaP*(iP - nP/2 - 0.50)
       do iQ = -1, nQ+2
          Q = DeltaQ*(iQ - nQ/2 - 0.50)
          VDF_G(     iQ,iP) = initial_cap(Q,P)
          VDFFinal_G(iQ,iP) = final_cap(  Q,P)
       end do
    end do
    Plot_VC(1,:,:) = VDF_G(     1:nQ,1:nP)
    Plot_VC(2,:,:) = VDFFinal_G(1:nQ,1:nP)
    Source_C = 0.0
    ! Computation
    Time = 0.0; iStep = 0
    do
       call explicit(nQ, nP, VDF_G, Volume_G, Source_C, &
            Hamiltonian_N,   &
            CFLIn=0.50, DtOut = Dt)
       iStep = iStep +1
       if(Time + Dt >= tFinal)then
          Dt = tFinal - Time
          call explicit(nQ, nP, VDF_G, Volume_G, Source_C, &
               Hamiltonian_N,   &
               DtIn = Dt)
          Time = tFinal
          write(*,*) Time
          VDF_G(1:nQ, 1:nP) = VDF_G(1:nQ, 1:nP) + Source_C
          EXIT
       else
          Time = Time + Dt
          VDF_G(1:nQ, 1:nP) = VDF_G(1:nQ, 1:nP) + Source_C
          write(*,*)Time
       end if
    end do
    Plot_VC(3,1:nQ,1:nP) = VDF_G(1:nQ, 1:nP)
    Error = sum(abs(VDF_G(1:nQ,1:nP) - VDFFinal_G(1:nQ,1:nP))*&
         Volume_G(1:nQ,1:nP))/(WidthX*WidthY)
    write(*,*)'Error=',Error
    call save_plot_file(NameFile='test_poisson_2d_smooth.out', &
         TypeFileIn='real8', TimeIn=tFinal, nStepIn = iStep, &
         NameVarIn='-Py    Px   VDFInit VDFAnal VDFFinal  Error', &
         CoordMinIn_D=[-12.0 + 0.50*DeltaQ, -12.0 +  0.50*DeltaP], &
         CoordMaxIn_D=[12.0 - 0.50*DeltaQ, 12.0 - 0.50*DeltaP], &
         ParamIn_I=[Error], &
         VarIn_VII = Plot_VC)
  contains
    !==========================================================================
    real function initial_cap(x,y)
      real, intent(in) :: x,y
      real, parameter :: xCenter = -6.0, yCenter = 0.0
      !------------------------------------------------------------------------
      initial_cap = 0.0
      if(abs(x  - xCenter)<=WidthX/2.and.abs(y - yCenter)<=WidthY/2)then
         initial_cap = cos(cPi*(x - xCenter)/WidthX)**4*&
              cos(cPi*(y - yCenter)/WidthY)**4
      end if
    end function initial_cap
    !==========================================================================
    real function final_cap(x,y)
      real, intent(in) :: x,y
      real :: CosPhi,SinPhi,Phi, P, P2, OmegaT, xInit, yInit
      !------------------------------------------------------------------------
      P2  = x**2 + y**2
      OmegaT = tFinal/hamiltonian(P2)
      P  = sqrt(P2)
      CosPhi = x/P; SinPhi = y/P
      if(SinPhi>=0.0)then
         Phi =  acos(CosPhi)
      else
         Phi = -acos(CosPhi)
      end if
      xInit = P*cos(Phi + OmegaT)
      yInit = P*sin(Phi + OmegaT)
      final_cap = initial_cap(xInit,yInit)
    end function final_cap
    !==========================================================================
  end subroutine test_poisson_2d_smooth
  !============================================================================
  subroutine test_energy_conservation(tFinal)

    real, intent(in) :: tFinal

    ! Misc:
    ! Harmonic oscillators, nQ is the number of meshes over coordinate,
    ! nP is number of meshes overr generalized momentum
    integer, parameter::  nQ = 360,  nP = 360

    ! Loop variables
    integer           ::  iQ, iP, iStep

    ! Mesh size
    real, parameter   :: DeltaQ = 24.0/nQ, DeltaP = 24.0/nP
    real :: VDF_G(-1:nQ+2, -1:nP+2),VDFInitial_C(nQ,nP),PlotVar_VC(2,nQ,nP)
    real :: Volume_G(0:nQ+1, 0:nP+1)
    real :: Hamiltonian_N(-1:nQ+1, -1:nP+1)
    real :: Energy_C(nQ, nP)
    real :: Time, Dt, Source_C(nQ,nP)
    real :: NormL1Init, NormL1, EnergyInit, Energy, qNode, pNode, Q, P
    real, parameter :: CFL =0.99
    real, parameter:: pWidth = 10.0, qWidth = 1.0
    integer, parameter:: nPower = 4
    logical, parameter:: IsSmooth = .false.
    ! Control volume, for a uniform rectangular grid
    !--------------------------------------------------------------------------
    Volume_G = DeltaQ*DeltaP
    ! Hamiltonian at the nodes
    do iP = -1, nP+1
       pNode = DeltaP*(iP - nP/2)
       do iQ = -1, nQ+1
          qNode = DeltaQ*(iQ - nQ/2)
          Hamiltonian_N(iQ,iP) = 0.5*(qNode**2 + pNode**2)
       end do
    end do
    ! Energy at the cell centers
    do iP = 1, nP
       P = DeltaP*(iP - nP/2 - 0.50)
       do iQ = 1, nQ
          Q = DeltaQ*(iQ - nQ/2  - 0.50)
          Energy_C(iQ,iP) = 0.5*(Q**2 + P**2)
       end do
    end do
    ! Initial distribution function
    VDFInitial_C = 0.0
    do iP = 1, nP; do iQ = 1, nQ
       P = DeltaP*(iP - nP/2 - 0.5)
       Q = DeltaQ*(iQ - nQ/2 - 0.5)
       if(abs(Q) < qWidth .and. abs(P) < pWidth)then
          if(IsSmooth)then
             VDFInitial_C(iQ,iP) = cos(0.5*cPi*Q/qWidth)**nPower &
                  *                cos(0.5*cPi*P/pWidth)**nPower
          else
             VDFInitial_C(iQ,iP) = 1.0
          endif
       end if
    end do; end do

    VDF_G = 0.0; VDF_G(1:nQ, 1:nP) = VDFInitial_C
    Source_C = 0.0
    ! Initial NormL2 and Energy
    NormL1Init = sum(VDFInitial_C)
    EnergyInit = sum(VDFInitial_C*Energy_C)
    ! Compiutation
    Time = 0.0; iStep = 0
    do
       call explicit(nQ, nP, VDF_G, Volume_G, Source_C, &
            Hamiltonian_N,   &
            CFLIn=CFL, DtOut = Dt)
       iStep = iStep +1
       if(Time + Dt >= tFinal)then
          Dt = tFinal - Time
          call explicit(nQ, nP, VDF_G, Volume_G, Source_C, &
               Hamiltonian_N,   &
               DtIn = Dt)
          VDF_G(1:nQ, 1:nP) = VDF_G(1:nQ, 1:nP) + Source_C
          NormL1 = sum( abs(VDF_G(1:nQ,1:nP) - VDFInitial_C) )
          Energy = sum(VDF_G(1:nQ,1:nP)*Energy_C)
          write(*,*)tFinal, &
               NormL1/NormL1Init,&
               Energy/EnergyInit - 1.0
          EXIT
       else
          Time = Time + Dt
          VDF_G(1:nQ, 1:nP) = VDF_G(1:nQ, 1:nP) + Source_C
          NormL1 = sum( abs(VDF_G(1:nQ,1:nP) - VDFInitial_C) )
          Energy = sum(VDF_G(1:nQ,1:nP)*Energy_C)
          write(*,*)Time, &
               NormL1/NormL1Init,&
               Energy/EnergyInit - 1.0
       end if
    end do
    PlotVar_VC(1,:,:) = VDF_G(1:nQ,1:nP)
    PlotVar_VC(2,:,:) = VDFInitial_C
    call save_plot_file(NameFile='test_energy.out', &
         TypeFileIn='real8', TimeIn=tFinal, nStepIn = iStep, &
         NameVarIn='Q    P VDF VDFIni ErrorL1 EnergyDefect', &
         ParamIn_I=[NormL1/NormL1Init, Energy/EnergyInit - 1.0],&
         CoordMinIn_D=[-12.0 + 0.5*DeltaQ, -12.0 + 0.5*DeltaP], &
         CoordMaxIn_D=[12.0 - 0.5*DeltaQ, 12.0 - 0.5*DeltaP], &
         VarIn_VII = PlotVar_VC)

  end subroutine test_energy_conservation
  !============================================================================
  subroutine test_in_action_angle(tFinal)

    real, intent(in) :: tFinal

    ! Misc:
    ! Harmonic oscillators, nQ is the number of meshes over coordinate,
    ! nP is number of meshes overr generalized momentum
    integer, parameter::  nJ = 300,  nPhi = 360

    ! Loop variables
    integer           ::  iJ, iPhi, iStep

    ! Mesh size
    ! Momentum max and min
    real, parameter   :: JMax = 12.0, JMin = 0.01

    ! Mesh size in \phi
    real, parameter   :: DeltaPhi = cTwoPi/nPhi
    real :: MomentumRatio, MomentumMin, MomentumMax
    real :: VDF_G(-1:nJ+2, -1:nPhi+2), Volume_G(0:nJ+1, 0:nPhi+1)
    real :: VDFInitial_C(nJ,nPhi), PlotVar_VC(2,nJ,nPhi)

    real :: Hamiltonian_N(-1:nJ+1, -1:nPhi+1)
    real :: LogMomentum_I(0:nJ+1), Momentum2_I(-1:nJ+1)

    real :: Energy_C(nJ, nPhi)
    real :: Time, Dt, Source_C(nJ,nPhi)
    real :: NormL2Init, NormL2, EnergyInit, Energy, J, Phi, Q, P
    !--------------------------------------------------------------------------
    MomentumRatio = exp(log(JMax/JMin)/nJ)
    MomentumMax =  JMin/MomentumRatio
    Momentum2_I(-1) = MomentumMax**2
    Hamiltonian_N(-1,:) = 0.50*Momentum2_I(-1)
    do iJ = 0, nJ+1
       MomentumMin = MomentumMax
       MomentumMax = MomentumMin*MomentumRatio
       Momentum2_I( iJ) = MomentumMax**2
       Hamiltonian_N(iJ,:) = 0.50*Momentum2_I(iJ)
       Volume_G(iJ,:) = 0.5*DeltaPhi*&
            (Momentum2_I(iJ) - Momentum2_I(iJ-1))
       LogMomentum_I(iJ)   = 0.50*log10(MomentumMin*MomentumMax)
    end do
    ! Energy at the cell centers,  and initial distribution
    VDFInitial_C = 0.0
    do iPhi = 1, nPhi
       Phi = DeltaPhi*(iPhi - 0.50)
       do iJ = 1, nJ
          J = (Momentum2_I(iJ)*Momentum2_I(iJ-1))**0.250
          Energy_C(iJ,iPhi) = 0.5*J**2
          Q = J*cos(Phi)
          P = J*sin(Phi)
          ! Initial distribution function
          if(Q < -1.0.or.Q > 1.0.or.P < -10.0.or.P > 10.0)CYCLE
          VDFInitial_C(iJ,iPhi) = 1.0
       end do
    end do
    VDF_G = 0.0; VDF_G(1:nJ, 1:nPhi) = VDFInitial_C
    ! Periodic boundary conditions
    VDF_G(1:nJ,-1:0 )   = VDF_G(1:nJ, nPhi-1:nPhi)
    VDF_G(1:nJ, nPhi+1:nPhi+2 ) = VDF_G(1:nJ, 1:2)

    Source_C = 0.0
    ! Initial NormL2 and Energy
    NormL2Init = sum(VDFInitial_C**2*Volume_G(1:nJ,1:nPhi))
    EnergyInit = sum(VDFInitial_C*Energy_C*Volume_G(1:nJ,1:nPhi))
    ! Compiutation
    Time = 0.0; iStep = 0
    ! write(*,*)'Time NormL2/NormL2Init Energy/EnergyInit-1'
    do
       call explicit(nJ, nPhi, VDF_G, Volume_G, Source_C, &
            Hamiltonian_N,   &
            CFLIn=0.99, DtOut = Dt)
       iStep = iStep +1
       if(Time + Dt >= tFinal)then
          Dt = tFinal - Time
          call explicit(nJ, nPhi, VDF_G, Volume_G, Source_C, &
               Hamiltonian_N,   &
               DtIn = Dt)
          VDF_G(1:nJ, 1:nPhi) = VDF_G(1:nJ, 1:nPhi) + Source_C
          NormL2 = sum( (VDF_G(1:nJ,1:nPhi) - VDFInitial_C)**2 &
               *Volume_G(1:nJ,1:nPhi))
          Energy = sum(VDF_G(1:nJ,1:nPhi)*Energy_C &
               *Volume_G(1:nJ,1:nPhi))
          write(*,*)tFinal, &
               NormL2/NormL2Init,&
               Energy/EnergyInit - 1.0
          EXIT
       else
          Time = Time + Dt
          VDF_G(1:nJ, 1:nPhi) = VDF_G(1:nJ, 1:nPhi) + Source_C
          ! Periodic boundary conditions
          VDF_G(1:nJ,-1:0 )   = VDF_G(1:nJ, nPhi-1:nPhi)
          VDF_G(1:nJ, nPhi+1:nPhi+2 ) = VDF_G(1:nJ, 1:2)
          NormL2 = sum( (VDF_G(1:nJ,1:nPhi) - VDFInitial_C)**2 &
               *Volume_G(1:nJ,1:nPhi))
          Energy = sum(VDF_G(1:nJ,1:nPhi)*Energy_C &
               *Volume_G(1:nJ,1:nPhi))
          write(*,*)Time, &
               NormL2/NormL2Init,&
               Energy/EnergyInit - 1.0
       end if
    end do
    PlotVar_VC(1,:,:) = VDF_G(1:nJ,1:nPhi)
    PlotVar_VC(2,:,:) = VDFInitial_C
    call save_plot_file(NameFile='test_action_angle.out', &
         TypeFileIn='real8', TimeIn=tFinal, nStepIn = iStep, &
         NameVarIn='Q    P VDF VDFIni ErrorL2 EnergyDefect', &
         CoordMinIn_D=[10.0**LogMomentum_I(1 ), 0.50*DeltaPhi],&
         CoordMaxIn_D=[10.0**LogMomentum_I(nJ), cTwoPi - 0.50*DeltaPhi], &
         ParamIn_I=[NormL2/NormL2Init, Energy/EnergyInit - 1.0],&
         Coord1In_I = 10.0**LogMomentum_I(1:nJ), &
         VarIn_VII = PlotVar_VC)
  end subroutine test_in_action_angle
  !============================================================================
  subroutine set_dsa_grid(nPIn, pMinIn, pMaxIn, FloorIn, LogPFitMaxIn)
    ! Momentum grid of the shock acceleration tests: nPIn bins from
    ! pMinIn to pMaxIn, the injection in the bin whose lower face is at
    ! p = 1, the floor FloorIn of the distribution function, and the
    ! upper end of the fit of the spectral index, log10(p) < LogPFitMaxIn
    integer, intent(in) :: nPIn
    real,    intent(in) :: pMinIn, pMaxIn, FloorIn, LogPFitMaxIn
    !--------------------------------------------------------------------------
    nP = nPIn; pMinDsa = pMinIn; pMaxDsa = pMaxIn; FloorDsa = FloorIn
    LogPFitMaxDsa = LogPFitMaxIn
    iPInj = 1 + nint(nP*log(1.0/pMinDsa)/log(pMaxDsa/pMinDsa))
  end subroutine set_dsa_grid
  !============================================================================
  real function vdf_inj(Momentum3_I)
    ! The distribution function of the injected particles (see LevelDsa)
    ! in the injection bin, whose faces have p^3/3 = Momentum3_I
    real, intent(in) :: Momentum3_I(2)
    real :: p0, p1
    !--------------------------------------------------------------------------
    if(LevelDsa > 0.0)then
       p0 = (3*Momentum3_I(1))**(1.0/3); p1 = (3*Momentum3_I(2))**(1.0/3)
       vdf_inj = LevelDsa*p0**4/(p1**4 - p0**4)
    else
       vdf_inj = 1/(Momentum3_I(2) - Momentum3_I(1))
    end if
  end function vdf_inj
  !============================================================================
  subroutine test_dsa_poisson

    use ModDiffusion

    real,    parameter :: tFinal = 600.00
    real,    parameter :: DtTrial = 1.0
    real ::  DtInv = 1/DtTrial
    real ::  MomentumRatio, MomentumMin, MomentumMax
    real ::  VDF_G(-1:nX+2, -1:nP+2)
    real ::  Volume_G(0:nX+1, 0:nP+1), VolumeNew_G(0:nX+1, 0:nP+1)
    real ::  VolumeX_I(0:nX+1), VolumeNewX_I(0:nX+1)
    real ::  VolumeP_I(0:nP+1)
    real ::  dVolumeDt_G(0:nX+1, 0:nP+1), dVolumeXDt_G(0:nX+1)
    real ::  DOuter_I(nX) = 1.0, dInner_I(0:nX+1) = 20.0 ! 20.0 For Fig5
    real ::  dHamiltonian02_FY(0:nX+1, -1:nP+1)
    real ::  LogMomentum_I(0:nP+1)     ! Cell centered, for plots
    real ::  Momentum3_I(-1:nP+1)
    real ::  Time = 0.0, Dt, DtNext, Source_C(nX,nP)
    real ::  Coord_I(-1:nX+2)  ! Time-dependent coordinate of a mesh
    real ::  Dist_I(-1:nX+1)   ! Distance from mesh i to mesh i+1

    ! Loop variables
    integer :: iX, iP, iStep
    !--------------------------------------------------------------------------
    Source_C(nX,nP) = 0.0

    ! Logarithmic grid in momentum. pMinDsa is the value at the left face
    ! of the first physical cell, pMaxDsa is the value at the right face
    ! off the last cell. The momentum ratio at the faces of each cell is:
    MomentumRatio   = exp(log(pMaxDsa/pMinDsa)/nP)

    ! The momentum at the left fface of 0th ghost cell is:
    MomentumMax     = pMinDsa/MomentumRatio

    ! Calculate the generalized variable p^3/3 at each face, starting from
    ! the left fface of 0th cell
    Momentum3_I(-1) = MomentumMax**3/3

    ! For cells from 0 to nP+1 calculate face values and volume factor:
    do iP = 0, nP + 1
       MomentumMin = MomentumMax
       MomentumMax = MomentumMin*MomentumRatio
       Momentum3_I(iP) = MomentumMax**3/3
       VolumeP_I(iP)   = Momentum3_I(iP) - Momentum3_I(iP-1)
       ! Only for visualization log10 of the cell centered momentum
       LogMomentum_I(iP) = 0.50*log10(MomentumMin*MomentumMax)
    end do
    VDF_G = FloorDsa; VDF_G(:,iPInj) = vdf_inj(Momentum3_I(iPInj-1:iPInj))

    Time = 0.0; iStep = 0; DtNext = DtTrial
    call update_coords(Time)
    do
       Volume_G = VolumeNew_G
       VolumeX_I= VolumeNewX_I
       Dt = min(DtNext, tFinal - Time); DtInv = 1/Dt
       call update_coords(Time + Dt)
       dVolumeDt_G  = DtInv*(VolumeNew_G  - Volume_G)
       dVolumeXDt_G = DtInv*(VolumeNewX_I - VolumeX_I)
       do iP = -1, nP+1
          dHamiltonian02_FY(:,iP) = - Momentum3_I(iP)*dVolumeXDt_G
       end do
       call explicit(nX, nP, VDF_G, Volume_G, Source_C, &
            dHamiltonian02_FY=dHamiltonian02_FY,        &
            dVolumeDt_G = dVolumeDt_G,                  &
            DtIn=Dt,           & ! Input time step, which may be reduced
            CFLIn=0.98,        & ! Input CFL to calculate next time step
            DtOut=DtNext)        ! Calculation of next time step
       ! Correct final volumes, if Dt had been reduced
       VolumeNew_G  = Volume_G  + Dt*dVolumeDt_G
       VolumeNewX_I = VolumeX_I + Dt*dVolumeXDt_G
       iStep = iStep +1
       VDF_G(1:nX, 1:nP) = VDF_G(1:nX, 1:nP) + Source_C
       do iP =1, nP
          call advance_diffusion1(Dt,0,nX,1,Dist_I(1:nX),VDF_G(1:nX+1,&
               iP),DOuter_I(1:nX),DInner_I(1:nX+1))
       end do
       Time = Time + Dt
       if(Time> tFinal - 1.0e-8*DtNext)EXIT
       VDF_G(1:nX,-1:0 ) = FloorDsa
       VDF_G(1:nX,nP+1:nP+2) = FloorDsa
       VDF_G( 0,      :) = VDF_G(1,       :)
       VDF_G(-1,      :) = VDF_G(1,       :)
       VDF_G(nX+1:nX+2,:) = FloorDsa
       VDF_G(nX+1:nX+2,iPInj) = vdf_inj(Momentum3_I(iPInj-1:iPInj))
    end do
    call save_plot_file(NameFile='test_dsa_poisson.out',       &
         TypeFileIn='ascii', TimeIn=tFinal, nStepIn = iStep,   &
         NameVarIn='LogMomentum VDF'  ,                        &
         CoordMinIn_D=[ LogMomentum_I(1 ) ],                   &
         CoordMaxIn_D=[ LogMomentum_I(nP) ],                   &
         StringFormatIn = '(2F16.9)',                          &
         Coord1In_I = LogMomentum_I(1:nP),                     &
         VarIn_I = alog10(VDF_G(500,1:nP)))
    VDFOld_G = VDF_G
  contains
    !==========================================================================
    subroutine update_coords(Time)
      real, intent(in):: Time

      ! Effective width of the shock wave, in cell size
      integer, parameter:: nWidth = 1

      ! terms in the functions, describing the
      ! lagrangian point motion
      real, parameter:: cQuad = 0.375/nWidth, cConst = 0.375*nWidth
      ! Loop variables
      integer :: iX, iP
      real    :: CoordLagr
      ! Shock wave with the copmtression ratio of 4, with width ~1
      !------------------------------------------------------------------------
      do iX = -1, nX +2
         CoordLagr = iX - 0.50 -Time
         if(CoordLagr>=0.0) then
            Coord_I(iX) = CoordLagr + Time
         elseif(CoordLagr>=-real(nWidth)) then
            Coord_I(iX) = CoordLagr + Time + cQuad*CoordLagr**2
         else
            Coord_I(iX) = 0.250*CoordLagr -cConst + Time
         end if
      end do
      Dist_I(-1:nX+1) = Coord_I(0:nX+2) - Coord_I(-1:nX+1)
      VolumeNewX_I(0:nX+1) = 0.50*(Dist_I(-1:nX) + Dist_I(0:nX+1))
      do iP = 0, nP+1
         VolumeNew_G(0:nX+1,iP) = VolumeP_I(iP)*VolumeNewX_I(0:nX+1)
      end do
    end subroutine update_coords
    !==========================================================================
  end subroutine test_dsa_poisson
  !============================================================================
  subroutine test_dsa_sa_mhd

    use ModDiffusion

    integer :: nStep
    real ::  MomentumRatio, MomentumMin, MomentumMax
    real ::  VDF_G(-1:nX+2, -1:nP+2), Dt_C(nX,nP)
    real ::  Hamiltonian_N(-1:nX+1, -1:nP+1)
    real ::  Volume_G(0:nX+1, 0:nP+1)
    real ::  VolumeX_I(-1:nX+2)
    real ::  VolumeP_I(0:nP+1)
    real ::  DOuter_I(nX) = 1.0, dInner_I(0:nX+1) = 20.0
    real ::  LogMomentum_I(0:nP+1)     ! Cell centered, for plots
    real ::  Momentum3_I(-1:nP+1)
    real ::  Source_C(nX,nP)
    real ::  Coord_I(-1:nX+2)  ! Time-dependent coordinate of a mesh
    real ::  Dist_I(-1:nX+1)   ! Distance from mesh i to mesh i+1
    ! Convergence history: f before the iteration, the unit of the file
    ! and the CPU time at the start
    real, allocatable :: VdfOld_C(:,:)
    real ::  CpuStart
    integer :: iUnit

    ! Loop variables
    integer :: iX, iP, iStep
    !--------------------------------------------------------------------------
    nStep = nStepSaMhd
    Source_C(nX,nP) = 0.0
    ! Logarithmic grid in momentum. pMinDsa is the value at the left face
    ! of the first physical cell, pMaxDsa is the value at the right face
    ! off the last cell. The momentum ratio at the faces of each cell is:
    MomentumRatio   = exp(log(pMaxDsa/pMinDsa)/nP)
    ! The momentum at the left fface of 0th ghost cell is:
    MomentumMax     = pMinDsa/MomentumRatio
    ! Calculate the generalized variable p^3/3 at each face, starting from
    ! the left fface of 0th cell
    Momentum3_I(-1) = MomentumMax**3/3
    ! For cells from 0 to nP+1 calculate face values and volume factor:
    do iP = 0, nP + 1
       MomentumMin = MomentumMax
       MomentumMax = MomentumMin*MomentumRatio
       Momentum3_I(iP) = MomentumMax**3/3
       VolumeP_I(iP)   = Momentum3_I(iP) - Momentum3_I(iP-1)
       ! Only for visualization log10 of the cell centered momentum
       LogMomentum_I(iP) = 0.50*log10(MomentumMin*MomentumMax)
    end do
    VDF_G = VDFOld_G

    call update_coords(600.0)
    if(DoHistoryDsa)then
       allocate(VdfOld_C(nX,nP))
       call open_history_dsa('test_dsa_sa_mhd', iUnit, CpuStart)
    end if
    do iStep = 1, nStep
       if(DoHistoryDsa) VdfOld_C = VDF_G(1:nX,1:nP)
       call explicit(nX, nP, VDF_G, Volume_G, Source_C,     &
            Hamiltonian12_N=Hamiltonian_N, CFLIn=CflSaMhd, &
            IsSteadyState=.true., DtOut_C=Dt_C)

       VDF_G(1:nX, 1:nP) = VDF_G(1:nX, 1:nP) + Source_C
       do iP =1, nP
          call advance_diffusion1(Dt_C(:,iP),0,nX,1,Dist_I(1:nX),VDF_G(1:nX+1,&
               iP),DOuter_I(1:nX),DInner_I(1:nX+1))
       end do
       ! Time = Time + Dt
       VDF_G(1:nX,-1:0 ) = FloorDsa
       VDF_G(1:nX,nP+1:nP+2) = FloorDsa
       VDF_G( 0,      :) = VDF_G(1,       :)
       VDF_G(-1,      :) = VDF_G(1,       :)
       VDF_G(nX+1:nX+2,iPInj) = vdf_inj(Momentum3_I(iPInj-1:iPInj))
       if(DoHistoryDsa) call save_history_dsa(iUnit, iStep, nStep, &
            VdfOld_C, VDF_G(1:nX,1:nP), LogMomentum_I(1:nP), CpuStart)
    end do
    if(DoHistoryDsa) close(iUnit)
    call save_plot_file(NameFile='test_dsa_sa_mhd.out', &
         TypeFileIn='ascii', nStepIn = nStep, &
         NameVarIn='LogMomentum VDF'  ,                  &
         CoordMinIn_D=[ LogMomentum_I(1 ) ],             &
         CoordMaxIn_D=[ LogMomentum_I(nP) ],             &
         StringFormatIn = '(2F16.9)',                    &
         Coord1In_I = LogMomentum_I(1:nP),               &
         VarIn_I = alog10(VDF_G(500,1:nP)))
  contains
    !==========================================================================
    subroutine update_coords(Time)

      real, intent(in):: Time

      ! Effective width of the shock wave, in cell size
      integer, parameter:: nWidth = 1

      ! terms in the functions, describing the
      ! lagrangian point motion
      real, parameter:: cQuad = 0.375/nWidth, cConst = 0.375*nWidth
      ! Loop variables
      integer :: iX, iP
      real    :: CoordLagr
      ! Shock wave with the copmtression ratio of 4, with width ~1
      !------------------------------------------------------------------------
      do iX = -1, nX +2
         CoordLagr = iX - 0.50 -Time
         if(CoordLagr>=0.0)then
            Coord_I(iX) = CoordLagr + Time
         elseif(CoordLagr>=-real(nWidth))then
            Coord_I(iX) = CoordLagr + Time + cQuad*CoordLagr**2
         else
            Coord_I(iX) = 0.250*CoordLagr -cConst + Time
         end if
      end do
      Dist_I(-1:nX+1) = Coord_I(0:nX+2) - Coord_I(-1:nX+1)
      VolumeX_I(0:nX+1) = 0.50*(Dist_I(-1:nX) + Dist_I(0:nX+1))
      VolumeX_I(-1)  = Dist_I(-1); VolumeX_I(nX+2) =  Dist_I(nX+1)
      do iP = 0, nP+1
         Volume_G(0:nX+1,iP) = VolumeP_I(iP)*VolumeX_I(0:nX+1)
      end do
      do iP = -1, nP+1
         Hamiltonian_N(-1:nX+1, iP)= -0.50*Momentum3_I(iP)*&
              (VolumeX_I(-1:nX+1) + VolumeX_I(0:nX+2))
      end do

    end subroutine update_coords
    !==========================================================================
  end subroutine test_dsa_sa_mhd
  !============================================================================
  subroutine open_history_dsa(NameTest, iUnit, CpuStart)
    ! Open the file of the convergence history of the test NameTest and
    ! start its clock
    character(len=*), intent(in) :: NameTest
    integer, intent(out) :: iUnit
    real,    intent(out) :: CpuStart
    !--------------------------------------------------------------------------
    iUnit = io_unit_new()
    open(iUnit, file=NameTest//'_history.txt', status='replace')
    write(iUnit,'(a)') &
         '#     step     change     index      p^4 f       cpu'
    call cpu_time(CpuStart)
  end subroutine open_history_dsa
  !============================================================================
  subroutine save_history_dsa(iUnit, iStep, nStep, VdfOld_C, Vdf_C, &
       LogP_I, CpuStart)
    ! A line of the convergence history, for each of the first 100 steps,
    ! every 10th step to 1000, every 100th to 10000 and so on, and for
    ! the last step: the step, the largest relative change of f in it (as
    ! in test_dsa_semi), the spectral index at iX = 500, fitted for
    ! 0.5 < log10(p) < LogPFitMaxDsa, the mean of p^4 f there and the CPU
    ! time of the steps so far
    integer, intent(in) :: iUnit, iStep, nStep
    real,    intent(in) :: VdfOld_C(nX,nP), Vdf_C(nX,nP), LogP_I(nP)
    real,    intent(in) :: CpuStart
    logical :: IsFit_I(nP)
    real :: LogF_I(nP), xMean, yMean, SpectralIndex, Cpu
    !--------------------------------------------------------------------------
    if(mod(iStep, 10**max(0, int(log10(real(iStep))) - 1)) /= 0 &
         .and. iStep /= nStep) RETURN
    IsFit_I = LogP_I > 0.5 .and. LogP_I < LogPFitMaxDsa
    LogF_I = log10(max(Vdf_C(500,:), 1e-30))
    xMean = sum(LogP_I, MASK=IsFit_I)/count(IsFit_I)
    yMean = sum(LogF_I, MASK=IsFit_I)/count(IsFit_I)
    SpectralIndex = sum((LogP_I - xMean)*(LogF_I - yMean), MASK=IsFit_I) &
         /sum((LogP_I - xMean)**2, MASK=IsFit_I)
    call cpu_time(Cpu)
    write(iUnit,'(i10,es11.3,f10.5,es11.4,f10.3)') iStep,          &
         maxval(abs(Vdf_C - VdfOld_C)/(abs(Vdf_C) + 1e-30)),       &
         SpectralIndex,                                            &
         sum(10**(4*LogP_I + LogF_I), MASK=IsFit_I)/count(IsFit_I), &
         Cpu - CpuStart
  end subroutine save_history_dsa
  !============================================================================
  subroutine test_dsa_impl
    real ::  MomentumRatio, MomentumMin, MomentumMax
    real ::  VDF_G(-1:nX+2, -1:nP+2), Dt_C(nX,nP)
    real ::  Hamiltonian_N(-1:nX+1, -1:nP+1)
    real ::  Volume_G(0:nX+1, 0:nP+1)
    real ::  VolumeX_I(-1:nX+2)
    real ::  VolumeP_I(0:nP+1)
    real ::  DOuter_I(nX) = 1.0, dInner_I(0:nX+1) = 20.0, Error
    real ::  Diffusion_FX(0:nX,1:nP)
    real ::  LogMomentum_I(0:nP+1)     ! Cell centered, for plots
    real ::  Momentum3_I(-1:nP+1)
    real ::  Coord_I(-1:nX+2)  ! Time-dependent coordinate of a mesh
    real ::  Dist_I(-1:nX+1)   ! Distance from mesh i to mesh i+1
    ! Convergence history: f before the step, the unit of the file and
    ! the CPU time at the start
    real, allocatable :: VdfOld_C(:,:)
    real :: CpuStart, Source_C(nX,nP)
    integer :: iUnit, nSweep
    logical :: UseKorenSemiOld

    ! Loop variables
    integer :: iX, iP, iStep
    !--------------------------------------------------------------------------

    ! Logarithmic grid in momentum. pMinDsa is the value at the left face
    ! of the first physical cell, pMaxDsa is the value at the right face
    ! off the last cell. The momentum ratio at the faces of each cell is:
    MomentumRatio   = exp(log(pMaxDsa/pMinDsa)/nP)
    ! The momentum at the left fface of 0th ghost cell is:
    MomentumMax     = pMinDsa/MomentumRatio
    ! Calculate the generalized variable p^3/3 at each face, starting from
    ! the left fface of 0th cell
    Momentum3_I(-1) = MomentumMax**3/3
    ! For cells from 0 to nP+1 calculate face values and volume factor:
    do iP = 0, nP + 1
       MomentumMin = MomentumMax
       MomentumMax = MomentumMin*MomentumRatio
       Momentum3_I(iP) = MomentumMax**3/3
       VolumeP_I(iP)   = Momentum3_I(iP) - Momentum3_I(iP-1)
       ! Only for visualization log10 of the cell centered momentum
       LogMomentum_I(iP) = 0.50*log10(MomentumMin*MomentumMax)
    end do
    ! Initial condition for the distribution function
    VDF_G = VDFOld_G
    call update_coords(600.0)
    if(DoHistoryDsa)then
       allocate(VdfOld_C(nX,nP))
       call open_history_dsa('test_dsa_impl', iUnit, CpuStart)
       UseKorenSemiOld = UseKorenSemi
       ! The shock acceleration test uses the Koren limiter, within the
       ! beta bound (TVD with the Koren limiter, as in the manuscript)
       UseKorenSemi = .true.
    end if
    do iStep = 1, nStepImpl
       if(DoHistoryDsa) VdfOld_C = VDF_G(1:nX,1:nP)
       if(DoHistoryDsa) then
          ! Paper benchmark: use the same backward-Euler increments and
          ! spatial line solves as the production steady Parker driver.
          call calc_implicit_source(nX, nP, VDF_G, Volume_G, Source_C, &
               Hamiltonian_N, AlphaIn=1.0, UseTvdIn=.true.,          &
               Diffusion_FX=Diffusion_FX, IsSteadyState=.true.,      &
               CFLIn=CflImpl, TypeSolverIn='lGS', nSweepMaxIn=2000,  &
               nSweepOut=nSweep, GsErrorOut=Error)
          if(nSweep < 0) call CON_stop('DSA benchmark: GS did not converge')
          VDF_G(1:nX,1:nP) = VDF_G(1:nX,1:nP) + Source_C
       else
          ! Keep the original nightly-regression path.
          call implicit2(nX, nP, VDF_G, Volume_G, &
               Hamiltonian12_N=Hamiltonian_N, Diffusion_FX=Diffusion_FX, &
               CFLIn=CflImpl, ErrorOut=Error)
       end if
       write(*,*)'iStep, Error=', iStep, Error
       VDF_G(1:nX,-1:0 ) = FloorDsa
       VDF_G(1:nX,nP+1:nP+2) = FloorDsa
       VDF_G( 0,      :) = VDF_G(1,       :)
       VDF_G(-1,      :) = VDF_G(1,       :)
       VDF_G(nX+1:nX+2,iPInj) = vdf_inj(Momentum3_I(iPInj-1:iPInj))
       if(DoHistoryDsa) call save_history_dsa(iUnit, iStep, nStepImpl, &
            VdfOld_C, VDF_G(1:nX,1:nP), LogMomentum_I(1:nP), CpuStart)
    end do
    if(DoHistoryDsa) then
       close(iUnit)
       UseKorenSemi = UseKorenSemiOld
    end if
    call save_plot_file(NameFile='test_dsa_impl.out', &
         TypeFileIn='ascii', nStepIn = nStepImpl, &
         NameVarIn='LogMomentum VDF'  ,                  &
         CoordMinIn_D=[ LogMomentum_I(1 ) ],             &
         CoordMaxIn_D=[ LogMomentum_I(nP) ],             &
         StringFormatIn = '(2F16.9)',                    &
         Coord1In_I = LogMomentum_I(1:nP),               &
         VarIn_I = alog10(VDF_G(500,1:nP)))
  contains
    !==========================================================================
    subroutine update_coords(Time)

      real, intent(in):: Time

      ! Effective width of the shock wave, in cell size
      integer, parameter:: nWidth = 1

      ! terms in the functions, describing the
      ! lagrangian point motion
      real, parameter:: cQuad = 0.375/nWidth, cConst = 0.375*nWidth
      ! Loop variables
      integer :: iX, iP
      real    :: CoordLagr, DiffusionX_I(0:nX)
      ! Shock wave with the copmtression ratio of 4, with width ~1
      !------------------------------------------------------------------------
      ! Coord_I is the cell-centered coordinate
      do iX = -1, nX +2
         CoordLagr = iX - 0.50 -Time
         if(CoordLagr>=0.0)then
            Coord_I(iX) = CoordLagr + Time
         elseif(CoordLagr>=-real(nWidth))then
            Coord_I(iX) = CoordLagr + Time + cQuad*CoordLagr**2
         else
            Coord_I(iX) = 0.250*CoordLagr -cConst + Time
         end if
      end do
      ! Dist is the distance to the next cell center
      Dist_I(-1:nX+1) = Coord_I(0:nX+2) - Coord_I(-1:nX+1)
      VolumeX_I(0:nX+1) = 0.50*(Dist_I(-1:nX) + Dist_I(0:nX+1))
      VolumeX_I(-1)  = Dist_I(-1); VolumeX_I(nX+2) =  Dist_I(nX+1)
      do iP = 0, nP+1
         Volume_G(0:nX+1,iP) = VolumeP_I(iP)*VolumeX_I(0:nX+1)
      end do
      do iP = -1, nP+1
         Hamiltonian_N(-1:nX+1, iP)= -0.50*Momentum3_I(iP)*&
              (VolumeX_I(-1:nX+1) + VolumeX_I(0:nX+2))
      end do
      ! Get diffusion coefficients
      ! Flot (no-flux) BC at the left boundary:
      DiffusionX_I(0) = 0.0
      do iX = 1, nX
         DiffusionX_I(iX) = &
              DOuter_I(iX)*0.50*(DInner_I(iX) + DInner_I(iX+1))/Dist_I(iX)
      end do
      do iP = 1, nP
         Diffusion_FX(:,iP) =  VolumeP_I(iP)*DiffusionX_I
      end do
    end subroutine update_coords
    !==========================================================================
  end subroutine test_dsa_impl
  !============================================================================
  subroutine test_semi_gyration(tFinal)
    ! The cosine cap of test_poisson_2d_smooth advanced with the
    ! semi-implicit scheme: forward (alpha = 0), central (1/2) and
    ! backward (1) time weights combined with the backward (upwind),
    ! central and forward (downwind) face values and with the limited
    ! correction, without (tvd) and with (tvd-b) the semi-implicit
    ! bound. Prints the L1 error, the extrema and the change of the
    ! particle number (through the boundary) at two Courant numbers and
    ! saves the solutions at the first one. Then the orders of accuracy
    ! as in Tables 1 and 2 of the semi-implicit manuscript (tFinal = 2*pi
    ! there) and their mean rates. The spatial order: the L1 error on
    ! three grids on the short arc tFinal/16. Over tFinal, the
    ! numerical diffusion of the first order schemes spreads the cap
    ! over more than its width on these grids, and the rate is about
    ! 0.5. The temporal order: the L1 difference from a time-exact
    ! reference on the fixed grid, with 25, 50, 100 and 200 steps to
    ! tFinal/8, also with the correction weighted with alpha = 1/2 by
    ! Picard iterations (tvd-P).

    real, intent(in) :: tFinal

    integer, parameter :: nQ = 100, nGrid = 3, nDt = 4, nStepRef = 800
    integer, parameter :: nQ_I(nGrid) = [50, 100, 200]
    real, parameter :: Cfl_I(2) = [0.5, 4.0]
    character(len=8), parameter :: NameTime_I(3) = &
         ['forward ', 'central ', 'backward']
    character(len=8), parameter :: NameSpace_I(5) = &
         ['backward', 'central ', 'forward ', 'tvd     ', 'tvd-b   ']
    real :: Plot_VC(16,nQ,nQ), Error_I(nGrid), Error, fMax, fMin, dMass
    real :: Sol_VC(nDt,nQ,nQ), Ref_VC(4,nQ,nQ), Cfl
    integer :: iCfl, iTime, iSpace, iGrid, iDt, nStep, iStepFail
    logical :: IsFailed
    character(len=8) :: NameRate
    !--------------------------------------------------------------------------
    do iCfl = 1, size(Cfl_I)
       write(*,'(a,f4.1,a,i4,a)')' test_semi_gyration: CFL =', &
            Cfl_I(iCfl), ',', nQ, '**2 cells'
       write(*,'(a)')' time     space     nStep   L1 error     max f'// &
            '    min f     dN/N'
       do iTime = 1, 3; do iSpace = 1, 5
          call run_gyration(nQ, tFinal, iTime, iSpace, Cfl_I(iCfl), &
               Error, fMax, fMin, dMass, nStep, iStepFail, Plot_VC( &
               5*(iTime-1)+iSpace,:,:), Exact_C=Plot_VC(16,:,:))
          if(iStepFail > 0)then
             write(*,'(1x,a,1x,a,i6,a,i5)')NameTime_I(iTime),        &
                  NameSpace_I(iSpace), nStep, '   blows up at step', &
                  iStepFail
          else
             write(*,'(1x,a,1x,a,i6,es11.3,2f9.4,es10.2)')              &
                  NameTime_I(iTime), NameSpace_I(iSpace), nStep, Error, &
                  fMax, fMin, dMass
          end if
       end do; end do
       if(iCfl == 1) call save_plot_file(                        &
            NameFile='test_semi_gyration.out',                   &
            TypeFileIn='real8', TimeIn=tFinal,                   &
            NameVarIn='-Py Px FB FC FF FT FTb CB CC CF CT CTb'// &
            ' BB BC BF BT BTb Exact Cfl', ParamIn_I=[Cfl_I(1)],  &
            CoordMinIn_D=[-12.0 + 12.0/nQ, -12.0 + 12.0/nQ],     &
            CoordMaxIn_D=[ 12.0 - 12.0/nQ,  12.0 - 12.0/nQ],     &
            VarIn_VII=Plot_VC)
    end do
    ! Spatial order on the arc tFinal/16 at CFL = 0.5. The backward time
    ! weight has Dt ~ DeltaQ**2 (CFL = 25/nQ), so that its first order
    ! temporal error stays below the spatial one
    write(*,'(a,3i4)')' test_semi_gyration: t = tFinal/16, CFL = 0.5'// &
         ' (25/nQ backward), L1 error for nQ =', nQ_I
    do iTime = 1, 3; do iSpace = 1, 5
       IsFailed = .false.
       do iGrid = 1, nGrid
          Cfl = Cfl_I(1)
          if(iTime == 3) Cfl = min(Cfl, 25.0/nQ_I(iGrid))
          call run_gyration(nQ_I(iGrid), tFinal/16, iTime, iSpace, Cfl, &
               Error_I(iGrid), fMax, fMin, dMass, nStep, iStepFail)
          IsFailed = IsFailed .or. iStepFail > 0
       end do
       call write_rate(Error_I)
    end do; end do
    ! Temporal order to tFinal/8 on the fixed grid. As Dt -> 0, all time
    ! weights converge to the same semi-discrete solution, since the
    ! Courant weight of the correction, max(0, 1 - 2*alpha)*CFL, vanishes.
    ! The reference is Crank-Nicolson with nStepRef steps, the correction
    ! weighted with alpha by Picard iterations. The bound of tvd-b,
    ! (1 - alpha)*|f_j - f|, does not depend on Dt, so the bounded
    ! schemes converge to a different solution for every time weight:
    ! tvd-b is left out
    do iSpace = 1, 4
       call run_gyration(nQ, tFinal/8, 2, iSpace, Cfl_I(1), Error, fMax, &
            fMin, dMass, nStep, iStepFail, Ref_VC(iSpace,:,:),          &
            nStepIn=nStepRef, nPicardIn=50)
    end do
    write(*,'(a,i4,a,i4,a)')' test_semi_gyration: t = tFinal/8,', nQ,  &
         '**2 cells, L1 difference from the reference with', nStepRef, &
         ' steps for 25, 50, 100, 200 steps'
    do iTime = 1, 3; do iSpace = 1, 4
       IsFailed = .false.
       do iDt = 1, nDt
          call run_gyration(nQ, tFinal/8, iTime, iSpace, Cfl_I(1), &
               Error, fMax, fMin, dMass, nStep, iStepFail,         &
               Sol_VC(iDt,:,:), nStepIn=25*2**(iDt-1))
          IsFailed = IsFailed .or. iStepFail > 0
       end do
       call write_rate(difference())
    end do; end do
    iTime = 2; iSpace = 4; IsFailed = .false.
    do iDt = 1, nDt
       call run_gyration(nQ, tFinal/8, iTime, iSpace, Cfl_I(1), &
            Error, fMax, fMin, dMass, nStep, iStepFail,         &
            Sol_VC(iDt,:,:), nStepIn=25*2**(iDt-1), nPicardIn=50)
       IsFailed = IsFailed .or. iStepFail > 0
    end do
    call write_rate(difference(), 'tvd-P')
  contains
    !==========================================================================
    function difference() result(Diff_I)
      ! L1 differences between the solutions and the reference
      real :: Diff_I(nDt)
      !------------------------------------------------------------------------
      do iDt = 1, nDt
         Diff_I(iDt) = sum(abs(Sol_VC(iDt,:,:) - Ref_VC(iSpace,:,:))) &
              *(24.0/nQ)**2/36.0
      end do
    end function difference
    !==========================================================================
    subroutine write_rate(Error_I, NameSpaceIn)
      ! Errors and their mean rate of decrease per halving
      real, intent(in) :: Error_I(:)
      character(len=*), optional, intent(in) :: NameSpaceIn
      character(len=8) :: NameSpace
      character(len=40) :: StringFormat
      integer :: n
      !------------------------------------------------------------------------
      NameSpace = NameSpace_I(iSpace)
      if(present(NameSpaceIn)) NameSpace = NameSpaceIn
      n = size(Error_I)
      if(IsFailed .or. any(Error_I <= 0.0))then
         write(*,'(1x,a,1x,a,a)')NameTime_I(iTime), NameSpace, &
              '   blows up'
         RETURN
      end if
      write(NameRate,'(f8.2)') log(Error_I(1)/Error_I(n))/log(2.0)/(n - 1)
      write(StringFormat,'(a,i0,a)') '(1x,a,1x,a,', n, 'es11.3,a,a)'
      write(*,StringFormat)NameTime_I(iTime), NameSpace, Error_I, &
           '  rate', NameRate
    end subroutine write_rate
    !==========================================================================
  end subroutine test_semi_gyration
  !============================================================================
  subroutine run_gyration(nQ, tFinal, iTime, iSpace, Cfl, Error, fMax, &
       fMin, dMass, nStep, iStepFail, Plot_C, nStepIn, nPicardIn, Exact_C)
    ! Advance the cosine cap on nQ**2 cells of [-12,12]**2 to tFinal,
    ! with the time weight 0.5*(iTime-1) and the face scheme iSpace of
    ! test_semi_gyration. The time step is fixed: tFinal is divided into
    ! equal steps at CFL <= Cfl, or into nStepIn steps if present.
    ! iStepFail > 0 is the step at which the solution blows up. The
    ! diagonally dominant matrices (upwind, or central at alpha*CFL < 1)
    ! are solved by Gauss-Seidel iterations, the others by the band LU.
    ! Optionally, nPicardIn iterations weight the limited correction
    ! with alpha. Plot_C and Exact_C are the numerical and the exact
    ! solutions at tFinal.

    integer, intent(in) :: nQ, iTime, iSpace
    real,    intent(in) :: tFinal, Cfl
    real,    intent(out):: Error, fMax, fMin, dMass
    integer, intent(out):: nStep, iStepFail
    real, optional, intent(out):: Plot_C(nQ,nQ), Exact_C(nQ,nQ)
    integer, optional, intent(in):: nStepIn, nPicardIn

    real, parameter :: WidthX = 6.0, WidthY = 6.0
    real, parameter :: Sigma_I(5) = [0.0, 0.5, 1.0, 0.0, 0.0]
    real :: DeltaQ, VDF_G(-1:nQ+2,-1:nQ+2), VDFInit_G(-1:nQ+2,-1:nQ+2)
    real :: VDFFinal_C(nQ,nQ), Source_C(nQ,nQ)
    real :: Volume_G(0:nQ+1,0:nQ+1), Hamiltonian_N(-1:nQ+1,-1:nQ+1)
    real :: Dt, DtCfl, Q, P, Alpha
    integer :: iQ, iP, iStep
    character(len=3) :: TypeSolver
    !--------------------------------------------------------------------------
    DeltaQ = 24.0/nQ
    Volume_G = DeltaQ**2
    do iP = -1, nQ+1; do iQ = -1, nQ+1
       Hamiltonian_N(iQ,iP) = hamiltonian( &
            (DeltaQ*(iQ - nQ/2))**2 + (DeltaQ*(iP - nQ/2))**2)
    end do; end do
    do iP = -1, nQ+2
       P = DeltaQ*(iP - nQ/2 - 0.50)
       do iQ = -1, nQ+2
          Q = DeltaQ*(iQ - nQ/2 - 0.50)
          VDFInit_G(iQ,iP) = initial_cap(Q,P)
          if(min(iQ,iP) >= 1 .and. max(iQ,iP) <= nQ) &
               VDFFinal_C(iQ,iP) = final_cap(Q,P)
       end do
    end do
    VDF_G = VDFInit_G
    ! The time step at CFL = Cfl, reduced to divide tFinal
    call calc_implicit_source(nQ, nQ, VDF_G, Volume_G, Source_C, &
         Hamiltonian_N, AlphaIn=0.0, CFLIn=Cfl, DtOut=DtCfl)
    nStep = ceiling(tFinal/DtCfl)
    if(present(nStepIn)) nStep = nStepIn
    Dt = tFinal/nStep
    Alpha = 0.5*(iTime - 1)
    TypeSolver = 'LU'
    if(Sigma_I(iSpace) == 0.0 .or. (Sigma_I(iSpace) == 0.5 .and. &
         Alpha*Cfl*Dt/DtCfl < 1.0)) TypeSolver = 'lGS'
    iStepFail = 0
    do iStep = 1, nStep
       call calc_implicit_source(nQ, nQ, VDF_G, Volume_G, Source_C,    &
            Hamiltonian_N, AlphaIn=Alpha, SigmaIn=Sigma_I(iSpace), &
            UseTvdIn=iSpace >= 4, UseSemiBoundIn=iSpace == 5,      &
            nPicardIn=nPicardIn, DtIn=Dt, TypeSolverIn=TypeSolver)
       VDF_G(1:nQ,1:nQ) = VDF_G(1:nQ,1:nQ) + Source_C
       if(.not.(maxval(abs(VDF_G(1:nQ,1:nQ))) < 1e3))then
          iStepFail = iStep
          EXIT
       end if
    end do
    Error = sum(abs(VDF_G(1:nQ,1:nQ) - VDFFinal_C) &
         *Volume_G(1:nQ,1:nQ))/(WidthX*WidthY)
    fMax = maxval(VDF_G(1:nQ,1:nQ)); fMin = minval(VDF_G(1:nQ,1:nQ))
    dMass = sum(VDF_G(1:nQ,1:nQ) - VDFInit_G(1:nQ,1:nQ)) &
         /sum(VDFInit_G(1:nQ,1:nQ))
    if(present(Plot_C))then
       where(abs(VDF_G(1:nQ,1:nQ)) < 1e3)
          Plot_C = VDF_G(1:nQ,1:nQ)
       elsewhere
          Plot_C = 1e3
       end where
    end if
    if(present(Exact_C)) Exact_C = VDFFinal_C
  contains
    !==========================================================================
    real function initial_cap(x,y)
      real, intent(in) :: x,y
      real, parameter :: xCenter = -6.0, yCenter = 0.0
      !------------------------------------------------------------------------
      initial_cap = 0.0
      if(abs(x  - xCenter)<=WidthX/2.and.abs(y - yCenter)<=WidthY/2)then
         initial_cap = cos(cPi*(x - xCenter)/WidthX)**4*&
              cos(cPi*(y - yCenter)/WidthY)**4
      end if
    end function initial_cap
    !==========================================================================
    real function final_cap(x,y)
      real, intent(in) :: x,y
      real :: CosPhi,SinPhi,Phi, P, P2, OmegaT, xInit, yInit
      !------------------------------------------------------------------------
      P2  = x**2 + y**2
      OmegaT = tFinal/hamiltonian(P2)
      P  = sqrt(P2)
      CosPhi = x/P; SinPhi = y/P
      if(SinPhi>=0.0)then
         Phi =  acos(CosPhi)
      else
         Phi = -acos(CosPhi)
      end if
      xInit = P*cos(Phi + OmegaT)
      yInit = P*sin(Phi + OmegaT)
      final_cap = initial_cap(xInit,yInit)
    end function final_cap
    !==========================================================================
  end subroutine run_gyration
  !============================================================================
  subroutine test_dsa_semi
    ! Steady shock acceleration of test_dsa_impl, started from the
    ! solution of test_dsa_poisson, with the semi-implicit scheme:
    ! forward, central and backward time weights combined with the
    ! face schemes of test_semi_gyration. The diffusion is implicit.
    ! Prints the number of steps (to a relative change below Tolerance,
    ! or nStepMax, 10 per momentum bin, or to a blow-up), the last
    ! relative change, the minimum of f and the spectral index at
    ! iX = 500, fitted for 0.5 < log10(p) < LogPFitMaxDsa (1.5 on the
    ! nightly grid); saves the spectra. The matrices of the upwind scheme
    ! or of the diffusion alone are solved by Gauss-Seidel iterations
    ! (TypeSolverDsa), the others by the band LU.

    integer :: nStepMax
    real,    parameter :: Tolerance = 1.0e-12
    ! CFL for each time weight: local time steps at the explicit limit
    ! for the forward weight, a global time step otherwise
    real,    parameter :: Cfl_I(3) = [0.98, 98.0, 98.0]
    logical, parameter :: UseLocalDt_I(3) = [.true., .false., .false.]
    character(len=8), parameter :: NameTime_I(3) = &
         ['forward ', 'central ', 'backward']
    character(len=8), parameter :: NameSpace_I(5) = &
         ['backward', 'central ', 'forward ', 'tvd     ', 'tvd-b   ']
    real,    parameter :: Sigma_I(5) = [0.0, 0.5, 1.0, 0.0, 0.0]
    real ::  MomentumRatio, MomentumMin, MomentumMax
    real ::  VDF_G(-1:nX+2, -1:nP+2)
    real ::  Hamiltonian_N(-1:nX+1, -1:nP+1)
    real ::  Volume_G(0:nX+1, 0:nP+1)
    real ::  VolumeX_I(-1:nX+2)
    real ::  VolumeP_I(0:nP+1)
    real ::  DOuter_I(nX) = 1.0, dInner_I(0:nX+1) = 20.0
    real ::  Diffusion_FX(0:nX,1:nP)
    real ::  LogMomentum_I(0:nP+1)     ! Cell centered, for plots
    real ::  Momentum3_I(-1:nP+1)
    real ::  Coord_I(-1:nX+2)  ! Time-dependent coordinate of a mesh
    real ::  Dist_I(-1:nX+1)   ! Distance from mesh i to mesh i+1
    real ::  Source_C(nX,nP), Plot_VI(15,nP), Change, SpectralIndex
    logical :: IsFit_I(nP), IsFailed
    integer :: iP, iStep, iTime, iSpace
    character(len=12) :: NameIndex
    character(len=5)  :: TypeSolver
    !--------------------------------------------------------------------------
    ! Logarithmic grid in momentum, as in test_dsa_impl
    MomentumRatio   = exp(log(pMaxDsa/pMinDsa)/nP)
    MomentumMax     = pMinDsa/MomentumRatio
    Momentum3_I(-1) = MomentumMax**3/3
    do iP = 0, nP + 1
       MomentumMin = MomentumMax
       MomentumMax = MomentumMin*MomentumRatio
       Momentum3_I(iP) = MomentumMax**3/3
       VolumeP_I(iP)   = Momentum3_I(iP) - Momentum3_I(iP-1)
       LogMomentum_I(iP) = 0.50*log10(MomentumMin*MomentumMax)
    end do
    IsFit_I = LogMomentum_I(1:nP) > 0.5 .and. &
         LogMomentum_I(1:nP) < LogPFitMaxDsa
    nStepMax = 10*nP
    call update_coords(600.0)
    write(*,'(a)')' test_dsa_semi: time     space     CFL   nStep'// &
         '  change    min f    index'
    do iTime = 1, 3; do iSpace = 1, 5
       VDF_G = VDFOld_G
       TypeSolver = 'LU'
       if(Sigma_I(iSpace) == 0.0 .or. iTime == 1) TypeSolver = TypeSolverDsa
       do iStep = 1, nStepMax
          if(IsTransposedDsa)then
             call step_transposed
          else
             call calc_implicit_source(nX, nP, VDF_G, Volume_G, Source_C,    &
                  Hamiltonian_N, AlphaIn=0.5*(iTime - 1),                &
                  SigmaIn=Sigma_I(iSpace), UseTvdIn=iSpace >= 4,         &
                  UseSemiBoundIn=iSpace == 5,                            &
                  Diffusion_FX=Diffusion_FX, AlphaDiffIn=1.0,            &
                  IsSteadyState=UseLocalDt_I(iTime), CFLIn=Cfl_I(iTime), &
                  TypeSolverIn=TypeSolver, nSweepMaxIn=nSweepDsa)
          end if
          VDF_G(1:nX,1:nP) = VDF_G(1:nX,1:nP) + Source_C
          Change = maxval(abs(Source_C)/(abs(VDF_G(1:nX,1:nP)) + 1e-30))
          VDF_G(1:nX,-1:0 ) = FloorDsa
          VDF_G(1:nX,nP+1:nP+2) = FloorDsa
          VDF_G( 0,      :) = VDF_G(1,       :)
          VDF_G(-1,      :) = VDF_G(1,       :)
          VDF_G(nX+1:nX+2,iPInj) = vdf_inj(Momentum3_I(iPInj-1:iPInj))
          IsFailed = .not.(maxval(abs(VDF_G(1:nX,1:nP))) < 1e6)
          if(Change < Tolerance .or. IsFailed) EXIT
       end do
       iStep = min(iStep, nStepMax)
       ! Least squares slope of log f versus log p
       if(all(VDF_G(500,1:nP) > 0.0 .or. .not.IsFit_I) .and. &
            .not.IsFailed)then
          SpectralIndex = slope(LogMomentum_I(1:nP), &
               log10(abs(VDF_G(500,1:nP))))
          write(NameIndex,'(f9.3)') SpectralIndex
       else
          NameIndex = '   none'
       end if
       if(IsFailed)then
          write(*,'(15x,a,1x,a,f6.2,i7,a)')NameTime_I(iTime), &
               NameSpace_I(iSpace), Cfl_I(iTime), iStep, '   blows up'
       else
          write(*,'(15x,a,1x,a,f6.2,i7,es10.2,es10.2,a)')            &
               NameTime_I(iTime), NameSpace_I(iSpace), Cfl_I(iTime), &
               iStep, Change, minval(VDF_G(1:nX,1:nP)), trim(NameIndex)
       end if
       Plot_VI(5*(iTime-1)+iSpace,:) = &
            log10(max(abs(VDF_G(500,1:nP)), 1e-30))
    end do; end do
    call save_plot_file(NameFile='test_dsa_semi.out',              &
         TypeFileIn='ascii', nStepIn = nStepMax,                   &
         NameVarIn='LogMomentum FB FC FF FT FTb CB CC CF CT CTb'// &
         ' BB BC BF BT BTb',                                       &
         CoordMinIn_D=[ LogMomentum_I(1 ) ],                       &
         CoordMaxIn_D=[ LogMomentum_I(nP) ],                       &
         Coord1In_I = LogMomentum_I(1:nP),                         &
         VarIn_VI = Plot_VI)
  contains
    !==========================================================================
    real function slope(x_I, y_I)
      real, intent(in) :: x_I(nP), y_I(nP)
      real :: xMean, yMean
      !------------------------------------------------------------------------
      xMean = sum(x_I, MASK=IsFit_I)/count(IsFit_I)
      yMean = sum(y_I, MASK=IsFit_I)/count(IsFit_I)
      slope = sum((x_I - xMean)*(y_I - yMean), MASK=IsFit_I) &
           /sum((x_I - xMean)**2, MASK=IsFit_I)
    end function slope
    !==========================================================================
    subroutine step_transposed
      ! The step with the coordinates in the other order, the lines of
      ! the solver along the momentum; the Poisson bracket changes sign
      real, allocatable :: F_G(:,:), H_N(:,:), V_G(:,:), D_FY(:,:), &
           S_C(:,:)
      !------------------------------------------------------------------------
      allocate(F_G(-1:nP+2,-1:nX+2), H_N(-1:nP+1,-1:nX+1), &
           V_G(0:nP+1,0:nX+1), D_FY(1:nP,0:nX), S_C(nP,nX))
      F_G = transpose(VDF_G); H_N = -transpose(Hamiltonian_N)
      V_G = transpose(Volume_G); D_FY = transpose(Diffusion_FX)
      call calc_implicit_source(nP, nX, F_G, V_G, S_C, H_N,             &
           AlphaIn=0.5*(iTime - 1), SigmaIn=Sigma_I(iSpace),        &
           UseTvdIn=iSpace >= 4, UseSemiBoundIn=iSpace == 5,        &
           Diffusion_FY=D_FY, AlphaDiffIn=1.0,                      &
           IsSteadyState=UseLocalDt_I(iTime), CFLIn=Cfl_I(iTime),   &
           TypeSolverIn=TypeSolver, nSweepMaxIn=nSweepDsa)
      Source_C = transpose(S_C)
    end subroutine step_transposed
    !==========================================================================
    subroutine update_coords(Time)

      real, intent(in):: Time

      ! Effective width of the shock wave, in cell size
      integer, parameter:: nWidth = 1

      ! terms in the functions, describing the
      ! lagrangian point motion
      real, parameter:: cQuad = 0.375/nWidth, cConst = 0.375*nWidth
      ! Loop variables
      integer :: iX, iP
      real    :: CoordLagr, DiffusionX_I(0:nX)
      ! Shock wave with the compression ratio of 4, with width ~1
      !------------------------------------------------------------------------
      ! Coord_I is the cell-centered coordinate
      do iX = -1, nX +2
         CoordLagr = iX - 0.50 -Time
         if(CoordLagr>=0.0)then
            Coord_I(iX) = CoordLagr + Time
         elseif(CoordLagr>=-real(nWidth))then
            Coord_I(iX) = CoordLagr + Time + cQuad*CoordLagr**2
         else
            Coord_I(iX) = 0.250*CoordLagr -cConst + Time
         end if
      end do
      ! Dist is the distance to the next cell center
      Dist_I(-1:nX+1) = Coord_I(0:nX+2) - Coord_I(-1:nX+1)
      VolumeX_I(0:nX+1) = 0.50*(Dist_I(-1:nX) + Dist_I(0:nX+1))
      VolumeX_I(-1)  = Dist_I(-1); VolumeX_I(nX+2) =  Dist_I(nX+1)
      do iP = 0, nP+1
         Volume_G(0:nX+1,iP) = VolumeP_I(iP)*VolumeX_I(0:nX+1)
      end do
      do iP = -1, nP+1
         Hamiltonian_N(-1:nX+1, iP)= -0.50*Momentum3_I(iP)*&
              (VolumeX_I(-1:nX+1) + VolumeX_I(0:nX+2))
      end do
      ! Get diffusion coefficients
      ! No-flux BC at the left boundary:
      DiffusionX_I(0) = 0.0
      do iX = 1, nX
         DiffusionX_I(iX) = &
              DOuter_I(iX)*0.50*(DInner_I(iX) + DInner_I(iX+1))/Dist_I(iX)
      end do
      do iP = 1, nP
         Diffusion_FX(:,iP) =  VolumeP_I(iP)*DiffusionX_I
      end do
    end subroutine update_coords
    !==========================================================================
  end subroutine test_dsa_semi
  !============================================================================
  subroutine test_semi_implicit_one_rotation(tFinal)
    ! Order of accuracy of the first order (upwind) scheme,
    ! semi_implicit_one, with the explicit (alpha = 0), Crank-Nicolson (1/2)
    ! and implicit (1) time weights, on the rotation of run_rotation.
    ! Expected: spatial order (CFL = 0.5) about 1 for all three; temporal
    ! order 1, 2 and 1.

    real, intent(in) :: tFinal

    integer :: iAlpha
    character(len=*), parameter:: NameSub = 'test_semi_implicit_one_rotation'
    !--------------------------------------------------------------------------
    write(*,*) '=== ', NameSub, ': SPATIAL order (CFL=0.5) ==='
    write(*,*) '   (first-order upwind -> order ~ 1 for all three schemes)'
    do iAlpha = 1, 3
       write(*,'(a)') ' '//NameAlpha_I(iAlpha)
       call write_spatial_order('semi_implicit_one', 0.5*(iAlpha - 1), &
            [40, 80, 160, 320], tFinal, Cfl=0.5)
    end do
    write(*,*) '=== ', NameSub, ': TEMPORAL order (nC=100, ref=tiny dt) ==='
    write(*,*) '   (explicit -> 1, implicit -> 1, Crank-Nicolson -> 2)'
    call write_temporal_order('semi_implicit_one', 100, 0.02, &
         [20, 40, 80, 160], AlphaRef=0.5)
    ! The finest Crank-Nicolson solution and the exact rotated blob
    call save_rotation('semi_implicit_one', 0.5, .true., 0.5, tFinal)
  end subroutine test_semi_implicit_one_rotation
  !============================================================================
  subroutine test_semi_implicit_two_rotation(tFinal)
    ! Order of accuracy of the second order scheme without limiter,
    ! semi_implicit_two, on the rotation of run_rotation. Expected: spatial
    ! order about 2 with Crank-Nicolson (CFL = 0.1); temporal order 2 for
    ! Crank-Nicolson, 1 for the implicit scheme.

    real, intent(in) :: tFinal

    character(len=*), parameter:: NameSub = 'test_semi_implicit_two_rotation'
    !--------------------------------------------------------------------------
    write(*,*) '=== ', NameSub, ': SPATIAL order (Crank-Nicolson, CFL=0.1) ==='
    write(*,*) '   (centered 2nd-order operator + CN time -> order ~ 2)'
    call write_spatial_order('semi_implicit_two', 0.5, [40, 80, 160, 320], &
         tFinal, Cfl=0.1)
    write(*,*) '=== ', NameSub, ': TEMPORAL order (nC=100, ref=tiny dt) ==='
    write(*,*) '   (Crank-Nicolson -> 2, backward Euler alpha=1 -> 1)'
    call write_temporal_order('semi_implicit_two', 100, 0.02, &
         [20, 40, 80, 160], AlphaRef=0.5)
    ! The finest Crank-Nicolson solution and the exact rotated blob
    call save_rotation('semi_implicit_two', 0.5, .true., 0.5, tFinal)
  end subroutine test_semi_implicit_two_rotation
  !============================================================================
  subroutine test_semi_implicit_tvd_rotation(tFinal)
    ! Order of accuracy of the TVD scheme, semi_implicit_tvd, on the
    ! rotation of run_rotation. On the smooth blob the limiters are
    ! mostly inactive away from the peak. The spatial order is measured
    ! with alpha = 1 and Dt = CsDt*Delta**2, so that the first order
    ! temporal error does not mask the spatial one. At alpha = 1 the
    ! bound of the correction at the extrema, (1 - alpha)*|f_j - f|, is
    ! zero, which switches the correction off at the extrema along the
    ! trajectories: the order is about 1 (about 2 without the bound).
    ! With the original CsDt = 0.5 the Picard iterations do not converge
    ! on 50**2 cells (CFL = 2.2), and 400**2 cells take more than 10
    ! minutes. The bound does not vanish as Dt -> 0, so the schemes with
    ! different alpha converge to different solutions: each is compared
    ! with its own reference.

    real, intent(in) :: tFinal

    character(len=*), parameter:: NameSub = &
         'test_semi_implicit_tvd_rotation'
    !--------------------------------------------------------------------------
    write(*,*) '=== ', NameSub, ': SPATIAL order (alpha=1, Dt=CsDt*dC^2) ==='
    write(*,*) '   (the bound at the extrema is zero at alpha=1 -> order ~ 1)'
    ! The finest solution and the exact rotated blob are saved
    call write_spatial_order('semi_implicit_tvd', 1.0, [50, 100, 200], &
         tFinal, CsDt=0.25, DoSave=.true.)
    write(*,*) '=== ', NameSub, ': TEMPORAL order (nC=50, ref=tiny dt) ==='
    write(*,*) '   (the reference has the same alpha)'
    call write_temporal_order('semi_implicit_tvd', 50, 0.01, &
         [40, 80, 160, 320])
  end subroutine test_semi_implicit_tvd_rotation
  !============================================================================
  subroutine write_spatial_order(NameScheme, Alpha, nC_I, tFinal, Cfl, &
       CsDt, DoSave)
    ! L1 errors of run_rotation on the grids nC_I and their orders. The
    ! time step is from the Courant number Cfl or, if CsDt is present,
    ! CsDt times the cell size squared. With DoSave, the solution on the
    ! finest grid is saved by write_rotation
    character(len=*), intent(in) :: NameScheme
    real,    intent(in) :: Alpha, tFinal
    integer, intent(in) :: nC_I(:)
    real,    optional, intent(in) :: Cfl, CsDt
    logical, optional, intent(in) :: DoSave

    real, allocatable :: Plot_VC(:,:,:)
    real    :: Error, ErrorPrev, CflOrDt
    integer :: iGrid, nC
    !--------------------------------------------------------------------------
    ErrorPrev = -1.0
    do iGrid = 1, size(nC_I)
       nC = nC_I(iGrid)
       if(present(CsDt))then
          CflOrDt = CsDt*(WidthRotation/nC)**2
       else
          CflOrDt = Cfl
       end if
       if(iGrid == size(nC_I) .and. present(DoSave))then
          allocate(Plot_VC(2,nC,nC))
          call run_rotation(NameScheme, nC, Alpha, .not.present(CsDt), &
               CflOrDt, tFinal, Error, fSave_C=Plot_VC(1,:,:),         &
               Exact_C=Plot_VC(2,:,:))
          if(DoSave) call write_rotation(NameScheme, tFinal, Plot_VC)
       else
          call run_rotation(NameScheme, nC, Alpha, .not.present(CsDt), &
               CflOrDt, tFinal, Error)
       end if
       call write_error('nC', nC, Error, ErrorPrev)
       ErrorPrev = Error
    end do
  end subroutine write_spatial_order
  !============================================================================
  subroutine write_temporal_order(NameScheme, nC, tFinal, nStep_I, &
       AlphaRef)
    ! L1 differences between the solutions of run_rotation on nC**2 cells
    ! with nStep_I steps and the reference with 1280 steps, for the time
    ! weights 0, 1/2 and 1. The reference has the time weight AlphaRef
    ! if present, else the same one
    character(len=*), intent(in) :: NameScheme
    integer, intent(in) :: nC, nStep_I(:)
    real,    intent(in) :: tFinal
    real, optional, intent(in) :: AlphaRef

    real    :: Ref_C(nC,nC), Error, ErrorPrev
    integer :: iAlpha, iStep
    !--------------------------------------------------------------------------
    if(present(AlphaRef)) call run_rotation(NameScheme, nC, AlphaRef, &
         .false., tFinal/1280.0, tFinal, Error, fSave_C=Ref_C)
    do iAlpha = 1, 3
       write(*,'(a)') ' '//NameAlpha_I(iAlpha)
       if(.not.present(AlphaRef)) call run_rotation(NameScheme, nC, &
            0.5*(iAlpha - 1), .false., tFinal/1280.0, tFinal, Error, &
            fSave_C=Ref_C)
       ErrorPrev = -1.0
       do iStep = 1, size(nStep_I)
          call run_rotation(NameScheme, nC, 0.5*(iAlpha - 1), .false., &
               tFinal/nStep_I(iStep), tFinal, Error, fRef_C=Ref_C)
          call write_error('nStep', nStep_I(iStep), Error, ErrorPrev)
          ErrorPrev = Error
       end do
    end do
  end subroutine write_temporal_order
  !============================================================================
  subroutine write_error(NameN, n, Error, ErrorPrev)
    ! A line of the tables: the order from the previous error, if any
    character(len=*), intent(in) :: NameN
    integer, intent(in) :: n
    real,    intent(in) :: Error, ErrorPrev
    !--------------------------------------------------------------------------
    if(ErrorPrev > 0.0)then
       write(*,'(a,i5,a,es11.3,a,f6.2)') '    '//NameN//'=', n, &
            '  L1=', Error, '   order=', log(ErrorPrev/Error)/log(2.0)
    else
       write(*,'(a,i5,a,es11.3)') '    '//NameN//'=', n, '  L1=', Error
    end if
  end subroutine write_error
  !============================================================================
  subroutine save_rotation(NameScheme, Alpha, UseCfl, CflOrDt, tFinal)
    ! Save the solution of run_rotation on 200**2 cells with
    ! write_rotation
    character(len=*), intent(in) :: NameScheme
    real,    intent(in) :: Alpha, CflOrDt, tFinal
    logical, intent(in) :: UseCfl

    integer, parameter :: nC = 200
    real :: Plot_VC(2,nC,nC), Error
    !--------------------------------------------------------------------------
    call run_rotation(NameScheme, nC, Alpha, UseCfl, CflOrDt, tFinal, &
         Error, fSave_C=Plot_VC(1,:,:), Exact_C=Plot_VC(2,:,:))
    call write_rotation(NameScheme, tFinal, Plot_VC)
  end subroutine save_rotation
  !============================================================================
  subroutine write_rotation(NameScheme, tFinal, Plot_VC)
    ! Write the numerical and the exact solutions of run_rotation to
    ! test_NameScheme_rotation.out
    character(len=*), intent(in) :: NameScheme
    real, intent(in) :: tFinal, Plot_VC(:,:,:)

    real :: dC
    !--------------------------------------------------------------------------
    dC = WidthRotation/size(Plot_VC, 2)
    call save_plot_file(NameFile='test_'//NameScheme//'_rotation.out', &
         TypeFileIn='real8', TimeIn=tFinal, nStepIn=0,                &
         NameVarIn='x y VDF Exact',                                   &
         CoordMinIn_D=[-0.5*WidthRotation + 0.5*dC,                   &
         -0.5*WidthRotation + 0.5*dC],                                &
         CoordMaxIn_D=[0.5*WidthRotation - 0.5*dC,                    &
         0.5*WidthRotation - 0.5*dC], VarIn_VII=Plot_VC)
  end subroutine write_rotation
  !============================================================================
  subroutine run_rotation(NameScheme, nC, Alpha, UseCfl, CflOrDt, tFinal, &
       Error, fSave_C, fRef_C, Exact_C)
    ! E x B rigid rotation of a Gaussian blob, the test of the schemes
    ! of the semi-implicit manuscript: the Hamiltonian is
    ! (Omega/2)*(x**2 + y**2), the period 2*pi/Omega = 1, the domain
    ! [-3,3]**2 with nC**2 cells, and the exact solution is the initial
    ! blob rotated rigidly. Advance it to tFinal with the scheme
    ! NameScheme (semi_implicit_one, semi_implicit_two or semi_implicit_tvd)
    ! and the time weight Alpha. With UseCfl the time step is from
    ! CFL = CflOrDt, the last one reduced to reach tFinal; otherwise it is
    ! CflOrDt, adjusted to divide tFinal. Error is the relative L1 error
    ! with respect to fRef_C if present, else to the exact solution.
    ! fSave_C and Exact_C are the numerical and exact solutions.

    character(len=*), intent(in) :: NameScheme
    integer, intent(in) :: nC
    real,    intent(in) :: Alpha, CflOrDt, tFinal
    logical, intent(in) :: UseCfl
    real,    intent(out):: Error
    real, optional, intent(out):: fSave_C(nC,nC), Exact_C(nC,nC)
    real, optional, intent(in) :: fRef_C(nC,nC)

    real, parameter :: Omega = cTwoPi
    ! Initial center and width of the blob
    real, parameter :: xBlob = -0.6, yBlob = 0.6, WidthBlob = 0.4
    real, allocatable :: VDF_G(:,:), Volume_G(:,:), Hamiltonian_N(:,:), &
         Source_C(:,:), Exact_II(:,:)
    real    :: dC, Time, Dt, DtLocal
    integer :: i, j, iStep, nStep
    !--------------------------------------------------------------------------
    dC = WidthRotation/nC
    allocate(VDF_G(-1:nC+2,-1:nC+2), Volume_G(0:nC+1,0:nC+1), &
         Hamiltonian_N(-1:nC+1,-1:nC+1), Source_C(nC,nC),    &
         Exact_II(nC,nC))
    Volume_G = dC*dC
    do j = -1, nC+1; do i = -1, nC+1
       Hamiltonian_N(i,j) = 0.5*Omega* &
            ((dC*(i - nC/2))**2 + (dC*(j - nC/2))**2)
    end do; end do
    do j = -1, nC+2; do i = -1, nC+2
       VDF_G(i,j) = blob(dC*(i - nC/2) - 0.5*dC, dC*(j - nC/2) - 0.5*dC, &
            0.0)
    end do; end do
    if(UseCfl)then
       Time = 0.0
       do
          call advance(CFLIn=CflOrDt, DtOut=Dt)
          if(Time + Dt >= tFinal)then
             Dt = tFinal - Time
             call advance(DtIn=Dt)
             VDF_G(1:nC,1:nC) = VDF_G(1:nC,1:nC) + Source_C
             EXIT
          end if
          VDF_G(1:nC,1:nC) = VDF_G(1:nC,1:nC) + Source_C
          Time = Time + Dt
       end do
    else
       nStep = max(1, nint(tFinal/CflOrDt))
       Dt = tFinal/nStep
       do iStep = 1, nStep
          DtLocal = Dt
          call advance(DtIn=DtLocal)
          VDF_G(1:nC,1:nC) = VDF_G(1:nC,1:nC) + Source_C
       end do
    end if
    do j = 1, nC; do i = 1, nC
       Exact_II(i,j) = blob(dC*(i - nC/2) - 0.5*dC, &
            dC*(j - nC/2) - 0.5*dC, tFinal)
    end do; end do
    if(present(fRef_C))then
       Error = sum(abs(VDF_G(1:nC,1:nC) - fRef_C))/sum(abs(fRef_C))
    else
       Error = sum(abs(VDF_G(1:nC,1:nC) - Exact_II))/sum(abs(Exact_II))
    end if
    if(present(fSave_C)) fSave_C = VDF_G(1:nC,1:nC)
    if(present(Exact_C)) Exact_C = Exact_II
  contains
    !==========================================================================
    real function blob(x, y, Time)
      ! The blob with its center rotated by the angle Omega*Time:
      ! dx/dt = Omega*y, dy/dt = -Omega*x
      real, intent(in) :: x, y, Time
      real :: xCenter, yCenter
      !------------------------------------------------------------------------
      xCenter =  xBlob*cos(Omega*Time) + yBlob*sin(Omega*Time)
      yCenter = -xBlob*sin(Omega*Time) + yBlob*cos(Omega*Time)
      blob = exp(-((x - xCenter)**2 + (y - yCenter)**2)/(2*WidthBlob**2))
    end function blob
    !==========================================================================
    subroutine advance(DtIn, CFLIn, DtOut)
      real, optional, intent(inout) :: DtIn
      real, optional, intent(in)    :: CFLIn
      real, optional, intent(out)   :: DtOut
      !------------------------------------------------------------------------
      select case(NameScheme)
      case('semi_implicit_one')
         call semi_implicit_one(nC, nC, VDF_G, Volume_G, Source_C, Alpha, &
              Hamiltonian12_N=Hamiltonian_N, DtIn=DtIn, CFLIn=CFLIn,  &
              DtOut=DtOut)
      case('semi_implicit_two')
         call semi_implicit_two(nC, nC, VDF_G, Volume_G, Source_C, Alpha, &
              Hamiltonian12_N=Hamiltonian_N, DtIn=DtIn, CFLIn=CFLIn,  &
              DtOut=DtOut)
      case('semi_implicit_tvd')
         call semi_implicit_tvd(nC, nC, VDF_G, Volume_G, Source_C, Alpha, &
              Hamiltonian12_N=Hamiltonian_N, DtIn=DtIn, CFLIn=CFLIn,     &
              DtOut=DtOut)
      case default
         call CON_stop('run_rotation: unknown scheme '//NameScheme)
      end select
    end subroutine advance
    !==========================================================================
  end subroutine run_rotation
  !============================================================================
end module ModTestPoissonBracket
!==============================================================================
module ModTestMultiPoisson
  ! Solve energetic particle transport (focused transport equation) via the
  ! numerical particle-conserved scheme with multi Poisson bracket:
  !   {f,H_1}_t,p^3/3 + {f,H_2}_s_L,\mu + {f,H_3}_p^3/3,\mu = 0
  ! Dimensionless parameters: Energy: GeV, Speed: c Mass: GeV/c^2 Length: c*sec
  ! In the code, label Q: s_L  P: \mu  R: ln(P^3/3)  NLoop: t
  ! this code consider only the most simple equation
  use ModPoissonBracket, ONLY: explicit
  use ModUtilities, ONLY: CON_stop
  use ModNumConst, ONLY: cTwoPi, cPi
  use ModPlotFile, ONLY: save_plot_file
  use ModConst, ONLY: cRmeProton, cGeV, &
       cMu, cAtomicMass, rSun, cLightSpeed
  implicit none
  ! Grid setups
  integer, parameter :: nQ = 160   ! Number of grid of s_L
  integer, parameter :: nP = 20    ! Number of grid of \mu = cos(pitch angle)
  integer, parameter :: nR = 60    ! Number of grid of ln(p^3/3) axis
  ! Unit setups
  real, parameter :: cRmeProtonGeV = cRmeProton/cGeV ! Proton mass: GeV/c^2
  real, parameter :: CoeffMuToxx = 27.0/14.0         ! DMuMu => Dxx constant
  ! Physical parameter reading setups
  real, parameter :: tEachFile = 100.0    ! Time range of each input file
  ! Number of variables in each line of the input files
  integer, parameter :: nVar = 14
  ! Position of the variables in the input line
  integer, parameter :: x_ = 2, y_ = 3, z_ = 4,    &
       rho_ = 5, ux_ = 7, uy_ = 8, uz_ = 9,        &
       Bx_ = 10, By_= 11, Bz_ = 12, Wave1_ = 13, Wave2_ = 14
  ! OLD and NEW variables from raw data files
  ! Raw data imported from MHD files
  real :: RawData1_II(nQ, nVar), RawData2_II(nQ, nVar)
  ! \Deltas/b, 1/(2B), ln(B\deltas^2) at Old time
  real :: DeltaSOverBOld_C(nQ), InvBOld_C(nQ), LnBDeltaS2Old_C(nQ)
  ! \Deltas/b, 1/(2B), ln(B\deltas^2) at New time
  real :: DeltaSOverBNew_C(nQ), InvBNew_C(nQ), LnBDeltaS2New_C(nQ)
  ! Midpoint for to consecutive points, \deltas
  real :: MidPoint_ID(nQ-1, x_:z_), DeltaS_I(nQ), DeltaSface_I(1:nQ-1)

contains
  !============================================================================
  subroutine test_multi_poisson(tFinal, TypeScatter)
    ! Test for multiple Poisson bracket

    real, intent(in) :: tFinal            ! Overall final Time
    ! No scatter, or advect and scatter, or only diffusion scatter
    character(LEN=*), intent(in) :: TypeScatter

    ! Local VARs:
    real, parameter  :: cMinVal = 1.0e-8  ! We set Min Value of VDF as 1.0e-8
    ! ------------ File ------------
    ! Unit for field line files: Old and new field lines
    integer, parameter :: nInFileOld = 8, nInFileNew = 9
    ! Stride of line of the data reading
    integer, parameter :: nStride = 25
    ! Lagr index that is gaped for field line at each time
    integer, parameter :: iLagrIdSkip = 20
    ! Read the redundant data
    real               :: GapData
    ! Heliocentric distance, for visualization
    real               :: HelioDist_I(nQ)
    ! Number of outputs:
    integer, parameter :: nOutputFile = 100
    ! The number index of output files in character
    character(LEN=3)   :: NameFileSuffix
    ! The data calculated directly from RawData (MHD data)
    real :: DeltaSOverB_C(nQ), dDeltaSOverBDt_C(nQ),  &
         InvB_C(nQ), bDuDt_C(nQ), dLnBDeltaS2Dt_C(nQ)

    ! ------------ Pitch angle ------------
    ! Now we have the parameters for "Q"-coordinate (s_L)
    ! We then consider the "P"-coordinate (mu)
    ! It is the cosine value of the pitch angle, varying from -1 to 1
    ! Here we calculate \delta\mu for each grid along the \mu axis
    real :: DeltaMu = 2.0/nP

    ! ------------ Energy and momentum ------------
    ! Last but not least, we focus on the "R"-coordinate (p^3/3)
    ! Min and Max value of particle kinetic energy, in units of GeV
    real, parameter :: EkMin = 1.0e-3, EkMax = 1.0
    ! Ln(P^3/3)_min, Ln(P^3/3)_max, \deltaLn(P^3/3)
    real            :: LnP3min, LnP3max, DeltaLnP3
    ! Momentum and velocity at the center of each cell
    real            :: Momentum_I(nR),  Velocity_I(nR)
    ! Momentum at the face of each cell and their distance \delta(P^3/3)
    real            :: MomentumFace_I(0:nR), DeltaP3_I(nR)
    ! Energy in unit of KeV, for output
    real            :: Energy_I(nR)

    ! ------------ Hamiltonian functions ------------
    ! Poisson bracket with regard to the first and second vars
    ! considering the case when there are more than one Poisson bracket
    ! and when there is a Poisson bracket with respect to the time
    real :: DeltaHamiltonian1_N(0:nQ+1, 0:nP+1, -1:nR+1)
    real :: Hamiltonian2_N(-1:nQ+1, -1:nP+1, 0:nR+1)
    real :: Hamiltonian3_N(0:nQ+1, -1:nP+1, -1:nR+1)

    ! ------------ Control volume ------------
    real :: Volume_G(0:nQ+1, 0:nP+1, 0:nR+1)    ! Control volume (phase)
    real :: dVolumeDt_G(0:nQ+1, 0:nP+1, 0:nR+1) ! Derivative regarding to time

    ! ------------ Distribution functions ------------
    real :: VDF_G(-1:nQ+2, -1:nP+2, -1:nR+2) ! Velocity distribution function
    real :: VDFOutput_II(nQ, nR)             ! VDF outputs
    real :: Source_C(nQ, nP, nR)             ! Source at each time step
    ! Parameters of the place of initial particles, for initial VDF_G in test
    real, parameter :: ParticleRange1 = 2.0, ParticleRange2 = 4.0

    ! ------------ Diffusion ------------
    real :: LambdaMuMu_II(nQ, nR)            ! lambda_\mu\mu from raw data

    ! ------------ Time variables ------------
    real            :: Time, Dt, DtNext      ! Simulation time, time step
    real            :: tOutput               ! Time of each output file
    logical         :: DoExit                ! Exit time marching
    real            :: Cfl = 0.99            ! CFL number

    ! ------------ Loop variables ------------
    integer :: iQ, iP, iR           ! For VDF variables
    integer :: iStep                ! For time marching
    integer :: iStride, iVar        ! For reading files
    integer :: iOutputFile          ! For output files
    !--------------------------------------------------------------------------

    call read_fieldline             ! Read MhData from saved files
    ! Calculate the Lagrangian particle positions
    HelioDist_I = sqrt(sum(RawData1_II(:, x_:z_)**2, DIM=2))*cLightSpeed/rSun

    ! Kinetic energy convert to momentum (1 MeV – 1 GeV) => (Pmin, Pmax)
    LnP3min = log(sqrt((EkMin + cRmeProtonGeV)**2 - cRmeProtonGeV**2)**3.0/3)
    LnP3max = log(sqrt((EkMax + cRmeProtonGeV)**2 - cRmeProtonGeV**2)**3.0/3)
    ! Get \DeltaLn(P^3/3)
    DeltaLnP3 = (LnP3max - LnP3min)/nR

    ! Then we can culculate P and deltaP^3/3
    MomentumFace_I(0) = (3.0*exp(LnP3min))**(1.0/3.0)
    do iR = 1, nR
       ! The momentum is calculated at each center of the cell
       Momentum_I(iR) = (3.0*exp(LnP3min + (iR-0.5)*DeltaLnP3))**(1.0/3.0)
       ! The momentum is calculated at the face of the cell
       MomentumFace_I(iR) = (3.0*exp(LnP3min + iR*DeltaLnP3))**(1.0/3.0)
    end do

    ! Calculate \deltaP^3/3 for each grid: Face - Face, like the Volume_P
    DeltaP3_I = MomentumFace_I(1:nR)**3.0/3 - MomentumFace_I(0:nR-1)**3.0/3
    ! Calculate the position for ln(p^3/3) axis, in the unit of KeV now
    Energy_I = 6.0+log10(sqrt(Momentum_I**2+cRmeProtonGeV**2) - cRmeProtonGeV)
    ! Considering the law of relativity, v=1/sqrt(1+m^2*c^2/p^2), v can be
    ! calculated as a function of p. Note that light speed is the unit of
    ! speed here, so we do not need to multiply c^2 in the following steps
    Velocity_I = 1.0/sqrt(1.0 + cRmeProtonGeV**2/Momentum_I**2)

    call init_test_VDF           ! Initialize the VDF
    ! Calculate VDF for output: integrate over the \mu axis
    VDFOutput_II = sum(VDF_G(1:nQ, 1:nP, 1:nR), dim=2)*DeltaMu
    ! Save the initial outpur of VDF
    call save_plot_file(NameFile='test_multipoisson_000.out',  &
         TypeFileIn = 'ascii', TimeIn = Time, nStepIn = iStep, &
         NameVarIn = 'sL logE VDF',                            &
         CoordMinIn_D = [HelioDist_I(1) ,    Energy_I(1)],     &
         CoordMaxIn_D = [HelioDist_I(nQ),    Energy_I(nR)],    &
         StringFormatIn = '(3F15.7)',                          &
         Coord1In_I = HelioDist_I,  Coord2In_I = Energy_I,     &
         VarIn_II = log10(VDFOutput_II(1:nQ, 1:nR)))

    ! Start the main simulation loop
    iStep = 0; Time = 0.0
    call calc_data_states
    call update_states(nQ, Time, dDeltaSOverBDt_C = dDeltaSOverBDt_C)
    ! Calculate the total control volume
    do iR = 1, nR
       do iQ = 1, nQ
          dVolumeDt_G(iQ, 0:nP+1, iR) = &
               dDeltaSOverBDt_C(iQ)*DeltaMu*DeltaP3_I(iR)
       end do
    end do
    ! Boundary conditions for total control volume time derivative
    dVolumeDt_G(0,    :,    :) = dVolumeDt_G(1 , :,  :)
    dVolumeDt_G(nQ+1, :,    :) = dVolumeDt_G(nQ, :,  :)
    dVolumeDt_G(:,    :,    0) = dVolumeDt_G(: , :,  1)
    dVolumeDt_G(:,    :, nR+1) = dVolumeDt_G(: , :, nR)

    ! Set up first trial to get DtNext
    call advance_advect(nP, nQ, nR, DeltaMu,                      &
         MomentumFace_I, dDeltaSOverBDt_C, DeltaHamiltonian1_N,   &
         DeltaP3_I, Velocity_I, InvB_C, Hamiltonian2_N,           &
         DeltaSOverB_C, dLnBDeltaS2Dt_C, bDuDt_C, Hamiltonian3_N, &
         iStep, Time, Dt, DtNext, Cfl,                            &
         Volume_G, dVolumeDt_G, VDF_G, Source_C)

    ! In the loop, we output several snapshots of VDF
    do iOutputFile = 1, nOutputFile
       tOutput = tFinal/real(nOutputFile)*real(iOutputFile)
       ! In this loop, we first update the RawData and calculate three
       ! Hamiltonian functions and the total control volume; then, we
       ! use subroutine explicit3 to calculate Source_C and update VDF
       do
          if(Time + DtNext >= tOutput) then
             Dt = tOutput - Time
             DoExit = .true.
          else
             DoExit = .false.
             Dt = DtNext
          end if

          if (trim(TypeScatter)=='no' .or.   &
               trim(TypeScatter)=='advect_scatter') then
             ! Use Poisson scheme for advection, and calculate Source_C
             call advance_advect(nP, nQ, nR, DeltaMu,                      &
                  MomentumFace_I, dDeltaSOverBDt_C, DeltaHamiltonian1_N,   &
                  DeltaP3_I, Velocity_I, InvB_C, Hamiltonian2_N,           &
                  DeltaSOverB_C, dLnBDeltaS2Dt_C, bDuDt_C, Hamiltonian3_N, &
                  iStep, Time, Dt, DtNext, Cfl,                            &
                  Volume_G, dVolumeDt_G, VDF_G, Source_C)

             if (DoExit) then
                ! This step is at tOutput
                VDF_G(1:nQ, 1:nP, 1:nR) = VDF_G(1:nQ, 1:nP, 1:nR) + Source_C
             else
                ! This step is not at tOutput
                ! Update VDF_G considering the time-dependent Volume_G

                Source_C = Source_C*Volume_G(1:nQ, 1:nP, 1:nR)
                ! Update DeltaSOverB_C for the calculation of volume
                call update_states(nQ, Time, DeltaSOverB_C = DeltaSOverB_C)
                ! Calculate the total control volume
                do iR = 1, nR
                   Volume_G(1:nQ, 0, iR) = DeltaSOverB_C*   &
                        DeltaMu*DeltaP3_I(iR)
                   do iP = 1, nP+1
                      Volume_G(1:nQ, iP, iR) = Volume_G(1:nQ, 0, iR)
                   end do
                end do

                ! Update VDF_G to CURRENT time: no BCs for (1:nQ, 1:nP, 1;nR)
                VDF_G(1:nQ, 1:nP, 1:nR) = VDF_G(1:nQ, 1:nP, 1:nR) +  &
                     Source_C/Volume_G(1:nQ, 1:nP, 1:nR)
             end if

             ! If there is no scatter, there is only advection
             ! If there is advect_scatter, we will include scattering effect
             if (trim(TypeScatter)=='advect_scatter') then
                call calc_scatter(nQ, nP, nR, Dt, DeltaMu,     &
                     Velocity_I, InvB_C, RawData1_II(:, rho_), &
                     LambdaMuMu_II, VDF_G(1:nQ, 1:nP, 1:nR))
             end if
          elseif (trim(TypeScatter)=='diffuse_scatter') then
             ! There is only scattering by diffusion
             call diffuse_scatter(nP, nQ, nR, Dt, DeltaMu,  &
                  Velocity_I, InvB_C, DeltaSface_I,         &
                  LambdaMuMu_II, VDF_G(1:nQ, 1:nP, 1:nR))
          else
             ! Unrecognized specified scatter type
             call CON_stop('test_multi_poisson: '//   &
                  'Unknown TypeScatter = '//TypeScatter)
          end if

          call set_VDF_Bc      ! Set VDF BCs (mu: reflection; others: outflow)
          iStep = iStep + 1    ! March step
          Time = Time + DtNext ! March time
          write(*,*) iStep, Time
          if(DoExit) EXIT      ! Reach tOutput, exit the current loop

       end do

       ! Set a min value for VDF_G
       VDF_G = max(VDF_G, cMinVal)
       ! Set the name of files, with the largest number being 999
       write(NameFileSuffix, '(I3.3)') iOutputFile
       ! Calculate VDF for output: integrate over the \mu axis
       VDFOutput_II = sum(VDF_G(1:nQ, 1:nP, 1:nR), dim=2)*DeltaMu

       ! Save the files
       call save_plot_file(NameFile='test_multipoisson_'          &
            //NameFileSuffix//'.out',   TypeFileIn = 'ascii',     &
            TimeIn=Time, nStepIn=iStep, NameVarIn='sL logE VDF',  &
            CoordMinIn_D = [HelioDist_I(1) ,    Energy_I(1)],     &
            CoordMaxIn_D = [HelioDist_I(nQ),    Energy_I(nR)],    &
            StringFormatIn = '(3F15.7)',                          &
            Coord1In_I = HelioDist_I,   Coord2In_I = Energy_I,    &
            VarIn_II = log10(VDFOutput_II(1:nQ, 1:nR)))
    end do

  contains
    !==========================================================================
    subroutine read_fieldline
      ! To read the raw data of fieldlines from the MHD simulations
      ! Import data from files
      !------------------------------------------------------------------------
      open(unit=nInFileOld, file='FieldLineOld.in')
      open(unit=nInFileNew, file='FieldLineNew.in')

      ! Get the lines that we need from the OLD data at OLD time
      GapData = 0.0
      do while(nint(GapData) /= iLagrIdSkip)
         read(nInFileOld, *) GapData
      end do
      do iQ = 1, nQ
         read(nInFileOld,*) (RawData1_II(iQ, iVar), iVar = 1, nVar)
         do iStride = 1, nStride-1
            read(nInFileOld,*)
         end do
      end do
      close(nInFileOld)

      ! Get the lines that we need from the NEW data at NEW time
      GapData = 0.0
      do while(nint(GapData) /= iLagrIdSkip)
         read(nInFileNew, *) GapData
      end do
      do iQ = 1, nQ
         read(nInFileNew,*) (RawData2_II(iQ, iVar), iVar = 1, nVar)
         do iStride = 1, nStride-1
            read(nInFileNew,*)
         end do
      end do
      close(nInFileNew)

      ! Normalization process: due to the differences of the units
      RawData1_II(:, x_ :z_ ) = RawData1_II(:, x_ :z_ )*rSun/cLightSpeed
      RawData1_II(:, ux_:uz_) = RawData1_II(:, ux_:uz_)/cLightSpeed
      RawData2_II(:, x_ :z_ ) = RawData2_II(:, x_ :z_ )*rSun/cLightSpeed
      RawData2_II(:, ux_:uz_) = RawData2_II(:, ux_:uz_)/cLightSpeed
    end subroutine read_fieldline
    !==========================================================================
    subroutine init_test_VDF
      ! Initialize the distribution function VDF_G
      !------------------------------------------------------------------------
      VDF_G = cMinVal

      ! Assign the VDF at positive \mu
      do iR = 1, nR
         do iQ = 1, nQ
            if(iQ <= ParticleRange1) then
               ! Assume momentum follow a power distribution: f ~ p^{-5}
               VDF_G(iQ, nP/2:nP, iR) = max(0.1*exp(-5.0/3.0* &
                    (LnP3min + iR*DeltaLnP3)), cMinVal)
            elseif(iQ <= ParticleRange2) then
               ! Use a linear slope for the middle part
               VDF_G(iQ, nP/2:nP, iR) = max(0.1*exp(-5.0/3.0* &
                    (LnP3min + iR*DeltaLnP3))*(ParticleRange2 - real(iQ))/ &
                    (ParticleRange2 - ParticleRange1), cMinVal)
            end if
         end do
      end do

      ! Set VDF Boundary Conditions (reflection for mu; outflow for others)
      call set_VDF_Bc
    end subroutine init_test_VDF
    !==========================================================================
    subroutine set_VDF_Bc
      ! Set Boundary Conditions for VDF

      ! Use reflective boundary condition for \mu axis
      !------------------------------------------------------------------------
      VDF_G(:,    0, :) = VDF_G(:,    1, :)
      VDF_G(:,   -1, :) = VDF_G(:,    2, :)
      VDF_G(:, nP+1, :) = VDF_G(:,   nP, :)
      VDF_G(:, nP+2, :) = VDF_G(:, nP-1, :)

      ! Use outflow boundary condition for s_L axis
      VDF_G(-1,   :, :) = VDF_G(1,    :, :)
      VDF_G(0,    :, :) = VDF_G(1,    :, :)
      VDF_G(nQ+1, :, :) = VDF_G(nQ,   :, :)
      VDF_G(nQ+2, :, :) = VDF_G(nQ,   :, :)

      ! Use outflow boundary condition for Ln(p^3/3) axis
      VDF_G(:, :,   -1) = VDF_G(:, :,    1)
      VDF_G(:, :,    0) = VDF_G(:, :,    1)
      VDF_G(:, :, nR+1) = VDF_G(:, :,   nR)
      VDF_G(:, :, nR+2) = VDF_G(:, :,   nR)
    end subroutine set_VDF_Bc
    !==========================================================================
    subroutine calc_data_states
      ! Calculate data states from the input files
      real :: MagneticOverWave_C(nQ), LarmorRadius13_C(nQ), Lmax_C(nQ)
      real, parameter :: CoeffMu = 6.0/(2.0**(2.0/3.0)*cPi**(5.0/3.0))

      ! Calculate values at OLD time
      ! Calculate midpoints
      !------------------------------------------------------------------------
      MidPoint_ID(1:nQ-1, x_:z_) = (RawData1_II(2:nQ,x_:z_)    &
           + RawData1_II(1:nQ-1,x_:z_))*0.5
      ! Calculate DeltaS
      DeltaS_I(2:nQ-1) = sqrt(sum((MidPoint_ID(2:nQ-1, x_:z_)  &
           - MidPoint_ID(1:nQ-2, x_:z_))**2, dim=2))
      ! Linear interpolate the deltas such that there will be nQ deltas
      DeltaS_I(1 ) = 2*sqrt(sum( &
           (MidPoint_ID(1, x_:z_) - RawData1_II(1, x_:z_))**2))
      DeltaS_I(nQ) = 2*sqrt(sum( &
           (MidPoint_ID(nQ-1, x_:z_) - RawData1_II(nQ, x_:z_))**2))
      ! Calculate 1/2B and \DeltaS/B at grid center
      InvBOld_C = 1.0/(2.0*sqrt(sum(RawData1_II(:, Bx_:Bz_)**2, dim=2)))
      DeltaSOverBOld_C = DeltaS_I*2.0*InvBOld_C
      LnBDeltaS2Old_C = log(1.0/(2.0*InvBOld_C)*DeltaS_I**2)

      ! At OLD time, we calculate LambdaMuMu_II
      ! Calculate DeltaSface_I: Distance of face at OLD time
      DeltaSface_I = sqrt(sum((RawData1_II(2:nQ, x_:z_)        &
           - RawData1_II(1:nQ-1, x_:z_))**2, dim=2))
      ! Calculate quantites related to magnetic field
      MagneticOverWave_C = (1.0/(2.0*InvBOld_C)**2 +           &
           cMu*sum(RawData1_II(:, wave1_:wave2_), dim=2))/     &
           (cMu*sum(RawData1_II(:, wave1_:wave2_), dim=2))
      LarmorRadius13_C = (1.0e9/(cLightSpeed**2/(2.0*InvBOld_C)))**(1.0/3.0)
      ! Mean free path in the unit of Rsun, i.e., Lmax = 0.03*Rsun
      Lmax_C = 0.03*sqrt(sum(RawData1_II(:, x_:z_)**2, dim=2))
      ! Calculate \lambda_mumu and \lambda_xx
      do iR = 1, nR
         LambdaMuMu_II(:, iR) = CoeffMu*MagneticOverWave_C*    &
              LarmorRadius13_C*Lmax_C**(2.0/3.0)*Momentum_I(iR)**(1.0/3.0)
      end do

      ! Calculate values at NEW time
      ! Calculate midpoints
      MidPoint_ID(1:nQ-1, x_:z_) = (RawData2_II(2:nQ,x_:z_)    &
           + RawData2_II(1:nQ-1,x_:z_))*0.5
      ! Calculate DeltaS
      DeltaS_I(2:nQ-1) = sqrt(sum((MidPoint_ID(2:nQ-1, x_:z_)  &
           - MidPoint_ID(1:nQ-2, x_:z_))**2, dim=2))
      ! Linear interpolate the deltas such that there will be nQ deltas
      DeltaS_I(1 ) = 2*sqrt(sum( &
           (MidPoint_ID(1, x_:z_) - RawData2_II(1, x_:z_))**2))
      DeltaS_I(nQ) = 2*sqrt(sum( &
           (MidPoint_ID(nQ-1, x_:z_) - RawData2_II(nQ, x_:z_))**2))
      ! Calculate 1/2B and \DeltaS/B at grid center
      InvBNew_C = 1.0/(2.0*sqrt(sum(RawData2_II(:, Bx_:Bz_)**2, dim=2)))
      DeltaSOverBNew_C = DeltaS_I*2.0*InvBNew_C
      LnBDeltaS2New_C = log(1.0/(2.0*InvBNew_C)*DeltaS_I**2)

    end subroutine calc_data_states
    !==========================================================================
  end subroutine test_multi_poisson
  !============================================================================
  subroutine update_states(nQ, Time, DeltaSOverB_C, &
       dDeltaSOverBDt_C, InvB_C, bDuDt_C, dLnBDeltaS2Dt_C)
    ! Update states according to the current time

    integer, intent(in) :: nQ     ! Number of grid along s_L axis
    real, intent(in)    :: Time   ! Current time
    real, optional, intent(inout) :: DeltaSOverB_C(nQ),   &
         dDeltaSOverBDt_C(nQ), InvB_C(nQ), bDuDt_C(nQ), dLnBdeltaS2Dt_C(nQ)
    ! Magnetic field strength B at cell center
    real :: B_C(nQ, 3)

    ! Calculate values for CURRENT time: here we use linear interpolation
    ! to get the data of each time step from every to consecutive files
    !--------------------------------------------------------------------------
    ! Calculate B at cell center
    if (present(bDuDt_C)) B_C = RawData1_II(:, Bx_:Bz_) + &
         (RawData2_II(:, Bx_:Bz_) - RawData1_II(:, Bx_:Bz_))/tEachFile*Time
    ! Calculate 1/2B
    if (present(InvB_C)) InvB_C = InvBOld_C + &
         (InvBNew_C - InvBOld_C)/tEachFile*Time
    ! Calculate \deltas/B
    if (present(DeltaSOverB_C)) DeltaSOverB_C = DeltaSOverBOld_C + &
         (DeltaSOverBNew_C - DeltaSOverBOld_C)/tEachFile*Time
    ! Calculate \deltas/B time derivative (for next time)
    if (present(dDeltaSOverBDt_C)) dDeltaSOverBDt_C =  &
         (DeltaSOverBNew_C - DeltaSOverBOld_C)/tEachFile
    ! Calculate Dln(B\deltas^2)/Dt
    if (present(dLnBdeltaS2Dt_C)) dLnBdeltaS2Dt_C = &
         (LnBDeltaS2New_C - LnBDeltaS2Old_C)/tEachfile
    ! Calculate b*Du/Dt
    if (present(bDuDt_C)) bDuDt_C = sum(B_C *                         &
         (RawData2_II(:, ux_:uz_) - RawData1_II(:, ux_:uz_)), dim=2)* &
         2.0*InvB_C/tEachFile
  end subroutine update_states
  !============================================================================
  subroutine calc_hamiltonian_1(nP, nQ, nR, DeltaMu,  &
       MomentumFace_I, dDeltaSOverBDt_C, DeltaHamiltonian1_N)
    ! Calculate the 1st Hamiltonian function with time:
    ! p^3/3*\deltas/B, regarding to tau and p^3/3

    integer, intent(in) :: nP, nQ, nR     ! Number of s_L, mu, ln(p^3/3) grids
    real, intent(in)    :: DeltaMu        ! Distance of adjacent two mu's
    real, intent(in)    :: MomentumFace_I(0:nR) ! Momentum at cell face
    real, intent(in)    :: dDeltaSOverBDt_C(nQ) ! D[delta(s_L)/B]/Dt
    real, intent(inout) :: DeltaHamiltonian1_N(0:nQ+1, 0:nP+1, -1:nR+1)
    integer             :: iP, iR         ! Loop variables
    ! Calculate the first Hamiltonian function
    !--------------------------------------------------------------------------

    do iR = 0, nR
       do iP = 0, nP+1
          DeltaHamiltonian1_N(1:nQ, iP, iR) =   &
               -MomentumFace_I(iR)**3.0/3.0*dDeltaSOverBDt_C
       end do
    end do

    ! Calculate the Hamiltonian function used actually：\tilde\deltaH
    DeltaHamiltonian1_N(1:nQ, 0:nP+1, 0:nR) =   &
         DeltaHamiltonian1_N(1:nQ, 0:nP+1, 0:nR)*DeltaMu

    ! Boundary condition of Hamiltonian function
    DeltaHamiltonian1_N(0   , :,    :) = DeltaHamiltonian1_N(1 , :,  :)
    DeltaHamiltonian1_N(nQ+1, :,    :) = DeltaHamiltonian1_N(nQ, :,  :)
    DeltaHamiltonian1_N(:   , :,   -1) = DeltaHamiltonian1_N(: , :,  0)
    DeltaHamiltonian1_N(:   , :, nR+1) = DeltaHamiltonian1_N(: , :, nR)
  end subroutine calc_hamiltonian_1
  !============================================================================
  subroutine calc_hamiltonian_2(nP, nQ, nR, DeltaMu,  &
       DeltaP3_I, Velocity_I, InvB_C, Hamiltonian2_N)
    ! Calculate the 2nd Hamiltonian function at each fixed time:
    ! (mu^2-1)*v/(2B), regarding to s_L and mu

    integer, intent(in) :: nP, nQ, nR     ! Number of s_L, mu, ln(p^3/3) grids
    real, intent(in)    :: DeltaMu        ! Distance of adjacent two mu
    real, intent(in)    :: DeltaP3_I(nR)  ! Distance of adjacent two p^3/3
    real, intent(in)    :: Velocity_I(nR) ! Particle velocity array
    real, intent(in)    :: InvB_C(nQ)     ! 1/(2B) array
    real, intent(inout) :: Hamiltonian2_N(-1:nQ+1, -1:nP+1, 0:nR+1)
    real                :: InvB_F(0:nQ)   ! Intermediate array
    integer             :: iP, iR         ! Loop variables
    !--------------------------------------------------------------------------
    ! Calculate 1/(2B) on the boundary of grid
    InvB_F(1:nQ-1) = (InvB_C(2:nQ) + InvB_C(1:nQ-1))*0.5
    InvB_F(0 )     = InvB_C(1 ) - (InvB_C(2 ) - InvB_C(1   ))*0.5
    InvB_F(nQ)     = InvB_C(nQ) + (InvB_C(nQ) - InvB_C(nQ-1))*0.5

    ! Calculate the second hamiltonian function = (mu^2-1)*v/(2B)
    do iR = 1, nR
       do iP = 0, nP
          Hamiltonian2_N(0:nQ, iP, iR) = ((iP*DeltaMu - 1.0)**2 &
               - 1.0)*Velocity_I(iR)*InvB_F
       end do
    end do

    ! Calculate the Hamiltonian function used actually: \tilde\deltaH
    do iR = 1, nR
       Hamiltonian2_N(0:nQ, 0:nP, iR) = &
            Hamiltonian2_N(0:nQ, 0:nP, iR)*DeltaP3_I(iR)
    end do

    ! Boundary condition of Hamiltonian function
    Hamiltonian2_N(-1  , :   , :   ) = Hamiltonian2_N(0 , :   , : )
    Hamiltonian2_N(nQ+1, :   , :   ) = Hamiltonian2_N(nQ, :   , : )
    Hamiltonian2_N(:   , -1  , :   ) = Hamiltonian2_N(: , 1   , : )
    Hamiltonian2_N(:   , nP+1, :   ) = Hamiltonian2_N(: , nP-1, : )
    Hamiltonian2_N(:   , :   , 0   ) = Hamiltonian2_N(: , :   , 1 )
    Hamiltonian2_N(:   , :   , nR+1) = Hamiltonian2_N(: , :   , nR)
  end subroutine calc_hamiltonian_2
  !============================================================================
  subroutine calc_hamiltonian_3(nP, nQ, nR, DeltaMu, MomentumFace_I,  &
       DeltaSOverB_C, dLnBDeltaS2Dt_C, bDuDt_C, Hamiltonian3_N)
    ! Calculate the 3rd Hamiltonian function at each fixed time:
    ! (1-mu^2)/2 * [ mu*(p^3/3)*(3\vec{b}\vec{b}:\nabla\vec{u} - div\vec{u})
    ! + ProtonMass*p^2*(\vec{b}*d\vec{u}/dt) ], regarding to p^3/3 and mu
    ! Here we list the variables in the analytical function:
    ! (3\vec{b}\vec{b}:\nabla\vec{u} - div\vec{u}) = d(ln(B*ds^2))/dt
    ! ProtonMass = cRmeProtonGeV, in the unit of GeV/c^2
    ! \vec{b}*d\vec{u}/dt = bDuDt_C

    integer, intent(in) :: nP, nQ, nR     ! Number of s_L, mu, ln(p^3/3) grids
    real, intent(in)    :: DeltaMu        ! Distance of adjacent two mu
    real, intent(in)    :: MomentumFace_I(0:nR) ! Momentum at cell face
    real, intent(in)    :: DeltaSOverB_C(nQ), dLnBDeltaS2Dt_C(nQ), bDuDt_C(nQ)
    real, intent(inout) :: Hamiltonian3_N(0:nQ+1, -1:nP+1, -1:nR+1)
    integer             :: iP, iQ, iR     ! Loop variables
    ! Calculate the third hamiltonian function
    !--------------------------------------------------------------------------
    do iR = 0, nR
       do iP = 0, nP
          Hamiltonian3_N(1:nQ, iP, iR) = 0.5*(1.0 - (iP*DeltaMu-1.0)**2)*   &
               ((iP*DeltaMu-1.0)*MomentumFace_I(iR)**3.0/3.0*               &
               dLnBDeltaS2Dt_C + cRmeProtonGeV*MomentumFace_I(iR)**2*bDuDt_C)
       end do
    end do

    ! Calculate the Hamiltonian function used actually：\tilde\deltaH
    do iQ = 1, nQ
       Hamiltonian3_N(iQ, 0:nP, 0:nR) = &
            Hamiltonian3_N(iQ,0:nP,0:nR)*DeltaSOverB_C(iQ)
    end do

    ! Boundary condition of Hamiltonian function
    Hamiltonian3_N(0   , :   , :   ) = Hamiltonian3_N(1 , :   , : )
    Hamiltonian3_N(nQ+1, :   , :   ) = Hamiltonian3_N(nQ, :   , : )
    Hamiltonian3_N(:   , -1  , :   ) = Hamiltonian3_N(: , 1   , : )
    Hamiltonian3_N(:   , nP+1, :   ) = Hamiltonian3_N(: , nP-1, : )
    Hamiltonian3_N(:   , :   , -1  ) = Hamiltonian3_N(: , :   , 0 )
    Hamiltonian3_N(:   , :   , nR+1) = Hamiltonian3_N(: , :   , nR)
  end subroutine calc_hamiltonian_3
  !============================================================================
  subroutine advance_advect(nP, nQ, nR, DeltaMu,               &
       MomentumFace_I, dDeltaSOverBDt_C, DeltaHamiltonian1_N,  &
       DeltaP3_I, Velocity_I, InvB_C, Hamiltonian2_N,          &
       DeltaSOverB_C, dLnBDeltaS2Dt_C, bDuDt_C, Hamiltonian3_N,&
       iStep, Time, Dt, DtNext, CflIn,                         &
       Volume_G, dVolumeDt_G, VDF_G, Source_C)
    ! Advect transport equation by multi_Poisson scheme

    integer, intent(in) :: nP, nQ, nR     ! Number of s_L, mu, ln(p^3/3) grids
    real, intent(in)    :: DeltaMu        ! Distance of adjacent two mu's
    real, intent(in)    :: MomentumFace_I(0:nR) ! Momentum at cell face
    real, intent(inout) :: dDeltaSOverBDt_C(nQ) ! D[delta(s_L)/B]/Dt
    real, intent(inout) :: DeltaHamiltonian1_N(0:nQ+1, 0:nP+1, -1:nR+1)
    real, intent(in)    :: DeltaP3_I(nR)  ! Distance of adjacent two p^3/3
    real, intent(in)    :: Velocity_I(nR) ! Particle velocity array
    real, intent(inout) :: InvB_C(nQ)     ! 1/(2B) array
    real, intent(inout) :: Hamiltonian2_N(-1:nQ+1, -1:nP+1, 0:nR+1)
    real, intent(inout) :: DeltaSOverB_C(nQ), dLnBDeltaS2Dt_C(nQ), bDuDt_C(nQ)
    real, intent(inout) :: Hamiltonian3_N(0:nQ+1, -1:nP+1, -1:nR+1)
    integer, intent(in) :: iStep          ! Current step index
    real, intent(in)    :: Time           ! Current physical time
    real, intent(inout) :: Dt             ! DtIn, can be reduced
    real, intent(inout) :: DtNext         ! DtOut, i.e., DtNext
    real, intent(in)    :: CflIn          ! Input CFL number
    real, intent(inout) :: Volume_G(0:nQ+1, 0:nP+1, 0:nR+1)    ! Control Volume
    real, intent(in)    :: dVolumeDt_G(0:nQ+1, 0:nP+1, 0:nR+1) ! dVolume/dt
    real, intent(inout) :: VDF_G(-1:nQ+2, -1:nP+2, -1:nR+2)    ! VDF
    real, intent(inout) :: Source_C(nQ, nP, nR) ! Source at each time step
    integer             :: iP, iR         ! Loop variables

    ! Advect by multiple Poisson bracket considering the pitch angle
    !--------------------------------------------------------------------------
    call update_states(nQ, Time,                &
         DeltaSOverB_C = DeltaSOverB_C,         &
         dDeltaSOverBDt_C = dDeltaSOverBDt_C,   &
         InvB_C = InvB_C, bDuDt_C = bDuDt_C,    &
         dLnBDeltaS2Dt_C = dLnBDeltaS2Dt_C)     ! Update states first
    Hamiltonian2_N = 0.0
    Hamiltonian3_N = 0.0

    ! Calculate 1st Hamiltonian function used in the time-dependent
    ! poisson bracket: {f_jk, (p^3/3)*(DeltaS/B)}_{tau, p^3/3}
    call calc_hamiltonian_1(nP, nQ, nR, DeltaMu,   &
         MomentumFace_I, dDeltaSOverBDt_C, DeltaHamiltonian1_N)
    ! Calculate 2nd Hamiltonian function used in the time-dependent
    ! poisson bracket: {f_jk, (mu^2-1)*v/(2B)}_{s_L, mu}
    call calc_hamiltonian_2(nP, nQ, nR, DeltaMu,   &
         DeltaP3_I, Velocity_I, InvB_C, Hamiltonian2_N)
    ! Calculate 3rd Hamiltonian function used in the time-dependent
    ! poisson bracket: {f_jk, (1-mu^2)/2 * [ mu*(p^3/3)*
    ! d(ln(B*ds^2))/dt + ProtonMass*p^2*bDuDt_C ]}_{p^3/3, mu}
    call calc_hamiltonian_3(nP, nQ, nR, DeltaMu, MomentumFace_I,  &
         DeltaSOverB_C, dLnBDeltaS2Dt_C, bDuDt_C, Hamiltonian3_N)

    ! Calculate the total control volume
    do iR = 1, nR
       Volume_G(1:nQ, 0, iR) = DeltaSOverB_C*DeltaMu*DeltaP3_I(iR)
       do iP = 1, nP+1
          Volume_G(1:nQ, iP, iR) = Volume_G(1:nQ, 0, iR)
       end do
    end do

    ! Boundary conditions for total control volume
    Volume_G(0,    :,    :) = Volume_G(1 , :,  :)
    Volume_G(nQ+1, :,    :) = Volume_G(nQ, :,  :)
    Volume_G(:,    :,    0) = Volume_G(: , :,  1)
    Volume_G(:,    :, nR+1) = Volume_G(: , :, nR)

    ! First trial: Need DtNext as DtOut
    if (iStep == 0) then
       call explicit(nQ, nP, nR, VDF_G, Volume_G, Source_C, &
            Hamiltonian12_N = Hamiltonian2_N,               &
            Hamiltonian23_N = -Hamiltonian3_N,              &
            dHamiltonian03_FZ = DeltaHamiltonian1_N,        &
            dVolumeDt_G = dVolumeDt_G,                      &
            DtOut = DtNext, CFLIn=CflIn)
       ! Subsequent trials: with DtIn=Dt and DtOut=DtNext
    else
       call explicit(nQ, nP, nR, VDF_G, Volume_G, Source_C, &
            Hamiltonian12_N = Hamiltonian2_N,               &
            Hamiltonian23_N = -Hamiltonian3_N,              &
            dHamiltonian03_FZ = DeltaHamiltonian1_N,        &
            dVolumeDt_G = dVolumeDt_G,                      &
            DtIn = Dt, DtOut = DtNext, CFLIn=CflIn)
    end if

  end subroutine advance_advect
  !============================================================================
  subroutine calc_scatter(nQ, nP, nR, Dt, DeltaMu,    &
       Velocity_I, InvB_C, n_I, LambdaMuMu_II, VDF_G)
    ! Calculate scatter: \deltaf/\deltat = (Dmumu*f_mu)_mu

    use ModDiffusion, ONLY: tridiag
    integer, intent(in) :: nQ, nP, nR     ! Number of s_L, mu, ln(p^3/3) grids
    real, intent(in) :: Dt                ! Time step
    real, intent(in) :: DeltaMu           ! Distance of adjacent two mu
    real, intent(in) :: Velocity_I(nR)    ! Particle velocity in units of c
    real, intent(in) :: InvB_C(nQ)        ! 1/(2B) along s_L axis
    real, intent(in) :: n_I(nQ)           ! Density in the unit of amu/m^3
    real, intent(in) :: LambdaMuMu_II(nQ, nR)   ! lambda_\mu\mu from data
    real, intent(inout) :: VDF_G(nQ, nP, nR)    ! VDF at s_L, mu, ln(p^3/3)

    ! LOCAL VARs
    ! Dmumu for each fixed iQ and iR. Dmumu is the coefficient
    ! of diffusion along the pitch angle, \mu
    real :: DMuMu_I(0:nP)
    ! Factorize DMuMu, each factor being only a function of s_L, \mu, ln(p^3/3)
    real :: FactorMu_F(0:nP)
    ! Lower, main, upper diagonal, output values for the scatter calculation
    real :: L_I(nP), M_I(nP), U_I(nP), W_I(nP)
    ! Physical VARs
    real :: LowerLimMu           ! Lower limit of Factor_mu
    real :: DtOverDMu2           ! (\Delta t) / (\Delta \mu^2)
    real :: Mu_I(0:nP)           ! Face centered value of the pitch angle
    real :: AlfvenSpeed_I(nQ)    ! Alfven wave speed
    integer :: iQ, iP, iR        ! Loop integers
    integer :: switch1, Pnumber  ! Control parameter
    !--------------------------------------------------------------------------

    DtOverDMu2 = Dt/DeltaMu**2   ! (\Delta t) / (\Delta \mu^2)

    ! Calculate factorized diffusion coefficient
    LowerLimMu = (1.0 - DeltaMu**2)*abs(DeltaMu)**(2.0/3.0)
    do iP = 0, nP
       Mu_I(iP) = -1.0 + real(iP)*DeltaMu
       FactorMu_F(iP) = (1.0 - Mu_I(iP)**2)*abs(Mu_I(iP))**(2.0/3.0)
    end do

    ! Get the Alfven wave speed for each s_L
    AlfvenSpeed_I = 0.5/(InvB_C*cLightSpeed*sqrt(n_I*cAtomicMass*cMu))
    ! Calculate the effect of scatter along \mu axis for each fixed iR and iQ
    do iR = 1, nR
       do iQ = 1, nQ
          DMuMu_I = 0.0
          switch1 = 1

          ! For each pitch angle, we will set DMuMu and solve VDF for scatter
          do iP = 0, nP
             ! Control whether we will floor the value of D_\mu\mu or not
             if(Velocity_I(iR)*abs(Mu_I(iP)) >= 10.0*AlfvenSpeed_I(iQ)) then
                DMuMu_I(iP) = Velocity_I(iR)/LambdaMuMu_II(iQ, iR)* &
                     FactorMu_F(iP)*DtOverDMu2
             else
                if(switch1==1) Pnumber = iP
                switch1 = 0
                ! Set the lower limit of D_\mu\mu when |\mu| is close to zero
                DMuMu_I(iP) = Velocity_I(iR)/LambdaMuMu_II(iQ,iR)* &
                     max(FactorMu_F(Pnumber), LowerLimMu)*DtOverDMu2
             end if
          end do

          ! Set up coefficients for solving the linear equation set
          L_I = -DMuMu_I(0:nP-1)
          U_I = -DMuMu_I(1:nP  )
          M_I = 1.0 - L_I - U_I
          W_I = VDF_G(iQ, 1:nP, iR)

          ! For each pitch angle, we finally solve VDF for scattering effect
          call tridiag(nP, L_I, M_I, U_I, W_I, VDF_G(iQ, 1:nP, iR))
       end do
    end do
  end subroutine calc_scatter
  !============================================================================
  subroutine diffuse_scatter(nP, nQ, nR, Dt, DeltaMu, &
       Velocity_I, InvB_C, DeltaSface_I, LambdaMuMu_II, VDF_G)
    ! Calculate scatter_diffusion: \deltaf/\deltat = (Dmumu*f_mu)_mu

    use ModDiffusion, ONLY: advance_diffusion1
    integer, intent(in) :: nQ, nP, nR     ! Number of s_L, mu, ln(p^3/3) grids
    real, intent(in) :: Dt                ! Time step
    real, intent(in) :: DeltaMu           ! Distance of adjacent two mu
    real, intent(in) :: Velocity_I(nR)    ! Particle velocity in units of c
    real, intent(in) :: InvB_C(nQ)        ! 1/(2B) along s_L axis
    real, intent(in) :: DeltaSface_I(nQ-1)! Distance between s_L centers
    real, intent(in) :: LambdaMuMu_II(nQ, nR)   ! lambda_\mu\mu from data
    real, intent(inout) :: VDF_G(nQ, nP, nR)    ! VDF at s_L, mu, ln(p^3/3)

    ! LOCAL VARs

    real :: DInner_I(nQ), DOuter_I(nQ)    ! Inner/Outer diffusion coffecients
    real :: Dist_I(nQ)                    ! Distance between s_L centers
    real :: dxx(nQ)                       ! Dxx diffusion coffecient
    integer :: iQ, iP, iR                 ! Loop integers

    !--------------------------------------------------------------------------
    Dist_I = 1.0                          ! Initial Dist_I
    Dist_I(1: nQ-1) = DeltaSface_I        ! Assign DeltaSface_I(1:nQ-1)

    ! Calculate the effect of scatter along s_L axis for each fixed iR and iP
    do iR = 1, nR
       do iP = 1, nP
          ! Set inner and outer coefficients for diffusion equation
          DInner_I = 2.0/3.0*Velocity_I(iR)*CoeffMuToxx* &
               LambdaMuMu_II(:, iR)*InvB_C
          DOuter_I = 0.50/InvB_C          ! There should be no Dt here!!!
          ! Finally we diffuse the VDF for scattering effect
          call advance_diffusion1(Dt, 0, nQ, 0, Dist_I,  &
               VDF_G(1:nQ, iP, iR), DOuter_I, DInner_I)
       end do
    end do
    dxx = DInner_I*0.5/InvB_C
  end subroutine diffuse_scatter
  !============================================================================
end module ModTestMultiPoisson
!==============================================================================
module ModHillVortex

  use ModPoissonBracket, ONLY: explicit
  use ModUtilities, ONLY: CON_stop
  use ModPlotFile, ONLY: save_plot_file
  use ModConst

  implicit none

  PRIVATE ! Except
  public :: test_hill_vortex

contains
  !============================================================================
  subroutine test_hill_vortex

    ! Grid in RSph-Theta variables
    ! Theta goes from 0 to 360 because it is mirrored over the pole for plots
    integer, parameter:: nR = 54, nTheta = 2*nR
    real,    parameter:: Dr = 4.50/nR, Cfl = 0.9

    ! Streamlines:
    integer, parameter:: r_ = 2, Ur_ = 2, z_ = 1, Uz_  = 1

    ! Position and width of the smooth bump or sharp ellipsoid
    real, parameter:: xCenter = -3.0, yCenter = 0.0, WidthX = 3.0, WidthY = 4.0
    logical:: IsSmooth = .true.

    real:: x, y, z, r, Uz, Ur, Vel_VC(Uz_:Ur_,-500:500,-500:500), RSph

    ! Loop variables
    integer:: i, j, iR, iTheta, iStep, iPlot

    ! Mesh size in Theta
    real, parameter   :: dTheta = cTwoPi/nTheta
    real :: VDF_G(-1:nR+2,-1:nTheta), Volume_G(0:nR+1, 0:nTheta/2+1)
    real :: R_I(-1:nR+1), Theta_I(-1:nTheta+1), CosTheta_I(-1:nTheta/2+1)
    real :: Hamiltonian_N(-1:nR+1, -1:nTheta/2+1)
    real :: Time, Dt, Source_C(nR,nTheta/2), tFinal
    real :: PlotVar_VC(5,nR,nTheta), DfDt_C(nR,nTheta), DfDt2_C(nR,nTheta)
    real :: Error, Error2, URSph, UTheta
    logical :: DoExit
    !--------------------------------------------------------------------------
    do j = -500, 500
       r = 0.01*j
       do i = -500, 500
          z = 0.01*i
          RSph = sqrt(z**2 + r**2)
          if(RSph < 1.0)then
             ! Inside the vortex
             Uz =  1.50*(-1.0 + r**2 + RSph**2)
             Ur = -1.50*r*z
          else
             Uz = 1.0 - 1.0/RSph**3 + 1.50*r**2/RSph**5
             Ur = -1.50*r*z/RSph**5
          end if
          Vel_VC(Uz_:Ur_, i, j) = [Uz, Ur]
       end do
    end do
    call save_plot_file( &
         NameFile = 'streamlines.out',                &
         StringHeaderIn='Hill vortex velocity field', &
         NameVarIn='z r Uz Ur'  ,                     &
         CoordMinIn_D=[-5.0, -5.0],                   &
         CoordMaxIn_D=[ 5.0,  5.0 ],                  &
         StringFormatIn = '(4F16.4)',                 &
         VarIn_VII = Vel_VC)

    ! Initialize coords
    ! These are cell face coordinates:
    ! r_I(iR) is the radius at iR+1/2, Theta_I(iTheta) is at iTheta+1/2.
    do iR = -1, nR + 1
       R_I(iR) = iR*Dr + 0.5
    end do
    do iTheta = -1, nTheta + 1
       Theta_I(iTheta) = iTheta*dTheta
    end do
    do iTheta = -1, nTheta/2+1
       CosTheta_I(iTheta) = cos(Theta_I(iTheta))
    end do
    ! Overwrite ghost cell values
    CosTheta_I(nTheta/2+1) = CosTheta_I(nTheta/2-1)
    CosTheta_I(-1) = CosTheta_I(1)

    ! Calculate Hamiltonian
    do iTheta = -1, nTheta/2+1
       ! Calculate theta-dependent factor
       Hamiltonian_N(:, iTheta) = cTwoPi*(1.0 - CosTheta_I(iTheta)**2)
       do iR = -1, nR/9 ! Unit radius
          Hamiltonian_N(iR,iTheta) = Hamiltonian_N(iR,iTheta)*0.750*&
               R_I(iR)**2*(R_I(iR)**2 - 1.0)
       end do
       do iR = nR/9+1, nR+1
          Hamiltonian_N(iR,iTheta) = Hamiltonian_N(iR,iTheta)*0.50*&
               R_I(iR)**2*(1.0 - 1.0/R_I(iR)**3)
       end do
    end do
    ! write(*,*)'Maxval(Hamiltonian)=',maxval(Hamiltonian_N)

    ! Calculate volume
    ! angle-dependent factor
    do iTheta = 1, nTheta/2
       Volume_G(0:nR+1,iTheta) = (CosTheta_I(iTheta-1) - CosTheta_I(iTheta))*&
            (R_I(0:nR+1)**3 - R_I(-1:nR)**3)*cTwoPi/3
    end do
    Volume_G(0:nR+1,0) = (CosTheta_I(0) - CosTheta_I(-1)) &
         *(R_I(0:nR+1)**3 - R_I(-1:nR)**3)*cTwoPi/3
    Volume_G(0:nR+1,nTheta/2+1) = &
         (CosTheta_I(nTheta/2+1) - CosTheta_I(nTheta/2)) &
         *(R_I(0:nR+1)**3 - R_I(-1:nR)**3)*cTwoPi/3
    ! write(*,*)'Minval(Volume)=',minval(Volume_G)
    ! Periodic boundary conditions
    Source_C = 0.0

    ! Initial conditions:
    VDF_G = 0.0; DfDt_C = 0.0; DfDt2_C = 0.0
    do iTheta = 1, nTheta; do iR = 1, nR
       x = cos(Theta_I(iTheta) - dTheta/2)*(R_I(iR) - dR/2) - xCenter
       y = sin(Theta_I(iTheta) - dTheta/2)*(R_I(iR) - dR/2) - yCenter
       if(IsSmooth)then
          ! inside a WidthX*WidthY rectangle centered on xCenter,yCenter
          ! cos^2(kx*(x-xCenter))*cos^2(ky*(y-yCenter))
          if(abs(x) < WidthX/2 .and. abs(y) < WidthY/2) then
             VDF_G(iR,iTheta) = cos(cPi*x/WidthX)**4 * cos(cPi*y/WidthY)**4
             z = x + xCenter; r = y + yCenter
             RSph = sqrt(z**2 + r**2)
             Uz = 1.0 - 1.0/RSph**3 + 1.50*r**2/RSph**5
             Ur = -1.50*r*z/RSph**5
             ! Anaylitical (df/dt=-u.grad f), to compare with numeric source
             DfDt_C(iR,iTheta) =cos(cPi*x/WidthX)**3 *cos(cPi*y/WidthY)**3&
                  *(Uz*cPi/WidthX*sin(cPi*x/WidthX)*cos(cPi*y/WidthY) + &
                  Ur*cPi/WidthY*cos(cPi*x/WidthX)*sin(cPi*y/WidthY))*4.0
          end if
       else
          if( (x/WidthX)**2 + (y/WidthY)**2 < 0.25) &
               VDF_G(iR,iTheta) = 1.0
       end if
    end do; end do
    if(IsSmooth)then
       do iTheta = 1, nTheta; do iR = 1, nR
          x = cos(Theta_I(iTheta) - dTheta/2)*(R_I(iR) - dR/2) - xCenter
          y = sin(Theta_I(iTheta) - dTheta/2)*(R_I(iR) - dR/2) - yCenter
          if(abs(x) < WidthX/2 .and. abs(y) < WidthY/2) then
             VDF_G(iR,iTheta) = cos(cPi*x/WidthX)**4 * cos(cPi*y/WidthY)**4
             z = x + xCenter; r = y + yCenter
             RSph = sqrt(z**2 + r**2)
             Uz = 1.0 - 1.0/RSph**3 + 1.50*r**2/RSph**5
             Ur = -1.50*r*z/RSph**5
             URsph = (Uz*z + Ur*r)/RSph
             UTheta = (-Uz*r + Ur*z)/RSph
             ! By differentiating  df/dt=-u.grad f, the second order derivative
             ! may be evaluated: d2f/dt2=-u.grad(df/dt), so that analytical
             ! df/dt may be differentiated numerically:
             DfDt2_C(iR,iTheta) =-0.50*(&
                  URsph*(DfDt_C(iR+1,iTheta) - DfDt_C(iR-1,iTheta))/Dr + &
                  UTheta*(DfDt_C(iR,iTheta+1) - DfDt_C(iR,iTheta-1))/&
                  (RSph*dTheta) )
          end if
       end do; end do
    end if

    Time = 0.0; iStep = 0

    ! Save initial conditions
    call save_plot_file('hill.outs', 'rewind', 'real8', &
         'Hill vortex', iStep, Time, &
         NameVarIn='r Theta Rho'  , &
         CoordMinIn_D=[0.5 + 0.5*Dr, 180.0/nTheta],&
         CoordMaxIn_D=[5.0 - 0.5*Dr, 360.0 - 180.0/nTheta],&
         VarIn_II = VDF_G(1:nR,1:nTheta) )

    call save_plot_file('hillcut.outs', 'rewind', 'real8', &
         'Hill vortex', iStep, Time, &
         NameVarIn='r Theta Rho'  , &
         CoordMinIn_D=[0.5 + 0.5*Dr, 180 + 180.0/nTheta],&
         CoordMaxIn_D=[5.0 - 0.5*Dr, 360 - 180.0/nTheta],&
         VarIn_II = VDF_G(1:nR,nTheta/2+1:nTheta) )
    if(IsSmooth)then
       ! Calculate source for small time step, to compare with
       ! analytical DfDt
       call explicit(nR, nTheta/2, VDF_G(-1:nR+2,-1:nTheta/2+2), Volume_G,&
            Source_C, Hamiltonian_N,   &
            CFLIn=1.0e-6, DtOut = Dt)
       PlotVar_VC(1,:,:) = VDF_G(1:nR,1:nTheta)
       PlotVar_VC(2,:,:) = DfDt_C(1:nR,1:nTheta)
       PlotVar_VC(3,:,1:nTheta/2) = Source_C/Dt
       Error = sum(abs(PlotVar_VC(2,:,1:nTheta/2) - &
            PlotVar_VC(3,:,1:nTheta/2))*&
            Volume_G(1:nR,1:nTheta/2))/(4/3.0*cPi*5.0**3)
       call explicit(nR, nTheta/2, VDF_G(-1:nR+2,-1:nTheta/2+2), Volume_G,&
            Source_C, Hamiltonian_N,   &
            CFLIn=CFL, DtOut = Dt)
       PlotVar_VC(4,:,:) = 0.50*Dt*DfDt2_C(1:nR,1:nTheta)
       PlotVar_VC(5,:,1:nTheta/2) = Source_C/Dt
       Error2 = sum(abs(&
            PlotVar_VC(2,:,1:nTheta/2) + PlotVar_VC(4,:,1:nTheta/2) - &
            PlotVar_VC(5,:,1:nTheta/2))*&
            Volume_G(1:nR,1:nTheta/2))/(4/3.0*cPi*5.0**3)
       do iTheta = nTheta/2 + 1, nTheta
          ! Symmetric prolongation for visualization
          PlotVar_VC(3,:, iTheta) = PlotVar_VC(3,:, 1+nTheta-iTheta)
          PlotVar_VC(5,:, iTheta) = PlotVar_VC(5,:, 1+nTheta-iTheta)
       end do
       call save_plot_file('InitialSource.out', 'rewind','real8', &
            'Hill vortex', iStep, Time, &
            NameVarIn=&
            'z r Rho DfDt Source DfDt2 Source2 Error Error2', &
            CoordMinIn_D=[0.5 + 0.5*Dr, 180.0/nTheta],&
            CoordMaxIn_D=[5.0 - 0.5*Dr, 360.0 - 180.0/nTheta],&
            ParamIn_I=[Error, Error2],&
            VarIn_VII = PlotVar_VC )
       write(*,*)'Error=', Error, ' Error2=', Error2
       ! stop
    end if
    ! Computation
    do iPlot = 0, 99
       tFinal = (iPlot + 1)*0.1
       DoExit = .false.
       PLOT:do
          call explicit(nR, nTheta/2, VDF_G(-1:nR+2,-1:nTheta/2+2), Volume_G,&
               Source_C, Hamiltonian_N,   &
               CFLIn=Cfl, DtOut = Dt)
          iStep = iStep +1
          if(Time + Dt >= tFinal)then
             Dt = tFinal - Time
             call explicit(nR, nTheta/2, VDF_G(-1:nR+2,-1:nTheta/2+2), &
                  Volume_G, Source_C, Hamiltonian_N, DtIn = Dt)
             VDF_G(1:nR, 1:nTheta/2) = VDF_G(1:nR, 1:nTheta/2) + Source_C
             Time = tFinal
             DoExit = .true.
          else
             Time = Time + Dt
             VDF_G(1:nR, 1:nTheta/2) = VDF_G(1:nR, 1:nTheta/2) + Source_C
          end if
          do iTheta = nTheta/2 + 1, nTheta
             ! Symmetric prolongation
             VDF_G(1:nR, iTheta) = VDF_G(1:nR, 1+nTheta-iTheta)
          end do
          ! Periodic prolongation
          VDF_G(1:nR,-1:0) = VDF_G(1:nR, nTheta-1:nTheta)
          write(*,*)'Time=',Time
          if(DoExit) EXIT PLOT
       end do PLOT
       call save_plot_file('hill.outs', 'append', 'real8', &
            'Hill vortex', iStep, tFinal, &
            NameVarIn='r Theta Rho'  , &
            CoordMinIn_D=[0.50 + 0.50*Dr, 180.0/nTheta],&
            CoordMaxIn_D=[5.0 - 0.50*Dr, 360.0 - 180.0/nTheta],&
            VarIn_II = VDF_G(1:nR,1:nTheta) )
       call save_plot_file('hillcut.outs', 'append', 'real8', &
            'Hill vortex', iStep, tFinal, &
            NameVarIn='r Theta Rho'  , &
            CoordMinIn_D=[0.5 + 0.5*Dr, 180 + 180.0/nTheta],&
            CoordMaxIn_D=[5.0 - 0.5*Dr, 360 - 180.0/nTheta],&
            VarIn_II = VDF_G(1:nR,nTheta/2+1:nTheta) )
    end do
  end subroutine test_hill_vortex
  !============================================================================
end module ModHillVortex
!==============================================================================
module ModStochastic
  use ModPoissonBracket, ONLY: explicit
  use ModUtilities, ONLY: CON_stop
  use ModPlotFile, ONLY: save_plot_file
  use ModConst
  implicit none
  PRIVATE ! Except
  public :: test_stochastic
  ! Streamlines:
  integer, parameter :: nPointPer2Pi = 180
  integer, parameter :: nJ = (4*nPointPer2Pi)/2, nTheta=nPointPer2Pi
  real, parameter    :: Delta = cTwoPi/nPointPer2Pi
  ! Loop variables
  integer           ::  iJ, iTheta, iStep, iPlot
  ! Mesh size in \Theta
  real :: VDF_G(  -1:nTheta+2, -1-nJ:nJ+2)
  real :: Volume_G(0:nTheta+1,   -nJ:nJ+1)
  real :: Theta_I(-1:nTheta+1), Action_I(-1-nJ:nJ+1), Action
  real :: HamiltonianFree_N(-1:nTheta+1, -1-nJ:nJ+1)
  real :: HamiltonianPush_N(-1:nTheta+1, -1-nJ:nJ+1)
  real :: Time, Dt, Source_C(nTheta,1-nJ:nJ), tFinal
  logical :: DoAgain
  character(LEN=11)::NameFile
  !----------------------------------------------------------------------------
contains
  !============================================================================
  subroutine test_stochastic(KChirikov)
    real, intent(in) :: KChirikov ! Left hand side in the Chhirikov criterion
    !--------------------------------------------------------------------------
    ! Uniform Delta*Delta grid:
    Volume_G = Delta**2
    do iJ = -1-nJ, nJ+1
       Action_I(iJ)  = iJ*Delta
    end do
    do iTheta = -1, nTheta+1
       Theta_I(iTheta) = iTheta*Delta
    end do
    ! Free rotator:
    do iJ = -1-nJ, nJ+1
       HamiltonianFree_N(-1:nTheta+1, iJ) = 0.50*Action_I(iJ)**2
    end do
    ! Sinusoidal perturbation
    do iTheta = -1, nTheta+1
       HamiltonianPush_N(iTheta, -1-nJ:nJ+1) = cos(Theta_I(iTheta))
    end do
    ! Account for:
    ! 1. Amplitude factor, KChirikov
    ! 2. Dirac delta-function, by dividing the amplitude by the short time
    !    interval and apply the perturbation only within this interval
    ! 3. Sum up both Hamiltonians:
    HamiltonianPush_N = (KChirikov/0.10)*HamiltonianPush_N + HamiltonianFree_N
    ! Initial distribution function: zero at high energies and in ghost cells:
    VDF_G(-1:nTheta+2,-1-nJ:nJ+2) =  exp(-12.50)
    do iJ  = 1, nJ
       Action = 0.50*(Action_I(iJ-1) + Action_I(iJ))
       if(Action < 1.0)then
          VDF_G(:,iJ) = exp(-12.50*Action**2)
       else
          EXIT
       end if
    end do
    do iJ  = 0, 1-nJ,-1
       Action = 0.50*(Action_I(iJ-1) + Action_I(iJ))
       if(abs(Action) < 1.0)then
          VDF_G(:,iJ) = exp(-12.50*Action**2)
       else
          EXIT
       end if
    end do
    Source_C = 0.0
    tFinal  = 0.0
    ! Compiutation
    Time = 0.0; iStep = 0
    call save_plot_file(NameFile='initial.out', &
         TypeFileIn='real8', TimeIn=tFinal, nStepIn = iStep, &
         NameVarIn='theta JT/2pi Log10VDF'  , &
         CoordMinIn_D=[         0.50*Delta, -2.0 + 0.50/nPointPer2Pi],&
         CoordMaxIn_D=[cTwoPi - 0.50*Delta,  2.0 - 0.50/nPointPer2Pi],&
         VarIn_II = log10(VDF_G(1:nTheta,1-nJ:nJ)) )
    do iPlot = 0, 999
       tFinal = tFinal + 0.1
       DoAgain = .true.
       do while(DoAgain)
          call explicit(nTheta,2*nJ, VDF_G, Volume_G,&
               Source_C, HamiltonianPush_N,   &
               CFLIn=0.5, DtOut = Dt,IsPeriodicIn_D=[.true.,.false.])
          iStep = iStep +1
          if(Time + Dt >= tFinal)then
             Dt = tFinal - Time
             call explicit(nTheta,2*nJ,VDF_G,Volume_G,&
                  Source_C, HamiltonianPush_N,   &
                  DtIn = Dt,IsPeriodicIn_D=[.true.,.false.])
             VDF_G(1:nTheta,1-nJ:nJ) = VDF_G(1:nTheta,1-nJ:nJ) + Source_C
             Time = tFinal
             DoAgain = .false.
          else
             Time = Time + Dt
             VDF_G(1:nTheta,1-nJ:nJ) = VDF_G(1:nTheta,1-nJ:nJ) + Source_C
          end if
          ! Periodic prolongation
          VDF_G(-1:0,:) = VDF_G(nTheta-1:nTheta,:)
          VDF_G(nTheta+1:nTheta+2,:) = VDF_G(1:2,:)
          write(*,*)'Time=',Time
       end do
       tFinal = tFinal + 0.9
       DoAgain = .true.
       do while(DoAgain)
          call explicit(nTheta,2*nJ, VDF_G, Volume_G,&
               Source_C, HamiltonianFree_N,   &
               CFLIn=0.99, DtOut = Dt,IsPeriodicIn_D=[.true.,.false.])
          iStep = iStep +1
          if(Time + Dt >= tFinal)then
             Dt = tFinal - Time
             call explicit(nTheta,2*nJ,VDF_G,Volume_G,&
                  Source_C, HamiltonianFree_N,   &
                  DtIn = Dt, IsPeriodicIn_D=[.true.,.false.])
             VDF_G(1:nTheta,1-nJ:nJ) = VDF_G(1:nTheta,1-nJ:nJ) + Source_C
             Time = tFinal
             DoAgain = .false.
          else
             Time = Time + Dt
             VDF_G(1:nTheta,1-nJ:nJ) = VDF_G(1:nTheta,1-nJ:nJ) + Source_C
          end if
          ! Periodic prolongation
          VDF_G(-1:0,:) = VDF_G(nTheta-1:nTheta,:)
          VDF_G(nTheta+1:nTheta+2,:) = VDF_G(1:2,:)
          write(*,*)'Time=',Time
       end do
       write(NameFile,'(a,i3.3,a)')'stoc',iPlot,'.out'
       call save_plot_file(NameFile=NameFile, &
            TypeFileIn='real8', TimeIn=tFinal, nStepIn = iStep, &
            NameVarIn='theta JT/2pi Log10VDF'  , &
            CoordMinIn_D=[         0.50*Delta, -2.0 + 0.50/nPointPer2Pi],&
            CoordMaxIn_D=[cTwoPi - 0.50*Delta,  2.0 - 0.50/nPointPer2Pi],&
            VarIn_II = log10(VDF_G(1:nTheta,1-nJ:nJ)) )
    end do
  end subroutine test_stochastic
  !============================================================================
end module ModStochastic
!==============================================================================
program test_program
  use ModTestPoissonBracket, ONLY: test_poisson_bracket, test_dsa_sa_mhd, &
       test_dsa_poisson, test_energy_conservation, test_in_action_angle,  &
       test_poisson_2d, test_poisson_2d_smooth, test_dsa_impl,            &
       test_semi_gyration, test_dsa_semi, test_semi_implicit_one_rotation,   &
       test_semi_implicit_two_rotation, test_semi_implicit_tvd_rotation,     &
       set_dsa_grid, TypeSolverDsa, nSweepDsa, IsTransposedDsa,           &
       nStepSaMhd, CflSaMhd, nStepImpl, CflImpl, DoHistoryDsa, LevelDsa
  use ModNumConst, ONLY: cTwoPi
  use ModUtilities, ONLY: CON_stop, upper_case, lower_case
  use ModHillVortex, ONLY: test_hill_vortex
  use ModStochastic, ONLY: test_stochastic
  use ModTestMultiPoisson, ONLY: test_multi_poisson

  implicit none
  ! ./test_poisson.exe [NameTest] [name=value ...]. Without NameTest or
  ! with explicitorig the nightly tests run; the tests and the options
  ! are described at their cases below
  character(len=20) :: NameTest = '', StringArg, NameGrid = 'nightly'
  integer :: iArg, iEqual, iError
  logical :: IsDsaVariant = .false.
  real    :: pMaxP0 = 1.0e4, FloorP0 = 1.0e-16
  integer :: nBinP0 = 20
  !----------------------------------------------------------------------------
  if(command_argument_count() > 0) call get_command_argument(1, NameTest)
  ! The test names are not case sensitive
  call lower_case(NameTest)
  do iArg = 2, command_argument_count()
     call get_command_argument(iArg, StringArg)
     iEqual = index(StringArg, '=')
     iError = 0
     select case(StringArg(1:max(iEqual-1,0)))
     case('solver')
        ! dsa: solver of test_dsa_semi, lGS (line Gauss-Seidel) or lGS1
        ! (the iterations of implicit2)
        TypeSolverDsa = StringArg(iEqual+1:)
        call upper_case(TypeSolverDsa)
        ! The line Gauss-Seidel names start with a lower-case l
        if(TypeSolverDsa(1:2) == 'LG') TypeSolverDsa(1:1) = 'l'
     case('sweeps')
        ! dsa: maximum number of the Gauss-Seidel sweeps
        read(StringArg(iEqual+1:), *, iostat=iError) nSweepDsa
     case('transposed')
        ! dsa: T puts the solver lines along the momentum
        read(StringArg(iEqual+1:), *, iostat=iError) IsTransposedDsa
     case('steps')
        ! Comprho4Shock: iterations of the explicit scheme (default 1000)
        read(StringArg(iEqual+1:), *, iostat=iError) nStepSaMhd
     case('cfl')
        ! Comprho4Shock: CFL of its local time steps (default 0.98)
        read(StringArg(iEqual+1:), *, iostat=iError) CflSaMhd
     case('cin')
        ! Comprho4Shock: CFL of the implicit steps (default 98)
        read(StringArg(iEqual+1:), *, iostat=iError) CflImpl
     case('isteps')
        ! Comprho4Shock: number of the implicit steps (default 1000)
        read(StringArg(iEqual+1:), *, iostat=iError) nStepImpl
     case('level')
        ! Comprho4Shock: the injection gives the analytic spectrum
        ! C (p0/p)^4 with C = level
        read(StringArg(iEqual+1:), *, iostat=iError) LevelDsa
     case('grid')
        ! Comprho4Shock: momentum grid nightly (default), wide or p0
        NameGrid = StringArg(iEqual+1:)
        if(all(NameGrid /= ['nightly', 'wide   ', 'p0     '])) iError = 1
     case('pmax')
        ! Comprho4Shock, grid=p0: upper end of the grid (default 1e4)
        read(StringArg(iEqual+1:), *, iostat=iError) pMaxP0
     case('bins')
        ! Comprho4Shock, grid=p0: bins per decade (default 20)
        read(StringArg(iEqual+1:), *, iostat=iError) nBinP0
     case('floor')
        ! Comprho4Shock, grid=p0: floor of f (default 1e-16, the analytic
        ! spectrum p^-4 at the default upper end)
        read(StringArg(iEqual+1:), *, iostat=iError) FloorP0
     case default
        iError = 1
     end select
     ! steps=, cfl=, cin=, isteps=, grid=, level=, pmax=, bins= and floor=
     ! belong to Comprho4Shock, the others to dsa
     if(any(StringArg(1:max(iEqual-1,0)) == ['steps ', 'cfl   ', &
          'cin   ', 'isteps', 'grid  ', 'level ', 'pmax  ', 'bins  ',  &
          'floor ']))then
        if(NameTest /= 'comprho4shock') iError = 1
     elseif(NameTest /= 'dsa')then
        iError = 1
     end if
     if(iError /= 0) call CON_stop( &
          'test_program: unknown or wrong argument '//trim(StringArg))
     IsDsaVariant = NameTest == 'dsa'
  end do
  write(*,*)' Start poisson bracket tests'

  select case(NameTest)
  case('rotation')
     ! Rotation of a Gaussian blob by a quarter turn with the schemes of
     ! the semi-implicit manuscript; takes a few minutes
     call test_semi_implicit_one_rotation(0.25)
     call test_semi_implicit_two_rotation(0.25)
     call test_semi_implicit_tvd_rotation(0.25)
  case('gyration')
     ! Gyration test of the semi-implicit scheme over one period
     call test_semi_gyration(cTwoPi)
  case('dsa')
     ! 100 momentum bins from 0.1 to 1e4 instead of 40 from 1 to 100, the
     ! injection still at p = 1; the floor is lowered to 1e-16, so that it
     ! is the same fraction of the p^-4 spectrum at the upper boundary.
     ! The index is fitted for 3 < p < 300 instead of 3 < p < 30.
     ! With solver=, sweeps= or transposed= only test_dsa_poisson, the
     ! initial state, runs before test_dsa_semi
     call set_dsa_grid(100, 0.1, 1.0e4, 1.0e-16, 2.5)
     call test_dsa_poisson
     if(.not.IsDsaVariant)then
        call test_dsa_sa_mhd
        call test_dsa_impl
     end if
     call test_dsa_semi
  case('comprho4shock')
     ! The steady shock acceleration of the manuscript (compression ratio
     ! 4): the explicit scheme with local time steps and implicit2, with
     ! their convergence histories. grid=nightly (default) keeps the 40
     ! bins from p = 1 to 100 of the nightly test
     select case(NameGrid)
     case('wide')
        ! 100 bins from p = 0.1 to 1e4, as in the dsa test
        call set_dsa_grid(100, 0.1, 1.0e4, 1.0e-16, 2.5)
     case('p0')
        ! bins= per decade from p = 1 to pmax=, the injection in the first
        ! bin, the floor floor=
        call set_dsa_grid(nint(nBinP0*log10(pMaxP0)), 1.0, pMaxP0, &
             FloorP0, 2.5)
     end select
     DoHistoryDsa = .true.
     call test_dsa_poisson
     call test_dsa_sa_mhd
     call test_dsa_impl
  case('', 'explicitorig')
     ! Nightly tests (make test_poisson_bracket), also explicitorig
     call test_poisson_bracket(cTwoPi)       ! nightly test1
     call test_poisson_2d(cTwoPi)            ! nightly test2
     call test_dsa_poisson                   ! nightly test3
     ! call test_poisson_2d_smooth(cTwoPi)
     ! call test_stochastic(1.2)
     ! call test_hill_vortex
     ! call test_energy_conservation(cTwoPi)
     ! call test_in_action_angle(cTwoPi)
     call test_dsa_sa_mhd                    ! for Fig 5.Right Panel
     call test_dsa_impl

     ! Tests of focused transport equation:
     ! Advection without any scattering
     ! call test_multi_poisson(100.0, TypeScatter='no')
     ! Advection and scattering
     ! call test_multi_poisson(100.0, TypeScatter='advect_scatter')
     ! Scatter only by diffusion
     ! call test_multi_poisson(100.0, TypeScatter='diffuse_scatter')
  case default
     call CON_stop('test_program: expected gyration, Comprho4Shock, '// &
          'dsa, explicitorig or rotation')
  end select
  write(*,*)' Finish poisson bracket tests'

end program test_program
!==============================================================================
