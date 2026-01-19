!  Copyright (C) 2002 Regents of the University of Michigan,
!  portions used with permission
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module SP_ModDiffusion

  ! Solve the diffusion term in the Parker equation
  ! Revision history:
  ! Prototype: SP/FLAMPA/src/ModDiffusion.f90
  ! Adapted for the use in MFLAMPA (Dist_I is an input paramater,
  ! fixed contributions to M_I in the end points)-D.Borovikov, 2017
  ! Updated (identation, comments): I.Sokolov, Dec.17, 2017

  use ModMpi
  use ModNumConst, ONLY: cPi, cTwoPi, cTiny
  use ModConst, ONLY: cAu, cLightSpeed, ckeV, cGeV, cMu, Rsun, cGyroRadius
  use SP_ModSize, ONLY: nVertexMax
  use SP_ModDistribution, ONLY: SpeedSi_G, Momentum_G, GammaLorentz_G, &
       MomentumInjSi, Mu_F, DeltaMu, Distribution_CB
  use SP_ModGrid, ONLY: nP, nMu, State_VIB, MHData_VIB, &
       D_, R_, X_, Z_, Wave1_, Wave2_
  use SP_ModProc, ONLY: nProc, iProc, iComm, iError
  use SP_ModUnit, ONLY: Io2Si_V, UnitX_
  use ModUtilities, ONLY: CON_stop

  implicit none

  PRIVATE ! Except

  SAVE

  ! Public members:
  public :: read_param, diffuse_distribution, scatter_distribution, &
       diffuseperp_distribution, set_diffusion_coef
  ! Diffusion in space
  interface diffuse_distribution
     module procedure diffuse_distribution_s   ! Global time step Dt
     module procedure diffuse_distribution_arr ! DtLocal for (nP, nX)
     module procedure diffuse_distribution_c   ! DtLocal for (nP, nMu, nX)
  end interface diffuse_distribution
  ! Diffusion by pitch angle scattering
  interface scatter_distribution
     module procedure scatter_distribution_s   ! Global time step Dt
     module procedure scatter_distribution_c   ! DtLocal for (nP, nMu, nX)
  end interface scatter_distribution

  ! Whether to include diffusion term (parallel comes first)
  logical, public :: UseDiffusion = .true.
  logical, public :: UseDiffusionPerp = .false. ! Perpendicular Diffusion
  ! Determine which type of the diffusion in upstream and from mhd turbulence
  character(len=15) :: TypeMhdDiffusion = 'sokolov2004'
  ! Set diffusion or scatter coefficients in the diffusion operator
  ! df/dt = DOuter * d(DInner * df/dx)/dx
  ! DOuter = BSi in the cell center
  ! DInner = Diffusion Coefficient/BSi at the face
  real, public, dimension(nVertexMax, nP) :: DInnerSi_II, & ! Dinner/BSi
       CoefLambdaMuMu_II ! mfp, mu>1
  real, public, dimension(nVertexMax) :: DOuterSi_I, CoefLambdaxx_I ! mfp, mu=1

  ! Perpendicular diffusion
  real,    public :: DPerpRatio = 0.065         ! Simple ratio of Dperp/Dpara
  ! Grid of the triangulated mesh for solving perp. term: nR * nTheta * nPhi
  integer         :: nRPerp, nThetaPerp, nPhiPerp
  real            :: dThetaPerp, dPhiPerp       ! dTheta and dPhi
  real, allocatable:: dRPerpMesh_I(:), dRPerpFace_I(:) ! dR of cell center+face
  real            :: RMinPerp = 1.2, RMaxPerp = 240.0  ! RMin and RMax of Grid
  real            :: dLogRFacePerp = 0.0        ! Geometric Sequence for RMesh
  character(len=6):: ScaleRPerp                 ! Scale (Linear/Log) along R_
  integer, allocatable :: iRPerpStart_I(:), iRPerpEnd_I(:)
  real, allocatable:: RPerp_C(:), ThetaPerp_C(:), PhiPerp_C(:), &
       RPerp_F(:), ThetaPerp_F(:), PhiPerp_F(:), XyzPerp_CB(:,:,:,:)

  ! Local parameters!
  ! Diffusion as in Li et al. 2003, doi:10.1029/2002JA009666
  logical         :: UseBtotal = .false., UseFixedUps = .false.
  real            :: MeanFreePath0InAu = 1.0
  ! Parameter characterizing cut-off wavenumber of turbulent spectrum:
  ! value of scale turbulence at 1 AU for any type (const or linear)
  real            :: ScaleTurbulenceSi = 0.03*cAu
  integer, parameter :: Const_ = 0, Linear_ = 1
  integer         :: iScaleTurbulenceType = Const_
contains
  !============================================================================
  subroutine read_param(NameCommand)
    use ModReadParam, ONLY: read_var
    use ModUtilities, ONLY: lower_case
    character(len=*), intent(in) :: NameCommand ! From PARAM.in
    character(len=8) :: TypeScaleTurbulence
    character(len=*), parameter:: NameSub = 'read_param'
    !--------------------------------------------------------------------------
    select case(NameCommand)
    case('#SCALETURBULENCE')
       ! cut-off wavenumber of turbulence spectrum: k0 = 2 cPi / Scale
       call read_var('ScaleTurbulenceType', TypeScaleTurbulence)
       call lower_case(TypeScaleTurbulence)
       ! determine the turbulence scale wrt distance in space
       select case(trim(TypeScaleTurbulence))
       case('const', 'constant')
          iScaleTurbulenceType = Const_
       case('linear')
          iScaleTurbulenceType = Linear_
       case default
          call CON_stop(NameSub//": Unknown scale turbulence type")
       end select
       ! specify the turbulence scale at 1 au
       call read_var('ScaleTurbulence1AU', ScaleTurbulenceSi)
       ScaleTurbulenceSi = ScaleTurbulenceSi * cAu
    case('#DIFFUSIONPARA')
       ! Use calculated from the mhd turbulence for downstream; if not using
       ! the fixed upstream diffusion coefficient, also use the mhd turbulence
       ! Determine whether or not to use Btotal=sqrt(B**2+dB**2) as B
       call read_var('UseBtotal', UseBtotal)
       call read_var("TypeMhdDiffusion", TypeMhdDiffusion)
       call lower_case(TypeMhdDiffusion)
    case('#DIFFUSIONPERP')
       call read_var('UseDiffusionPerp', UseDiffusionPerp)
       ! Setup perpendicular diffusion coefficients only when using it
       if(UseDiffusionPerp) then
          ! Simply assume DPerp/DPara = DPerpRatio
          call read_var('DPerpRatio', DPerpRatio)

          ! Handle triangulated mesh in radial direction
          call read_var('nRPerp', nRPerp)
          call read_var('RMinPerp', RMinPerp)
          call read_var('RMaxPerp', RMaxPerp)
          if(nRPerp<=1 .or. RMinPerp<=0.0 .or. RMaxPerp<RMinPerp) RETURN
          if(nRPerp<nProc) then
             write(*,*) "nRPerp = ", nRPerp, "< nProc = ", nProc, &
                  "now set nRPerp = nProc = ", nProc
          end if

          ! From RMinPerp to RMaxPerp: Distribution of multiple spheres
          call read_var('ScaleRPerp', ScaleRPerp)

          ! Handle lon-lat grid in the triangulated mesh
          call read_var('nThetaPerp', nThetaPerp)
          call read_var('nPhiPerp', nPhiPerp)
          if(nThetaPerp<=0 .or. nPhiPerp<=0) RETURN
          ! Theta: nThetaPerp, dThetaPerp and ThetaPerp_C
          if(mod(nThetaPerp, 2) > 0) then
             write(*,*) "nThetaPerp reads", nThetaPerp, "and becomes", &
                  nThetaPerp+1, "to avoid singularity when theta=0"
             nThetaPerp = nThetaPerp + 1
          end if

          ! Setup the mesh used for triangulation in perpendicular diffusion
          call setup_multi_uniform_spheres
       end if
    case('#USEFIXEDMFPUPSTREAM')
       ! Determine whether or not to fix upstream diffusion coefficient
       call read_var('UseFixedMeanFreePathUpstream',&
            UseFixedUps)
       if(UseFixedUps)then
          ! see Li et al. (2003), doi:10.1029/2002JA009666
          call read_var('MeanFreePath0InAu', MeanFreePath0InAu)
       end if
    case('#DIFFUSION')
       ! To be used only for testing/developing, to switch off diffusion
       call read_var('UseDiffusion', UseDiffusion)
    case default
       call CON_stop('SP:'//NameSub//': unknown command '//NameCommand)
    end select
  end subroutine read_param
  !============================================================================
  subroutine diffuse_distribution_s(iLine, nX, iShock, nSi_I, &
       BSi_I, Dt, LowerEndSpectrumIn_I, UpperEndSpectrumIn_I, &
       LowerEndSpectrumIn_II, UpperEndSpectrumIn_II)
    ! diffuse the distribution function with scalar/global Dt, with
    ! the pitch-angle-averaged or dependent lower and upper spectra

    ! Variables as inputs
    ! input Line, End (for how many particles), and Shock indices
    integer, intent(in) :: iLine, nX, iShock
    ! Number density and magnetic field strength at the end of this iteration
    real, intent(in) :: nSi_I(nX), BSi_I(nX)
    ! Scalar time step for diffusion
    real, intent(in) :: Dt
    ! Given spectrum at low end (flare acceleration) and upper end (GCRs)
    real, intent(in), optional, dimension(nP) :: &
         LowerEndSpectrumIn_I, UpperEndSpectrumIn_I   ! Mu-averaged
    real, intent(in), optional, dimension(nP, nMu) :: &
         LowerEndSpectrumIn_II, UpperEndSpectrumIn_II ! Mu-dependent
    ! LOCAL VARS
    real :: DtFake_C(nP, nMu, nX)
    !--------------------------------------------------------------------------
    DtFake_C = Dt

    call diffuse_distribution_c(iLine, nX, &
         iShock, nSi_I, BSi_I, DtFake_C, &
         LowerEndSpectrumIn_I, UpperEndSpectrumIn_I, &
         LowerEndSpectrumIn_II, UpperEndSpectrumIn_II)

  end subroutine diffuse_distribution_s
  !============================================================================
  subroutine diffuse_distribution_arr(iLine, nX, iShock, nSi_I, &
       BSi_I, DtLocal_II, LowerEndSpectrumIn_I, UpperEndSpectrumIn_I, &
       LowerEndSpectrumIn_II, UpperEndSpectrumIn_II)
    ! diffuse the distribution function with vector/local Dt, with the
    ! pitch-angle-averaged local time step, and lower and upper spectra

    ! Variables as inputs
    ! input Line, End (for how many particles), and Shock indices
    integer, intent(in) :: iLine, nX, iShock
    ! Number density and magnetic field strength at the end of this iteration
    real, intent(in) :: nSi_I(nX), BSi_I(nX)
    ! Local time step for diffusion: Only with nP and nX
    real, intent(in) :: DtLocal_II(nP, nX)
    ! Given spectrum at low end (flare acceleration) and upper end (GCRs)
    real, intent(in), optional, dimension(nP) :: &
         LowerEndSpectrumIn_I, UpperEndSpectrumIn_I   ! Mu-averaged
    real, intent(in), optional, dimension(nP, nMu) :: &
         LowerEndSpectrumIn_II, UpperEndSpectrumIn_II ! Mu-dependent
    ! LOCAL VARS
    real :: DtFake_C(nP, nMu, nX)
    !--------------------------------------------------------------------------
    DtFake_C = spread(DtLocal_II, DIM=2, NCOPIES=nMu)

    call diffuse_distribution_c(iLine, nX, &
         iShock, nSi_I, BSi_I, DtFake_C, &
         LowerEndSpectrumIn_I, UpperEndSpectrumIn_I, &
         LowerEndSpectrumIn_II, UpperEndSpectrumIn_II)

  end subroutine diffuse_distribution_arr
  !============================================================================
  subroutine diffuse_distribution_c(iLine, nX, iShock, nSi_I, BSi_I, &
       DtLocalIn_C, LowerEndSpectrumIn_I, UpperEndSpectrumIn_I, &
       LowerEndSpectrumIn_II, UpperEndSpectrumIn_II)
    ! diffuse the distribution function with vector/local Dt, with the
    ! pitch-angle-dependent local time step, and lower and upper spectra

    use ModMpi
    use SP_ModTurbulence, ONLY: UseTurbulentSpectrum, set_dxx, Dxx

    ! Variables as inputs
    ! Input Line, End (for how many particles), and Shock indices
    integer, intent(in) :: iLine, nX, iShock
    ! Number density and magnetic field strength at the end of this iteration
    real, intent(in) :: nSi_I(nX), BSi_I(nX)
    ! Local time step for diffusion: with nP, nMu and nX
    real, intent(in) :: DtLocalIn_C(nP, nMu, nX)
    ! Given spectrum at low end (flare acceleration) and upper end (GCRs)
    real, optional, intent(in), dimension(nP) :: &
         LowerEndSpectrumIn_I, UpperEndSpectrumIn_I   ! Mu-averaged
    real, optional, intent(in), dimension(nP, nMu) :: &
         LowerEndSpectrumIn_II, UpperEndSpectrumIn_II ! Mu-dependent
    ! Variables declared in this subroutine
    integer :: iP, iMu                                ! Loop variables
    ! For an optimized loop, we need to change the index order by
    ! (nP, nMu, nX) into by (nX, nP, nMu) for the
    ! three-dimensional local time step DtLocal_III array
    real :: DtLocal_C(nX, nP, nMu)                    ! Local time stepping
    ! Same reason for LowerEndSpectrumIn_II and UpperEndSpectrumIn_II: keep
    ! (nP, nMu) as a good looping order, outermost=iMu, innermost=iP
    real, dimension(nP, nMu) :: LowerEndSpectrum_II, UpperEndSpectrum_II
    ! Logical variable for whether or not to use given lower and upper spectra
    logical :: UseLowerSpectrum = .false., UseUpperSpectrum = .false.
    real, parameter :: DiffCoeffMinSi = 1.0E+04*Rsun
    ! Mesh spacing and face spacing.
    real :: DsMeshSi_I(1:nX-1), DsFaceSi_I(2:nX-1)
    ! Main, upper, and lower diagonals, source
    real :: Main_I(nX), Upper_I(nX), Lower_I(nX), Res_I(nX)
    ! Auxiliary arrays for the tri-diagonal arrays
    real :: Aux1_I(nX), Aux2_I(nX)
    !--------------------------------------------------------------------------
    ! diffusion along the field line
    ! Set up the local time step: the reason is that, we loop through each
    ! iX at fixed iP and iMu, so we will visit each DtLocal_C(iX, iP, iMu)
    ! in the loop with the inner-most index iX, intermediate loop index iP,
    ! and the outer-most index iMu. It will visit and put the data in cache
    ! that are stored adjacently.
    do iMu = 1, nMu
       DtLocal_C(:, :, iMu) = transpose(DtLocalIn_C(:, iMu, :))
    end do

    ! Get LowerEndSpectrum_II and/or UpperEndSpectrum_II if given
    ! Given Mu-averaged arrays: spread into 2 dimensions with nMu
    if(present(LowerEndSpectrumIn_I)) then
       UseLowerSpectrum = .true.
       LowerEndSpectrum_II = spread(LowerEndSpectrumIn_I, DIM=2, NCOPIES=nMu)
    end if
    if(present(UpperEndSpectrumIn_I)) then
       UseUpperSpectrum = .true.
       UpperEndSpectrum_II = spread(UpperEndSpectrumIn_I, DIM=2, NCOPIES=nMu)
    end if
    ! Given Mu-dependent arrays: transpose with a better index for the loop
    if(present(LowerEndSpectrumIn_II)) then
       UseLowerSpectrum = .true.
       LowerEndSpectrum_II = LowerEndSpectrumIn_II
    end if
    if(present(UpperEndSpectrumIn_II)) then
       UseUpperSpectrum = .true.
       UpperEndSpectrum_II = UpperEndSpectrumIn_II
    end if

    ! if using turbulent spectrum: set_dxx for diffusion along the field line
    if(UseTurbulentSpectrum) call set_dxx(nX, BSi_I(1:nX))

    ! In M-FLAMPA DsMeshSi_I(i) is the distance between centers of meshes
    ! i and i+1. Therefore,
    DsMeshSi_I(1:nX-1) = max(State_VIB(D_,1:nX-1,iLine)*Io2Si_V(UnitX_), cTiny)
    ! Within the framework of finite volume method, the cell volume
    ! is used, which is proportional to the distance between the faces
    ! bounding the volume with an index, i, which is half of sum of
    ! distance between meshes i-1 and i, i.e. DsMeshSi_I(i-1), and that
    ! between meshes i and i+1, which is DsMeshSi_I(i):
    DsFaceSi_I(2:nX-1) = 0.5*(DsMeshSi_I(1:nX-2)+DsMeshSi_I(2:nX-1))

    ! In flux coordinates, the control volume associated with the given
    ! cell has a cross-section equal to (Magnetic Flux)/B, where the flux
    ! is a constant along the magnetic field line, set to one hereafter.
    ! Therefore, the diffusion equation has a following form:
    ! (DsFaceSi_I/BSi_I)*(f^(n+1) - f^n) = Flux_(i-1/2) - Flux_(i+1/2),
    ! where the particle density flux should be multiplied by the
    ! cross-section area too (magnetic flux factor is one!):
    !  Flux_(i-1/2) = (diffusion coefficient)_(i-1/2)/B_(i-1/2)*&
    !                 (f^(n+1)_(i-1) - f^(n+1)_i),
    !  The multiplier, DsFaceSi_I/BSi_I, is denoted as DsFaceSi_I/DOuter_I
    !  The face-centered combination,

    ! For each momentum, the dependence of the diffusion coefficient
    ! on momentum is D \propto r_L*v \propto Momentum**2/TotalEnergy
    if(UseTurbulentSpectrum) then
       DInnerSi_II(1:nX, 1:nP) = Dxx(nX, BSi_I(1:nX))/ &
            spread(BSi_I(1:nX), DIM=2, NCOPIES=nP)
    else
       ! Lower limit to floor the spatial diffusion coefficient For a
       ! given spatial and temporal resolution, the value of the
       ! diffusion coefficient should be artificially increased to get
       ! the diffusion length to be larger than the shock front width,
       ! which even for the steepened shock is as wide as a mesh size
       ! of the Largangian grid, State_VIB(D_,:,:)*RSun. In this way,
       ! the Lagrangian grid resolution is sufficient to resolve a
       ! precursor in the upstream distribution of DiffCoeffMin=0
       ! (default value), we do NOT enhance the diffusion coefficient!
       ! Physically, DiffCoeffMin should be given by the product of
       ! shock wave speed and local grid spacing.
       DInnerSi_II(1:nX, 1:nP) = max(DInnerSi_II(1:nX, 1:nP), &
            DiffCoeffMinSi/spread(DOuterSi_I(1:nX), DIM=2, NCOPIES=nP))
    end if

    MU:do iMu = 1, nMu
       MOMENTUM:do iP = 1, nP
          ! Now, we solve the matrix equation
          ! f^(n+1)_i-Dt*DOuter_I/DsFaceSi_I*(&
          !     DInner_(i+1/2)*(f^(n+1)_(i+1)-f^(n+1)_i)/DsMesh_(i+1)-&
          !     DInner_(i-1/2)*(f^(n+1)_i -f^(n+1)_(i-1)/DsMesh_(i ))=f^n_i
          ! Set source term in the RHS:
          Res_I = Distribution_CB(iP, iMu, 1:nX, iLine)
          ! Set elements of tri-diagonal matrix in the LHS
          Main_I = 1.0; Lower_I = 0.0; Upper_I = 0.0

          ! For i=1:
          Aux1_I(1)  = DtLocal_C(1,iP,iMu)*DOuterSi_I(1)*   &
               0.5*(DInnerSi_II(1,iP) + DInnerSi_II(2,iP))/DsMeshSi_I(1)**2
          Main_I(1)  = Main_I(1) + Aux1_I(1)
          Upper_I(1) = -Aux1_I(1)
          if(UseLowerSpectrum) then
             ! With the given value for f_0 behind the boundary,
             ! the above scheme reads:
             ! f^(n+1)_1-Dt*DOuter_I/DsFaceSi_I*(&
             !     DInner_(3/2)*(f^(n+1)_2-f^(n+1)_1)/DsMesh_(2)-&
             !     DInner_(1/2)*(f^(n+1)_1 -f_0/DsMesh_(1))=f^n_1
             Aux2_I(1) = DtLocal_C(1,iP,iMu)*DOuterSi_I(1) &
                  *DInnerSi_II(1,iP)/DsMeshSi_I(1)**2
             ! With these regards, Aux2 is added to Main_I(1)...
             Main_I(1) = Main_I(1) + Aux2_I(1)
             ! ...while the given Aux2*f_0 is moved to the RHS and summed up
             ! with the source:
             Res_I(1)  = Res_I(1) + Aux2_I(1)*LowerEndSpectrum_II(iP, iMu)
          end if

          ! For i=2, n-1:
          Aux1_I(2:nX-1) = DtLocal_C(2:nX-1,iP,iMu)*DOuterSi_I(2:nX-1)* &
               0.5*(DInnerSi_II(2:nX-1,iP) + DInnerSi_II(3:nX, iP))/ &
               (DsMeshSi_I(2:nX-1)*DsFaceSi_I(2:nX-1))
          Aux2_I(2:nX-1) = DtLocal_C(2:nX-1,iP,iMu)*DOuterSi_I(2:nX-1)* &
               0.5*(DInnerSi_II(1:nX-2,iP) + DInnerSi_II(2:nX-1,iP))/&
               (DsMeshSi_I(1:nX-2)*DsFaceSi_I(2:nX-1))
          Main_I(2:nX-1)  = Main_I(2:nX-1) + Aux1_I(2:nX-1) + Aux2_I(2:nX-1)
          Upper_I(2:nX-1) = -Aux1_I(2:nX-1)
          Lower_I(2:nX-1) = -Aux2_I(2:nX-1)

          ! For i=n:
          ! Version before Nov. 2023:
          ! Aux2 = Dt*DOuter_I(n)*0.5*(DInner_I(n-1) + DInner_I(n))/ &
          !     DsMeshSi_I(n-1)**2
          !
          ! After Nov. 2023: set free escaping at outer boundary for now
          ! Aux2 = 0.0
          ! In both these versions:
          ! Main_I( n) = Main_I(n) + Aux2
          ! Lower_I(n) = -Aux2
          ! So, effectively for the version after Nov. 2023
          ! Main_I(n) = 1; Lower_I(n) = 0 (equivalently to doing nothing)
          ! For backward compatibility, please keep UseUpperEndBc=.false.
          if(UseUpperSpectrum) then
             Aux2_I(nX)  = DtLocal_C(nX,iP,iMu)*DOuterSi_I(nX)*0.5* &
                  (DInnerSi_II(nX-1,iP)+DInnerSi_II(nX,iP))/DsMeshSi_I(nX-1)**2
             Main_I(nX)  = Main_I(nX) + Aux2_I(nX)
             Lower_I(nX) = -Aux2_I(nX)
             Aux1_I(nX)  = DtLocal_C(nX,iP,iMu)*DOuterSi_I(nX)* &
                  DInnerSi_II(nX,iP)/DsMeshSi_I(nX-1)**2
             Main_I(nX)  = Main_I(nX) + Aux1_I(nX)
             Res_I(nX)   = Res_I(nX) + Aux1_I(nX)*UpperEndSpectrum_II(iP, iMu)
          end if

          ! Update the solution from f^(n) to f^(n+1):
          call tridiag(nX, Lower_I, Main_I, Upper_I, Res_I, &
               Distribution_CB(iP, iMu, 1:nX, iLine))
       end do MOMENTUM
    end do MU

  end subroutine diffuse_distribution_c
  !============================================================================
  subroutine setup_multi_uniform_spheres

    use ModCoordTransform, ONLY: sph_to_xyz
    integer :: iRPerp, iThetaPerp, iPhiPerp ! loop variables
    integer :: iProcChunk, iRPerpRemainder
    character(len=*), parameter:: NameSub = 'setup_multi_uniform_spheres'
    !--------------------------------------------------------------------------
    ! setup the mesh arrays for triangulation used in perpendicular diffusion
    if(allocated(RPerp_C)) deallocate(RPerp_C)
    allocate(RPerp_C(nRPerp))
    if(allocated(RPerp_F)) deallocate(RPerp_F)
    allocate(RPerp_F(nRPerp+1))
    if(allocated(dRPerpMesh_I)) deallocate(dRPerpMesh_I)
    allocate(dRPerpMesh_I(nRPerp-1))
    if(allocated(dRPerpFace_I)) deallocate(dRPerpFace_I)
    allocate(dRPerpFace_I(nRPerp))
    if(allocated(XyzPerp_CB)) deallocate(XyzPerp_CB)
    allocate(XyzPerp_CB(X_:Z_, nPhiPerp, nThetaPerp, nRPerp))
    XyzPerp_CB = 0.0

    ! Determine chunk size, start and end indices for each processor
    if(allocated(iRPerpStart_I)) deallocate(iRPerpStart_I)
    if(allocated(iRPerpEnd_I)) deallocate(iRPerpEnd_I)
    allocate(iRPerpStart_I(0:nProc-1), iRPerpEnd_I(0:nProc-1))
    iProcChunk = nRPerp/nProc             ! Base chunk size
    iRPerpRemainder = mod(nRPerp, nProc)  ! Leftover elements
    ! Loop over processors to compute start and end indices
    do iProc = 0, nProc-1
       if(iProc < iRPerpRemainder) then
          iRPerpStart_I(iProc) = iProc*(iProcChunk+1) + 1
          iRPerpEnd_I(iProc)   = iRPerpStart_I(iProc) + iProcChunk
       else
          iRPerpStart_I(iProc) = iProc*iProcChunk + iRPerpRemainder+1
          iRPerpEnd_I(iProc)   = iRPerpStart_I(iProc) + iProcChunk-1
       endif
    end do

    ! R: radial direction
    select case(trim(ScaleRPerp))
    case("Linear", "linear")
       dRPerpMesh_I = (RMaxPerp - RMinPerp)/real(nRPerp) ! Same=Const.
       dRPerpFace_I = dRPerpMesh_I(1) ! Same=Const.
       do iRPerp = 1, nRPerp
          RPerp_F(iRPerp) = RMinPerp + (real(iRPerp)-1.0)*dRPerpMesh_I(1)
       end do
       RPerp_F(nRPerp+1) = RMaxPerp
       RPerp_C = 0.5*(RPerp_F(1:nRPerp) + RPerp_F(2:nRPerp+1))
    case("Exp", "exp", "Exponential", "exponential")
       dLogRFacePerp = log(RMaxPerp/RMinPerp)/real(nRPerp)
       do iRPerp = 1, nRPerp+1
          RPerp_F(iRPerp) = RMinPerp * exp(iRPerp*dLogRFacePerp)
       end do
       RPerp_C = 0.5*(RPerp_F(1:nRPerp) + RPerp_F(2:nRPerp+1))
       dRPerpFace_I = RPerp_F(2:nRPerp+1) - RPerp_F(1:nRPerp)
       dRPerpMesh_I = RPerp_C(2:nRPerp) - RPerp_C(1:nRPerp-1)
    case default
       call CON_stop(NameSub// &
            'Error in ScaleRPerp for Perpendicular Diffusion.')
    end select

    ! Theta: zenith/latitudinal direction
    dThetaPerp = cPi/real(nThetaPerp)
    if(allocated(ThetaPerp_C)) deallocate(ThetaPerp_C)
    allocate(ThetaPerp_C(nThetaPerp))
    do iThetaPerp = 1, nThetaPerp
       ThetaPerp_C(iThetaPerp) = (real(iThetaPerp)-0.5)*dThetaPerp
    end do

    ! Phi: azimuthal/longitudinal direction
    dPhiPerp = cTwoPi/real(nPhiPerp)
    if(allocated(PhiPerp_C)) deallocate(PhiPerp_C)
    allocate(PhiPerp_C(nPhiPerp))
    do iPhiPerp = 1, nPhiPerp
       PhiPerp_C(iPhiPerp) = (real(iPhiPerp)-0.5)*dPhiPerp
    end do

    ! (R, Theta, Phi) => (X, Y, Z) + parallelization
    do iRPerp = iRPerpStart_I(iProc), iRPerpEnd_I(iProc)
       do iThetaPerp = 1, nThetaPerp
          do iPhiPerp = 1, nPhiPerp
             call sph_to_xyz(RPerp_C(iRPerp), &
                  ThetaPerp_C(iThetaPerp), PhiPerp_C(iPhiPerp), &
                  XyzPerp_CB(:, iPhiPerp, iThetaPerp, iRPerp))
          end do
       end do
    end do
    if(nProc > 1) call MPI_ALLREDUCE(MPI_IN_PLACE, XyzPerp_CB, &
         3*nPhiPerp*nThetaPerp*nRPerp, MPI_REAL, MPI_SUM, iComm, iError)

  end subroutine setup_multi_uniform_spheres
  !============================================================================
  subroutine diffuseperp_distribution()

    use SP_ModGrid, ONLY: nLineAll
    use SP_ModTriangulate, ONLY: intersect_surf, build_trmesh, &
         interpolate_trmesh

    integer :: iRPerp, iThetaPerp, iPhiPerp, iPE ! loop variables
    real    :: DistrRPerp_5D(0:nP+1, 1:nMu, nPhiPerp, nThetaPerp, nRPerp)
    character(len=*), parameter:: NameSub = 'diffuseperp_distribution'
    !--------------------------------------------------------------------------
    ! cross-field (perpendicular) diffusion
    DistrRPerp_5D = 0.0

    ! step 1: field lines => intersection points on multiple uniform layers
    ! here, iProc is for field lines, not for sub-slices/layers
    do iPE = 0, nProc-1
       do iRPerp = iRPerpStart_I(iPE), iRPerpEnd_I(iPE)
          call intersect_surf(RPerp_C(iRPerp), iPE, iRperp)
       end do
    end do

    ! then, iProc is for sub-slices/layers now
    do iRPerp = iRPerpStart_I(iProc), iRPerpEnd_I(iProc)
       ! step 2: intersection points => construct the triangulation skeleton
       call build_trmesh(iRPerp)
       do iThetaPerp = 1, nThetaPerp
          do iPhiPerp = 1, nPhiPerp
             call interpolate_trmesh(XyzPerp_CB(:,iPhiPerp,iThetaPerp,iRPerp),&
                  iRIn = iRPerp, Log10DistrInterp_II=                         &
                  DistrRPerp_5D(:,:,iPhiPerp,iThetaPerp,iRPerp))
          end do
       end do
    end do

    ! Broadcast results to all processes
    if(nProc > 1) then
      do iPE = 0, nProc-1
       call MPI_BCAST(DistrRPerp_5D(:,:,:,:,             &
            iRPerpStart_I(iPE):iRPerpEnd_I(iPE)),    &
            (iRPerpEnd_I(iPE)-iRPerpStart_I(iPE)+1)* &
            (nP+2)*nMu*nPhiPerp*nThetaPerp, MPI_REAL, iPE, iComm, iError)
      end do
    end if
    DistrRPerp_5D = exp(DistrRPerp_5D*log(10.0)) ! log10(VDF) => VDF
    ! Check the error message from after interpolate_trmesh
    if(iError /= 0) then
       write(*,*) 'iProc = ', iProc, &
            NameSub//': interpolate_trmesh failed, Dperp stopped'
       RETURN
    end if

    ! step 4: solve the Dperp distribution Eq in multiple uniform layers

    ! step 5: interpolate back to the intersection points along field lines

  end subroutine diffuseperp_distribution
  !============================================================================
  subroutine scatter_distribution_s(iLine, nX, nSi_I, BSi_I, Dt)
    ! Calculate scattering: \deltaf/\deltat = (Dmumu*f_mu)_mu with global Dt

    ! Variables as inputs
    ! input Line and End index (for how many particles)
    integer, intent(in) :: iLine, nX
    ! Number density and magnetic field strength at the end of this iteration
    real, intent(in) :: nSi_I(nX), BSi_I(nX)
    ! Scalar time step for diffusion
    real, intent(in) :: Dt
    ! LOCAL VARs
    real :: DtFake_C(nP, nMu, nX)
    !--------------------------------------------------------------------------
    DtFake_C = Dt

    call scatter_distribution_c(iLine, nX, nSi_I, BSi_I, DtFake_C)

  end subroutine scatter_distribution_s
  !============================================================================
  subroutine scatter_distribution_c(iLine, nX, nSi_I, BSi_I, DtLocal_C)
    ! Calculate scattering: \deltaf/\deltat = (Dmumu*f_mu)_mu with local Dt

    use ModMpi
    use ModConst, ONLY: cAtomicMass

    ! Variables as inputs
    ! input Line and End index (for how many particles)
    integer, intent(in) :: iLine, nX
    ! Number density and magnetic field strength at the end of this iteration
    real, intent(in) :: nSi_I(nX), BSi_I(nX)
    ! Scalar time step for diffusion
    real, intent(in) :: DtLocal_C(nP, nMu, nX)
    ! LOCAL VARs
    ! For each fixed iMu and iX, we calculate Dmumu, the coefficient of
    ! diffusion about the pitch angle = v/lambda_mumu*(1-mu**2)*mu**(2.0/3.0)
    real    :: DMuMu_I(0:nMu)
    ! Factorize DMuMu, each factor being only a function of p**3/3, mu, s_L
    real    :: FactorMu_F(0:nMu)
    ! Lower, main, upper diagonals, and source for the scatter calculation
    real    :: Main_I(nMu), Upper_I(nMu), Lower_I(nMu), Res_I(nMu)
    ! Physical VARs
    real    :: LowerLimMu           ! Lower limit of Factor_mu
    real    :: DtOverDMu2_C(nP, 0:nMu, nX) ! Dt/DeltaMu**2
    real    :: AlfvenSpeed_I(nX)    ! Alfven wave speed
    integer :: iP, iX               ! Loop variables
    integer :: iMuSwitch            ! Control parameter for mu
    !--------------------------------------------------------------------------
    DtOverDMu2_C(:, 1:nMu, :) = DtLocal_C/DeltaMu**2 ! Dt/DeltaMu**2
    DtOverDMu2_C(:, 0, :) = DtOverDMu2_C(:, 1, :)    ! Bc

    ! Calculate factorized diffusion coefficient: (1-mu**2) * mu**(2.0/3.0)
    LowerLimMu = (1.0 - DeltaMu**2)*abs(DeltaMu)**(2.0/3.0)
    FactorMu_F = (1.0 - Mu_F**2)*abs(Mu_F)**(2.0/3.0)

    ! Get the Alfven wave speed for each s_L
    AlfvenSpeed_I = BSi_I/sqrt(cMu*nSi_I*cAtomicMass)

    ! Calculate the effect of scatter along mu axis for each fixed iX and iP
    SPACELOC:do iX = 1, nX
       MOMENTUM:do iP = 1, nP
          ! For each pitch angle, set DMuMu and solve VDF for mu scattering
          ! Meanwhile, we control whether we floor the value of D_mumu or not
          iMuSwitch = 0
          if(.not.any(SpeedSi_G(iP)*abs(Mu_F) >= 10.0*AlfvenSpeed_I(iX))) &
               iMuSwitch = minloc(Mu_F, DIM=1, MASK= &
               SpeedSi_G(iP)*abs(Mu_F) < 10.0*AlfvenSpeed_I(iX))

          where(SpeedSi_G(iP)*abs(Mu_F) >= 10.0*AlfvenSpeed_I(iX))
             ! Set Dmumu where mu is large enough
             DMuMu_I = SpeedSi_G(iP)/(CoefLambdaMuMu_II(iX, iP)* &
                  Momentum_G(iP))*FactorMu_F*DtOverDMu2_C(iP, :, iX)
          elsewhere
             ! Set the lower limit of D_mumu when |mu| is close to zero
             DMuMu_I = SpeedSi_G(iP)/(CoefLambdaMuMu_II(iX, iP)* &
                  Momentum_G(iP))*max(FactorMu_F(iMuSwitch), LowerLimMu)* &
                  DtOverDMu2_C(iP,:,iX)
          end where

          ! Set up coefficients for solving the linear equation set
          Lower_I = -DMuMu_I(0:nMu-1)
          Upper_I = -DMuMu_I(1:nMu  )
          Main_I  = 1.0 - Lower_I - Upper_I
          Res_I   = Distribution_CB(iP, 1:nMu, iX, iLine)

          ! Update the solution for mu scattering
          call tridiag(nMu, Lower_I, Main_I, Upper_I, Res_I, &
               Distribution_CB(iP, 1:nMu, iX, iLine))
       end do MOMENTUM
    end do SPACELOC

  end subroutine scatter_distribution_c
  !============================================================================
  subroutine set_diffusion_coef(iLine, nX, iShock, BSi_I)
    ! Set spatial diffusion and mu scattering coefficients for given field line

    use ModConst, ONLY: cElectronCharge, cProtonMass
    integer, intent(in):: iLine, nX, iShock
    integer            :: iP
    real, intent(in)   :: BSi_I(1:nX)
    real               :: ScaleSi_I(1:nX), RadiusSi_I(1:nX), OmegaCyclo_I(1:nX)
    real               :: dBSi_I(1:nX), BtotalSi_I(1:nX)
    real, parameter    :: cCoef = 81.0/(7.0*cPi*cTwoPi**(2.0/3))
    real, parameter    :: cCoefxx_to_mumu = 14.0/27.0
    character(len=*), parameter:: NameSub = 'set_diffusion_coef'
    !--------------------------------------------------------------------------
    DOuterSi_I(1:nX) = BSi_I(1:nX)
    RadiusSi_I(1:nX) = State_VIB(R_,1:nX,iLine)*Io2Si_V(UnitX_)
    ! if(UseTurbulentSpectrum) RETURN

    ! precompute scale of turbulence along the line
    select case(iScaleTurbulenceType)
    case(Const_)
       ScaleSi_I(1:nX) = ScaleTurbulenceSi
    case(Linear_)
       ScaleSi_I(1:nX) = ScaleTurbulenceSi*RadiusSi_I(1:nX)/cAu
    end select

    ! Handle MHD turbulence => Parallel diffusion coefficient of particles
    dBSi_I(1:nX) = sqrt(cMu*sum(MHData_VIB(Wave1_:Wave2_,1:nX,iLine),1))
    if(UseBtotal) then
       BtotalSi_I(1:nX) = sqrt(BSi_I(1:nX)**2 + dBSi_I(1:nX)**2)
    else
       BtotalSi_I(1:nX) = BSi_I(1:nX)
    end if
    ! Determine what eqn. is used for calculation
    select case(trim(TypeMhdDiffusion))
    case('giacalone1999')
       ! see Giacalone & Jokipii 1999, doi:10.1086/307452
       OmegaCyclo_I(1:nX) = cElectronCharge*BtotalSi_I(1:nX)/cProtonMass
       do iP = 1, nP
          DInnerSi_II(1:nX, iP) = (3.0*SpeedSi_G(iP)**3)/ &
               (40.0*(OmegaCyclo_I(1:nX)/GammaLorentz_G(iP))**2* &
               ScaleSi_I(1:nX)*sin(3.0*cPi/5.0))* (BtotalSi_I/dBSi_I)**2* &
               (1.0 + 72.0/7.0*((OmegaCyclo_I(1:nX)*ScaleSi_I(1:nX))/ &
               (GammaLorentz_G(iP)*SpeedSi_G(iP)))**(5.0/3))/BSi_I(1:nX)
          ! End: Add 1/B as the actual diffusion is D/B
       end do
    case('sokolov2004')
       ! see Sokolov et al., 2004: eq (4),
       ! also see Borovikov 2019, doi:10.48550/arXiv.1911.10165
       ! Gyroradius = cGyroRadius * momentum / |B|
       ! DInner \propto (B/\delta B)**2*Gyroradius*Vel/|B|
       ! ------------------------------------------------------
       ! effective level of turbulence is different for different momenta:
       ! (\delta B)**2 \propto Gyroradius**(1/3)
       CoefLambdaxx_I(1:nX) = (cCoef/3.0)*BSi_I(1:nX)**2 /        &
            (cMu*sum(MHData_VIB(Wave1_:Wave2_,1:nX,iLine),1))*    &
            (ScaleSi_I(1:nX)**2*cGyroRadius*MomentumInjSi/BSi_I(1:nX))&
            **(1.0/3)
       do iP = 1, nP
          DInnerSi_II(1:nX, iP) = CoefLambdaxx_I(1:nX)* &
               Momentum_G(iP)**(1.0/3)*SpeedSi_G(iP)/BSi_I(1:nX)
       end do
    case default
       call CON_Stop(NameSub//": Unknown TypeMhdDiffusion"//TypeMhdDiffusion)
    end select

    ! diffusion is different up- and down-stream
    ! Sokolov et al. 2004, paragraphs before and after eq (4)
    if(UseFixedUps) then
       ! Compute diffusion coefficient without the contribution of v (velocity)
       ! and p (momentum), as v and p are different for different iP
       where(RadiusSi_I(1:nX) > 0.9*RadiusSi_I(iShock))
          ! upstream: reset the diffusion coefficient to
          ! (1/3)*MeanFreePath0InAu[AU]*(R/1AU)*v*(pc/1GeV)**(1/3)
          ! see Li et al. (2003), doi:10.1029/2002JA009666
          ! ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
          ! 1/AU cancels with unit of Lambda0, no special attention needed;
          ! v (velocity) and (p)**(1/3) are calculated in momentum do loop
          CoefLambdaxx_I(1:nX) = (1.0/3.0)*MeanFreePath0InAu* &
               RadiusSi_I(1:nX)*(cLightSpeed*MomentumInjSi/cGeV)**(1.0/3)
       endwhere
       ! update the DInner diffusion coefficient
       do iP = 1, nP
          DInnerSi_II(1:nX, iP) = CoefLambdaxx_I(1:nX)* &
               Momentum_G(iP)**(1.0/3)*SpeedSi_G(iP)/BSi_I(1:nX)
       end do
    end if

    ! for perpendicular diffusion: it becomes "Dpara - Dperp" here
    if(UseDiffusionPerp) then
       DInnerSi_II = DInnerSi_II * (1.0-DPerpRatio)
    end if

    if(nMu > 1) then
       ! for focused transport eq.: CoefLambda_xx => CoefLambda_mumu
       do iP = 1, nP
          CoefLambdaMuMu_II(1:nX, iP) = cCoefxx_to_mumu * &
               DInnerSi_II(1:nX, iP)*BSi_I(1:nX)*3.0/SpeedSi_G(iP)
       end do
    end if

  end subroutine set_diffusion_coef
  !============================================================================
  subroutine tridiag(n, Lower_I, Main_I, Upper_I, Res_I, W_I)
    ! Solve tri-diagonal system of equations:
    !  ||m_1 u_1  0....        || ||w_1|| ||r_1||
    !  ||l_2 m_2 u_2...        || ||w_2|| ||r_2||
    !  || 0  l_3 m_3 u_3       ||.||w_3||=||r_3||
    !  ||...                   || ||...|| ||...||
    !  ||.............0 l_n m_n|| ||w_n|| ||r_n||
    ! From: Numerical Recipes, Chapter 2.6, p.40.

    ! Input parameters
    integer, intent(in) :: n
    real,    intent(in) :: Lower_I(n), Main_I(n), Upper_I(n), Res_I(n)
    ! Output parameters
    real,    intent(out):: W_I(n)
    ! Misc
    integer:: j
    real:: Aux,Aux_I(2:n)
    !--------------------------------------------------------------------------
    if(Main_I(1) == 0.0) call CON_stop('Error in tridiag: Main_I(1)=0')
    Aux = Main_I(1)
    W_I(1) = Res_I(1)/Aux
    do j = 2, n
       Aux_I(j) = Upper_I(j-1)/Aux
       Aux = Main_I(j) - Lower_I(j)*Aux_I(j)
       if(Aux == 0.0) then
          write(*,*) 'Main_I(j), Lower_I(j), Aux_I(j) = ',&
               Main_I(j),Lower_I(j),Aux_I(j)
          write(*,*) ' For j=',j
          call CON_stop('Tridiag failed')
       end if
       W_I(j) = (Res_I(j) - Lower_I(j)*W_I(j-1))/Aux
    end do
    do j = n-1, 1, -1
       W_I(j) = W_I(j) - Aux_I(j+1)*W_I(j+1)
    end do

  end subroutine tridiag
  !============================================================================
end module SP_ModDiffusion
!==============================================================================
