!  Copyright (C) 2002 Regents of the University of Michigan,
!  portions used with permission
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module SP_ModPerpDiffusion

  ! Solve the perpendicular diffusion term in the Parker equation
  ! Dperp acts on the mu-averaged f and its term is added to every mu
  use ModMpi
  use ModNumConst, ONLY: cPi, cTwoPi, cTiny
  use ModUtilities, ONLY: CON_stop
  use SP_ModDistribution, ONLY: Distribution_CB, IsDistNeg, check_dist_neg
  use SP_ModGrid, ONLY: nP, nMu, nLine, nLineAll, State_VIB, R_, X_, Z_
  use SP_ModSize, ONLY: nVertexMax
  use SP_ModTime, ONLY: IsSteadyState
  use SP_ModProc, ONLY: nProc, iProc, iComm, iError
  use SP_ModUnit, ONLY: Io2Si_V, Si2Io_V, UnitX_
  use SP_ModTurbulence, ONLY: DPerp_CB

  implicit none

  PRIVATE ! Except

  SAVE

  ! Public members:
  public :: read_param, init, diffuseperp_distribution, &
       interp_source_intersect, interp_source_linenode
  ! Diffusion in space
  interface diffuseperp_distribution
     ! Global time step Dt
     module procedure diffuseperp_distribution_s
     ! DtLocal for (nP, nVertexMax, nLine)
     module procedure diffuseperp_distribution_c
  end interface diffuseperp_distribution

  ! Whether to include Perpendicular Diffusion
  logical, public :: UseDiffusionPerp = .false.
  character(len=15), public :: TypeDPerp = 'UseDPerpRatio'
  real, public      :: DPerpRatio = 0.01        ! Simple ratio of Dperp/Dpara
  logical, public :: DoTestDPerp = .false.

  ! Grid of the triangulated mesh for solving perp. term: nR * nTheta * nPhi
  integer, public :: nRPerp, nThetaPerp, nPhiPerp
  real            :: dThetaPerp, dPhiPerp       ! dTheta and dPhi
  real, allocatable:: dRPerpMesh_I(:), dRPerpFace_I(:) ! dR of cell center+face
  real            :: RMinPerp = 1.2, RMaxPerp = 240.0  ! RMin and RMax of Grid
  ! Line nodes below RMinGhostPerp (in Rs) are not changed by the perp.
  ! update; nodes from there up to RPerp_C(1) take the shell-1 values
  real, parameter :: RMinGhostPerp = 1.105
  real            :: dLogRFacePerp = 0.0        ! Geometric Sequence for RMesh
  character(len=15):: ScaleRPerp                ! Scale (Linear/Log) along R_
  integer, allocatable :: iRPerpStart_I(:), iRPerpEnd_I(:)
  real, public, allocatable:: RPerp_C(:), ThetaPerp_C(:), PhiPerp_C(:), &
       RPerp_F(:), ThetaPerp_F(:), PhiPerp_F(:), XyzPerp_CB(:,:,:,:) ! SI units

  ! For the FVM solver on the uniform grid: face areas and cell volumes
  real, allocatable:: dVolumeR_I(:), dRPerpFace2_I(:)   ! radial direction
  real, allocatable:: ThetaPerpSin_C(:), dCosTheta_I(:) ! theta direction
  real, allocatable:: &
       AreaR_IIF(:,:,:), AreaTheta_IFI(:,:,:), AreaPhi_FII(:,:,:)    ! faces
  real, public, allocatable :: Volume_CB(:,:,:)                      ! volume
  real, public, allocatable :: Dt_CB(:,:,:)

  ! TEST_PERPDIFF/pos prototype (positivity-preserving line coupling):
  ! a grid node with f at or below DistrTiny gets no relative loss rate
  real, parameter :: DistrTiny = tiny(1.0)
contains
  !============================================================================
  subroutine read_param(NameCommand)
    use ModReadParam, ONLY: read_var
    character(len=*), intent(in) :: NameCommand ! From PARAM.in
    character(len=*), parameter:: NameSub = 'read_param'
    !--------------------------------------------------------------------------
    select case(NameCommand)
    case('#DIFFUSIONPERP')
       call read_var('UseDiffusionPerp', UseDiffusionPerp)
       ! Read in perpendicular diffusion information only when using it
       if(UseDiffusionPerp) then
          call read_var('TypeDPerp', TypeDPerp)
          select case(trim(TypeDPerp))
          case('UseDPerpRatio')
             ! Simply assume DPerp/DPara = DPerpRatio
             call read_var('DPerpRatio', DPerpRatio)
          case default
             call CON_stop(NameSub//": Unknown type of DPerp")
          end select

          ! Handle triangulated mesh in radial direction
          call read_var('nRPerp', nRPerp)
          call read_var('RMinPerp', RMinPerp)
          call read_var('RMaxPerp', RMaxPerp)
          if(nRPerp<=1 .or. RMinPerp<=0.0 .or. RMaxPerp<RMinPerp) RETURN
          if(nRPerp<nProc) then
             write(*,*) "nRPerp = ", nRPerp, "< nProc = ", nProc, &
                  "now set nRPerp = nProc = ", nProc
             ! Ensure every processor owns at least one radial shell
             nRPerp = nProc
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

          ! Print out test information
          call read_var('DoTestDPerp', DoTestDPerp)
       end if
    case default
       call CON_stop('SP:'//NameSub//': unknown command '//NameCommand)
    end select
  end subroutine read_param
  !============================================================================
  subroutine init
    ! setup the uniform spheres and time stepping in perpendicular diffusion

    integer :: iProcChunk, iRPerpRemainder, iPE
    character(len=*), parameter:: NameSub = 'init'
    !--------------------------------------------------------------------------
    ! setup the mesh arrays for triangulation used in perpendicular diffusion
    ! Grid
    if(allocated(RPerp_C)) &
         deallocate(RPerp_C, RPerp_F, &
         dRPerpMesh_I, dRPerpFace_I, XyzPerp_CB)
    allocate(RPerp_C(nRPerp), RPerp_F(nRPerp+1), &
         dRPerpMesh_I(nRPerp-1), dRPerpFace_I(nRPerp), &
         XyzPerp_CB(X_:Z_, nPhiPerp, nThetaPerp, nRPerp))
    XyzPerp_CB = 0.0

    ! Faces and volume used in the FVM solver
    if(allocated(Volume_CB)) &
         deallocate(dVolumeR_I, dRPerpFace2_I, &
         ThetaPerpSin_C, dCosTheta_I, &
         AreaR_IIF, AreaTheta_IFI, AreaPhi_FII, Volume_CB)
    allocate(dVolumeR_I(nRPerp), dRPerpFace2_I(nRPerp), &
         ThetaPerpSin_C(nThetaPerp), dCosTheta_I(nThetaPerp), &
         AreaR_IIF(nPhiPerp, nThetaPerp, nRPerp+1), &
         AreaTheta_IFI(nPhiPerp, nThetaPerp+1, nRPerp), &
         AreaPhi_FII(nPhiPerp+1, nThetaPerp, nRPerp), &
         Volume_CB(nPhiPerp, nThetaPerp, nRPerp))

    ! DPerp coefficient
    if(allocated(DPerp_CB)) deallocate(DPerp_CB)
    allocate(DPerp_CB(nP, nVertexMax, nLine), stat=iError)
    DPerp_CB = 0.0

    ! Determine chunk size, start and end indices for each processor
    if(allocated(iRPerpStart_I)) deallocate(iRPerpStart_I)
    if(allocated(iRPerpEnd_I)) deallocate(iRPerpEnd_I)
    allocate(iRPerpStart_I(0:nProc-1), iRPerpEnd_I(0:nProc-1))
    iProcChunk = nRPerp/nProc             ! Base chunk size
    iRPerpRemainder = mod(nRPerp, nProc)  ! Leftover elements
    ! Loop over processors to compute start and end indices
    do iPE = 0, nProc-1
       if(iPE < iRPerpRemainder) then
          iRPerpStart_I(iPE) = iPE*(iProcChunk+1) + 1
          iRPerpEnd_I(iPE)   = iRPerpStart_I(iPE) + iProcChunk
       else
          iRPerpStart_I(iPE) = iPE*iProcChunk + iRPerpRemainder+1
          iRPerpEnd_I(iPE)   = iRPerpStart_I(iPE) + iProcChunk-1
       end if
    end do

    ! Initialize the uniform grid
    call init_uniform_grid
    ! Precompute the face area and volume of the uniform grid
    call calc_uniform_grid_area_volume

    ! Allocate Dt_CB
    if(IsSteadyState) allocate(Dt_CB(nP, nVertexMax, nLine), stat=iError)

  contains
    !==========================================================================
    subroutine init_uniform_grid
      ! Initialize the uniform grid

      use ModCoordTransform, ONLY: sph_to_xyz
      integer :: iRPerp, iThetaPerp, iPhiPerp ! loop variables
      !------------------------------------------------------------------------
      ! R: radial direction
      select case(trim(ScaleRPerp))
      case("Linear", "linear")
         ! Linear (i.e., Uniform) Grid
         dRPerpMesh_I = (RMaxPerp - RMinPerp)/real(nRPerp) ! Same=Const.
         dRPerpFace_I = dRPerpMesh_I(1) ! Same=Const.
         do iRPerp = 1, nRPerp
            RPerp_F(iRPerp) = RMinPerp + (real(iRPerp)-1.0)*dRPerpMesh_I(1)
         end do
         RPerp_F(nRPerp+1) = RMaxPerp
         RPerp_C = 0.5*(RPerp_F(1:nRPerp) + RPerp_F(2:nRPerp+1))
      case("Exp", "exp", "Exponential", "exponential")
         ! Exponential Grid
         dLogRFacePerp = log(RMaxPerp/RMinPerp)/real(nRPerp)
         do iRPerp = 1, nRPerp+1
            RPerp_F(iRPerp) = RMinPerp * exp((real(iRPerp)-1.0)*dLogRFacePerp)
         end do
         RPerp_C = 0.5*(RPerp_F(1:nRPerp) + RPerp_F(2:nRPerp+1))
         dRPerpFace_I = RPerp_F(2:nRPerp+1) - RPerp_F(1:nRPerp)
         dRPerpMesh_I = RPerp_C(2:nRPerp) - RPerp_C(1:nRPerp-1)
      case("HalfExp", "halfexp", "HalfLinear", "halflinear")
         ! Half Exponential + Half Linear Grid -- From Experience
         dLogRFacePerp = log(RMaxPerp/RMinPerp)/real(nRPerp)
         dRPerpMesh_I = (RMaxPerp - RMinPerp)/real(nRPerp) ! Same=Const.
         do iRPerp = 1, nRPerp+1
            RPerp_F(iRPerp) = &
                 0.5*RMinPerp * exp((real(iRPerp)-1.0)*dLogRFacePerp) + &
                 0.5*(RMinPerp + (real(iRPerp)-1.0)*dRPerpMesh_I(1))
         end do
         RPerp_C = 0.5*(RPerp_F(1:nRPerp) + RPerp_F(2:nRPerp+1))
         dRPerpFace_I = RPerp_F(2:nRPerp+1) - RPerp_F(1:nRPerp)
         dRPerpMesh_I = RPerp_C(2:nRPerp) - RPerp_C(1:nRPerp-1)
      case default
         call CON_stop(NameSub// &
              'Error in ScaleRPerp for Perpendicular Diffusion.')
      end select
      dRPerpMesh_I = dRPerpMesh_I*Io2Si_V(UnitX_)
      dRPerpFace_I = dRPerpFace_I*Io2Si_V(UnitX_)
      RPerp_C = RPerp_C*Io2Si_V(UnitX_)
      RPerp_F = RPerp_F*Io2Si_V(UnitX_)

      ! Theta: zenith/latitudinal direction
      dThetaPerp = cPi/real(nThetaPerp)
      if(allocated(ThetaPerp_C)) deallocate(ThetaPerp_C, ThetaPerp_F)
      allocate(ThetaPerp_C(nThetaPerp), ThetaPerp_F(nThetaPerp+1))
      do iThetaPerp = 1, nThetaPerp
         ThetaPerp_C(iThetaPerp) = (real(iThetaPerp)-0.5)*dThetaPerp
      end do
      ThetaPerp_F(2:nThetaPerp) = 0.5*(ThetaPerp_C(1:nThetaPerp-1) + &
           ThetaPerp_C(2:nThetaPerp))
      ThetaPerp_F(1) = ThetaPerp_C(1) - 0.5*dThetaPerp
      ThetaPerp_F(nThetaPerp+1) = ThetaPerp_C(nThetaPerp) + 0.5*dThetaPerp

      ! Phi: azimuthal/longitudinal direction
      dPhiPerp = cTwoPi/real(nPhiPerp)
      if(allocated(PhiPerp_C)) deallocate(PhiPerp_C, PhiPerp_F)
      allocate(PhiPerp_C(nPhiPerp), PhiPerp_F(nPhiPerp+1))
      do iPhiPerp = 1, nPhiPerp
         PhiPerp_C(iPhiPerp) = (real(iPhiPerp)-0.5)*dPhiPerp
      end do
      PhiPerp_F(2:nPhiPerp) = 0.5*(PhiPerp_C(1:nPhiPerp-1) + &
           PhiPerp_C(2:nPhiPerp))
      PhiPerp_F(1) = PhiPerp_C(1) - 0.5*dPhiPerp
      PhiPerp_F(nPhiPerp+1) = PhiPerp_C(nPhiPerp) + 0.5*dPhiPerp

      ! (R, Theta, Phi) => (X, Y, Z) + Parallelization
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

    end subroutine init_uniform_grid
    !==========================================================================
    subroutine calc_uniform_grid_area_volume
      ! Pre-calculate the face areas and cell volumes

      integer :: i, j
      !------------------------------------------------------------------------
      ! Grids
      ! Radial: RPerp_C, RPerp_F, dRPerpMesh_I, dRPerpFace_I, dRPerpFace2_I
      dRPerpFace2_I = 0.5*(RPerp_F(2:nRPerp+1)**2 - RPerp_F(1:nRPerp)**2)
      ! Theta: ThetaPerp_C, ThetaPerp_F, dThetaPerp
      ThetaPerpSin_C = sin(ThetaPerp_C) ! needed also in FVM internal fluxes
      ! Phi: PhiPerp_C, PhiPerp_F, dPhiPerp -- Ready

      ! Precompute: Face area and cell volume
      dVolumeR_I = (RPerp_F(2:nRPerp+1)**3 - RPerp_F(1:nRPerp)**3)/3.0
      dCosTheta_I = cos(ThetaPerp_F(1:nThetaPerp)) - &
           cos(ThetaPerp_F(2:nThetaPerp+1))
      do i = 1, nRPerp+1
         do j = 1, nThetaPerp
            AreaR_IIF(:,j,i) = RPerp_F(i)**2*dCosTheta_I(j)*dPhiPerp
         end do
      end do
      do i = 1, nRPerp
         do j = 1, nThetaPerp
            Volume_CB(:,j,i) = dVolumeR_I(i)*dCosTheta_I(j)*dPhiPerp
            AreaTheta_IFI(:,j,i) = dRPerpFace2_I(i)* &
                 sin(ThetaPerp_F(j))*dPhiPerp
         end do
         AreaTheta_IFI(:,nThetaPerp+1,i) = dRPerpFace2_I(i)* &
              sin(ThetaPerp_F(nThetaPerp+1))*dPhiPerp
         AreaPhi_FII(:,:,i) = dRPerpFace2_I(i)*dThetaPerp
      end do

    end subroutine calc_uniform_grid_area_volume
    !==========================================================================
  end subroutine init
  !============================================================================
  subroutine diffuseperp_distribution_s(IsFirstCall, dtIn)
    logical, intent(in) :: IsFirstCall
    real, intent(in)    :: dtIn
    real :: dtFake_C(nP, nVertexMax, nLine)
    !--------------------------------------------------------------------------
    dtFake_C = dtIn
    call diffuseperp_distribution_c(IsFirstCall, dtFake_C)
  end subroutine diffuseperp_distribution_s
  !============================================================================
  subroutine diffuseperp_distribution_c(IsFirstCall, dtIn_CB)
    ! cross-field (perpendicular) diffusion related steps and manipulations

    use SP_ModTriangulate, ONLY: reset_intersect_surf, intersect_surf, &
         build_trmesh, interpolate_trmesh, nTriMesh_I, XyzOrig_DII

    ! inputs
    logical, intent(in) :: IsFirstCall
    real, intent(in)    :: dtIn_CB(nP, nVertexMax, nLine)
    ! loop variables
    integer :: iRPerp, iThetaPerp, iPhiPerp, iPE
    ! mu-averaged VDF and the source term (4D)
    !--------------------------------------------------------------------------
    real, dimension(0:nP+1, nPhiPerp, nThetaPerp, nRPerp)&
         :: DistrPerp_4D, source_4D
    ! DPerp coefficients (4D)
    real :: DPerp_4D(nP, nPhiPerp, nThetaPerp, nRPerp)
    ! Positivity-preserving coupling: grid f before the sub-steps, and the
    ! split of the averaged df/dt (source_4D) into an absolute gain rate
    ! Gain >= 0 and a relative loss rate LossRel <= 0 (per unit f), on the
    ! uniform grid (4D) and at the line-shell intersections (IV)
    real, dimension(nP, nPhiPerp, nThetaPerp, nRPerp) :: &
         DistrOld_4D, Gain_4D, LossRel_4D
    real, dimension(nP, nLineAll, nRPerp) :: Gain_IV, LossRel_IV
    ! VDF interpolated to a grid cell for every mu (before the mu-average)
    real :: DistrMu_II(0:nP+1, nMu)
    character(len=*), parameter:: NameSub = 'diffuseperp_distribution_c'
    !--------------------------------------------------------------------------
    if(IsFirstCall) then
       ! In the 1st call of DPerp: we set up the skeleton, coefficient, and VDF
       DistrPerp_4D = 0.0; DPerp_4D = 0.0

       ! Step 1: field lines => intersection points on multiple uniform layers
       ! here, iProc is for field lines, not for sub-slices/layers.
       ! In fact, this step of getting the intersection points is needed ONLY
       ! the 1st time, since MHData would be changed every coupling interval.
       call reset_intersect_surf(nRPerp)
       do iRPerp = 1, nRPerp
          call intersect_surf(RPerp_C(iRPerp)*Si2Io_V(UnitX_), &
               0, iRPerp, IsIncludeDPerpIn=.true.)
       end do
       if(nProc > 1) then
          do iPE = 0, nProc-1
             call MPI_BCAST(nTriMesh_I, nRPerp, &
                  MPI_INTEGER, iPE, iComm, iError)
          end do
          ! Intersection points of all lines are complete only on the
          ! root (0) of intersect_surf; interp_source_intersect needs
          ! them on every processor
          call MPI_BCAST(XyzOrig_DII, size(XyzOrig_DII), &
               MPI_REAL, 0, iComm, iError)
       end if

       if(iProc == 0) then
          ! Step 2: Get DPerp at the cell center (uniform grid; ONLY 1st time)
          ! Now, iProc is for sub-slices/layers in the `DistrPerp_4D` array
          do iRPerp = 1, nRPerp
             ! intersection points => construct the triangulation skeleton
             ! (also fix DPerp at the polar nodes before interpolating it)
             call build_trmesh(iRPerp, IsIncludeDPerpIn=.true.)
             do iThetaPerp = 1, nThetaPerp
                do iPhiPerp = 1, nPhiPerp
                   ! Get DPerp coefficient at the cell center (uniform grid)
                   call interpolate_trmesh(XyzPerp_CB(:,             &
                        iPhiPerp,iThetaPerp,iRPerp)*Si2Io_V(UnitX_), &
                        iRIn = iRPerp, DPerp_I =                     &
                        DPerp_4D(:,iPhiPerp,iThetaPerp,iRPerp))
                end do
             end do
          end do
       end if
       ! Broadcast DPerp coefficient to all processes
       if(nProc > 1) then
          do iPE = 0, nProc-1
             call MPI_BCAST(DPerp_4D, nPhiPerp*nThetaPerp*nRPerp* &
                  nP, MPI_REAL, iPE, iComm, iError)
          end do
       end if
       ! Check the error message from after interpolate_trmesh
       if(iError /= 0) then
          write(*,*) 'iProc = ', iProc, NameSub//&
               ': interpolate_trmesh for DPerp_4D failed, DPerp stopped'
          RETURN
       end if

       ! Output test information
       if(DoTestDPerp) write(*,*) NameSub//'Finish FirstCall: '// &
            'Triangulation skeleton and DPerp Coef in the uniform grid'
    end if

    ! Step 3: Triangulation skeleton => Interpolate the VDF
    if(iProc == 0) then
       do iRPerp = 1, nRPerp
          call build_trmesh(iRPerp)
          do iThetaPerp = 1, nThetaPerp
             do iPhiPerp = 1, nPhiPerp
                ! Get VDF at the cell center in the uniform grid
                call interpolate_trmesh(XyzPerp_CB(:,             &
                     iPhiPerp,iThetaPerp,iRPerp)*Si2Io_V(UnitX_), &
                     iRIn = iRPerp, DistrInterp_II = DistrMu_II)
                ! mu-average of f: the mu grid is uniform (DeltaMu), so
                ! all weights are equal; f itself for nMu = 1 (bitwise)
                if(nMu == 1) then
                   DistrPerp_4D(:,iPhiPerp,iThetaPerp,iRPerp) = &
                        DistrMu_II(:,1)
                else
                   DistrPerp_4D(:,iPhiPerp,iThetaPerp,iRPerp) = &
                        sum(DistrMu_II, DIM=2)/real(nMu)
                end if
             end do
          end do
       end do
    end if
    ! Broadcast Dperp coefficient to all processes
    if(nProc > 1) then
       do iPE = 0, nProc-1
          call MPI_BCAST(DistrPerp_4D, nPhiPerp*nThetaPerp*nRPerp* &
               (nP+2), MPI_REAL, iPE, iComm, iError)
       end do
    end if

    ! Save the grid f before the sub-steps (the same on all processors)
    DistrOld_4D = DistrPerp_4D(1:nP,:,:,:)

    ! Step 4: Solve the DPerp equation in multiple uniform layers by FVM
    source_4D = 0.0
    ! Pass the shells of this processor only, so that the dummy shell
    ! range iRStart:iREnd refers to the same global shells
    call solver_fvm_diffuse3d(DPerp_4D,                 &
         iRPerpStart_I(iProc), iRPerpEnd_I(iProc),      &
         IsSteadyState, dtIn_CB(1,1,1),                 &
         DistrPerp_4D(:,:,:,iRPerpStart_I(iProc):iRPerpEnd_I(iProc)), &
         source_4D(:,:,:,iRPerpStart_I(iProc):iRPerpEnd_I(iProc)))
    ! Broadcast source_4D (uniform grid) to all processes
    if(nProc > 1) then
       do iPE = 0, nProc-1
          call MPI_BCAST(source_4D(:,:,:,               &
               iRPerpStart_I(iPE):iRPerpEnd_I(iPE)),    &
               (iRPerpEnd_I(iPE)-iRPerpStart_I(iPE)+1)* &
               (nP+2)*nPhiPerp*nThetaPerp, MPI_REAL, iPE, iComm, iError)
       end do
    end if
    ! Check the error message from after interpolate_trmesh
    if(iError /= 0) then
       write(*,*) 'iProc = ', iProc, NameSub//&
            ': interpolate_trmesh for source_4D failed, DPerp stopped'
       RETURN
    end if

    ! Step 4b: split df/dt at every grid node into an absolute gain rate
    ! and a relative loss rate. Every processor holds the full source_4D
    ! and DistrOld_4D, so this element-wise split is identical on all.
    Gain_4D = max(source_4D(1:nP,:,:,:), 0.0)
    where(DistrOld_4D > DistrTiny)
       LossRel_4D = min(source_4D(1:nP,:,:,:), 0.0)/DistrOld_4D
    elsewhere
       LossRel_4D = 0.0
    end where
    ! Time-accurate: a node cannot lose more than it had, dt*LossRel >= -1
    ! (interp_source_linenode also clips 1 + dt*LossRel at 0 on every node)
    if(.not.IsSteadyState .and. dtIn_CB(1,1,1) > 0.0) &
         LossRel_4D = max(LossRel_4D, -1.0/dtIn_CB(1,1,1))

    ! Step 5: Interpolate gain and relative loss back to the line nodes
    ! uniform grid => intersection points
    ! (pass the shells of this processor, matching the dummy shell range)
    call interp_source_intersect(                                     &
         Gain_4D(:,:,:,iRPerpStart_I(iProc):iRPerpEnd_I(iProc)),      &
         LossRel_4D(:,:,:,iRPerpStart_I(iProc):iRPerpEnd_I(iProc)),   &
         Gain_IV, LossRel_IV)
    ! intersection points => field line nodes
    call interp_source_linenode(Gain_IV, LossRel_IV, dtIn_CB)

  end subroutine diffuseperp_distribution_c
  !============================================================================
  subroutine solver_fvm_diffuse3d(DPerp_4D, iRStart, iREnd, &
       IsSteadyStateIn, dtIn, DistrPerp_4D, source_4D)
    ! FVM solver for pure 3d diffusion equation: df/dt = div (DPerp grad_f)

    integer, intent(in) :: iRStart, iREnd
    logical, intent(in) :: IsSteadyStateIn
    real,    intent(in) :: dtIn
    real,    intent(in) :: DPerp_4D(1:nP,nPhiPerp,nThetaPerp,nRPerp)
    real, intent(inout) :: DistrPerp_4D(0:nP+1, &
         nPhiPerp, nThetaPerp, iRStart:iREnd)
    real, intent(out):: source_4D(0:nP+1,nPhiPerp,nThetaPerp,iRStart:iREnd)

    ! Proc index used for exchanging ghost VDF
    integer :: iProcPrev, iProcNext
    ! loop variables
    integer :: i, j, k, iStep
    ! time step and iterations
    real :: dt
    real :: dtLocal_III(nPhiPerp,nThetaPerp,nRPerp)
    integer :: nStep
    ! smallest radial spacing in the radial stencil of a cell
    real :: dRMin
    ! largest Dperp in a cell; largest sum of the FVM face coefficients of
    ! a cell (the update keeps f >= 0 if dt*LambdaMax <= 1); smallest grid
    ! f over all sub-steps (diagnostic)
    real :: DPerpMax, LambdaMax, DistrMinGrid

    ! ghost VDF when sendrecv across prev-current-next processors
    real :: DistrPrevProc_G(0:nP+1, nPhiPerp, nThetaPerp), &
         DistrNextProc_G(0:nP+1, nPhiPerp, nThetaPerp)
    real :: sourceIncrement_4D(nP,nPhiPerp,nThetaPerp,iRStart:iREnd)
    character(len=*), parameter:: NameSub = 'solver_fvm_diffuse3d'
    !--------------------------------------------------------------------------
    iProcPrev = iProc - 1
    iProcNext = iProc + 1
    source_4D = 0.0

    if(IsSteadyStateIn) then
       ! Steady-state mode: Only 1 step, no use of dt indeed
       nStep = 1
       dt = 0.0
    else
       ! Time-accurate mode: Use exact dt; otherwise, no time marching
       ! Stability check: Determine maximal allowed global time step:
       ! dt <= 0.5*min_{i,j,k} [D_{ijk} * (
       !       1/DeltaR_i**2 + 1/(R_i*Delta(theta_j))**2 +
       !       1/(R_i*sin(theta_j)*Delta(phi_k))**2 )]^{-1}
       dtLocal_III = 0.0
       do i = 1, nRPerp
          ! DeltaR_i: smallest of the cell width and the centre-to-centre
          ! spacings of the inner radial faces of cell i (as used in
          ! calc_fvm_netflux); it equals dRPerpFace_I(i) on Linear grid
          dRMin = dRPerpFace_I(i)
          if(i > 1) dRMin = min(dRMin, dRPerpMesh_I(i-1))
          if(i < nRPerp) dRMin = min(dRMin, dRPerpMesh_I(i))
          do j = 1, nThetaPerp
             do k = 1, nPhiPerp
                ! A cell with no Dperp > 0 does not limit dt (no division
                ! by Dperp = 0, e.g. for DPerpRatio = 0)
                DPerpMax = maxval(DPerp_4D(:,k,j,i))
                if(DPerpMax > 0.0) dtLocal_III(k,j,i) = 0.5/(1.0/dRMin**2 + &
                     1.0/(RPerp_C(i)*dThetaPerp)**2 + &
                     1.0/(RPerp_C(i)*ThetaPerpSin_C(j)*dPhiPerp)**2)/DPerpMax
             end do
          end do
       end do
       if(any(dtLocal_III > 0.0) .and. dtIn > 0.0) then
          dt = minval(dtLocal_III, mask=(dtLocal_III>0))
          nStep = int(dtIn/dt) + 1
       else
          ! Dperp = 0 everywhere, or dtIn <= 0: one sub-step, which does
          ! not change DistrPerp_4D (never nStep = 0, so no 0/0 below)
          nStep = 1
       end if
       dt = dtIn/real(nStep)

       ! Positivity check of the explicit update (all processors compute
       ! the same LambdaMax from the full DPerp_4D)
       call calc_lambda_max
       if(DoTestDPerp .and. iProc == 0) write(*,*) 'DPerp pos: LambdaMax=', &
            LambdaMax, ' dt*LambdaMax (stability dt)=', dt*LambdaMax
       ! The stability dt above can exceed 1/LambdaMax slightly on the
       ! spherical grid; if so, use the fewest extra sub-steps that give
       ! dt*LambdaMax < 1 (no change if the bound already holds)
       if(dtIn*LambdaMax >= real(nStep)) then
          nStep = int(dtIn*LambdaMax) + 1
          dt = dtIn/real(nStep)
       end if

       ! Output test information
       if(DoTestDPerp) &
            write(*,*) 'DPerp: Stable (Used) dt=', dt, 'with nStep=', nStep
    end if
    ! Smallest grid f before and during the sub-steps (diagnostic)
    DistrMinGrid = minval(DistrPerp_4D(1:nP,:,:,:))

    ! Time stepping: Diffusion in 3D by FVM
    DIFFUSE3D:do iStep = 1, nStep
       ! Exchange ghostVDF: iProcPrev |<=| (Start) iProc (End) |=>| iProcNext
       call exchange_ghostVDF
       ! Compute internal fluxes of the FVM solver = source term of df/dt
       call calc_fvm_netflux
       ! Update source_4D (= df/dt) and DistrPerp_4D
       source_4D(1:nP,:,:,:) = source_4D(1:nP,:,:,:) + sourceIncrement_4D
       DistrPerp_4D(1:nP,:,:,:) = DistrPerp_4D(1:nP,:,:,:) + &
            sourceIncrement_4D*dt   ! f^{n+1} = f^{n} + (df/dt)*delta_t
       DistrMinGrid = min(DistrMinGrid, minval(DistrPerp_4D(1:nP,:,:,:)))
    end do DIFFUSE3D
    ! Average df/dt over the nStep sub-steps, so that source_4D*dtIn
    ! equals the change of DistrPerp_4D over the full time step dtIn
    source_4D(1:nP,:,:,:) = source_4D(1:nP,:,:,:)/real(nStep)

    ! Report the smallest grid f over all shells and sub-steps
    if(DoTestDPerp) then
       if(nProc > 1) call MPI_ALLREDUCE(MPI_IN_PLACE, DistrMinGrid, 1, &
            MPI_REAL, MPI_MIN, iComm, iError)
       if(iProc == 0) write(*,*) &
            'DPerp pos: grid min f over all sub-steps=', DistrMinGrid
    end if

  contains
    !==========================================================================
    subroutine calc_lambda_max
      ! LambdaMax = max over cells and p of the sum, over the faces
      ! with a nonzero gradient, of DPerp_face*Area_face/(dX_face*Volume),
      ! exactly as used in calc_fvm_netflux. The update
      ! f_c + dt*sum_f c_f*(f_nb - f_c) = (1 - dt*Lambda_c)*f_c + dt*sum ...
      ! keeps f >= 0 if dt*LambdaMax <= 1 (all c_f >= 0 since Dperp >= 0).
      ! Zero-gradient faces (radial ends, theta=0 and pi) are skipped.

      integer :: kP, kM
      real    :: Lambda_I(nP)
      !------------------------------------------------------------------------
      LambdaMax = 0.0
      do i = 1, nRPerp
         do j = 1, nThetaPerp
            do k = 1, nPhiPerp
               Lambda_I = 0.0
               ! Radial faces (i, i+1) and (i-1, i)
               if(i < nRPerp) Lambda_I = Lambda_I + &
                    0.5*(DPerp_4D(:,k,j,i) + DPerp_4D(:,k,j,i+1)) &
                    /dRPerpMesh_I(i)*AreaR_IIF(k,j,i+1)
               if(i > 1) Lambda_I = Lambda_I + &
                    0.5*(DPerp_4D(:,k,j,i) + DPerp_4D(:,k,j,i-1)) &
                    /dRPerpMesh_I(i-1)*AreaR_IIF(k,j,i)
               ! Theta faces (j, j+1) and (j-1, j)
               if(j < nThetaPerp) Lambda_I = Lambda_I + &
                    0.5*(DPerp_4D(:,k,j,i) + DPerp_4D(:,k,j+1,i)) &
                    /(RPerp_C(i)*dThetaPerp)*AreaTheta_IFI(k,j+1,i)
               if(j > 1) Lambda_I = Lambda_I + &
                    0.5*(DPerp_4D(:,k,j-1,i) + DPerp_4D(:,k,j,i)) &
                    /(RPerp_C(i)*dThetaPerp)*AreaTheta_IFI(k,j,i)
               ! Phi faces (periodic)
               kP = mod(k,nPhiPerp)+1
               kM = mod(k-2+nPhiPerp,nPhiPerp)+1
               Lambda_I = Lambda_I + &
                    (0.5*(DPerp_4D(:,k,j,i) + DPerp_4D(:,kP,j,i)) &
                    *AreaPhi_FII(kP,j,i) &
                    + 0.5*(DPerp_4D(:,kM,j,i) + DPerp_4D(:,k,j,i)) &
                    *AreaPhi_FII(kM,j,i)) &
                    /(RPerp_C(i)*ThetaPerpSin_C(j)*dPhiPerp)
               LambdaMax = max(LambdaMax, maxval(Lambda_I)/Volume_CB(k,j,i))
            end do
         end do
      end do
    end subroutine calc_lambda_max
    !==========================================================================
    subroutine exchange_ghostVDF
      ! Exchange ghostVDF: iProcPrev |<=| (Start) iProc (End) |=>| iProcNext
      !------------------------------------------------------------------------
      ! Exchange VDF at ghost cells 1: iProc End |=>| iProcNext Start
      if(iProcNext < nProc) then
         call MPI_SENDRECV(DistrPerp_4D(:,:,:,iREnd), &
              (nP+2)*nPhiPerp*nThetaPerp, &
              MPI_REAL, iProcNext, 0, &
              DistrNextProc_G, (nP+2)*nPhiPerp*nThetaPerp, &
              MPI_REAL, iProcNext, 1, &
              iComm, MPI_STATUS_IGNORE, iError)
      else
         ! zero-gradient (Neumann BC) at outer boundary
         DistrNextProc_G = DistrPerp_4D(:,:,:,iREnd)
      end if

      ! Exchange VDF at ghost cells 2: iProcPrev End |<=| iProc Start
      if(iProcPrev >= 0) then
         call MPI_SENDRECV(DistrPerp_4D(:,:,:,iRStart), &
              (nP+2)*nPhiPerp*nThetaPerp, &
              MPI_REAL, iProcPrev, 1, &
              DistrPrevProc_G, (nP+2)*nPhiPerp*nThetaPerp, &
              MPI_REAL, iProcPrev, 0, &
              iComm, MPI_STATUS_IGNORE, iError)
      else
         ! zero-gradient (Neumann BC) at inner boundary
         DistrPrevProc_G = DistrPerp_4D(:,:,:,iRStart)
      end if

    end subroutine exchange_ghostVDF
    !==========================================================================
    subroutine calc_fvm_netflux
      ! calculate internal net fluxes in the FVM solver

      integer :: jPlus, jMinus, kPlus, kMinus
      real    :: VolumeInv
      ! radial distance used in gradF at face(i, ip) and face(im, i)
      real    :: dRPlus, dRMinus
      real, dimension(1:nP) :: DPerpPlus_I, DPerpMinus_I
      real, dimension(0:nP+1) :: DistrPlus_I, DistrMinus_I
      !------------------------------------------------------------------------
      sourceIncrement_4D = 0.0
      do i = iRStart, iREnd
         do j = 1, nThetaPerp
            ! Ready for Theta fluxes
            jPlus = min(j+1, nThetaPerp) ! j or nThetaPerp
            jMinus = max(j-1, 1) ! j or 1

            do k = 1, nPhiPerp
               VolumeInv = 1.0/Volume_CB(k,j,i)

               ! Part 1 -- Face_r: Radial fluxes
               ! VDF at i+1 (ip)
               if(i < iREnd) then
                  DistrPlus_I = DistrPerp_4D(:,k,j,i+1)
               else
                  DistrPlus_I = DistrNextProc_G(:,k,j)
               end if
               ! VDF at i-1 (im)
               if(i > iRStart) then
                  DistrMinus_I = DistrPerp_4D(:,k,j,i-1)
               else
                  DistrMinus_I = DistrPrevProc_G(:,k,j)
               end if
               ! DPerp coefficient at face(i, i+1)
               ! gradF at an inner face uses the centre-to-centre spacing
               ! (conservative flux); at the outer boundary face the
               ! gradient is zero (Neumann BC), so keep the cell width
               if(i < nRPerp) then
                  DPerpPlus_I = 0.5*(Dperp_4D(:,k,j,i) &
                       + Dperp_4D(:,k,j,i+1))
                  dRPlus = dRPerpMesh_I(i)
               else
                  DPerpPlus_I = Dperp_4D(:,k,j,i) + &
                       0.5*(Dperp_4D(:,k,j,i) - Dperp_4D(:,k,j,i-1))
                  dRPlus = dRPerpFace_I(i)
               end if
               ! DPerp coefficient at face(i-1, i)
               ! (same spacing choice; zero gradient at the inner boundary)
               if(i > 1) then
                  DPerpMinus_I = 0.5*(Dperp_4D(:,k,j,i) &
                       + Dperp_4D(:,k,j,i-1))
                  dRMinus = dRPerpMesh_I(i-1)
               else
                  DPerpMinus_I = Dperp_4D(:,k,j,i) - &
                       0.5*(Dperp_4D(:,k,j,i+1) - Dperp_4D(:,k,j,i))
                  dRMinus = dRPerpFace_I(i)
               end if
               ! FVM net flux: sourceIncrement_4D += (flux_ip - flux_im)
               sourceIncrement_4D(:,k,j,i) = sourceIncrement_4D(:,k,j,i)+ &
                    DPerpPlus_I &         ! + DPerp at face(i, ip)
                    *(DistrPlus_I(1:nP) - DistrPerp_4D(1:nP,k,j,i))&
                    /dRPlus &             ! * gradF at face(i, ip)
                    *AreaR_IIF(k,j,i+1) & ! * area of face(i, ip)
                    *VolumeInv            ! / volume of (k, j, i)
               sourceIncrement_4D(:,k,j,i) = sourceIncrement_4D(:,k,j,i)- &
                    DPerpMinus_I &        ! - DPerp at face(im, i)
                    *(DistrPerp_4D(1:nP,k,j,i) - DistrMinus_I(1:nP))&
                    /dRMinus &            ! * gradF at face(im, i)
                    *AreaR_IIF(k,j,i) &   ! * area of face(im, i)
                    *VolumeInv            ! / volume of (k, j, i)

               ! Part 2 -- Face_theta: Theta fluxes
               ! DPerp coefficient at face(j, j+1)
               if(j < nThetaPerp) then
                  DPerpPlus_I = 0.5*(Dperp_4D(:,k,j,i) &
                       + Dperp_4D(:,k,j+1,i))
               else
                  DPerpPlus_I = Dperp_4D(:,k,j,i) + &
                       0.5*(Dperp_4D(:,k,j,i) - Dperp_4D(:,k,j-1,i))
               end if
               ! DPerp coefficient at face(j-1, j)
               if(j > 1) then
                  DPerpMinus_I = 0.5*(Dperp_4D(:,k,j-1,i) &
                       + Dperp_4D(:,k,j,i))
               else
                  DPerpMinus_I = Dperp_4D(:,k,j,i) - &
                       0.5*(Dperp_4D(:,k,j+1,i) - Dperp_4D(:,k,j,i))
               end if
               ! FVM net flux: sourceIncrement_4D += (flux_jp - flux_jm)
               sourceIncrement_4D(:,k,j,i) = sourceIncrement_4D(:,k,j,i)+ &
                    DPerpPlus_I &               ! + DPerp at face(j, jp)
                    *(DistrPerp_4D(1:nP,k,jPlus,i)&
                    -DistrPerp_4D(1:nP,k,j,i))&
                    /(RPerp_C(i)*dThetaPerp) &  ! * gradF at face(j, jp)
                    *AreaTheta_IFI(k,j+1,i) &   ! * area of face(j, jp)
                    *VolumeInv                  ! / volume of (k, j, i)
               sourceIncrement_4D(:,k,j,i) = sourceIncrement_4D(:,k,j,i)- &
                    DPerpMinus_I &              ! - DPerp at face(jm, j)
                    *(DistrPerp_4D(1:nP,k,j,i)&
                    -DistrPerp_4D(1:nP,k,jMinus,i))&
                    /(RPerp_C(i)*dThetaPerp) &  ! * gradF at face(jm, j)
                    *AreaTheta_IFI(k,j,i) &     ! * area of face(jm, j)
                    *VolumeInv                  ! / volume of (k, j, i)

               ! Part 3 -- Face_phi: Phi fluxes
               kPlus = mod(k,nPhiPerp)+1              ! next index with wrap
               kMinus = mod(k-2+nPhiPerp,nPhiPerp)+1  ! prev index with wrap
               ! DPerp coefficient at face(k, k+1)
               DPerpPlus_I = 0.5*(Dperp_4D(:,k,j,i) &
                    + Dperp_4D(:,kPlus,j,i))
               ! DPerp coefficient at face(k-1, k)
               DPerpMinus_I = 0.5*(Dperp_4D(:,kMinus,j,i) &
                    + Dperp_4D(:,k,j,i))
               ! FVM net flux: sourceIncrement_4D += (flux_kp - flux_km)
               sourceIncrement_4D(:,k,j,i) = sourceIncrement_4D(:,k,j,i)+ &
                    DPerpPlus_I &               ! + DPerp at face(k, kp)
                    *(DistrPerp_4D(1:nP,kPlus,j,i)& ! * gradF at face(k, kp)
                    -DistrPerp_4D(1:nP,k,j,i))&
                    /(RPerp_C(i)*ThetaPerpSin_C(j)*dPhiPerp) &
                    *AreaPhi_FII(kPlus,j,i) &   ! * area of face(k, kp)
                    *VolumeInv                  ! / volume of (k, j, i)
               sourceIncrement_4D(:,k,j,i) = sourceIncrement_4D(:,k,j,i)- &
                    DPerpMinus_I &              ! - DPerp at face(km, k)
                    *(DistrPerp_4D(1:nP,k,j,i)& ! * gradF at face(km, k)
                    -DistrPerp_4D(1:nP,kMinus,j,i))&
                    /(RPerp_C(i)*ThetaPerpSin_C(j)*dPhiPerp) &
                    *AreaPhi_FII(kMinus,j,i) &  ! * area of face(km, k)
                    *VolumeInv                  ! / volume of (k, j, i)
            end do
         end do
      end do

    end subroutine calc_fvm_netflux
    !==========================================================================
  end subroutine solver_fvm_diffuse3d
  !============================================================================
  subroutine interp_source_intersect(Gain_4D, LossRel_4D, Gain_IV, LossRel_IV)
    ! Interpolate the gain and relative loss rates from the uniform grid to
    ! all intersection points of the lines and multiple slices, bilinearly
    ! in (lon, colat) with convex weights: each in [0,1], summing to 1

    use ModCoordTransform, ONLY: xyz_to_rlonlat
    use SP_ModTriangulate, ONLY: XyzOrig_DII
    real, intent(in), dimension(nP, nPhiPerp, nThetaPerp, &
         iRPerpStart_I(iProc):iRPerpEnd_I(iProc)) :: Gain_4D, LossRel_4D
    real, intent(out), dimension(nP, nLineAll, nRPerp) :: &
         Gain_IV, LossRel_IV
    integer :: iRPerp, iLine, iPE, nSize
    integer :: iLon1, iLon2, iColat1, iColat2
    real    :: XyzPoint_D(X_:Z_), rPoint, lonPoint, latPoint, colatPoint
    real    :: lonFrac, colatFrac, Weight_I(4)
    character(len=*), parameter:: NameSub = 'interp_source_intersect'
    !--------------------------------------------------------------------------
    do iRPerp = iRPerpStart_I(iProc), iRPerpEnd_I(iProc)
       do iLine = 1, nLineAll
          ! Extract the intersection point coordinates
          XyzPoint_D = XyzOrig_DII(:, iLine, iRPerp)
          if(sum(XyzPoint_D**2) > cTiny) then
             ! CoordTransform: XyzPoint_D => (rPoint, lonPoint, latPoint)
             call xyz_to_rlonlat(XyzPoint_D, rPoint, lonPoint, latPoint)
             colatPoint = 0.5*cPi - latPoint

             ! Longitude (periodic; cell centres PhiPerp_C): measure from
             ! PhiPerp_C(1) with modulo, so that a point at
             ! lon < PhiPerp_C(1) lies between nPhiPerp and 1 (lon seam)
             lonFrac = modulo(lonPoint - PhiPerp_C(1), cTwoPi)/dPhiPerp
             iLon1 = min(max(floor(lonFrac), 0), nPhiPerp-1)
             lonFrac = min(max(lonFrac - real(iLon1), 0.0), 1.0)
             iLon1 = iLon1 + 1
             iLon2 = mod(iLon1, nPhiPerp) + 1

             ! Colatitude (not periodic): in the polar caps outside
             ! ThetaPerp_C(1) ... ThetaPerp_C(nThetaPerp) the fraction is
             ! clamped to [0,1], i.e. the nearest ring is used
             colatFrac = (colatPoint - ThetaPerp_C(1))/dThetaPerp
             iColat1 = min(max(floor(colatFrac), 0), nThetaPerp-2)
             colatFrac = min(max(colatFrac - real(iColat1), 0.0), 1.0)
             iColat1 = iColat1 + 1
             iColat2 = iColat1 + 1

             ! Convex bilinear weights on the uniform LON-COLAT grid
             Weight_I = [(1.0-lonFrac)*(1.0-colatFrac), &
                  (1.0-lonFrac)*colatFrac, lonFrac*(1.0-colatFrac), &
                  lonFrac*colatFrac]
             Gain_IV(:,iLine,iRPerp) = &
                  Weight_I(1)*Gain_4D(:,iLon1,iColat1,iRPerp) + &
                  Weight_I(2)*Gain_4D(:,iLon1,iColat2,iRPerp) + &
                  Weight_I(3)*Gain_4D(:,iLon2,iColat1,iRPerp) + &
                  Weight_I(4)*Gain_4D(:,iLon2,iColat2,iRPerp)
             LossRel_IV(:,iLine,iRPerp) = &
                  Weight_I(1)*LossRel_4D(:,iLon1,iColat1,iRPerp) + &
                  Weight_I(2)*LossRel_4D(:,iLon1,iColat2,iRPerp) + &
                  Weight_I(3)*LossRel_4D(:,iLon2,iColat1,iRPerp) + &
                  Weight_I(4)*LossRel_4D(:,iLon2,iColat2,iRPerp)
          else
             Gain_IV(:,iLine,iRPerp) = 0.0
             LossRel_IV(:,iLine,iRPerp) = 0.0
          end if
       end do
    end do

    ! Broadcast Gain_IV and LossRel_IV (line intersections) to all processes
    if(nProc > 1) then
       do iPE = 0, nProc-1
          nSize = (iRPerpEnd_I(iPE)-iRPerpStart_I(iPE)+1)*nP*nLineAll
          call MPI_BCAST(Gain_IV(:,:,                   &
               iRPerpStart_I(iPE):iRPerpEnd_I(iPE)),    &
               nSize, MPI_REAL, iPE, iComm, iError)
          call MPI_BCAST(LossRel_IV(:,:,                &
               iRPerpStart_I(iPE):iRPerpEnd_I(iPE)),    &
               nSize, MPI_REAL, iPE, iComm, iError)
       end do
    end if
    ! Check the error message from after interpolate_trmesh
    if(iError /= 0) then
       write(*,*) 'iProc = ', iProc, NameSub//&
            ': interpolate_trmesh for source_IV failed, DPerp stopped'
       RETURN
    end if
  end subroutine interp_source_intersect
  !============================================================================
  subroutine interp_source_linenode(Gain_IV, LossRel_IV, dtIn_CB)
    ! Interpolate the gain and relative loss rates from the intersection
    ! points of the lines and multiple slices to the original points along
    ! the field line (convex radial weights), and update the VDF there as
    !    f = f*max(0, 1 + dt*LossRel) + dt*Gain
    ! which keeps f >= 0 for any Dperp, and leaves f unchanged for Dperp=0

    use SP_ModGrid, ONLY: Used_B, nVertex_B, iLineAll0
    real, intent(in), dimension(nP, nLineAll, nRPerp) :: &
         Gain_IV, LossRel_IV
    real, intent(in):: dtIn_CB(nP, nVertexMax, nLine)
    integer :: iLine, iX, iPoint, iLineAll, iMu
    real :: rNode, r1, r2, Weight
    real :: Gain_I(nP), LossRel_I(nP)
    ! Diagnostics: smallest f right after the update (before check_dist_neg
    ! drops a line); number of f < 0 values, flagged lines, unused lines
    real :: DistrMin
    integer :: nCount_I(3)
    character(len=*), parameter:: NameSub = 'interp_source_linenode'
    !--------------------------------------------------------------------------
    DistrMin = huge(1.0); nCount_I = 0
    ! Loop over all field lines on this processor
    LINE:do iLine = 1, nLine
       if(.not.Used_B(iLine)) CYCLE LINE
       iLineAll = iLine + iLineAll0
       ! Loop over all nodes along the valid field line
       do iX = 1, nVertex_B(iLine)
          ! Extract the radial coordinate of the current node
          rNode = State_VIB(R_, iX, iLine)*Io2Si_V(UnitX_)

          ! Nodes above the last shell centre are not changed
          if(rNode > RPerp_C(nRPerp)) CYCLE
          ! Nodes below RMinGhostPerp (inner boundary) are not changed either
          if(rNode < RMinGhostPerp*Io2Si_V(UnitX_)) CYCLE

          if(rNode <= RPerp_C(1)) then
             ! rNode is at or below the first shell centre: use shell 1
             Gain_I    = Gain_IV(:, iLineAll, 1)
             LossRel_I = LossRel_IV(:, iLineAll, 1)
          else
             ! RPerp_C(iPoint) < rNode <= RPerp_C(iPoint+1); this also
             ! covers rNode = RPerp_C(nRPerp), which used to get shell 1
             iPoint = count(RPerp_C(1:nRPerp) < rNode)
             r1 = RPerp_C(iPoint)
             r2 = RPerp_C(iPoint+1)
             ! Linear interpolation between r1 and r2, weight in [0,1]
             Weight = min(max((rNode - r1)/(r2 - r1), 0.0), 1.0)
             Gain_I = (1.0 - Weight)*Gain_IV(:, iLineAll, iPoint) &
                  + Weight*Gain_IV(:, iLineAll, iPoint+1)
             LossRel_I = (1.0 - Weight)*LossRel_IV(:, iLineAll, iPoint) &
                  + Weight*LossRel_IV(:, iLineAll, iPoint+1)
          end if
          ! Positivity-preserving update: relative loss, then absolute gain,
          ! with the same Gain_I and LossRel_I (from the mu-averaged f) for
          ! every mu; it keeps f >= 0 and, while 1 + dt*LossRel_I >= 0,
          ! changes the mu-average of f by dt*(Gain_I + LossRel_I*fAvg).
          ! Alternative (same mu-averaged change), same term for every mu:
          ! f(1:nP,iMu) = f(1:nP,iMu) + dt*(Gain_I + LossRel_I*fAvg)
          do iMu = 1, nMu
             Distribution_CB(1:nP, iMu, iX, iLine) = &
                  Distribution_CB(1:nP, iMu, iX, iLine)* &
                  max(0.0, 1.0 + dtIn_CB(:, iX, iLine)*LossRel_I) + &
                  dtIn_CB(:, iX, iLine)*Gain_I
          end do
       end do

       if(DoTestDPerp) then
          ! Own check right after the perp update of this line
          DistrMin = min(DistrMin, &
               minval(Distribution_CB(:, :, 1:nVertex_B(iLine), iLine)))
          nCount_I(1) = nCount_I(1) + &
               count(Distribution_CB(:, :, 1:nVertex_B(iLine), iLine) < 0.0)
       end if
       ! Check if the VDF includes negative values after Dperp
       call check_dist_neg(NameSub// &
            ' after perpendicular diffusion', 1, nVertex_B(iLine), iLine)
       if(.not.Used_B(iLine)) nCount_I(2) = nCount_I(2) + 1
       ! Was RETURN, which skipped the perp update of all remaining lines
       ! on this processor (rank-dependent results)
       if(IsDistNeg) CYCLE LINE
    end do LINE

    if(DoTestDPerp) then
       ! Lines not in use (dropped now or earlier) on this processor
       nCount_I(3) = count(.not.Used_B(1:nLine))
       if(nProc > 1) then
          call MPI_ALLREDUCE(MPI_IN_PLACE, DistrMin, 1, MPI_REAL, MPI_MIN, &
               iComm, iError)
          call MPI_ALLREDUCE(MPI_IN_PLACE, nCount_I, 3, MPI_INTEGER, &
               MPI_SUM, iComm, iError)
       end if
       if(iProc == 0) write(*,*) 'DPerp pos: after perp update: min f=', &
            DistrMin, ' nNegative=', nCount_I(1), ' nFlaggedLine=', &
            nCount_I(2), ' nUnusedLine=', nCount_I(3)
    end if
  end subroutine interp_source_linenode
  !============================================================================
end module SP_ModPerpDiffusion
!==============================================================================
