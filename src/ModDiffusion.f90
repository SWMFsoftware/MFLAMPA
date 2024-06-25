!  Copyright (C) 2002 Regents of the University of Michigan,
!  portions used with permission
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module SP_ModDiffusion

  ! Solve the diffusion term in the Parker equation
  ! Revision history:
  ! Prototype:Sp/FLAMPA/src/SP_ModDiffusion.f90
  ! Adapted for the use in MFLAMPA (Dist_I is an input paramater,
  ! fixed contributions to M_I in the end points)-D.Borovikov, 2017
  ! Updated (identation, comments):  I.Sokolov, Dec.17, 2017

  use ModNumConst,  ONLY: cPi, cTwoPi, cTiny
  use ModConst,     ONLY: cAu, cLightSpeed, cGeV, cMu, Rsun, cGyroRadius
  use SP_ModSize,   ONLY: nVertexMax
  use SP_ModDistribution, ONLY: MomentumInjSi, nP
  use SP_ModGrid,   ONLY: State_VIB, MHData_VIB, D_, R_, Wave1_, Wave2_
  use SP_ModUnit,   ONLY: UnitX_, Io2Si_V
  use ModUtilities, ONLY: CON_stop

  implicit none

  SAVE

  PRIVATE
  logical, public :: UseDiffusion = .true.
  real,    public :: DOuterSi_I(1:nVertexMax), CoefDInnerSi_I(1:nVertexMax)
  ! Local parameters!
  ! Diffusion as in Li et al. (2003), doi:10.1029/2002JA009666
  logical, public :: UseFixedMeanFreePathUpstream = .false.
  real            :: MeanFreePath0InAu = 1.0

  ! Parameter characterizing cut-off wavenumber of turbulent spectrum:
  ! value of scale turbulence at 1 AU for any type (const or linear)
  real            :: ScaleTurbulenceSi = 0.03*cAu
  integer, parameter :: Const_ = 0, Linear_ = 1
  integer         :: iScaleTurbulenceType = Const_

  ! Public members:
  public          :: read_param, diffuse_distribution, set_diffusion_coef
  interface diffuse_distribution
     module procedure diffuse_distribution_s    ! Global time step Dt
     module procedure diffuse_distribution_arr  ! DtLocal_I array
  end interface diffuse_distribution

contains
  !============================================================================
  subroutine read_param(NameCommand)
    use ModReadParam, ONLY: read_var
    character(len=*), intent(in) :: NameCommand ! From PARAM.in
    character(len=8) :: StringScaleTurbulenceType
    character(len=*), parameter:: NameSub = 'read_param'
    !--------------------------------------------------------------------------
    select case(NameCommand)
    case('#USEFIXEDMFPUPSTREAM')
       call read_var('UseFixedMeanFreePathUpstream',&
            UseFixedMeanFreePathUpstream)
       if(UseFixedMeanFreePathUpstream)then
          ! see Li et al. (2003), doi:10.1029/2002JA009666
          call read_var('MeanFreePath0InAu', MeanFreePath0InAu)
       end if
    case('#SCALETURBULENCE')
       ! cut-off wavenumber of turbulence spectrum: k0 = 2 cPi / Scale
       call read_var('ScaleTurbulenceType', StringScaleTurbulenceType)
       select case(StringScaleTurbulenceType)
       case('const', 'constant')
          iScaleTurbulenceType = Const_
       case('linear')
          iScaleTurbulenceType = Linear_
       case default
          call CON_stop(NameSub//": unknown scale turbulence type")
       end select
       call read_var('ScaleTurbulence1AU', ScaleTurbulenceSi)
       ScaleTurbulenceSi = ScaleTurbulenceSi * cAu
    case('#DIFFUSION')
       call read_var('UseDiffusion', UseDiffusion)
    end select
  end subroutine read_param
  !============================================================================
  subroutine diffuse_distribution_s(iLine, nX, iShock, Dt, &
       nSi_I, BSi_I, LowerEndSpectrum_I, UpperEndSpectrum_I)
    ! diffuse the distribution function with scalar/global Dt

    ! Variables as inputs
    ! input Line, End (for how many particles), and Shock indices
    integer, intent(in) :: iLine, nX, iShock
    real, intent(in) :: Dt              ! Time step for diffusion
    real, intent(in) :: nSi_I(1:nX), BSi_I(1:nX)
    ! Given spectrum of particles at low end (flare acceleration)
    real, intent(in), optional :: LowerEndSpectrum_I(nP)
    ! Given spectrum of particles at upper end (GCRs)
    real, intent(in), optional :: UpperEndSpectrum_I(nP)
    ! LOCAL VARS
    real :: DtFake_C(nP, nX)
    !--------------------------------------------------------------------------
    DtFake_C = Dt
    call diffuse_distribution_arr(iLine, nX, iShock, DtFake_C, &
         nSi_I, BSi_I, LowerEndSpectrum_I, UpperEndSpectrum_I)
  end subroutine diffuse_distribution_s
  !============================================================================
  subroutine diffuse_distribution_arr(iLine, nX, iShock, DtLocalIn_II,  &
       nSi_I, BSi_I, LowerEndSpectrum_I, UpperEndSpectrum_I)
    ! diffuse the distribution function with vector/local Dt
    use SP_ModDistribution, ONLY: SpeedSi_I, Momentum_I, Distribution_CB
    use SP_ModTurbulence,   ONLY: UseTurbulentSpectrum, set_dxx, Dxx

    ! Variables as inputs
    ! input Line, End (for how many particles), and Shock indices
    integer, intent(in) :: iLine, nX, iShock
    real, intent(in) :: DtLocalIn_II(nP, nX)  ! Local time step for diffusion
    real, intent(in) :: nSi_I(1:nX), BSi_I(1:nX)
    ! Given spectrum of particles at low end (flare acceleration)
    real, intent(in), optional :: LowerEndSpectrum_I(nP)
    ! Given spectrum of particles at upper end (GCRs)
    real, intent(in), optional :: UpperEndSpectrum_I(nP)
    ! Variables declared in this subroutine
    integer :: iP                ! loop variables
    ! For an optimized loop, we need to change the (nP, nX) into (nX, nP)
    ! for the two-dimensional local time step DtLocal_II array
    real :: DtLocal_II(nX, nP)   ! Local time step for diffusion
    ! Coefficients in the diffusion operator
    ! df/dt = DOuter * d(DInner * df/dx)/dx
    ! DOuter = BSi in the cell center
    ! DInner = DiffusionCoefficient/BSi at the face
    real :: DInnerSi_II(nX, nP)  ! Calculate once and use later
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
    real, parameter :: DiffCoeffMinSi = 1.0E+04*Rsun
    ! Mesh spacing and face spacing.
    real :: DsMesh_I(2:nX), DsFace_I(2:nX-1)
    ! Main, upper, and lower diagonals, source
    real :: Main_I(nX), Upper_I(nX), Lower_I(nX), Res_I(nX)
    ! Auxiliary arrays for the tri-diagonal arrays
    real :: Aux1_I(nX), Aux2_I(nX)
    !--------------------------------------------------------------------------
    ! diffusion along the field line
    ! Set up the local time step: the reason is that, we loop through
    ! each iX at the fixed iP, so we will visit each DtLocal_II(iX, iP)
    ! in the loop with the inner-most iX and outer-most iP. It will
    ! visit and put the data in cache that are stored adjacently.
    DtLocal_II = transpose(DtLocalIn_II)

    ! if using turbulent spectrum: set_dxx for diffusion along the field line
    if(UseTurbulentSpectrum) call set_dxx(nX, BSi_I(1:nX))

    ! In M-FLAMPA DsMesh_I(i) is the distance between centers of meshes
    ! i-1 and i. Therefore,
    DsMesh_I(2:nX) = max(State_VIB(D_,1:nX-1,iLine)*Io2Si_V(UnitX_), cTiny)
    ! Note: The max(..., cTiny) in DsMesh_I and DsFace_I will be removed later

    ! Within the framework of finite volume method, the cell volume
    ! is used, which is proportional to the distance between the faces
    ! bounding the volume with an index, i, which is half of sum of
    ! distance between meshes i-1 and i (i.e. DsMesh_I(i) and that
    ! between meshes i and i+1 (which is DsMesh_I(i+1)):
    DsFace_I(2:nX-1) = max(0.5*(DsMesh_I(3:nX)+DsMesh_I(2:nX-1)), cTiny)
    ! In flux coordinates, the control volume associated with the given
    ! cell has a cross-section equal to (Magnetic Flux)/B, where the flux
    ! is a constant along the magnetic field line, set to one hereafter.
    ! Therefore, the diffusion equation has a following form:
    ! (DsFace_i/B_i)(f^(n+1) - f^n) = Flux_(i-1/2) - Flux_(i+1/2),
    ! where the particle density flux should be multiplied by the
    ! cross-section area too (magnetic flux factor is one!):
    !  Flux_(i-1/2) = (diffusion coefficient)_(i-1/2)/B_(i-1/2)*&
    !                 (f^(n+1)_(i-1) - f^(n+1)_i),
    !  The multiplier, DsFace_i/B_i, is denoted as DsFace_i/DOuter_i
    !  The face-centered combination,

    ! For each momentum, the dependence of the diffusion coefficient
    ! on momentum is D \propto r_L*v \propto Momentum**2/TotalEnergy
    if(UseTurbulentSpectrum) then
       DInnerSi_II = Dxx(nX, BSi_I(1:nX))/ &
            spread(BSi_I(1:nX), DIM=2, NCOPIES=nP)
    else
       ! Add v (= p*c^2/E_total in the relativistic case) and p^(1/3)
       DInnerSi_II = spread(CoefDInnerSi_I(1:nX), DIM=2, NCOPIES=NP) &
            * spread(SpeedSi_I(1:nP)*Momentum_I(1:nP)**(1.0/3),      &
            DIM=1, NCOPIES=NX)
       DInnerSi_II = max(DInnerSi_II, &
            DiffCoeffMinSi/spread(DOuterSi_I(1:nX), DIM=2, NCOPIES=NP))
    end if

    MOMENTUM:do iP = 1, nP
       ! Now, we solve the matrix equation
       ! f^(n+1)_i-Dt*DOuter_I/DsFace_I*(&
       !     DInner_(i+1/2)*(f^(n+1)_(i+1)-f^(n+1)_i)/DsMesh_(i+1)-&
       !     DInner_(i-1/2)*(f^(n+1)_i -f^(n+1)_(i-1)/DsMesh_(i ))=f^n_i
       ! Set source term in the RHS:
       Res_I = Distribution_CB(iP, 1, 1:nX, iLine)
       ! Set elements of tri-diagonal matrix in the LHS
       Main_I = 1.0; Lower_I = 0.0; Upper_I = 0.0

       ! For i=1:
       Aux1_I(1)  = DtLocal_II(1,iP)*DOuterSi_I(1)*   &
            0.5*(DInnerSi_II(1,iP) + DInnerSi_II(2,iP))/DsMesh_I(2)**2
       Main_I(1)  = Main_I(1) + Aux1_I(1)
       Upper_I(1) = -Aux1_I(1)
       if(present(LowerEndSpectrum_I)) then
          ! With the given value for f_0 behind the boundary,
          ! the above scheme reads:
          ! f^(n+1)_1-Dt*DOuter_I/DsFace_I*(&
          !     DInner_(3/2)*(f^(n+1)_2-f^(n+1)_1)/DsMesh_(2)-&
          !     DInner_(1/2)*(f^(n+1)_1 -f_0/DsMesh_(1))=f^n_1
          Aux2_I(1) = DtLocal_II(1,iP)*DOuterSi_I(1)  &
               *DInnerSi_II(1,iP)/DsMesh_I(2)**2
          ! With these regards, Aux2 is added to Main_I(1)...
          Main_I(1) = Main_I(1) + Aux2_I(1)
          ! ...while the given Aux2*f_0 should moved to the RHS and summed up
          ! with the source:
          Res_I(1)  = Res_I(1) + Aux2_I(1)*LowerEndSpectrum_I(iP)
       end if

       ! For i=2, n-1:
       Aux1_I(2:nX-1) = DtLocal_II(2:nX-1,iP)*DOuterSi_I(2:nX-1)* &
            0.5*(DInnerSi_II(2:nX-1,iP) + DInnerSi_II(3:nX, iP))/ &
            (DsMesh_I(3:nX)*DsFace_I(2:nX-1))
       Aux2_I(2:nX-1) = DtLocal_II(2:nX-1,iP)*DOuterSi_I(2:nX-1)* &
            0.5*(DInnerSi_II(1:nX-2,iP) + DInnerSi_II(2:nX-1,iP))/&
            (DsMesh_I(2:nX-1)*DsFace_I(2:nX-1))
       Main_I(2:nX-1)  = Main_I(2:nX-1) + Aux1_I(2:nX-1) + Aux2_I(2:nX-1)
       Upper_I(2:nX-1) = -Aux1_I(2:nX-1)
       Lower_I(2:nX-1) = -Aux2_I(2:nX-1)

       ! For i=n:
       ! Version before Nov. 2023:
       ! Aux2 = Dt*DOuter_I(n)*0.50*(DInner_I(n-1) + DInner_I(n))/&
       !     DsMesh_I(n)**2
       !
       ! After Nov. 2023: set free escaping at outer boundary for now
       ! Aux2 = 0.0
       ! In both these versions:
       ! Main_I( n) = Main_I(n) + Aux2
       ! Lower_I(n) = -Aux2
       ! So, effectively for the version after Nov. 2023
       ! Main_I(n) = 1; Lower_I(n) = 0 (equivalently to doing nothing)
       ! For backward compatibility, keep this option for UseUpperEndBc=.false.
       if(present(UpperEndSpectrum_I)) then
          Aux2_I(nX)  = DtLocal_II(nX,iP)*DOuterSi_I(nX)* &
               0.5*(DInnerSi_II(nX-1,iP) + DInnerSi_II(nX,iP))/DsMesh_I(nX)**2
          Main_I(nX)  = Main_I(nX) + Aux2_I(nX)
          Lower_I(nX) = -Aux2_I(nX)
          Aux1_I(nX)  = DtLocal_II(nX,iP)*DOuterSi_I(nX)* &
               DInnerSi_II(nX,iP)/DsMesh_I(nX)**2
          Main_I(nX)  = Main_I(nX) + Aux1_I(nX)
          Res_I(nX)   = Res_I(nX) + Aux1_I(nX)*UpperEndSpectrum_I(iP)
       end if

       ! Update the solution from f^(n) to f^(n+1):
       call tridiag(nX, Lower_I, Main_I, Upper_I, Res_I, &
            Distribution_CB(iP, 1, 1:nX, iLine))
    end do MOMENTUM
  end subroutine diffuse_distribution_arr
  !============================================================================
  subroutine set_diffusion_coef(iLine, nX, iShock, BSi_I)
    ! set diffusion coefficient for the current line
    use SP_ModBc, ONLY: SpectralIndex

    integer, intent(in) :: iLine, nX, iShock
    real, intent(in)    :: BSi_I(1:nX)
    real                :: ScaleSi_I(1:nX), RadiusSi_I(1:nX)
    real, parameter     :: cCoef = 81.0/7/cPi/cTwoPi**(2.0/3)
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
    ! Compute diffusion coefficient without the contribution of v (velocity)
    ! and p (momentum), as v and p are different for different iP
    if(UseFixedMeanFreePathUpstream) then
       ! diffusion is different up- and down-stream
       ! Sokolov et al. 2004, paragraphs before and after eq (4)
       where(RadiusSi_I(1:nX) > 0.9*RadiusSi_I(iShock))
          ! upstream: reset the diffusion coefficient to
          ! (1/3)*MeanFreePath0InAu[AU]*(R/1AU)*v*(pc/1GeV)^(1/3)
          ! see Li et al. (2003), doi:10.1029/2002JA009666
          ! ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
          ! 1/AU cancels with unit of Lambda0, no special attention needed;
          ! v (velocity) and (p)^(1/3) are calculated in momentum do loop
          CoefDInnerSi_I(1:nX) = (1.0/3)*MeanFreePath0InAu *      &
               RadiusSi_I(1:nX)*(cLightSpeed*MomentumInjSi/cGeV)**(1.0/3)
       elsewhere
          CoefDInnerSi_I(1:nX) = (cCoef/3)*BSi_I(1:nX)**2 /       &
               (cMu*sum(MHData_VIB(Wave1_:Wave2_,1:nX,iLine),1))* &
               (ScaleSi_I(1:nX)**2*cGyroRadius*MomentumInjSi/BSi_I(1:nX))&
               **(1.0/3)
       end where
    else
       ! Sokolov et al., 2004: eq (4),
       ! Gyroradius = cGyroRadius * momentum / |B|
       ! DInner \propto (B/\delta B)**2*Gyroradius*Vel/|B|
       ! ------------------------------------------------------
       ! effective level of turbulence is different for different momenta:
       ! (\delta B)**2 \propto Gyroradius^(1/3)
       CoefDInnerSi_I(1:nX) = (cCoef/3)*BSi_I(1:nX)**2 /          &
            (cMu*sum(MHData_VIB(Wave1_:Wave2_,1:nX,iLine),1))*    &
            (ScaleSi_I(1:nX)**2*cGyroRadius*MomentumInjSi/BSi_I(1:nX))&
            **(1.0/3)
    end if

    ! Add 1/B as the actual diffusion is D/B
    CoefDInnerSi_I(1:nX) = CoefDInnerSi_I(1:nX)/BSi_I(1:nX)

  end subroutine set_diffusion_coef
  !============================================================================
  subroutine tridiag(n, Lower_I, Mean_I, Upper_I, Res_I, W_I)

    ! Solve tri-diagonal system of equations:
    !  ||m_1 u_1  0....        || ||w_1|| ||r_1||
    !  ||l_2 m_2 u_2...        || ||w_2|| ||r_2||
    !  || 0  l_3 m_3 u_3       ||.||w_3||=||r_3||
    !  ||...                   || ||...|| ||...||
    !  ||.............0 l_n m_n|| ||w_n|| ||r_n||
    ! From: Numerical Recipes, Chapter 2.6, p.40.

    ! input parameters
    integer, intent(in) :: n
    real,    intent(in) :: Lower_I(n), Mean_I(n), Upper_I(n), Res_I(n)
    ! Output parameters
    real,    intent(out):: W_I(n)
    ! Misc
    integer:: j
    real:: Aux,Aux_I(2:n)
    !--------------------------------------------------------------------------
    if(Mean_I(1) == 0.0) call CON_stop('Error in tridiag: Mean_I(1)=0')
    Aux = Mean_I(1)
    W_I(1) = Res_I(1)/Aux
    do j = 2, n
       Aux_I(j) = Upper_I(j-1)/Aux
       Aux = Mean_I(j) - Lower_I(j)*Aux_I(j)
       if(Aux == 0.0) then
          write(*,*) 'M_I(j), L_I(j), Aux_I(j) = ',&
               Mean_I(j),Lower_I(j),Aux_I(j)
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
