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

  use ModNumConst, ONLY: cPi, cTiny
  use ModConst,   ONLY: cAu, cLightSpeed, cGEV, cMu,     &
       Rsun, cProtonMass, cGyroradius
  use SP_ModSize, ONLY: nVertexMax
  use SP_ModGrid, ONLY: State_VIB, MHData_VIB, D_,       &
       R_, x_, y_, z_, Wave1_, Wave2_
  use SP_ModUnit, ONLY: UnitX_, Io2Si_V
  ! use SP_ModTurbulence, ONLY: update_spectrum
  use ModUtilities, ONLY: CON_stop

  implicit none

  SAVE

  PRIVATE
  logical, public :: UseDiffusion = .true.
  real, public :: DOuterSi_I(1:nVertexMax), CoefDInnerSi_I(1:nVertexMax)
  ! Local parameters!
  ! Diffusion as in Li et al. (2003), doi:10.1029/2002JA009666
  logical, public :: UseFixedMFPUpstream = .false.
  real    :: MeanFreePath0InAu = 1.0

  ! Parameter characterizing cut-off wavenumber of turbulent spectrum:
  ! value of scale turbulence at 1 AU for any type (const or linear)
  real :: ScaleTurbulenceSi = 0.03 * cAu
  integer, parameter :: Const_ = 0, Linear_ = 1
  integer :: iScaleTurbulenceType = Const_

  ! Public members:
  public :: read_param, diffuse_distribution, set_diffusion_coef

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
       call read_var('UseFixedMFPUpstream',UseFixedMFPUpstream)
       if(UseFixedMFPUpstream)then
          ! see Li et al. (2003), doi:10.1029/2002JA009666
          call read_var('MeanFreePath0 [AU]', MeanFreePath0InAu)
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
       call read_var('ScaleTurbulence [AU] at 1 AU', ScaleTurbulenceSi)
       ScaleTurbulenceSi = ScaleTurbulenceSi * cAu
    case('#DIFFUSION')
       call read_var('UseDiffusion', UseDiffusion)
    end select
  end subroutine read_param
  !============================================================================
  subroutine diffuse_distribution(iLine, nX, iShock, Dt, &
       nSi_I, BSi_I, LowerEndSpectrum_I, UpperEndSpectrum_I)
    ! set up the diffusion coefficients
    ! diffuse the distribution function
    !  use ModLinearSolver, ONLY: tridiag
    use SP_ModDistribution, ONLY: nP, SpeedSi_I, DLogP,  &
         MomentumSi_I, Distribution_IIB
    !  use SP_ModTurbulence, ONLY: UseTurbulentSpectrum, set_dxx, Dxx

    ! Variables as inputs
    ! input Line, End (for how many particles), and Shock indices
    integer, intent(in) :: iLine, nX, iShock
    real, intent(in) :: Dt              ! Time step for diffusion
    real, intent(in) :: nSi_I(1:nX), BSi_I(1:nX)
    ! Given spectrum of particles at low end (flare acceleration)
    real, intent(in), optional :: LowerEndSpectrum_I(nP)
    ! Given spectrum of particles at upper end (GCRs)
    real, intent(in), optional :: UpperEndSpectrum_I(nP)
    ! Variables declared in this subroutine
    real :: XyzSi_DI(3, 1:nX), DsSi_I(1:nX)
    integer :: iP, iVertex              ! loop variables
    ! Coefficients in the diffusion operator
    ! df/dt = DOuter * d(DInner * df/dx)/dx
    ! DOuter = BSi in the cell center
    ! DInner = DiffusionCoefficient/BSi at the face
    real  :: DInnerSi_I(1:nX)
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
    real   :: DsMesh_I(2:nX), DsFace_I(2:nX-1)
    ! Main, upper, and lower diagonals, source
    real   :: Main_I(nX), Upper_I(nX), Lower_I(nX), R_I(nX)
    real   :: Aux1, Aux2
    ! diffusion along the field line

    !--------------------------------------------------------------------------
    ! if using turbulent spectrum:
    ! set_dxx for diffusion along the field line
    ! if(UseTurbulentSpectrum) then
    !   call set_dxx(nX, nP, BSi_I(1:nX))
    ! end if

    ! In M-FLAMPA DsSi_I(i) is the distance between meshes i and i+1
    ! while DsMesh_I(i) is the distance between centers of meshes
    ! i-1 and i. Therefore,
    DsSi_I(1:nX) = State_VIB(D_,1:nX,iLine)*Io2Si_V(UnitX_)

    ! Within the framework of finite volume method, the cell
    ! volume is used, which is proportional to  the distance between
    ! the faces bounding the volume with an index, i, which is half of
    ! sum of distance between meshes i-1 and i (i.e. D_I(i-1) and that
    ! between meshes i and i+1 (which is D_I(i)):
    DsMesh_I(2:nX) = max(DsSi_I(1:nX-1), cTiny)

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
    DsFace_I(2:nX-1) = max(0.5*(DsSi_I(2:nX-1)+DsSi_I(1:nX-2)), cTiny)

    MOMENTUM:do iP = 1, nP
       ! For each momentum account for dependence
       ! of the diffusion coefficient on momentum
       ! D\propto r_L*v\propto Momentum**2/TotalEnergy
       ! if (UseTurbulentSpectrum) then
       !   do iVertex=1, nX
       !      DInnerSi_I(iVertex) = Dxx(iVertex, iP,       &
       !           MomentumSi_I(iP), SpeedSi_I(iP),        &
       !           BSi_I(iVertex)) / BSi_I(iVertex)
       !   end do
       ! else
       ! Add v (= p*c^2/E_total in the relativistic case)
       ! and (p)^(1/3)
       DInnerSi_I(1:nX) = CoefDInnerSi_I(1:nX)     &
            *SpeedSi_I(iP)*(MomentumSi_I(iP))**(1.0/3)

       DInnerSi_I(1:nX) = max(DInnerSi_I(1:nX),    &
            DiffCoeffMinSi/DOuterSi_I(1:nX))
       ! end if
       ! Now, we solve the matrix equation
       ! f^(n+1)_i-Dt*DOuter_I/DsFace_I*(&
       !     DInner_(i+1/2)*(f^(n+1)_(i+1)-f^(n+1)_i)/DsMesh_(i+1)-&
       !     DInner_(i-1/2)*(f^(n+1)_i -f^(n+1)_(i-1)/DsMesh_(i ))=f^n_i

       ! Set source term in the RHS:
       R_I = Distribution_IIB(iP, 1:nX, iLine)
       ! Set elements of tri-diagonal matrix in the LHS
       Main_I = 1.0

       ! For i=1:
       Aux1 = Dt*DOuterSi_I(1)*0.5*(DInnerSi_I(1)+DInnerSi_I(2))/&
            DsMesh_I(2)**2
       Main_I( 1) = Main_I(1) + Aux1
       Upper_I(1) = -Aux1
       if(present(LowerEndSpectrum_I))then
          Aux2 = Dt*DOuterSi_I(1)*DInnerSi_I(1)/DsMesh_I(2)**2
          Main_I(1) = Main_I(1) + Aux2
          R_I(1) = R_I(1) + Aux2*LowerEndSpectrum_I(iP)
       end if
       ! For i=2,n-1:
       do iVertex = 2, nX-1
          Aux1 = Dt*DOuterSi_I(iVertex)*0.5*(DInnerSi_I(iVertex  ) + &
               DInnerSi_I(iVertex+1))/(DsMesh_I(iVertex+1)*DsFace_I(iVertex))
          Aux2 = Dt*DOuterSi_I(iVertex)*0.5*(DInnerSi_I(iVertex-1) + &
               DInnerSi_I(iVertex  ))/(DsMesh_I(iVertex)*DsFace_I(iVertex))
          Main_I(iVertex)  = Main_I(iVertex) + Aux1 + Aux2
          Upper_I(iVertex) = -Aux1
          Lower_I(iVertex) = -Aux2
       end do

       ! For i=n:
       Aux2 = Dt*DOuterSi_I(nX)*0.5*(DInnerSi_I(nX-1) + DInnerSi_I(nX))/&
            DsMesh_I(nX)**2
       Main_I( nX) = Main_I(nX) + Aux2
       Lower_I(nX) = -Aux2
       if(present(UpperEndSpectrum_I))then
          Aux1 = Dt*DOuterSi_I(nX)*DInnerSi_I(nX)/DsMesh_I(nX)**2
          Main_I( nX) = Main_I(nX) + Aux1
          R_I(nX) = R_I(nX) + Aux1*UpperEndSpectrum_I(iP)
       end if
       ! Update the solution from f^(n) to f^(n+1):
       call tridiag(nX, Lower_I, Main_I, Upper_I, R_I,   &
            Distribution_IIB(iP, 1:nX, iLine))
    end do MOMENTUM

    ! if (UseTurbulentSpectrum) then
    !    XyzSi_DI(x_:z_,1:nX) = MhData_VIB(x_:z_,1:nX,iLine)*IO2Si_V(UnitX_)
    !    call update_spectrum(nX,nP,MomentumSi_I,DLogP, &
    !       XyzSi_DI(:,1:nX), DsSi_I(1:nX),             &
    !       Distribution_IIB(:,1:nX), BSi_I(1:nX),      &
    !       nSi_I(1:nX)*cProtonMass, Dt)
    ! end if

  end subroutine diffuse_distribution
  !============================================================================
  subroutine set_diffusion_coef(iLine, nX, iShock, BSi_I)
    ! set diffusion coefficient for the current line

    integer, intent(in) :: iLine, nX, iShock
    real, intent(in)    :: BSI_I(1:nX)
    real                :: ScaleSi_I(1:nX), RadiusSi_I(1:nX)
    real, parameter     :: cCoef = 81.0/7/cPi/(2*cPi)**(2.0/3)
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
    ! Compute the diffusion coefficient without the contribution of
    ! v (velocity) and p (momentum), as v and p are different for
    ! different iP
    if(UseFixedMFPUpstream) then
       ! diffusion is different up- and down-stream
       ! Sokolov et al. 2004, paragraphs before and after eq (4)
       where(RadiusSi_I(1:nX) > 0.9 * RadiusSi_I(iShock))
          ! upstream: reset the diffusion coefficient to
          ! (1/3)*MeanFreePath0InAu[AU]*(R/1AU)*v*(pc/1GeV)^(1/3)
          ! see Li et al. (2003), doi:10.1029/2002JA009666
          ! ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
          ! 1/AU cancels with unit of Lambda0,no special attention needed;
          ! v (velocity) and (p)^(1/3) are calculated in momentum do loop
          CoefDInnerSi_I(1:nX) = (1.0/3)*MeanFreePath0InAu *      &
               RadiusSi_I(1:nX) * (cLightSpeed/cGEV)**(1.0/3)
       elsewhere
          CoefDInnerSi_I(1:nX) = (cCoef/3)*BSi_I(1:nX)**2 /       &
               (cMu*sum(MHData_VIB(Wave1_:Wave2_,1:nX,iLine),1))* &
               (ScaleSi_I(1:nX)**2*cGyroRadius/BSi_I(1:nX))**(1.0/3)
       end where
    else  ! .not.UseFixedMFPUpstream
       ! Sokolov et al., 2004: eq (4),
       ! note: Momentum = TotalEnergy * Vel / C**2
       ! Gyroradius = cGyroRadius * momentum / |B|
       ! DInner \propto (B/\delta B)**2*Gyroradius*Vel/|B|
       ! ------------------------------------------------------
       ! effective level of turbulence is different for different momenta:
       ! (\delta B)**2 \propto Gyroradius^(1/3)
       CoefDInnerSi_I(1:nX) = (cCoef/3)*BSi_I(1:nX)**2 /          &
            (cMu*sum(MHData_VIB(Wave1_:Wave2_,1:nX,iLine),1))*    &
            (ScaleSi_I(1:nX)**2*cGyroRadius/BSi_I(1:nX))**(1.0/3)
    end if

    ! Add 1/B as the actual diffusion is D/B
    CoefDInnerSi_I(1:nX) = CoefDInnerSi_I(1:nX) / BSi_I(1:nX)

  end subroutine set_diffusion_coef
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
    real:: Aux,Aux_I(2:n)
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
end module SP_ModDiffusion
!==============================================================================
