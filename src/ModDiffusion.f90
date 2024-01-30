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

  use ModNumConst, ONLY: cPi
  use ModConst,   ONLY: cAu, cLightSpeed, cGEV, cMu
  use SP_ModGrid, ONLY: Wave1_, Wave2_
  ! use SP_ModTurbulence, ONLY: update_spectrum
  use ModUtilities, ONLY: CON_stop

  implicit none

  SAVE

  PRIVATE

!!!!!!!!!!!!!!!!!!!!!!!!! Local parameters!!!!!!!!!!!!!!!
  ! Diffusion as in Li et al. (2003), doi:10.1029/2002JA009666
  logical, public :: UseFixedMFPUpstream = .false.
  real    :: MeanFreePath0InAu = 1.0

  ! Parameter characterizing cut-off wavenumber of turbulent spectrum:
  ! value of scale turbulence at 1 AU for any type (const or linear)
  real :: ScaleTurbulenceSI = 0.03 * cAu
  integer :: iScaleTurbulenceType
  integer, parameter :: Const_ = 0, Linear_ = 1

  ! Public members:
  public :: read_param, diffuse_distribution, advance_diffusion

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
       call read_var('ScaleTurbulence [AU] at 1 AU', ScaleTurbulenceSI)
       ScaleTurbulenceSI = ScaleTurbulenceSI * cAu
    end select
  end subroutine read_param
  !============================================================================
  subroutine diffuse_distribution(iLine, iEnd, iShock, Dt,        &
       Distribution_IIB, XyzSI_DI, nSI_I, BSI_I,             &
       DsSI_I, RadiusSi_I)!, DOuterSI_I, CoefDInnerSI_I)
    ! set up the diffusion coefficients
    ! diffuse the distribution function

    use ModConst, ONLY: cProtonMass, cGyroradius
    use SP_ModSize, ONLY: nVertexMax
    use SP_ModGrid, ONLY: MHData_VIB
    use SP_ModDistribution, ONLY: nP, SpeedSI_I, MomentumSI_I, DLogP
    ! use SP_ModTurbulence, ONLY: UseTurbulentSpectrum, set_dxx, Dxx

    ! Variables as inputs
    ! input Line, End (for how many particles), and Shock indices
    integer, intent(in) :: iLine, iEnd, iShock
    real, intent(in) :: Dt              ! Time step for diffusion
    real, intent(inout) :: Distribution_IIB(0:nP+1, 1:iEnd)
    real, intent(in) :: XyzSI_DI(3, 1:nVertexMax)
    real, intent(in), dimension(1:nVertexMax) :: nSI_I, BSI_I,  &
         DsSI_I, RadiusSi_I
    ! Variables declared in this subroutine
    integer :: iP, iVertex              ! loop variables
    ! Coefficients in the diffusion operator
    ! df/dt = DOuter * d(DInner * df/dx)/dx
    ! DOuter =BSI in the cell center
    ! DInner = DiffusionCoefficient/BSI at the face
    real, dimension(1:nVertexMax) :: DOuterSI_I, &
         DInnerSI_I, CoefDInnerSI_I
    ! real, dimension(1:nVertexMax) :: DInnerSI_I
    ! real, intent(in), dimension(1:nVertexMax) :: DOuterSI_I, CoefDInnerSI_I
    ! Full difference between DataInputTime and SPTime
    real, parameter :: DiffCoeffMinSI = 1.0E+04

    ! diffusion along the field line

    !--------------------------------------------------------------------------
    call set_coef_diffusion

    ! if using turbulent spectrum:
    ! set_dxx for diffusion along the field line
    ! if(UseTurbulentSpectrum) then
    !   call set_dxx(iEnd, nP, BSI_I(1:iEnd))
    ! end if

    MOMENTUM:do iP = 1, nP
       ! For each momentum account for dependence
       ! of the diffusion coefficient on momentum
       ! D\propto r_L*v\propto Momentum**2/TotalEnergy
       ! if (UseTurbulentSpectrum) then
       !   do iVertex=1, iEnd
       !      DInnerSI_I(iVertex) = Dxx(iVertex, iP,       &
       !           MomentumSI_I(iP), SpeedSI_I(iP),          &
       !           BSI_I(iVertex)) / BSI_I(iVertex)
       !   end do
       ! else
       ! Add v (= p*c^2/E_total in the relativistic case)
       ! and (p)^(1/3)
       DInnerSI_I(1:iEnd) = CoefDInnerSI_I(1:iEnd)     &
            *SpeedSI_I(iP)*(MomentumSI_I(iP))**(1.0/3)

       DInnerSI_I(1:iEnd) = max(DInnerSI_I(1:iEnd),    &
            DiffCoeffMinSI/DOuterSI_I(1:iEnd))
       ! end if

       call advance_diffusion(Dt, iEnd, DsSI_I(1:iEnd),   &
            Distribution_IIB(iP, 1:iEnd),                   &
            DOuterSI_I(1:iEnd), DInnerSI_I(1:iEnd))
    end do MOMENTUM

    ! if (UseTurbulentSpectrum) then
    !    call update_spectrum(iEnd,nP,MomentumSI_I,DLogP,   &
    !       XyzSI_DI(:,1:iEnd), DsSI_I(1:iEnd),             &
    !       Distribution_IIB(:,1:iEnd), BSI_I(1:iEnd),      &
    !       nSI_I(1:iEnd)*cProtonMass, Dt)
    ! end if
  contains
    !==========================================================================
    subroutine set_coef_diffusion
      ! set diffusion coefficient for the current line
      real, dimension(1:nVertexMax) :: ScaleSI_I
      real, parameter :: cCoef = 81./7/cPi/(2*cPi)**(2.0/3)

      !------------------------------------------------------------------------
      DOuterSI_I(1:iEnd) = BSI_I(1:iEnd)

      ! if(UseTurbulentSpectrum) RETURN

      ! precompute scale of turbulence along the line
      select case(iScaleTurbulenceType)
      case(Const_)
         ScaleSI_I(1:iEnd) = ScaleTurbulenceSI
      case(Linear_)
         ScaleSI_I(1:iEnd) = ScaleTurbulenceSI*RadiusSI_I(1:iEnd)/cAU
      end select
      ! Compute the diffusion coefficient without the contribution of
      ! v (velocity) and p (momentum), as v and p are different for
      ! different iP
      if(UseFixedMFPUpstream) then
         ! diffusion is different up- and down-stream
         ! Sokolov et al. 2004, paragraphs before and after eq (4)
         where(RadiusSI_I(1:iEnd) > 0.9 * RadiusSI_I(iShock))
            ! upstream: reset the diffusion coefficient to
            ! (1/3)*MeanFreePath0InAu[AU]*(R/1AU)*v*(pc/1GeV)^(1/3)
            ! see Li et al. (2003), doi:10.1029/2002JA009666
            ! ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
            ! 1/AU cancels with unit of Lambda0,no special attention needed;
            ! v (velocity) and (p)^(1/3) are calculated in momentum do loop
            CoefDInnerSI_I(1:iEnd) =                     &
                 (1.0/3)*MeanFreePath0InAU * RadiusSI_I(1:iEnd)&
                 *(cLightSpeed/cGEV)**(1.0/3)
         elsewhere
            CoefDInnerSI_I(1:iEnd) =  (cCoef/3)*BSI_I(1:iEnd)**2 /       &
                 (cMu*sum(MHData_VIB(Wave1_:Wave2_,1:iEnd,iLine),1))*    &
                 (ScaleSI_I(1:iEnd)**2*cGyroRadius/BSI_I(1:iEnd))**(1.0/3)
         end where
      else    ! .not.UseFixedMFPUpstream
         ! Sokolov et al., 2004: eq (4),
         ! note: Momentum = TotalEnergy * Vel / C**2
         ! Gyroradius = cGyroRadius * momentum / |B|
         ! DInner \propto (B/\delta B)**2*Gyroradius*Vel/|B|
         ! ------------------------------------------------------
         ! effective level of turbulence is different for different momenta:
         ! (\delta B)**2 \propto Gyroradius^(1/3)
         CoefDInnerSI_I(1:iEnd) =  (cCoef/3)*BSI_I(1:iEnd)**2 /       &
              (cMu*sum(MHData_VIB(Wave1_:Wave2_,1:iEnd,iLine),1))*    &
              (ScaleSI_I(1:iEnd)**2*cGyroRadius/BSI_I(1:iEnd))**(1.0/3)
      end if

      ! Add 1/B as the actual diffusion is D/B
      CoefDInnerSI_I(1:iEnd) = CoefDInnerSI_I(1:iEnd) / BSI_I(1:iEnd)

    end subroutine set_coef_diffusion
    !==========================================================================
  end subroutine diffuse_distribution
  !============================================================================
  subroutine advance_diffusion(Dt,n,Dist_I,F_I,DOuter_I,DInner_I)

    ! This routine solves the diffusion equation:
    !         f_t-D_outer(D_inner*f_x)_x=0,
    ! with zero Neumann boundary condition. The solution is advanced in time
    ! using fully implicit scheme.

    use ModNumConst, ONLY: cTiny

    real,   intent(in   ):: Dt     ! Time step
    integer,intent(in   ):: n      ! Number of meshes along the x-coordinate
    real,   intent(in   ):: Dist_I(n) ! Distance to the next mesh
    real,   intent(inout):: F_I(n) ! In:sol.to be advanced; Out:advanced sol
    ! Laplace multiplier and diffusion coefficient.
    real,   intent(in   ):: DOuter_I(n), DInner_I(n)

    ! Mesh spacing and face spacing.
    real                 :: DsMesh_I(2:n), DsFace_I(2:n-1)
    ! Main, upper, and lower diagonals.
    real, dimension(n)   :: Main_I,Upper_I,Lower_I, R_I
    integer:: i
    real:: Aux1,Aux2
    !--------------------------------------------------------------------------

    ! In M-FLAMPA D_I(i) is the distance between meshes i   and i+1
    ! while DsMesh_I(i) is the distance between centers of meshes
    ! i-1 and i. Therefore,
    do i=2,n
       DsMesh_I(i) = max(Dist_I(i-1),cTiny)
    end do

    ! Within the framework of finite volume method, the cell
    ! volume is used, which is proportional to  the distance between
    ! the faces bounding the volume with an index, i, which is half of
    ! sum of distance between meshes i-1 and i (i.e. D_I(i-1) and that
    ! between meshes i and i+1 (which is D_I(i)):
    do i=2,n-1
       DsFace_I(i) = max(0.5*(Dist_I(i) + Dist_I(i-1)),cTiny)
    end do
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
    Main_I = 1.0

    ! For i=1:
    Aux1 = Dt*DOuter_I(1)*0.50*(DInner_I(1)+DInner_I(2))/&
         DsMesh_I(2)**2
    Main_I( 1) = Main_I(1)+Aux1
    Upper_I(1) = -Aux1

    ! For i=2,n-1:
    do i=2,n-1
       Aux1 = Dt*DOuter_I(i)*0.50*(DInner_I(i  ) + DInner_I(i+1))/&
            (DsMesh_I(i+1)*DsFace_I(i))
       Aux2 = Dt*DOuter_I(i)*0.50*(DInner_I(i-1) + DInner_I(i  ))/&
            (DsMesh_I(i  )*DsFace_I(i))
       Main_I(i)  = Main_I(i) + Aux1 + Aux2
       Upper_I(i) = -Aux1
       Lower_I(i) = -Aux2
    end do

    ! For i=n:
    ! Aux2 = Dt*DOuter_I(n)*0.50*(DInner_I(n-1) + DInner_I(n))/&
    !     DsMesh_I(n)**2
    ! set free escaping at outerboundary for now
    Aux2=0.
    Main_I( n) = Main_I(n) + Aux2
    Lower_I(n) = -Aux2

    ! Update the solution from f^(n) to f^(n+1):
    R_I = F_I
    call tridiag(n,Lower_I,Main_I,Upper_I,R_I,F_I)

  end subroutine advance_diffusion
  !============================================================================
  subroutine tridiag(n, L_I, M_I, U_I, R_I, W_I)

    ! This routine solves three-diagonal system of equations:
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
    if (M_I(1)==0.0)&
         call CON_stop(' Error in tridiag: M_I(1)=0')
    Aux = M_I(1)
    W_I(1) = R_I(1)/Aux
    do j=2,n
       Aux_I(j) = U_I(j-1)/Aux
       Aux = M_I(j)-L_I(j)*Aux_I(j)
       if (Aux == 0.0) then
          write(*,*)'M_I(j), L_I(j), Aux_I(j) = ',&
               M_I(j),L_I(j),Aux_I(j)
          write(*,*)'  For j=',j
          call CON_stop('Tridiag failed')
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
