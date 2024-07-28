!  Copyright (C) 2002 Regents of the University of Michigan,
!  portions used with permission
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module SP_ModBc

  ! The module sets up boundary conditions
  use ModNumConst,        ONLY: cPi
  use ModCosmicRay,       ONLY: local_interstellar_spectrum,          &
       TypeLisBc, UseModulationPot, ModulationPot
  use SP_ModDistribution, ONLY: nP, nMu, Distribution_CB, Momentum_I, &
       MomentumInjSi, Background_I
  use SP_ModGrid,         ONLY: MhData_VIB, NoShock_, nWidth, T_, X_, Z_
  use SP_ModUnit,         ONLY: kinetic_energy_to_momentum,           &
       UnitEnergy_, UnitX_, Io2Si_V
  use ModUtilities,       ONLY: CON_stop
  implicit none

  SAVE

  PRIVATE ! except

  ! Public members:
  public:: read_param
  public:: set_momentum_bc   ! Set the Distribution_CB Bc at min/max energy
  public:: set_upper_end_bc  ! Set the UpperEndBc_I Bc at the far end point
  public:: set_lower_end_bc  ! Set the LowerEndBc_I Bc at the near-Sun end
  public:: set_lower_end_vdf ! Set the Distribution_CB Bc at the near-Sun end
  public:: set_VDF           ! Set the VDF Bc, typically with 2 ghost cells
  interface set_VDF
     module procedure set_VDF2 ! VDF2: along P^3/3 and s_L
     module procedure set_VDF3 ! VDF3: along P^3/3, mu and s_L
  end interface set_VDF

  ! Boundary condition at the injection energy
  ! Injection efficiency and assumed spectral index with the energy
  ! range k_B*T_i < Energy < EnergyInjection, to be read from PARAM.in
  real, public      :: CoefInj = 0.25, SpectralIndex = 5.0
  real, parameter   :: CoefInjTiny = 2.5E-11 ! Before Nov 2023: set to be 0.0
  ! Injection coefficient at lower end boundary, should set 0.0 for GCRs
  real, public      :: CoefInjLowBc = 0.25

  ! Lower end BC, set at the firsr point along the field line
  logical, public   :: UseLowerEndBc = .true.
  ! Type of lower end BC: float, escape, inject
  character(LEN=6)  :: TypeLowerEndBc = 'inject'
  ! Index of the left lower end BC
  integer, public   :: iStart = 1
  integer, parameter:: &
       iStartUseLeft_   = 1, & ! when UseLowerEndBc = .true.
       iStartNoUseLeft_ = 2    ! when UseLowerEndBc = .false.

  ! Upper end BC, set at the last point along the field line
  logical, public   :: UseUpperEndBc = .false.
  ! Type of upper end BC: none, float, escape, lism
  character(LEN=6)  :: TypeUpperEndBc = 'none'
  real, public      :: UpperEndBc_I(1:nP), LowerEndBc_I(0:nP+1)
contains
  !============================================================================
  subroutine read_param(NameCommand)

    use ModReadParam, ONLY: read_var
    use ModUtilities, ONLY: lower_case

    character(len=*), intent(in):: NameCommand ! From PARAM.in
    character(len=*), parameter:: NameSub = 'read_param'
    !--------------------------------------------------------------------------
    select case(NameCommand)
    case('#INJECTION')
       call read_var('SpectralIndex',   SpectralIndex)
       call read_var('EfficiencyGlob',  CoefInj)
       call read_var('EfficiencyLowBc', CoefInjLowBc)
    case('#LOWERENDBC')
       ! Read whether to use LowerLowBc
       call read_var('UseLowerEndBc', UseLowerEndBc)
       if(UseLowerEndBc) then
          iStart = iStartUseLeft_
          ! Read the type of LowerEndBc and Select
          call read_var('TypeLowerEndBc', TypeLowerEndBc)
          call lower_case(TypeLowerEndBc)
       else
          iStart = iStartNoUseLeft_
       end if
    case('#UPPERENDBC')
       ! Read whether to use UpperLowBc
       call read_var('UseUpperEndBc', UseUpperEndBc)
       if(UseUpperEndBc) then
          ! Read the type of UppenEndBc and Select
          call read_var('TypeUpperEndBc', TypeUpperEndBc)
          call lower_case(TypeUpperEndBc)

          select case(trim(TypeUpperEndBc))
          case('none', 'f', 'false')
             ! Reset UseUpperEndBc
             UseUpperEndBc = .false.
          case('float', 'escape')
             ! Do nothing
          case('lism')
             ! We want to read the type of LIS here
             call read_var('TypeLisBc', TypeLisBc)
             call lower_case(TypeLisBc)
             ! Read whether using ModulationPot to get GCR spectrum at ~1 AU
             call read_var('UseModulationPot', UseModulationPot)
             if(UseModulationPot) call read_var('ModulationPot', ModulationPot)
          case default
             call CON_stop(NameSub//&
                  ': Unknown type of upper end BC '//TypeUpperEndBc)
          end select
       end if
    case default
       call CON_stop(NameSub//': Unknown command '//NameCommand)
    end select
  end subroutine read_param
  !============================================================================
  subroutine set_momentum_bc(iLine, iEnd, nSi_I, iShock)
    ! Set boundary conditions for Distribution_CB on grid point,
    ! on the given field line.

    integer, intent(in) :: iLine, iEnd, iShock
    real,    intent(in) :: nSi_I(1:iEnd)
    ! local variables
    integer :: iVertex     ! loop variable
    real    :: MomentumSi  ! Momentum for the thermal energy k_BTi
    real    :: CoefInjLocal, DistributionBc
    !--------------------------------------------------------------------------

    do iVertex = 1, iEnd
       ! injection(Ti, Rho), see Sokolov et al., 2004, eq (3)
       ! f = CoefInj/2/pi * N / (2*m*T_p)^(3/2) * ((2*m*T_p)^(1/2)/p_inj)^5
       !   = CoefInj/2/pi * N / p^3 * (p/p_inj)^5
       ! where p = sqrt(2*m*T_p) is the momentum of thermal ion
       CoefInjLocal = CoefInjTiny
       MomentumSi   = kinetic_energy_to_momentum( &
            MhData_VIB(T_,iVertex,iLine)*Io2Si_V(UnitEnergy_))

       DistributionBc = (SpectralIndex-3)/(4*cPi)   &
            * MomentumInjSi**2*Io2Si_V(UnitEnergy_) &
            * nSi_I(iVertex)/MomentumSi**3          &
            * (MomentumSi/MomentumInjSi)**SpectralIndex

       if(iShock /= NoShock_ .and. iVertex <= iShock + nWidth .and.  &
            iVertex >= iShock - nWidth) CoefInjLocal = CoefInj

       Distribution_CB(0, :, iVertex, iLine) = DistributionBc*CoefInjLocal
    end do
    ! Modified for lower end boundary
    Distribution_CB(0, :, 1, iLine) = &
         Distribution_CB(0, :, 1, iLine)/CoefInj*CoefInjLowBc

  end subroutine set_momentum_bc
  !============================================================================
  subroutine set_upper_end_bc(iLine, iEnd)
    ! Set boundary condition at the last grid point on the given field line.
    ! Assign the calculated BC to UpperEndBc_I.

    integer, intent(in) :: iLine, iEnd
    real :: XyzSi_D(3)                          ! Where to set BC
    !--------------------------------------------------------------------------
    select case(trim(TypeUpperEndBc))
    case('float')
       UpperEndBc_I = Distribution_CB(1:nP,nMu,iEnd,iLine)
    case('escape')
       UpperEndBc_I = Background_I(1:nP)
    case('lism')
       XyzSi_D = MhData_VIB(X_:Z_,iEnd,iLine)*Io2Si_V(UnitX_)
       call local_interstellar_spectrum(&
            nP = nP,                         &  ! # of grid points
            MomentumSi_I = Momentum_I(1:nP)* &
            MomentumInjSi,                   &  ! momentum (SI) in grid points
            XyzSi_D = XyzSi_D,               &  ! Coords
            DistTimesP2Si_I = UpperEndBc_I)
       ! Now, in UpperEndBc_I there is Distribution[Si]*Momentum[Si]**2
       ! Our Momentum_I is MomentumSi_I/MomentumInjSi
       ! So, UpperEndBc_I is Distribution[Si]*MomentumInjSi**2*Momentum_I**2
       ! The distribution used in our code is
       ! Distribution[Si]*MomentumInjSi**2*Io2Si_V(UnitEnergy_)
       UpperEndBc_I = (UpperEndBc_I/Momentum_I(1:nP)**2)*Io2Si_V(UnitEnergy_)
    end select
  end subroutine set_upper_end_bc
  !============================================================================
  subroutine set_lower_end_bc(iLine)
    ! Set boundary condition at the zeroth grid point on the given field line.
    ! Assign the calculated BC to LowerEndBc_I.

    integer, intent(in) :: iLine
    character(len=*), parameter:: NameSub = 'set_lower_end_bc'
    !--------------------------------------------------------------------------
    select case(trim(TypeLowerEndBc))
    case('float')
       LowerEndBc_I = Distribution_CB(0:nP+1, nMu, 1, iLine)
    case('escape')
       LowerEndBc_I = Background_I
    case('inject')
       LowerEndBc_I = Distribution_CB(0, 1, 1, iLine) &
            /Momentum_I(0:nP+1)**SpectralIndex
    case default
       call CON_stop(NameSub//&
            ': Unknown type of lower end BC '//TypeLowerEndBc)
    end select
  end subroutine set_lower_end_bc
  !============================================================================
  subroutine set_lower_end_vdf(iLine)
    ! Set boundary condition for Distribution_CB at the zeroth grid point,
    ! namely the footpoint on the solar surface, on the current field line.

    integer, intent(in) :: iLine
    !--------------------------------------------------------------------------
    ! Set the left boundary condition of VDF for diffusion
    Distribution_CB(1:nP+1, 1, 1, iLine) = &
         Distribution_CB(0, 1, 1, iLine)/Momentum_I(1:nP+1)**SpectralIndex
  end subroutine set_lower_end_vdf
  !============================================================================
  subroutine set_VDF2(iLine, nX, VDF_G)
    ! We need the VDF on the extended grid with two layers of ghost cells, to
    ! solve the second-order accurate scheme. Add solution in physical cells
    ! and in a single layer of the ghost cells along the momentum coordinate,
    ! plus two layers of ghost cells along the pitch-angle coordinates.
    ! This function is for the 2-dimensional VDF, with momentum (P^3/3)
    ! and line coordinates (s_L).

    integer, intent(in) :: iLine     ! Indices of line and shock
    integer, intent(in) :: nX        ! Number of meshes along s_L axis
    real, intent(inout) :: VDF_G(-1:nP+2, -1:nX+2)
    real :: VDF3_G(-1:nP+2, -1:nMu+2, -1:nX+2) ! nMu = 1 for VDF2
    !--------------------------------------------------------------------------

    call set_VDF3(iLine, nX, VDF3_G)
    VDF_G = VDF3_G(:, nMu, :) ! nMu = 1 for VDF2

  end subroutine set_VDF2
  !============================================================================
  subroutine set_VDF3(iLine, nX, VDF_G)
    ! We need the VDF on the extended grid with two layers of ghost cells, to
    ! solve the second-order accurate scheme. Add solution in physical cells
    ! and in a single layer of the ghost cells along the momentum coordinate,
    ! plus two layers of ghost cells along the pitch-angle coordinates.
    ! This function is for the 3-dimensional VDF, with momentum (P^3/3),
    ! pitch angle (mu), and line coordinates (s_L).

    integer, intent(in) :: iLine     ! Indices of line and shock
    integer, intent(in) :: nX        ! Number of meshes along s_L axis
    real, intent(inout) :: VDF_G(-1:nP+2, -1:nMu+2, -1:nX+2)
    !--------------------------------------------------------------------------

    VDF_G(0:nP+1, 1:nMu, iStart:nX) = Distribution_CB(:, :, iStart:nX, iLine)
    VDF_G(0:nP+1, 1:nMu, 1) = Distribution_CB(:, :, iStart, iLine)

    ! Manipulate the LowerEndBc along the line coordinate:
    if(UseLowerEndBc) then
       call set_lower_end_bc(iLine)
       VDF_G(0:nP+1, 1:nMu, 0) = spread(max(LowerEndBc_I, &
            Background_I), DIM=2, NCOPIES=nMu)
    else
       call set_lower_end_vdf(iLine)
       VDF_G(0:nP+1, 1:nMu, 0) = spread(Background_I, DIM=2, NCOPIES=nMu)
    end if

    ! Manipulate the UpperEndBc along the line coordinate:
    if(UseUpperEndBc) then
       call set_upper_end_bc(iLine, nX)
       VDF_G(1:nP, 1:nMu, nX+1) = spread(max(UpperEndBc_I, &
            Background_I(1:nP)), DIM=2, NCOPIES=nMu)
       VDF_G(0   , 1:nMu, nX+1) = max(VDF_G(0, 1:nMu, nX), Background_I(0))
       VDF_G(nP+1, 1:nMu, nX+1) = max(VDF_G(nP+1,1:nMu,nX), Background_I(nP+1))
    else
       VDF_G(0:nP+1, 1:nMu, nX+1) = spread(Background_I, DIM=2, NCOPIES=nMu)
    end if

    ! Add a second layer of the ghost cells along the line coordinate:
    VDF_G(0:nP+1, 1:nMu,   -1) = VDF_G(0:nP+1, 1:nMu,    0)
    VDF_G(0:nP+1, 1:nMu, nX+2) = VDF_G(0:nP+1, 1:nMu, nX+1)
    ! Add two layer of the ghost cells along the line coordinate:
    VDF_G(0:nP+1,     0, :) = VDF_G(0:nP+1,    1,  :)
    VDF_G(0:nP+1,    -1, :) = VDF_G(0:nP+1,    0,  :)
    VDF_G(0:nP+1, nMu+1, :) = VDF_G(0:nP+1,   nMu, :)
    VDF_G(0:nP+1, nMu+2, :) = VDF_G(0:nP+1, nMu+1, :)
    ! Add a second layer of the ghost cells along the momentum coordinate:
    VDF_G(  -1, :, :) = VDF_G(   0, :, :)
    VDF_G(nP+2, :, :) = VDF_G(nP+1, :, :)

  end subroutine set_VDF3
  !============================================================================
end module SP_ModBc
!==============================================================================
