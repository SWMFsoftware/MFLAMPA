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
  public:: set_momentum_bc  ! Set the boundary condition at min/max energy
  public:: set_upper_end_bc ! Set the boundary condition at the far end point
  public:: set_lower_end_bc ! Set the boundary condition at the near-Sun end
  ! Boundary condition at the injection energy
  ! Injection efficiency and assumed spectral index with the energy
  ! range k_B*T_i < Energy < EnergyInjection, to be read from PARAM.in
  real, public     :: CoefInj = 0.25, SpectralIndex = 5.0
  real, parameter  :: CoefInjTiny = 2.5E-11 ! Before Nov 2023: set to be 0.0
  ! Type of lower end BC: float, escape, inject
  character(LEN=6) :: TypeLowerEndBc = 'inject'
  ! Upper end BC, set at the last point along the field line
  logical, public  :: UseUpperEndBc = .false.
  ! Type of upper end BC: none, float, escape, lism
  character(LEN=6) :: TypeUpperEndBc = 'none'
  real, public     :: UpperEndBc_I(1:nP), LowerEndBc_I(0:nP+1)
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
       call read_var('SpectralIndex', SpectralIndex)
       call read_var('Efficiency',    CoefInj)
    case('#ENDBC')
       ! Read the type of LowerEndBc and Select
       call read_var('TypeLowerEndBc', TypeLowerEndBc)

       ! Read whether to use UpperEndBc
       call read_var('UseUpperEndBc', UseUpperEndBc)
       if(UseUpperEndBc) then
          ! Read the type of UppenEndBc and Select
          call read_var('TypeUpperEndBc', TypeUpperEndBc)
          select case(trim(TypeUpperEndBc))
          case('none')
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
    ! set boundary conditions on grid point on the current line

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
       MomentumSi   = kinetic_energy_to_momentum(  &
            MhData_VIB(T_,iVertex,iLine)*Io2Si_V(UnitEnergy_))

       DistributionBc = (SpectralIndex-3)/(4*cPi)     &
            * MomentumInjSi**2*Io2Si_V(UnitEnergy_)   &
            * nSi_I(iVertex)/MomentumSi**3            &
            * (MomentumSi/MomentumInjSi)**SpectralIndex

       if(iShock /= NoShock_ .and. iVertex <= iShock + nWidth .and.  &
            iVertex >= iShock - nWidth) CoefInjLocal = CoefInj

       Distribution_CB(0, :, iVertex, iLine) = DistributionBc*CoefInjLocal
    end do

  end subroutine set_momentum_bc
  !============================================================================
  subroutine set_upper_end_bc(iLine, iEnd)
    ! set boundary condition at the last grid point on the given line
    ! assign the calculated BC to UpperEndBc_I

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
    ! set boundary condition at the first grid point on the given line
    ! assign the calculated BC to LowerEndBc_I

    integer, intent(in) :: iLine
    character(len=*), parameter:: NameSub = 'set_lower_end_bc'
    !--------------------------------------------------------------------------
    select case(trim(TypeLowerEndBc))
    case('float')
       LowerEndBc_I = Distribution_CB(0:nP+1, nMu, 1, iLine)
    case('escape')
       LowerEndBc_I = Background_I
    case('inject')
       LowerEndBc_I = Distribution_CB(0, 1, 1, iLine)  &
            /Momentum_I(0:nP+1)**SpectralIndex
    case default
       call CON_stop(NameSub//&
            ': Unknown type of lower end BC '//TypeLowerEndBc)
    end select
  end subroutine set_lower_end_bc
  !============================================================================
end module SP_ModBc
!==============================================================================
