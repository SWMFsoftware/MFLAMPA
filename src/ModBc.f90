!  Copyright (C) 2002 Regents of the University of Michigan,
!  portions used with permission
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module SP_ModBc

  ! The module sets up boundary conditions
  use ModNumConst, ONLY: cPi
  use SP_ModDistribution, ONLY: nP, Distribution_IIB, MomentumSi_I,   &
       MomentumInjSi
  use SP_ModGrid, ONLY: MHData_VIB, NoShock_, nWidth, T_
  use SP_ModUnit, ONLY: kinetic_energy_to_momentum, UnitEnergy_, Io2Si_V
  use ModUtilities, ONLY: CON_stop
  implicit none

  SAVE

  PRIVATE ! except

  ! Public members:
  public:: read_param
  public:: set_momentum_bc  ! Sets the boundary condition at min/max energy
  ! Boundary condition at the injection energy
  ! Injection efficiency and assumed spectral index with the energy
  ! range k_BT_i< Energy < EnergyInjection, to be read from PARAM.in
  real, public :: CoefInj = 0.25, SpectralIndex = 5.0
  real, parameter :: CoefInjTiny = 2.5E-11
contains
  !============================================================================
  subroutine read_param(NameCommand)
    use ModReadParam, ONLY: read_var
    character(len=*), intent(in):: NameCommand ! From PARAM.in
    character(len=8) :: StringScaleTurbulenceType
    character(len=*), parameter:: NameSub = 'read_param'
    !--------------------------------------------------------------------------
    select case(NameCommand)
    case('#INJECTION')
       call read_var('Efficiency',   CoefInj)
       call read_var('SpectralIndex',SpectralIndex)
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

    integer:: iVertex     ! loop variable
    real   :: MomentumSi    ! Momentum for the thermal energy k_BTi
    real   :: CoefInjLocal, DistributionBc
    !--------------------------------------------------------------------------
    do iVertex = 1, iEnd
       ! injection(Ti, Rho), see Sokolov et al., 2004, eq (3)
       ! f = CoefInj/2/pi * N / (2*m*T_p)^(3/2) * ((2*m*T_p)^(3/2)/p_inj)^5
       !   = CoefInj/2/pi * N / p^3 * (p/p_inj)^5
       ! where p = sqrt(2*m*T_p) is the momentum of thermal ion
       CoefInjLocal = CoefInjTiny
       MomentumSi   = kinetic_energy_to_momentum(                     &
            MHData_VIB(T_,iVertex,iLine)*Io2Si_V(UnitEnergy_))

       DistributionBc = (SpectralIndex-3)/(4*cPi)                     &
            * nSi_I(iVertex)/MomentumSi**3                            &
            * (MomentumSi/MomentumInjSi)**SpectralIndex

       if (iShock /= NoShock_ .and. iVertex <= iShock + nWidth .and.  &
            iVertex >= iShock - nWidth) CoefInjLocal = CoefInj

       Distribution_IIB(0,iVertex,iLine) = DistributionBc*CoefInjLocal
    end do

  end subroutine set_momentum_bc
  !============================================================================
end module SP_ModBc
!==============================================================================
