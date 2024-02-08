!  Copyright (C) 2002 Regents of the University of Michigan,
!  portions used
!  with permission
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module SP_ModTurbulence

  use ModConst
  use SP_ModDistribution, ONLY: nP
  use SP_ModGrid,         ONLY: iPTest, iParticleTest
  implicit none
  SAVE

  private

  public :: init, finalize, UseTurbulentSpectrum, set_dxx, &
       read_param, dxx, DoTraceShock, advance_log_advection

  logical:: UseTurbulentSpectrum        = .false.
  logical:: DoTraceShock                = .true.

  integer, parameter :: nK = nP
  real    :: dLogK
  real    :: DtReduced

  real, allocatable :: Gamma_I(:,:)
  real, allocatable :: IPlusSi_IX(:,:),IMinusSi_IX(:,:),ICSi_X(:)
  real, allocatable :: vAlfvenSi_I(:)
  real, allocatable :: kOverBSi_I(:)
  real, allocatable :: kSi_I(:)

  ! Rate of advection in k space, neglecting
  ! a spacial advection with the Alfven speed
  real, allocatable:: DispersionA_I(:)

  ! Rate of advection in k space, for I_+/I_- wave
  real, allocatable:: DispersionPlus_I(:)
  real, allocatable:: DispersionMinus_I(:)

  ! This is the ratio of densities powered 3/2
  real, allocatable:: RhoCompression_I(:)

  !------------------------------------------------------------------------!
  !          Grid in the momentum space                                    !
  ! iP     0     1                         nP   nP+1                       !
  !       |     |    ....                 |     |                          !
  ! P      P_inj P_inj*exp(\Delta (Ln P))  P_Max P_Max*exp(\Delta (Ln P))  !
  !             |    Grid in k-space      |     |                          !
  ! K/B         KMax                      KMin                             !
  ! ik     0     1                         nP   nP+1                       !
  !------------------------------------------------------------------------!

  real,allocatable,private:: AK_II(:,:)
  real,allocatable,private:: BK_II(:,:)
  real,allocatable,private:: CFL_I(:)
  integer,allocatable::      CorrectionMode_X(:)

  ! the intensity of the back travelling wave in the initial condition
  real :: Alpha       = 1.0/10
  real :: Lambda0InAu = 4.0/10  ! [AU]
contains
  !============================================================================
  subroutine read_param(NameCommand)
    use ModReadParam, ONLY: read_var
    character(len=*), intent(in):: NameCommand
    character(len=*), parameter:: NameSub = 'read_param'
    !--------------------------------------------------------------------------
    select case(NameCommand)
    case('#TURBULENTSPECTRUM')
       call read_var('UseTurbulentSpectrum', UseTurbulentSpectrum)
    case('#TRACESHOCK')
       call read_var('DoTraceShock', DoTraceShock)
    case default
       call CON_stop(NameSub//' Unknown command '//NameCommand)
    end select
  end subroutine read_param
  !============================================================================
  subroutine init
    ! Init all the allocatable vars
    use SP_ModSize, ONLY: nVertexMax
    use SP_ModProc, ONLY: iProc
    character(len=*), parameter:: NameSub = 'init'
    !--------------------------------------------------------------------------

    allocate(IPlusSi_IX(0:nP+1,1:nVertexMax), &
         IMinusSi_IX(0:nP+1,1:nVertexMax))
    allocate(kOverBSi_I(0:nP+1), kSi_I(0:nP+1))
    allocate(ICSi_X(1:nVertexMax),CorrectionMode_X(1:nVertexMax))
    allocate(vAlfvenSi_I(1:nVertexMax))
    allocate(DispersionA_I(1:nVertexMax))

    allocate(RhoCompression_I(1:nVertexMax))

    ! if(DoOutputGamma) &
    !     allocate(Gamma_I(nP,(iXOutputLast-iXOutputStart)/iXOutputStride+1))

    allocate(DispersionPlus_I(1:nVertexMax))
    allocate(DispersionMinus_I(1:nVertexMax))
    allocate(CFL_I(1:nVertexMax))

    allocate(AK_II(nP,nVertexMax),BK_II(nP,nVertexMax))
  end subroutine init
  !============================================================================
  subroutine finalize

    character(len=*), parameter:: NameSub = 'finalize'
    !--------------------------------------------------------------------------
    deallocate(IPlusSi_IX,IMinusSi_IX)
    deallocate(kOverBSi_I, kSi_I)
    deallocate(ICSi_X,CorrectionMode_X)
    deallocate(vAlfvenSi_I)
    deallocate(DispersionA_I)

    deallocate(RhoCompression_I)

    ! if(DoOutputGamma) &
    !     allocate(Gamma_I(nP,(iXOutputLast-iXOutputStart)/iXOutputStride+1))

    deallocate(DispersionPlus_I)
    deallocate(DispersionMinus_I)
    deallocate(CFL_I)

    deallocate(AK_II,BK_II)

  end subroutine finalize
  !============================================================================
  subroutine set_dxx(iEnd, nP, BSi_I)
    integer,intent(in) :: iEnd,nP
    real,intent(in)    :: BSi_I(iEnd)

    integer:: iVertex,iK
    real :: F01,F02,F11,F12
    real :: k0Si,k1Si
    real :: SpectralIndexAtKMax, ISumSi

    logical :: DoTestMe = .false.

    ! The coefficient of a spatial diffusion along the magnetic field is
    ! given by a formula in the SI unit:
    ! D_{xx}=v*B^2/(cPi*cMu)*(AK_I(k_{res})-k_{res}^2*BK_I(k_{res})),

    ! where AK_I=\int_{k_{res}}^\infty{d(\log k)/(k^2*(I_{+}(k) + I_{-}(k) )) }
    ! and   BK_I=\int_{k_{res}}^\infty{d(\log k)/(k^4*(I_{+}(k) + I_{-}(k) )) }

    !          Grid in the momentum space                                     !
    ! iP     0     1                         nP   nP+1                         !
    !       |     |    ....                 |     |                           !
    ! P      P_inj P_inj*exp(\Delta (Ln P))  P_Max P_Max*exp(\Delta (Ln P))    !
    !             |    Grid in k-space      |     |                           !
    ! K/B         kMax(k0)                  kMin                               !
    ! ik     0     1                         nP   nP+1                         !
    !--------------------------------------------------------------------------

    do iVertex=1,iEnd
       select case(CorrectionMode_X(iVertex))
       case(1)
          SpectralIndexAtKMax = 5.0/3
       case(2)
          SpectralIndexAtKMax = 1.0
       end select

       kSi_I = BSi_I(iVertex)*kOverBSi_I

       ! Initially from kMax
       k0Si = kSi_I(1)

       ! The sum of I_{plus}+I_{minus} at P_max
       ISumSi = IPlusSi_IX(1,iVertex)+IMinusSi_IX(1,iVertex)

       ! The integrand for AK_I, BK_I
       F01 = 1.0/(k0Si**2)/ISumSi
       F02 = 1.0/(k0Si**4)/ISumSi

       ! As the starting values for AK_I and BK_I at the minimum momentum,
       ! solve the integrals from K_{max} up to \infty, assuming the
       ! power law spectrum of turbulence at K>K_{max}
       AK_II(1,iVertex)=F01/(2.0-SpectralIndexAtKMax)
       BK_II(1,iVertex)=F02/(4.0-SpectralIndexAtKMax)

       do iK=2,nP
          ! We calculate the partial sums for a set of the wave number values.
          ! The integral is taken from KRes up to infinity, so we start from
          ! the maximal wave number and add the contributions from each of
          ! the wave number intervals.

          ! The current value of the wave number
          k1Si = kSi_I(iK)

          ! The sum of I_{plus}+I_{minus} at P
          ISumSi = IPlusSi_IX(iK,iVertex)+IMinusSi_IX(iK,iVertex)

          ! The integrands at the lower value of the wave number
          F11 = 1.0/(k1Si**2)/ ISumSi
          F12 = 1.0/(k1Si**4)/ ISumSi

          ! Calculate the new partial sums

          AK_II(iK,iVertex)=AK_II(iK-1,iVertex)+0.5*(F01+F11)*dLogK
          BK_II(iK,iVertex)=BK_II(iK-1,iVertex)+0.5*(F02+F12)*dLogK

          ! current values saved as the initial values for the next step in
          ! the loop
          k0Si = k1Si; F01=F11; F02=F12
       end do
    end do
  end subroutine set_dxx
  !============================================================================
  real function Dxx(iX, iP, MomentumSi, SpeedSi, BSi)
    integer,intent(in) :: iX, iP
    real,   intent(in) :: MomentumSi, SpeedSi, BSi

    real    :: kRSi

    logical :: DoTestMe =.false.

    !          Grid in the momentum space                                     !
    ! iP     0     1                         nP   nP+1                         !
    !       |     |    ....                 |     |                           !
    ! P      P_inj P_inj*exp(\Delta (Ln P))  P_Max P_Max*exp(\Delta (Ln P))    !
    !             |    Grid in k-space      |     |                           !
    ! K/B         KMax                      KMin                               !
    ! ik     0     1                         nP   nP+1                         !
    !--------------------------------------------------------------------------

    ! The coefficient of a spatial diffusion along the magnetic field is
    ! given by a formula in the SI unit:
    ! D_{xx}=v*B^2/(cPi*cMu)*(AK_I(k_{res})-k_{res}^2*BK_I(k_{res})),

    ! where AK_I=\int_{k_{res}}^\infty{d(\log k)/(k^2*(I_{+}(k) + I_{-}(k) )) }
    ! and   BK_I=\int_{k_{res}}^\infty{d(\log k)/(k^4*(I_{+}(k) + I_{-}(k) )) }

    ! The resonant wave number, kr = e*B/p in the SI unit
    kRSi = cElectronCharge*BSi/MomentumSi

    ! Calculate D_{xx}: KRes-dependent part
    Dxx = BSi**2*SpeedSi/(cMu*cPi) * (AK_II(iP,iX)-BK_II(iP,iX)*kRSi**2)

    if (iX == iParticleTest .and. iP == iPTest .and. DoTestMe) then
       write(*,*) 'AK_II(iP,iX), BK_II(iP,iX) =', &
            AK_II(iP,iX), BK_II(iP,iX)
       write(*,*) 'Dxx =', Dxx
       write(*,*) 'D from Li =', 1./3*Lambda0InAu*9.3286125374064124E+08 &
            *(MomentumSi*cLightSpeed/cGeV)**(1./3)*SpeedSi
    end if

  end function Dxx
  !============================================================================
  subroutine advance_log_advection(CFLIn, nGCLeft, nGCRight,        &
       FInOut_I, IsConservative, DeltaLnP)
    ! The procedure integrates the log-advection equation, in
    ! the conservative or non-conservative formulation, at a logarithmic grid,
    ! using a single-stage second order scheme

    real,   intent(in):: CFLIn        ! Time step * acceleration rate/(Dlnp)
    integer,intent(in):: nGCLeft      ! The solution in the ghost cells is not
    integer,intent(in):: nGCRight     ! advanced in time but used as the BCs
    real,intent(inout):: FInOut_I(1-nGCLeft:nP+nGCRight)  ! Solution
    ! sol. index 0 is to set boundary condition at the injection energy
    logical,intent(in):: IsConservative  ! Solve (C) if .true. (NC) otherwise
    real,optional,intent(in)::DeltaLnP   ! Used only for NonConservative=.true.!
    ! Subcycling to achieve the condition CFL<1, if nStep>1
    integer :: nStep
    ! Loop variables
    integer :: iStep
    ! Extended version of the sulution array to implement BCc
    real    :: F_I(1 - max(nGCLeft,2):nP+max(nGCRight,2))
    real    :: FSemiintUp_I(0:nP+1), FSemiintDown_I(0:nP+1)
    real    :: CFL, HalfADtIfNeeded
    !-------------------------NonConservative---------------------------
    ! This is a one-stage second-order scheme for the advection equation:
    !             f_t+A f_{ln p}=0,
    !
    ! Herewith CFL = A*Delta t/Delta(ln p):
    ! f^(n+1)_i-CFL*(f^n_(i-1/2)-f^n_(i+1/2)=f^n_i, where
    ! For CFL>0:
    ! f^n_(i-1/2)=f^n_(i-1)+cHalf*(cOne-CFL)*df_lim^n_(i-1)
    ! For CFL<0:
    ! f^n_(i-1/2)=f^n_(i  )-cHalf*(cOne+CFL)*df_lim^n_(i  )
    !-------------------------   Conservative---------------------------
    ! This is a one-stage second-order scheme for the advection equation:
    !             f_t+A (f p}_p=0,
    ! We now need to use both
    ! A*Delta t = CFLIn*Delta (LnP)
    ! and
    ! CFL = A*Delta t/(2 tanh(Delta (ln p)/2)
    ! f^(n+1)_i = CFL*(f^n_(i-1/2)-f^n_(i+1/2))-
    !        (1/2)(f^n_(i-1/2)+f^n_(i+1/2))A*Delta t + f^n_i,
    ! where
    ! For CFL>0:
    ! f^n_(i-1/2)=f^n_(i-1)*(1-(1/2)*A*Delta t)+cHalf*(cOne-CFL)*df_lim^n_(i-1)
    ! For CFL<0:
    ! f^n_(i-1/2)=f^n_(i  )*(1-(1/2)*A*Delta t)-cHalf*(cOne+CFL)*df_lim^n_(i  )

    !--------------------------------------------------------------------------
    if(IsConservative)then
       HalfADtIfNeeded = 0.50*CFLIn*DeltaLnP
       CFL = HalfADtIfneeded/tanh(0.50 * DeltaLnP)
       HalfADtIfNeeded = -HalfADtIfNeeded
    else
       HalfADtIfNeeded = 0.0; CFL = CFLIn
    end if

    nStep = 1 + int(abs(CFL)); CFL = CFL/real(nStep)
    HalfADtIfNeeded=HalfADtIfNeeded/real(nStep)
    F_I(1-nGCLeft:nP+nGCRight) = FInOut_I(1 - nGCLeft:nP+nGCRight)

    ! Check for positivity
    if(any(F_I(1-nGCLeft:nP+nGCRight)<=0.0))then
       write(*,*)'Before advection F_I <=0'
       write(*,*)IsConservative
       write(*,*)F_I
       stop
    end if

    ! One stage second order upwind scheme

    if (CFL>0.0) then
       do iStep = 1, nStep
          ! Boundary condition at the left boundary
          if(nGCLeft<2)F_I(            -1:0-nGCLeft) = F_I( 1-nGCLeft )
          ! Boundary condition at the right boundary
          if(nGCRight<2)F_I(nP+1-nGCRight:nP+2     ) = F_I(nP+nGCRight)

          ! f_(i+1/2):
          FSemiintUp_I(0:nP) = F_I(0:nP)*(1.0 + HalfADtIfNeeded)&
               + 0.5*(1.0 - CFL)*df_lim_array(0, nP)
          ! f_(i-1/2):
          FSemiintDown_I(1:nP) = FSemiintUp_I(0:nP-1)

          ! Update the solution from f^(n) to f^(n+1):
          F_I(1:nP) = F_I(1:nP)+CFL*(FSemiintDown_I(1:nP)-FSemiintUp_I(1:nP))+&
               HalfADtIfNeeded*(FSemiintDown_I(1:nP)+FSemiintUp_I(1:nP))
       end do
    else
       do iStep = 1, nStep
          ! Boundary condition at the left boundary
          if(nGCLeft<2)F_I(            -1:0-nGCLeft) = F_I( 1-nGCLeft)
          ! Boundary condition at the right boundary
          if(nGCRight<2)F_I(nP+1-nGCRight:nP+2     ) = F_I(nP+nGCRight)

          ! f_(i-1/2):
          FSemiintDown_I(1:nP+1) = F_I(1:nP+1)*(1.0 + HalfADtIfNeeded)&
               - 0.5*(1.0 + CFL)*df_lim_array(1, nP+1)
          ! f_(i+1/2):
          FSemiintUp_I(1:nP) = FSemiintDown_I(2:nP+1)

          ! Update the solution from f^(n) to f^(n+1):
          F_I(1:nP) = F_I(1:nP)+CFL*(FSemiintDown_I(1:nP)-FSemiintUp_I(1:nP))+&
               HalfADtIfNeeded*(FSemiintDown_I(1:nP)+FSemiintUp_I(1:nP))
       end do
    end if

    FInOut_I(1:nP) = F_I(1:nP)
    if(any(FInOut_I(1-nGCLeft:nP+nGCRight)<=0.0))then
       write(*,*)'After advection F_I <=0, for CFLFermi= ',CFL
       write(*,*)F_I
       call CON_stop('Negative distribution function')
    end if
  contains
    !==========================================================================
    function df_lim_array(iLeft, iRight) result(df_lim_arr)
      integer, intent(in):: iLeft, iRight ! start/left and end/right indices
      real :: df_lim_arr(iLeft:iRight)    ! output results
      real :: dF1(iLeft:iRight), dF2(iLeft:iRight)

      !------------------------------------------------------------------------
      dF1 = F_I(iLeft+1:iRight+1)-F_I(iLeft:iRight)
      dF2 = F_I(iLeft:iRight)-F_I(iLeft-1:iRight-1)

      ! df_lim=0 if dF1*dF2<0, sign(dF1) otherwise:
      df_lim_arr = sign(0.50,dF1)+sign(0.50,dF2)
      dF1 = abs(dF1)
      dF2 = abs(dF2)
      df_lim_arr = df_lim_arr*min(max(dF1,dF2),2.0*dF1,2.0*dF2)
    end function df_lim_array
    !==========================================================================
  end subroutine advance_log_advection
  !============================================================================
end module SP_ModTurbulence
!==============================================================================
