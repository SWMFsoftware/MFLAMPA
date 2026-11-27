!  Copyright (C) 2002 Regents of the University of Michigan,
!  portions used with permission
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!===========================TESTS==============================================
! Unit and functional tests of the focused transport equation in M-FLAMPA.
!
! Purpose: each test turns on one term (or a few), runs the numerical
! solver on a simple setup with an analytic solution, and checks the
! accuracy against the exact values. These benchmarks show that each
! term does the right physics and that the solvers are implemented right.
!
! Framework: one module per term (or group of terms) with its own tests,
! and one program at the end that runs the tests by name:
!   ModTestDperp    perpendicular diffusion across the field lines
!   ModTestFocused  focused transport along the field lines (later)
! Each test is one subroutine that contains its own helper routines.
! Usage: ./test_transport.exe [NameTest]; without NameTest all tests run.
! The program runs on 1 or more MPI processors, e.g.
!   mpirun -np 4 ./test_transport.exe dperp
! and the results do not depend on the number of processors.
! Each test prints its error norms and writes test_<name>.out.
!==============================================================================
module ModTestDperp

  ! Tests of the perpendicular diffusion solver (SP_ModPerpDiffusion in
  ! src/ModDiffusionPerp.f90) against analytic solutions:
  ! - test_dperpUniSphere3d: a point source with a uniform Dperp inside a
  !   sphere spreads as the 3D Gaussian of free space, and the number of
  !   particles is conserved
  ! Later, e.g.: a Y_lm on a sphere of radius r decays as
  ! exp(-l(l+1)*Dperp*t/r^2).
  ! The tests call the solver of SP_ModPerpDiffusion on its own mesh; they
  ! do not have a copy of it.

  use ModMpi
  use ModNumConst, ONLY: cPi, cRadToDeg
  use ModUtilities, ONLY: CON_stop, open_file, close_file
  use ModIoUnit, ONLY: UnitTmp_
  use ModPlotFile, ONLY: save_plot_file
  use ModReadParam, ONLY: read_file, read_init, read_line, &
       read_command, n_line_read, lStringLine
  use SP_ModProc, ONLY: iComm, iProc, nProc, iError
  use SP_ModGrid, ONLY: nP, nLine
  use SP_ModUnit, ONLY: init_unit => init, Io2Si_V, UnitX_
  use SP_ModPerpDiffusion, ONLY: read_param_dperp => read_param, &
       init_dperp => init, solver_fvm_diffuse3d, nRPerp, nThetaPerp, &
       nPhiPerp, RPerp_C, ThetaPerp_C, PhiPerp_C, XyzPerp_CB, Volume_CB

  implicit none

  private ! except

  public:: test_dperpUniSphere3d ! point source, uniform Dperp in a sphere

  ! Setup of test_dperpUniSphere3d. Units: lengths in Rs (the unit of the
  ! mesh in the #DIFFUSIONPERP command), time in s and Dperp in Rs^2/s.
  ! The solver works in SI units (m, s, m^2/s): the test converts with
  ! Io2Si_V(UnitX_).
  !
  ! Final time and Dperp: the Gaussian has the width
  ! sqrt(4*DperpTest*tFinalTest) = 2 Rs at tFinalTest = 5 s
  real, parameter:: tFinalTest = 5.0, DperpTest = 0.2
  ! Base mesh of the solver: nRTest shells between RMinTest and RMaxTest
  ! (Linear) and nThetaTest x nPhiTest cells on each shell. The source
  ! cell (iR = nR/2+1) is at r = 11.5 Rs, where the cells are about 1 Rs
  ! (half of the width) wide in all 3 directions. The test is repeated
  ! on the mesh refined by 2 in all 3 directions.
  real, parameter:: RMinTest = 2.0, RMaxTest = 20.0
  integer, parameter:: nRTest = 18, nThetaTest = 36, nPhiTest = 72
  ! Dperp = DperpTest in the cells closer to the source than RadiusDperp
  ! times the width, where the exact f is above exp(-RadiusDperp**2) =
  ! 2.6e-5 of its peak; Dperp = 0 in the other cells. The solver takes its
  ! stable time step from the cells with Dperp > 0 only, so the tiny cells
  ! near the poles and RMin do not slow down the test: with Dperp = const
  ! everywhere the base mesh would need 22129 sub-steps instead of 21.
  real, parameter:: RadiusDperp = 3.25

contains
  !============================================================================
  subroutine test_dperpUniSphere3d

    ! Point source test of the 3D finite volume Dperp solver
    ! solver_fvm_diffuse3d on its own r-theta-phi mesh (from init).
    !
    ! Setup: f = 1 in the middle cell (iR, iTheta, iPhi) =
    ! (nR/2+1, nTheta/2+1, nPhi/2+1) and f = 0 in all other cells, the
    ! same at all momenta. Dperp is constant around the source (see
    ! RadiusDperp). One call of the solver advances f from t = 0 to
    ! tFinalTest with its own stable sub-steps.
    !
    ! Analytic solution: far from the walls and the poles, f spreads as
    ! a point source in free space,
    !   f(x,t) = M/(4*pi*Dperp*t)**1.5 * exp(-|x - x0|**2/(4*Dperp*t)),
    ! where x0 is the center of the source cell and M = its volume (the
    ! number of particles, since f = 1 there). The free-space solution
    ! holds: Dperp > 0 only within RadiusDperp widths of x0, and the walls
    ! and the poles are farther away. Representing the source cell as a
    ! point is a part of the discretization error, which is of 2nd order.
    ! The errors include the error of the explicit sub-steps, which partly
    ! cancels the spatial error (smaller sub-steps give larger errors).
    !
    ! Outputs (on processor 0): for the base mesh and the 2x refined mesh
    ! the number of particles sum(f*V) at t = 0 and tFinalTest, min f,
    ! two checks of the solver output (the same f at all momenta, and its
    ! df/dt times t = the change of f; both 0 up to round-off) and the
    ! relative max, L1 and L2 errors, then a table of the errors and the
    ! order of convergence (2 is expected). Files (f, exact f and error):
    !   test_dperpUniSphere3d.out     phi-theta map on the source shell,
    !   test_dperpUniSphere3d_r.out   radial cut through the source,
    !                                 both on the base mesh
    !   test_dperpUniSphere3d_x2*.out the same for the refined mesh
    ! Any number of MPI processors: the solver updates the shells of each
    ! processor; the sums are made per shell and the plots are collected
    ! on processor 0, so the results are the same for any nProc.

    integer, parameter:: nMesh = 2
    ! Relative max, L1 and L2 errors on the base and the refined meshes
    real:: Error_VI(3,nMesh)
    integer:: iMesh
    character(len=*), parameter:: NameSub = 'test_dperpUniSphere3d'
    !--------------------------------------------------------------------------
    ! Unit conversion: Io2Si_V(UnitX_) = Rs [m]
    call init_unit
    ! No field lines in this test (init_dperp allocates DPerp_CB with nLine)
    nLine = 0

    do iMesh = 1, nMesh
       call run_dperpUniSphere3d(2**(iMesh-1), Error_VI(:,iMesh))
    end do

    if(iProc /= 0) RETURN
    write(*,'(/,a)') NameSub//': relative errors at t = tFinalTest'
    write(*,'(a)') '   nR nTheta  nPhi   max error    L1 error    L2 error'
    do iMesh = 1, nMesh
       write(*,'(i5,i7,i6,3es12.4)') nRTest*2**(iMesh-1), &
            nThetaTest*2**(iMesh-1), nPhiTest*2**(iMesh-1), Error_VI(:,iMesh)
    end do
    write(*,'(a,3f12.3)') ' order of convergence', &
         log(Error_VI(:,1)/Error_VI(:,2))/log(2.0)

  contains
    !==========================================================================
    subroutine run_dperpUniSphere3d(nRefine, Error_V)

      ! One mesh of the test: the base mesh refined by nRefine in all 3
      ! directions. Set the mesh, f and Dperp, call the solver to
      ! tFinalTest, compare with the exact solution, print the results and
      ! save the plot files

      integer, intent(in):: nRefine
      real,    intent(out):: Error_V(3) ! relative max, L1 and L2 errors

      ! Shells of this processor: the same split as in init_dperp
      integer:: iRStart, iREnd, nRChunk, nRExtra
      ! Indexes of the source cell
      integer:: iR0, iTheta0, iPhi0
      ! Dperp [m^2/s], 4*Dperp*t [m^2], center of the source cell [m], and
      ! the number of particles in it (f = 1 times its volume [m^3])
      real:: DperpSi, Width2, Xyz0_D(3), Mass0
      ! f on the shells of this processor at all momenta (with the ghost
      ! momenta 0 and nP+1), df/dt from the solver (averaged over its
      ! sub-steps), and Dperp on all shells (input of the solver)
      real, allocatable:: Distr_4D(:,:,:,:), DfDt_4D(:,:,:,:), &
           Dperp_4D(:,:,:,:)
      ! Sums over the cells of each shell, zero for the shells of the other
      ! processors (so the sums do not depend on nProc): sum(f*V) at t = 0
      ! and t = tFinalTest, sum(|f - fExact|*V), sum(fExact*V),
      ! sum((f - fExact)**2*V) and sum(fExact**2*V)
      integer, parameter:: nSum = 6
      real, allocatable:: Sum_VI(:,:)
      real:: Sum_V(nSum)
      ! Max |f - fExact|, max fExact, max |f - f(p=1)| over the momenta,
      ! max |DfDt*t - (f - f(t=0))| (DfDt*t is the change of f), and min f
      ! over all cells and momenta
      real:: Max_I(4), DistrMin
      ! Plot data: f, fExact and f - fExact on the source shell and along
      ! the radial line through the source
      real, allocatable:: Plot_VII(:,:,:), Plot_VI(:,:)
      ! Exact, numerical and initial f and volume of a cell, Rs [m],
      ! wall time [s]
      real:: Exact, Distr, Distr0, Volume, Rs, CpuTime
      integer:: i, j, k
      character(len=40):: NameFile
      character(len=*), parameter:: NameSub = 'run_dperpUniSphere3d'
      !------------------------------------------------------------------------
      call set_mesh_dperp(nRTest*nRefine, nThetaTest*nRefine, &
           nPhiTest*nRefine)

      ! Shells of this processor: the first mod(nRPerp,nProc) processors
      ! have one more (as iRPerpStart_I and iRPerpEnd_I in init_dperp)
      nRChunk = nRPerp/nProc
      nRExtra = mod(nRPerp, nProc)
      iRStart = iProc*nRChunk + min(iProc, nRExtra) + 1
      iREnd   = iRStart + nRChunk - 1
      if(iProc < nRExtra) iREnd = iREnd + 1

      ! The source cell (middle of the mesh) and the exact solution
      iR0 = nRPerp/2 + 1; iTheta0 = nThetaPerp/2 + 1; iPhi0 = nPhiPerp/2 + 1
      Xyz0_D = XyzPerp_CB(:,iPhi0,iTheta0,iR0)
      Mass0  = Volume_CB(iPhi0,iTheta0,iR0)
      Rs      = Io2Si_V(UnitX_)
      DperpSi = DperpTest*Rs**2
      Width2  = 4*DperpSi*tFinalTest

      ! Dperp on all shells (the same at all momenta)
      allocate(Dperp_4D(nP,nPhiPerp,nThetaPerp,nRPerp))
      do i = 1, nRPerp; do j = 1, nThetaPerp; do k = 1, nPhiPerp
         if(sum((XyzPerp_CB(:,k,j,i) - Xyz0_D)**2) &
              <= RadiusDperp**2*Width2)then
            Dperp_4D(:,k,j,i) = DperpSi
         else
            Dperp_4D(:,k,j,i) = 0.0
         end if
      end do; end do; end do

      ! Initial f on the shells of this processor (all momenta)
      allocate(Distr_4D(0:nP+1,nPhiPerp,nThetaPerp,iRStart:iREnd), &
           DfDt_4D(0:nP+1,nPhiPerp,nThetaPerp,iRStart:iREnd))
      Distr_4D = 0.0
      if(iR0 >= iRStart .and. iR0 <= iREnd) &
           Distr_4D(:,iPhi0,iTheta0,iR0) = 1.0

      ! Number of particles at t = 0
      allocate(Sum_VI(nSum,nRPerp))
      Sum_VI = 0.0
      do i = iRStart, iREnd
         Sum_VI(1,i) = sum(Distr_4D(1,:,:,i)*Volume_CB(:,:,i))
      end do

      ! Advance f to tFinalTest with the solver (time accurate mode)
      call MPI_barrier(iComm, iError)
      CpuTime = MPI_wtime()
      call solver_fvm_diffuse3d(Dperp_4D, iRStart, iREnd, .false., &
           tFinalTest, Distr_4D, DfDt_4D)
      CpuTime = MPI_wtime() - CpuTime

      ! Compare with the exact solution on the shells of this processor.
      ! All momenta have the same f; check it and use p = 1. Check that
      ! the df/dt returned by the solver (averaged over its sub-steps)
      ! times tFinalTest is the change of f.
      Max_I = 0.0
      DistrMin = huge(1.0)
      do i = iRStart, iREnd; do j = 1, nThetaPerp; do k = 1, nPhiPerp
         Exact  = exact_distr(XyzPerp_CB(:,k,j,i), Xyz0_D, Mass0, Width2)
         Distr  = Distr_4D(1,k,j,i)
         Volume = Volume_CB(k,j,i)
         Distr0 = 0.0
         if(i == iR0 .and. j == iTheta0 .and. k == iPhi0) Distr0 = 1.0
         Sum_VI(2:nSum,i) = Sum_VI(2:nSum,i) + Volume*[Distr, &
              abs(Distr - Exact), Exact, (Distr - Exact)**2, Exact**2]
         Max_I(1) = max(Max_I(1), abs(Distr - Exact))
         Max_I(2) = max(Max_I(2), Exact)
         Max_I(3) = max(Max_I(3), maxval(abs(Distr_4D(1:nP,k,j,i) - Distr)))
         Max_I(4) = max(Max_I(4), maxval(abs(DfDt_4D(1:nP,k,j,i)*tFinalTest &
              - (Distr_4D(1:nP,k,j,i) - Distr0))))
         DistrMin = min(DistrMin, minval(Distr_4D(1:nP,k,j,i)))
      end do; end do; end do
      ! Collect the sums of all shells and add them up in the same order on
      ! every processor; collect the extrema
      call MPI_allreduce(MPI_IN_PLACE, Sum_VI, nSum*nRPerp, MPI_REAL, &
           MPI_SUM, iComm, iError)
      Sum_V = 0.0
      do i = 1, nRPerp
         Sum_V = Sum_V + Sum_VI(:,i)
      end do
      call MPI_allreduce(MPI_IN_PLACE, Max_I, 4, MPI_REAL, MPI_MAX, &
           iComm, iError)
      call MPI_allreduce(MPI_IN_PLACE, DistrMin, 1, MPI_REAL, MPI_MIN, &
           iComm, iError)
      ! Relative errors: max |f - fExact|/max fExact, the L1 error over the
      ! L1 norm of fExact and the same for L2 (sums weighted by the volume)
      Error_V(1) = Max_I(1)/Max_I(2)
      Error_V(2) = Sum_V(3)/Sum_V(4)
      Error_V(3) = sqrt(Sum_V(5)/Sum_V(6))

      if(iProc == 0)then
         write(*,'(/,a,i0,a,i0,a,i0,a,i0)') NameSub// &
              ': mesh nR x nTheta x nPhi = ', nRPerp, ' x ', nThetaPerp, &
              ' x ', nPhiPerp, ', nProc = ', nProc
         write(*,'(a,f7.3,a,2f8.3,a)') ' source cell at r =', &
              RPerp_C(iR0)/Rs, ' Rs, (theta, phi) =', &
              ThetaPerp_C(iTheta0)*cRadToDeg, PhiPerp_C(iPhi0)*cRadToDeg, &
              ' deg'
         write(*,'(a,f6.3,a,f6.2,a,f6.3,a)') ' Dperp =', DperpTest, &
              ' Rs^2/s, t =', tFinalTest, ' s, width sqrt(4*Dperp*t) =', &
              sqrt(Width2)/Rs, ' Rs'
         write(*,'(a,f9.3,a)') ' wall time of the solver   =', CpuTime, ' s'
         write(*,'(a,2es22.14)') ' sum(f*V) [Rs^3] at t=0, t =', &
              Sum_V(1:2)/Rs**3
         write(*,'(a,es11.3)') ' relative change of sum(f*V)   =', &
              (Sum_V(2) - Sum_V(1))/Sum_V(1)
         write(*,'(a,es11.3)') ' min f                         =', DistrMin
         write(*,'(a,es11.3)') ' max |f(p) - f(p=1)|, p=1..nP  =', Max_I(3)
         write(*,'(a,es11.3)') ' max |DfDt*t - (f - f(t=0))|   =', Max_I(4)
         write(*,'(a,3es11.3)') ' relative errors: max, L1, L2  =', Error_V
      end if

      ! Plots of f, exact f and error: phi-theta map on the source shell
      ! (from the processor of that shell) and radial cut through the source
      ! (from all processors), collected (summed with zeros) on all
      ! processors and written by processor 0
      allocate(Plot_VII(3,nPhiPerp,nThetaPerp), Plot_VI(3,nRPerp))
      Plot_VII = 0.0; Plot_VI = 0.0
      if(iR0 >= iRStart .and. iR0 <= iREnd)then
         do j = 1, nThetaPerp; do k = 1, nPhiPerp
            Distr = Distr_4D(1,k,j,iR0)
            Exact = exact_distr(XyzPerp_CB(:,k,j,iR0), Xyz0_D, Mass0, Width2)
            Plot_VII(:,k,j) = [Distr, Exact, Distr - Exact]
         end do; end do
      end if
      do i = iRStart, iREnd
         Distr = Distr_4D(1,iPhi0,iTheta0,i)
         Exact = exact_distr(XyzPerp_CB(:,iPhi0,iTheta0,i), &
              Xyz0_D, Mass0, Width2)
         Plot_VI(:,i) = [Distr, Exact, Distr - Exact]
      end do
      call MPI_allreduce(MPI_IN_PLACE, Plot_VII, size(Plot_VII), MPI_REAL, &
           MPI_SUM, iComm, iError)
      call MPI_allreduce(MPI_IN_PLACE, Plot_VI, size(Plot_VI), MPI_REAL, &
           MPI_SUM, iComm, iError)
      if(iProc == 0)then
         NameFile = 'test_dperpUniSphere3d'
         if(nRefine > 1) write(NameFile,'(a,i0)') &
              'test_dperpUniSphere3d_x', nRefine
         call save_plot_file(NameFile=trim(NameFile)//'.out', &
              TypeFileIn='ascii', TimeIn=tFinalTest, nStepIn=1, &
              StringHeaderIn='test_dperpUniSphere3d: f, exact f and '// &
              'f - exact f on the source shell, angles in deg', &
              NameVarIn='Phi Theta f fExact Error rShell Dperp', &
              ParamIn_I=[RPerp_C(iR0)/Rs, DperpTest], &
              Coord1In_I=PhiPerp_C*cRadToDeg, &
              Coord2In_I=ThetaPerp_C*cRadToDeg, VarIn_VII=Plot_VII)
         call save_plot_file(NameFile=trim(NameFile)//'_r.out', &
              TypeFileIn='ascii', TimeIn=tFinalTest, nStepIn=1, &
              StringHeaderIn='test_dperpUniSphere3d: f, exact f and '// &
              'f - exact f along r through the source, r in Rs', &
              NameVarIn='r f fExact Error Theta Phi Dperp', &
              ParamIn_I=[ThetaPerp_C(iTheta0)*cRadToDeg, &
              PhiPerp_C(iPhi0)*cRadToDeg, DperpTest], &
              Coord1In_I=RPerp_C/Rs, VarIn_VI=Plot_VI)
      end if

      deallocate(Dperp_4D, Distr_4D, DfDt_4D, Sum_VI, Plot_VII, Plot_VI)

    end subroutine run_dperpUniSphere3d
    !==========================================================================
    real function exact_distr(Xyz_D, Xyz0_D, Mass0, Width2)

      ! Exact f at Xyz_D [m]: the point source of Mass0 particles at
      ! Xyz0_D [m] after time tFinalTest; Width2 = 4*Dperp*tFinalTest [m^2]

      real, intent(in):: Xyz_D(3), Xyz0_D(3), Mass0, Width2
      !------------------------------------------------------------------------
      exact_distr = Mass0/(cPi*Width2)**1.5 &
           *exp(-sum((Xyz_D - Xyz0_D)**2)/Width2)

    end function exact_distr
    !==========================================================================
    subroutine set_mesh_dperp(nR, nTheta, nPhi)

      ! Make the mesh of the solver with its own interfaces: the
      ! #DIFFUSIONPERP command (see PARAM.XML) is written into a file, read
      ! and broadcast by ModReadParam and passed to read_param of
      ! SP_ModPerpDiffusion; then init_dperp makes the mesh, its face areas
      ! and cell volumes. The file is deleted after reading.

      integer, intent(in):: nR, nTheta, nPhi

      character(len=*), parameter:: NameParamFile = 'test_dperp_PARAM.in'
      character(len=*), parameter:: StringTab = char(9)//char(9)//char(9)
      character(len=lStringLine):: NameCommand
      character(len=*), parameter:: NameSub = 'set_mesh_dperp'
      !------------------------------------------------------------------------
      if(iProc == 0)then
         call open_file(File=NameParamFile, NameCaller=NameSub)
         write(UnitTmp_,'(a)')      '#DIFFUSIONPERP'
         write(UnitTmp_,'(a)')      'T'//StringTab//'UseDiffusionPerp'
         write(UnitTmp_,'(a)')      'UseDPerpRatio'//StringTab//'TypeDPerp'
         write(UnitTmp_,'(a)')      '0.01'//StringTab//'DPerpRatio'
         write(UnitTmp_,'(i0,a)')   nR, StringTab//'nRPerp'
         write(UnitTmp_,'(f0.2,a)') RMinTest, StringTab//'RMinPerp'
         write(UnitTmp_,'(f0.2,a)') RMaxTest, StringTab//'RMaxPerp'
         write(UnitTmp_,'(a)')      'Linear'//StringTab//'ScaleRPerp'
         write(UnitTmp_,'(i0,a)')   nTheta, StringTab//'nThetaPerp'
         write(UnitTmp_,'(i0,a)')   nPhi, StringTab//'nPhiPerp'
         write(UnitTmp_,'(a)')      'T'//StringTab//'DoTestDPerp'
         call close_file(NameCaller=NameSub)
      end if
      call read_file(NameParamFile, iComm)
      if(iProc == 0)then
         call open_file(File=NameParamFile, Status='old', NameCaller=NameSub)
         call close_file(Status='delete', NameCaller=NameSub)
      end if
      call read_init(nLineIn=n_line_read())
      do
         if(.not.read_line()) EXIT
         if(.not.read_command(NameCommand)) CYCLE
         call read_param_dperp(NameCommand)
      end do
      call init_dperp
      ! read_param changes an odd nTheta, and nR if it is less than nProc
      if(nRPerp /= nR .or. nThetaPerp /= nTheta .or. nPhiPerp /= nPhi) &
           call CON_stop(NameSub//': mesh changed by read_param, nR < nProc?')

    end subroutine set_mesh_dperp
    !==========================================================================
  end subroutine test_dperpUniSphere3d
  !============================================================================
end module ModTestDperp
!==============================================================================
program test_transport

  use ModTestDperp, ONLY: test_dperpUniSphere3d
  use ModUtilities, ONLY: CON_stop, lower_case
  use ModMpi
  use SP_ModProc, ONLY: iComm, iProc, nProc, iError

  implicit none

  character(len=20):: NameTest = ''
  !----------------------------------------------------------------------------
  ! MPI on 1 or more processors, as in SP_stand_alone
  call MPI_init(iError)
  iComm = MPI_COMM_WORLD
  call MPI_comm_rank(iComm, iProc, iError)
  call MPI_comm_size(iComm, nProc, iError)

  if(command_argument_count() > 0) call get_command_argument(1, NameTest)
  ! The test names are not case sensitive
  call lower_case(NameTest)
  select case(NameTest)
  case('', 'dperp', 'dperpunisphere3d')
     call test_dperpUniSphere3d
  case default
     call CON_stop('test_transport: unknown test '//trim(NameTest))
  end select

  call MPI_finalize(iError)

end program test_transport
!==============================================================================
