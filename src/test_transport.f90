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
!   ModTestDperp    perpendicular diffusion across the field lines:
!                   test_UniSphere3d (uniform Dperp, Linear mesh)
!                   test_R2Sphere3d  (Dperp = k*r**2, Exp mesh)
!                   test_Line2Sphere3d (field lines coupled with the
!                   Dperp mesh: Dpara along the lines, Dperp on the mesh)
!   ModTestFocused  focused transport along the field lines (later)
! Each test is one subroutine that contains its own helper routines.
! Usage: ./test_transport.exe [NameTest]; without NameTest all tests run.
! NameTest is not case sensitive: the name of one test without test_
! (e.g. R2Sphere3d), or dperp for all tests of ModTestDperp.
! The program runs on 1 or more MPI processors, e.g.
!   mpirun -np 4 ./test_transport.exe dperp
! and the results do not depend on the number of processors.
! Each test prints its error norms and writes test_<name>*.out files.
!==============================================================================
module ModTestDperp

  ! Tests of the perpendicular diffusion solver (SP_ModPerpDiffusion in
  ! src/ModDiffusionPerp.f90) against analytic solutions:
  ! - test_UniSphere3d: a point source with a uniform Dperp inside a
  !   sphere spreads as the 3D Gaussian of free space, and the number of
  !   particles is conserved
  ! - test_R2Sphere3d: a point source with Dperp = k*r**2 on the
  !   whole sphere spreads as a Gaussian in ln(r) times the heat kernel of
  !   the unit sphere (a Y_lm decays as exp(-l(l+1)*k*t)), and the number
  !   of particles is conserved
  ! - test_Line2Sphere3d: the coupled cycle of a run on radial field
  !   lines: diffusion along the lines (Dpara), then Dperp on the mesh
  !   (line -> mesh -> line), with Dpara and Dperp proportional to r**2;
  !   the solution is that of test_R2Sphere3d with different k in r and
  !   in angle
  ! The tests call the solvers of SP_ModPerpDiffusion (and the line solver
  ! of SP_ModAdvancePoisson); they do not have a copy of them.

  use ModMpi
  use ModNumConst, ONLY: cPi, cRadToDeg, cDegToRad
  use ModUtilities, ONLY: CON_stop, open_file, close_file
  use ModIoUnit, ONLY: UnitTmp_
  use ModPlotFile, ONLY: save_plot_file
  use ModReadParam, ONLY: read_file, read_init, read_line, &
       read_command, n_line_read, lStringLine, read_echo_set
  use SP_ModProc, ONLY: iComm, iProc, nProc, iError
  use SP_ModGrid, ONLY: nP, nLine
  use SP_ModUnit, ONLY: init_unit => init, Io2Si_V, UnitX_
  use SP_ModPerpDiffusion, ONLY: read_param_dperp => read_param, &
       init_dperp => init, solver_fvm_diffuse3d, nRPerp, nThetaPerp, &
       nPhiPerp, RPerp_C, ThetaPerp_C, PhiPerp_C, XyzPerp_CB, Volume_CB, &
       RPerp_F, ThetaPerp_F, PhiPerp_F

  implicit none

  private ! except

  public:: test_UniSphere3d ! point source, uniform Dperp in a sphere
  public:: test_R2Sphere3d  ! point source, Dperp = k*r**2, whole sphere
  public:: test_Line2Sphere3d ! lines + mesh, Dpara and Dperp ~ r**2

  ! Setup of test_UniSphere3d. Units: lengths in Rs (the unit of the
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
  subroutine test_UniSphere3d

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
    ! Outputs (on processor 0): a banner, then one block: the setup
    ! (Dperp, t, the width, nProc and the source cell), a table with one
    ! row per mesh (the base mesh and the 2x refined mesh): the wall time
    ! of the solver, the relative change of the number of particles
    ! sum(f*V) from t = 0 to tFinalTest, min f and two checks of the
    ! solver output, both 0 up to round-off: the same f at all momenta
    ! (max |f - f(p=1)|) and its df/dt times t = the change of f (max
    ! |DfDt*t - Df|, Df = f - f(t=0)); then the table of the relative
    ! max, L1 and L2 errors and the order of convergence (2 is expected).
    ! Files (f, exact f and error):
    !   test_UniSphere3d.out     phi-theta map on the source shell,
    !   test_UniSphere3d_r.out   radial cut through the source,
    !                                 both on the base mesh
    !   test_UniSphere3d_x2*.out the same for the refined mesh
    ! Any number of MPI processors: the solver updates the shells of each
    ! processor; the sums are made per shell and the plots are collected
    ! on processor 0, so the results are the same for any nProc.

    integer, parameter:: nMesh = 2
    ! Relative max, L1 and L2 errors on the base and the refined meshes
    real:: Error_VI(3,nMesh)
    integer:: iMesh
    ! Title of the banner
    character(len=*), parameter:: NameTitle = 'test_UniSphere3d: '// &
         'point source, Dperp = const, uniform spherical mesh'
    character(len=*), parameter:: NameSub = 'test_UniSphere3d'
    !--------------------------------------------------------------------------
    ! Banner (as in SP_stand_alone), the setup (the same for all meshes)
    ! and the header of the table of checks (run_UniSphere3d writes one
    ! row per mesh)
    if(iProc == 0)then
       write(*,'(/,4x,a,/,4x,a)') NameTitle, repeat('-', len(NameTitle))
       write(*,'(a,f6.3,a,f6.2,a,f6.3,a,i0)') ' Dperp =', DperpTest, &
            ' Rs^2/s, t =', tFinalTest, ' s, width sqrt(4*Dperp*t) =', &
            sqrt(4*DperpTest*tFinalTest), ' Rs, nProc = ', nProc
       write(*,'(a)') ' source cell (nR/2+1, nTheta/2+1, nPhi/2+1): f = 1'
       write(*,'(a)') '   nR nTheta  nPhi wall[s] dsum(fV)/sum      min f'// &
            ' max|f-f(p=1)| max|DfDt*t-Df|'
    end if

    ! Unit conversion: Io2Si_V(UnitX_) = Rs [m]
    call init_unit
    ! No field lines in this test (init_dperp allocates DPerp_CB with nLine)
    nLine = 0

    do iMesh = 1, nMesh
       call run_UniSphere3d(2**(iMesh-1), Error_VI(:,iMesh))
    end do

    if(iProc /= 0) RETURN
    write(*,'(a)') NameSub//': relative errors at t = tFinalTest'
    write(*,'(a)') '   nR nTheta  nPhi   max error    L1 error    L2 error'
    do iMesh = 1, nMesh
       write(*,'(i5,i7,i6,3es12.4)') nRTest*2**(iMesh-1), &
            nThetaTest*2**(iMesh-1), nPhiTest*2**(iMesh-1), Error_VI(:,iMesh)
    end do
    write(*,'(a,3f12.3)') ' order of convergence', &
         log(Error_VI(:,1)/Error_VI(:,2))/log(2.0)

  contains
    !==========================================================================
    subroutine run_UniSphere3d(nRefine, Error_V)

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
      ! Only the commented-out mesh line used NameSub (unused parameter)
      ! character(len=*), parameter:: NameSub = 'run_UniSphere3d'
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
         ! One row of the table of checks (header in test_UniSphere3d):
         ! the mesh, the wall time of the solver, the relative change of
         ! sum(f*V), min f, max |f - f(p=1)| and max |DfDt*t - Df|
         write(*,'(i5,i7,i6,f8.3,es13.3,es11.3,es14.3,es15.3)') nRPerp, &
              nThetaPerp, nPhiPerp, CpuTime, (Sum_V(2) - Sum_V(1))/Sum_V(1), &
              DistrMin, Max_I(3:4)
         ! The lines below are commented out: the row above has their
         ! values, or they are redundant.
         ! The mesh is in the row, nProc in the setup (test_UniSphere3d)
         ! write(*,'(/,a,i0,a,i0,a,i0,a,i0)') NameSub// &
         !      ': mesh nR x nTheta x nPhi = ', nRPerp, ' x ', nThetaPerp, &
         !      ' x ', nPhiPerp, ', nProc = ', nProc
         ! The setup gives the source cell; its center is in the .out files
         ! (rShell of the map, Theta and Phi of the _r cut)
         ! write(*,'(a,f7.3,a,2f8.3,a)') ' source cell at r =', &
         !      RPerp_C(iR0)/Rs, ' Rs, (theta, phi) =', &
         !      ThetaPerp_C(iTheta0)*cRadToDeg, PhiPerp_C(iPhi0)*cRadToDeg, &
         !      ' deg'
         ! The same on all meshes: in the setup (test_UniSphere3d)
         ! write(*,'(a,f6.3,a,f6.2,a,f6.3,a)') ' Dperp =', DperpTest, &
         !      ' Rs^2/s, t =', tFinalTest, ' s, width sqrt(4*Dperp*t) =', &
         !      sqrt(Width2)/Rs, ' Rs'
         ! In the row (wall[s])
         ! write(*,'(a,f9.3,a)') ' wall time of the solver   =', CpuTime, ' s'
         ! Absolute sums: the row has their relative change (dsum(fV)/sum)
         ! write(*,'(a,2es22.14)') ' sum(f*V) [Rs^3] at t=0, t =', &
         !      Sum_V(1:2)/Rs**3
         ! In the row (dsum(fV)/sum, min f, max|f-f(p=1)|, max|DfDt*t-Df|)
         ! write(*,'(a,es11.3)') ' relative change of sum(f*V)   =', &
         !      (Sum_V(2) - Sum_V(1))/Sum_V(1)
         ! write(*,'(a,es11.3)') ' min f                         =', DistrMin
         ! write(*,'(a,es11.3)') ' max |f(p) - f(p=1)|, p=1..nP  =', Max_I(3)
         ! write(*,'(a,es11.3)') ' max |DfDt*t - (f - f(t=0))|   =', Max_I(4)
         ! The same as the table of errors (test_UniSphere3d)
         ! write(*,'(a,3es11.3)') ' relative errors: max, L1, L2  =', Error_V
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
         NameFile = 'test_UniSphere3d'
         if(nRefine > 1) write(NameFile,'(a,i0)') &
              'test_UniSphere3d_x', nRefine
         call save_plot_file(NameFile=trim(NameFile)//'.out', &
              TypeFileIn='ascii', TimeIn=tFinalTest, nStepIn=1, &
              StringHeaderIn='test_UniSphere3d: f, exact f and '// &
              'f - exact f on the source shell, angles in deg', &
              NameVarIn='Phi Theta f fExact Error rShell Dperp', &
              ParamIn_I=[RPerp_C(iR0)/Rs, DperpTest], &
              Coord1In_I=PhiPerp_C*cRadToDeg, &
              Coord2In_I=ThetaPerp_C*cRadToDeg, VarIn_VII=Plot_VII)
         call save_plot_file(NameFile=trim(NameFile)//'_r.out', &
              TypeFileIn='ascii', TimeIn=tFinalTest, nStepIn=1, &
              StringHeaderIn='test_UniSphere3d: f, exact f and '// &
              'f - exact f along r through the source, r in Rs', &
              NameVarIn='r f fExact Error Theta Phi Dperp', &
              ParamIn_I=[ThetaPerp_C(iTheta0)*cRadToDeg, &
              PhiPerp_C(iPhi0)*cRadToDeg, DperpTest], &
              Coord1In_I=RPerp_C/Rs, VarIn_VI=Plot_VI)
      end if

      deallocate(Dperp_4D, Distr_4D, DfDt_4D, Sum_VI, Plot_VII, Plot_VI)

    end subroutine run_UniSphere3d
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
         ! F: T prints the sub-step and positivity diagnostics of
         ! solver_fvm_diffuse3d at every call; set T to debug
         write(UnitTmp_,'(a)')      'F'//StringTab//'DoTestDPerp'
         call close_file(NameCaller=NameSub)
      end if
      ! Quiet: no 'read and broadcast' line, and no echo of the
      ! parameters (read_file turns the echo on, on processor 0)
      call read_file(NameParamFile, iComm, IsVerbose=.false.)
      call read_echo_set(.false.)
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
  end subroutine test_UniSphere3d
  !============================================================================
  subroutine test_R2Sphere3d

    ! Point source test of the 3D finite volume Dperp solver
    ! solver_fvm_diffuse3d with Dperp = k*r**2 on the whole sphere, on its
    ! own r-theta-phi mesh (from init) with the exponential radial spacing
    ! (ScaleRPerp = Exp).
    !
    ! Setup: f = 1 in the cell that contains the point (rSource,
    ! ThetaSource, PhiSource) and f = 0 in all other cells, the same at
    ! all momenta. Dperp = kDperp*r**2 at every cell center and at every
    ! momentum (SI, r in m). One call of the solver advances f from t = 0
    ! to tFinal with its own stable sub-steps (time accurate mode).
    !
    ! Analytic solution: with D = k*r**2 and s = ln(r), d/dr = (1/r)d/ds,
    ! the radial part of div(D grad f) is
    !   (1/r**2) d/dr(r**2*D df/dr) = (k/r**3) d/ds(r**3 df/ds)
    !                               = k*(f_ss + 3*f_s)
    ! and the angular part is (D/r**2)*Lap_S f = k*Lap_S f, where Lap_S
    ! is the Laplacian on the unit sphere, so
    !   df/dt = k*(f_ss + 3*f_s) + k*Lap_S f.
    ! The radial and the angular operators commute, so the point source
    ! solution is the product of their point source solutions:
    !   f = M/rc**3 * G(s) * K(gamma),
    !   G = exp(-(s - ln(rc) + 3*k*t)**2/(4*k*t))/sqrt(4*pi*k*t),
    !   K = sum_{l=0..lmax} (2l+1)/(4*pi)*P_l(cos(gamma))*exp(-l(l+1)*k*t),
    ! where rc is the radius of the source point and gamma is the angle
    ! between the cell center and the source point (seen from r = 0).
    ! G solves G_t = k*(G_ss + 3*G_s): a Gaussian in s with the variance
    ! 2*k*t whose center moves to ln(rc) - 3*k*t. K is the heat kernel of
    ! the unit sphere: the Y_lm decay as exp(-l(l+1)*k*t), and the
    ! addition theorem sums the Y_lm of degree l into the Legendre
    ! polynomial P_l (here from the three-term recurrence). The sum stops
    ! at the first l with exp(-l(l+1)*k*t) < 1e-16 (l = 14 at k*t = 0.2).
    ! The number of particles is conserved: dV = r**3 ds dOmega, K
    ! integrates to 1 over the sphere, and the integral of G*exp(3*s) over
    ! s is the mean of exp(3*s) for the Gaussian G,
    ! exp(3*(ln(rc) - 3*k*t) + 9*k*t) = rc**3, so the integral of f dV is
    ! M at all times.
    ! The exact solution uses the center (rc, theta_c, phi_c) of the
    ! source cell as the source point and M = its volume (f = 1 there).
    ! Representing the source cell as a point is a part of the
    ! discretization error, which is of 2nd order. The exact solution is
    ! for all s, while the solver has no flux through the radial walls at
    ! RMin and RMax. The mass density f*r**3 of the exact solution is a
    ! Gaussian in s centered at ln(rc) + 3*k*t, about 1.9*sqrt(2*k*t)
    ! below ln(RMax) for rc = 40 Rs, so about 3% of M is beyond RMax in the
    ! exact solution and inside the mesh in the numerical one. The test
    ! prints sum(fExact*V)/M - 1 (the exact mass outside the mesh).
    !
    ! Outputs (on processor 0), as in test_UniSphere3d: a banner, then one
    ! block: the setup (k, t, k*t, nProc, the source point, M = sum(f*V)
    ! at t = 0), a table with one row per mesh: the wall time of the
    ! solver, sum(fExact*V)/M - 1, min f and the two checks of the solver
    ! output (max |f - f(p=1)| and max |DfDt*t - Df|, both 0 up to
    ! round-off); then the table of the relative max, L1 and L2 errors
    ! (the same norms) and the mass error (sum(f*V) - M)/M at tFinal, and
    ! the order of convergence (2 is expected) for each pair of meshes.
    ! Files for each mesh (f, exact f and error):
    !   test_R2Sphere3d_<nTheta>x<nPhi>.out       phi-theta map on
    !                                             the source shell,
    !   test_R2Sphere3d_<nTheta>x<nPhi>_r.out     cut along r and
    !   test_R2Sphere3d_<nTheta>x<nPhi>_theta.out cut along theta
    !                                             through the source.
    ! Any number of MPI processors, as in test_UniSphere3d.

    ! Units: lengths in Rs (the unit of the mesh in the #DIFFUSIONPERP
    ! command), angles in deg, time in s and kDperp in 1/s. The solver
    ! works in SI units: Dperp = kDperp*r**2 [m^2/s] with r in m.
    ! Final time and k: k*t = 0.2, so G has the width sqrt(2*k*t) = 0.63
    ! in ln(r), and K is about 36 deg wide.
    real, parameter:: tFinal = 5.0, kDperp = 0.04
    ! Meshes: whole sphere, Exp radial spacing from RMin to RMax. The base
    ! mesh nR x nTheta x nPhi is run and then refined by 2 in all 3
    ! directions. read_param makes an odd nTheta even, so the base mesh is
    ! 9 x 10 x 18 (not 9 x 9 x 18) and its 2x is 18 x 20 x 36.
    real, parameter:: RMin = 2.0, RMax = 240.0
    integer, parameter:: nPair = 1
    integer, parameter:: nCellBase_DI(3,nPair) = reshape( &
         [9, 10, 18], [3, nPair])
    ! The source cell contains this point: r [Rs], theta (colatitude) and
    ! phi [deg]. No face of the meshes is there.
    real, parameter:: rSource = 40.0, ThetaSource = 52.0, PhiSource = 185.0

    ! Relative max, L1, L2 and mass errors on the base and refined meshes
    real:: Error_VII(4,2,nPair)
    integer:: iPair, iLevel
    ! Title of the banner
    character(len=*), parameter:: NameTitle = 'test_R2Sphere3d: '// &
         'point source, Dperp = k*r**2, Exp spherical mesh'
    character(len=*), parameter:: NameSub = 'test_R2Sphere3d'
    !--------------------------------------------------------------------------
    ! Banner, the setup (the same for all meshes) and the header of the
    ! table of checks (run_R2Sphere3d writes one row per mesh)
    if(iProc == 0)then
       write(*,'(/,4x,a,/,4x,a)') NameTitle, repeat('-', len(NameTitle))
       write(*,'(a,f7.4,a,f6.2,a,f6.3,a,i0)') ' Dperp = k*r**2 with k =', &
            kDperp, ' 1/s, t =', tFinal, ' s, k*t =', kDperp*tFinal, &
            ', nProc = ', nProc
       write(*,'(3(a,f0.1),a)') ' source cell of r = ', rSource, &
            ' Rs, theta = ', ThetaSource, ', phi = ', PhiSource, &
            ' deg: f = 1, M = sum(f*V)'
       write(*,'(a)') '   nR nTheta  nPhi wall[s] sum(fExV)/M-1      min f'// &
            ' max|f-f(p=1)| max|DfDt*t-Df|'
    end if

    ! Unit conversion: Io2Si_V(UnitX_) = Rs [m]
    call init_unit
    ! No field lines in this test (init_dperp allocates DPerp_CB with nLine)
    nLine = 0

    do iPair = 1, nPair; do iLevel = 1, 2
       call run_R2Sphere3d(iLevel*nCellBase_DI(:,iPair), &
            Error_VII(:,iLevel,iPair))
    end do; end do

    if(iProc /= 0) RETURN
    write(*,'(a)') NameSub//': relative errors at t = tFinal'
    write(*,'(a)') '   nR nTheta  nPhi   max error    L1 error    L2 error'// &
         '  mass error'
    do iPair = 1, nPair
       do iLevel = 1, 2
          write(*,'(i5,i7,i6,4es12.4)') iLevel*nCellBase_DI(:,iPair), &
               Error_VII(:,iLevel,iPair)
       end do
       write(*,'(a,3f12.3)') ' order of convergence', &
            log(Error_VII(1:3,1,iPair)/Error_VII(1:3,2,iPair))/log(2.0)
    end do

  contains
    !==========================================================================
    subroutine run_R2Sphere3d(nCell_D, Error_V)

      ! One mesh of the test with nCell_D = (nR, nTheta, nPhi) cells. Set
      ! the mesh, f and Dperp, call the solver to tFinal, compare with the
      ! exact solution, print the results and save the plot files

      integer, intent(in):: nCell_D(3)
      ! Relative max, L1 and L2 errors and (sum(f*V) - M)/M at tFinal
      real,    intent(out):: Error_V(4)

      ! Shells of this processor: the same split as in init_dperp
      integer:: iRStart, iREnd, nRChunk, nRExtra
      ! Indexes of the source cell
      integer:: iR0, iTheta0, iPhi0
      ! Center of the source cell: r [m], theta and phi [rad], and the
      ! number of particles M in it (f = 1 times its volume [m^3])
      real:: r0, Theta0, Phi0, Mass0
      ! f on the shells of this processor at all momenta (with the ghost
      ! momenta 0 and nP+1), df/dt from the solver (averaged over its
      ! sub-steps), and Dperp on all shells (input of the solver)
      real, allocatable:: Distr_4D(:,:,:,:), DfDt_4D(:,:,:,:), &
           Dperp_4D(:,:,:,:)
      ! Sums over the cells of each shell, zero for the shells of the other
      ! processors (so the sums do not depend on nProc): sum(f*V) at t = 0
      ! and t = tFinal, sum(|f - fExact|*V), sum(fExact*V),
      ! sum((f - fExact)**2*V) and sum(fExact**2*V)
      integer, parameter:: nSum = 6
      real, allocatable:: Sum_VI(:,:)
      real:: Sum_V(nSum)
      ! Max |f - fExact|, max fExact, max |f - f(p=1)| over the momenta,
      ! max |DfDt*t - (f - f(t=0))| (DfDt*t is the change of f), and min f
      ! over all cells and momenta
      real:: Max_I(4), DistrMin
      ! Plot data: f, fExact and f - fExact on the source shell, along r
      ! and along theta through the source cell
      real, allocatable:: Plot_VII(:,:,:), PlotR_VI(:,:), PlotTheta_VI(:,:)
      ! Exact, numerical and initial f and volume of a cell, Rs [m],
      ! wall time [s]
      real:: Exact, Distr, Distr0, Volume, Rs, CpuTime
      integer:: i, j, k
      character(len=60):: NameFile
      ! Only the commented-out mesh line used NameSub (unused parameter)
      ! character(len=*), parameter:: NameSub = 'run_R2Sphere3d'
      !------------------------------------------------------------------------
      call set_mesh_dperp_r2(nCell_D(1), nCell_D(2), nCell_D(3))

      ! Shells of this processor: the first mod(nRPerp,nProc) processors
      ! have one more (as iRPerpStart_I and iRPerpEnd_I in init_dperp)
      nRChunk = nRPerp/nProc
      nRExtra = mod(nRPerp, nProc)
      iRStart = iProc*nRChunk + min(iProc, nRExtra) + 1
      iREnd   = iRStart + nRChunk - 1
      if(iProc < nRExtra) iREnd = iREnd + 1

      ! The source cell: count the lower faces at or below the source point
      Rs      = Io2Si_V(UnitX_)
      iR0     = count(RPerp_F(1:nRPerp) <= rSource*Rs)
      iTheta0 = count(ThetaPerp_F(1:nThetaPerp) <= ThetaSource*cDegToRad)
      iPhi0   = count(PhiPerp_F(1:nPhiPerp) <= PhiSource*cDegToRad)
      r0      = RPerp_C(iR0)
      Theta0  = ThetaPerp_C(iTheta0)
      Phi0    = PhiPerp_C(iPhi0)
      Mass0   = Volume_CB(iPhi0,iTheta0,iR0)

      ! Dperp = k*r**2 [m^2/s] at the cell centers of all shells (input of
      ! the solver), the same at all momenta
      allocate(Dperp_4D(nP,nPhiPerp,nThetaPerp,nRPerp))
      do i = 1, nRPerp
         Dperp_4D(:,:,:,i) = kDperp*RPerp_C(i)**2
      end do

      ! Initial f on the shells of this processor (all momenta)
      allocate(Distr_4D(0:nP+1,nPhiPerp,nThetaPerp,iRStart:iREnd), &
           DfDt_4D(0:nP+1,nPhiPerp,nThetaPerp,iRStart:iREnd))
      Distr_4D = 0.0
      if(iR0 >= iRStart .and. iR0 <= iREnd) &
           Distr_4D(:,iPhi0,iTheta0,iR0) = 1.0

      ! Number of particles at t = 0 (= Mass0)
      allocate(Sum_VI(nSum,nRPerp))
      Sum_VI = 0.0
      do i = iRStart, iREnd
         Sum_VI(1,i) = sum(Distr_4D(1,:,:,i)*Volume_CB(:,:,i))
      end do

      ! Advance f to tFinal with the solver (time accurate mode)
      call MPI_barrier(iComm, iError)
      CpuTime = MPI_wtime()
      call solver_fvm_diffuse3d(Dperp_4D, iRStart, iREnd, .false., &
           tFinal, Distr_4D, DfDt_4D)
      CpuTime = MPI_wtime() - CpuTime

      ! Compare with the exact solution on the shells of this processor.
      ! All momenta have the same f; check it and use p = 1. Check that
      ! the df/dt returned by the solver (averaged over its sub-steps)
      ! times tFinal is the change of f.
      Max_I = 0.0
      DistrMin = huge(1.0)
      do i = iRStart, iREnd; do j = 1, nThetaPerp; do k = 1, nPhiPerp
         Exact  = exact_distr_r2(RPerp_C(i), ThetaPerp_C(j), PhiPerp_C(k), &
              r0, Theta0, Phi0, Mass0)
         Distr  = Distr_4D(1,k,j,i)
         Volume = Volume_CB(k,j,i)
         Distr0 = 0.0
         if(i == iR0 .and. j == iTheta0 .and. k == iPhi0) Distr0 = 1.0
         Sum_VI(2:nSum,i) = Sum_VI(2:nSum,i) + Volume*[Distr, &
              abs(Distr - Exact), Exact, (Distr - Exact)**2, Exact**2]
         Max_I(1) = max(Max_I(1), abs(Distr - Exact))
         Max_I(2) = max(Max_I(2), Exact)
         Max_I(3) = max(Max_I(3), maxval(abs(Distr_4D(1:nP,k,j,i) - Distr)))
         Max_I(4) = max(Max_I(4), maxval(abs(DfDt_4D(1:nP,k,j,i)*tFinal &
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
      ! Relative errors as in test_UniSphere3d: max |f - fExact|/max
      ! fExact, the L1 error over the L1 norm of fExact and the same for L2
      ! (sums weighted by the volume); relative mass error at tFinal
      Error_V(1) = Max_I(1)/Max_I(2)
      Error_V(2) = Sum_V(3)/Sum_V(4)
      Error_V(3) = sqrt(Sum_V(5)/Sum_V(6))
      Error_V(4) = (Sum_V(2) - Mass0)/Mass0

      if(iProc == 0)then
         ! One row of the table of checks (header in test_R2Sphere3d): the
         ! mesh, the wall time of the solver, sum(fExact*V)/M - 1, min f,
         ! max |f - f(p=1)| and max |DfDt*t - Df|
         write(*,'(i5,i7,i6,f8.3,es14.3,es11.3,es14.3,es15.3)') nRPerp, &
              nThetaPerp, nPhiPerp, CpuTime, Sum_V(4)/Mass0 - 1, DistrMin, &
              Max_I(3:4)
         ! The lines below are commented out: the row above has their
         ! values, or they are redundant.
         ! The mesh is in the row, nProc in the setup (test_R2Sphere3d)
         ! write(*,'(/,a,i0,a,i0,a,i0,a,i0)') NameSub// &
         !      ': mesh nR x nTheta x nPhi = ', nRPerp, ' x ', nThetaPerp, &
         !      ' x ', nPhiPerp, ', nProc = ', nProc
         ! The setup gives the source point; the center of its cell is in
         ! the .out files (rShell of the map, Theta and Phi of the _r cut)
         ! write(*,'(a,f8.3,a,2f8.3,a)') ' source cell at r =', r0/Rs, &
         !      ' Rs, (theta, phi) =', Theta0*cRadToDeg, Phi0*cRadToDeg, ' deg'
         ! The same on all meshes: in the setup (test_R2Sphere3d)
         ! write(*,'(a,f7.4,a,f6.2,a,f6.3)') ' Dperp = k*r**2 with k =', &
         !      kDperp, ' 1/s, t =', tFinal, ' s, k*t =', kDperp*tFinal
         ! In the row (wall[s])
         ! write(*,'(a,f9.3,a)') ' wall time of the solver   =', CpuTime, ' s'
         ! Absolute sums: the table of errors has the mass error
         ! write(*,'(a,2es22.14)') ' sum(f*V) [Rs^3] at t=0 (=M), t =', &
         !      Sum_V(1:2)/Rs**3
         ! The same as the mass error column of the table of errors
         ! write(*,'(a,es11.3)') ' mass error (sum(f*V) - M)/M   =', Error_V(4)
         ! In the row (sum(fExV)/M-1, min f, max|f-f(p=1)|, max|DfDt*t-Df|)
         ! write(*,'(a,es11.3)') ' sum(fExact*V)/M - 1           =', &
         !      Sum_V(4)/Mass0 - 1
         ! write(*,'(a,es11.3)') ' min f                         =', DistrMin
         ! write(*,'(a,es11.3)') ' max |f(p) - f(p=1)|, p=1..nP  =', Max_I(3)
         ! write(*,'(a,es11.3)') ' max |DfDt*t - (f - f(t=0))|   =', Max_I(4)
         ! The same as the table of errors (test_R2Sphere3d)
         ! write(*,'(a,3es11.3)') ' relative errors: max, L1, L2  =', &
         !      Error_V(1:3)
      end if

      ! Plots of f, exact f and error: phi-theta map on the source shell
      ! (from the processor of that shell) and radial cut through the
      ! source (from all processors), collected (summed with zeros) on all
      ! processors; the theta cut through the source is taken from the map
      ! (the meridian of the source cell). Processor 0 writes the files.
      allocate(Plot_VII(3,nPhiPerp,nThetaPerp), PlotR_VI(3,nRPerp), &
           PlotTheta_VI(3,nThetaPerp))
      Plot_VII = 0.0; PlotR_VI = 0.0
      if(iR0 >= iRStart .and. iR0 <= iREnd)then
         do j = 1, nThetaPerp; do k = 1, nPhiPerp
            Distr = Distr_4D(1,k,j,iR0)
            Exact = exact_distr_r2(r0, ThetaPerp_C(j), PhiPerp_C(k), &
                 r0, Theta0, Phi0, Mass0)
            Plot_VII(:,k,j) = [Distr, Exact, Distr - Exact]
         end do; end do
      end if
      do i = iRStart, iREnd
         Distr = Distr_4D(1,iPhi0,iTheta0,i)
         Exact = exact_distr_r2(RPerp_C(i), Theta0, Phi0, &
              r0, Theta0, Phi0, Mass0)
         PlotR_VI(:,i) = [Distr, Exact, Distr - Exact]
      end do
      call MPI_allreduce(MPI_IN_PLACE, Plot_VII, size(Plot_VII), MPI_REAL, &
           MPI_SUM, iComm, iError)
      call MPI_allreduce(MPI_IN_PLACE, PlotR_VI, size(PlotR_VI), MPI_REAL, &
           MPI_SUM, iComm, iError)
      PlotTheta_VI = Plot_VII(:,iPhi0,:)
      if(iProc == 0)then
         write(NameFile,'(a,i0,a,i0)') 'test_R2Sphere3d_', &
              nThetaPerp, 'x', nPhiPerp
         call save_plot_file(NameFile=trim(NameFile)//'.out', &
              TypeFileIn='ascii', TimeIn=tFinal, nStepIn=1, &
              StringHeaderIn='test_R2Sphere3d: f, exact f and '// &
              'f - exact f on the source shell, angles in deg', &
              NameVarIn='Phi Theta f fExact Error rShell kDperp', &
              ParamIn_I=[r0/Rs, kDperp], &
              Coord1In_I=PhiPerp_C*cRadToDeg, &
              Coord2In_I=ThetaPerp_C*cRadToDeg, VarIn_VII=Plot_VII)
         call save_plot_file(NameFile=trim(NameFile)//'_r.out', &
              TypeFileIn='ascii', TimeIn=tFinal, nStepIn=1, &
              StringHeaderIn='test_R2Sphere3d: f, exact f and '// &
              'f - exact f along r through the source, r in Rs', &
              NameVarIn='r f fExact Error Theta Phi kDperp', &
              ParamIn_I=[Theta0*cRadToDeg, Phi0*cRadToDeg, kDperp], &
              Coord1In_I=RPerp_C/Rs, VarIn_VI=PlotR_VI)
         call save_plot_file(NameFile=trim(NameFile)//'_theta.out', &
              TypeFileIn='ascii', TimeIn=tFinal, nStepIn=1, &
              StringHeaderIn='test_R2Sphere3d: f, exact f and '// &
              'f - exact f along theta through the source, theta in deg', &
              NameVarIn='Theta f fExact Error rShell Phi kDperp', &
              ParamIn_I=[r0/Rs, Phi0*cRadToDeg, kDperp], &
              Coord1In_I=ThetaPerp_C*cRadToDeg, VarIn_VI=PlotTheta_VI)
      end if

      deallocate(Dperp_4D, Distr_4D, DfDt_4D, Sum_VI, Plot_VII, PlotR_VI, &
           PlotTheta_VI)

    end subroutine run_R2Sphere3d
    !==========================================================================
    real function exact_distr_r2(r, Theta, Phi, r0, Theta0, Phi0, Mass0)

      ! Exact f at the point (r [m], Theta, Phi [rad]) at t = tFinal for
      ! the point source of Mass0 particles at (r0, Theta0, Phi0):
      ! Mass0/r0**3*G(ln(r))*K(gamma), see test_R2Sphere3d

      real, intent(in):: r, Theta, Phi, r0, Theta0, Phi0, Mass0

      ! k*t, cos(gamma), Legendre polynomials P_(l-1), P_l and P_(l+1),
      ! exp(-l(l+1)*k*t), and K(gamma) times 4*pi
      real:: kt, CosGamma, LegendreMinus, Legendre, LegendrePlus, Decay, &
           SumK
      integer:: l
      !------------------------------------------------------------------------
      kt = kDperp*tFinal
      ! Angle gamma between the point and the source point (from r = 0)
      CosGamma = max(-1.0, min(1.0, cos(Theta)*cos(Theta0) &
           + sin(Theta)*sin(Theta0)*cos(Phi - Phi0)))
      ! Sum of (2l+1)*P_l(cos gamma)*exp(-l(l+1)*k*t) over l, with P_0 = 1,
      ! P_1 = x and (l+1)*P_(l+1) = (2l+1)*x*P_l - l*P_(l-1), up to the
      ! first l with exp(-l(l+1)*k*t) < 1e-16
      SumK = 0.0
      LegendreMinus = 0.0
      Legendre = 1.0
      l = 0
      do
         Decay = exp(-l*(l+1)*kt)
         if(Decay < 1e-16) EXIT
         SumK = SumK + (2*l+1)*Legendre*Decay
         LegendrePlus = ((2*l+1)*CosGamma*Legendre - l*LegendreMinus)/(l+1)
         LegendreMinus = Legendre
         Legendre = LegendrePlus
         l = l + 1
      end do
      exact_distr_r2 = Mass0/r0**3 &
           *exp(-(log(r/r0) + 3*kt)**2/(4*kt))/sqrt(4*cPi*kt) &
           *SumK/(4*cPi)

    end function exact_distr_r2
    !==========================================================================
    subroutine set_mesh_dperp_r2(nR, nTheta, nPhi)

      ! Make the mesh of the solver with its own interfaces, as
      ! set_mesh_dperp of test_UniSphere3d, but with the exponential
      ! radial spacing (ScaleRPerp = Exp) from RMin to RMax: the
      ! #DIFFUSIONPERP command (see PARAM.XML) is written into a file, read
      ! and broadcast by ModReadParam and passed to read_param of
      ! SP_ModPerpDiffusion; then init_dperp makes the mesh, its face areas
      ! and cell volumes. The file is deleted after reading.

      integer, intent(in):: nR, nTheta, nPhi

      character(len=*), parameter:: NameParamFile = 'test_dperp_PARAM.in'
      character(len=*), parameter:: StringTab = char(9)//char(9)//char(9)
      character(len=lStringLine):: NameCommand
      character(len=*), parameter:: NameSub = 'set_mesh_dperp_r2'
      !------------------------------------------------------------------------
      if(iProc == 0)then
         call open_file(File=NameParamFile, NameCaller=NameSub)
         write(UnitTmp_,'(a)')      '#DIFFUSIONPERP'
         write(UnitTmp_,'(a)')      'T'//StringTab//'UseDiffusionPerp'
         write(UnitTmp_,'(a)')      'UseDPerpRatio'//StringTab//'TypeDPerp'
         write(UnitTmp_,'(a)')      '0.01'//StringTab//'DPerpRatio'
         write(UnitTmp_,'(i0,a)')   nR, StringTab//'nRPerp'
         write(UnitTmp_,'(f0.2,a)') RMin, StringTab//'RMinPerp'
         write(UnitTmp_,'(f0.2,a)') RMax, StringTab//'RMaxPerp'
         write(UnitTmp_,'(a)')      'Exp'//StringTab//'ScaleRPerp'
         write(UnitTmp_,'(i0,a)')   nTheta, StringTab//'nThetaPerp'
         write(UnitTmp_,'(i0,a)')   nPhi, StringTab//'nPhiPerp'
         ! F: T prints the sub-step and positivity diagnostics of
         ! solver_fvm_diffuse3d at every call; set T to debug
         write(UnitTmp_,'(a)')      'F'//StringTab//'DoTestDPerp'
         call close_file(NameCaller=NameSub)
      end if
      ! Quiet: no 'read and broadcast' line, and no echo of the
      ! parameters (read_file turns the echo on, on processor 0)
      call read_file(NameParamFile, iComm, IsVerbose=.false.)
      call read_echo_set(.false.)
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
           call CON_stop(NameSub// &
           ': mesh changed by read_param (odd nTheta or nR < nProc?)')

    end subroutine set_mesh_dperp_r2
    !==========================================================================
  end subroutine test_R2Sphere3d
  !============================================================================
  subroutine test_Line2Sphere3d

    ! Test of the coupled solvers of a run: diffusion along the field lines
    ! (Dpara) and the 3D finite volume Dperp solver on the uniform
    ! r-theta-phi mesh, coupled in each time step as in the time-accurate
    ! advance of SP_ModAdvance (src/ModAdvance.f90):
    !   1. on each line of this processor: advect_via_poisson_parker
    !      (#ADVECTION T; here with no advection), which calls
    !      diffuse_distribution (implicit Dpara along the line),
    !   2. diffuseperp_distribution: f on the lines is interpolated to the
    !      mesh (intersections of the lines with the shells, triangulation
    !      of each shell), solver_fvm_diffuse3d advances the mesh f with
    !      Dperp, and the change of f is interpolated back and applied on
    !      the line nodes.
    ! The test calls these routines in this order; it has no copy of them.
    ! It does not call advance itself: advance gets Dpara from
    ! set_diffusion_coef, which computes it from the turbulence (growing as
    ! p**(1/3)*v with the momentum), so the test sets the outputs of
    ! set_diffusion_coef on each line before the line solver instead.
    !
    ! Setup (lengths in Rs, angles in deg, time in s, k in 1/s; the
    ! solvers use SI units, D = k*r**2 [m^2/s] with r in m):
    ! - Lines: straight radial lines from RMin = 2 to RMax = 240 through
    !   the cell centers of a lat-lon partition of the sphere, nThetaLine
    !   rows in theta (colatitude, row 1 at the north pole) and nPhiLine
    !   columns in phi, numbered and split between the processors as in a
    !   run (ilineall_to_lon_lat and init of SP_ModGrid), with nNode nodes
    !   per line uniform in ln(r) and B = B0*(RMin/r)**2, so that the
    !   diffusion along a line is the radial part of the 3D diffusion.
    !   Base: 18 x 36 lines (10 deg) with 36 nodes; refined: 36 x 72 lines
    !   (5 deg) with 72 nodes. The background (MhData_VIB) is constant.
    ! - Mesh: as in test_R2Sphere3d, 9 x 10 x 18 cells for the base lines
    !   and 18 x 20 x 36 for the refined lines (Exp in r, 2 to 240 Rs).
    ! - Coefficients, the same at all momenta: Dpara = kDpara*r**2 is
    !   split with DPerpRatio (#DIFFUSIONPERP) as in set_diffusion_coef:
    !   the line solver gets Dpara*(1 - DPerpRatio) (DInnerSi_II times
    !   DOuterSi_I = B) and the mesh gets Dperp = DPerpRatio*Dpara
    !   (DPerp_CB on the line nodes), so the total radial diffusion is
    !   Dpara. Three cases:
    !     coupled:    kDpara = kPara = 0.04, DPerpRatio = kPerp/kPara = 0.5
    !                 (Dperp = kPerp*r**2 with kPerp = 0.02),
    !     Dpara only: kDpara = kPara, DPerpRatio = 0 (Dperp = 0),
    !     Dperp only: kDpara = kPerp, DPerpRatio = 1 (the line Dpara is 0,
    !                 raised by diffuse_distribution to its floor of
    !                 1e4 Rs [m^2/s], 2e-4 of kPerp*r**2 at r = 2 Rs).
    !   No other terms: the background does not change (no advection, no
    !   adiabatic terms) and nMu = 1 (no focusing, no scattering).
    ! - Boundary conditions (#LOWERENDBC, #UPPERENDBC, #MOMENTUMBC): float
    !   at both ends of the lines and of the momentum grid. The ghost f at
    !   the line ends is at least the background f, which #FLUXINITIAL
    !   1e-200 makes negligible. The mesh has zero gradients at its radial
    !   walls.
    ! - Time: nStep coupled steps of Dt = tFinal/nStep with tFinal = 5 s:
    !   5 steps (Dt = 1 s) with the base lines and 20 steps with the
    !   refined lines (the line solver is first order in time). The line
    !   solver makes one step of Dt per coupled step (no advection, so no
    !   CFL limit); solver_fvm_diffuse3d makes its own stable sub-steps
    !   (printed if DoTestDPerp is T in set_param_line).
    !
    ! Exact solution: in s = ln(r) the 3D operator is
    !   kR*(f_ss + 3*f_s) + kA*Lap_S f
    ! (see test_R2Sphere3d) with kR = kDpara (lines and mesh) and
    ! kA = DPerpRatio*kDpara (mesh), so the point source solution is
    !   f = M/rc**3 * G(s; TauR) * K(gamma; TauA),
    ! the G and K of test_R2Sphere3d with k*t replaced by TauR = kR*t and
    ! TauA = kA*t. All cases start at t0 from the exact f of the coupled
    ! case with kPara*t0 = 0.05 (TauR = 0.05, TauA = 0.025), so after the
    ! run TauR = 0.05 + kR*tFinal and TauA = 0.025 + kA*tFinal (TauA stays
    ! 0.025 for Dpara only). Source rc = 40 Rs, theta = 52 deg, phi = 185
    ! deg; M makes the peak of f at t0 equal to 1. The initial f is the
    ! exact f at the nodes at all momenta, at least DistrFloor = 1e-100
    ! (where the Legendre sum is round-off, below 1e-15 and possibly
    ! negative). The exact solution is for
    ! all s, while the lines and the mesh end at 2 and 240 Rs: at tFinal
    ! the exact f of the coupled case at r = 2 (240) Rs is 6.4e-3 (1.6e-3)
    ! of its maximum along the same line, so the ends add errors.
    !
    ! Outputs (on processor 0), as in test_UniSphere3d: a banner, then one
    ! block: the setup (the base lines and mesh, nProc, the source point,
    ! tFinal, TauR and TauA of f(t0), M0 = sum(f*V) at t0, and for each
    ! case the Dpara of the lines and Dperp over r**2 and TauR and TauA at
    ! t0 + tFinal), a table with one row per case and set of lines: nNode,
    ! the number of steps, the wall time of the steps, the number of lines
    ! in use, sum(fExact*V)/M0 - 1, min f and max |f(p) - f(p=1)| over the
    ! momenta; then the table of the relative max, L1 and L2 errors over
    ! all line nodes at t0 + tFinal: max |f - fExact|/max fExact and the
    ! sums of |f - fExact| and (f - fExact)**2 over those of fExact and
    ! fExact**2, weighted by the volume V = r**3*dln(r)*dOmega of each node
    ! (half at the line ends; dOmega of the lat-lon cell of the line), and
    ! the relative change of sum(f*V), with the orders of convergence
    ! (base -> refined). Files for the coupled case (f, exact f and error):
    !   test_Line2Sphere3d_<nThetaLine>x<nPhiLine>.out       phi-theta map
    !                                        at the node nearest to rc,
    !   test_Line2Sphere3d_<nThetaLine>x<nPhiLine>_r.out     along the line
    !                                        nearest to the source,
    !   test_Line2Sphere3d_<nThetaLine>x<nPhiLine>_theta.out along theta at
    !                                        the phi and node of the above.
    ! Any number of MPI processors up to the number of lines: the sums are
    ! made per line, the plots are collected on processor 0, and f does
    ! not depend on nProc.
    !
    ! Memory: as in a run, the line arrays (MhData_VIB, State_VIB,
    ! Distribution_CB, DPerp_CB and the time step of
    ! diffuseperp_distribution) have nVertexMax = 20000 nodes per line,
    ! 52 MB per line. init_stand_alone and the init routines of
    ! SP_ModDistribution and SP_ModPerpDiffusion would fill all of them
    ! (135 GB for 36 x 72 lines), so the test calls these init routines
    ! with nLine = 0 and allocates the line arrays with the same shapes. It
    ! sets the nodes in use only; the memory of the others is not touched.

    use SP_ModSize, ONLY: nVertexMax
    use SP_ModGrid, ONLY: nLat, nLon, nLineAll, iLineAll0, nMu, &
         nVertex_B, Used_B, MhData_VIB, State_VIB, LagrID_, X_, Z_, Rho_, &
         Bx_, Bz_, nMhData, nVar, R_, D_, B_, NoShock_, &
         ilineall_to_lon_lat, get_other_state_var
    use SP_ModDistribution, ONLY: init_dist => init, &
         read_param_dist => read_param, Distribution_CB, IsDistNeg
    use SP_ModBc, ONLY: read_param_bc => read_param
    use SP_ModDiffusion, ONLY: DInnerSi_II, DOuterSi_I
    use SP_ModTurbulence, ONLY: DPerp_CB
    use SP_ModPerpDiffusion, ONLY: diffuseperp_distribution, DPerpRatio
    use SP_ModAdvancePoisson, ONLY: advect_via_poisson_parker

    ! Coefficients [1/s] of the coupled case, start time t0 [s] with
    ! kPara*t0 = 0.05, and the run time [s]
    real, parameter:: kPara = 0.04, kPerp = 0.02, tStart = 0.05/kPara, &
         tFinal = 5.0
    ! Lines from RMin to RMax [Rs]
    real, parameter:: RMin = 2.0, RMax = 240.0
    ! Base lines: nodes per line, rows in theta, columns in phi; the
    ! refined lines have twice as many of each
    integer, parameter:: nNodeBase = 36, nThetaLineBase = 18, &
         nPhiLineBase = 36
    ! Base mesh nR x nTheta x nPhi (as test_R2Sphere3d), refined by 2
    integer, parameter:: nCellBase_D(3) = [9, 10, 18]
    ! Coupled steps for the base and the refined lines
    integer, parameter:: nStep_I(2) = [5, 20]
    ! Source point: r [Rs], theta (colatitude) and phi [deg]
    real, parameter:: rSource = 40.0, ThetaSource = 52.0, PhiSource = 185.0
    ! Cases: names, kDpara and DPerpRatio
    integer, parameter:: nCase = 3
    character(len=10), parameter:: NameCase_I(nCase) = &
         ['coupled   ', 'Dpara only', 'Dperp only']
    real, parameter:: kDpara_I(nCase) = [kPara, kPara, kPerp]
    real, parameter:: DPerpRatio_I(nCase) = [kPerp/kPara, 0.0, 1.0]
    ! Background along the lines at RMin: B [T] and density [m^-3]
    real, parameter:: B0 = 1e-5, Rho0 = 1e10
    ! CFL number of the line solver (default of #CFL in SP_ModAdvance)
    real, parameter:: Cfl = 0.9
    ! Floor of the initial f: it replaces the round-off of the Legendre
    ! sum (below 1e-15, also negative) where the exact f is tiny, and
    ! keeps the ghost f at the line ends above the background f (1e-203),
    ! so f stays the same at all momenta
    real, parameter:: DistrFloor = 1e-100

    ! Number of particles M of the exact solution
    real:: Mass
    ! Lines of the current run: nodes per line, rows and columns, and the
    ! distance of the nodes in ln(r)
    integer:: nNode, nThetaLine, nPhiLine
    real:: dLogNode
    ! The test has allocated the line arrays (pointers of SP_ModGrid)
    logical:: IsLineAllocated
    ! Relative max, L1, L2 errors and the relative change of sum(f*V) for
    ! each case on the base and refined lines
    real:: Error_VII(4,nCase,2)
    integer:: iCase, iLevel
    ! Title of the banner
    character(len=*), parameter:: NameTitle = 'test_Line2Sphere3d: '// &
         'radial lines (Dpara) coupled with Exp mesh (Dperp)'
    character(len=*), parameter:: NameSub = 'test_Line2Sphere3d'
    !--------------------------------------------------------------------------
    ! Banner, the setup (the same for all runs), the cases and the header
    ! of the table of checks (run_Line2Sphere3d writes one row per run)
    if(iProc == 0)then
       write(*,'(/,4x,a,/,4x,a)') NameTitle, repeat('-', len(NameTitle))
       write(*,'(7(a,i0))') ' lines ', nThetaLineBase, ' x ', &
            nPhiLineBase, ' with ', nNodeBase, ' nodes, mesh ', &
            nCellBase_D(1), ' x ', nCellBase_D(2), ' x ', nCellBase_D(3), &
            '; refined 2x in all; nProc = ', nProc
       write(*,'(3(a,f0.1),a,f0.2,a)') ' source point r = ', rSource, &
            ' Rs, theta = ', ThetaSource, ', phi = ', PhiSource, &
            ' deg, tFinal = ', tFinal, ' s'
       write(*,'(a,2f8.4,a)') ' f(t0) = exact f (peak 1) at TauR, TauA =', &
            kPara*tStart, kPerp*tStart, '; M0 = sum(f(t0)*V)'
       ! Each case: the Dpara of the lines and Dperp (as in
       ! run_Line2Sphere3d) over r**2, TauR and TauA at t0 + tFinal
       write(*,'(a)') ' case       line Dpara/r**2, Dperp/r**2 [1/s], '// &
            'TauR, TauA at t0 + tFinal'
       do iCase = 1, nCase
          write(*,'(1x,a,4f8.4)') NameCase_I(iCase), &
               (1 - DPerpRatio_I(iCase))*kDpara_I(iCase), &
               DPerpRatio_I(iCase)*kDpara_I(iCase), &
               kPara*tStart + kDpara_I(iCase)*tFinal, &
               kPerp*tStart + DPerpRatio_I(iCase)*kDpara_I(iCase)*tFinal
       end do
       write(*,'(a)') ' case       nNode nStep wall[s] lines '// &
            'sum(fExV)/M0-1       min f max|f-f(p=1)|'
    end if

    ! Unit conversion: Io2Si_V(UnitX_) = Rs [m]
    call init_unit
    IsLineAllocated = .false.
    ! M such that the peak of the exact f at t0 (at r = rc*exp(-3*TauR),
    ! gamma = 0) is 1
    Mass = 1.0
    Mass = 1.0/exact_line(rSource*exp(-3*kPara*tStart), &
         ThetaSource*cDegToRad, PhiSource*cDegToRad, kPara*tStart, &
         kPerp*tStart)

    do iLevel = 1, 2; do iCase = 1, nCase
       call run_Line2Sphere3d(iLevel, iCase, Error_VII(:,iCase,iLevel))
    end do; end do

    ! Free the line arrays
    deallocate(MhData_VIB, State_VIB, nVertex_B, Used_B, Distribution_CB, &
         DPerp_CB)
    nLine = 0

    if(iProc /= 0) RETURN
    write(*,'(a)') NameSub//': relative errors at t0 + tFinal on the '// &
         'line nodes (L1, L2 weighted by V = r**3*dln(r)*dOmega)'
    write(*,'(a)') ' case        nNode nTheta  nPhi   max error'// &
         '    L1 error    L2 error dsum(fV)/sum'
    do iCase = 1, nCase
       do iLevel = 1, 2
          write(*,'(1x,a,i6,i7,i6,4es12.4)') NameCase_I(iCase), &
               iLevel*nNodeBase, iLevel*nThetaLineBase, &
               iLevel*nPhiLineBase, Error_VII(:,iCase,iLevel)
       end do
       write(*,'(a,3f12.3)') ' order of convergence      ', &
            log(Error_VII(1:3,iCase,1)/Error_VII(1:3,iCase,2))/log(2.0)
    end do

  contains
    !==========================================================================
    subroutine run_Line2Sphere3d(iLevel, iCase, Error_V)

      ! One case (iCase) on one set of lines (iLevel = 1: base, 2:
      ! refined): read the parameters, set the lines and the initial f, make
      ! the coupled steps to tFinal, compare with the exact solution, print
      ! the results and save the plot files (coupled case only)

      integer, intent(in):: iLevel, iCase
      ! Relative max, L1 and L2 errors and the relative change of sum(f*V)
      real,    intent(out):: Error_V(4)

      ! kDpara of the case [1/s], time step [s], radial and angular k*t
      ! accumulated at t0 + tFinal
      real:: kDpara, Dt, TauR, TauA
      integer:: nStep, iStep
      ! Data along one line in SI units, as in advance: distance to the
      ! next node, density times the length of a node over B, density, B
      ! and r
      real, allocatable:: DsSi_I(:), Mass_C(:), nSi_I(:), BSi_I(:), rSi_I(:)
      ! Time step of diffuseperp_distribution at the line nodes
      real, allocatable:: DtLine_CB(:,:,:)
      ! Sums over the nodes of each line, zero for the lines of the other
      ! processors (so the sums do not depend on nProc): sum(f*V) at t0
      ! and t0 + tFinal, sum(|f - fExact|*V), sum(fExact*V),
      ! sum((f - fExact)**2*V) and sum(fExact**2*V)
      integer, parameter:: nSum = 6
      real, allocatable:: Sum_VI(:,:)
      real:: Sum_V(nSum)
      ! Max |f - fExact|, max fExact and max |f(p) - f(p=1)| over the
      ! momenta, min f over all nodes and momenta, and the lines in use
      real:: Max_I(3), DistrMin
      integer:: nLineUsed
      ! Plot data: f, fExact and f - fExact at the node nearest to rc on
      ! all lines, along the line nearest to the source, and along theta
      real, allocatable:: Plot_VII(:,:,:), PlotR_VI(:,:), PlotTheta_VI(:,:)
      ! Row, column and index of the line nearest to the source, and its
      ! node nearest to rc
      integer:: iLat0, iLon0, iLineSource, iX0
      ! Exact and numerical f, volume of a node [Rs^3], angles of a line
      ! [rad], Rs [m], wall time [s]
      real:: Exact, Distr, Volume, Theta, Phi, Rs, CpuTime
      integer:: iLine, iLineAll, iLat, iLon, iX, iP, iEnd, i
      character(len=60):: NameFile
      character(len=*), parameter:: NameSub = 'run_Line2Sphere3d'
      !------------------------------------------------------------------------
      kDpara = kDpara_I(iCase)
      nNode = iLevel*nNodeBase
      nThetaLine = iLevel*nThetaLineBase
      nPhiLine = iLevel*nPhiLineBase
      nStep = nStep_I(iLevel)
      Dt = tFinal/nStep
      Rs = Io2Si_V(UnitX_)

      ! Parameters, momentum grid and mesh; then the lines and f at t0
      call set_param_line(iLevel*nCellBase_D, DPerpRatio_I(iCase))
      call set_lines

      ! Number of particles at t0 (sum over the nodes of each line)
      allocate(Sum_VI(nSum,nLineAll))
      Sum_VI = 0.0
      do iLine = 1, nLine
         iLineAll = iLineAll0 + iLine
         do iX = 1, nNode
            Sum_VI(1,iLineAll) = Sum_VI(1,iLineAll) &
                 + volume_node(iLine, iX)*Distribution_CB(1,1,iX,iLine)
         end do
      end do

      ! Time step of the Dperp update on the nodes in use
      allocate(DtLine_CB(nP,nVertexMax,nLine))
      DtLine_CB(:,1:nNode,:) = Dt
      allocate(DsSi_I(nNode), Mass_C(nNode), nSi_I(nNode), BSi_I(nNode), &
           rSi_I(nNode))

      ! Coupled steps, as advance does each time step
      call MPI_barrier(iComm, iError)
      CpuTime = MPI_wtime()
      do iStep = 1, nStep
         ! 1. Line solver on the lines of this processor
         do iLine = 1, nLine
            if(.not.Used_B(iLine)) CYCLE
            iEnd = nVertex_B(iLine)
            ! Data along the line in SI units, as in advance (the
            ! background does not change: nOldSi_I = nSi_I, BOld = B)
            DsSi_I(1:iEnd-1) = State_VIB(D_,1:iEnd-1,iLine)*Rs
            Mass_C(2:iEnd-1) = 0.5*(DsSi_I(2:iEnd-1) + DsSi_I(1:iEnd-2)) &
                 *MhData_VIB(Rho_,2:iEnd-1,iLine)/State_VIB(B_,2:iEnd-1,iLine)
            Mass_C(1) = DsSi_I(1)*MhData_VIB(Rho_,1,iLine) &
                 /State_VIB(B_,1,iLine)
            Mass_C(iEnd) = DsSi_I(iEnd-1)*MhData_VIB(Rho_,iEnd,iLine) &
                 /State_VIB(B_,iEnd,iLine)
            nSi_I = MhData_VIB(Rho_,1:iEnd,iLine)
            BSi_I = State_VIB(B_,1:iEnd,iLine)
            rSi_I = State_VIB(R_,1:iEnd,iLine)*Rs
            ! The outputs of set_diffusion_coef for Dpara = kDpara*r**2:
            ! DOuter = B, DInner = Dpara*(1 - DPerpRatio)/B on the line,
            ! and Dperp = DPerpRatio*Dpara for the mesh
            DOuterSi_I(1:iEnd) = BSi_I
            do iP = 1, nP
               DInnerSi_II(1:iEnd,iP) = &
                    (1.0 - DPerpRatio)*kDpara*rSi_I**2/BSi_I
               DPerp_CB(iP,1:iEnd,iLine) = DPerpRatio*kDpara*rSi_I**2
            end do
            ! Line solver of advance (#ADVECTION T, IsMuAvg)
            call advect_via_poisson_parker(iLine, iEnd, NoShock_, Dt, Cfl, &
                 nSi_I, nSi_I, BSi_I, Mass_C)
            if(IsDistNeg) write(*,*) NameSub//': f < 0 on line ', &
                 iLineAll0 + iLine
         end do
         ! 2. Dperp on the mesh and back to the line nodes (all processors)
         call diffuseperp_distribution(IsFirstCall=.true., &
              dtIn_CB=DtLine_CB)
      end do
      CpuTime = MPI_wtime() - CpuTime

      ! Compare with the exact solution on the lines of this processor.
      ! All momenta have the same f; check it and use p = 1.
      TauR = kPara*tStart + kDpara*tFinal
      TauA = kPerp*tStart + DPerpRatio*kDpara*tFinal
      Max_I = 0.0
      DistrMin = huge(1.0)
      do iLine = 1, nLine
         if(.not.Used_B(iLine)) CYCLE
         iLineAll = iLineAll0 + iLine
         call line_angles(iLine, iLat, iLon, Theta, Phi)
         do iX = 1, nVertex_B(iLine)
            Exact  = exact_line(State_VIB(R_,iX,iLine), Theta, Phi, TauR, &
                 TauA)
            Distr  = Distribution_CB(1,1,iX,iLine)
            Volume = volume_node(iLine, iX)
            Sum_VI(2:nSum,iLineAll) = Sum_VI(2:nSum,iLineAll) + Volume* &
                 [Distr, abs(Distr - Exact), Exact, (Distr - Exact)**2, &
                 Exact**2]
            Max_I(1) = max(Max_I(1), abs(Distr - Exact))
            Max_I(2) = max(Max_I(2), Exact)
            Max_I(3) = max(Max_I(3), &
                 maxval(abs(Distribution_CB(1:nP,1,iX,iLine) - Distr)))
            DistrMin = min(DistrMin, &
                 minval(Distribution_CB(1:nP,1:nMu,iX,iLine)))
         end do
      end do
      ! Collect the sums of all lines and add them up in the same order on
      ! every processor; collect the extrema and the lines in use
      call MPI_allreduce(MPI_IN_PLACE, Sum_VI, nSum*nLineAll, MPI_REAL, &
           MPI_SUM, iComm, iError)
      Sum_V = 0.0
      do iLineAll = 1, nLineAll
         Sum_V = Sum_V + Sum_VI(:,iLineAll)
      end do
      call MPI_allreduce(MPI_IN_PLACE, Max_I, 3, MPI_REAL, MPI_MAX, &
           iComm, iError)
      call MPI_allreduce(MPI_IN_PLACE, DistrMin, 1, MPI_REAL, MPI_MIN, &
           iComm, iError)
      nLineUsed = count(Used_B(1:nLine))
      call MPI_allreduce(MPI_IN_PLACE, nLineUsed, 1, MPI_INTEGER, MPI_SUM, &
           iComm, iError)
      ! Relative errors as in test_R2Sphere3d (weighted by the volume of
      ! the nodes) and the relative change of sum(f*V)
      Error_V(1) = Max_I(1)/Max_I(2)
      Error_V(2) = Sum_V(3)/Sum_V(4)
      Error_V(3) = sqrt(Sum_V(5)/Sum_V(6))
      Error_V(4) = (Sum_V(2) - Sum_V(1))/Sum_V(1)

      if(iProc == 0)then
         ! One row of the table of checks (header in test_Line2Sphere3d):
         ! the case, nNode, the coupled steps, the wall time of the steps,
         ! the lines in use, sum(fExact*V)/sum(f(t0)*V) - 1, min f and
         ! max |f(p) - f(p=1)|
         write(*,'(1x,a,2i6,f8.3,i6,es15.3,es12.3e3,es14.3e3)') &
              NameCase_I(iCase), nNode, nStep, CpuTime, nLineUsed, &
              Sum_V(4)/Sum_V(1) - 1, DistrMin, Max_I(3)
         ! The lines below are commented out: the row above has their
         ! values, or they are redundant.
         ! The case and nNode are in the row, the base lines (nThetaLine x
         ! nPhiLine, refined 2x) and nProc in the setup (test_Line2Sphere3d)
         ! write(*,'(/,a,i0,a,i0,a,i0,a,i0)') NameSub//' ('// &
         !      trim(NameCase_I(iCase))//'): lines nNode x nTheta x nPhi = ', &
         !      nNode, ' x ', nThetaLine, ' x ', nPhiLine, ', nProc = ', nProc
         ! The same for all cases: the setup gives the base mesh (refined 2x)
         ! write(*,'(a,i0,a,i0,a,i0)') ' mesh nR x nTheta x nPhi = ', &
         !      nRPerp, ' x ', nThetaPerp, ' x ', nPhiPerp
         ! The same for both sets of lines: the setup has them (the table
         ! of the cases, and TauR, TauA at t0 in the f(t0) line)
         ! write(*,'(a,f7.4,a,f7.4,a)') ' line Dpara = ', &
         !      (1 - DPerpRatio)*kDpara, '*r**2, Dperp = ', &
         !      DPerpRatio*kDpara, '*r**2 [m^2/s], r in m'
         ! write(*,'(a,2f8.4,a,2f8.4)') ' TauR, TauA at t0 =', &
         !      kPara*tStart, kPerp*tStart, ', at t0 + tFinal =', TauR, TauA
         ! nStep is in the row and tFinal (Dt = tFinal/nStep) in the setup
         ! write(*,'(a,i0,a,f6.3,a)') ' coupled steps: ', nStep, &
         !      ' of Dt = ', Dt, ' s (one line solver step each)'
         ! In the row (wall[s], and lines: nLineAll = nThetaLine*nPhiLine)
         ! write(*,'(a,f9.3,a)') ' wall time of the steps    =', CpuTime, ' s'
         ! write(*,'(a,i0,a,i0)') ' lines in use              = ', &
         !      nLineUsed, ' of ', nLineAll
         ! Absolute sums: the table of errors has their relative change
         ! write(*,'(a,2es22.14)') ' sum(f*V) [Rs^3] at t0, t0 + tFinal =', &
         !      Sum_V(1:2)
         ! The same as the dsum(fV)/sum column of the table of errors
         ! write(*,'(a,es11.3)') ' relative change of sum(f*V)   =', &
         !      Error_V(4)
         ! In the row (sum(fExV)/M0-1, min f, max|f-f(p=1)|)
         ! write(*,'(a,es11.3)') ' sum(fExact*V)/sum(f(t0)*V) - 1 =', &
         !      Sum_V(4)/Sum_V(1) - 1
         ! write(*,'(a,es12.3e3)') ' min f                         =', &
         !      DistrMin
         ! write(*,'(a,es12.3e3)') ' max |f(p) - f(p=1)|, p=1..nP  =', &
         !      Max_I(3)
         ! The same as the table of errors (test_Line2Sphere3d)
         ! write(*,'(a,3es11.3)') ' relative errors: max, L1, L2  =', &
         !      Error_V(1:3)
      end if

      if(iCase == 1)then
         ! Plots of f, exact f and error: phi-theta map at the node nearest
         ! to rc, along the line nearest to the source (its lat-lon cell
         ! contains the source direction) and along theta at its phi,
         ! collected (summed with zeros) on all processors and written by
         ! processor 0
         iLat0 = int(ThetaSource/180.0*nThetaLine) + 1
         iLon0 = int(PhiSource/360.0*nPhiLine) + 1
         iLineSource = (iLat0 - 1)*nPhiLine + iLon0
         iX0 = nint(log(rSource/RMin)/dLogNode) + 1
         allocate(Plot_VII(3,nPhiLine,nThetaLine), PlotR_VI(3,nNode), &
              PlotTheta_VI(3,nThetaLine))
         Plot_VII = 0.0; PlotR_VI = 0.0
         do iLine = 1, nLine
            iLineAll = iLineAll0 + iLine
            call line_angles(iLine, iLat, iLon, Theta, Phi)
            Distr = Distribution_CB(1,1,iX0,iLine)
            Exact = exact_line(State_VIB(R_,iX0,iLine), Theta, Phi, TauR, &
                 TauA)
            Plot_VII(:,iLon,iLat) = [Distr, Exact, Distr - Exact]
            if(iLineAll /= iLineSource) CYCLE
            do iX = 1, nNode
               Distr = Distribution_CB(1,1,iX,iLine)
               Exact = exact_line(State_VIB(R_,iX,iLine), Theta, Phi, &
                    TauR, TauA)
               PlotR_VI(:,iX) = [Distr, Exact, Distr - Exact]
            end do
         end do
         call MPI_allreduce(MPI_IN_PLACE, Plot_VII, size(Plot_VII), &
              MPI_REAL, MPI_SUM, iComm, iError)
         call MPI_allreduce(MPI_IN_PLACE, PlotR_VI, size(PlotR_VI), &
              MPI_REAL, MPI_SUM, iComm, iError)
         PlotTheta_VI = Plot_VII(:,iLon0,:)
         if(iProc == 0)then
            write(NameFile,'(a,i0,a,i0)') 'test_Line2Sphere3d_', &
                 nThetaLine, 'x', nPhiLine
            call save_plot_file(NameFile=trim(NameFile)//'.out', &
                 TypeFileIn='ascii', TimeIn=tStart+tFinal, nStepIn=nStep, &
                 StringHeaderIn='test_Line2Sphere3d: f, exact f and '// &
                 'f - exact f at the nodes nearest to rc, angles in deg', &
                 NameVarIn='Phi Theta f fExact Error rNode kPara kPerp', &
                 ParamIn_I=[RMin*exp((iX0 - 1)*dLogNode), kPara, kPerp], &
                 Coord1In_I=[((i - 0.5)*360.0/nPhiLine, i = 1, nPhiLine)], &
                 Coord2In_I=[((i - 0.5)*180.0/nThetaLine, i = 1, &
                 nThetaLine)], VarIn_VII=Plot_VII)
            call save_plot_file(NameFile=trim(NameFile)//'_r.out', &
                 TypeFileIn='ascii', TimeIn=tStart+tFinal, nStepIn=nStep, &
                 StringHeaderIn='test_Line2Sphere3d: f, exact f and '// &
                 'f - exact f along the line nearest to the source, r in Rs', &
                 NameVarIn='r f fExact Error Theta Phi kPara kPerp', &
                 ParamIn_I=[(iLat0 - 0.5)*180.0/nThetaLine, &
                 (iLon0 - 0.5)*360.0/nPhiLine, kPara, kPerp], &
                 Coord1In_I=[(RMin*exp((i - 1)*dLogNode), i = 1, nNode)], &
                 VarIn_VI=PlotR_VI)
            call save_plot_file(NameFile=trim(NameFile)//'_theta.out', &
                 TypeFileIn='ascii', TimeIn=tStart+tFinal, nStepIn=nStep, &
                 StringHeaderIn='test_Line2Sphere3d: f, exact f and '// &
                 'f - exact f along theta at the node nearest to rc, '// &
                 'theta in deg', &
                 NameVarIn='Theta f fExact Error rNode Phi kPara kPerp', &
                 ParamIn_I=[RMin*exp((iX0 - 1)*dLogNode), &
                 (iLon0 - 0.5)*360.0/nPhiLine, kPara, kPerp], &
                 Coord1In_I=[((i - 0.5)*180.0/nThetaLine, i = 1, &
                 nThetaLine)], VarIn_VI=PlotTheta_VI)
         end if
         deallocate(Plot_VII, PlotR_VI, PlotTheta_VI)
      end if

      deallocate(Sum_VI, DtLine_CB, DsSi_I, Mass_C, nSi_I, BSi_I, rSi_I)

    end subroutine run_Line2Sphere3d
    !==========================================================================
    subroutine set_lines

      ! Lines of this processor (nNode nodes, nThetaLine x nPhiLine lines),
      ! the background along them and f at t0 at all momenta. The lines are
      ! split between the processors as in init of SP_ModGrid; the line
      ! arrays get the shapes of a run (see the memory note above).

      integer:: iLine, iX, iLat, iLon
      ! Angles of the line [rad], its direction, r of a node [Rs]
      real:: Theta, Phi, Dir_D(3), r
      character(len=*), parameter:: NameSub = 'set_lines'
      !------------------------------------------------------------------------
      ! Lines of this processor: iLineAll0+1 to iLineAll0+nLine
      nLat = nThetaLine; nLon = nPhiLine; nLineAll = nLat*nLon
      if(nLineAll < nProc) call CON_stop(NameSub//': nProc > nLineAll')
      iLineAll0 = (iProc*nLineAll)/nProc
      nLine = ((iProc + 1)*nLineAll)/nProc - iLineAll0

      ! Line arrays with nVertexMax nodes per line, as in a run
      if(IsLineAllocated) deallocate(MhData_VIB, State_VIB, nVertex_B, &
           Used_B)
      allocate(MhData_VIB(LagrID_:nMhData,nVertexMax,nLine), &
           State_VIB(nMhData+1:nVar,nVertexMax,nLine), nVertex_B(nLine), &
           Used_B(nLine))
      IsLineAllocated = .true.
      if(allocated(Distribution_CB)) deallocate(Distribution_CB)
      allocate(Distribution_CB(0:nP+1,nMu,nVertexMax,nLine))
      if(allocated(DPerp_CB)) deallocate(DPerp_CB)
      allocate(DPerp_CB(nP,nVertexMax,nLine))
      nVertex_B = nNode
      Used_B = .true.

      ! Background: radial lines, nodes uniform in ln(r), B ~ 1/r**2, and
      ! the other MHD data zero as in init_stand_alone
      dLogNode = log(RMax/RMin)/(nNode - 1)
      do iLine = 1, nLine
         call line_angles(iLine, iLat, iLon, Theta, Phi)
         Dir_D = [sin(Theta)*cos(Phi), sin(Theta)*sin(Phi), cos(Theta)]
         do iX = 1, nNode
            r = RMin*exp((iX - 1)*dLogNode)
            MhData_VIB(:,iX,iLine) = 0.0
            MhData_VIB(LagrID_,iX,iLine) = real(iX)
            MhData_VIB(X_:Z_,iX,iLine) = r*Dir_D
            MhData_VIB(Rho_,iX,iLine) = Rho0*(RMin/r)**2
            MhData_VIB(Bx_:Bz_,iX,iLine) = B0*(RMin/r)**2*Dir_D
            State_VIB(:,iX,iLine) = 0.0
         end do
      end do
      ! B, distances, r etc. along the lines (State_VIB), as in a run
      call get_other_state_var

      ! f at t0: the exact f at all momenta, at least DistrFloor (and
      ! Dperp, set each step)
      do iLine = 1, nLine
         call line_angles(iLine, iLat, iLon, Theta, Phi)
         do iX = 1, nNode
            Distribution_CB(:,:,iX,iLine) = max(DistrFloor, exact_line( &
                 State_VIB(R_,iX,iLine), Theta, Phi, kPara*tStart, &
                 kPerp*tStart))
            DPerp_CB(:,iX,iLine) = 0.0
         end do
      end do

    end subroutine set_lines
    !==========================================================================
    subroutine line_angles(iLine, iLat, iLon, Theta, Phi)

      ! Row iLat, column iLon, colatitude Theta and longitude Phi [rad] of
      ! the line iLine of this processor: the center of its lat-lon cell

      integer, intent(in) :: iLine
      integer, intent(out):: iLat, iLon
      real,    intent(out):: Theta, Phi
      !------------------------------------------------------------------------
      call ilineall_to_lon_lat(iLineAll0 + iLine, iLon, iLat)
      Theta = (iLat - 0.5)*cPi/nThetaLine
      Phi   = (iLon - 0.5)*2*cPi/nPhiLine

    end subroutine line_angles
    !==========================================================================
    real function volume_node(iLine, iX)

      ! Volume [Rs^3] of node iX of line iLine of this processor:
      ! r**3*dln(r)*dOmega (half at the line ends), dOmega of the lat-lon
      ! cell of the line

      integer, intent(in):: iLine, iX

      integer:: iLat, iLon
      real:: Theta, Phi, dTheta
      !------------------------------------------------------------------------
      call line_angles(iLine, iLat, iLon, Theta, Phi)
      dTheta = cPi/nThetaLine
      volume_node = State_VIB(R_,iX,iLine)**3*dLogNode*2*cPi/nPhiLine &
           *(cos(Theta - 0.5*dTheta) - cos(Theta + 0.5*dTheta))
      if(iX == 1 .or. iX == nVertex_B(iLine)) volume_node = 0.5*volume_node

    end function volume_node
    !==========================================================================
    real function exact_line(r, Theta, Phi, TauR, TauA)

      ! Exact f at the point (r [Rs], Theta, Phi [rad]) for the radial and
      ! angular k*t accumulated since the point source, TauR and TauA:
      ! Mass/rc**3*G(ln(r); TauR)*K(gamma; TauA), the G and K of
      ! test_R2Sphere3d

      real, intent(in):: r, Theta, Phi, TauR, TauA

      ! Source angles [rad], cos(gamma), Legendre polynomials P_(l-1), P_l
      ! and P_(l+1), exp(-l(l+1)*TauA), and K(gamma) times 4*pi
      real:: Theta0, Phi0, CosGamma, LegendreMinus, Legendre, &
           LegendrePlus, Decay, SumK
      integer:: l
      !------------------------------------------------------------------------
      Theta0 = ThetaSource*cDegToRad
      Phi0   = PhiSource*cDegToRad
      ! Angle gamma between the point and the source point (from r = 0)
      CosGamma = max(-1.0, min(1.0, cos(Theta)*cos(Theta0) &
           + sin(Theta)*sin(Theta0)*cos(Phi - Phi0)))
      ! Sum of (2l+1)*P_l(cos gamma)*exp(-l(l+1)*TauA) over l, with P_0 = 1,
      ! P_1 = x and (l+1)*P_(l+1) = (2l+1)*x*P_l - l*P_(l-1), up to the
      ! first l with exp(-l(l+1)*TauA) < 1e-16
      SumK = 0.0
      LegendreMinus = 0.0
      Legendre = 1.0
      l = 0
      do
         Decay = exp(-l*(l+1)*TauA)
         if(Decay < 1e-16) EXIT
         SumK = SumK + (2*l+1)*Legendre*Decay
         LegendrePlus = ((2*l+1)*CosGamma*Legendre - l*LegendreMinus)/(l+1)
         LegendreMinus = Legendre
         Legendre = LegendrePlus
         l = l + 1
      end do
      exact_line = Mass/rSource**3 &
           *exp(-(log(r/rSource) + 3*TauR)**2/(4*TauR))/sqrt(4*cPi*TauR) &
           *SumK/(4*cPi)

    end function exact_line
    !==========================================================================
    subroutine set_param_line(nCell_D, DPerpRatioIn)

      ! Read the parameters with the read_param routines of the modules, as
      ! a run does: #DIFFUSIONPERP (the mesh as in test_R2Sphere3d and
      ! DPerpRatioIn), #MOMENTUMBC, #LOWERENDBC and #UPPERENDBC (float),
      ! #FLUXINITIAL (a negligible background f). The commands are written
      ! into a file, read and broadcast by ModReadParam; the file is
      ! deleted after reading. Then init of SP_ModDistribution (momentum
      ! grid, first call only) and init of SP_ModPerpDiffusion (mesh, face
      ! areas and volumes) are called with nLine = 0 (see the memory note).

      integer, intent(in):: nCell_D(3)
      real,    intent(in):: DPerpRatioIn

      character(len=*), parameter:: NameParamFile = 'test_line_PARAM.in'
      character(len=*), parameter:: StringTab = char(9)//char(9)//char(9)
      character(len=lStringLine):: NameCommand
      character(len=*), parameter:: NameSub = 'set_param_line'
      !------------------------------------------------------------------------
      if(iProc == 0)then
         call open_file(File=NameParamFile, NameCaller=NameSub)
         write(UnitTmp_,'(a)')      '#DIFFUSIONPERP'
         write(UnitTmp_,'(a)')      'T'//StringTab//'UseDiffusionPerp'
         write(UnitTmp_,'(a)')      'UseDPerpRatio'//StringTab//'TypeDPerp'
         write(UnitTmp_,'(es22.15,a)') DPerpRatioIn, StringTab//'DPerpRatio'
         write(UnitTmp_,'(i0,a)')   nCell_D(1), StringTab//'nRPerp'
         write(UnitTmp_,'(f0.2,a)') RMin, StringTab//'RMinPerp'
         write(UnitTmp_,'(f0.2,a)') RMax, StringTab//'RMaxPerp'
         write(UnitTmp_,'(a)')      'Exp'//StringTab//'ScaleRPerp'
         write(UnitTmp_,'(i0,a)')   nCell_D(2), StringTab//'nThetaPerp'
         write(UnitTmp_,'(i0,a)')   nCell_D(3), StringTab//'nPhiPerp'
         ! F: T prints the sub-step and positivity diagnostics of
         ! solver_fvm_diffuse3d and diffuseperp_distribution at every
         ! call; set T to debug
         write(UnitTmp_,'(a)')      'F'//StringTab//'DoTestDPerp'
         write(UnitTmp_,'(/,a)')    '#MOMENTUMBC'
         write(UnitTmp_,'(a)')      'float'//StringTab//'TypeMomentumMinBc'
         write(UnitTmp_,'(a)')      'float'//StringTab//'TypeMomentumMaxBc'
         write(UnitTmp_,'(/,a)')    '#LOWERENDBC'
         write(UnitTmp_,'(a)')      'float'//StringTab//'TypeLowerEndBc'
         write(UnitTmp_,'(/,a)')    '#UPPERENDBC'
         write(UnitTmp_,'(a)')      'float'//StringTab//'TypeUpperEndBc'
         write(UnitTmp_,'(/,a)')    '#FLUXINITIAL'
         write(UnitTmp_,'(a)')      '1e-200'//StringTab//'FluxInitIo'
         call close_file(NameCaller=NameSub)
      end if
      ! Quiet: no 'read and broadcast' line, and no echo of the
      ! parameters (read_file turns the echo on, on processor 0)
      call read_file(NameParamFile, iComm, IsVerbose=.false.)
      call read_echo_set(.false.)
      if(iProc == 0)then
         call open_file(File=NameParamFile, Status='old', NameCaller=NameSub)
         call close_file(Status='delete', NameCaller=NameSub)
      end if
      call read_init(nLineIn=n_line_read())
      do
         if(.not.read_line()) EXIT
         if(.not.read_command(NameCommand)) CYCLE
         select case(NameCommand)
         case('#DIFFUSIONPERP')
            call read_param_dperp(NameCommand)
         case('#MOMENTUMBC', '#LOWERENDBC', '#UPPERENDBC')
            call read_param_bc(NameCommand)
         case('#FLUXINITIAL')
            call read_param_dist(NameCommand)
         case default
            call CON_stop(NameSub//': unknown command '//NameCommand)
         end select
      end do
      ! Momentum grid and mesh, with no lines
      nLine = 0
      call init_dist
      call init_dperp
      ! read_param changes an odd nTheta, and nR if it is less than nProc
      if(nRPerp /= nCell_D(1) .or. nThetaPerp /= nCell_D(2) .or. &
           nPhiPerp /= nCell_D(3)) call CON_stop(NameSub// &
           ': mesh changed by read_param (odd nTheta or nR < nProc?)')

    end subroutine set_param_line
    !==========================================================================
  end subroutine test_Line2Sphere3d
  !============================================================================
end module ModTestDperp
!==============================================================================
program test_transport

  use ModTestDperp, ONLY: test_UniSphere3d, test_R2Sphere3d, &
       test_Line2Sphere3d
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
  case('', 'dperp')
     ! All tests of ModTestDperp
     call test_UniSphere3d
     call test_R2Sphere3d
     call test_Line2Sphere3d
  case('unisphere3d')
     call test_UniSphere3d
  case('r2sphere3d')
     call test_R2Sphere3d
  case('line2sphere3d')
     call test_Line2Sphere3d
  case default
     call CON_stop('test_transport: unknown test '//trim(NameTest))
  end select

  call MPI_finalize(iError)

end program test_transport
!==============================================================================
