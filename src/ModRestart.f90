!  Copyright (C) 2002 Regents of the University of Michigan,
!  portions used with permission
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module SP_ModRestart

  ! This module contains methods for writing output files
  use SP_ModSize,   ONLY: nVertexMax
  use SP_ModGrid,   ONLY: iblock_to_lon_lat, LagrID_, Z_,&
       nLine, MhData_VIB, iShock_IB,  Used_B,&
       FootPoint_VB, nVertex_B, nShockParam, &
       nLon, nLat
  use SP_ModDistribution, ONLY: Distribution_IIB
  use SP_ModTime,   ONLY: SPTime, iIter
  use SP_ModUnit,   ONLY: NameEnergyUnit
  use ModPlotFile,  ONLY: save_plot_file, read_plot_file
  use ModUtilities, ONLY: open_file, close_file, CON_stop
  use ModIoUnit,    ONLY: UnitTmp_

  implicit none

  SAVE

  private ! except

  ! Public members
  public:: save_restart, read_restart, stand_alone_save_restart, read_param
  public:: NameRestartInDir, NameRestartOutDir

  ! the restart directory
  character (len=100) :: NameRestartOutDir="SP/restartOUT/"
  character (len=100) :: NameRestartInDir ="SP/restartIN/"
  ! name of the header file
  character (len=100) :: NameHeaderFile   ="restart.H"

  ! Whether and when to save the restart file (STAND ALONE ONLY)
  logical:: DoSaveRestart=.false.
  integer:: DnSaveRestart=-1,  nIterSinceRestart = 0
  real   :: DtSaveRestart=-1.0, TimeSinceRestart = 0.0

contains
  !============================================================================
  subroutine read_param
    use ModReadParam, ONLY: read_var

    character(len=*), parameter:: NameSub = 'read_param'
    !--------------------------------------------------------------------------
    call read_var("DoSaveRestart",DoSaveRestart)
    call read_var("DnSaveRestart",DnSaveRestart)
    call read_var("DtSaveRestart",DtSaveRestart)
    if(DtSaveRestart < 0. .and. DnSaveRestart < 0)&
         call CON_stop(NameSub//': incorrectly set ')

  end subroutine read_param
  !============================================================================
  subroutine stand_alone_save_restart(Dt)

    real, intent(in) :: Dt
    !--------------------------------------------------------------------------
    if(.not.DoSaveRestart) RETURN
    nIterSinceRestart = nIterSinceRestart + 1
    TimeSinceRestart  = TimeSinceRestart  + Dt
    if(  DtSaveRestart > 0..and.  TimeSinceRestart >= DtSaveRestart .or. &
         DnSaveRestart > 0 .and. nIterSinceRestart == DnSaveRestart)then
       call save_restart
       if(DtSaveRestart > 0.)then
          TimeSinceRestart  = modulo(TimeSinceRestart, DtSaveRestart)
       else
          TimeSinceRestart = 0.
       end if
       nIterSinceRestart = 0
    end if

  end subroutine stand_alone_save_restart
  !============================================================================
  subroutine save_restart

    use ModIoUnit,     ONLY: UnitTmp_
    use ModUtilities,  ONLY: open_file, close_file
    ! write the restart data

    ! name of the output file
    character(len=100):: NameFile
    ! loop variable
    integer:: iLine
    ! indexes of corresponding node, latitude and longitude
    integer:: iLat, iLon
    character(len=*), parameter:: NameSub = 'save_restart'
    !--------------------------------------------------------------------------

    call write_restart_header

    do iLine = 1, nLine
       if(.not.Used_B(iLine))CYCLE
       call iblock_to_lon_lat(iLine, iLon, iLat)

       ! set the file name
       write(NameFile,'(a,i3.3,a,i3.3,a)') &
            trim(NameRestartOutDir)//'data_',iLon,'_',iLat,&
            '.rst'
       call open_file(file=NameFile, form='UNFORMATTED',&
            NameCaller=NameSub)
       write(UnitTmp_)real(nVertex_B(iLine)),&
            real(iShock_IB(:, iLine))
       write(UnitTmp_)&
            FootPoint_VB(:, iLine),&
            MhData_VIB(LagrID_:Z_,1:nVertex_B(iLine), iLine),&
            Distribution_IIB(:,1:nVertex_B(iLine), iLine)
       call close_file
    end do

  end subroutine save_restart
  !============================================================================
  subroutine read_restart

    ! read the restart data

    ! name of the input file
    character(len=100):: NameFile
    ! loop variables
    integer:: iLine
    ! indexes of corresponding node, latitude and longitude
    integer:: iLat, iLon
    real   :: Aux, Aux_I(nShockParam) ! For reading integers
    integer:: iError
    character(len=*), parameter:: NameSub = 'read_restart'
    !--------------------------------------------------------------------------

    do iLine = 1, nLine
       call iBlock_to_lon_lat(iLine, iLon, iLat)
       ! set the file name
       write(NameFile,'(a,i3.3,a,i3.3,a)') &
            trim(NameRestartInDir)//'data_',iLon,'_',iLat,&
            '.rst'
       inquire(file=NameFile,exist=Used_B(iLine))
       if(.not.Used_B(iLine))then
          write(*,'(a)')NameSub//': the restart file '//NameFile//' lacks'
          write(*,'(a)')NameSub//': line is marked as unused'
          nVertex_B(iLine) = 0
          CYCLE
       end if
       call open_file(file=NameFile,  status='old',&
            form='UNFORMATTED', NameCaller=NameSub)
       read(UnitTmp_,iostat = iError)Aux, Aux_I
       if(iError>0)then
          write(*,*)'Error in reading nPoint in line=', iLine
          call close_file
          call CON_stop('Run stops')
       end if
       ! process buffer
       nVertex_B(iLine) = nint(Aux)
       ! general parameters
       iShock_IB(:, iLine) = nint(Aux_I)
       read(UnitTmp_, iostat = iError) &
            FootPoint_VB(:, iLine),&
            MhData_VIB(LagrID_:Z_,1:nVertex_B(iLine), iLine),&
            Distribution_IIB(:,1:nVertex_B(iLine), iLine)
       call close_file
    end do

  end subroutine read_restart
  !============================================================================
  subroutine write_restart_header

    use SP_ModPlot,         ONLY: nTag
    use SP_ModProc,         ONLY: iProc
    use ModUtilities,       ONLY: cTab
    use SP_ModDistribution, ONLY: nP, EnergyInjIo, EnergyMaxIo
    ! full name of the header file
    character(len=100):: NameFile

    character(len=*), parameter:: NameSub = 'write_restart_header'
    !--------------------------------------------------------------------------
    if (iProc/=0) RETURN
    NameFile = trim(NameRestartOutDir)//trim(NameHeaderFile)

    call open_file(file=NameFile, NameCaller=NameSub)
    write(UnitTmp_,*)
    write(UnitTmp_,'(a)')'#RESTART'
    write(UnitTmp_,'(a)')'T'//cTab//cTab//cTab//'DoRestart'
    write(UnitTmp_,*)
    write(UnitTmp_,'(a)')'#CHECKGRIDSIZE'
    write(UnitTmp_,'(i8,a)') nVertexMax,cTab//cTab//'nVertexMax'
    write(UnitTmp_,'(i8,a)') nLon,     cTab//cTab//'nLon'
    write(UnitTmp_,'(i8,a)') nLat,     cTab//cTab//'nLat'
    write(UnitTmp_,*)
    write(UnitTmp_,'(a)')'#NSTEP'
    write(UnitTmp_,'(i8,a)')iIter,cTab//cTab//'nStep'
    write(UnitTmp_,*)
    write(UnitTmp_,'(a)')'#TIMESIMULATION'
    write(UnitTmp_,'(es22.15,a)')SPTime,cTab//cTab//'tSimulation'
    write(UnitTmp_,*)
    write(UnitTmp_,'(a)')'#NTAG'
    write(UnitTmp_,'(i8,a)')nTag,cTab//cTab//'nTag'
    write(UnitTmp_,*)
    write(UnitTmp_,'(a)')'#MOMENTUMGRID'
    write(UnitTmp_,'(es22.15,a)')EnergyInjIo,cTab//cTab//'EnergyMin'
    write(UnitTmp_,'(es22.15,a)')EnergyMaxIo,cTab//cTab//'EnergyMax'
    write(UnitTmp_,'(i8,a)')nP,cTab//cTab//'nP'
    write(UnitTmp_,'(a8,a)')NameEnergyUnit,cTab//cTab//'NameEnergyUnit'
    write(UnitTMP_,'(a)')'#END'
    write(UnitTmp_,*)
    call close_file

  end subroutine write_restart_header
  !============================================================================
end module SP_ModRestart
!==============================================================================
