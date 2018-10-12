!  Copyright (C) 2002 Regents of the University of Michigan, 
!  portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!=================================================================!
module SP_ModPlot
  !\
  ! Methods for saving plots
  !/
  use SP_ModSize, ONLY: nParticleMax
  use SP_ModGrid, ONLY: get_node_indexes, nVar, nMHData, nBlock,   &
       State_VIB, iShock_IB, iNode_B, nParticle_B, Shock_, X_, Z_, &
       R_, NameVar_V, TypeCoordSystem, LagrID_, nNode
  use SP_ModDistribution, ONLY: nP, KinEnergySI_I, MomentumSI_I,      &
       Distribution_IIB,                                           &
       Flux_VIB, Flux0_, FluxMax_, NameFluxChannel_I, nFluxChannel
  use SP_ModTime, ONLY: SPTime, iIter, StartTime, StartTimeJulian
  use SP_ModProc, ONLY: iProc
  use ModPlotFile, ONLY: save_plot_file, read_plot_file
  use ModUtilities, ONLY: open_file, close_file, remove_file
  use ModIoUnit, ONLY: UnitTmp_
  implicit none
  SAVE
  private ! except
  !\
  !Public members
  public:: init         ! Initialize module parameters
  public:: read_param   ! Read module parameters
  public:: save_plot_all! Save output (plot) files 
  public:: finalize     ! Save final list  
  ! the output directory
  character(len=*), public, parameter :: NamePlotDir="SP/IO2/"
  !/
  !\
  ! Number of plot file types
  integer:: nFileOut = 0
  ! Types of output files in terms of output dataa
  integer, parameter:: &
       ! Background mhd data
       MH1D_      = 0, & ! along each line
       MH2D_      = 1, & ! at given radius as Lon-Lat plot
       MHTime_    = 2, & ! at given radius as time series for lines 
       ! Distribution
       Distr1D_   = 3    ! along each line
  ! Momentum or energy axis to use for 2D plots 
  integer, parameter:: & 
       Momentum_= 1,   &
       Energy_  = 2
  ! Plot types for distribution function
  integer, parameter:: &
       CDF_ = 1,       &
       DEF_ = 2
  !/
  !\
  type TypePlotFile
     !\
     ! General information
     !/
     ! kind of data printed to a file
     integer:: iKindData   !=MH1D_or MH2D_ or MHTime_ or Distr1D_ 
     ! file name extension
     character(len=4 ):: NameFileExtension  !.out, .tec etc
     character(len=20):: TypeFile           !tec, idl, ascii, real8
     ! whether it is the first call
     ! USED ONLY IN write_mh_time  FOR NOW!!!
     logical:: IsFirstCall
     ! names of variables to be written
     character(len=300):: NameVarPlot
     ! output buffer
     real, pointer:: Buffer_II(:,:)
     !\
     ! MH data
     !/
     ! variables from the state vector to be written
     logical:: DoPlot_V(X_:nVar)
     ! variables from the flux vector to be written
     logical, allocatable :: DoPlotFlux_I(:)
     ! total numbers of variables to be written
     integer:: nVarPlot, nFluxPlot
     ! their indices in the state vectors
     integer, pointer:: iVarPlot_V(:), iFluxPlot_V(:)
     !\
     ! Distribution
     !/
     ! Momentum or energy axis to use for 2D plots 
     integer:: iScale      ! =Momentum_ or Energy_
     ! type out output (CDF or differential energy flow)
     integer:: iTypeDistr  ! =CDF_ or EDF_
     !\
     ! Data on the sphere
     !/
     ! radius of the sphere the data to be written at
     real:: Radius
  end type TypePlotFile
  !/
  !\
  ! All plot files
  type(TypePlotFile), allocatable:: File_I(:)
  !/
  !\
  ! Arrays used to visualize the distribution function
  real, dimension(0:nP+1) :: Log10MomentumSI_I, Log10KinEnergySI_I,      &
       DMomentumOverDEnergy_I
  !/
  !\
  ! auxilary array, used to write data on a sphere
  integer, allocatable:: iNodeIndex_I(:)
  !/
  !\
  ! info for MH1D header and tag list
  logical:: DoWriteHeader = .false.
  character(len=*), parameter, public :: NameMHData = "MH_data"
  ! name of the header file
  character(len=*), parameter :: NameHeaderFile = NameMHData//'.H'
  ! name of the tag list file
  character(len=*), parameter :: NameTagFile  = NameMHData//'.lst'
  ! number of different output file tags
  integer,  public            :: nTag = 0 
  !/
  !\ 
  character(len=20)           :: TypeMHDataFile 
  !/
  !\
  !Fromat for saving time tag. If .true. the time tag format is
  !YYYYMMDDHHMMSS 
  logical, public:: UseDateTime = .false.
  !/
  !\
  ! If DoSaveInitial=.false.,the initial files are not saved 
  logical :: DoSaveInitial = .true.
  !/
  logical :: DoInit        = .true.
contains
  subroutine init
    use ModConst,           ONLY: cLightSpeed
    use SP_ModDistribution, ONLY: EnergySI_I
    ! storage for existing tags (possible during restart
    character(len=50),allocatable:: StringTag_I(:)
    ! full tag file name
    character(len=100)::NameFile
    ! loop variable
    integer:: iTag
    
    character(len=*),parameter:: NameSub='SP:ModPlot:init'
    !--------------------------------------------------
    if(.not.DoInit) RETURN
    DoInit = .false.
    ! Array for plotting distribution function
    Log10MomentumSI_I      = log10(MomentumSI_I) 
    Log10KinEnergySI_I     = log10(KinEnergySI_I)
    DMomentumOverDEnergy_I = EnergySI_I/(MomentumSI_I*cLightSpeed**2)
    
    !\
    ! Reset/trim NameTagFile if nTag==0/nTag>0; the latter happens at restart.
    !
    ! During the run new tags are continuously appended to NameTagFile,
    ! however, only when the run is succesfully finalized the list of tags is
    ! considered to be valid, i.e. new #NTAG is written to the header file.
    !
    ! If previous run hasn't been properly finalized, NameTagFile may be
    ! inconsistent with nTag, therefore the file is trimmed according to nTag
    !/
    if(iProc/=0) RETURN ! done only by the root
    ! full file name
    NameFile = trim(NamePlotDir)//trim(NameTagFile)
    if(nTag>0)then 
       allocate(StringTag_I(nTag))
       call open_file(file=NameFile, status='old', NameCaller=NameSub)
       do iTag = 1, nTag
          read(UnitTmp_,'(a)') StringTag_I(iTag)
       end do
       call close_file
    end if
    call open_file(file=NameFile, status='replace', NameCaller=NameSub)
    if(nTag > 0)then
       do iTag = 1, nTag
          write(UnitTmp_,'(a)') StringTag_I(iTag)
       end do
       deallocate(StringTag_I)
    end if
    call close_file
  end subroutine init
  !===============================================================
  subroutine read_param(NameCommand)
    use ModUtilities, ONLY: split_string, lower_case
    use ModReadParam, ONLY: read_var
    character(len=*), intent(in):: NameCommand
    ! set parameters of output files: file format, kind of output etc.
    character(len=300):: StringPlot
    ! loop variables
    integer:: iFile, iNode, iVar
    integer:: nStringPlot
    character(len=20):: TypeFile, KindData, StringPlot_I(2*nVar)
    character(len=*), parameter :: NameSub='SP:set_write_param'
    !--------------------------------------------------------------------------
    select case(NameCommand)
    case("#SAVEPLOT")
       ! initialize auxilary arrays

       ! GOES by default, #SAVEPLOT must be after #FLUXCHANNEL
       if (.not. allocated(NameFluxChannel_I)) then
          if (nFluxChannel == 6) then
             allocate(NameFluxChannel_I(0:nFluxChannel+1))
             NameFluxChannel_I = (/'flux_total', 'flux_00005', 'flux_00010', &
                  'flux_00030', 'flux_00050', 'flux_00060', 'flux_00100', &
                  'eflux     '/)
          else
             call CON_stop(NameSub//' check nFluxChannel ')
          end if
       end if

       if (.not. allocated(iNodeIndex_I)) &
            allocate(iNodeIndex_I(nNode))

       do iNode = 1, nNode
          iNodeIndex_I(iNode) = iNode
       end do
       ! number of output files
       call read_var('nFileOut', nFileOut)
       ! check correctness
       if(nFileOut == 0) RETURN ! no output file requested
       if(nFileOut  < 0) call CON_stop(NameSub//': incorrect SAVEPLOT section')
       
       ! allocate the storage for file info
       allocate(File_I(nFileOut))
       
       ! read info about each file
       do iFile = 1, nFileOut
          ! reset and read the file info
          StringPlot = ''
          call read_var('StringPlot', StringPlot)
          
          ! make comparison case insensitive: convert strings to lower case
          call lower_case(StringPlot)
          
          ! put individual variables' and format names in separate array entries
          call split_string(StringPlot, StringPlot_I, nStringPlot)
          
          ! data kind is the first entry
          KindData = StringPlot_I(1)
          ! check whether set properly
          select case(KindData)
          case('mh1d')
             File_I(iFile) % iKindData = MH1D_
          case('mh2d')
             File_I(iFile) % iKindData = MH2D_
          case('mhtime')
             File_I(iFile) % iKindData = MHTime_
          case('distr1d')
             File_I(iFile) % iKindData = Distr1D_
          case default
             call CON_stop(NameSub//&
                  ": kind of data isn't properly set in PARAM.in")
          end select
          
          ! format of output is the last entry
          TypeFile = StringPlot_I(nStringPlot)
          ! check whether set properly
          select case(TypeFile)
          case('tec')
             File_I(iFile) % NameFileExtension='.dat'
             File_I(iFile) % TypeFile  ='tec'
          case('idl','ascii')
             File_I(iFile) % NameFileExtension='.out'
             File_I(iFile) % TypeFile  ='ascii'
          case('real4','real8')
             File_I(iFile) % NameFileExtension='.out'
             File_I(iFile) % TypeFile  = TypeFile
          case default
             call CON_stop(NameSub//&
                  ": output format isn't properly set in PARAM.in")
          end select
          
          ! reset variables' names
          File_I(iFile) % NameVarPlot = ''
          
          ! based on kind of data process the requested output
          select case(File_I(iFile) % iKindData)
          case(MH1D_)
             call process_mh
             ! prepare the output data container
             allocate(File_I(iFile) % Buffer_II(&
                  File_I(iFile)%nVarPlot + File_I(iFile)%nFluxPlot,&
                  1:nParticleMax))
             ! add particle index to variable names
             File_I(iFile) % NameVarPlot = &
                  NameVar_V(LagrID_)//' '//&
                  trim(File_I(iFile) % NameVarPlot)
             do iVar = LagrID_,Z_
                File_I(iFile)%NameVarPlot = &
                     trim(File_I(iFile)%NameVarPlot)//&
                     ' '//NameVar_V(iVar)  
             end do
             File_I(iFile) % NameVarPlot = &
                  trim(File_I(iFile) % NameVarPlot)//&
                  ' iShock RShock StartTime StartTimeJulian'
             TypeMHDataFile = File_I(iFile) % TypeFile
             DoWriteHeader = .true.
          case(MH2D_)
             call process_mh
             ! prepare the output data container
             allocate(File_I(iFile) % Buffer_II(&
                  File_I(iFile)%nVarPlot + File_I(iFile)%nFluxPlot,&
                  nNode))
             ! add line index to variable names
             File_I(iFile) % NameVarPlot = &
                  'LineIndex '//trim(File_I(iFile) % NameVarPlot)//&
                  ' StartTime StartTimeJulian'
             ! get radius
             call read_var('Radius [Rs]', File_I(iFile) % Radius)
          case(MHTime_)
             call process_mh
             ! prepare the output data container;
             ! note extra space reserved for time of the output
             allocate(File_I(iFile) % Buffer_II(&
                  1 + File_I(iFile)%nVarPlot + File_I(iFile)%nFluxPlot, 1))
             ! add time interval index to variable names
             File_I(iFile) % NameVarPlot = &
                  'Time '//trim(File_I(iFile) % NameVarPlot)//&
                  ' StartTime StartTimeJulian'
             ! get radius
             call read_var('Radius [Rs]', File_I(iFile) % Radius)
             ! reset indicator of the first call
             File_I(iFile) % IsFirstCall = .true.
          case(Distr1D_)
             call process_distr
             ! prepare the output data container
             allocate(File_I(iFile) % &
                  Buffer_II(nP,1:nParticleMax))
          end select
       end do
       !\
       ! Check consistency
       !/
       ! only 1 MH1D file can be requested
       if(count(File_I(1:nFileOut)%iKindData == MH1D_,1) > 1)&
            call CON_stop(NameSub//&
            ": only one MH1D output file can be requested")
    case("#USEDATETIME")
       call read_var('UseDateTime',UseDateTime)
    case('#SAVEINITIAL')
       call read_var('DoSaveInitial',DoSaveInitial)
    case('#NTAG')
       call read_var('nTag', nTag)
    case default
       call CON_stop('Unknown command '//NameCommand//' in '//&
            NameSub)
    end select

  contains
    !=============================================================
    subroutine process_mh
      ! process variables to plot
      ! NOTE: for iKindData == MH1D_ certain variables are always printed:
      !       Rho_, T_, Ux_:Uz_, Bx_:Bz_, Wave1_, Wave2_
      integer:: iVar, iVarPlot, iStringPlot
      character(len=10) ::  NameVarLowerCase
      !-----------------------------------------------

      if (allocated(File_I(iFile)%DoPlotFlux_I)) &
           deallocate(File_I(iFile)%DoPlotFlux_I)
      allocate(File_I(iFile)%DoPlotFlux_I(0:nFluxChannel+1))

      ! reset
      File_I(iFile) % DoPlot_V     = .false.
      File_I(iFile) % DoPlotFlux_I = .false.

      ! for MH1D_ minimal set of variables is printed
      if(File_I(iFile)%iKindData == MH1D_)then
         ! for MH1D_ minimal set of variables is printed
         File_I(iFile) % DoPlot_V(1:nMHData) = .true.
      else
         ! coordinates are always printed
         File_I(iFile) % DoPlot_V(X_:Z_) = .true.
      end if
      !\
      ! determine, which variables were requested to be in the output file
      !/
      do iStringPlot = 2, nStringPlot - 1
         ! if the string is flux, then save ALL the fluxes
         if (trim(StringPlot_I(iStringPlot)) == 'flux') then
            File_I(iFile)%DoPlotFlux_I = .true.
            CYCLE
         end if

         ! check names of individual variables
         do iVar = 1, nVar
            NameVarLowerCase = NameVar_V(iVar)
            call lower_case(NameVarLowerCase)
            if(StringPlot_I(iStringPlot) == NameVarLowerCase)&
                 File_I(iFile)%DoPlot_V(iVar) = .true.
            CYCLE
         end do

         do iVar = 0, nFluxChannel+1
            NameVarLowerCase = NameFluxChannel_I(iVar)
            call lower_case(NameVarLowerCase)
            if(StringPlot_I(iStringPlot) == NameVarLowerCase)&
                 File_I(iFile)%DoPlotFlux_I(iVar) = .true.
            CYCLE
         end do
      end do

      File_I(iFile)%nVarPlot  = count(File_I(iFile)%DoPlot_V(1:nVar))
      File_I(iFile)%nFluxPlot = count(File_I(iFile)%DoPlotFlux_I)
      ! indices in the state vector
      allocate(File_I(iFile)%iVarPlot_V(File_I(iFile)%nVarPlot))
      ! determine indices and names of variables
      iVarPlot = 0
      do iVar = 1, nVar
         if(.not.File_I(iFile)%DoPlot_V(iVar)) CYCLE
         File_I(iFile)%NameVarPlot = &
              trim(File_I(iFile)%NameVarPlot)//' '//&
              trim(NameVar_V(iVar))
         iVarPlot = iVarPlot + 1
         File_I(iFile)%iVarPlot_V(iVarPlot) = iVar
      end do
      if(File_I(iFile)%nFluxPlot > 0)then
         allocate(File_I(iFile)%iFluxPlot_V(File_I(iFile)%nFluxPlot))
         iVarPlot = 0
         do iVar = Flux0_, FluxMax_
             if(.not.File_I(iFile)%DoPlotFlux_I(iVar)) CYCLE
            File_I(iFile)%NameVarPlot = &
                 trim(File_I(iFile)%NameVarPlot)//' '//&
                 trim(NameFluxChannel_I(iVar))
            iVarPlot = iVarPlot + 1
            File_I(iFile) % iFluxPlot_V(iVarPlot) = iVar
         end do
      end if

      deallocate(File_I(iFile)%DoPlotFlux_I)
    end subroutine process_mh
    !=============================================================
    subroutine process_distr
      ! process output parameters for distribution output
      integer:: iStringPlot
      character(len=20):: NameVar, NameScale
      !-------------------------------------------------
      ! only 1 variable is printed
      File_I(iFile) % nVarPlot = 1
      ! reset string with variables' names and put defaults
      NameScale = 'Log10Momentum'
      NameVar   = 'Log10DiffEnergyFlux'
      File_I(iFile) % iScale     = Momentum_
      File_I(iFile) % iTypeDistr = DEF_
      do iStringPlot = 2, nStringPlot - 1
         ! may contain type of output scale (momentum/energy)
         ! or utput variable (canonical distr func/differential energy flux)
         select case(StringPlot_I(iStringPlot))
         case('momentum')
            File_I(iFile) % iScale = Momentum_
            NameScale = 'Log10Momentum'
         case('energy')
            File_I(iFile) % iScale = Energy_
            NameScale = 'Log10Energy'
         case('cdf')
            ! canonical distribution function
            File_I(iFile) % iTypeDistr = CDF_
            NameVar = 'Log10Distribution'
         case('def')
            ! differential energy flux
            File_I(iFile) % iTypeDistr = DEF_
            NameVar = 'Log10DiffEnergyFlux'
         end select
      end do
      ! form the name with variables' names'
      File_I(iFile) % NameVarPlot = &
           trim(NameScale)//' Distance '//trim(NameVar)
    end subroutine process_distr
  end subroutine read_param
  !===============================================================
  subroutine save_plot_all(IsInitialOutputIn)
    use SP_ModDistribution, ONLY: get_integral_flux
    ! write the output data
    logical, intent(in), optional:: IsInitialOutputIn

    ! loop variables
    integer:: iFile
    integer:: iKindData

    logical:: IsInitialOutput

    character(len=*), parameter:: NameSub = 'SP:write_output'
    !------------------------------------------------------------
    ! check whether this is a call for initial output
    if(present(IsInitialOutputIn))then
       IsInitialOutput = IsInitialOutputIn
    else
       IsInitialOutput = .false.
    end if
    if(nFileOut == 0) RETURN
    if(IsInitialOutput)then
       if(.not.DoSaveInitial)RETURN
    else
       call get_integral_flux
    end if
    do iFile = 1, nFileOut
       iKindData = File_I(iFile) % iKindData

       ! during initial call only background 1D data is printed
       if(IsInitialOutput .and. iKindData /= MH1D_)&
            iKindData = -1

       select case(iKindData)
       case(MH1D_)
          call write_mh_1d
       case(MH2D_)
          call write_mh_2d
       case(MHTime_)
          call write_mh_time
       case(Distr1D_)
          call write_distr_1d
       end select
    end do
  contains
    !=============================================================    
    subroutine write_mh_1d
      use SP_ModGrid,    ONLY: FootPoint_VB, NoShock_
      !\
      ! write output with 1D MH data in the format to be read by 
      ! IDL/TECPLOT; separate file is created for each field line, 
      ! name format is:
      ! MH_data_<iLon>_<iLat>_n<ddhhmmss>_n<iIter>.{out/dat}
      !/
      ! name of the output file
      character(len=100):: NameFile
      ! header for the file
      character(len=200):: StringHeader
      ! loop variables
      integer:: iBlock, iParticle, iVarPlot
      ! indexes of corresponding node, latitude and longitude
      integer:: iNode, iLat, iLon
      ! index of last particle on the field line
      integer:: iLast
      ! for better readability
      integer:: nVarPlot, nFluxPlot
      ! shock location
      integer:: iShock
      real   :: RShock
      integer, parameter:: RShock_     = Z_ + 2
      integer, parameter:: StartTime_  = RShock_ + 1
      integer, parameter:: StartJulian_= StartTime_ + 1
      real :: Param_I(LagrID_:StartJulian_)
      ! timetag
      character(len=15):: StringTime
      character(len=*), parameter:: NameSub='write_mh_1d'
      !------------------------------------------------------------------------
      !\
      !  Update number of time tags and write to tag list file
      !/
      if(iProc==0)then
         ! increase the file counter
         nTag = nTag + 1
         ! add to the tag list file
         NameFile = trim(NamePlotDir)//trim(NameTagFile)
         call open_file(file=NameFile, position='append', status='unknown', &
              NameCaller=NameSub)
         if(UseDateTime)then
            ! create date_time-iteration tag
            call get_date_time_string(SPTime, StringTime)
            write(UnitTmp_,'(a,i6.6)') 'e'//StringTime//'_n',iIter
         else
            ! create time-iteration tag
            call get_time_string(SPTime, StringTime(1:8))
            write(UnitTmp_,'(a,i6.6)') 't'//StringTime(1:8)//'_n',iIter
         end if
         call close_file
      end if
      !\
      ! Write ouput files themselves
      !/
      nVarPlot  = File_I(iFile)%nVarPlot
      nFluxPlot = File_I(iFile)%nFluxPlot
      StringHeader = &
           'MFLAMPA: data along a field line; '//&
           'Coordindate system: '//trim(TypeCoordSystem)           
      do iBlock = 1, nBlock
         iNode = iNode_B(iBlock)
         call get_node_indexes(iNode, iLon, iLat)
         ! set the file name
         if(UseDateTime)then
            call get_date_time_string(SPTime, StringTime)
            write(NameFile,'(a,i3.3,a,i3.3,a,i6.6,a)') &
                 trim(NamePlotDir)//NameMHData//'_',iLon,'_',iLat,&
                 '_e'//StringTime//'_n',iIter,&
                 File_I(iFile)%NameFileExtension
         else
            call get_time_string(SPTime, StringTime(1:8))
            write(NameFile,'(a,i3.3,a,i3.3,a,i6.6,a)') &
                 trim(NamePlotDir)//NameMHData//'_',iLon,'_',iLat,&
                 '_t'//StringTime(1:8)//'_n',iIter,&
                 File_I(iFile)%NameFileExtension
         end if
         ! get min and max particle indexes on this field line
         iLast  = nParticle_B(   iBlock)
         ! fill the output buffer
         File_I(iFile) % Buffer_II(1:nVarPlot, 1:iLast) = &
              State_VIB(File_I(iFile) % iVarPlot_V(1:nVarPlot), &
              1:iLast, iBlock)
         if(nFluxPlot > 0)File_I(iFile)%Buffer_II(&
              nVarPlot + 1:nVarPlot + nFluxPlot, 1:iLast) = &
              Flux_VIB(File_I(iFile)%iFluxPlot_V(1:nFluxPlot), &
              1:iLast, iBlock)
              
         !Parameters
         Param_I(LagrID_:Z_) = FootPoint_VB(LagrID_:Z_,iBlock)
         ! shock location
         if(iShock_IB(Shock_,iBlock)/=NoShock_)then
            iShock             = iShock_IB(Shock_,iBlock)
            Param_I(RShock_)   = State_VIB(R_,iShock,iBlock)
            Param_I(RShock_-1) = real(iShock)
         else
            Param_I(RShock_-1:RShock_) = -1.0
         end if
         ! start time in seconds from base year
         Param_I(StartTime_)  = StartTime
         ! start time in Julian date
         Param_I(StartJulian_)= StartTimeJulian
         ! print data to file
         call save_plot_file(&
              NameFile      = NameFile, &
              StringHeaderIn= StringHeader, &
              TypeFileIn    = File_I(iFile) % TypeFile, &
              nDimIn        = 1, &
              TimeIn        = SPTime, &
              nStepIn       = iIter, &
              CoordMinIn_D  = (/State_VIB(LagrID_,1,iBlock)/), &
              CoordMaxIn_D  = (/State_VIB(LagrID_,iLast,iBlock)/), &
              NameVarIn     = File_I(iFile) % NameVarPlot, &
              VarIn_VI      = &
              File_I(iFile) % Buffer_II(1:nVarPlot + nFluxPlot,1:iLast),&
              ParamIn_I    = Param_I(LagrID_:StartJulian_))
      end do
    end subroutine write_mh_1d
    !=============================================================
    subroutine write_mh_2d
      use SP_ModProc, ONLY: iComm, nProc
      use ModMpi
      ! write output with 2D MH data in the format to be read by IDL/TECPLOT;
      ! single file is created for all field lines, name format is
      ! MH_data_R=<Radius [AU]>_t<ddhhmmss>_n<iIter>.{out/dat}
      !-----------------------------------------------------------
      ! name of the output file
      character(len=100):: NameFile
      ! header of the output file
      character(len=200):: StringHeader
      ! loop variables
      integer:: iBlock, iParticle, iVarPlot, iVarIndex
      ! indexes of corresponding node, latitude and longitude
      integer:: iNode, iLat, iLon
      ! index of first/last particle on the field line
      integer:: iLast
      ! index of particle just above the radius
      integer:: iAbove
      ! radii of particles, added for readability
      real:: Radius0, Radius1
      ! interpolation weight
      real:: Weight
      ! for better readability
      integer:: nVarPlot, nFluxPlot
      ! MPI error
      integer:: iError
      ! skip a field line not reaching radius of output sphere
      logical:: DoPrint_I(nNode)
      ! additional parameters
      integer, parameter:: StartTime_  = 1
      integer, parameter:: StartJulian_= StartTime_ + 1
      real :: Param_I(1:StartJulian_)
      ! timetag
      character(len=15):: StringTime
      character(len=*), parameter:: NameSub='write_mh_2d'
      !-----------------------------------------------------------
      nVarPlot = File_I(iFile)%nVarPlot
      nFluxPlot= File_I(iFile)%nFluxPlot
      ! header for the output file
      StringHeader = &
           'MFLAMPA: data on a sphere at fixed heliocentric distance; '//&
           'Coordindate system: '//trim(TypeCoordSystem)           

      ! reset the output buffer
      File_I(iFile) % Buffer_II = 0

      ! reset, all field lines are printed reaching output sphere
      DoPrint_I = .true.

      ! go over all lines on the processor and find the point of 
      ! intersection with output sphere if present
      do iBlock = 1, nBlock
         iNode = iNode_B(iBlock)

         ! set the file name
         if(UseDateTime)then
            call get_date_time_string(SPTime, StringTime)
            write(NameFile,'(a,i4.4,f0.2,a,i6.6,a)')  &
                 trim(NamePlotDir)//NameMHData//'_R=',     &
                 int(File_I(iFile)%Radius),           &
                 File_I(iFile)%Radius - &
                 int(File_I(iFile)%Radius),           &
                 '_e'//StringTime//'_n', iIter, &
                 File_I(iFile) % NameFileExtension
         else
            call get_time_string(SPTime, StringTime(1:8))
            write(NameFile,'(a,i4.4,f0.2,a,i6.6,a)')  &
                 trim(NamePlotDir)//NameMHData//'_R=',     &
                 int(File_I(iFile) % Radius),         &
                 File_I(iFile) % Radius -             &
                 int(File_I(iFile)%Radius),           &
                 '_t'//StringTime(1:8)//'_n',         &
                 iIter, File_I(iFile)%NameFileExtension
         end if
         ! get max particle indexes on this field line
         iLast  = nParticle_B(iBlock)

         ! find the particle just above the given radius
         do iParticle = 1 , iLast
            Radius0 = State_VIB(R_, iParticle, iBlock)
            if( Radius0 > File_I(iFile)%Radius) then
               iAbove = iParticle
               !check if line started above output sphere, i.e. no intersection
               DoPrint_I(iNode) = iAbove /= 1
               EXIT
            end if
            ! check if reached the end, i.e. there is no intersection
            if(iParticle == iLast) &
                 DoPrint_I(iNode) = .false.
         end do

         ! if no intersection found -> proceed to the next line
         if(.not.DoPrint_I(iNode)) CYCLE

         ! intersection is found -> get data at that location;
         ! interpolate data and fill buffer
         Radius0 = sum(State_VIB(X_:Z_, iAbove-1, iBlock)**2)**0.5
         Radius1 = sum(State_VIB(X_:Z_, iAbove,   iBlock)**2)**0.5
         Weight  = (File_I(iFile)%Radius - Radius0) / (Radius1 - Radius0)
         ! interpolate each requested variable
         do iVarPlot = 1, nVarPlot
            iVarIndex = File_I(iFile) % iVarPlot_V(iVarPlot)
            File_I(iFile) % Buffer_II(iVarPlot, iNode) = &
                 State_VIB(iVarIndex, iAbove-1, iBlock) * (1-Weight) + &
                 State_VIB(iVarIndex, iAbove,   iBlock) *    Weight
         end do
         do iVarPlot = 1, nFluxPlot
            iVarIndex = File_I(iFile)%iFluxPlot_V(iVarPlot)
            File_I(iFile) % Buffer_II(iVarPlot + nVarPlot, iNode) = &
                 Flux_VIB(iVarIndex, iAbove-1, iBlock) * (1-Weight) + &
                 Flux_VIB(iVarIndex, iAbove,   iBlock) *    Weight
         end do

      end do

      ! gather interpolated data on the source processor
      if(nProc > 1)then
         if(iProc==0)then
            call MPI_Reduce(MPI_IN_PLACE, File_I(iFile) % Buffer_II, &
                 nNode * (nVarPlot + nFluxPlot), MPI_REAL, MPI_Sum, &
                 0, iComm, iError)
            call MPI_Reduce(MPI_IN_PLACE, DoPrint_I, &
                 nNode, MPI_Logical, MPI_Land, &
                 0, iComm, iError)
         else
            call MPI_Reduce(File_I(iFile)%Buffer_II, &
                 File_I(iFile)%Buffer_II,&
                 nNode * (nVarPlot + nFluxPlot), MPI_REAL, MPI_Sum,&
                 0, iComm, iError)
            call MPI_Reduce(DoPrint_I, DoPrint_I, &
                 nNode, MPI_Logical, MPI_Land, &
                 0, iComm, iError)
         end if
      end if

      ! start time in seconds from base year
      Param_I(StartTime_)  = StartTime
      ! start time in Julian date
      Param_I(StartJulian_)= StartTimeJulian

      if(iProc==0)&
           ! print data to file
           call save_plot_file(&
           NameFile      = NameFile, &
           StringHeaderIn= StringHeader, &
           TypeFileIn    = File_I(iFile) % TypeFile, &
           nDimIn        = 1, &
           TimeIn        = SPTime, &
           nStepIn       = iIter, &
           ParamIn_I     = Param_I, &
           Coord1In_I    = real(pack(iNodeIndex_I, MASK=DoPrint_I)),&
           NameVarIn     = File_I(iFile) % NameVarPlot, &
           VarIn_VI      = &
           reshape(&
           pack(File_I(iFile) % Buffer_II(1:nVarPlot+nFluxPlot,1:nNode),&
           MASK = spread(DoPrint_I, 1, nVarPlot+nFluxPlot)), &
           (/nVarPlot+nFluxPlot, count(DoPrint_I)/)))
    end subroutine write_mh_2d
    !=============================================================
    subroutine write_mh_time
      ! write output w/time series MH data in format to be read by IDL/TECPLOT;
      ! a file is created for each field lines, name format is
      ! MH_data_R=<Radius [AU]>_<iLon>_<iLat>.{out/dat}
      ! the file has no timetag as it is updated during the run
      !-----------------------------------------------------------
      ! name of the output file
      character(len=100):: NameFile
      ! loop variables
      integer:: iBlock, iParticle, iVarPlot, iVarIndex
      ! indexes of corresponding node, latitude and longitude
      integer:: iNode, iLat, iLon
      ! index of first/last particle on the field line
      integer:: iLast
      ! index of particle just above the radius
      integer:: iAbove
      ! radii of particles, added for readability
      real:: Radius0, Radius1
      ! interpolation weight
      real:: Weight
      ! for better readability
      integer:: nVarPlot, nFluxPlot
      ! skip a field line if it fails to reach radius of output sphere
      logical:: DoPrint
      ! size of the already written data
      integer:: nDataLine
      ! current size of the buffer
      integer:: nBufferSize
      ! whether the output file already exists
      logical:: IsPresent
      ! additional parameters
      integer, parameter:: StartTime_  = 1
      integer, parameter:: StartJulian_= StartTime_ + 1
      real :: Param_I(1:StartJulian_)
      ! header for the file
      character(len=200):: StringHeader
      character(len=*), parameter:: NameSub='write_mh_time'
      !------------------------------------------------------------------------
      nVarPlot = File_I(iFile)%nVarPlot
      nFluxPlot= File_I(iFile)%nFluxPlot
      ! set header
      StringHeader = &
           'MFLAMPA: data on a field line at fixed heliocentric distance; '//&
           'Coordindate system: '//trim(TypeCoordSystem)           

      !\
      ! at first call, remove files if they exist to reset time series output
      !/
      if(File_I(iFile) % IsFirstCall)then
         ! mark that the 1st call has alredy happened
         File_I(iFile) % IsFirstCall = .false.
         ! go over list of lines and remove file for each one
         do iBlock = 1, nBlock
            iNode = iNode_B(iBlock)
            call get_node_indexes(iNode, iLon, iLat)
            ! set the file name
            write(NameFile,'(a,i4.4,f0.2,a,i3.3,a,i3.3,a)') &
                 trim(NamePlotDir)//NameMHData//'_R=', &
                 int(File_I(iFile)%Radius),&
                 File_I(iFile) % Radius - int(File_I(iFile) % Radius), &
                 '_', iLon, '_', iLat, File_I(iFile) % NameFileExtension
            call remove_file(NameFile)
         end do
      end if

      ! reset the output buffer
      File_I(iFile) % Buffer_II = 0

      ! go over all lines on the processor and find the point of intersection
      ! with output sphere if present
      do iBlock = 1, nBlock

         ! reset, the field line is printed unless fail to reach output sphere
         DoPrint = .true.

         iNode = iNode_B(iBlock)
         call get_node_indexes(iNode, iLon, iLat)

         ! get max particle indexes on this field line
         iLast  = nParticle_B(   iBlock)

         ! find the particle just above the given radius
         do iParticle = 1 , iLast
            Radius0 = sum(State_VIB(X_:Z_, iParticle, iBlock)**2)**0.5
            if( Radius0 > File_I(iFile) % Radius)then
               iAbove = iParticle
               !check if line started above output sphere, i.e. no intersection
               DoPrint = iAbove /= 1
               EXIT
            end if
            ! check if reached the end, i.e. there is no intersection
            if(iParticle == iLast) &
                 DoPrint = .false.
         end do
         ! if no intersection found -> proceed to the next line
         if(.not.DoPrint) CYCLE

         ! set the file name
         write(NameFile,'(a,i4.4,f0.2,a,i3.3,a,i3.3,a)') &
              trim(NamePlotDir)//NameMHData//'_R=', &
              int(File_I(iFile)%Radius),&
              File_I(iFile) % Radius - int(File_I(iFile) % Radius), &
              '_', iLon, '_', iLat, File_I(iFile) % NameFileExtension
         !\
         ! if file already exists -> read its content
         nDataLine = 0
         inquire(FILE=NameFile, EXIST=IsPresent)
         if(IsPresent)then

            ! first, determine its size
            call read_plot_file(&
                 NameFile   = NameFile, &
                 TypeFileIn = File_I(iFile) % TypeFile, &
                 n1Out      = nDataLine)
            ! if buffer is too small then reallocate it
            nBufferSize = ubound(File_I(iFile)%Buffer_II, 2)
            if(nBufferSize < nDataLine + 1)then
               deallocate(File_I(iFile)%Buffer_II)
               allocate(File_I(iFile)%Buffer_II(&
                    nVarPlot + nFluxPlot + 1, 2*nBufferSize))
            end if

            ! read the data itself
            call read_plot_file(&
                 NameFile   = NameFile, &
                 TypeFileIn = File_I(iFile)%TypeFile,&
                 Coord1Out_I= File_I(iFile)%Buffer_II(1+nVarPlot+nFluxPlot,:),&
                 VarOut_VI  = File_I(iFile)%Buffer_II(1:nVarPlot+nFluxPlot,:))
         end if

         !\
         ! add new data
         nDataLine = nDataLine + 1
         ! put time into buffer
         File_I(iFile)%Buffer_II(nVarPlot+nFluxPlot+1,nDataLine) = SPTime
         ! interpolate data and fill buffer
         Radius0 = sum(State_VIB(X_:Z_, iAbove-1, iBlock)**2)**0.5
         Radius1 = sum(State_VIB(X_:Z_, iAbove,   iBlock)**2)**0.5
         Weight  = (File_I(iFile)%Radius - Radius0) / (Radius1 - Radius0)
         ! interpolate each requested variable
         do iVarPlot = 1, nVarPlot
            iVarIndex = File_I(iFile)%iVarPlot_V(iVarPlot)
            File_I(iFile)%Buffer_II(iVarPlot, nDataLine) = &
                 State_VIB(iVarIndex, iAbove-1, iBlock) * (1-Weight) + &
                 State_VIB(iVarIndex, iAbove,   iBlock) *    Weight
         end do
         do iVarPlot = 1, nFluxPlot
            iVarIndex = File_I(iFile)%iFluxPlot_V(iVarPlot)
            File_I(iFile)%Buffer_II(iVarPlot + nVarPlot, nDataLine) = &
                 Flux_VIB(iVarIndex, iAbove-1, iBlock) * (1-Weight) + &
                 Flux_VIB(iVarIndex, iAbove,   iBlock) *    Weight
         end do

         ! start time in seconds from base year
         Param_I(StartTime_)  = StartTime
         ! start time in Julian date
         Param_I(StartJulian_)= StartTimeJulian
         ! reprint data to file
         call save_plot_file(&
              NameFile      = NameFile, &
              StringHeaderIn= StringHeader, &
              TypeFileIn    = File_I(iFile) % TypeFile, &
              nDimIn        = 1, &
              TimeIn        = SPTime, &
              nStepIn       = iIter, &
              ParamIn_I     = Param_I, &
              Coord1In_I    = &
              File_I(iFile) % Buffer_II(1+nVarPlot+nFluxPlot,1:nDataLine), &
              NameVarIn     = File_I(iFile) % NameVarPlot, &
              VarIn_VI      = &
              File_I(iFile) % Buffer_II(1:nVarPlot+nFluxPlot,1:nDataLine))
      end do

    end subroutine write_mh_time
    !=============================================================
    subroutine write_distr_1d
      use SP_ModGrid, ONLY: S_
      !\
      ! write file with distribution in the format to be read by IDL/TECPLOT;
      ! separate file is created for each field line, name format is
      ! Distribution_<iLon>_<iLat>_t<ddhhmmss>_n<iIter>.{out/dat}
      !/
      ! name of the output file
      character(len=100):: NameFile
      ! loop variables
      integer:: iBlock, iParticle, iVarPlot
      ! indexes of corresponding node, latitude and longitude
      integer:: iNode, iLat, iLon
      ! index of first/last particle on the field line
      integer:: iLast
      ! scale and conversion factor
      real:: Scale_I(0:nP+1), Factor_I(0:nP+1)
      real:: Unity_I(0:nP+1) = 1.0
      ! timetag
      character(len=15):: StringTime
      character(len=*), parameter:: NameSub='write_distr_1d'
      !-----------------------------------------------------------
      select case(File_I(iFile)%iScale)
      case(Momentum_)
         Scale_I = Log10MomentumSI_I
         Factor_I= Unity_I
      case(Energy_)
         Scale_I = Log10KinEnergySI_I
         Factor_I= DMomentumOverDEnergy_I
      end select
      do iBlock = 1, nBlock
         iNode = iNode_B(iBlock)
         call get_node_indexes(iNode, iLon, iLat)

         ! set the file name
         if(UseDateTime)then
            call get_date_time_string(SPTime, StringTime)
            write(NameFile,'(a,i3.3,a,i3.3,a,i6.6,a)') &
                 trim(NamePlotDir)//'Distribution_',iLon,'_',iLat,&
                 '_e'//StringTime//'_n',iIter,&
                 File_I(iFile) % NameFileExtension
         else
            call get_time_string(SPTime, StringTime(1:8))
            write(NameFile,'(a,i3.3,a,i3.3,a,i6.6,a)') &
                 trim(NamePlotDir)//'Distribution_',iLon,'_',iLat,&
                 '_t'//StringTime(1:8)//'_n',iIter,&
                 File_I(iFile) % NameFileExtension
         end if
         ! get max particle indexes on this field line
         iLast  = nParticle_B(   iBlock)

         do iParticle = 1, nParticleMax
            ! reset values outside the line's range
            if(iParticle > iLast)then
               File_I(iFile) % Buffer_II(:,iParticle) = 0.0
               CYCLE
            end if
            ! the actual distribution
            File_I(iFile) % Buffer_II(:,iParticle) = &
                 log10(Distribution_IIB(1:nP,iParticle,iBlock)*Factor_I(1:nP))
            ! account for the requested output
            select case(File_I(iFile) % iTypeDistr)
            case(CDF_)
               ! do nothing
            case(DEF_)
               File_I(iFile) % Buffer_II(:,iParticle) = &
                    File_I(iFile) % Buffer_II(:,iParticle) + &
                    2*Log10MomentumSI_I
            end select
         end do
         ! print data to file
         call save_plot_file(&
              NameFile   = NameFile, &
              TypeFileIn = File_I(iFile) % TypeFile, &
              nDimIn     = 2, &
              TimeIn     = SPTime, &
              nStepIn    = iIter, &
              Coord1In_I = Scale_I, &
              Coord2In_I = State_VIB(S_,1:iLast,iBlock), &
              NameVarIn  = File_I(iFile) % NameVarPlot, &
              VarIn_II   = File_I(iFile) % Buffer_II(:,1:iLast))
      end do
    end subroutine write_distr_1d
  end subroutine save_plot_all
  !===============================================================
  subroutine finalize
     use SP_ModGrid, ONLY: nLat, nLon
    ! write the header file that contains necessary information 
    ! for reading input files in a separate run
    character(len=*), parameter:: NameSub='write_mh_1d_header'
    !-------------------------------------------------------------
    if(.not.DoWriteHeader)RETURN
    ! performed on root proc only
    if (iProc/=0) RETURN
    ! write the header file
    call open_file(file= trim(NamePlotDir)//trim(NameHeaderFile),&
         NameCaller=NameSub)
    write(UnitTmp_,*)
    write(UnitTmp_,'(a)')'#CHECKGRIDSIZE'
    write(UnitTmp_,'(i8,a32)') nParticleMax,'nParticleMax'
    write(UnitTmp_,'(i8,a32)') nLon,     'nLon'
    write(UnitTmp_,'(i8,a32)') nLat,     'nLat'
    write(UnitTmp_,*)
    write(UnitTmp_,'(a)')'#MHDATA'
    write(UnitTmp_,'(a,a35)')trim(TypeMHDataFile),'  TypeFile'
    write(UnitTmp_,'(i8,a32)')nTag,'nFileRead'
    write(UnitTmp_,'(a,a29)')trim(NameTagFile), '  NameTagFile'
    write(UnitTmp_,*)
    write(UnitTmp_,'(a)')'#END'
    write(UnitTmp_,*)
    call close_file
  end subroutine finalize
  !===============================================================
  subroutine get_time_string(Time, StringTime)
    ! the subroutine converts real variable Time into a string,
    ! the structure of the string is 'ddhhmmss', 
    ! i.e shows number of days, hours, minutes and seconds 
    ! after the beginning of the simulation
    real,             intent(in) :: Time
    character(len=8), intent(out):: StringTime
    !-------------------------------------------------------------
    ! This is the value if the time is too large
    StringTime = '99999999'
    if(Time < 100.0*86400) &
         write(StringTime,'(i2.2,i2.2,i2.2,i2.2)') &
         int(                  Time          /86400.), & ! # days
         int((Time-(86400.*int(Time/86400.)))/ 3600.), & ! # hours
         int((Time-( 3600.*int(Time/ 3600.)))/   60.), & ! # minutes
         int( Time-(   60.*int(Time/   60.)))            ! # seconds
  end subroutine get_time_string
  !===============================================================
  subroutine get_date_time_string(Time, StringTime)
    use SP_ModTime,    ONLY: StartTime
    use ModTimeConvert,ONLY: time_real_to_int
    ! the subroutine converts real variable Time into a string,
    ! the structure of the string is 'ddhhmmss', 
    ! i.e shows number of days, hours, minutes and seconds 
    ! after the beginning of the simulation
    real,              intent(in) :: Time
    character(len=15), intent(out):: StringTime
    integer :: iTime_I(7)
    !-------------------------------------------------------------
    call time_real_to_int(StartTime+Time, iTime_I)
    write(StringTime,'(i4.4,i2.2,i2.2,a,i2.2,i2.2,i2.2)')&
         iTime_I(1:3),'_',iTime_I(4:6)
  end subroutine get_date_time_string
end module SP_ModPlot
