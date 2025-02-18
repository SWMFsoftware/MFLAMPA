!  Copyright (C) 2002 Regents of the University of Michigan,
!  portions used with permission
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module SP_ModTriangulate

  use ModUtilities,       ONLY: CON_stop

  implicit none

  PRIVATE ! Except

  SAVE

  public:: read_param    ! Read module parameters

  ! Test triangulation
  logical, public :: DoTestTri = .false.
  integer, public :: nLocTestTri = 1
  real, allocatable, public :: XyzLocTestTri_II(:,:)

  ! If we use poles in triangulation
  logical, public :: UsePoleTri   = .false.
  logical, public :: UsePlanarTri = .true.
  real,    public, parameter :: iSouthPoleTri_ = -1, iNorthPoleTri_ = -2
contains
  !============================================================================
  subroutine read_param(NameCommand)

    use ModReadParam, ONLY: read_var
    character(len=*), intent(in) :: NameCommand ! From PARAM.in
    character(len=*), parameter:: NameSub = 'read_param'
    !--------------------------------------------------------------------------
    select case(NameCommand)
    case("#TESTTRIANGULATE")
       ! To do
       call read_var('TestTriangulate', DoTestTri)
       if(DoTestTri) then
          call read_var('nLocTestTri', nLocTestTri)
          allocate(XyzLocTestTri_II(3, nLocTestTri))
       end if
    case('#TRIANGULATION')
       ! get pole triangulartion flag
       call read_var('UsePoleTriangulation', UsePoleTri)
       ! get the triangulation approach flag
       call read_var('UsePlanarTriangles', UsePlanarTri)
    case default
       call CON_stop(NameSub//': Unknown command '//NameCommand)
    end select

  end subroutine read_param
  !============================================================================
end module SP_ModTriangulate
!==============================================================================
