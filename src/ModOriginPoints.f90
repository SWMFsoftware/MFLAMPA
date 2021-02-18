!  Copyright (C) 2002 Regents of the University of Michigan
! portions used with permission
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module SP_ModOriginPoints
  implicit none
  PRIVATE  ! Except
  public :: read_param
  public :: get_origin_points
  ! Grid size, boundaries, coordinates
  ! Starting position of field lines in Rs
  real         :: ROrigin = 2.5
  ! Size of angular grid, latitude and longitude, at origin
  ! surface R=ROrigin
  real         :: LonMin = 0.0, LonMax = 360.0
  real         :: LatMin = -70.0, LatMax = 70.0
contains
  !============================================================================

  subroutine read_param
    ! Read input parameters for SP component
    use ModReadParam, ONLY: read_var

    character(len=*), parameter:: NameSub = 'read_param'
    !--------------------------------------------------------------------------
    call read_var('ROrigin', ROrigin)
    call read_var('LonMin', LonMin)
    call read_var('LatMin', LatMin)
    call read_var('LonMax', LonMax)
    call read_var('LatMax', LatMax)
  end subroutine read_param
  !============================================================================
  subroutine get_origin_points
    use ModNumConst,       ONLY: cDegToRad
    use ModCoordTransform, ONLY: rlonlat_to_xyz
    use SP_ModGrid,        ONLY: nLat, nLon, X_, Z_, MHData_VIB, nParticle_B,&
         nBlock, iNode_B, get_node_indexes
   
    integer:: iLat, iLon, iNode, iBlock
    ! Sell size on the origin surface, per line
    real         ::  DLon, DLat

    ! convert angels from degrees to radians

    !--------------------------------------------------------------------------
    LonMax = LonMax*cDegToRad
    LonMin = LonMin*cDegToRad

    ! angular grid's step
    DLon = (LonMax - LonMin)/nLon

    ! convert angels from degrees to radians
    LatMax = LatMax*cDegToRad
    LatMin = LatMin*cDegToRad

    ! angular grid's step
    DLat = (LatMax - LatMin)/nLat
    do iBlock = 1, nBlock
       call get_node_indexes(iNode_B(iBlock), iLon, iLat)
       nParticle_B(iBlock) = 1
       call rlonlat_to_xyz([ROrigin, LonMin + (iLon - 0.5)*DLon, &
            LatMin + (iLat - 0.5)*DLat], MHData_VIB(X_:Z_,1,iBlock))
    end do
  end subroutine get_origin_points
  !============================================================================
end module SP_ModOriginPoints
