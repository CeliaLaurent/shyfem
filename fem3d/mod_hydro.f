
!--------------------------------------------------------------------------
!
!    Copyright (C) 2015-2017,2019  Georg Umgiesser
!
!    This file is part of SHYFEM.
!
!    SHYFEM is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    SHYFEM is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with SHYFEM. Please see the file COPYING in the main directory.
!    If not, see <http://www.gnu.org/licenses/>.
!
!    Contributions to this file can be found below in the revision log.
!
!--------------------------------------------------------------------------

! revision log :
!
! 05.06.2015	ggu	changed VERS_7_1_12
! 10.07.2015	ggu	changed VERS_7_1_50
! 28.04.2016	ggu	changed VERS_7_5_9
! 16.05.2017	ggu	changed VERS_7_5_27
! 16.02.2019	ggu	changed VERS_7_5_60
! 21.05.2019	ggu	changed VERS_7_5_62

!**************************************************************************

!==================================================================
	module mod_hydro
!==================================================================

	implicit none

	integer, private, save :: nkn_hydro = 0
	integer, private, save :: nel_hydro = 0
	integer, private, save :: nlv_hydro = 0

	real, allocatable, save :: zov(:), znv(:)
	real, allocatable, save :: zeov(:,:), zenv(:,:)
	real, allocatable, save :: utlov(:,:)
	real, allocatable, save :: utlnv(:,:)
	real, allocatable, save :: vtlov(:,:)
	real, allocatable, save :: vtlnv(:,:)
c-----------------------------------------
        ! creating next arrays used in new3di.f instead of allocating
        ! them at every iteration and for every element solves an
        ! optimization issue encountered on the m100 cluster.
	real, allocatable, save ::   hact(:) 
	real, allocatable, save ::   rhact(:)
	real, allocatable, save ::   alev(:) 
	double precision, allocatable, save ::  rmat(:)
	double precision, allocatable, save ::  smat(:,:)
	double precision, allocatable, save ::  s2dmat(:,:)	!for 2D
	double precision, allocatable, save ::  rvec(:)  	!ASYM (3 systems to solve)
	double precision, allocatable, save ::  rvecp(:) 	!ASYM (3 systems to solve)
	double precision, allocatable, save ::  solv(:)  	!ASYM (3 systems to solve)

!==================================================================
	contains
!==================================================================

        subroutine mod_hydro_init(nkn,nel,nlv)
        
	use levels, only : nlvdi
        implicit none

        integer nkn, nel, nlv
        
        if( nkn == nkn_hydro .and. nel == nel_hydro .and.
     +      nlv == nlv_hydro ) return

        if( nel > 0 .or. nkn > 0 .or. nlv > 0 ) then
          if( nel == 0 .or. nkn == 0 .or. nlv == 0 ) then
            write(6,*) 'nel,nkn,nlv: ',nel,nkn,nlv
            stop 'error stop mod_hydro_init: incompatible parameters'
          end if
        end if

        if( nkn_hydro > 0 ) then
          deallocate(zov)
          deallocate(znv)
        
          deallocate(zeov)
          deallocate(zenv)

          deallocate(utlov)
          deallocate(utlnv)
          deallocate(vtlov)
          deallocate(vtlnv)
        end if

        nkn_hydro = nkn
        nel_hydro = nel
        nlv_hydro = nlv

        if( nkn == 0 ) return

        allocate(zov(nkn))
        allocate(znv(nkn))

        allocate(zeov(3,nel))
        allocate(zenv(3,nel))

        allocate(utlov(nlv,nel))
        allocate(utlnv(nlv,nel))
        allocate(vtlov(nlv,nel))
        allocate(vtlnv(nlv,nel))

	allocate(hact(0:nlvdi+1))
	allocate(rhact(0:nlvdi+1))
	allocate(alev(0:nlvdi))
	allocate(rmat(10*nlvdi))
	allocate(smat(-2:2,2*nlvdi))
	allocate(s2dmat(-1:1,2))		!for 2D
	allocate(rvec(6*nlvdi))		!ASYM (3 systems to solve)
	allocate(rvecp(6*nlvdi))		!ASYM (3 systems to solve)
	allocate(solv(6*nlvdi))		!ASYM (3 systems to solve)

	zov = 0.
	znv = 0.
	zeov = 0.
	zenv = 0.
	utlov = 0.
	vtlov = 0.
	utlnv = 0.
	vtlnv = 0.

        end subroutine mod_hydro_init

!==================================================================
        end module mod_hydro
!==================================================================

