c
c $Id: ousinf.f,v 1.15 2009-11-18 16:50:37 georg Exp $
c
c info on OUS files
c
c revision log :
c
c 02.09.2003	ggu	adapted to new OUS format
c 24.01.2005	ggu	computes maximum velocities for 3D (only first level)
c 16.10.2007	ggu	new debug routine
c 27.10.2009    ggu     include evmain.h, compute volume
c 23.03.2011    ggu     compute real u/v-min/max of first level
c 16.12.2011    ggu     bug fix: call to init_sigma_info, makehev (common hev)
c 25.01.2013    ggu     utility routines to ousutil.f
c
c***************************************************************

	program ousinf

c reads ous file and writes info to terminal
c
c we would not even need to read basin

	use mod_depth !COMMON_GGU_SUBST
	use mod_hydro !COMMON_GGU_SUBST
	use evgeom !COMMON_GGU_SUBST
	use levels !COMMON_GGU_SUBST
	use basin !COMMON_GGU_SUBST

	implicit none

        include 'param.h'
COMMON_GGU_DELETED	include 'evmain.h'

COMMON_GGU_DELETED	include 'basin.h'
	include 'simul.h'


COMMON_GGU_DELETED	include 'nlevel.h'
COMMON_GGU_DELETED	include 'levels.h'

COMMON_GGU_DELETED	include 'depth.h'

COMMON_GGU_DELETED	include 'hydro.h'

        real ut2v(neldim)
        real vt2v(neldim)
        real u2v(neldim)
        real v2v(neldim)

        integer nvers,nin,nlev
        integer itanf,itend,idt,idtous
	integer it,ie,i
        integer ierr,nread,ndry
        integer nknous,nelous,nlvous
        real href,hzoff,hlvmin
	real volume
	real zmin,zmax
	real umin,umax
	real vmin,vmax

c	integer rdous,rfous
	integer iapini,ideffi

c-----------------------------------------------------------------

	nread=0

c-----------------------------------------------------------------
c open basin and simulation
c-----------------------------------------------------------------

	call shyfem_copyright('ousinf - Info on OUS files')

	if(iapini(3,nkndim,neldim,0).eq.0) then
		stop 'error stop : iapini'
	end if

	call open_ous_type('.ous','old',nin)

	call read_ous_header(nin,nkndim,neldim,nlvdim,ilhv,hlv,hev)
	call ous_get_params(nin,nknous,nelous,nlvous)
	nlv = nlvous

	call set_ev
        call init_sigma_info(nlv,hlv)

c-----------------------------------------------------------------
c loop on data of simulation
c-----------------------------------------------------------------

	do while(.true.)

	  call ous_read_record(nin,it,nlvdim,ilhv,znv,zenv
     +				,utlnv,vtlnv,ierr)

          if(ierr.gt.0) write(6,*) 'error in reading file : ',ierr
          if(ierr.ne.0) goto 100

	  nread=nread+1

	  call mima(znv,nknous,zmin,zmax)
          call comp_barotropic(nel,nlvdim,ilhv,utlnv,vtlnv,ut2v,vt2v)
          call comp_vel2d(nel,hev,zenv,ut2v,vt2v,u2v,v2v
     +                          ,umin,vmin,umax,vmax)
	  call compute_volume(nel,zenv,hev,volume)

c          call debug_write_node(it,nread,nkndim,neldim,nlvdim,nkn,nel,nlv
c     +          ,nen3v,zenv,znv,utlnv,vtlnv)

	  write(6,*) 
	  write(6,*) 'time : ',it
	  write(6,*) 
	  write(6,*) 'zmin/zmax : ',zmin,zmax
	  write(6,*) 'umin/umax : ',umin,umax
	  write(6,*) 'vmin/vmax : ',vmin,vmax
	  write(6,*) 'volume    : ',volume

	end do

c-----------------------------------------------------------------
c end of loop
c-----------------------------------------------------------------

  100	continue

	write(6,*)
	write(6,*) nread,' records read'
	write(6,*)

c-----------------------------------------------------------------
c end of routine
c-----------------------------------------------------------------

	end

c******************************************************************

	subroutine write_zeta(nkn,znv)

	implicit none

	integer nkn
	real znv(nkn)

	integer k

	open(1,file='zinit.dat',form='unformatted',status='unknown')
	write(1) (znv(k),k=1,nkn)
	close(1)

	end

c******************************************************************
