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
c 16.12.2011    ggu     bug fix: call to init_sigma_info and makehev (common hev)
c
c***************************************************************

	program ousinf

c reads ous file and writes info to terminal
c
c we would not even need to read basin

	implicit none

        include 'param.h'
	include 'evmain.h'

	character*80 descrr,descrp
	common /descrr/ descrr
	common /descrp/ descrp
	integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw

	real xgv(nkndim), ygv(nkndim)
	real hm3v(3,neldim)
	integer nen3v(3,neldim)
	integer ipev(neldim), ipv(nkndim)
	integer iarv(neldim)
	common /xgv/xgv, /ygv/ygv
	common /hm3v/hm3v
	common /nen3v/nen3v
	common /ipev/ipev, /ipv/ipv
	common /iarv/iarv

	integer ilhv(neldim)
	real hlv(nlvdim)
        real utlnv(nlvdim,neldim)
        real vtlnv(nlvdim,neldim)
	common /ilhv/ilhv
	common /hlv/hlv
        common /utlnv/utlnv
        common /vtlnv/vtlnv

	real hev(neldim)
	common /hev/hev

	real znv(nkndim)
	real zenv(3,neldim)

        integer nvers,nin,nlv
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
c initialize basin and simulation
c-----------------------------------------------------------------

	nread=0

	if(iapini(3,nkndim,neldim,0).eq.0) then
		stop 'error stop : iapini'
	end if

	call set_ev

	nin=ideffi('datdir','runnam','.ous','unform','old')
	if(nin.le.0) goto 100

c-----------------------------------------------------------------
c read header of simulation
c-----------------------------------------------------------------

	nvers=1
        call rfous(nin
     +			,nvers
     +			,nknous,nelous,nlvous
     +			,href,hzoff
     +			,descrp
     +			,ierr)

	nlv=nlvous
	call dimous(nin,nkndim,neldim,nlvdim)

        write(6,*)
        write(6,*)   descrp
        write(6,*)
        write(6,*) ' nvers        : ',nvers
        write(6,*) ' href,hzoff   : ',href,hzoff
        write(6,*) ' nkn,nel      : ',nknous,nelous
        write(6,*) ' nlv          : ',nlvous
        write(6,*)

	call rsous(nin,ilhv,hlv,hev,ierr)

        call init_sigma_info(nlv,hlv)
	call makehev(hev)

	call init_date()

c-----------------------------------------------------------------
c loop on data of simulation
c-----------------------------------------------------------------

  300   continue

        call rdous(nin,it,nlvdim,ilhv,znv,zenv,utlnv,vtlnv,ierr)

        if(ierr.gt.0) then
		write(6,*) 'error in reading file : ',ierr
		goto 100
        else if(ierr.lt.0) then
		goto 100
	end if

	nread=nread+1

	call mima(znv,nknous,zmin,zmax)
        call comp_vel(1,nel,hev,zenv,nlvdim,utlnv,vtlnv
     +			,umin,vmin,umax,vmax)
	call compute_volume(nel,zenv,hev,volume)

	if( it .eq. 410400 ) then
	  write(6,*)
	  write(6,*) '****** elaborating data for it = ',it
	  write(6,*)
	  call elab_z(nkn,znv)
	  write(6,*)
	  stop
	end if

c        call debug_write_node(it,nread,nkndim,neldim,nlvdim,nkn,nel,nlv
c     +          ,nen3v,zenv,znv,utlnv,vtlnv)

	write(6,*) 
	write(6,*) 'time : ',it
	call get_date(it)
	write(6,*) 
	write(6,*) 'zmin/zmax : ',zmin,zmax
	write(6,*) 'umin/umax : ',umin,umax
	write(6,*) 'vmin/vmax : ',vmin,vmax
	write(6,*) 'volume    : ',volume

	goto 300

  100	continue

c-----------------------------------------------------------------
c end of loop
c-----------------------------------------------------------------

	write(6,*)
	write(6,*) nread,' records read'
	write(6,*)

c-----------------------------------------------------------------
c end of routine
c-----------------------------------------------------------------

	stop
	end

c******************************************************************

	subroutine compute_volume(nel,zenv,hev,volume)

	implicit none

	include 'param.h'
	include 'evmain.h'

	integer nel
	real zenv(3,neldim)
	real hev(neldim)
	real volume

	integer ie,ii
	real zav,area
	double precision vol,voltot,areatot

	voltot = 0.
	areatot = 0.

	do ie=1,nel
	  zav = 0.
	  do ii=1,3
	    zav = zav + zenv(ii,ie)
	  end do
	  area = 12. * ev(10,ie)
	  vol = area * (hev(ie) + zav/3.)
	  voltot = voltot + vol
	  !areatot = areatot + area
	end do

	volume = voltot

	end

c******************************************************************

        subroutine comp_vel(level,nel,hev,zenv,nlvdim,utlnv,vtlnv
     +			,umin,vmin,umax,vmax)

        implicit none

        integer level
        integer nel
        real hev(1)
        real zenv(3,1)
        integer nlvdim
        real utlnv(nlvdim,1)
        real vtlnv(nlvdim,1)
        real umin,vmin
        real umax,vmax

        integer ie,ii
        real zmed,hmed,u,v

	umin =  1.e+30
	vmin =  1.e+30
        umax = -1.e+30
        vmax = -1.e+30

        do ie=1,nel
          zmed = 0.
          do ii=1,3
            zmed = zmed + zenv(ii,ie)
          end do
          zmed = zmed / 3.
          hmed = hev(ie) + zmed
          if( hmed .le. 0. ) stop 'error stop hmed...'

          u = utlnv(level,ie) / hmed
          v = vtlnv(level,ie) / hmed

          umin = min(umin,u)
          vmin = min(vmin,v)
          umax = max(umax,u)
          vmax = max(vmax,v)
        end do

        end

c******************************************************************

        subroutine debug_write_node(it,nrec
     +		,nkndim,neldim,nlvdim,nkn,nel,nlv
     +          ,nen3v,zenv,znv,utlnv,vtlnv)

c debug write

        implicit none

        integer it,nrec
        integer nkndim,neldim,nlvdim,nkn,nel,nlv
        integer nen3v(3,neldim)
        real znv(nkndim)
        real zenv(3,neldim)
        real utlnv(nlvdim,neldim)
        real vtlnv(nlvdim,neldim)

        integer ie,ii,k,l,ks
        logical bk

        ks = 6068

        write(66,*) 'time: ',it,nrec
        write(66,*) 'kkk: ',znv(ks)

        do ie=1,nel
          bk = .false.
          do ii=1,3
            k = nen3v(ii,ie)
            if( k .eq. ks ) then
              write(66,*) 'ii: ',ii,ie,zenv(ii,ie)
              bk = .true.
            end if
          end do
          if( bk ) then
          do l=1,nlv
            write(66,*) 'ie: ',ie,l,utlnv(l,ie),vtlnv(l,ie)
          end do
          end if
        end do

        end

c******************************************************************

	subroutine get_date(it)

	implicit none

	integer it
	integer year,month,day,hour,min,sec
	character*20 line

	call dts2dt(it,year,month,day,hour,min,sec)
	call dtsform(year,month,day,hour,min,sec,line)

	write(6,*) line

	end

c******************************************************************

	subroutine init_date()

	implicit none

	character*60 name
	integer date,time

	name='date0'

	open(1,file=name,status='old',form='formatted')
	read(1,*) date
	close(1)

	time = 0
	call dtsini(date,time)

	end

c******************************************************************

	subroutine elab_z(nkn,znv)

	implicit none

	integer nkn
	real znv(nkn)

	character*60 name
	real x,y,aux,ssl
	real zeta
	real z(3)
	integer ie,i

	name='sev23122012_ore18.txt'
	open(1,file=name,status='old',form='formatted')

	i = 0
	ie = 0
    1	continue
	  read(1,*,end=2) y,x,aux,ssl
	  call find_elem_from_old(ie,x,y,ie)
	  zeta = -999.
	  if( ie .gt. 0 .and. ssl .lt. 10. .and. y .ge. 40. ) then
	    call get_scal_elem(ie,znv,z)
	    call femintp(ie,z,x,y,zeta)
	    i = i + 1
	    write(66,*) i,ssl,zeta
	    write(67,*) x,y,ssl,ie,zeta
	  end if
	  write(6,*) x,y,ssl,ie,zeta
	goto 1
    2	continue

	close(1)

	end

c******************************************************************

