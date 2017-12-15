!
! tree utility routines
!
! still to do: implement tree for string
!
!===============================================================
	module tree
!===============================================================

	implicit none

	private

        type :: node
          integer :: parent
          integer :: left
          integer :: right
          integer :: key
          integer :: info
        end type node

        type :: entry
          integer :: root
          integer :: ndim
          integer :: nfill
          integer :: nfree
          type(node), allocatable :: data(:)
          integer, allocatable :: free(:)
        end type entry

        integer, save :: idlast = 0
        integer, save :: ndim = 0
        integer, parameter :: ndim_first = 5
        integer, parameter :: ndim_data = 32
        type(entry), save, allocatable, target :: pentry(:)

	public :: tree_init		!call tree_init(id)
	public :: tree_delete		!call tree_delete(id)

	public :: tree_insert		!call tree_push(id,value)
	public :: tree_walk		!logical tree_pop(id,value)
	public :: tree_get_data		!logical tree_pop(id,value)
	public :: tree_check		!logical tree_pop(id,value)

	public :: tree_is_empty		!logical tree_is_empty(id)
	public :: tree_info		!call tree_info(id)

!        INTERFACE tree_push
!        MODULE PROCEDURE         tree_push_d
!     +                          ,tree_push_r
!     +                          ,tree_push_i
!        END INTERFACE


!===============================================================
	contains
!===============================================================

        subroutine tree_init_alloc

        type(entry), allocatable :: paux(:)

        if( ndim == 0 ) then
          ndim = ndim_first
          allocate(pentry(ndim))
        else
          ndim = ndim*2
          allocate(paux(ndim))
          paux(1:ndim/2) = pentry(1:ndim/2)
          call move_alloc(paux,pentry)
        end if

        end subroutine tree_init_alloc

!******************************************************************

        subroutine tree_init_new_id(id)

        integer id

        idlast = idlast + 1
        if( idlast > ndim ) then
          call tree_init_alloc
        end if
        id = idlast

        call tree_init_id(id)

        end subroutine tree_init_new_id

!******************************************************************

        subroutine tree_init_id(id)

        integer id

        if( id > ndim ) then
          stop 'error stop tree_init_id: ndim'
        end if

        pentry(id)%root = 0
        pentry(id)%ndim = ndim_data
        pentry(id)%nfill = 0
        pentry(id)%nfree = 0

	if( allocated(pentry(id)%data) ) deallocate(pentry(id)%data)
	if( allocated(pentry(id)%free) ) deallocate(pentry(id)%free)

	allocate(pentry(id)%data(0:ndim_data))
	allocate(pentry(id)%free(ndim_data))

	pentry(id)%data(:) = node(0,0,0,0,0)
	pentry(id)%free(:) = 0

        end subroutine tree_init_id

!******************************************************************

	subroutine realloc_data(id)

	integer id

        type(node), allocatable :: daux(:)
        integer, allocatable :: faux(:)

	logical bmove,bdebug
	integer nsize,nsize2

	bmove = .true.
	bmove = .false.
	bdebug = .false.
	nsize = pentry(id)%ndim
	nsize2 = 2*nsize

	if( bdebug ) write(6,*) 'alloc: ',nsize,nsize2

	if( bmove ) then

        allocate(daux(0:nsize2))
        allocate(faux(nsize2))
        daux(0:nsize) = pentry(id)%data(0:nsize)
        faux(1:nsize) = pentry(id)%free(1:nsize)
        call move_alloc(daux,pentry(id)%data)
        call move_alloc(faux,pentry(id)%free)

	else

        allocate(daux(0:nsize))
        allocate(faux(nsize))
        daux(0:nsize) = pentry(id)%data(0:nsize)
        faux(1:nsize) = pentry(id)%free(1:nsize)
	deallocate(pentry(id)%data)
	deallocate(pentry(id)%free)
	allocate(pentry(id)%data(0:nsize2))
	allocate(pentry(id)%free(1:nsize2))
	pentry(id)%data(0:nsize2) = node(0,0,0,0,0)
	pentry(id)%free(1:nsize2) = 0
        pentry(id)%data(0:nsize) = daux(0:nsize)
        pentry(id)%free(1:nsize) = faux(1:nsize)
	deallocate(daux)
	deallocate(faux)

	end if

	pentry(id)%ndim = nsize2

	if( bdebug ) then
	write(6,*) '------------ alloc --------------'
	write(6,*) pentry(id)%data
	write(6,*) pentry(id)%free
	write(6,*) '------------ alloc --------------'
	end if

	end subroutine realloc_data

!******************************************************************

	integer function tree_get_free_node(id)

	integer id

	if( pentry(id)%nfree > 0 ) then
	  tree_get_free_node = pentry(id)%free(pentry(id)%nfree)
	  pentry(id)%nfree = pentry(id)%nfree - 1
	else 
	  if( pentry(id)%nfill == pentry(id)%ndim ) then
	    call realloc_data(id)
	  end if
	  pentry(id)%nfill = pentry(id)%nfill + 1
	  tree_get_free_node = pentry(id)%nfill
	end if

	end function tree_get_free_node

!******************************************************************

	subroutine tree_error(id,error)

	integer id,error

	!if( error == empty_error ) then
	!  write(6,*) 'tree: ',id
	!  stop 'error stop tree: tree is empty'
	!else if( error == type_error ) then
	!  write(6,*) 'tree: ',id
	!  write(6,*) 'type: ',pentry(id)%type
	!  stop 'error stop tree: variable is of wrong type'
	!else
	!  stop 'error stop tree: internal error (1)'
	!end if

	end subroutine tree_error

!******************************************************************
!******************************************************************
!******************************************************************

	subroutine tree_init(id)
	integer id
        call tree_init_new_id(id)
	end subroutine tree_init

	subroutine tree_delete(id)
	integer id
        call tree_init_id(id)
	if( id == idlast ) idlast = idlast - 1
	end subroutine tree_delete

!--------------------

	subroutine tree_get_data(id,x,key,info)

	integer id,x,key,info

	key = pentry(id)%data(x)%key
	info = pentry(id)%data(x)%info

	end subroutine tree_get_data

!--------------------

	subroutine tree_nodes_info(id,i)

	integer id,i

	write(6,*) 'tree_nodes_info: node'
	call tree_node_info(id,i)
	write(6,*) 'tree_nodes_info: parent'
	call tree_node_info(id,pentry(id)%data(i)%parent)
	write(6,*) 'tree_nodes_info: left'
	call tree_node_info(id,pentry(id)%data(i)%left)
	write(6,*) 'tree_nodes_info: right'
	call tree_node_info(id,pentry(id)%data(i)%right)

	end subroutine tree_nodes_info

	subroutine tree_node_info(id,i)

	integer id,i

	write(6,*) 'tree_node_info: ',i
	write(6,*) 'parent: ',pentry(id)%data(i)%parent
	write(6,*) 'left  : ',pentry(id)%data(i)%left
	write(6,*) 'right : ',pentry(id)%data(i)%right
	write(6,*) 'key   : ',pentry(id)%data(i)%key

	end subroutine tree_node_info

!--------------------

	subroutine tree_check(id)

	integer id
	integer i,n,key,p,l,r,x
        type(node), pointer :: data(:)

	data => pentry(id)%data

	n = pentry(id)%nfill

	do i=1,n
	  key = data(i)%key
	  p = data(i)%parent
	  r = data(i)%right
	  l = data(i)%left
	  if( p == 0 .and. r == 0 .and. l == 0 ) cycle
	  if( p > 0 ) then
	    if( data(p)%key > key ) then
	      if( data(p)%left /= i ) goto 99
	    else
	      if( data(p)%right /= i ) goto 99
	    end if
	  end if
	  if( l > 0 .and. data(l)%parent /= i ) goto 99
	  if( r > 0 .and. data(r)%parent /= i ) goto 99
	end do

	x = tree_walk(id,0)
	key = data(x)%key
	do
	  if( x == 0 ) exit
	  key = data(x)%key
	  p = data(x)%parent
	  l = data(x)%left
	  r = data(x)%right
	  if( l /= 0 .and. data(l)%parent /= x ) goto 98
	  if( l /= 0 .and. data(l)%key >= key ) goto 98
	  if( r /= 0 .and. data(r)%parent /= x ) goto 98
	  if( r /= 0 .and. data(r)%key <= key ) goto 98
	  x = tree_walk(id,x)
	end do

	return
   98	continue
	call tree_nodes_info(id,x)
	stop 'error stop tree_check: inconsistency 2'
   99	continue
	call tree_nodes_info(id,i)
	stop 'error stop tree_check: inconsistency 1'
	end subroutine tree_check

!--------------------

	integer function tree_walk(id,x)

	integer id,x

	if( x == 0 ) then
	  tree_walk = tree_minimum(id,pentry(id)%root)
	else
	  tree_walk = tree_successor(id,x)
	end if

	end function tree_walk

!--------------------

	integer function tree_successor(id,z)

	integer id,z

	integer x,y
        type(node), pointer :: data(:)
	
	data => pentry(id)%data
	x = z

	if( data(x)%right /= 0 ) then
	  tree_successor = tree_minimum(id,data(x)%right)
	else
	  y = data(x)%parent
	  do
	    if( y == 0 .or. x /= data(y)%right ) exit
	    x = y
	    y = data(y)%parent
	  end do
	  tree_successor = y
	end if

	end function tree_successor

!--------------------

	integer function tree_search(id,y,key)

	integer id,y,key
        type(node), pointer :: data(:)
	integer x
	
	x = y
	if( x == 0 ) x = pentry(id)%root
	data => pentry(id)%data

	do
	  if( x == 0 .or. data(x)%key == key ) exit
	  if( key < data(x)%key ) then
	    x = data(x)%left
	  else
	    x = data(x)%right
	  end if
	end do

	tree_search = x

	end function tree_search

!--------------------

	integer function tree_insert(id,key,info)

	integer id,key,info
        type(node), pointer :: data(:)
	integer x,y,z
	logical bdebug
	
	bdebug = .false.
	bdebug = (key == 5 )

	data => pentry(id)%data

	if( bdebug ) write(6,*) 'looking for: ',key
	y = 0
	x = pentry(id)%root
	do
	  if( x == 0 ) exit
	  y = x
	if( bdebug ) write(6,*) 'iter: ',x,data(x)%key
     +			,data(x)%left,data(x)%right
	  if( key == data(x)%key ) then
	    tree_insert = 0
	    return
	  else if( key < data(x)%key ) then
	    x = data(x)%left
	  else
	    x = data(x)%right
	  end if
	end do

	if( bdebug ) write(6,*) 'debug: ',x,y,data(y)%key,key

	z = tree_get_free_node(id)

	if( bdebug ) then
	  write(6,*) 'new node: ',z
	  call tree_node_info(id,z)
	end if

	data(z)%parent = y
	data(z)%key = key
	data(z)%info = info

	if( y == 0 ) then
	  pentry(id)%root = z
	else 
	  if( key < data(y)%key ) then
	    data(y)%left = z
	  else
	    data(y)%right = z
	  end if
	end if

	if( bdebug ) then
	  write(6,*) 'debug: '
	  write(6,*) 'node y: '
	  call tree_node_info(id,y)
	  write(6,*) 'node z: '
	  call tree_node_info(id,z)
	end if

	tree_insert = z

	end function tree_insert

!--------------------

	integer function tree_minimum(id,y)

	integer id,y
        type(node), pointer :: data(:)
	integer x

	x = y
	if( x == 0 ) x = pentry(id)%root
	data => pentry(id)%data

	do
	  if( data(x)%left == 0 ) exit
	  x = data(x)%left
	end do

	tree_minimum = x

	end function tree_minimum

!--------------------

	integer function tree_maximum(id,y)

	integer id,y
        type(node), pointer :: data(:)
	integer x

	x = y
	if( x == 0 ) x = pentry(id)%root
	data => pentry(id)%data

	do
	  if( data(x)%right == 0 ) exit
	  x = data(x)%right
	end do

	tree_maximum = x

	end function tree_maximum

!--------------------

	logical function tree_is_empty(id)
	integer id
	tree_is_empty = ( pentry(id)%root == 0 )
	end function tree_is_empty

!--------------------

	subroutine tree_info(id)
	integer id
	write(6,*) 'tree_info: ',id
     +			,pentry(id)%root
     +			,pentry(id)%ndim
     +			,pentry(id)%nfill
     +			,pentry(id)%nfree
	end subroutine tree_info

!===============================================================
	end module
!===============================================================

	subroutine tree_test

	use tree

	implicit none

	integer, parameter :: max = 100
	!integer, parameter :: nloop = 10
	integer, parameter :: nloop = 10000
	integer nin,nl,x,key,id
	logical bdebug
	real r

	bdebug = .true.
	bdebug = .false.

	call tree_init(id)

	call random_seed
	nin = 0

	do nl=1,nloop
	  if( nin > (max*9)/10 ) exit
	  call tree_rand_int(1,max,key)
	  !write(6,*) 'new key to insert: ',key
	  x = tree_insert(id,key,2*key)
	  if( x /= 0 ) nin = nin + 1
	  !write(6,*) 'inserted: ',nl,nin,x,key
	  !if( mod(nl,10) == 0 ) call tree_inorder(id)
	  call tree_check(id)
	end do

	call tree_inorder(id)
	write(6,*) nl,nin
	call tree_delete(id)

	write(6,*) 'tree test successfully finished: ',nloop,nl

	end

!******************************************************************

	subroutine tree_inorder(id)

	use tree

	implicit none

	integer id

	integer x,key,info

	write(6,*) 'inorder walk start'

	x = 0
	do
	  x = tree_walk(id,x)
	  if( x == 0 ) exit
	  call tree_get_data(id,x,key,info)
	  write(6,*) x,key,info
	end do

	write(6,*) 'inorder walk end'

	end

!******************************************************************

        subroutine tree_assert(bcheck,text,id)
        use tree
        implicit none
        logical bcheck
        character*(*) text
        integer id
        if( .not. bcheck ) then
          write(6,*) 'tree_assertion: ',trim(text)
          call tree_info(id)
          stop 'assertion failed'
        end if
        end

	subroutine tree_rand_int(min,max,irand)

	implicit none
	integer min,max
	integer irand
	real r

	call random_number(r)
	irand = min + (max-min+1)*r

	end

!******************************************************************

	programme tree_main
	call tree_test
	end programme tree_main

!******************************************************************
