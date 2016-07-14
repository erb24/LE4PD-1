	program resarea
	implicit none
	integer i,j,k,n,nres
	real r,aa,sdev
	real, dimension(1000) :: ra
	n=0
	r=0.0
	aa=0.0
	sdev=0.0
	ra=0.0
	open(unit=23,file='protname.txt',status='old')
	read(23,*)
	read(23,*)n
	close(23)
	open(unit=1,file="resarea.xvg")
	do i=1,20
	read(1,*)
	end do
	do i=1,n
	read(1,*)nres,aa,sdev
	ra(i)=((aa/(4*3.1415927))**.5)*10.0
	end do
	close(1)
	open(unit=2,file="avresrad")
	do i=1,n
	write(2,*)ra(i)
	end do 
	end program
