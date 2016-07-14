	program inputreader
	integer nfrs,n
	open(unit=5,file="protname.txt",status='old')
	read(5,*)
	read(5,*)n
	read(5,*)nfrs
	close(5)
	write(*,*)n,nfrs
!	nfrs=500 !testing!
	call umatrix(n,nfrs)
	End program inputreader

	subroutine umatrix(n,nfrs)
	integer i,j,k,l,m,imin,jmin,a,nca,n,nfrs
	real, dimension(n) :: rx,ry,rz,lavmsq
	real, dimension(n) :: lx,ly,lz
	real, dimension(n,nfrs) :: lix,liy,liz
	real, dimension(n,0:nfrs) :: dipcorr
	real, dimension(0:nfrs) :: avdipcorr
	character(32)protname
	character(16)aa,ii
	lix=0.0
	liy=0.0
	liz=0.0
	limag=0.0
	open(unit=5,file="protname.txt",status='old')
	read(5,'(A)')protname
	close(5)
	protname=adjustl(protname)
	rx=0.0
	ry=0.0
	rz=0.0
	lx=0.0
	ly=0.0
	lz=0.0
	lavmsq=0.0
	!read from trajectory
	open(unit=13,file=trim(protname)//'.g96',status='old')
	!skip first 7,now read and calculate stuff
	do i=1,7
	read(13,*)
	end do

	do k=1,nfrs
	rx=0.0
	ry=0.0
	rz=0.0
	lx=0.0
	ly=0.0
	lz=0.0
	if(mod(k,500).eq.0)write(*,*)"frame",k
	do j=1,n !ca
	read(13,*)rx(j),ry(j),rz(j)
	end do

	do j=1,n-1 !ca
	lx(j)=rx(j+1)-rx(j)
	ly(j)=ry(j+1)-ry(j)
	lz(j)=rz(j+1)-rz(j)
	lavmsq(j)=lavmsq(j)+lx(j)**2+ly(j)**2+lz(j)**2
	end do

	!skip 8 lines
	do j=1,8
	read(13,*)
	end do
!	put into array
	do j=1,n-1 !ca
	lix(j,k)=lx(j)
	liy(j,k)=ly(j)
	liz(j,k)=lz(j)
	end do


	!come out of time loop
	end do
	close(11)
	close(12)
	close(13)

	do j=1,n-1 !ca
	lavmsq(j)=lavmsq(j)/real(nfrs)
	end do

	dipcorr=0.0
	avdipcorr=0.0
	!calculate tcf of bond autocorrelation
	do l=1,n-1
c	i=BONDNUMBER
	write(ii,*)l
	ii=adjustl(ii)
	open(unit=200+l,file='m1CAsim_'//trim(ii))
	end do
	open(unit=100,file='avm1CAsim')
	do j=0,nfrs/4
	do l=1,n-1 !bond loop
	if(mod(j,100).eq.0)write(*,*)"tcf frame:",j
	do k=1,nfrs-j
!and now the tcf
	dipcorr(l,j)=dipcorr(l,j)+((lix(l,k)*lix(l,j+k)
     &+liy(l,k)*liy(l,j+k)+liz(l,k)*liz(l,j+k)))
     &/(lavmsq(l))
	end do
	dipcorr(l,j)=dipcorr(l,j)/real(nfrs-j)
	avdipcorr(j)=dipcorr(l,j)+avdipcorr(j)
	write(200+l,*)j*.2,dipcorr(l,j)
	end do !end bond loop
	avdipcorr(j)=avdipcorr(j)/real(n-1)
	write(100,*)j*.2,avdipcorr(j)
	end do !end j time loop
	close(100)
	do l=1,n-1
	close(200+l)
	end do


	end subroutine

