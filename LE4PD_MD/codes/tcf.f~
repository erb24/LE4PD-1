	program inputreader
	integer nfrs,n
	open(unit=5,file="protname.txt",status='old')
	read(5,*)
	read(5,*)n
	read(5,*)nfrs
	close(5)
	nca=n
	n=2*n
	write(*,*)n,nca,nfrs
!	nfrs=500 !testing!
	call umatrix(n,nca,nfrs)
	End program inputreader

	subroutine umatrix(n,nca,nfrs)
	integer i,j,k,l,m,imin,jmin,a,nca,n,nfrs
	real, dimension(n) :: rx,ry,rz,lavmsq,avlx,avly,avlz,avlmag
	real, dimension(n) :: lx,ly,lz,lmag,avlxca,avlyca,avlzca,avlmagca
	real, dimension(nca) :: rxca,ryca,rzca,lavmca
	real, dimension(nca) :: lxca,lyca,lzca,lmagca,lavmsqca
	real,dimension(n) :: s1NH
	real,dimension(3,3) :: RM
	real,dimension(3,1) :: vc1,vc2,vdip1,vdip2
	real rnx,rny,rnz,dipl,rnl,theta,thetadot,dot,cross,rlx,rly,rlz
	real lcadot1,lcadot2
	real, dimension(n,n) :: qinvm,qm
	real, dimension(n,nfrs) :: lix,liy,liz,limag
     &,dlnhx,dlnhy,dlnhz,dlnhmag
	real, dimension(nca,nfrs) :: lixca,liyca,lizca,limagca,
     &rlixca,rliyca,rlizca
	real, dimension(n,0:nfrs) :: dipcorr,dipcorr2,dipcorr3
	character(32)protname
	character(16)aa,ii
	nca=n/2 !number of calphas
	write(*,*)"no. of ca:",nca
	lix=0.0
	liy=0.0
	liz=0.0
	limag=0.0
	lixca=0.0
	liyca=0.0
	lizca=0.0
	limagca=0.0
	qinvm=0.0
	qm=0.0
	dipcorr=0.0
	open(unit=5,file="protname.txt",status='old')
	read(5,'(A)')protname
	close(5)
	protname=adjustl(protname)
	rx=0.0
	ry=0.0
	rz=0.0
	rxca=0.0
	ryca=0.0
	rzca=0.0
	lx=0.0
	ly=0.0
	lz=0.0
	lxca=0.0
	lyca=0.0
	lzca=0.0
	lmag=0.0
	lmagca=0.0
	lavm=0.0
	s1NH=0.0
	rnx=0.0
	rny=0.0
	rnz=0.0
	rnl=0.0
	RM=0.0
	vcc=0.0
	vdip=0.0
	avrnhin3=0.0
	thetadot=0.0
	avlx=0.0
	avly=0.0
	avlz=0.0
	avlmag=0.0
	avlxca=0.0
	avlyca=0.0
	avlzca=0.0
	avlmagca=0.0
	dlnhx=0.0
	dlnhy=0.0
	dlnhz=0.0
	dlnhmag=0.0
	lavmsq=0.0
	rlixca=0.0
	rliyca=0.0
	rlizca=0.0
	lavmsqca=0.0
	!read from trajectory
	open(unit=11,file='nitro.g96',status='old')
	open(unit=12,file='hydro.g96',status='old')
	open(unit=13,file=trim(protname)//'_rot.g96',status='old')
	!skip first 7,now read and calculate stuff
	do i=1,7
	read(11,*)
	read(12,*)
	read(13,*)
	end do

	do k=1,nfrs
	rx=0.0
	ry=0.0
	rz=0.0
	rxca=0.0
	ryca=0.0
	rzca=0.0
	lx=0.0
	ly=0.0
	lz=0.0
	lxca=0.0
	lyca=0.0
	lzca=0.0
	lmag=0.0
	lmagca=0.0
	if(mod(k,100).eq.0)write(*,*)"frame",k
	do j=1,n
c	write(*,*)j
	if(mod(j,2).eq.1)then
	read(11,*)rx(j),ry(j),rz(j)
	end if
	if(mod(j,2).eq.0)then
	read(12,*)rx(j),ry(j),rz(j)
	end if
	end do

	do j=1,nca !ca
	read(13,*)rxca(j),ryca(j),rzca(j)
	end do

	do j=1,n-1
	if(mod(j,2).eq.0)then
	lx(j)=rx(j+1)-rx(j-1)
	ly(j)=ry(j+1)-ry(j-1)
	lz(j)=rz(j+1)-rz(j-1)
	lmag(j)=(lx(j)**2+ly(j)**2+lz(j)**2)**.5
c	write(*,*)"NN:",lmag(j)
	end if
	if(mod(j,2).eq.1)then
	lx(j)=rx(j)-rx(j+1)
	ly(j)=ry(j)-ry(j+1)
	lz(j)=rz(j)-rz(j+1)
	lmag(j)=(lx(j)**2+ly(j)**2+lz(j)**2)**.5
c	write(*,*)"NH:",lmag(j)
	avlx(j)=avlx(j)+lx(j)
	avly(j)=avly(j)+ly(j)
	avlz(j)=avlz(j)+lz(j)
	lavmsq(j)=lavmsq(j)+lmag(j)**2
	end if	
	end do

	do j=1,nca-1 !ca
	lxca(j)=rxca(j+1)-rxca(j)
	lyca(j)=ryca(j+1)-ryca(j)
	lzca(j)=rzca(j+1)-rzca(j)
	lmagca(j)=(lxca(j)**2+lyca(j)**2+lzca(j)**2)**.5
c	write(*,*)"ca:",lmagca(j)
	avlxca(j)=avlxca(j)+lxca(j)
	avlyca(j)=avlyca(j)+lyca(j)
	avlzca(j)=avlzca(j)+lzca(j)
	lavmsqca(j)=lavmsqca(j)+lmagca(j)**2
	end do

	!skip 8 lines
	do j=1,8
	read(11,*)
	read(12,*)
	read(13,*)
	end do
!	put into array
	do j=1,n-1 !bond loop
	lix(j,k)=lx(j)
	liy(j,k)=ly(j)
	liz(j,k)=lz(j)
	limag(j,k)=lmag(j)
	end do
	do j=1,nca-1 !ca
	lixca(j,k)=lxca(j)
	liyca(j,k)=lyca(j)
	lizca(j,k)=lzca(j)
	limagca(j,k)=lmagca(j)
	end do


	!come out of time loop
	end do
	close(11)
	close(12)
	close(13)

	!normalize
	do i=1,n-1
	l=(i+1)/2-1
	if(mod(i,2).eq.1)then
	avlx(i)=avlx(i)/real(nfrs)
	avly(i)=avly(i)/real(nfrs)
	avlz(i)=avlz(i)/real(nfrs) 
	lavmsq(i)=lavmsq(i)/real(nfrs) 
	avlmag(i)=(avlx(i)**2+avly(i)**2+avlz(i)**2)**.5
	write(*,*)l,lavmsq(i)
	end if
	end do

	open(unit=20,file='bondorder_sim.dat')
	do j=1,nca-1 !ca
	avlxca(j)=avlxca(j)/real(nfrs)
	avlyca(j)=avlyca(j)/real(nfrs)
	avlzca(j)=avlzca(j)/real(nfrs) 
	lavmsqca(j)=lavmsqca(j)/real(nfrs)
	avlmagca(j)=(avlxca(j)**2+avlyca(j)**2+avlzca(j)**2)**.5
	write(20,*)j,(avlmagca(j)**2)/(lavmsqca(j))
	end do
	close(20)

	!calculate nh-nh angle
c	do i=3,n-2,2
c	do i=1,1
c	l=(i-1)/2
c	write(*,*)"bond",i,l
c	write(ii,*)l
c	ii=adjustl(ii)
c	open(unit=200,file='NHNH_'//trim(ii))
c	do k=1,nfrs
c	write(200,*)k*.2,acos((lix(i,k)*lix(i+1,k)+liy(i+1,k)*liy(i,k)
c     &+liz(i,k)*liz(i+1,k))/(limag(i,k)*limag(i+1,k)))*(360./6.28)
c	end do
c	end do
c	close(200)

	!calculate nh-ca angle
c	do i=3,n-1,2
c	do i=1,1
c	l=(i-1)/2
c	write(*,*)"bond",i,l
c	write(ii,*)l
c	ii=adjustl(ii)
c	open(unit=200,file='NHCA_'//trim(ii))
c	do k=1,nfrs
c	write(200,*)k*.2,acos((lix(i,k)*lixca(l,k)+liy(i,k)*liyca(l,k)
c     &+liz(i,k)*lizca(l,k))/(limag(i,k)*limagca(l,k)))*(360./6.28)
c	end do
c	end do
c	close(200)

	dipcorr=0.0
	dipcorr2=0.0
	dipcorr3=0.0
	s1NH=0.0
	!calculate tcf of bond autocorrelation
	do i=3,n-1,2
c	i=BONDNUMBER
	l=(i-1)/2
	write(*,*)"bond",i,l
	write(ii,*)l
	ii=adjustl(ii)
	open(unit=200,file='m1CAsim_'//trim(ii))
c	open(unit=201,file='dip123_'//trim(ii))
!	open(unit=202,file='m1NH_'//trim(ii))
!	do j=0,nfrs/4
	do j=0,20000
	if(mod(j,500).eq.0)write(*,*)i,"tcf frame:",j
	m=0
	do k=1,nfrs-j
c	write(*,*)i,"tcf frame:",j
!and now the tcf of the NH and the CA
	dipcorr(i,j)=dipcorr(i,j)+((lixca(l,k)*lixca(l,j+k)
     &+liyca(l,k)*liyca(l,j+k)+lizca(l,k)*lizca(l,j+k)))
     &/(lavmsqca(l))
!	dipcorr2(i,j)=dipcorr2(i,j)+((lix(i,k)*lix(i,j+k)
!     &+liy(i,k)*liy(i,j+k)+liz(i,k)*liz(i,j+k)))
!     &/(limag(i,k)*limag(i,j+k))
	m=m+1

	end do
	dipcorr(i,j)=dipcorr(i,j)/real(nfrs-j)
!	dipcorr2(i,j)=dipcorr2(i,j)/real(nfrs-j)
c	write(*,*)nfrs-j,m
	write(200,*)j*.2,dipcorr(i,j)
!	write(202,*)j*.2,dipcorr2(i,j)
c	write(*,*)j*.2,dipcorr(i,j)
	end do
	close(200)
!	close(202)
	end do

	end subroutine

