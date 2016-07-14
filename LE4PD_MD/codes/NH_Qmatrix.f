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
c	nfrs=500 !testing!
	call umatrix(n,nca,nfrs)
	End program inputreader

	subroutine umatrix(n,nca,nfrs)
	integer i,j,k,l,m,imin,jmin,a,nca,n,nfrs
	real, dimension(n) :: rx,ry,rz,lavmsq,avlx,avly,avlz,avlmag
	real, dimension(n) :: lx,ly,lz,lmag,avlxca,avlyca,avlzca,avlmagca
	real, dimension(nca) :: rxca,ryca,rzca,lavmca
	real, dimension(nca) :: lxca,lyca,lzca,lmagca
	real,dimension(n) :: s1NH,avximsq,NHorder
	real,dimension(3,3) :: RM
	real,dimension(3,1) :: vc1,vc2,vdip1,vdip2
	real rnx,rny,rnz,dipl,rnl,theta,thetadot,dot,cross,rlx,rly,rlz
	real lcadot1,lcadot2,avblnhsq
	real, dimension(n,n) :: qinvm,qm,Aia,AiaCA,QmatrixNH
	real, dimension(n,nfrs) :: lix,liy,liz,limag
     &,dlnhx,dlnhy,dlnhz,dlnhmag
	real, dimension(nca,nfrs) :: lixca,liyca,lizca,limagca,
     &rlixca,rliyca,rlizca,xix,xiy,xiz,xim
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
	xix=0.0
	xiy=0.0
	xiz=0.0
	xim=0.0
	qinvm=0.0
	qm=0.0
	dipcorr=0.0
	Aia=0.0
	AiaCA=0.0
	QmatrixNH=0.0
	open(unit=21,file='QINVmatrix',status='old')
	do i=1,nca-1
	do j=1,nca-1
	read(21,*)qinvm(i,j)
	end do
	end do
	close(21)
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
	avximsq=0.0
	!read from trajectory
	open(unit=11,file='nitro.g96',status='old')
	open(unit=12,file='hydro.g96',status='old')
	open(unit=13,file=trim(protname)//'.g96',status='old')
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
	if(mod(k,10000).eq.0)write(*,*)"reading frame",k
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

!	calculating instantaneous mode vector
	do a=1,nca-1 !mode loop
	do j=1,nca-1 !residue loop
	xix(a,k)=qinvm(a,j)*lxca(j)+xix(a,k)
	xiy(a,k)=qinvm(a,j)*lyca(j)+xiy(a,k)
	xiz(a,k)=qinvm(a,j)*lzca(j)+xiz(a,k)
	xim(a,k)=((xix(a,k)**2+xiy(a,k)**2+xiz(a,k)**2)**.5)
	end do
!	xix(a,k)=xix(a,k)/(xim(a,k)) !unit vector needed
!	xiy(a,k)=xiy(a,k)/(xim(a,k))
!	xiz(a,k)=xiz(a,k)/(xim(a,k))
c	write(*,*)"norm mode:",((xix(a,k)**2+xiy(a,k)**2+xiz(a,k)**2)**.5)
	avximsq(a)=avximsq(a)+xim(a,k)**2
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

	do j=1,nca-1 !ca
	avlxca(j)=avlxca(j)/real(nfrs)
	avlyca(j)=avlyca(j)/real(nfrs)
	avlzca(j)=avlzca(j)/real(nfrs) 
	avlmagca(j)=(avlxca(j)**2+avlyca(j)**2+avlzca(j)**2)**.5
c	write(*,*)j,avlmagca(j)
	avximsq(j)=avximsq(j)/real(nfrs)
	end do

	QmatrixNH=0.0
	do i=3,n-1,2
	l=(i-1)/2
	write(*,*)"bond",i,l
	write(ii,*)l
	ii=adjustl(ii)
	do k=1,nfrs !time loop

	do a=1,nca-1 !loop over modes to obtain QmatrixNH as <lNH dot xi_a>/<xi_a^2>
	QmatrixNH(l,a)=QmatrixNH(l,a)+((lix(i,k))*xix(a,k)+
     &(liy(i,k))*xiy(a,k)
     &+(liz(i,k))*xiz(a,k))
	end do !end mode loop
	
	end do !end time loop
	end do !end residue loop

	do i=1,nca-1 !normalize
	do a=1,nca-1
	QmatrixNH(i,a)=QmatrixNH(i,a)/real(nfrs)
	QmatrixNH(i,a)=QmatrixNH(i,a)/avximsq(a)
	end do
	end do

	open(unit=1,file="Qmatrix_NH")
	do i=1,nca-1
	do a=1,nca-1
	write(1,*)QmatrixNH(i,a)
	end do
	end do
	close(1)

	open(unit=1,file="ximsq")
	do i=1,nca-1
	write(1,*)avximsq(i)
	end do
	close(1)

	avblsqNH=0.0
	k=0
	open(unit=2,file='blsqNH')
	open(unit=3,file="NHorder")
	do j=3,n-1
	if(mod(j,2).eq.1)then
	avblsqNH=lavmsq(j)+avblsqNH
	NHorder(j)=(avlmag(j)**2)/lavmsq(j)
	k=k+1
	write(2,*)lavmsq(j)
	write(3,*)NHorder(j)
	end if	
	end do
	close(2)
	close(3)
	avblsqNH=avblsqNH/real(k)
	open(unit=1,file="avblsqNH")
	write(1,*)avblsqNH
	close(1)

	end subroutine

