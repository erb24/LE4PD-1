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
!	nfrs=1000 !testing!
	call umatrix(n,nca,nfrs)
	End program inputreader

	subroutine umatrix(n,nca,nfrs)
	integer i,j,k,l,m,imin,jmin,a,nca,n,nfrs,itime(11),ip,nsims
	real, dimension(nca) :: rxca,ryca,rzca,lavmca
	real, dimension(n,0:nfrs) :: dipcorr,rixca,riyca,rizca
	real, dimension(nfrs) :: comx,comy,comz
	character(32)protname
	character(16)aa,ii
	nca=n/2 !number of calphas
	write(*,*)"no. of ca:",nca
	lix=0.0
	liy=0.0
	liz=0.0
	limag=0.0
	rixca=0.0
	riyca=0.0
	rizca=0.0
	limagca=0.0
	qinvm=0.0
	qm=0.0
	dipcorr=0.0
	nsims=10
c	open(unit=5,file="frames.dat")
c	read(5,*)
c	do i=1,nsims
c	read(5,*)itime(i),a
c	end do
c	close(5)
c	itime(nsims+1)=nfrs
c	write(*,*)itime(1),itime(2),itime(3),itime(4),itime(5),
c    &itime(6),itime(7),itime(8),itime(9),itime(10)
	open(unit=5,file="protname.txt",status='old')
	read(5,'(A)')protname
	close(5)
	rxca=0.0
	ryca=0.0
	rzca=0.0
	comx=0.0
	comy=0.0
	comz=0.0
	!read from trajectory
	open(unit=13,file=trim(protname)//'.g96',status='old')
	!skip first 7,now read and calculate stuff
	do i=1,7
	read(13,*)
	end do

	do k=1,nfrs
	rxca=0.0
	ryca=0.0
	rzca=0.0
	if(mod(k,100).eq.0)write(*,*)"frame",k

	do j=1,nca !ca
	read(13,*)rixca(j,k),riyca(j,k),rizca(j,k)
	end do

!	do i=1,nca
!	comx(k)=comx(k)+rxca(i)
!	comy(k)=comy(k)+ryca(i)
!	comz(k)=comz(k)+rzca(i)
!	end do
!	comx(k)=comx(k)/real(nca)
!	comy(k)=comy(k)/real(nca)
!	comz(k)=comz(k)/real(nca)

	!skip 8 lines
	do j=1,8
	read(13,*)
	end do


	!come out of time loop
	end do
	close(11)
	close(12)
	close(13)


	dipcorr=0.0
	open(unit=202,file='smsd')
	ip=1
c	do j=0,nfrs/4
	do j=0,nfrs/2
c	do ip=1,nsims
	if(mod(j,500).eq.0)write(*,*)i,"tcf frame:",j
	m=0
c	do k=1,nfrs-j
c	do k=itime(ip),itime(ip+1)-j
	do k=1,nfrs-j
c	write(*,*)i,"tcf frame:",j
!and now the tcf
	do l=1,nca
	dipcorr(ip,j)=dipcorr(ip,j)+((rixca(l,k)-rixca(l,j+k))**2.
     &+(riyca(l,k)-riyca(l,j+k))**2.+(rizca(l,k)-rizca(l,j+k))**2.)
	end do
	end do !end inner time for k
	dipcorr(ip,j)=dipcorr(ip,j)/real(nfrs-j)
	dipcorr(ip,j)=dipcorr(ip,j)/real(nca)
c	end do !end ip loop
c	do ip=1,nsims
c	dipcorr(20,j)=dipcorr(20,j)+dipcorr(ip,j)
c	end do
c	dipcorr(20,j)=dipcorr(20,j)/real(nsims)
	write(202,*)j*.2,dipcorr(ip,j)
c	write(*,*)j*.2,dipcorr(i,j)
	end do !end j loop
	close(202)

	end subroutine

