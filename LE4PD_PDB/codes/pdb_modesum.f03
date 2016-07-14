	program pdbcalc !calculates umatrix, hmatrix, input
	integer i,j,k,anum,rnum,ntot,ios,nlines,maxlines,nres,natoms
	integer start,stop,nr,is,imin,jmin,nc,nm1,a,icrap,nn
	real occ,bfac,um,dotpij,rij,mrad,bl,avfrw,avfr
	real Rb,T,cc,gc,blsq,sg,eps,eigmuscale,ap,amax
	character(16)crap,atype,restype,chain,aat,protname,cnc,cmamp,ctau,canum
	character(128)line
	real, dimension(:), allocatable :: x,y,z,rx,ry,rz
	real, dimension(:,:), allocatable :: qm,ma
	real, dimension(:), allocatable :: tau,ml,eigmu,s1,Am,eiglam,barrier,fri,friw
	integer, dimension(:), allocatable :: imax
	DATA P/3.141592654/
	!get number of lines, starting point
	open(unit=2,file="prot",status='old')
	protname=adjustl(protname)
	read(2,*)protname
	close(2)
	open(unit=2,file="nres.dat",status='old')
	read(2,*)nres
	close(2)
	open(unit=2,file="natoms.dat",status='old')
	read(2,*)natoms
	close(2)
	open(unit=2,file="ncopies.dat",status='old')
	read(2,*)nc
	close(2)
	nm1=nres-1
	Rb=.00198 !(boltzmanns constant in kcal/mol*K)
	eps=6.42 !energy per unit length in kcal/mol/nm
	open(unit=19,file='temp',status='old')
	read(19,*)T
	close(19)

	ALLOCATE ( x(natoms),y(natoms),z(natoms),rx(natoms),ry(natoms),rz(natoms) )
	ALLOCATE ( tau(nres-1),Am(nres),eiglam(nres),eigmu(nres),ml(nres),s1(nres),barrier(nres),imax(nres),fri(nres),friw(nres) )
	ALLOCATE ( ma(nres,nres),qm(nres,nres) )
	rx=0.0
	ry=0.0
	rz=0.0
	tau=0.0
	ma=0.0
	x=0.0
	y=0.0
	z=0.0
	imax=0
	open(unit=19,file='avbl',status='old')
	read(19,*)blsq
	close(19)
	blsq=blsq*blsq*1.002

	do k=1,nc !open copy loop
	write(cnc,*)k
	cnc=adjustl(cnc)
	open(unit=13,file=trim(protname)//trim(cnc)//".pdb")
	nlines=0
	ios=0
	maxlines=150000
	start=0
	is=0
	stop=0
	nr=0
	do i=1,maxlines
	Read(13,*,IOSTAT=ios)line
	If(ios /= 0) EXIT
	If (J == maxlines) then
	write(*,*) "Maxlines exceeded"
	STOP
	End If
	if(line(1:4) .eq. "ATOM" .AND. is .eq. 0)then
	start=nlines
	is=1
	end if
	if(is .eq. 1 .AND. line(1:3) .eq. "TER") EXIT
	nlines=nlines+1
	end do
	stop=nlines
	nr=stop-start

	write(*,*)start,stop,nr,nlines
	close(13)

!	now read in to grab values
	write(*,'(A)')trim(protname)//trim(cnc)//".pdb"
	open(unit=13,file=trim(protname)//trim(cnc)//".pdb")
	ntot=0
	mrad=0.0
	do i=1,start
	read(13,*)
	end do
	do i=1,nr
!	read(13,*)crap,anum,atype,restype,chain,rnum,x(i),y(i),z(i),
!     &occ,bfac,aat
	read(13,'(A)')line
	read(line(14:16),'(A)')atype
	read(line(18:20),'(A)')restype
	atype=adjustl(atype)
	restype=adjustl(restype)
!	write(*,*)atype,restype,line(14:16),line(18:20)
	read(line(32:38),'(F8.3)')x(i)
	read(line(40:46),'(F8.3)')y(i)
	read(line(48:54),'(F8.3)')z(i)
!	write(*,*)line(32:38),x(i),line(40:46),y(i),line(48:54),z(i)
	ntot=ntot+1
	rx(ntot)=rx(ntot)+x(i)
	ry(ntot)=ry(ntot)+y(i)
	rz(ntot)=rz(ntot)+z(i)
!	write(*,*)bfac

	end do
	close(13)

	end do !end copy loop

	do i=1,natoms !normalize average structure
	rx(i)=rx(i)/real(nc)
	ry(i)=ry(i)/real(nc)
	rz(i)=rz(i)/real(nc)
!	write(*,*)i,rx(i),ry(i),rz(i)
	end do

	!calculate mode amplitudes as integral over normalized ma and tau_a
	ma=0.0
	ml=0.0
	qm=0.0
	eigmu=0.0
	tau=0.0

	open(unit=99,file="Qmatrix",status='old')
	do i=1,nm1
	do j=1,nm1
	read(99,*)qm(i,j)
!	write(*,*)qm(i,j)
	end do
	end do
	close(99)
	open(unit=99,file="mu_eig",status='old')
	do i=1,nm1
	read(99,*)eigmu(i)
!	write(*,*)eigmu(i)
	end do
	open(unit=99,file="lambda_eig",status='old') !read in all values
	do i=1,nm1
	read(99,*)eiglam(i)
!	write(*,*)eiglam(i)
	end do
	open(unit=747,file='sigma',status='old')
	read(747,*) 
	read(747,*) sg ! sigma
	close(747)
	open(unit=10,file='barriers.dat')
	do i=1,nm1
	read(10,*)icrap,barrier(i)
	end do
	close(10)



	do i=1,nm1
	ap=0.0
	do a=4,nm1
	ma(i,a)=(qm(i,a)**2.)/eigmu(a)
	s1(i)=s1(i)+ma(i,a)
	if(ma(i,a).ge.ap)then
	amax=ma(i,a)
	imax(i)=a
	ap=amax
	end if
	end do
	ml(i)=(s1(i)*blsq)**.5
!	write(*,*)s1(i)
!	write(*,*)imax(i)
	end do

	do i=1,nm1
	do a=4,nm1
	tau(i)=tau(i)+(ma(i,a)/s1(i))*((1./(eiglam(a)*sg))*(exp(barrier(a)/(Rb*T))))
	end do
	end do

	open(unit=10,file='ml')
	open(unit=11,file='taum1')
	do i=1,nm1
	write(10,*)ml(i)
	write(11,*)i,tau(i)
	end do

	open(unit=40,file='fric') !friction coefficients
	read(40,*)avfr,avfrw
	do i=1,nres
	read(40,*)friw(i),fri(i)
	end do
	close(40)

	open(unit=13,file=trim(protname)//trim(cnc)//".pdb")
	open(unit=14,file=trim(protname)//"_ms.pdb")
	nn=0
	do i=1,start
	Read(13,*)
	end do
	do i=1,nr-1
	Read(13,'(A)')line
!	write(*,'(A)')line
!	canum=line(8:12)
!	read(canum,'(I4)')anum
!	write(*,*)anum
	if(line(13:15) .eq. " N ")then
	nn=nn+1
	end if
	if(nn.eq.1)write(cmamp,'(F4.2)')ml(nn)*10.0
	if(nn.ne.1)write(cmamp,'(F4.2)')ml(nn-1)*10.0
	line(62:66)=cmamp(1:4)
	if(nn.eq.1)write(cmamp,'(F4.1)')(tau(nn)/maxval(tau))
	if(nn.ne.1)write(cmamp,'(F4.1)')(tau(nn-1)/maxval(tau))
!	write(cmamp,'(F4.2)')friw(nres)/fri(nres)
	write(line(32:38),'(F7.3)')rx(i)
	write(line(40:46),'(F7.3)')ry(i)
	write(line(48:54),'(F7.3)')rz(i)
	line(57:60)=cmamp(1:4)
	write(14,'(A)')line
	end do

	close(13)
	close(14)

	end program
	
