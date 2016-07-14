	program pdbcalc !calculates umatrix, hmatrix, input
	integer i,j,k,anum,rnum,ntot,ios,nlines,maxlines,nres,natoms
	integer start,stop,nr,is,imin,jmin
	real occ,bfac,um,dotpij,rij,mrad,bl
	real Rb,T,cc,gc
	character(16)crap,atype,restype,chain,aat,protname
	character(128)line
	real, dimension(:), allocatable :: x,y,z,rx,ry,rz,lx,ly,lz,lmag,nx,ny,nz,hx,hy,hz,nhx,nhy,nhz,nhmag
	real, dimension(:,:), allocatable :: ldot
	double precision, dimension(:,:), allocatable :: qm,msf
	double precision, dimension(:), allocatable :: eig
	DATA P/3.141592654/
	!get number of lines, starting point
	open(unit=2,file="protname.txt",status='old')
	protname=adjustl(protname)
	read(2,*)protname
	close(2)
	open(unit=2,file="nres.dat",status='old')
	read(2,*)nres
	close(2)
	open(unit=2,file="natoms.dat",status='old')
	read(2,*)natoms
	close(2)

	ALLOCATE ( x(natoms),y(natoms),z(natoms),rx(nres),ry(nres),rz(nres),lx(nres),ly(nres),lz(nres),lmag(nres) )
	ALLOCATE ( eig(nres) )
	ALLOCATE ( nx(nres),ny(nres),nz(nres),hx(nres),hy(nres),hz(nres),nhx(nres),nhy(nres),nhz(nres),nhmag(nres) )
	ALLOCATE ( ldot(nres,nres) )
	ALLOCATE ( qm(nres,nres),msf(nres,nres) )

	open(unit=13,file=trim(protname)//".pdb")
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
	open(unit=13,file=trim(protname)//".pdb")
	open(unit=23,file="mrad.dat")
	open(unit=103,file="length")
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
	if(atype(1:2).eq."CA")then
!	write(*,*)line(32:38),x(i),line(40:46),y(i),line(48:54),z(i)
	ntot=ntot+1
	rx(ntot)=x(i)
	ry(ntot)=y(i)
	rz(ntot)=z(i)
	nx(ntot)=x(i-2)
	ny(ntot)=y(i-2)
	nz(ntot)=z(i-2)
	hx(ntot)=x(i-1)
	hy(ntot)=y(i-1)
	hz(ntot)=z(i-1)
	if(restype(1:3).eq."PRO")then
	nx(ntot)=x(i-10)
	ny(ntot)=y(i-10)
	nz(ntot)=z(i-10)
	hx(ntot)=x(i-9)
	hy(ntot)=y(i-9)
	hz(ntot)=z(i-9)
	end if
!	write(*,*)bfac
!	assign Van Der Waals radii from miller values
	if(restype(1:3).eq."ALA")mrad=(113.0/(4*3.14))**.5
	if(restype(1:3).eq."ARG")mrad=(241.0/(4*3.14))**.5
	if(restype(1:3).eq."ASN")mrad=(158.0/(4*3.14))**.5
	if(restype(1:3).eq."ASP")mrad=(151.0/(4*3.14))**.5
	if(restype(1:3).eq."CYS")mrad=(140.0/(4*3.14))**.5
	if(restype(1:3).eq."GLN")mrad=(189.0/(4*3.14))**.5
	if(restype(1:3).eq."GLU")mrad=(183.0/(4*3.14))**.5
	if(restype(1:3).eq."GLY")mrad=( 85.0/(4*3.14))**.5
	if(restype(1:3).eq."HIS")mrad=(194.0/(4*3.14))**.5
	if(restype(1:3).eq."ILE")mrad=(182.0/(4*3.14))**.5
	if(restype(1:3).eq."LEU")mrad=(180.0/(4*3.14))**.5
	if(restype(1:3).eq."LYS")mrad=(211.0/(4*3.14))**.5
	if(restype(1:3).eq."MET")mrad=(204.0/(4*3.14))**.5
	if(restype(1:3).eq."PHE")mrad=(218.0/(4*3.14))**.5
	if(restype(1:3).eq."PRO")mrad=(143.0/(4*3.14))**.5
	if(restype(1:3).eq."SER")mrad=(122.0/(4*3.14))**.5
	if(restype(1:3).eq."THR")mrad=(146.0/(4*3.14))**.5
	if(restype(1:3).eq."TRP")mrad=(259.0/(4*3.14))**.5
	if(restype(1:3).eq."TYR")mrad=(229.0/(4*3.14))**.5
	if(restype(1:3).eq."VAL")mrad=(160.0/(4*3.14))**.5
	write(23,*)mrad
	end if
	end do
	close(13)
!	calculate bond lengths
	lmag=0.0
	do j=1,ntot-1
	lx(j)=rx(j+1)-rx(j)
	ly(j)=ry(j+1)-ry(j)
	lz(j)=rz(j+1)-rz(j)
	lmag(j)=(lx(j)**2+ly(j)**2+lz(j)**2)**.5
	nhx(j)=hx(j+1)-nx(j+1)
	nhy(j)=hy(j+1)-ny(j+1)
	nhz(j)=hz(j+1)-nz(j+1)
	nhmag(j)=(nhx(j)**2+nhy(j)**2+nhz(j)**2)**.5
	end do
	do i=1,ntot-1
	write(103,*)lmag(i)*.1
	end do
	close(103)

	ldot=0.0
	!calculate lNH dot lCA	
	do i=1, ntot-1
	do j=1, ntot-1
	ldot(i,j)=(nhx(i)*lx(j)+nhy(i)*ly(j)+nhz(i)*lz(j))/(lmag(j)*nhmag(i))
	end do
	end do

	open(unit=36,file='ldot.dat')
	do i=1,ntot-1
	do j=1,ntot-1
	write(36,*)ldot(i,j)
	end do
	end do

	close(36)

	open(unit=2,file="protname.txt",status='old',position='append')
	write(2,*)ntot
	close(2)

	!assign Rij
	rijmin=0.0
	rijp=100.0
	
	open(unit=15,file="Rij")
	do i=1,ntot
	do j=1,ntot
	if(i.ne.j)then
	rij=(rx(i)-rx(j))**2+(ry(i)-ry(j))**2+(rz(i)-rz(j))**2
	rij=rij**.5
	rij=rij*.1
	if(rij.le.rijp.and.rij.gt.0.0)then
	rijmin=rij
	imin=i
	jmin=j
	end if
	rijp=rij
	rij=1.0/rij
	end if
	if(i.eq.j)then
	rij=0.0
	end if
	write(15,*)rij
	end do
	end do
	close(15)
	open(unit=21,file="Rij_min")
	write(21,*)rijmin,imin,jmin
	close(21)
	bl=0.0
	!average bond length
	do i=1,ntot-1
	bl=bl+lmag(i)*.1
	end do
	close(18)
	bl=bl/real(ntot-1)
	open(unit=19,file='avbl')
	write(19,*)bl
	close(19)
	close(23)

	!calculate neighbor list
	call contactm(ntot)
!	cmi=0.0
!	open(unit=10,file='gammaI.dat',status='old')
!	do i=1,ntot
!	do j=1,ntot
!	read(10,*)cmi(i,j)
!	end do
!	end do
!	close(10)
	open(unit=10,file='Qgamma.dat')
	do i=1,ntot
	do j=1,ntot
	read(10,*)qm(i,j)
	end do
	end do
	close(10)

	open(unit=10,file='gammaeig.dat')
	do i=1,ntot
	read(10,*)eig(i)
	end do
	close(10)
	Rb=.00198 !(boltzmanns constant in kcal/mol*K)
	T=300.0
	gc=6. !potential strength adjustable 100-150 from bahar 2001
	cc=(3*Rb*T)/(gc)
!	write(*,*)"constant:",cc
!	write(*,*)(8./3.)*(P**2)*cc

	msf=0.0
	do i=1,ntot
	do j=1,ntot
	do k=2,ntot !throw out zero eigenvalue
	msf(i,j)=msf(i,j)+(qm(i,k)*qm(j,k))/eig(k)
	end do
	end do	
	end do

	!bfactors
	open(unit=11,file='bfactors.dat')
	do i=1,ntot
	write(11,*)i,(8./3.)*(P**2)*cc*msf(i,i)*100.
	end do
	close(11)
	

!	calculate Umatrix
	open(unit=16,file='Umatrix')
	um=0.0
	dotpij=0.0
	write(16,*)ntot-1
	do i=1,ntot-1
	do j=1,ntot-1
	dotpij=lx(i)*lx(j)+ly(i)*ly(j)+lz(i)*lz(j)
	dotpij=(dotpij+cc*msf(i,j)+cc*msf(i+1,j+1)-cc*msf(i,j+1)-cc*msf(i+1,j))/(lmag(i)*lmag(j))
	um=dotpij
	write(16,*)um
!	write(*,*)cc*msf(i,j), cc*msf(i,j)+cc*msf(i+1,j+1)-cc*msf(i,j+1)-cc*msf(i+1,j)
	end do
	end do
	close(16)



	

	end program
	
	subroutine contactm(ntot)
	integer i,j,k,io
	integer,dimension(ntot) :: csum,ipiv
	integer,dimension(ntot,ntot) :: neighlist
	double precision,dimension(ntot,ntot) :: rij,cm,cmi,ctest
	double precision rc,work(ntot*ntot),eigc(ntot)
	rc=.7 !cutoff in nm
	write(*,*)ntot
	write(*,*)"using cutoff of:",rc,"nm"
	open(unit=15,file='Rij',status='old')
!	read(15,*)
	nloc=0
	rij=0.0
	cm=0.0
	cmi=0.0
	do i=1,ntot
	do j=1,ntot
	read(15,*)rij(i,j)
	if(rij(i,j).eq.0.0)rij(i,j)=.001
	rij(i,j)=1.0/rij(i,j)
	end do	
	end do

	!calculate gamma contact matrix
	!start with contacts for off-diagonals
	do i=1,ntot
	do j=1,ntot
	if(i.ne.j)then
	 if(rij(i,j).le.rc)then
	 cm(i,j)=-1.0
!	 write(*,*)"contact!",i,j
	 end if
	end if
	end do
	end do

	csum=0.0
	!get sums for diagonals
	do i=1,ntot
	do j=1,ntot
	if(i.ne.j)csum(i)=csum(i)-cm(i,j)
	end do
	end do
	!assign diagonals
	do i=1,ntot
	cm(i,i)=csum(i)
	end do

!	cmi=cm

	!check
!	do i=1,ntot
!	write(*,*)
!	do j=1,ntot
!	if(cm(i,j).ne.0.0)write(*,*)cm(i,j)
!	end do
!	end do

!	call DGETRF(ntot, ntot, cmi, ntot, ipiv, io)
!	write(*,*)"sgetrf (0 is success):",io
!	call DGETRI(ntot, cmi, ntot, ipiv, work, ntot*ntot, io)
!	write(*,*)"sgetri: (0 is success)",io

	!diagonalize for eigen expansion
	eig=0.0
	call DSYEV("V","U",ntot,cm,ntot,eigc,work,ntot*ntot,io)
	write(*,*)"ssyev: (0 is success)",io,work(1)

	!check
!	ctest=0.0
!	ctest=matmul(cm,cmi)
!	do i=1,ntot
!	write(*,*)
!	do j=1,ntot
!	if(cmi(i,j).ne.0.0)write(*,*)cmi(i,j)
!	write(*,*)ctest(i,j)
!	end do
!	end do
	
	open(unit=10,file='Qgamma.dat')
	do i=1,ntot
	do j=1,ntot
	write(10,*)cm(i,j)
	end do
	end do
	close(10)

	open(unit=10,file='gammaeig.dat')
	do i=1,ntot
	write(10,*)eigc(i)
	end do
	close(10)

	end subroutine
	
	


