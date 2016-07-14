	program pdbcalc !calculates umatrix, hmatrix, input
	integer i,j,k,anum,rnum,ios,nlines,maxlines,ntot,HYDtot,NITtot,nt,nat
	integer start,stop,nr,is,imin,jmin,nmol,nres,natomstot
	real occ,bfac,um,dotpij,rij,mrad,bl
	real Rb,T,cc,gc
	character(16)crap,atype,restype,chain,aat,protname,cnmol
	character(128)line
	real, dimension(:), allocatable :: x,y,z,rx,ry,rz,lx,ly,lz,lmag,nx,ny,nz,hx,hy,hz,nhx,nhy,nhz,nhmag
	real, dimension(:,:), allocatable :: ldot
	integer, dimension(:), allocatable :: n,natoms
	double precision, dimension(:,:), allocatable :: qm,msf
	double precision, dimension(:), allocatable :: eig
	DATA P/3.141592654/
	!get number of lines, starting point
	open(unit=2,file="protname.txt",status='old')
	protname=adjustl(protname)
	read(2,*)protname
	close(2)
	open(unit=2,file="nmol.dat",status='old')
	read(2,*)nmol
	close(2)
	ALLOCATE ( natoms(nmol),n(nmol) )
	n=0
	nres=0
	natoms=0
	natomstot=0
	do i=1,nmol
	write(cnmol,*)i
	cnmol=adjustl(cnmol)
	open(unit=2,file="nres"//trim(cnmol)//".dat",status='old')
	read(2,*)n(i)
	close(2)
	nres=nres+n(i)
	open(unit=2,file="natoms"//trim(cnmol)//".dat",status='old')
	read(2,*)natoms(i)
	close(2)
	natomstot=natomstot+natoms(i)
	end do

	ALLOCATE ( x(natomstot),y(natomstot),z(natomstot),rx(nres),ry(nres),rz(nres),lx(nres),ly(nres),lz(nres),lmag(nres) )
	ALLOCATE ( eig(nres) )
	ALLOCATE ( nx(nres),ny(nres),nz(nres),hx(nres),hy(nres),hz(nres),nhx(nres),nhy(nres),nhz(nres),nhmag(nres) )
	ALLOCATE ( ldot(nres,nres) )
	ALLOCATE ( qm(nres,nres),msf(nres,nres) )

	ntot=0
	HYDtot=0
	NITtot=0
	mrad=0.0
	nat=0

! Edited by @pgromano to include explicit search or atomtypes, and not assume
! consistent order to input topology file.

	open(unit=23,file="mrad.dat")
	do k=1,nmol
		write(cnmol,*)k
		cnmol=adjustl(cnmol)
		open(unit=13,file=trim(protname)//"_"//trim(cnmol)//".pdb")
	!	now read in to grab values
		do i=1,natoms(k)
			nat=nat+1
			read(13,'(A)')line
			read(line(14:16),'(A)')atype
			read(line(18:20),'(A)')restype
			atype=adjustl(atype)
			restype=adjustl(restype)
			read(line(32:38),'(F8.3)')x(nat)
			read(line(40:46),'(F8.3)')y(nat)
			read(line(48:54),'(F8.3)')z(nat)
			IF(atype(1:2).eq."CA")THEN
				ntot=ntot+1
				rx(ntot)=x(nat)
				ry(ntot)=y(nat)
				rz(ntot)=z(nat)
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
			ELSE IF(atype(1:2).eq."N ")THEN
				NITtot=NITtot+1
				nx(NITtot)=x(nat)
				ny(NITtot)=y(nat)
				nz(NITtot)=z(nat)
			ELSE IF(restype(1:3).ne."PRO")THEN
				IF(atype(1:2).eq."H ")THEN
					HYDtot = HYDtot+1
					IF(HYDtot.lt.ntot)THEN
						HYDtot = ntot
					END IF
					hx(HYDtot)=x(nat)
					hy(HYDtot)=y(nat)
					hz(HYDtot)=z(nat)
				END IF
			ELSE IF(restype(1:3).eq."PRO")THEN
				IF(atype(1:2).eq."CD")THEN
					HYDtot = HYDtot+1
					IF(HYDtot.lt.ntot)THEN
						HYDtot = ntot
					END IF
					hx(HYDtot)=x(nat)
					hy(HYDtot)=y(nat)
					hz(HYDtot)=z(nat)
				END IF
			END IF
		end do
		close(13)
	end do !end loop over molecules
	close(23)

!	calculate bond lengths
	open(unit=103,file="length")
	open(unit=12345, file='NH_length')
	lmag=0.0
	j=1
	i=1 !mol number
	nt=n(1)
	do k=1,ntot-1
!	write(*,*)j
	lx(j)=rx(k+1)-rx(k)
	ly(j)=ry(k+1)-ry(k)
	lz(j)=rz(k+1)-rz(k)
	lmag(j)=(lx(j)**2+ly(j)**2+lz(j)**2)**.5
	nhx(j)=hx(k+1)-nx(k+1)
	nhy(j)=hy(k+1)-ny(k+1)
	nhz(j)=hz(k+1)-nz(k+1)
	nhmag(j)=(nhx(j)**2+nhy(j)**2+nhz(j)**2)**.5
!	if(lmag(j).gt.4.5)write(*,*)k,j,lmag(j),n(i)
	if(j.eq.nt)then !drop bonds between molecules
	j=j-1
	i=i+1
	nt=nt+n(i)
	end if
	j=j+1
	end do
	do i=1,ntot-nmol
	write(103,*)lmag(i)*.1
	write(12345,*)nhmag(i)*.1
	end do
	close(103)
	close(12345)


	ldot=0.0
	!calculate lNH dot lCA
	do i=1, ntot-1
	do j=1, ntot-1
	ldot(i,j)=(nhx(i)*lx(j)+nhy(i)*ly(j)+nhz(i)*lz(j))/(lmag(j)*nhmag(i))
	end do
	end do

	open(unit=36,file='ldot.dat')
	do i=1,ntot-nmol
	do j=1,ntot-nmol
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
	do i=1,ntot-nmol
	bl=bl+lmag(i)*.1
	end do
	close(18)
	bl=bl/real(ntot-nmol)
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
	write(16,*)ntot-nmol
	do i=1,ntot-nmol
	do j=1,ntot-nmol
	dotpij=lx(i)*lx(j)+ly(i)*ly(j)+lz(i)*lz(j)
	dotpij=(dotpij+cc*msf(i,j)+cc*msf(i+1,j+1)-cc*msf(i,j+1)-cc*msf(i+1,j))
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

	open(unit=100,file='cm_pre',status='unknown')
	do i=1,ntot
	do j=1,ntot
	write(100,*)cm(i,j)
	end do
	end do

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
