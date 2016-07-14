	program pdbcalc !calculates umatrix, hmatrix, input
	integer i,j,k,anum,rnum,ntot,ios,nlines,maxlines
	integer start,stop,nr,is
	real occ,bfac,pi,um,dotpij,rij,mrad
	realcomx,comy,comz,asph,proj
	real, dimension(15000) :: x,y,z
	character(16)crap,atype,restype,chain,aat,protname
	character(128)line
	real, dimension(500) :: rx,ry,rz,var,lx,ly,lz,lmag
	real, dimension(500) :: nx,ny,nz,hx,hy,hz,nhx,nhy,nhz,nhmag
	real, dimension(3,3) :: gtens,dgtens
	real, dimension(12) :: work
	real, dimension(3) :: eig,px,py,pz

	
	!get number of lines, starting point
	open(unit=2,file="protname.txt",status='old')
	read(2,*)protname
	close(2)
	protname=adjustl(protname)
	open(unit=13,file=trim(protname)//"_first.pdb")
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
	nlines=nlines+1
	if(is .eq. 1 .AND. line(1:3) .eq. "TER") EXIT
	end do
	stop=nlines
	nr=stop-start

	write(*,*)start,stop,nr,nlines
	close(13)

c	now read in to grab values
	open(unit=13,file=trim(protname)//"_first.pdb")
	open(unit=23,file="mrad.dat")
	ntot=0
	mrad=0.0
	x=0.0
	y=0.0
	z=0.0
	do i=1,start
	read(13,*)
	end do
	do i=1,nr-1
!	read(13,*)crap,anum,atype,restype,chain,rnum,x,y,z,
!     &occ,bfac,aat
	read(13,*)crap,anum,atype,restype,rnum,x(i),y(i),z(i),
     &occ,bfac
	if(atype(1:2).eq."CA")then
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
    	var(ntot)=bfac
c	write(*,*)bfac
c	assign Van Der Waals radii from miller values
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
c	calculate bond vectors
	pi=3.14159
	do j=1,ntot-1
	lx(j)=rx(j+1)-rx(j)
	ly(j)=ry(j+1)-ry(j)
	lz(j)=rz(j+1)-rz(j)
	lmag(j)=(lx(j)**2+ly(j)**2+lz(j)**2)**.5
	nhx(j)=hx(j+1)-nx(j+1)
	nhy(j)=hy(j+1)-ny(j+1)
	nhz(j)=hz(j+1)-nz(j+1)
	nhmag(j)=(nhx(j)**2+nhy(j)**2+nhz(j)**2)**.5
	write(*,*)j,"NH",nhmag(j)
	end do
	!calculate COM, gyration tensor, asphericity
	comx=0.0
	comy=0.0
	comz=0.0
	gtens=0.0
	eig=0.0
	asph=0.0
	px=0.0
	py=0.0
	pz=0.0
	proj=0.0
	do i=1,ntot
	comx=comx+rx(i)
	comy=comy+ry(i)
	comz=comz+rz(i)
	end do
	comx=comx/real(ntot)
	comy=comy/real(ntot)
	comz=comz/real(ntot)
	write(*,*)"COM",comx,comy,comz
	!assign new bead positions relative to COM
	do i=1,ntot
	rx(i)=rx(i)-comx
	ry(i)=ry(i)-comy
	rz(i)=rz(i)-comz
c	write(*,*)"coord",rx(i),ry(i),rz(i)
	end do
	!calculate gyration tensor
	do i=1,ntot
	gtens(1,1)=gtens(1,1)+rx(i)*rx(i)
	gtens(1,2)=gtens(1,2)+rx(i)*ry(i)
	gtens(1,3)=gtens(1,3)+rx(i)*rz(i)
	gtens(2,1)=gtens(2,1)+ry(i)*rx(i)
	gtens(2,2)=gtens(2,2)+ry(i)*ry(i)
	gtens(2,3)=gtens(2,3)+ry(i)*rz(i)
	gtens(3,1)=gtens(3,1)+rz(i)*rx(i)
	gtens(3,2)=gtens(3,2)+rz(i)*ry(i)
	gtens(3,3)=gtens(3,3)+rz(i)*rz(i)
	end do
	gtens=gtens/(real(ntot))
	do i=1,3
	write(*,*)gtens(i,1),gtens(i,2),gtens(i,3)
	end do
	!diagonalize gyration tensor
	call SSYEV("V","U",3,gtens,3,eig,work,12,io)
	write(*,*)"eig:",eig(1),eig(2),eig(3)
	do i=1,3
	write(*,*)gtens(i,1),gtens(i,2),gtens(i,3)
	end do
	!aspherocity
	open(unit=41,file="asphericity.dat")
	asph=eig(3)-.5*(eig(1)+eig(2))
	write(*,*)"asphericity:",asph
	write(41,*)"asphericity:",asph
	close(41)
	!projection of bond vectors
	px(1)=gtens(1,3)
	px(2)=gtens(1,2)
	px(3)=gtens(1,1)
	py(1)=gtens(2,3)
	py(2)=gtens(2,2)
	py(3)=gtens(2,1)
	pz(1)=gtens(3,3)
	pz(2)=gtens(3,2)
	pz(3)=gtens(3,1)
	open(unit=31,file="ca_proj1.dat")
	open(unit=32,file="ca_proj2.dat")
	open(unit=33,file="ca_proj3.dat")
	do i=1,ntot-1
	proj=px(1)*lx(i)+py(1)*ly(i)+pz(1)*lz(i)
	proj=proj/lmag(i)
	write(31,*)i,abs(proj)
	proj=px(2)*lx(i)+py(2)*ly(i)+pz(2)*lz(i)
	proj=proj/lmag(i)
	write(32,*)i,abs(proj)
	proj=px(3)*lx(i)+py(3)*ly(i)+pz(3)*lz(i)
	proj=proj/lmag(i)
	write(33,*)i,abs(proj)
	end do
	close(31)
	close(32)
	close(33)
	open(unit=31,file="NH_proj1.dat")
	open(unit=32,file="NH_proj2.dat")
	open(unit=33,file="NH_proj3.dat")
	do i=1,ntot-1
	proj=px(1)*nhx(i)+py(1)*nhy(i)+pz(1)*nhz(i)
	proj=proj/nhmag(i)
	write(31,*)i,abs(proj)
	proj=px(2)*nhx(i)+py(2)*nhy(i)+pz(2)*nhz(i)
	proj=proj/nhmag(i)
	write(32,*)i,abs(proj)
	proj=px(3)*nhx(i)+py(3)*nhy(i)+pz(3)*nhz(i)
	proj=proj/nhmag(i)
	write(33,*)i,abs(proj)
	end do
	close(31)
	close(32)
	close(33)
	open(unit=16,file='protname.txt',status='old',
     &position='append')
	write(16,*)nr-1
	close(16)

	end program
	

	


