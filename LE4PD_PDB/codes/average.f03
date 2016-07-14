	program average !calculates averages over conformers
	integer i,j,k,ntot,nc
	real avbl,bl
	character(16)crap,protname,cnc
	character(128)line
	real, dimension(:), allocatable :: resrad,avresrad
	double precision, dimension(:,:), allocatable :: um,avum,rij,avrij,ldot,avldot
	open(unit=2,file="protname.txt",status='old')
	read(2,*)protname
	read(2,*)ntot
	close(2)
	open(unit=2,file="ncopies.dat",status='old')
	read(2,*)nc
	close(2)
	ALLOCATE ( resrad(ntot),avresrad(ntot) )
	ALLOCATE ( um(ntot,ntot),avum(ntot,ntot),rij(ntot,ntot),avrij(ntot,ntot),ldot(ntot,ntot),avldot(ntot,ntot) )
	avum=0.0
	avresrad=0.0
	avrij=0.0
	avldot=0.0
	avbl=0.0

	protname=adjustl(protname)
	write(cnc,*)nc
	cnc=adjustl(cnc)
	do i=1,nc
	um=0.0
	resrad=0.0
	rij=0.0
	ldot=0.0
	bl=0.0
	write(cnc,*)i
	cnc=adjustl(cnc)
	open(unit=10,file='avbl'//trim(cnc),status='old')
	open(unit=11,file='resrad'//trim(cnc),status='old')
	open(unit=12,file='Umatrix'//trim(cnc),status='old')
	open(unit=13,file='Rij'//trim(cnc),status='old')
	open(unit=14,file='ldot'//trim(cnc)//'.dat',status='old')
	
	read(10,*)bl
	avbl=avbl+bl

	read(12,*)
	do j=1,ntot-1
	do k=1,ntot-1
	read(12,*)um(j,k)
	read(14,*)ldot(j,k)
	avum(j,k)=avum(j,k)+um(j,k)
	avldot(j,k)=avldot(j,k)+ldot(j,k)
	end do
	end do
	
	do j=1,ntot
	do k=1,ntot
	read(13,*)rij(j,k)
	avrij(j,k)=avrij(j,k)+rij(j,k)
	end do
	end do
	
	do j=1,ntot
	read(11,*)resrad(j)
	avresrad(j)=avresrad(j)+resrad(j)
	end do
	
	close(10)
	close(11)
	close(12)
	close(13)
	close(14)
	end do !end copy loop

	avbl=avbl/real(nc)

	do j=1,ntot-1
	do k=1,ntot-1
	avum(j,k)=avum(j,k)/real(nc)
	avldot(j,k)=avldot(j,k)/real(nc)
	end do
	end do
	
	do j=1,ntot
	do k=1,ntot
	avrij(j,k)=avrij(j,k)/real(nc)
	end do
	end do
	
	do j=1,ntot
	avresrad(j)=avresrad(j)/real(nc)
	end do

	open(unit=10,file='avbl')
	open(unit=11,file='avresrad')
	open(unit=12,file='Umatrix')
	open(unit=13,file='Rij')
	open(unit=14,file='ldot.dat')
	write(10,*)avbl

	write(12,*)
	do j=1,ntot-1
	do k=1,ntot-1
	write(12,*)avum(j,k)
	write(14,*)avldot(j,k)
	end do
	end do
	
	do j=1,ntot
	do k=1,ntot
	write(13,*)avrij(j,k)
	end do
	end do
	
	do j=1,ntot
	write(11,*)avresrad(j)
	end do
	
	close(10)
	close(11)
	close(12)
	close(13)
	close(14)
	end program
