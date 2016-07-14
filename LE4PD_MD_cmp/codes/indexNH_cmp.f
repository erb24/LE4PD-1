	program input
	integer ntot,nmol
	ntot=0
	open(unit=2,file="protname.txt",status='old')
	read(2,*)
	read(2,*)ntot
	close(2)
	open(unit=2,file="nmol.dat",status='old')
	read(2,*)nmol
	close(2)
	call index(ntot,nmol)
	end program

	subroutine index(ntot,nmol) !gets atomindex for C2 and C4' from .gro file
	integer i,j,k,anum,rnum,ntot,ios,nlines,maxlines,nmol
	integer start,stop,nr,nres(nmol),nt,l
	real x,y,z
	character(16)crap,atype,restype,protname,cnmol
	character(128)line
	integer, dimension(ntot) :: bbN,bbH
	
	!get number of lines, starting point
	write(*,*)"ntot: ", ntot
	open(unit=2,file="protname.txt",status='old')
	read(2,*)protname
	close(2)
	do i=1,nmol
	write(cnmol,*)i
	cnmol=adjustl(cnmol)	
	open(unit=2,file="nres"//trim(cnmol)//".dat",status='old')
	read(2,*)nres(i)
	close(2)
	end do
	protname=adjustl(protname)
	open(unit=13,file=trim(protname)//".gro")
	nlines=0
	ios=0
	maxlines=150000
	start=3
	is=0
	stop=0
	nr=0
	anum=0
	rnum=0
	x=0.0
	y=0.0
	z=0.0
	do i=1,maxlines
	Read(13,*,IOSTAT=ios)line
	If(ios /= 0) EXIT
	If (J == maxlines) then
	write(*,*) "Maxlines exceeded"
	STOP
	End If
	nlines=nlines+1
	end do
	stop=nlines-1
	nr=stop-start
	write(*,*)start,stop,nr,nlines
	close(13)
	!read in C4' indexes for chain
	open(unit=13,file=trim(protname)//".gro")
	bbN=0
c	c4B=0
	bbH=0
c	c2B=0
	rnum=1
	l=0 !mol number
	nt=1
	do i=1,start-1
	read(13,*)
	end do
	do i=1,nr
	read(13,'(A)')line
	if(line(14:16).eq.' N ')then
	write(*,*)"HEREA"
	read(line(17:20),*)anum
	bbN(rnum)=anum
	write(*,'(A)')line(14:16),line(17:20),anum
c	if(line(14:16).eq.' N ')write(*,'(A)')"yes!",line(14:16)
	rnum=rnum+1
	end if
	end do
	open(unit=10,file="nitro.ndx")
	write(10,'(A)')"[ nitro ]"
	do i=1,ntot
	write(10,*)bbN(i)
	end do
	close(10)
	close(13)
	!read in C2 index for chain
	open(unit=13,file=trim(protname)//".gro")
	rnum=1
	do i=1,start-1
	read(13,*)
	end do
	do i=1,nr
	read(13,'(A)')line
	if(line(6:8).ne.'PRO')then
	if(line(14:16).eq.' H ')then
	write(*,*)"HEREA"
	read(line(17:20),*)anum
	bbH(rnum)=anum
	write(*,*)line(14:16),line(17:20),anum
	rnum=rnum+1
	end if
	end if
	if(line(6:8).eq.'PRO')then
	if(line(13:16).eq.' CD ')then
	write(*,*)"HEREA"
	read(line(17:20),*)anum
	bbH(rnum)=anum
	write(*,*)line(14:16),line(17:20),anum
	rnum=rnum+1
	end if
	end if
	if(rnum.eq.nt)then
	if(line(13:16).eq.' H1 ')then
	write(*,*)"HEREA"
	read(line(17:20),*)anum
	bbH(rnum)=anum
	write(*,*)line(14:16),line(17:20),anum
	rnum=rnum+1
	l=l+1
	nt=nt+nres(l)
	end if
	end if
	end do
	open(unit=10,file="hydro.ndx")
	write(10,'(A)')"[ hydro ]"
	do i=1,ntot
	write(10,*)bbH(i)
	end do
	close(10)
	close(13)
	end subroutine
