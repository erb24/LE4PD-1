	program input
	integer i,j,k,anum,rnum,ntot,ios,nlines,maxlines,start,nr,stop
	integer natoms
	character(128)posx,posy,posz,line
	character(32)protname
	!get number of lines, starting point
	open(unit=2,file="../../protname.txt",status='old')
	read(2,*)protname
	read(2,*)ntot
	close(2)
	open(unit=2,file="../../natoms.dat",status='old')
	read(2,*)natoms
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

	call pdbmaker(start,ntot,nr,natoms)
	end program

	subroutine pdbmaker(start,ntot,nr,natoms)
	integer i,j,k,anum,rnum,ntot,ios,nlines,maxlines,nr,start
	integer c2(natoms),c4(natoms),nres,frame,mode,natoms
	character(128)posx,posy,posz,line,filename
	character(32)protname,canum,cframe,cmode
	real c4x(natoms),c4y(natoms),c4z(natoms),c2x(ntot),c2y(ntot)
     &,c2z(ntot)
	open(unit=1,file="framefile",status='old')
	read(1,'(A)')filename
	filename=adjustl(filename)
	close(1)
	c4x=0.0
	c4y=0.0
	c4z=0.0
	open(unit=1,file=trim(filename),status='old')
	do i=1,natoms
	read(1,*)c4x(i),c4y(i),c4z(i)
	c4x(i)=c4x(i)*10.0
	c4y(i)=c4y(i)*10.0
	c4z(i)=c4z(i)*10.0
	end do
	close(1)
	
	open(unit=1,file="frame",status='old')
	read(1,*)frame
	close(1)
	open(unit=1,file="mode",status='old')
	read(1,*)mode
	close(1)
c	open(unit=1,file="CA.ndx",status='old')
c	read(1,*)
c	do i=1,ntot
c	read(1,*)c4(i)
c	end do
	close(1)
	open(unit=2,file="protname.txt",status='old')
	read(2,*)protname
	close(2)
	protname=adjustl(protname)
	open(unit=13,file=trim(protname)//"_first.pdb")
	write(cframe,*)frame
	cframe=adjustl(cframe)
	write(cmode,*)mode
	cmode=adjustl(cmode)
	open(unit=14,file="struct2_m"//trim(cmode)//
     &'_'//trim(cframe)//".pdb")
	nres=0
	do i=1,start
	Read(13,*)
	end do
	do i=1,nr-1
	Read(13,'(A)')line
	line(63:66)="1.00"
c	canum=line(8:12)
c	read(canum,'(I4)')anum
c	write(*,*)anum
!	if(line(13:15) .eq. " P ")then
	nres=nres+1
	read(line(8:12),'(I4)')c4(nres)
!	write(posx,'(F6.2)')c4x(nres)
!	write(posy,'(F6.2)')c4y(nres)
!	write(posz,'(F6.2)')c4z(nres)
	write(posx,*)c4x(nres)
	write(posy,*)c4y(nres)
	write(posz,*)c4z(nres)
	posx=adjustl(posx)
	posy=adjustl(posy)
	posz=adjustl(posz)
	line(32:37)=posx(1:6)
	line(40:45)=posy(1:6)
	line(48:53)=posz(1:6)
	write(14,'(A)')line
	write(*,'(A)')line
!	end if
	end do
!	do i=1,ntot-1
!	line="                                                "
!	line(1:7)="CONECT "
!	write(line(8:11),'(I4)')c4(i)
!	write(line(13:16),'(I4)')c4(i+1)
!	write(*,*)line
!	write(14,'(A)')line
!	end do
	close(13)
	close(14)

	end subroutine
