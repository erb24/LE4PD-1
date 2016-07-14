	program input
	integer i,j,k,anum,rnum,ntot,ios,nlines,maxlines,start,nr,stop
	integer nc
	character(128)posx,posy,posz,line
	character(32)protname,cnc
	!get number of lines, starting point
	open(unit=2,file="protname.txt",status='old')
	read(2,*)protname
	close(2)
	protname=adjustl(protname)
	nc=1
	do k=1,nc
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
	nlines=nlines+1
	if(is .eq. 1 .AND. line(1:3) .eq. "TER") EXIT
	end do
	stop=nlines
	nr=stop-start
	write(*,*)start,stop,nr,nlines
	close(13)
	open(unit=2,file="protname.txt",status='old')
	read(2,*)
	read(2,*)ntot
	close(2)

	call pdbmaker(start,ntot,nr,cnc)
	end do
	end program

	subroutine pdbmaker(start,ntot,nr,cnc)
	integer i,j,k,anum,rnum,ntot,ios,nlines,maxlines,nr,start
	integer c2(ntot),c4(ntot),nres,frame,mode
	character(128)posx,posy,posz,line,filename
	character(32)protname,canum,cframe,cmode,cnc,cmamp
	real amp(ntot),c4y(ntot),c4z(ntot),c2x(ntot),c2y(ntot),c2z(ntot)
	amp=0.0
	do mode=4,94
	write(cmode,*)mode
	cmode=adjustl(cmode)
	open(unit=1,file="qm_"//trim(cmode),status='old')
	do i=1,ntot-1
	read(1,*)amp(i)
!	amp(i)=(-1.0*log(amp(i)))**(-.6)
	amp(i)=(amp(i)**.5)*.385*10.0
	end do
	close(1)
	open(unit=2,file="protname.txt",status='old')
	read(2,*)protname
	close(2)
	protname=adjustl(protname)
	open(unit=13,file=trim(protname)//trim(cnc)//".pdb")
	open(unit=14,file=trim(protname)//trim(cnc)//'_m'//trim(cmode)//
     &".pdb")
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
	if(line(13:15) .eq. " N ")then
	nres=nres+1
	end if
	write(cmamp,'(F4.2)')amp(nres)
	line(63:67)=cmamp(1:4)
	write(14,'(A)')line
	end do

	close(13)
	close(14)
	end do !close mode loop
	end subroutine
