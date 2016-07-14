	program input !calculates input from g96
	integer i,j,k,anum,rnum,ntot,ios,nlines,maxlines
	integer start,stop,n,is,nfrs
	real x,y,z,t
	character(16)crap,protname
	character(128)line
	!get number of lines, starting point
	open(unit=2,file="protname.txt")
	read(2,*)protname
	close(2)
	protname=adjustl(protname)
	open(unit=13,file=trim(protname)//".g96")

	nlines=0
	ios=0
	maxlines=150000000
	start=0
	is=0
	stop=0
	n=0
	nfrs=0
	do i=1,maxlines
	Read(13,*,IOSTAT=ios)line
	If(ios /= 0) EXIT
	If (J == maxlines) then
	write(*,*) "Maxlines exceeded"
	STOP
	End If
	if(line(1:4) .eq. "POSI")then
	start=nlines
	is=1
	nfrs=nfrs+1
	end if
	if(line(1:3) .eq. "END" .and. is.eq.1)then
	stop=nlines-1
	is=0
	n=stop-start
	end if
	nlines=nlines+1
	end do
	nfrs=nfrs-1
	write(*,*)n,nfrs
	close(13)
	open(unit=2,file="protname.txt",status='old',
     &position='append')
	write(2,*)n
	write(2,*)nfrs
	close(2)
	end program
