	program pdbsplit !splits pdb into conformer structures
	integer i,j,k,anum,rnum,ntot,ios,nlines,maxlines,nc,l,m,natoms
	integer start,stop,nr,is,imin,jmin
	character(16)protname,cnc
	character(128)line
	nc=1
	write(cnc,*)nc
	cnc=adjustl(cnc)

	!get number of lines, starting point
	open(unit=2,file="protname.txt",status='old')
	read(2,*)protname
	close(2)
	protname=adjustl(protname)
	open(unit=13,file=trim(protname)//".pdb")
	open(unit=14,file=trim(protname)//trim(cnc)//".pdb")
	nlines=0
	ios=0
	maxlines=150000
	start=0
	is=0
	stop=0
	nr=0
	ntot=0
	natoms=0
	do i=1,maxlines
	Read(13,'(A)',IOSTAT=ios)line
	If(ios /= 0) EXIT
	If (J == maxlines) then
	write(*,*) "Maxlines exceeded"
	STOP
	End If
	if(line(1:4) .eq. "ATOM" .AND. is .eq. 0)then
	start=nlines
	is=1
	end if
	if(line(1:4) .eq. "ATOM")then
	write(14,'(A)')line
	natoms=natoms+1
	end if
	if(line(14:15).eq."CA")then
	ntot=ntot+1
	end if
	if(is .eq. 1 .AND. line(1:3) .eq. "TER") then
	open(unit=1,file="nres.dat",position='append')
	write(1,*)ntot
	close(1)
	open(unit=1,file="natoms.dat",position='append')
	write(1,*)natoms
	close(1)
	ntot=0
	close(14)
	nc=nc+1
	write(cnc,*)nc
	cnc=adjustl(cnc)
	open(unit=14,file=trim(protname)//trim(cnc)//".pdb")
	end if	
	if(is .eq. 1 .AND. line(1:6) .eq. "MASTER") EXIT
	nlines=nlines+1
	end do
	stop=nlines
	nr=stop-start

	write(*,*)start,stop,nr,nlines
	close(13)

	open(unit=1,file="ncopies.dat")
	write(1,*)nc-1
	close(1)

	end program
	
