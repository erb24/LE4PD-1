	program pdbsplit !splits pdb into conformer structures
	integer i,j,k,anum,rnum,nmol,ios,nlines,natoms(10),ntot(10),maxlines,nc,l,m
	integer start,stop,nr,is,imin,jmin,stat
	character(16)protname,cnc,cnmol
	character(128)line
	nc=1
	write(cnc,*)nc
	cnc=adjustl(cnc)

	!get number of lines, starting point
	open(unit=2,file="protname.txt",status='old')
	read(2,*)protname
	close(2)
	protname=adjustl(protname)
	nlines=0
	ios=0
	maxlines=1500000
	start=0
	is=0
	stop=0
	nr=0
	ntot=0
	natoms=0
	nmol=1
	write(cnmol,*)nmol
	cnmol=adjustl(cnmol)
	open(unit=13,file=trim(protname)//".pdb")
	open(unit=14,file=trim(protname)//trim(cnc)//"_"//trim(cnmol)//".pdb")
	do i=1,maxlines
	write(*,*)"nc:",nc
	Read(13,'(A)',IOSTAT=ios)line
	If(ios /= 0) EXIT
	If (J == maxlines) then
	write(*,*) "Maxlines exceeded"
	STOP
	End If
	if(line(1:4) .eq. "ATOM" .AND. is .eq. 0)then
	start=nlines1
	is=1
	end if
	if(line(1:4) .eq. "ATOM")then
	write(14,'(A)')line
	natoms(nmol)=natoms(nmol)+1
	end if
	if(line(14:15).eq."CA")then
	ntot(nmol)=ntot(nmol)+1
	end if
	if(is .eq. 1 .AND. line(1:3) .eq. "TER") then
	close(14)
	open(unit=1,file="nres"//trim(cnmol)//".dat",position='append')
	write(1,*)ntot(nmol)
	close(1)
	open(unit=1,file="natoms"//trim(cnmol)//".dat",position='append')
	write(1,*)natoms(nmol)
	close(1)
	nmol=nmol+1
	write(cnmol,*)nmol
	cnmol=adjustl(cnmol)
	close(14)
	open(unit=14,file=trim(protname)//trim(cnc)//"_"//trim(cnmol)//".pdb")
	write(cnmol,*)nmol
	cnmol=adjustl(cnmol)
	end if	
	if(is .eq. 1 .AND. line(1:6) .eq. "ENDMDL") then
	close(14, status='delete')
	open(unit=1,file="nmol.dat",position='append')
	write(1,*)nmol-1
	close(1)
	ntot=0
	natoms=0
	nc=nc+1
	write(*,*)nc
	write(cnc,*)nc
	cnc=adjustl(cnc)
	nmol=1
	write(cnmol,*)nmol
	cnmol=adjustl(cnmol)
	close(14)
	open(unit=14,file=trim(protname)//trim(cnc)//"_"//trim(cnmol)//".pdb")
	end if
	write(cnmol,*)nmol
	cnmol=adjustl(cnmol)
	if(is .eq. 1 .AND. line(1:6) .eq. "MASTER")then
	close(14, status='delete')
	EXIT
	end if
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
	
