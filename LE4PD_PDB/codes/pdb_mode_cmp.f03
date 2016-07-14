	program inputreader
	integer n,nmol
	open(unit=114,file='protname.txt',status='old')
	read (114,*)
	read(114,*)n
	close(114)
	open(unit=114,file='nmol.dat',status='old')
	read(114,*)nmol
	close(114)
	call pdb(n,nmol)
	end program inputreader

	subroutine pdb(n1,nmol1) !splits pdb into conformer structures
	integer i,j,k,anum,rnum,nmol,ios,nlines,natoms(nmol1),ntot(nmol1),maxlines,nc,l,m,a
	integer start,stop,nr,is,imin,jmin,stat,n,n1,nmol1,ncmp
	character(16)protname,cnc,cnmol,cmode,cfml,atype
	character(128)line
	real fml,afml(n1)

	open(unit=2,file="prot",status='old')
	read(2,*)protname
	close(2)
	protname=adjustl(protname)
	
	afml=0.0
	do a=4,10
	write(cmode,*)a
	cmode=adjustl(cmode)
	open(unit=1,file="mlength_"//trim(cmode),status='old')
	do i=1,n1-nmol1
	read(1,*)l,fml
	afml(l)=fml
	end do
	close(1)

	nc=1
	write(cnc,*)nc
	cnc=adjustl(cnc)

	open(unit=15,file=trim(protname)//"_m"//trim(cmode)//".pdb")

	!get number of lines, starting point
	nlines=0
	ios=0
	maxlines=150000
	start=0
	is=0
	stop=0
	nr=0
	ntot=0
	natoms=0
	nmol=1
	ncmp=0
	write(cnmol,*)nmol
	cnmol=adjustl(cnmol)
	open(unit=13,file=trim(protname)//".pdb")
!	open(unit=14,file=trim(protname)//trim(cnc)//"_"//trim(cnmol)//".pdb")
	do i=1,maxlines
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
!	write(14,'(A)')line
	read(line(14:16),'(A)')atype
	atype=adjustl(atype)
	if(atype(1:3).eq."N  ")then
	ncmp=ncmp+1
	write(*,*)atype,ncmp
	end if
	write(cfml,'(F5.2)')afml(ncmp)
	line(62:67)=cfml(1:5)
	write(15,'(A)')line
	natoms(nmol)=natoms(nmol)+1
	end if
	if(line(14:15).eq."CA")then
	ntot(nmol)=ntot(nmol)+1
	end if
	if(is .eq. 1 .AND. line(1:3) .eq. "TER") then
	write(15,'(A)')line	
	close(14)
!	open(unit=1,file="nres"//trim(cnmol)//".dat",position='append')
!	write(1,*)ntot(nmol)
!	close(1)
!	open(unit=1,file="natoms"//trim(cnmol)//".dat",position='append')
!	write(1,*)natoms(nmol)
!	close(1)
	nmol=nmol+1
	write(cnmol,*)nmol
	cnmol=adjustl(cnmol)
!	close(14)
!	open(unit=14,file=trim(protname)//trim(cnc)//"_"//trim(cnmol)//".pdb")
	write(cnmol,*)nmol
	cnmol=adjustl(cnmol)
	end if	
	if(is .eq. 1 .AND. line(1:6) .eq. "ENDMDL") then
!	close(14, status='delete')
	write(15,'(A)')line
	open(unit=1,file="nmol.dat",position='append')
	write(1,*)nmol-1
	close(1)
	ntot=0
	ncmp=0
	natoms=0
	nc=nc+1
	write(cnc,*)nc
	cnc=adjustl(cnc)
	nmol=1
	write(cnmol,*)nmol
	cnmol=adjustl(cnmol)
	close(14)
!	open(unit=14,file=trim(protname)//trim(cnc)//"_"//trim(cnmol)//".pdb")
	end if
	write(cnmol,*)nmol
	cnmol=adjustl(cnmol)
	if(is .eq. 1 .AND. line(1:6) .eq. "MASTER")then
!	close(14, status='delete')
	write(15,'(A)')line
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

	end do !close loop over modes
	close(15)

	end subroutine
	
