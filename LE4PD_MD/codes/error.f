	program inputreader
	integer n
	character(32)protname
	open(unit=114,file='protname.txt',status='old')
	read (114,*)protname
	read(114,*)n
	close(114)
	call errors(n)
	end program inputreader

	subroutine errors(n)
	integer i,j,k,bn,nlines,maxlines,ios, resn,nexp
	real, dimension(n) :: expres,T1,T2,NOE,xT1,xt2,xNOE
	real xx,err,tt1,tt2,ttn,totalerr,fratio,vv
	real sigxT1,sigxT2,sigxNOE,sigT1,sigT2,sigNOE,avxT1,
     &avxT2,avxNOE,avT1,avT2,avNOE,corrT1,corrT2,corrNOE
	character(32)protname
	open(unit=114,file='protname.txt',status='old')
	read (114,*)protname
	close(114)
	protname=adjustl(protname)
c	write(*,*)n
	vv=0.0
	nexp=0
	
	open(unit=10,file="T1")
	open(unit=11,file="T2")
	open(unit=12,file="NOE")

	open(unit=13,file=trim(protname)//"_T1_exp.dat",status='old')
	open(unit=14,file=trim(protname)//"_T2_exp.dat",status='old')
	open(unit=15,file=trim(protname)//"_NOE_exp.dat",status='old')
	ttn=0.0
	tt2=0.0
	tt1=0.0
	resn=0
	T1=0.0
	T2=0.0
	NOE=0.0
	!read in data from files
	do i=1,n-1
	read(10,*)resn,tt1
	read(11,*)resn,tt2
	read(12,*)resn,ttn
	T1(resn)=tt1
	T2(resn)=tt2
	NOE(resn)=ttn
	end do
	nlines=0
	ios=0
	maxlines=150000
	xT1=0.0
	xT2=0.0
	xNOE=0.0
	xx=0.0

	do i=1,maxlines
	Read(13,*,IOSTAT=ios)
	If(ios /= 0) EXIT
	If (J == maxlines) then
	write(*,*) "Maxlines exceeded"
	STOP
	End If
	nlines=nlines+1
	end do
c	nlines=140
	nlines=nlines-1
c	write(*,*)"nlines",nlines
	close(13)
	open(unit=13,file=trim(protname)//"_T1_exp.dat",status='old')
	do i=1,nlines
c	write(*,*)i
	read(13,*)resn,xx
	xT1(resn)=xx
	end do

	nlines=0
	do i=1,maxlines
	Read(14,*,IOSTAT=ios)
	If(ios /= 0) EXIT
	If (J == maxlines) then
	write(*,*) "Maxlines exceeded"
	STOP
	End If
	nlines=nlines+1
	end do
c	nlines=140
	nlines=nlines-1
c	write(*,*)"nlines",nlines
	close(14)
	open(unit=14,file=trim(protname)//"_T2_exp.dat",status='old')
	do i=1,nlines
c	write(*,*)i
	read(14,*)resn,xx
	xT2(resn)=xx
	end do

	nlines=0
	do i=1,maxlines
	Read(15,*,IOSTAT=ios)
	If(ios /= 0) EXIT
	If (J == maxlines) then
	write(*,*) "Maxlines exceeded"
	STOP
	End If
	nlines=nlines+1
	end do
c	nlines=140
	nlines=nlines-1
c	write(*,*)"nlines",nlines
	close(15)
	open(unit=15,file=trim(protname)//"_NOE_exp.dat",status='old')
	do i=1,nlines
c	write(*,*)i
	read(15,*)resn,xx
	xNOE(resn)=xx
	end do

	!write scatter data sets
	open(unit=1,file='T1_scatter.dat')
	open(unit=2,file='T2_scatter.dat')
	open(unit=3,file='NOE_scatter.dat')
	do i=2,n
	if(xT1(i).ne.0.0)write(1,*)xT1(i),T1(i)
	if(xT2(i).ne.0.0)write(2,*)xT2(i),T2(i)
	if(xNOE(i).ne.0.0)write(3,*)xNOE(i),NOE(i)
	end do
	close(1)
	close(2)
	close(3)

	totalerr=0.0
	do i=2,n
c	write(*,*)"theory: ",T1(i),T2(i),NOE(i)
c	write(*,*)"exp: ",xT1(i),xT2(i),xNOE(i)
! if experimental data is missing, set equal to theoretical
	if(xT1(i).ne.0.0.AND.xT2(i).ne.0.0.AND.xNOE(i).ne.0.0)then
! calculate sum of relative errors of T1, T2, NOE
	nexp=nexp+1
	totalerr=totalerr+abs(T1(i)-xT1(i))/abs(xT1(i))
     &+abs(T2(i)-xT2(i))/abs(xT2(i))+abs(NOE(i)-xNOE(i))/abs(xNOE(i))
	end if
	end do
	totalerr=totalerr/real(3*nexp)

!obtain correlation coefficient
	j=0
	avxT1=0.0
	avxT2=0.0
	avxNOE=0.0
	avT1=0.0
	avT2=0.0
	avNOE=0.0
	do i=2,n !first get averages
	if(xT1(i).ne.0.0)then
	j=j+1
	avxT1=avxT1+xT1(i)
	avxT2=avxT2+xT2(i)
	avxNOE=avxNOE+xNOE(i)
	avT1=avT1+T1(i)
	avT2=avT2+T2(i)
	avNOE=avNOE+NOE(i)
	end if
	end do
	avxT1=avxT1/Real(j)
	avxT2=avxT2/Real(j)
	avxNOE=avxNOE/Real(j)
	avT1=avT1/Real(j)
	avT2=avT2/Real(j)
	avNOE=avNOE/Real(j)
	sigxT1=0.0
	sigxT2=0.0
	sigxNOE=0.0
	sigT1=0.0
	sigT2=0.0
	sigNOE=0.0
	corrT1=0.0
	corrT2=0.0
	corrNOE=0.0

	do i=2,n !now get correlations
	if(xT1(i).ne.0.0)then
	j=j+1
	sigxT1=sigxT1+(xT1(i)-avxT1)**2	
	sigxT2=sigxT2+(xT2(i)-avxT2)**2
	sigxNOE=sigxNOE+(xNOE(i)-avxNOE)**2
	sigT1=sigT1+(T1(i)-avT1)**2	
	sigT2=sigT2+(T2(i)-avT2)**2
	sigNOE=sigNOE+(NOE(i)-avNOE)**2
	corrT1=corrT1+(T1(i)-avT1)*(xT1(i)-avxT1)
	corrT2=corrT2+(T2(i)-avT2)*(xT2(i)-avxT2)
	corrNOE=corrNOE+(NOE(i)-avNOE)*(xNOE(i)-avxNOE)
	end if
	end do
	corrT1=corrT1/((sigxT1*sigT1)**.5)
	corrT2=corrT2/((sigxT2*sigT2)**.5)
	corrNOE=corrNOE/((sigxNOE*sigNOE)**.5)
	write(*,*)(corrT1+corrT2+corrNOE)*.33
	write(*,*)"cT1:",corrT1,"cT2:",corrT2,"cNOE:",corrNOE

	open(unit=20,file='error.dat')
	write(*,*)"mean relative error from exp:",totalerr
	write(20,*)totalerr
	close(20)
	open(unit=20,file='correlation.dat')
	write(20,*)'#total'
	write(20,*)(corrT1+corrT2+corrNOE)*.33
	write(20,*)'#T1,T2,NOE'
	write(20,*)corrT1,corrT2,corrNOE
	close(20)
	end subroutine
	
	
	
