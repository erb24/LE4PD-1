c******************************************************************************************
c p2m1.f : program to calculate P2(t) from simulations and from theory
C ******************************************************************************************
	program inputreader
	integer n
	open(unit=114,file='protname.txt',status='old')
	read (114,*)
	read(114,*)n
	close(114)
	call p2m1(n)
	end program inputreader
c*********************************************************
c  catena lineare per una proteina. E' il programma linearp1p2.f
C     CALCULATION OF FILE P2(T),TAU        
	subroutine p2m1(n)
	implicit double precision(a-h,o-z)
	integer ai1,ai2,ai3,ai4,fintime,k,l,i,j
	integer,dimension(n,4) :: nlist
	real factor,fricorr(n),eiglam1,eiglam2,eiglam3
	DATA P/3.141592654/
	dimension cc(n,n),tc(5000),tt(5000),avm1(5000)
        dimension dt(5000),rm(n,3),rnorm(5000),rmodenorm(5000)
	dimension fint0(5000,n),fints0(5000,n),fj0(5000)
	dimension epi0(5000,n),epis0(5000,n),pj0(5000),nbg(n)
	dimension ned(n),mexp(n),difm1(n),difp2(n)
        dimension al0(n),avmode(n),itcr(10)
	dimension qinvm(n,n),qm(n,n),eigmu(n),eiglam(n),dot(n,n)
	character(32)protname,ii

c
!        Open(unit=63,file='tau',status='unknown')
!	open(unit=67,file='taum1',status='unknown')
!	open(unit=65,file='p2',status='unknown')
c	open(unit=66,file='m1',status='unknown')
	open(unit=747,file='sigma',status='old')
	open(unit=789,file='length',status='old')
!	open(unit=790, file='rotamp.dat')
!	open(unit=791, file='intamp.dat')
	open(unit=114,file='protname.txt',status='old')
	read(114,*)protname
	close(114)
	protname=adjustl(protname)

c ***************************************************
c ************now the calculation from the theory****
c****************************************************
	fricorr=0.0
	Rb=.00198 !(boltzmanns constant in kcal/mol*K)
	T=300.0
	open(unit=10,file='fmad_mp_60.dat',status='old')
	do i=1,n-1
	read(10,*)fricorr(i)
	fricorr(i)=exp(fricorr(i)/(Rb*T))
	end do
	qinvm=0.0
	qm=0.0
	eigmu=0.0
	eiglam=0.0
	dot=0.0
	nm1=n-1
	cc=0.0
	rnorm=0.0
	rmodenorm=0.0
	eiglam1=0.0
	eiglam2=0.0
	eiglam3=0.0
	open(unit=99,file="lambda_eig",status='old') !read in all values
	do i=1,nm1
	read(99,*)eiglam(i)
c	write(*,*)eiglam(i)
	end do
	close(99)
	open(unit=99,file="QINVmatrix",status='old')
	do i=1,nm1
	do j=1,nm1
	read(99,*)qinvm(i,j)
c	write(*,*)qinvm(i,j)
	end do
	end do
	close(99)
	open(unit=99,file="Qmatrix",status='old')
	do i=1,nm1
	do j=1,nm1
	read(99,*)qm(i,j)
c	write(*,*)qm(i,j)
	end do
	end do
	close(99)
	open(unit=99,file="mu_eig",status='old')
	do i=1,nm1
	read(99,*)eigmu(i)
c	write(*,*)eigmu(i)
	end do
	close(99)
!	open(unit=13,file=trim(protname)//'_ldot.av',status='old')
!	do i=1,nm1
!	do j=1,nm1
!	read(13,*)dot(i,j)
c	write(*,*)dot(i,j)
!	end do
!	end do
!	close(13)
!	do i=1,nm1
!	write(ii,*)i
!	ii=adjustl(ii)
!	open(unit=110,file='qm_'//trim(ii))
!	do j=1,nm1
!	write(110,*)(qm(j,i)**2)/eigmu(i)
!	end do
!	close(110)
!	end do

!now calculate the needed cc(i,j) (mode amplitudes in m1 relaxation)
!first 3 modes, use NH orientation rexpressed from CA basis
!	do i=1,nm1 !residue loop
!	do j=1,3 !mode loop
!	do l=1,nm1 !dummy loop over CA vector basis set
!	cc(i,j)=cc(i,j)+qinvm(j,l)*dot(i,l)
!	end do !close summation
!	cc(i,j)=(cc(i,j)**2)*eigmu(j)
c	if(i.eq.10)write(*,*)cc(i,j)
!	end do !end mode loop
!	end do !end residue loop
!All modes, use regular
	do i=1,nm1 !residue loop
	do j=1,nm1 !mode loop
	cc(i,j)=(qm(i,j)**2)/eigmu(j)
c	if(i.eq.10)write(*,*)cc(i,j)
	end do !end mode loop
	end do !end residue loop
!	do i=1,nm1
!	do j=1,3
!	rnorm(i)=rnorm(i)+(qm(i,j)**2)/eigmu(j)
!	rmodenorm(i)=rmodenorm(i)+cc(i,j)
!	end do
!	do j=1,3 !need to normalize the sum of the first three modes to what it was in CA-CA representation
!	cc(i,j)=cc(i,j)*(rnorm(i)/rmodenorm(i))
!	end do !end mode loop
c	write(*,*)"after",cc(i,1)+cc(i,2)+cc(i,3),"before",
c     &rmodenorm(i)
!	end do
	tua=0.d0
	taum=0.d0
	k=0
	avmode=0.0
	rm=0.0
	rnorm=0.0
	rmodenorm=0.0
! put in rigid body eig combinations
	eiglam1=.5*(eiglam(1)+eiglam(2))
	eiglam2=.5*(eiglam(1)+eiglam(3))
	eiglam3=.5*(eiglam(2)+eiglam(3))
	eiglam(1)=eiglam1
	eiglam(2)=eiglam2
	eiglam(3)=eiglam3
!	open(unit=30,file='roteig.dat')
!	do j=1,3
!	write(30,*)eiglam(j)
!	end do
!	close(30)

        do ia=1,nm1
	som=0.d0
	do i=1,nm1
	som=som+cc(ia,i)
	end do
c	difsum=1-som
c	if(difsum.gt.0.5) then
c	write(*,*)"eigenvalue/vector are wrong"
c	goto 111
c	end if
c	write(64,*)ia,som
c	write(*,*)"unnormalized",som
c	renormalize
	do i=1,nm1
	cc(ia,i)=cc(ia,i)/som
	end do
	som2=0.0
	do i=1,nm1
	som2=som2+cc(ia,i)
	end do
c	write(*,*)"normalized",som2
c	all added to renormalize
	end do
        close(70)

c	open(unit=31,file='mode1_amp.dat')
c	open(unit=32,file='mode2_amp.dat')
c	open(unit=33,file='mode3_amp.dat')
	do i=1,nm1 !write first mode amplitude
c	write(31,*)i,cc(i,1)
c	write(32,*)i,cc(i,2)
c	write(33,*)i,cc(i,3)
	end do
c	close(31)
c	close(32)
c	close(33)
	do i=1,nm1
	do j=1,nm1
	avmode(j)=avmode(j)+cc(i,j)
	end do
c	write(*,*)"norm",rnorm(i)
	end do
	do j=1,nm1
	avmode(j)=avmode(j)/real(nm1)
	end do
	do j=1,nm1
c	if(j.le.3)write(790,*)avmode(j)
c	if(j.gt.3)write(791,*)avmode(j)
	end do
c	close(790)
c	close(791)
	nit=nm1

	apg=2.d0/p
	read(747,*) 
	read(747,*) sg ! sigma
	close(747)
	test=0.255d0

c this loop is to have a good representation of the statistics (many points at the beginning, less at the end)
	index=1
	amint=0.001
	tc(1)=0.d0
	amint=amint*1.d0
	tc(2)=amint
	do jt=1,4 
	jt1=jt-1
	jt2=jt
	jt3=jt+1
	jt4=jt+2
	ai1=10**jt1
	ai2=10**jt2
	ai3=10**jt3
	ai4=10**jt4

	do itm=ai1,ai4,ai1
	index=index+1
	aincrement=ai1*amint
	Tc(index+1)=tc(index)+aincrement
	end do
	end do
	fintime=index

	avm1=0.0
      SOM=0. 
      ISOM=0 
      NME=Nm1
      DO 27  Ires=1,NME ! big residue loop
	write(ii,*)ires
	ii=adjustl(ii)
!	open(unit=66,file='p2CArot_'//trim(ii),status='unknown')
	open(unit=60,file='m1CArot_'//trim(ii),status='unknown')
	fints0(1,ires)=0.d0
	fint0(1,ires)=0.d0
	fj0(ires)=0.d0
	epis0(1,ires)=0.d0
	epi0(1,ires)=0.d0
	pj0(ires)=0.d0
      ISOM=ISOM+1  
      SAM=0.   
      DO 28  L=1,nit 
c      WRITE(64,*)L,Ires,CC(Ires,L)  
      AF1=CC(ires,L)/eiglam(L) 
   28 SAM =SAM+AF1  
      blcheck=sam-alb

      DO 29  I=1,index !time loop
      SUM=0. 
      IU=0    
      DO 30  L=1,nit !mode loop
      ES=eiglam(L)*Tc(I) 
c      IF (ES.GT.fintime)GO TO 32 
      ESP1=DEXP(-ES)  
	if(L.ge.4)ESP1=DEXP(-ES/(fricorr(L)))
!	if(L.le.3)ESP1=1.0
      GO TO 31  
   32 ESP1=0.   
      IU=IU+1   
   31 B=CC(Ires,L)*ESP1  
      SUM=SUM+B 
   30 CONTINUE  !end mode
!?
c      cte=(alb/alen(ires))**2
c      sum=sum*cte
!?
cccccccccccccccccccccccccccccccccccccccccccccccccccxxxx
cccccccccccccccccccccccccccccccccccccccccccccccccccxxxx
c	write(60,*)tt(I),sum
      AM2=SUM**2  
      AM3=SUM**3 
	p2=0.d0
      IF(SUM.GT.TEST) then
! added
	If(SUM.GT..99999)then
	X2=0.0
!
	else
	X2=(1.-AM2)/AM2
	end if
c      write (*,*) AM2
      X1=DSQRT(X2)
      X3=X1**3.
      AT=DATAN(X1)
      P1=(1.-X2)*(1.-(APG*AT))+(APG*X1)
      P2=1.-(3.*X2)+(3.*X3*(1./APG-AT))
! ??
        if(p2.le.10E-9) then
	p2=1E-9
	endif
! ??
	else
      P1=APG*((4./3.*SUM)+(2./15.*(SUM**3.))+(3./70.*(SUM**5))+ 
     *(5./252.*(SUM**7))) 
      P2=(3./5.*AM2)+(6./35.*(AM2**2))+(8./105.*(AM2**3))+  
     *(16./385.*(AM2**4)) 
	endif
c	if (i.eq.1257) then
c	difp2(ires)=p2t(int(t500),ires)-(p2*cte)
c	write(668,*)ires,difp2(ires)*difp2(ires)
c	end if
	
	fint0(i+1,ires)=p2
	tt(i)=tc(i)/sg
	tintv=tc(i+1)-tc(i)
	addt=0.5d0*(fint0(i+1,ires)+fint0(i,ires))
	addt=addt*tintv
	fints0(i+1,ires)=fints0(i,ires)+addt 
	epi0(i+1,ires)=sum
	addm=0.5d0*(epi0(i+1,ires)+epi0(i,ires))
	addm=addm*tintv
	epis0(i+1,ires)=epis0(i,ires)+addm 
	if(p2.le.1E-8)p2=1E-8
!	write(65,*)tt(I),P2  ! note that we need a*here if we want matlab 
!	write(66,*)tt(I),P2
	write(60,*)tt(I),sum
	avm1(I)=avm1(I)+sum
c to readit
	! since this is here and not before writing p2, we have 
c to make sure that we multiply tau by sigma
	tua=fints0(fintime,ires)/sg
	taum=epis0(fintime,ires)/sg
   29 CONTINUE ! end of i (time) loop 
c	if(p2.ge.1E-7)write(*,*)"p2 of",ires,"hasnt decayed,=",p2
c	fj0(ires)=ptch

!	write(63,*)ires,tua
!	write(67,*)ires,taum
	close(60)
!	close(66)
   27 continue ! end if (residue) loop
   55 continue

	open(unit=61,file='avm1rot',status='unknown')
	do l=1,index
	tt(l)=tc(l)/sg
	avm1(l)=avm1(l)/real(n-1)
	write(61,*)tt(l),avm1(l)
	end do
	close(61)
c
!	open(unit=231,file='eigenvalues',status='unknown')
!	do ia=1,nm1
!	write(231,*)eiglam(ia)
!	end do
c	close(64)
c	close(65)
!	close(63)
c      RETURN
111    continue !for errors
	return

	end
