c******************************************************************************************
c p2m1.f : program to calculate P2(t) from simulations and from theory
C ******************************************************************************************
	program inputreader
	integer n,nmol
	open(unit=114,file='protname.txt',status='old')
	read (114,*)
	read(114,*)n
	close(114)
	open(unit=2,file="nmol.dat",status='old')
	read(2,*)nmol
	close(2)
	call p2m1(n,nmol)
	end program inputreader
c*********************************************************
c  catena lineare per una proteina. E' il programma linearp1p2.f
C     CALCULATION OF FILE P2(T),TAU        
	subroutine p2m1(n,nmol)
	implicit double precision(a-h,o-z)
	integer ai1,ai2,ai3,ai4,fintime,k,l,i,j,mcut,nmol
	integer,dimension(n,4) :: nlist
	real factor,fricorr(n),eiglam1,eiglam2,eiglam3,sumamp(n)
	real rot,mint,ti,tip1,g,avblsq,avblsqNH,test,cutint(n)
	real intamp(n),NHorder(n),testi,testr
	DATA P/3.141592654/
	dimension cc(n,n),tc(5000),tt(5000),avximsq(n),qmNH(n,n)
        dimension dt(5000),rm(n,3),rnorm(5000),rmodenorm(5000)
	dimension fint0(5000,n),fints0(5000,n),fj0(5000)
	dimension epi0(5000,n),epis0(5000,n),pj0(5000),nbg(n)
	dimension ned(n),mexp(n),difm1(n),difp2(n),blsqNH(n)
        dimension al0(n),avmode(n),itcr(10),ampNH(n),qmrot(n,n)
	dimension qinvm(n,n),qm(n,n),eigmu(n),eiglam(n),dot(n,n)
	character(32)protname,ii

c
        Open(unit=63,file='tau',status='unknown')
	open(unit=67,file='taum1',status='unknown')
	open(unit=65,file='p2',status='unknown')
	open(unit=66,file='m1',status='unknown')
	open(unit=747,file='sigma',status='old')
	open(unit=789,file='length',status='old')
	open(unit=790, file='rotamp.dat')
	open(unit=791, file='intamp.dat')
	open(unit=114,file='protname.txt',status='old')
	read(114,*)protname
	close(114)
	protname=adjustl(protname)

c ***************************************************
c ************now the calculation from the theory****
c****************************************************
	ampNH=0.0
	Rb=.00198 !(boltzmanns constant in kcal/mol*K)
	T=300.0
	open(unit=10,file='fmad_mp_60.dat',status='old')
	do i=1,n-nmol
	read(10,*)fricorr(i)
	fricorr(i)=exp(fricorr(i)/(Rb*T))
	end do
	do i=1,n-nmol
	AmpNH(i)=.02 !from Bax JACS 1998
!	write(*,*)"ampNH:",ampNH(i)
	end do
	qinvm=0.0
	qm=0.0
	qmNH=0.0
	qmrot=0.0
	eigmu=0.0
	eiglam=0.0
	dot=0.0
	nm1=n-nmol
	cc=0.0
	rnorm=0.0
	rmodenorm=0.0
	eiglam1=0.0
	eiglam2=0.0
	eiglam3=0.0
	cut=.1
	small=.00000001
	rot=0.0
	mint=0.0
	ti=0.0
	tip1=0.0
	blsqNH=0.0
	avximsq=0.0
	mcut=75
	cutint=0.0
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
	open(unit=99,file="Qmatrix_NH",status='old')
	do i=1,nm1
	do j=1,nm1
	read(99,*)qmNH(i,j)
c	write(*,*)qm(i,j)
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
	open(unit=110,file='avblsqNH')
	read(110,*)avblsqNH
	close(110)
	open(unit=110,file='avblsq')
	read(110,*)avblsq
	close(110)
	open(unit=110,file='blsqNH')
	do j=1,nm1
	read(110,*)blsqNH(j)
	end do
	close(110)
	open(unit=110,file='ximsq')
	do j=1,nm1
	read(110,*)avximsq(j)
	end do
	close(110)
	NHorder=0.0




!now calculate the needed cc(i,j) (mode amplitudes in m1 relaxation)
	sumamp=0.0
	intamp=0.0
	do i=1,nm1 !residue loop
	do j=1,3 !mode loop
	cc(i,j)=((qmNH(i,j)**2)*avximsq(j))/(blsqNH(i))
	end do !end mode loop
	end do !end residue loop

	do i=1,nm1 !residue loop
	do j=4,nm1 !mode loop
	cc(i,j)=((qm(i,j)**2)/(eigmu(j)))
	intamp(i)=intamp(i)+cc(i,j)
	end do !end mode loop
	end do !end residue loop

	rmodenorm=0.0
	do i=1,nm1
	testi=0.0
	testr=0.0
	do j=1,3 !have to add in NH librational
	rnorm(i)=rnorm(i)+cc(i,j)
	rmodenorm(i)=rmodenorm(i)+(qm(i,j)**2)/eigmu(j)
	end do
	rmodenorm(i)=rmodenorm(i)-ampNH(i)
	if(rnorm(i).lt.0.0)rnorm(i)=0.0
	if(rnorm(i).lt.0.0)write(*,*)"rnorm:",i,"less than zero"
c	write(*,*)"S1NH:",rnorm(i)
	do j=1,3
	cc(i,j)=cc(i,j)*(rmodenorm(i)/rnorm(i))
	testr=testr+cc(i,j)
	end do !end mode loop
!	write(*,*)"before:",rnorm(i),"after",
!     &cc(i,1)+cc(i,2)+cc(i,3)
	end do

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

        do ia=1,nm1
	som=0.d0
	rot=0.
	mint=0. 
	do i=1,nm1
	som=som+cc(ia,i)
	end do
	som=som+ampNH(ia)
c	difsum=1-som
c	if(difsum.gt.0.5) then
c	write(*,*)"eigenvalue/vector are wrong"
c	goto 111
c	end if
c	write(64,*)ia,som
!	write(*,*)"unnormalized",som
c	renormalize
	do i=1,nm1
	cc(ia,i)=cc(ia,i)/som
	end do
	som2=0.0
	do i=1,nm1
	som2=som2+cc(ia,i)
	end do
	som2=som2+ampNH(ia)
!	write(*,*)"normalized",som2
c	all added to renormalize
	end do
        close(70)

c	open(unit=31,file='mode1_amp.dat')
c	open(unit=32,file='mode2_amp.dat')
c	open(unit=33,file='mode3_amp.dat')
	open(unit=31,file='S2gNHCAint.dat')
	do i=1,nm1 !write first mode amplitude
	g=cc(i,1)+cc(i,2)+cc(i,3)
	write(31,*)i+1,1.5*(g**2)-.5
c	write(31,*)i,cc(i,1)
c	write(32,*)i,cc(i,2)
c	write(33,*)i,cc(i,3)
	end do
	close(31)
c	close(31)
c	close(32)
c	close(33)

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

	
      SOM=0. 
      ISOM=0 
      NME=Nm1
      DO 27  Ires=1,NME ! big residue loop
	write(ii,*)ires
	ii=adjustl(ii)
c	open(unit=66,file='m1qmNHint_'//trim(ii),status='unknown')
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
c	if(L.le.3)ESP1=1.0
      GO TO 31  
   32 ESP1=0.   
      IU=IU+1   
   31 B=CC(Ires,L)*ESP1  
      SUM=SUM+B 
      if(L.le.3)rot=rot+B/rmodenorm(ires)
      if(L.gt.3)mint=mint+B/(1.-rmodenorm(ires))
   30 CONTINUE  !end mode
!?
c      cte=(alb/alen(ires))**2
c      sum=sum*cte
!?
cccccccccccccccccccccccccccccccccccccccccccccccccccxxxx
cccccccccccccccccccccccccccccccccccccccccccccccccccxxxx

	SUM=SUM+AmpNH(Ires)*DEXP(-((Tc(I))/(sg*.2)))
	if(I.eq.1)SUM=1.0
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
	write(65,*)tt(I),P2  ! note that we need a*here if we want matlab 
c	write(66,*)tt(I),SUM
c	write(66,*)tt(I),P2
c to readit
	! since this is here and not before writing p2, we have 
c to make sure that we multiply tau by sigma
	tua=fints0(fintime,ires)/sg
	taum=epis0(fintime,ires)/sg
   29 CONTINUE ! end of i (time) loop 
c	if(p2.ge.1E-7)write(*,*)"p2 of",ires,"hasnt decayed,=",p2
c	fj0(ires)=ptch

	write(63,*)ires,tua
	write(67,*)ires,taum
	close(60)
	close(66)
	close(61)
   27 continue ! end if (residue) loop
   55 continue
c
	open(unit=231,file='eigenvalues',status='unknown')
	do ia=1,nm1
	write(231,*)eiglam(ia)
	end do
c	close(64)
c	close(65)
	close(63)
c      RETURN
111    continue !for errors
	return

	end
