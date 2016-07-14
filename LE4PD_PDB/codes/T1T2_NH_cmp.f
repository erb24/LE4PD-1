c	Program to calculate T1, T2 and NOE from P2(t)
c t1t2a.f has been modified to be used with automater
	program inputreader
	integer n,nmol
	open(unit=114,file='protname.txt',status='old')
	read (114,*)
	read(114,*)n
	close(114)
	open(unit=114,file='nmol.dat',status='old')
	read(114,*)nmol
	close(114)
	call t1t2(n,nmol)
	end program inputreader

	subroutine t1t2(n,nmol)
	implicit double precision(a-h,o-z)
	integer n,nmol
c
	dimension tt(0:5000,n),p2(0:5000,n)
	dimension fint0(0:5000,n)
	dimension ttb(0:5000,n),fints0(0:5000,n)
	dimension dt(0:5000,n),fj0(n)
	dimension fint1(0:5000,n),fints1(0:5000,n)
	dimension fint2(0:5000,n)
	dimension fints2(0:5000,n),fint3(0:5000,n)
	dimension fints3(0:5000,n)
	dimension fint4(0:5000,n),fints4(0:5000,n)
	dimension fj1(n),fj2(n),fj3(n),fj4(n)
	dimension t1a(n),t1b(n),t1c(n),t1(n)
	dimension t2a(n),t2b(n),t2c(n),t2(n),xnoe(n)
c
	open(unit=1,file="p2",status='old')
	open(unit=2,file="t1t2",status='unknown')
c	open(unit=3,file="check.txt",status='unknown')
! all the parameters below (save for pi) should be changeable. How could we do that?
	data uo/1.256637D-06/
	data hp/6.62608D-34/
	data wh/599.98D+06/ !this will be default, can we change this?
	data wn/60.8D+06/   !same
	data gh/26.7519D+07/
	data gn/-2.7126D+07/
	data dn/160D-06/
	rnhin3=9.90688385D+29
c	data rnhin3/9.34056542D+29/ !corrected Bax 2008
	data pi/3.141592654d0/
	factor=0.0
	open(unit=9,file='NHfactor.dat')
	read(9,*)factor
	close(9)
	rnhin3=rnhin3*factor
	write(*,*)rnhin3
	nstep=4001
	do ib=1,n-nmol
c	read(1,*)
	tt(0,ib)=0.d0
	p2(0,ib)=1.d0
	do ia=1,nstep
c	write(*,*) ia
	read(1,*)ttb(ia,ib),p2(ia,ib)
	tt(ia,ib)=ttb(ia,ib)*1D-12
	end do
	end do
c	
	w0=0.d0
	w1=wn
	w2=wh
	w3=wh+wn
	w4=wh-wn
c
	do ib=1,n-nmol
	fint0(0,ib)=1.d0
	fints0(0,ib)=0
	fj0(ib)=0
	do ia=1,nstep
	fint0(ia,ib)=p2(ia,ib)*dcos(w0*tt(ia,ib)*2*pi)
	dt(ia,ib)=tt(ia,ib)-tt(ia-1,ib)
	fints0(ia,ib)=fints0(ia-1,ib)+
     &  0.5*dt(ia,ib)*(fint0(ia,ib)+fint0(ia-1,ib))
	end do
	fj0(ib)=0.4*fints0(nstep,ib)
c	write(3,*)ib,fints0(nstep,ib)
	end do
c
	do ib=1,n-nmol
	fint1(0,ib)=1.d0
	fints1(0,ib)=0
	fj1(ib)=0
	do ia=1,nstep
	fint1(ia,ib)=p2(ia,ib)*dcos(w1*tt(ia,ib)*2*pi)
	dt(ia,ib)=tt(ia,ib)-tt(ia-1,ib)
	fints1(ia,ib)=fints1(ia-1,ib)+
     &  0.5*dt(ia,ib)*(fint1(ia,ib)+fint1(ia-1,ib))
	end do
	fj1(ib)=0.4*fints1(nstep,ib)
	end do
c	
	do ib=1,n-nmol
	fint2(0,ib)=1.d0
	fints2(0,ib)=0
	fj2(ib)=0
	do ia=1,nstep
	fint2(ia,ib)=p2(ia,ib)*dcos(w2*tt(ia,ib)*2*pi)
	dt(ia,ib)=tt(ia,ib)-tt(ia-1,ib)
	fints2(ia,ib)=fints2(ia-1,ib)+
     &  0.5*dt(ia,ib)*(fint2(ia,ib)+fint2(ia-1,ib))
	end do
	fj2(ib)=0.4*fints2(nstep,ib)
	end do
c
	do ib=1,n-nmol
	fint3(0,ib)=1.d0
	fints3(0,ib)=0
	fj3(ib)=0
	do ia=1,nstep
	fint3(ia,ib)=p2(ia,ib)*dcos(w3*tt(ia,ib)*2*pi)
	dt(ia,ib)=tt(ia,ib)-tt(ia-1,ib)
	fints3(ia,ib)=fints3(ia-1,ib)+
     &  0.5*dt(ia,ib)*(fint3(ia,ib)+fint3(ia-1,ib))
	end do
	fj3(ib)=0.4*fints3(nstep,ib)
	end do
c	
	do ib=1,n-nmol
	fint4(0,ib)=1.d0
	fints4(0,ib)=0
	fj4(ib)=0
	do ia=1,nstep
	fint4(ia,ib)=p2(ia,ib)*dcos(w4*tt(ia,ib)*2*pi)
c	write(3,*)tt(ia,ib),w4*tt(ia,ib)*2*pi,fint4(ia,ib)
	dt(ia,ib)=tt(ia,ib)-tt(ia-1,ib)
	fints4(ia,ib)=fints4(ia-1,ib)+
     &  0.5*dt(ia,ib)*(fint4(ia,ib)+fint4(ia-1,ib))
	end do
	fj4(ib)=0.4*fints4(nstep,ib)
	end do
c	
	di=uo*uo*hp*hp*gh*gh*gn*gn
	d2=di*rnhin3*rnhin3/(64*pi*pi*pi*pi)
	c2=4*pi*pi*wn*wn*dn*dn/3
c
	do ib=1,n-nmol
	t1a(ib)=fj4(ib)+3*fj1(ib)+6*fj3(ib)
	t1b(ib)=0.25*d2*t1a(ib)
	t1c(ib)=t1b(ib)+c2*fj1(ib)
	t1(ib)=1/t1c(ib)
	t2a(ib)=4*fj0(ib)+fj4(ib)+3*fj1(ib)+6*fj2(ib)+6*fj3(ib)
	t2b(ib)=d2*t2a(ib)/8
	t2c(ib)=t2b(ib)+(c2*(3*fj1(ib)+4*fj0(ib))/6)
	t2(ib)=1/t2c(ib)
	xnoe(ib)=1+0.25*d2*(gh/gn)*(6*fj3(ib)-fj4(ib))*t1(ib)
c	write(3,*)ib,xnoe(ib),6*fj3(ib),4*fj4(ib)
	end do
c
	do ib=1,n-nmol
	write(2,400)ib,t1(ib),t2(ib),xnoe(ib)
400	format(1x,i3,2x,f15.8,2x,f15.8,2x,f15.8)
	end do
c	
	close(1)
	close(2)
	end
