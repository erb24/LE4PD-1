c	Program to calculate T1, T2 and NOE from P2(t)
c t1t2a.f has been modified to be used with automater
	implicit double precision(a-h,o-z)
c
	dimension tt(0:5000,700),p2(0:5000,700)
	dimension fint0(0:5000,700)
	dimension ttb(0:5000,700),fints0(0:5000,700)
	dimension dt(0:5000,700),fj0(700)
	dimension fint1(0:5000,700),fints1(0:5000,700)
	dimension fint2(0:5000,700)
	dimension fints2(0:5000,700),fint3(0:5000,700)
	dimension fints3(0:5000,700)
	dimension fint4(0:5000,700),fints4(0:5000,700)
	dimension fj1(700),fj2(700),fj3(700),fj4(700)
	dimension t1a(700),t1b(700),t1c(700),t1(700)
	dimension t2a(700),t2b(700),t2c(700),t2(700),xnoe(700)
c
	open(unit=1,file="p2",status='old')
	open(unit=2,file="t1t2",status='unknown')
c	open(unit=3,file="check.txt",status='unknown')
	open(unit=8,file='protname.txt',status='unknown')
! all the parameters below (save for pi) should be changeable. How could we do that?
	data uo/1.256637D-06/
	data hp/6.62608D-34/
	data wh/599.98D+06/ !this will be default, can we change this?
	data wn/60.8D+06/   !same
	data gh/26.7519D+07/
	data gn/-2.7126D+07/
	data dn/170D-06/
!	data dn/160D-06/
	rnhin3=9.90688385D+29
c	data rnhin3/9.34056542D+29/ !corrected Bax 2008
	data pi/3.141592654d0/
	read(8,*)
	read(8,*)ntot
	if(ntot.gt.700)write(*,*)'arrays too small for t1t2.f'
	factor=0.0
	open(unit=9,file='NHfactor.dat')
	read(9,*)factor
	close(9)
	rnhin3=rnhin3*factor
	write(*,*)rnhin3
	nstep=4001
	do ib=1,ntot-1
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
	do ib=1,ntot-1
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
	do ib=1,ntot-1
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
	do ib=1,ntot-1
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
	do ib=1,ntot-1
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
	do ib=1,ntot-1
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
	do ib=1,ntot-1
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
	do ib=1,ntot-1
	write(2,400)ib,t1(ib),t2(ib),xnoe(ib)
400	format(1x,i3,2x,f15.8,2x,f15.8,2x,f15.8)
	end do
c	
	close(1)
	close(2)
	end
