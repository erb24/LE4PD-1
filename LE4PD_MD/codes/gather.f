	implicit double precision(a-h,o-z)
	CHARACTER AMP*1
	integer ntot, crap

c	open(unit=1,file='newdata',status='old')
c	open(unit=2,file='T1them',status='unknown')
c	open(unit=3,file='T2them',status='unknown')
c	open(unit=4,file='NOEthem',status='unknown')
	open(unit=5,file='t1t2',status='old')
	open(unit=6,file='T1',status='unknown')
	open(unit=7,file='T2',status='unknown')
	open(unit=8,file='NOE',status='unknown')
	open(unit=9,file='protname.txt',status='old')
	read(9,*)
	read(9,*)ntot

c	AMP='&'
c	
c	do i=1,ntot-1
c	read(1,*)inx,b,c,d,e,f,g,h,o,p
c	IF(c.NE.99) THEN
c		write(2,*)inx,c
c	ELSE 
c		write(2,*)inx,AMP
c	END IF
c	IF(f.NE.99) THEN
c		write(3,*)inx,f
c	ELSE
c		write(3,*)inx,AMP
c	END IF
c	IF(o.NE.99) THEN
c		write(4,*)inx,o
c	ELSE
c		write(4,*)inx,AMP
c	END IF
c	end do

	do j=1,ntot-1
	read(5,*)indx,t1,t2,oen
	write(6,*)indx+1,t1
	write(7,*)indx+1,t2
	write(8,*)indx+1,oen
	end do
	
	close(1)
	close(2)
	close(3)
	close(4)
	close(5)
	close(6)
	close(7)
	close(8)
	END
