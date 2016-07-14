!	This programs reads naccess output	
	implicit double precision(a-h,o-z)
	Dimension ax(1000), ay(1000), atot(1000)
	Dimension ras(1000), tot(1000)
	character(32)protname   
	open(unit=23,file='protname.txt',status='old')
	read(23,*)protname
	read(23,*)ntot
	protname=adjustl(protname)
	pi=3.14159
        open(unit=1,file=trim(protname)//".rsa",status='old')
        open(unit=2,file="resrad", status='unknown',position='append')
	open(unit=3,file="resvan", status='unknown',position='append')
	do i=1,4
	read(1,*)
	end do       
        do ia=1,ntot
	read(1,'(14x,2E8.2)')ax(ia),ay(ia)
	if (ay(ia) .GT. 0) then
	atot(ia)=ax(ia)/(ay(ia)/100)
	else
	atot(ia)=0
	end if
	end do
	do ia=1,ntot
	if (0 .lt. ax(ia)) then
	ras(ia)= sqrt(ax(ia)/(4*pi))
	tot(ia)=sqrt(atot(ia)/(4*pi))
	else
	ras(ia)=0
	tot(ia)=0
	endif 
	write(2,*)ras(ia)
	write(3,*)tot(ia)
	end do

	close(1)
	close(44)
	close(2)
	close(3)
	stop
	end
