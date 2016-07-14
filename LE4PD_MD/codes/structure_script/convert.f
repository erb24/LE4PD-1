	program clickfinder
	integer i,j,k,nbins,ix,iy,itheta,iphi,phibins
	real x,y,theta,phi,delha,delphi,xmin,xmax,ymin,ymax
	real xrange,yrange,delx,dely,x2theta,y2phi,pi,degr,otheta,ophi
	character(32)ctheta,cphi
	nbins=90
	delha=360./real(nbins)
	pi=3.1415927
	degr=(2.*pi)/360.
	open(unit=1,file="xmax")	
	read(1,*)xmax
	close(1)
	open(unit=1,file="ymax")	
	read(1,*)ymax
	close(1)
	open(unit=1,file="xmin")	
	read(1,*)xmin
	close(1)
	open(unit=1,file="ymin")	
	read(1,*)ymin
	close(1)
	xrange=abs(xmax-xmin)
	yrange=abs(ymax-ymin)
	x2theta=(180.-2.*delha)/xrange
	y2phi=(360.-2.*delha)/yrange
	open(unit=1,file="x")
	read(1,*)x
	close(1)
	open(unit=1,file="y")
	read(1,*)y
	close(1)
	x=abs(x-xmin)
	y=abs(y-ymin)
	theta=x*x2theta
	phi=y*y2phi
	if(theta.lt.delha)theta=delha
	write(*,*)theta,phi
	itheta=nint(theta/delha)
	otheta=itheta*delha
	phibins=nint(nbins*sin(otheta*degr))
	delphi=360./real(phibins)
	iphi=nint(phi/delphi)
	ophi=iphi*delphi
	write(ctheta,'(F4.0)')otheta
	ctheta=adjustl(ctheta)
	write(cphi,'(F4.0)')ophi
	cphi=adjustl(cphi)
	write(*,*)"theta: "//ctheta//"  phi: "//cphi
	open(unit=1,file="theta")
	write(1,'(A)')trim(ctheta)
	open(unit=2,file="phi")
	write(2,'(A)')trim(cphi)
	close(1)
	close(2)
	open(unit=1,file="feline")
	write(1,*)(itheta-1)*(nbins-1)+(itheta-1)+nint(phi/delha)
	close(1)
	end program
	
