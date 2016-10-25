	module declare
	implicit none
	double precision :: b0,br,bz,bth,ga,r0,a,psi,q,dq,rb,pi
	double precision :: psi0,psimin,dpsi,r,z,rho,th,dth,dl
	double precision :: gradpsi,f,fp,rho0,rho1,diff,mu0,ip
	double precision :: psit,psi1,xi,de,area,darea
	integer :: i,j,m,n,mm,nn
	character(len=20) :: psich,qch,mch,psitch
	character(len=6) :: ich
	end module declare 

        program tearing_mode_profile
        use declare
        implicit none
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc                                                                     ccc
ccc This program computes the q-profile & flux surface area profile     ccc
ccc for a Solev'ev-type equilibrium                                     ccc
ccc                                                                     ccc
ccc Last modified 19/7/07.                                              ccc
ccc                                                                     ccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c Initiate output file
c
        open(8,file="qprofile.dat",status='unknown')
	psich="    Psi(Tm^2) "
	psitch="  area(m^2)"
	qch="     q(Psi)  "
	ich=" Ip = "
	mch="MA"
	write(8,fmt='(3a14)')psich,qch,psitch
c
c Input equilibrium parameters
c
	write(6,fmt='("Input B0, psi0, gamma, R0 & Rb")')
	read(5,*)b0,psi0,ga,r0,rb
	psi1=-psi0*ga/8.0*rb**4
c
c Input tearing mode parameters
c
	write(6,fmt='("Input m, n, xi & delta")')
	read(5,*)mm,nn,xi,de
c
c Input numerical parameters
c
	write(6,fmt='("Input no. of q-values & abscissae")')
	read(5,*)n,m
c
	pi=4.0*atan(1.0)
	mu0=4.0*pi*1.e-7
	psimin=-psi0*ga/8.0*rb**4
	dpsi=-psimin/(n-1.0)
c
c Compute q, area for magnetic axis
c
	psi=psimin/1.0001
	dth=pi/(m-1.0)
	th=0.0
	r=sqrt(r0**2+sqrt(rb**4+8.0*psi/ga/psi0))
	br=0.0
	bz=psi0*ga/2.0*(r**2-r0**2)
	bth=sqrt(br**2+bz**2)
	rho=r-r0	
	gradpsi=abs(bz)*r
	dl=rho*bth*dth/abs(bz*cos(th)-br*sin(th))
	q=dl/r/gradpsi/2.0
	do j=2,m-1
	   th=th+dth
c
c Use Newton-Raphson iteration to compute rho for given psi & th
c at magnetic axis 
c
	   rho0=rho
 20	   r=r0+rho0*cos(th)
	   z=rho0*sin(th)
	   f=ga/8.0*((r**2-r0**2)**2-rb**4)+(1.0-ga)/2.0*(r*z)**2
	   f=psi-psi0*f
	   fp=(1.0-ga)*z**2/rho0*r*(r+rho0*cos(th))
	   fp=-psi0*(ga/2.0*(r**2-r0**2)*r*cos(th)+fp)
	   rho1=rho0-f/fp
	   diff=abs((rho1-rho0)/rho0)
	   if(diff.lt.1.e-4)goto 10
	   rho0=rho1
	   goto 20
 10	   rho=rho1
c
	   br=-psi0*(1.0-ga)*r*z
	   bz=psi0*(ga/2.0*(r**2-r0**2)+(1.0-ga)*z**2)
	   bth=sqrt(br**2+bz**2)
     	   gradpsi=sqrt(br**2+bz**2)*r
	   dl=rho*bth*dth/abs(bz*cos(th)-br*sin(th))
	   dq=dl/r/gradpsi
	   q=q+dq
	   enddo
	th=pi
	r=sqrt(r0**2-sqrt(rb**4+8.0*psi/ga/psi0))
	br=0.0
	bz=psi0*ga/2.0*(r**2-r0**2)
	bth=sqrt(br**2+bz**2)
	rho=r0-r	
	gradpsi=abs(bz)*r
	dl=rho*bth*dth/abs(bz*cos(th)-br*sin(th))
	dq=dl/r/gradpsi/2.0 
	q=q+dq
	q=r0*b0/pi*q
	area=0.0
	write(8,fmt='(1p,3d14.5)')psimin,q,area
	psi=psimin
c
c Compute q, area for the rest of the plasma
c
	do i=2,n
	   psi=psi+dpsi
	   ip=0.0
	   th=0.0
	   r=sqrt(r0**2+sqrt(rb**4+8.0*psi/ga/psi0))
	   br=0.0
	   bz=psi0*ga/2.0*(r**2-r0**2)
	   bth=sqrt(br**2+bz**2)
	   rho=r-r0	
	   gradpsi=abs(bz)*r
	   dl=rho*bth*dth/abs(bz*cos(th)-br*sin(th))
	   q=dl/r/gradpsi/2.0
	   area=r*dl/2.0
	   ip=ip+bth*dl/2.0
	   do j=2,m-1
	      th=th+dth
c
c Use Newton-Raphson iteration to compute rho for given psi & th 
c
	      rho0=rho
 40	      r=r0+rho0*cos(th)
	      z=rho0*sin(th)
	      f=ga/8.0*((r**2-r0**2)**2-rb**4)+(1.0-ga)/2.0*(r*z)**2
	      f=psi-psi0*f
	      fp=(1.0-ga)*z**2/rho*r*(r+rho*cos(th))
	      fp=-psi0*(ga/2.0*(r**2-r0**2)*r*cos(th)+fp)
	      rho1=rho0-f/fp
	      diff=abs((rho1-rho0)/rho0)
	      if(diff.lt.1.e-4)goto 30
	      rho0=rho1
	      goto 40
 30	      rho=rho1
c
	      br=-psi0*(1.0-ga)*r*z
	      bz=psi0*(ga/2.0*(r**2-r0**2)+(1.0-ga)*z**2)
	      bth=sqrt(br**2+bz**2)
	      gradpsi=sqrt(br**2+bz**2)*r
	      dl=rho*bth*dth/abs(bz*cos(th)-br*sin(th))
	      dq=dl/r/gradpsi
	      darea=r*dl
	      q=q+dq
	      area=area+darea
	      ip=ip+bth*dl
	      enddo
	   th=pi
	   r=sqrt(r0**2-sqrt(rb**4+8.0*psi/ga/psi0))
	   br=0.0
	   bz=psi0*ga/2.0*(r**2-r0**2)
	   bth=sqrt(br**2+bz**2)
	   rho=r0-r	
	   gradpsi=abs(bz)*r
	   dl=rho*bth*dth/abs(bz*cos(th)-br*sin(th))
	   dq=dl/r/gradpsi/2.0 
	   darea=r*dl/2.0
	   q=q+dq
	   area=area+darea
	   q=r0*b0/pi*q
	   area=4.0*pi*area
	   ip=ip+bth*dl/2.0
	   write(8,fmt='(1p,3d14.5)')psi,q,area
	   enddo
	   psi=9.0
	   q=9.0
	   area=9.0
	   write(8,fmt='(1p,3d14.5)')psi,q,area
	   ip=ip/mu0/1.e6*2.0
	   write(8,fmt='()')
	   write(8,fmt='(a12,1p,d14.5,a6)')ich,ip,mch
        end program tearing_mode_profile































