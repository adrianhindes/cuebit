       program CUEBIT_MAST_ripple
        implicit none
        integer :: i,j,k,nt,ntt,n,idum,zzz,iseed
        double precision :: e,rm,pi,tpi,za,aa,r0,b0,th2
        double precision :: dt,dt1,det,phi2,psi0,tau,epsilon0
        double precision :: x,y,z,timp,vimp,bx,by,bz,ax,ay,az
        double precision :: r,phi,bphi,br,t,vx,vy,vz,vrotx,vroty
        double precision :: as,bs,vth,g(3),rco,rr,xcen,ycen
        double precision :: vxm,vym,vzm,tx,ty,tz,nop,rstart
        double precision,dimension(:),allocatable::x1,x2,y1,y2,z1,z2
        double precision,dimension(:),allocatable::vx1,vx2,vy1,vy2
        double precision,dimension(:),allocatable::vz1,vz2
        double precision,dimension(:),allocatable::phi1,r1,r2
        character(20) :: tch,rch,zch,pphich
        double precision :: vphi,pphi1,pphi2,errpphi
        double precision :: fi,gamma,err0,err1,derr,psi
        double precision :: rb,omega,er,ephi,ez,ex,ey,psi1
        double precision :: temppsi,npsi,temp0,temp1,n0,n1
        double precision :: coeffa,coeffb,coeffc,xroot1,xroot2
        double precision :: nnought,iii,r11,ph11,mu11,th11
        double precision :: t9,r9,z9,deltabph,deltabr,phi1d
        integer :: ntactual,i1,i2,radius,nco
        integer :: ione,itwo,nrem,nlost
        integer,dimension(:),allocatable::nflag
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc                                                                        ccc
ccc This version of CUEBIT solves the Lorentz force equation for           ccc
ccc collisionless charged particles in a MAST-like tokamak with Freidberg  ccc 
ccc poloidal flux surfaces. This version is used to model ripple transport ccc
ccc of beam ions. The beam ions are born with a uniform distribution       ccc
ccc within a sphere of radius 15cm & centre prescribed in input file.      ccc
ccc                                                                        ccc
ccc An implicit energy-conserving scheme is used. Iteration is used        ccc
ccc to make the overall scheme 2nd order accurate.                         ccc
ccc                                                                        ccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c Last modified on April 7 2011 by Ken McClements
c
c Initiate output files
c
        open(7,file="data_beam.dat",status='unknown')
c        open(15,file="orbit.dat",status='unknown')
c        open(16,file="full_orbit.dat",status='unknown')
c
        tch = "    t(s) "
        rch = "    R(m) "
        zch = "    Z(m) "
        pphich = " delta-Pphi"
c
c Physical & mathematical constants
c
        e = 1.602d-19
        rm = 1.672d-27
        epsilon0 = 8.85d-12
        pi = 4.d0*atan(1.d0)
        tpi = 2.d0*pi
c
c Specify parameters of Freidberg equilibrium
c
        r0 = 0.964
        rb = 0.93
        b0 = 0.4
        psi0 = 0.9
        gamma = 0.8
c
c Input particle mass number & charge state, time step (in units of 
c Larmor period at magnetic axis), number of
c time steps, and fraction of time steps output
c 
        read(7,fmt='(2f5.1,f6.3,i9,f7.4)')aa,za,dt1,nt,fi
        ntt = fi*(nt-2)
c
c Input random number seed & number of simulated particles
c
        read(7,fmt='(i5,i6)')iseed,n
        call srand(iseed)
        idum = 0
c
c Input coordinates of centre of beam deposition profile
c
        read(7,fmt='(2f6.2)')xcen,ycen
c
c Input number & major radius of coils 
c
        read(7,fmt='(i5,f6.2)')nco,rco
c
c Evaluate time step in seconds
c 
        dt = 2.d0*pi*aa*rm/za/e/b0*dt1
c
        allocate(x1(n),x2(n),y1(n),y2(n),z1(n),z2(n))
        allocate(vx1(n),vx2(n),vy1(n),vy2(n),vz1(n),vz2(n))
        allocate(phi1(n),r1(n),r2(n),nflag(n))
c
c Input initial position & velocity of particles   
c       
        nrem = 0
        nlost = 0
        do k = 1,n
           r11 = rand(idum)/3.0
           ph11 = tpi*rand(idum)
           mu11 = 2.0*(rand(idum)-0.5)
           th11 = acos(mu11)
	   x1(k)=xcen+0.15*(3.0*r11)**0.333*sin(th11)*cos(ph11)
	   y1(k)=ycen+0.15*(3.0*r11)**0.333*sin(th11)*sin(ph11)
	   z1(k)=0.15*(3.0*r11)**0.333*cos(ph11)
	   r1(k)=sqrt(x1(k)**2+y1(k)**2)
           vx1(k) = 0.0
           vy1(k) = -2.59e6
           vz1(k) = 0.0
           nrem = nrem+1
        enddo
c
c Open output file
c
      open(12,file="MAST_loss.dat",status='unknown')
	write(12,fmt='("Z =",f5.1)')za
	write(12,fmt='("A =",f5.1)')aa
	write(12,fmt='("delta-t =",1pd10.3)')dt
	write(12,fmt='("No. of time steps =",i10)')nt
	write(12,fmt='("No. of coils =",i4)')nco
	write(12,fmt='("Major radius of coils =",f6.2)')rco
	write(12,fmt='()')
c	write(12,fmt='("Initial particle positions")')
c	write(12,fmt='()')
c	write(12,fmt='("    R(m)     Z(m)")')
c
c Output initial particle positions
c
c        do k = 1,nrem
c           write(12,fmt='(1p2e10.3)')r1(k),z1(k)
c        enddo 
c
c Commence time loop
c
        t = 0.d0
        vphi = -vy1(1)*cos(phi1(1))
        psi = psi0*((gamma/8.0*((r1(1)**2-r0**2)**2-rb**4))
     &         +((1-gamma)/2.0*r1(1)**2*z1(1)**2))
        pphi1 = aa*rm*r1(1)*vphi+za*e*psi
        do i = 2,nt
        t = t+dt
c
c Commence particle loop (n > 1 only in multi-particle versions of code)
c
        do k = 1,nrem
           if(nflag(k).eq.1)cycle
c
c Evaluate ripple amplitude
c
        deltabph = (b0*r0/r1(k))*(r1(k)/rco)**nco*cos(nco*phi1(k))
	deltabr = (b0*r0/r1(k))*(r1(k)/rco)**nco*sin(nco*phi1(k)) 
c
c Evaluate magnetic field
c
        bphi = b0*r0/r1(k)
        bphi = bphi+deltabph
        br = -psi0*z1(k)*r1(k)*(1-gamma)
        br = br+deltabr
        bx = br*dcos(phi1(k))-bphi*dsin(phi1(k))
        by = br*dsin(phi1(k))+bphi*dcos(phi1(k))
        bz = psi0*(gamma/2.0*(r1(k)**2-r0**2)+(1-gamma)*z1(k)**2)
c
c Approximate new position & velocity of particle with old values 
c
        x2(k) = x1(k)
        y2(k) = y1(k)
        z2(k) = z1(k)
        vx2(k) = vx1(k)
        vy2(k) = vy1(k)
        vz2(k) = vz1(k)
c
c Commence iteration loop
c
        do j = 1,3
        x = (x1(k)+x2(k))/2.d0
        y = (y1(k)+y2(k))/2.d0
        z = (z1(k)+z2(k))/2.d0
        vx = (vx1(k)+vx2(k))/2.d0
        vy = (vy1(k)+vy2(k))/2.d0
        vz = (vz1(k)+vz2(k))/2.d0
        r = sqrt(x**2+y**2)
        phi = atan(y/x)
        if(x.lt.0.d0)phi = phi+pi
        if(x.gt.0.d0.and.y.lt.0.d0)phi = phi+tpi
        if(phi.gt.tpi)phi = phi-tpi
c
c Evaluate ripple amplitude at current time
c
        deltabph = (b0*r0/r)*(r/rco)**nco*cos(nco*phi)
	deltabr = (b0*r0/r)*(r/rco)**nco*sin(nco*phi) 
c
c Evaluate magnetic field components
c
        bphi = b0*r0/r
        bphi = bphi+deltabph
        br = -psi0*z*r*(1-gamma)
        br = br+deltabr
        bx = br*dcos(phi)-bphi*dsin(phi)
        by = br*dsin(phi)+bphi*dcos(phi)
        bz = psi0*(gamma/2.0*(r**2-r0**2)+(1-gamma)*z**2)
c
c Evaluate quantities appearing in matrix relating vi and vi+1
c
        ax = za*e*bx/2.d0/aa/rm*dt
        ay = za*e*by/2.d0/aa/rm*dt
        az = za*e*bz/2.d0/aa/rm*dt
        det = 1.d0+ax**2+ay**2+az**2
c
c Advance equation dvx/dt = Ze/m(vy*Bz-vZ*By)
c
        vx2(k) = (1.d0+ax**2-ay**2-az**2)*vx1(k)
        vx2(k) = vx2(k)+2.d0*(ay*ax+az)*vy1(k)
        vx2(k) = vx2(k)+2.d0*(az*ax-ay)*vz1(k)
        vx2(k) = vx2(k)/det
c
c Advance equation dvy/dt = Ze/m(vz*Bx-vx*By)
c
        vy2(k) = 2.d0*(ay*ax-az)*vx1(k)
        vy2(k) = vy2(k)+(1.d0-ax**2+ay**2-az**2)*vy1(k)
        vy2(k) = vy2(k)+2.d0*(ay*az+ax)*vz1(k)
        vy2(k) = vy2(k)/det
c
c Advance equation dvz/dt = Ze/m(vx*By-vy*Bx)
c
        vz2(k) = 2.d0*(az*ax+ay)*vx1(k)
        vz2(k) = vz2(k)+2.d0*(ay*az-ax)*vy1(k)
        vz2(k) = vz2(k)+(1.d0-ax**2-ay**2+az**2)*vz1(k)
        vz2(k) = vz2(k)/det
c
c Advance equation dx/dt = vx
c
        x2(k) = x1(k)+(vx2(k)+vx1(k))/2.d0*dt
c
c Advance equation dy/dt = vy
c
        y2(k) = y1(k)+(vy2(k)+vy1(k))/2.d0*dt
c
c Advance equation dz/dt = vz
c
        z2(k) = z1(k)+(vz2(k)+vz1(k))/2.d0*dt
c
        r2(k)=sqrt(x2(k)**2+y2(k)**2)
c
c Terminate iteration loop
c
        enddo
        x1(k) = x2(k)
        y1(k) = y2(k)
        z1(k) = z2(k)
        r1(k) = r2(k)
        phi1(k) = atan(y1(k)/x1(k))
        if(x1(k).lt.0.d0)phi1(k) = phi1(k)+pi
        if(x1(k).gt.0.d0.and.y1(k).lt.0.d0)phi1(k) = phi1(k)+tpi
        if(phi1(k).gt.tpi)phi1(k) = phi1(k)-tpi
        psi = psi0*((gamma/8.0*((r1(k)**2-r0**2)**2-rb**4))
     &  +((1-gamma)/2.0*r1(k)**2*z1(k)**2))
        vx1(k) = vx2(k)
        vy1(k) = vy2(k)
        vz1(k) = vz2(k)
        vphi = -vx1(k)*dsin(phi1(k))+vy1(k)*dcos(phi1(k))
        pphi2 = aa*rm*r1(k)*vphi+za*e*psi
        errpphi = (pphi2-pphi1)/pphi1
c
c Test to determine whether particle is lost
c
        if(psi > 0.0)nlost = nlost+1
        if(psi > 0.0)nflag(k) = 1
c
c End particle loop
c
        enddo
c
c End time loop
c
        enddo
	write(12,fmt='()')
	write(12,fmt='("Total no. of particles =",i7)')nrem
	write(12,fmt='("No. of lost particles =",i8)')nlost
        END
