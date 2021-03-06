        program CUEBIT_linear
        implicit none
        integer :: i,j,k,nt,ntt,n,idum,zzz
        double precision :: e,rm,pi,tpi,za,aa,r0,b0,th2
        double precision :: dt,dt1,det,phi2,psi0,tau,epsilon0
        double precision :: x,y,z,timp,vimp,bx,by,bz,ax,ay,az
        double precision :: r,phi,bphi,br,t,vx,vy,vz,vrotx,vroty
        double precision :: as,bs,vth,g(3)
        double precision :: vxm,vym,vzm,tx,ty,tz,nop,rstart
        double precision,dimension(:),allocatable::x1,x2,y1,y2,z1,z2
        double precision,dimension(:),allocatable::vx1,vx2,vy1,vy2
        double precision,dimension(:),allocatable::vz1,vz2
        double precision,dimension(:),allocatable::phi1,r1,r2
        character(20) :: tch,xch,ych,rch,zch,pphich
        double precision :: vphi,pphi1,pphi2,errpphi
        double precision :: fi,gamma,err0,err1,derr,psi
        double precision :: rb,omega,er,ephi,ez,ex,ey,psi1
        double precision :: temppsi,npsi,temp0,temp1,n0,n1
        double precision :: coeffa,coeffb,coeffc,xroot1,xroot2
        double precision :: nnought,iii
        double precision :: t9,r9,z9
        integer :: ntactual,i1,i2,radius
        integer :: ione,itwo
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc                                                                        ccc
ccc This version of CUEBIT solves the Lorentz force equation for a         ccc
ccc collisionless charged particle in a MAST-like tokamak with Freidberg   ccc
ccc poloidal flux surfaces. It is intended that this version will          ccc
ccc eventually be used to model ripple transport of beam ions.             ccc
ccc                                                                        ccc
ccc An implicit energy-conserving scheme is used. Iteration is used        ccc
ccc to make the overall scheme 2nd order accurate.                         ccc
ccc                                                                        ccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c Last modified on March 1 2011 by Ken McClements
c
c Initiate output files
c
        open(15,file="orbit.dat",status='unknown')
        open(16,file="full_orbit.dat",status='unknown')
c
        tch = "    t(s) "
        xch = "    X(m) "
        ych = "    Y(m) "
        rch = "    R(m) "
        zch = "    Z(m) "
        pphich = " delta-Pphi"
        write(16,fmt='(6a11)')tch,xch,ych,rch,zch,pphich
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
c Input initial position & velocity of particle
c
        n = 1
c
        allocate(x1(n),x2(n),y1(n),y2(n),z1(n),z2(n))
        allocate(vx1(n),vx2(n),vy1(n),vy2(n),vz1(n),vz2(n))
        allocate(phi1(n),r1(n),r2(n))
c
c Set all positions at origin, and speeds to zero except vx at sound speed
c
        r1(1) = 1.1
        z1(1) = 1.00
        phi1(1) = 0.00
        vx1(1) = 1.8d5
        vy1(1) = 1.8d5
        vz1(1) = 1.8d5
        phi1(1) = phi1(1)*pi/180.0
        x1(1) = r1(1)*dcos(phi1(1))
        y1(1) = r1(1)*dsin(phi1(1))
        vphi = -vx1(1)*dsin(phi1(1))+vy1(1)*dcos(phi1(1))
        psi = psi0*((gamma/8.0*((r1(1)**2-r0**2)**2-rb**4))
     &  +((1-gamma)/2.0*r1(1)**2*z1(1)**2))
c
c Input particle mass number & charge state, time step (in units of
c Larmor period at magnetic axis), number of
c time steps, and fraction of time steps output
c
        aa = 2.0
        za = 1.0
        dt1 = 0.01
        nt = 1000000
        fi = 1.000
        ntt = nt-2
        pphi1 = aa*rm*r1(1)*vphi+za*e*psi
c
c Evaluate time step in seconds
c
        dt = 2.d0*pi*aa*rm/za/e/b0*dt1
c
c Commence time loop
c
        t = 0.d0
        do i = 2,nt
        t = t+dt
c
c Commence particle loop (n > 1 only in multi-particle versions of code)
c
        do k = 1,n
c
c Evaluate magnetic field

           bz = b0 + 0.01*z1(1)*b0/r1(1)**2
           br = -0.01*bz/r1(1)**2/2
           bx = br*dcos(phi1(1))
           by = br*dsin(phi1(1))
           bphi = b0*0.01*phi1(1)/r1(1)
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
        enddo
c
c Evaluate magnetic field components
c
        bz = b0+0.01*z*b0/r**2
        br = -0.01*bz/r**2/2
        bphi = b0*0.01*phi/r
        bx = br*dcos(phi)
        by = br*dsin(phi)
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
c Use this section to plot orbit of single particle
c
        ione = idint(1.d0/fi+0.5d0)
        itwo = (i/ione)*ione
        if(itwo.eq.i)write(15,fmt='(1p2d11.3)')r2(k),z2(k)
        if(itwo.eq.i)write(16,fmt='(1p6d11.3)')t,x2(k),y2(k),r2(k),z2(k)
     &  ,errpphi

c        particle loop
c
        enddo
c
c End time loop
c
        enddo
c
c Output final positions & velocities
c
        t=9.0
        x2(1) = 9.0
        y2(1) = 9.0
        r2(1) = 9.0
        z2(1) = 9.0
        errpphi = 9.0
c        write(15,fmt='(1p4d11.3)')t,r2(1),z2(1),errpphi
        END
