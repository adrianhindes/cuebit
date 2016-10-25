        program CUEBIT
        implicit none
        integer :: i,j,k,nt,n,out1,out2,idum,zzz,mm
        double precision :: e,rm,pi,tpi,za,aa,r0,b0,th2
        double precision :: dt,dt1,det,phi2,psi0,tau,epsilon0
        double precision :: x,y,z,timp,vimp,bx,by,bz,ax,ay,az
        double precision :: r,phi,bphi,br,t,vx,vy,vz,vrotx,vroty
        double precision :: as,bs,vth,g(3),cx,cy,cz
        double precision :: vxm,vym,vzm,tx,ty,tz,nop,rstart
        double precision,dimension(:),allocatable::x1,x2,y1,y2,z1,z2
        double precision,dimension(:),allocatable::vx1,vx2,vy1,vy2
        double precision,dimension(:),allocatable::vz1,vz2,vsqa,vsq
        double precision,dimension(:),allocatable::phi1,r1,r2
        character(20) :: rch,zch,thch,phich,vxch,vych,vzch,psich
        double precision :: rval1,rval2,loglambda
        double precision :: fi,gamma,err0,err1,derr,psi
        double precision :: rb,omega,er,ephi,ez,ex,ey,psi1
        double precision :: temppsi,npsi,temp0,temp1,n0,n1
        double precision :: coeffa,coeffb,coeffc,xroot1,xroot2
        double precision :: nnought,iii,ntt
        double precision :: t9,r9,z9,psi9,rbar,drs,drbar
        double precision :: ziprime,miprime,fitwo
        double precision,dimension(:),allocatable::kkk
        integer :: ntactual,i1,i2,radius
        integer :: pchoice,profilechoice
        integer :: ione,itwo
        double precision g05ddf
        external function g05ddf,g05caf
        external g05cbf
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc                                                                        ccc
ccc This program solves the Lorentz force equation for charged particles   ccc
ccc in a MAST-like tokamak with Freidberg poloidal flux surfaces. The      ccc
ccc plasma is specified to be rotating co-currently or counter-currently   ccc
ccc in the toroidal direction.                                             ccc
ccc                                                                        ccc
ccc An implicit energy-conserving scheme is used. Iteration is used        ccc
ccc to make the overall scheme 2nd order accurate.                         ccc
ccc                                                                        ccc
ccc Collisions are taken into account by calculating the collision time    ccc
ccc from radial profiles for T and n. The test particles are initially in  ccc
ccc the outer midplane of the plasma with a Maxwellian temp. distribution. ccc
ccc                                                                        ccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c Last modified 12/03/10 by Ken McClements
c
c Initiate input & output files
c
        open(7,file="inputdata.dat",status='unknown')
        open(8,file="moments.dat",status='unknown')
c
        rch="    r(m) "
        zch="    z(m) "
        thch="    theta"
        phich="     phi"
        psich="     psi"
        vxch="      vx"
        vych="      vy "
        vzch="      vz"
c
c Electron charge, proton mass, pi & mu0
c
        e=1.602d-19
        rm=1.672d-27
        epsilon0=8.85d-12
        pi=4.d0*atan(1.d0)
        tpi=2.d0*pi
        out2=0
c
        vxm=0.0
        vym=0.0
        vzm=0.0
        tx=0.0
        ty=0.0
        tz=0.0       
c
c Read in parameters
c
c        write(6,fmt='("MAST-like parameters? Y=1 N=2")')
        read(7,*)pchoice
        if(pchoice.eq.2)goto 15
        r0=0.964
        rb=0.93
        b0=0.4
        psi0=0.9
        gamma=0.8
        ziprime=1.0
        miprime=2.0*rm
c
c Density/temperature parameters
c
        temp0=1
        temp1=0.1
        n0=5d19
        n1=1d19
c
c Coulomb logarithm value
c
        loglambda=15
        ephi=-0.3d0
c
c        write(6,fmt='("Omega, no. ions, impurity T (keV), Z and A?")')
        read(7,*)omega,n,timp,za,aa
        nop=n
        goto 20
c 15     write(6,fmt='("Values: Z,A,R0,Rb,B0,psi0,gamma,omega")')
 15     read(7,*)za,aa,r0,rb,b0,psi0,gamma,omega
c        write(6,fmt='("Timp,temp0,temp1,n,n0,n1,loglambda")')
        read(7,*)timp,temp0,temp1,n,n0,n1,loglambda
c
 20     coeffa=(psi0*gamma)/8
        coeffb=-(psi0*gamma*r0**2)/4
        coeffc=(psi0*gamma)*(r0**4-rb**4)/8
        xroot1=(-coeffb+sqrt(coeffb**2-4*coeffa*coeffc))/(2*coeffa)
        xroot2=(-coeffb-sqrt(coeffb**2-4*coeffa*coeffc))/(2*coeffa)
        rval1=sqrt(xroot1)
        rval2=sqrt(xroot2)
        timp=timp*e*1.e3
        vimp=sqrt(2.0*timp/aa/rm)
        allocate(x1(n),x2(n),y1(n),y2(n),z1(n),z2(n))
        allocate(vx1(n),vx2(n),vy1(n),vy2(n),vz1(n),vz2(n))
        allocate(phi1(n),r1(n),r2(n),kkk(n),vsqa(n),vsq(n))
c
c Input time step (in units of Larmor period at magnetic axis) & number of
c time steps. Fraction output refers to "boundary.dat" file, which
c tracks the number of particles still confined at each timestep: long
c simulations should have this value set to e.g. 0.1, 0.01 etc.
c
c        write(6,30)
c 30     format('Timestep,no. of steps,frac. output,seed')
        read(7,*)dt1,nt,fi,idum
        dt=2.d0*pi*aa*rm/za/e/b0*dt1
        mm=nt*fi
        write(8,36)mm
 36     format(i6)
c
c Select position to launch particles from
c
        read(7,*)radius
        if(radius.eq.1)rstart=rval1
        if(radius.eq.2)rstart=r0
c        if(radius.eq.3)rstart=r0+((rval1-r0)/2.0)       
        if(radius.eq.3)rstart=r0+(4.0*(1.3395-r0)/8.0)       
c
c Select temp/density profiles
c
c        write(6,fmt='("Choose temp./density profile (1-6!)")')
        read(7,*)profilechoice
c
c Initialise seed
c
        call g05cbf(idum)
c
c Initialise particle positions & velocities using NAG routines
c
        as=0.0
        bs=vimp/sqrt(2.0)
        do k=1,n
        r1(k)=rstart
        z1(k)=0.0
        phi1(k)=0.0
        vx1(k)=g05ddf(as,bs)
        vy1(k)=g05ddf(as,bs)+omega*r1(k)
c        vy1(k)=g05ddf(as,bs)
        vz1(k)=g05ddf(as,bs)
c        vx1(k)=0.0
c        vy1(k)=-2.8d4
c        vz1(k)=2.3d4
        vsqa(k)=vx1(k)**2+vy1(k)**2+vz1(k)**2
        psi=psi0*((gamma/8.0*((r1(k)**2-r0**2)**2-rb**4))
     &  +((1-gamma)/2.0*r1(k)**2*z1(k)**2))
        err0=za*e*omega*psi+0.5*aa*rm*vsqa(k)
        enddo
c
c Commence time loop
c
        t=0.d0
        do i=2,nt
        zzz=0
        out1=0
        t=t+dt
c        psi0=psi0-dt*r0*ephi
        iii=dble(i)
        ntt=dble(nt)
c
c Commence particle loop
c
        rbar=0.0
        drs=0.0
        do k=1,n
        if(kkk(k).eq.1)out1=out1+1
        if(kkk(k).eq.1)cycle
        x1(k)=r1(k)*cos(phi1(k))
        y1(k)=r1(k)*sin(phi1(k))
c
c Evaluate magnetic field
c
        bphi=b0*r0/r1(k)
        br=-psi0*z1(k)*r1(k)*(1-gamma)
        bx=br*dcos(phi1(k))-bphi*dsin(phi1(k))
        by=br*dsin(phi1(k))+bphi*dcos(phi1(k))
        bz=psi0/r1(k)*((gamma/8.0*(4.0*r1(k)**3-4.0*r0**2*r1(k)))
     &  +(1-gamma)*r1(k)*z1(k)**2)
c
c Evaluate electric field
c
        er=-omega*r1(k)*bz-2.0*rm*omega**2*r1(k)/(2.0*e)
        ez=omega*r1(k)*br
        ex=er*dcos(phi1(k))-ephi*sin(phi1(k))
        ey=er*dsin(phi1(k))+ephi*cos(phi1(k))
c
c Evaluate collision time via temp./density gradients
c
        psi=psi0*((gamma/8.0*((r1(k)**2-r0**2)**2-rb**4))
     &  +((1-gamma)/2.0*r1(k)**2*z1(k)**2))
        psi1=-(gamma*psi0*rb**4)/8.0
c
c Select desired temperature/density profile: 1-6
c
        if(profilechoice.eq.1)goto 41
        if(profilechoice.eq.2)goto 42
        if(profilechoice.eq.3)goto 43
        if(profilechoice.eq.4)goto 44
        if(profilechoice.eq.5)goto 45
        if(profilechoice.eq.6)goto 46 

 41     temppsi=(temp0-temp1)*(psi/psi1)+temp1
        npsi=(n0-n1)*(psi/psi1)+n1
        goto 47
c
 42     temppsi=(temp0-temp1)*(psi/psi1)+temp1
        npsi=(n0-n1)*sqrt(psi/psi1)+n1
        goto 47
c
 43     temppsi=(temp0-temp1)*sqrt(psi/psi1)+temp1
        npsi=(n0-n1)*(psi/psi1)+n1
        goto 47
c
 44     temppsi=(temp0-temp1)*(psi/psi1)+temp1
        nnought=(n0-n1)*(psi/psi1)+n1
        npsi=nnought*exp(2*rm*omega**2*(r1(k)**2-r0**2)/(4*temppsi))
        goto 47
c
 45     temppsi=(temp0-temp1)*(psi/psi1)+temp1
        nnought=(n0-n1)*sqrt(psi/psi1)+n1
        npsi=nnought*exp(2*rm*omega**2*(r1(k)**2-r0**2)/(4*temppsi))
        goto 47
c 
 46     temppsi=(temp0-temp1)*sqrt(psi/psi1)+temp1
        nnought=(n0-n1)*sqrt(psi/psi1)+n1
        npsi=nnought*exp(2*rm*omega**2*(r1(k)**2-r0**2)/(4*temppsi))
c
 47     temppsi=temppsi*e*1.e3
c
        tau=(aa*rm*6*sqrt(2.0)*pi**1.5*epsilon0**2*temppsi**1.5)/
     &  (sqrt(miprime)*za**2*ziprime**2*e**4*npsi*loglambda)
        vth=sqrt(2.0*temppsi/aa/rm)
c
c Select random numbers for collision term using NAG routines
c
        do j = 1,3
        as=0.d0
        bs=vth*sqrt(1.0/tau/dt)
        g(j)=g05ddf(as,bs)
        enddo
        x2(k)=x1(k)
        y2(k)=y1(k)
        z2(k)=z1(k)
        vx2(k)=vx1(k)
        vy2(k)=vy1(k)
        vz2(k)=vz1(k)
c
c Commence iteration loop
c
        do j=1,3
        x=(x1(k)+x2(k))/2.d0
        y=(y1(k)+y2(k))/2.d0
        z=(z1(k)+z2(k))/2.d0
        vx=(vx1(k)+vx2(k))/2.d0
        vy=(vy1(k)+vy2(k))/2.d0
        vz=(vz1(k)+vz2(k))/2.d0
        r=sqrt(x**2+y**2)
        phi=atan(y/x)
        if(x.lt.0.d0)phi=phi+pi
        if(x.gt.0.d0.and.y.lt.0.d0)phi=phi+tpi
        if(phi.gt.tpi)phi=phi-tpi
c
c Evaluate magnetic field components
c
        bphi=b0*r0/r
        br=-psi0*z*r*(1-gamma)
        bx=br*dcos(phi)-bphi*dsin(phi)
        by=br*dsin(phi)+bphi*dcos(phi)
        bz=psi0/r*((gamma/8.0*(4.0*r**3-4.0*r0**2*r))
     &  +(1-gamma)*r*z**2)
c
c Evaluate electric field components
c
        er=-omega*r*bz-2.0*rm*omega**2*r/(2.0*e)
        ephi=-0.3d0
        ez=omega*r*br
        ex=er*dcos(phi)-ephi*sin(phi)
        ey=er*dsin(phi)+ephi*cos(phi)
c
c Rotational velocities here
c
        vrotx=-omega*r*dsin(phi)
        vroty=omega*r*dcos(phi)
c
c Evaluate quantities appearing in matrix relating vi and vi+1
c
        ax=za*e*bx/2.d0/aa/rm*dt
        ay=za*e*by/2.d0/aa/rm*dt
        az=za*e*bz/2.d0/aa/rm*dt
        det=1.d0+ax**2+ay**2+az**2
        cx=-dt*(vx-vrotx)/tau+dt*g(1)+za*e*dt/aa/rm*ex
        cy=-dt*(vy-vroty)/tau+dt*g(2)+za*e*dt/aa/rm*ey
        cz=-dt*vz/tau+dt*g(3)+za*e*dt/aa/rm*ez
c        cx=za*e*dt/aa/rm*ex
c        cy=za*e*dt/aa/rm*ey
c        cz=za*e*dt/aa/rm*ez
c
c Advance equation dvx/dt = Ze/m(vy*Bz-vZ*By)
c
        vx2(k)=(1.d0+ax**2-ay**2-az**2)*vx1(k)
        vx2(k)=vx2(k)+2.d0*(ay*ax+az)*vy1(k)
        vx2(k)=vx2(k)+2.d0*(az*ax-ay)*vz1(k)
        vx2(k)=vx2(k)+(1.0+ax**2)*cx
        vx2(k)=vx2(k)+(ax*ay+az)*cy
        vx2(k)=vx2(k)+(ax*az-ay)*cz
        vx2(k)=vx2(k)/det
c
c Advance equation dvy/dt = Ze/m(vz*Bx-vx*By)
c
        vy2(k)=2.d0*(ay*ax-az)*vx1(k)
        vy2(k)=vy2(k)+(1.d0-ax**2+ay**2-az**2)*vy1(k)
        vy2(k)=vy2(k)+2.d0*(ay*az+ax)*vz1(k)
        vy2(k)=vy2(k)+(ax*ay-az)*cx
        vy2(k)=vy2(k)+(1.0+ay**2)*cy
        vy2(k)=vy2(k)+(ay*az+ax)*cz
        vy2(k)=vy2(k)/det
c
c Advance equation dvz/dt = Ze/m(vx*By-vy*Bx)
c
        vz2(k)=2.d0*(az*ax+ay)*vx1(k)
        vz2(k)=vz2(k)+2.d0*(ay*az-ax)*vy1(k)
        vz2(k)=vz2(k)+(1.d0-ax**2-ay**2+az**2)*vz1(k)
        vz2(k)=vz2(k)+(ax*az+ay)*cx
        vz2(k)=vz2(k)+(ay*az-ax)*cy
        vz2(k)=vz2(k)+(1.0+az**2)*cz
        vz2(k)=vz2(k)/det
c
c Advance equation dx/dt = vx
c
        x2(k)=x1(k)+(vx2(k)+vx1(k))/2.d0*dt
c
c Advance equation dy/dt = vy
c
        y2(k)=y1(k)+(vy2(k)+vy1(k))/2.d0*dt
c
c Advance equation dz/dt = vz
c
        z2(k)=z1(k)+(vz2(k)+vz1(k))/2.d0*dt
c
        r2(k)=sqrt(x2(k)**2+y2(k)**2)
        enddo
c
        x1(k)=x2(k)
        y1(k)=y2(k)
        z1(k)=z2(k)
        r1(k)=r2(k)
        phi1(k)=atan(y1(k)/x1(k))
        if(x1(k).lt.0.d0)phi1(k)=phi1(k)+pi
        if(x1(k).gt.0.d0.and.y1(k).lt.0.d0)phi1(k)=phi1(k)+tpi
        if(phi1(k).gt.tpi)phi1(k)=phi1(k)-tpi
        psi=psi0*((gamma/8.0*((r1(k)**2-r0**2)**2-rb**4))
     &  +((1-gamma)/2.0*r1(k)**2*z1(k)**2))
        if(psi.gt.0)out1=out1+1
        if(psi.gt.0)kkk(k)=1
        vx1(k)=vx2(k)
        vy1(k)=vy2(k)
        vz1(k)=vz2(k)
        vsq(k)=vx2(k)**2+vy2(k)**2+vz2(k)**2
        err1=za*e*omega*psi+0.5*aa*rm*vsq(k)
        derr=(err1-err0)
        rbar=rbar+r1(k)
c
c End particle loop
c
        x1(k)=x2(k)
        y1(k)=y2(k)
        z1(k)=z2(k)
        r1(k)=r2(k)
        vx1(k)=vx2(k)
        vy1(k)=vy2(k)
        vz1(k)=vz2(k)
c
        enddo
        rbar=rbar/float(n)
        drbar=rbar-rstart
        do k=1,n
           drs=drs+(r1(k)-rbar)**2
        enddo
        drs=drs/float(n)
        i1=idint(1.D0/fi+0.5d0)
        i2=(i/i1)*i1
        if(i2.eq.i)goto 152
        goto 153
 152     write(8,10)t,drbar,drs
 10     format(1p3d12.5)
c
c End time loop
c
 153    i1=idint(1.D0/fi+0.5d0)
        i2=(i/i1)*i1
        ntactual=n-out1
        enddo
        t9=9.0
        r9=9.0
        z9=9.0
        psi9=9.0
c
c Output final positions & velocities
c
        vxm=0.0
        vym=0.0
        vzm=0.0
        tx=0.0
        ty=0.0
        tz=0.0
        do k=1,n  
        vxm=vxm+vx2(k)
        vym=vym+vy2(k)
        vzm=vzm+vz2(k)
        tx=tx+vx2(k)**2
        ty=ty+vy2(k)**2
        tz=tz+vz2(k)**2
        phi2=atan(y2(k)/x2(k))
        if(x2(k).lt.0.d0)phi2=phi2+pi
        if(x2(k).gt.0.d0.and.y2(k).lt.0.d0)phi2=phi2+tpi
        if(phi2.gt.tpi)phi2=phi2-tpi
        phi2=180.0/pi*phi2
        th2=atan(z2(k)/(r2(k)-r0))
        if(r2(k).lt.r0)th2=th2+pi
        if(r2(k).gt.r0.and.z2(k).lt.0.d0)th2=th2+tpi
        if(th2.gt.tpi)th2=th2-tpi
        th2=180.d0/pi*th2
        psi=psi0*((gamma/8*((r1(k)**2-r0**2)**2-rb**4))
     &  +((1-gamma)/2*r1(k)**2*z1(k)**2))
        if(psi.gt.0)out2=out2+1
c        write(8,60)r2(k),z2(k),th2,phi2,psi,vx2(k),vy2(k),vz2(k)
c 60     format(1p8d11.3)
        enddo
        t=9.0
        r2(k)=9.0
        z2(k)=9.0
        psi=9.0
        t=9.0
        drbar=9.0
        drs=9.0
        write(8,10)t,drbar,drs
c
        END
