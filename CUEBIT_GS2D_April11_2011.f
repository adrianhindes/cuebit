	module declare
	implicit none
	double precision :: e,rm,tpi,za,aa,vr1,vphi1,xcen
	double precision :: vx,vy,vz,psipa,ph11,mu11,ycen,pphi
	double precision :: vz0,dt,fi,vsq,t,err,errphi,b1,det,phi2
	double precision :: x,y,z,phideg,psi1,pphi0,vr2,dt1
	double precision :: bx,by,ax,ay,az,r,phi,errpphi,phid
	double precision :: vphi2,vsq2,psi2,pphi2,rco,deltabph,deltabr
	double precision :: rg1,rg2,zg1,zg2,deltar,deltaz,vphi
	double precision :: psirg1,psirg2,psizg1,psizg2,r11,th11
        double precision,dimension(:),allocatable::x1,x2,y1,y2,z1,z2
        double precision,dimension(:),allocatable::vx1,vx2,vy1,vy2
        double precision,dimension(:),allocatable::vz1,vz2
        double precision,dimension(:),allocatable::phi1,r1,r2	
	integer :: i,j,nt,i1,i2,nco,n,iseed,idum,k
        integer :: ione,itwo,nrem,nlost
        integer,dimension(:),allocatable::nflag
	character(len=20) :: tch,rch,phich,zch,deech,depphich
	end module declare 

      program CUEBIT_GS2D_April11_2011
      use declare
      implicit none
      real :: terp2
      external coeff1,terp1,coeff2,terp2
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc                                                                      ccc
ccc This program reads in a MAST equilibrium reconstruction generated by ccc
ccc the code GS2D. A TF ripple is superimposed on this equilibrium.      ccc
ccc A subroutine package splines.f is used to interpolate the fields     ccc
ccc using bicubic splines.                                               ccc 
ccc                                                                      ccc
ccc This is a multi-particle version of the code.                        ccc  
ccc                                                                      ccc
ccc To compile this code type:                                           ccc
ccc                                                                      ccc 
ccc gfortran CUEBIT_GS2D_April11_2011.f splines.f -fdefault-real-8       ccc
ccc    -o CUEBIT_GS2D_April11_2011                                       ccc
ccc                                                                      ccc
ccc Last modified by Ken McClements 11/4/11.                             ccc
ccc                                                                      ccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      double precision,dimension(:),allocatable::ps,amin,q,f,p,pp
      double precision,dimension(:),allocatable::rsep,zsep,rgrid,zgrid
      double precision,dimension(:,:),allocatable:: pxx,pyy,pxxyy
      double precision,dimension(:,:),allocatable:: bxx,byy,bxxyy
      double precision,dimension(:,:),allocatable:: br,bphi,bz,psi
      double precision :: r0,a,rmag,zmag,pi,dz,dr
      double precision :: psmin,psip,b0,ippsi
      double precision :: brr,bphp,bzz
      integer :: nsu,nsep,nr,nz,nm
      character(len=80) :: line
      integer :: iop(2),isrch,int1
      integer :: itab(3),lenw,lenf,lenwk
      integer :: ibd(4),ixd,iyd,isav,jsav,idm
      double precision :: w(100),tab(3),wk(301)
c
c This part of the program reads in the MAST equilibrium
c
      open(1,file="gs2.dat",status='unknown')
      do i=1,6
         read(1,fmt='(a80)') line
      enddo
      read(1,*) r0,a,rmag,zmag
      read(1,fmt='(a80)') line ; read(1,*) psmin,psip,b0,ippsi
      read(1,fmt='(a80)') line ; read(1,*) nsu
      allocate(ps(nsu),amin(nsu),q(nsu),f(nsu),p(nsu),pp(nsu))
      read(1,fmt='(a80)') line ; read(1,fmt='(1p,8e12.4)') ps
      read(1,fmt='(a80)') line ; read(1,fmt='(1p,8e12.4)') amin
      read(1,fmt='(a80)') line ; read(1,fmt='(1p,8e12.4)') q
      read(1,fmt='(a80)') line ; read(1,fmt='(1p,8e12.4)') f
      read(1,fmt='(a80)') line ; read(1,fmt='(1p,8e12.4)') p
      read(1,fmt='(a80)') line ; read(1,fmt='(1p,8e12.4)') pp
      read(1,fmt='(T22,I6)') nsep
      allocate(rsep(nsep),zsep(nsep))
      read(1,fmt='(a80)') line ; read(1,fmt='(1p,8e12.4)') rsep
      read(1,fmt='(a80)') line ; read(1,fmt='(1p,8e12.4)') zsep
      read(1,fmt='(a80)') line ; read(1,*) nr,nz
      allocate(rgrid(nr),zgrid(nz),psi(nr,nz),bphi(nr,nz))
      allocate(br(nr,nz),bz(nr,nz))
      allocate(pxx(nr,nz),pyy(nr,nz),pxxyy(nr,nz))
      allocate(bxx(nr,nz),byy(nr,nz),bxxyy(nr,nz))
      read(1,fmt='(a80)') line ; read(1,*) rgrid
      read(1,fmt='(a80)') line ; read(1,*) zgrid
      read(1,fmt='(a80)') line ; read(1,fmt='(1p,8e12.4)') psi
      close(1)
c
c Evaluate BR 
c
      pi=4.d0*atan(1.d0)
      do i=1,nr
         do j=1,nz
            if (j.eq.1) then 
	      dz=zgrid(j+1)-zgrid(j)   
              br(i,j)=(psi(i,j+1)-psi(i,j))/dz/2.d0/pi/rgrid(i)
	    else if (j.eq.nz) then
	      dz=zgrid(j)-zgrid(j-1)
              br(i,j)=(psi(i,j)-psi(i,j-1))/dz/2.d0/pi/rgrid(i)
            else       
              dz=zgrid(j+1)-zgrid(j-1)
              br(i,j)=(psi(i,j+1)-psi(i,j-1))/dz/2.d0/pi/rgrid(i)
	    end if
         enddo
      enddo
c
c
c Evaluate BZ 
c
      do j=1,nz
         do i=1,nr
            if (i.eq.1) then
              dr=rgrid(i+1)-rgrid(i)
              bz(i,j)=-(psi(i+1,j)-psi(i,j))/dr/2.d0/pi/rgrid(i)
            else if (i.eq.nr) then
              dr=rgrid(i)-rgrid(i-1)
              bz(i,j)=-(psi(i,j)-psi(i-1,j))/dr/2.d0/pi/rgrid(i)
	    else
              dr=rgrid(i+1)-rgrid(i-1)
              bz(i,j)=-(psi(i+1,j)-psi(i-1,j))/dr/2.d0/pi/rgrid(i)
	    end if
         enddo
      enddo
c
c
c Evaluate Bphi 
c
      do i=1,nr
         do j=1,nz
         iop(1)=2
         iop(2)=2
         w(1)=0.d0
         w(nsu)=0.d0
         itab(1)=1
         itab(2)=1
         itab(3)=1
         isrch=1
         int1=1
         call coeff1(nsu,ps,lenf,f,lenw,w,iop,int1,wk)
         call terp1(nsu,ps,lenf,f,lenw,w,psi(i,j),int1,tab,itab,isrch)
         bphi(i,j)=tab(1)/rgrid(i)
         enddo
      enddo
c
      open(1,file="gs2_copy.dat",status='unknown')
      write(1,fmt='("GS2 input file")')
      write(1,fmt='("r0,a,rmag,zmag")')
      write(1,fmt='(1p,4e12.4)')r0,a,rmag,zmag
      write(1,fmt='("psmin,psip,b0,ippsi")')
      write(1,fmt='(1p,4e12.4)')psmin,psip,b0,ippsi
      write(1,fmt='("nfs"/I6)') nsu
      write(1,fmt='("Psi on 1d grid (for FS quantities) (T)")')
      write(1,fmt='(1p,8e12.4)') ps
      write(1,fmt='("amin (m)")')
      write(1,fmt='(1p,8e12.4)') amin
      write(1,fmt='("q")')
      write(1,fmt='(1p,8e12.4)') q
      write(1,fmt='("f =r B_phi (Tm)")')
      write(1,fmt='(1p,8e12.4)') f
      write(1,fmt='("p (Pa)")')
      write(1,fmt='(1p,8e12.4)') p
      write(1,fmt='("dp/dpsi (Pa/Wb)")')
      write(1,fmt='(1p,8e12.4)') pp
      write(1,fmt='("No of points on LCFS=",I6)') nsep
      write(1,fmt='("r(j) (m) on LCFS")')
      write(1,fmt='(1p,8e12.4)') rsep
      write(1,fmt='("z(j) (m) on LCFS")')
      write(1,fmt='(1p,8e12.4)') zsep
      write(1,fmt='("NR",T14,"NZ"/2I6)') NR, NZ
      write(1,fmt='("rgrid (m)")')
      write(1,fmt='(1p,8e12.4)') rgrid
      write(1,fmt='("zgrid (m)")')
      write(1,fmt='(1p,8e12.4)') zgrid
      write(1,fmt='("Psi on grid (Wb) : NB Br=(1/2pi r)*dpsi/dz")')
      write(1,fmt='(1p,8e12.4)') psi
      close(1)
c
c Evaluate coefficients for interpolation of fields
c
	deltar=1.d-6
	deltaz=1.d-6
        idm=nr
	lenwk=301
	ibd(1)=4
	ibd(2)=4
	ibd(3)=4
	ibd(4)=4
        ixd=0
        iyd=0
	isrch=1
        call coeff2(nr,rgrid,nz,zgrid,psi,pxx,pyy,pxxyy,idm,ibd,
     1      lenwk,wk)
        idm=nr
	lenwk=196
	ibd(1)=4
	ibd(2)=4
	ibd(3)=4
	ibd(4)=4
        call coeff2(nr,rgrid,nz,zgrid,bphi,bxx,byy,bxxyy,idm,ibd,
     1      lenwk,wk)
c
c This part of the program tracks the orbits of a population 
c of charged particles
c
c Physical & mathematical constants
c
	e=1.602d-19
	rm=1.672d-27
	pi=4.d0*atan(1.d0)
	tpi=2.d0*pi
c
c Initiate particle input file
c
      open(7,file="GS2_particle_data.dat",status='unknown')
      open(8,file="gs2_orbit.dat",status='unknown')
      open(9,file="gs2_RZ_orbit.dat",status='unknown')
c
c Input mass number, charge state, time step in cyclotron periods,
c number of time steps, & fraction of time steps output. 
c
	read(7,fmt='(2f5.1,f6.3,i9,f7.4)')aa,za,dt1,nt,fi
c
c Evaluate time step in seconds
c
	dt=2.d0*pi*aa*rm/za/e/abs(b0)*dt1
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
        allocate(x1(n),x2(n),y1(n),y2(n),z1(n),z2(n))
        allocate(vx1(n),vx2(n),vy1(n),vy2(n),vz1(n),vz2(n))
        allocate(phi1(n),r1(n),r2(n),nflag(n))
c
c Initialise positions & velocities of particles   
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
           vy1(k) = 2.4e6
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
	write(12,fmt='("Initial particle positions")')
	write(12,fmt='()')
	write(12,fmt='("    R(m)     Z(m)")')
c
c Output initial particle positions
c
        do k = 1,nrem
           write(12,fmt='(1p2e10.3)')r1(k),z1(k)
        enddo 
c
c Commence time loop
c
	t=0.d0
        vphi = -vy1(1)*cos(phi1(1))
        do i = 2,nt
        t = t+dt
c
c Commence particle loop (n > 1 only in multi-particle versions of code)
c
        do k = 1,nrem
           if(nflag(k).eq.1)cycle
c
c Evaluate magnetic field
c
	rg1=r1(k)-deltar
	rg2=r1(k)+deltar
	zg1=z1(k)-deltaz
	zg2=z1(k)+deltaz
c 
	psirg1=terp2(rg1,z1(k),nr,rgrid,nz,zgrid,psi,pxx,pyy,pxxyy,
     1      idm,ixd,iyd,isrch,isav,jsav)
	psirg2=terp2(rg2,z1(k),nr,rgrid,nz,zgrid,psi,pxx,pyy,pxxyy,
     1      idm,ixd,iyd,isrch,isav,jsav)
	psizg1=terp2(r1(k),zg1,nr,rgrid,nz,zgrid,psi,pxx,pyy,pxxyy,
     1      idm,ixd,iyd,isrch,isav,jsav)
	psizg2=terp2(r1(k),zg2,nr,rgrid,nz,zgrid,psi,pxx,pyy,pxxyy,
     1      idm,ixd,iyd,isrch,isav,jsav)
	brr=(psizg2-psizg1)/2.d0/deltaz/r1(k)/2.d0/pi
	bzz=-(psirg2-psirg1)/2.d0/deltar/r1(k)/2.d0/pi
	bphp=terp2(r1(k),z1(k),nr,rgrid,nz,zgrid,bphi,bxx,byy,bxxyy,
     1      idm,ixd,iyd,isrch,isav,jsav)
c
c Evaluate ripple amplitude
c
        deltabph = (b0*rmag/r1(k))*(r1(k)/rco)**nco*cos(nco*phi1(k))
	deltabr = (b0*rmag/r1(k))*(r1(k)/rco)**nco*sin(nco*phi1(k)) 
c
	bphp = bphp+deltabph
	brr = brr+deltabr
	bx=brr*cos(phi1(k))-bphp*sin(phi1(k))
	by=brr*sin(phi1(k))+bphp*cos(phi1(k))
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
c Evaluate magnetic field components
c
	rg1=r-deltar
	rg2=r+deltar
	zg1=z-deltaz
	zg2=z+deltaz
c 
	psirg1=terp2(rg1,z,nr,rgrid,nz,zgrid,psi,pxx,pyy,pxxyy,idm,
     1      ixd,iyd,isrch,isav,jsav)
	psirg2=terp2(rg2,z,nr,rgrid,nz,zgrid,psi,pxx,pyy,pxxyy,idm,
     1      ixd,iyd,isrch,isav,jsav)
	psizg1=terp2(r,zg1,nr,rgrid,nz,zgrid,psi,pxx,pyy,pxxyy,idm,
     1      ixd,iyd,isrch,isav,jsav)
	psizg2=terp2(r,zg2,nr,rgrid,nz,zgrid,psi,pxx,pyy,pxxyy,idm,
     1      ixd,iyd,isrch,isav,jsav)
	brr=(psizg2-psizg1)/2.d0/deltaz/r/2.d0/pi
	bzz=-(psirg2-psirg1)/2.d0/deltar/r/2.d0/pi
	bphp=terp2(r,z,nr,rgrid,nz,zgrid,bphi,bxx,byy,bxxyy,idm,
     1      ixd,iyd,isrch,isav,jsav)
c
c Evaluate ripple amplitude for current iteration
c
        deltabph = (b0*rmag/r)*(r/rco)**nco*cos(nco*phi)
	deltabr = (b0*rmag/r)*(r/rco)**nco*sin(nco*phi)
	bphp = bphp+deltabph
	brr = brr+deltabr
c
	bx=brr*cos(phi)-bphp*sin(phi)
	by=brr*sin(phi)+bphp*cos(phi)
c
c Evaluate quantities appearing in matrix relating vi and vi+1
c
	ax=za*e*bx/2.d0/aa/rm*dt
	ay=za*e*by/2.d0/aa/rm*dt
	az=za*e*bzz/2.d0/aa/rm*dt
	det=1.d0+ax**2+ay**2+az**2
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
        vx1(k) = vx2(k)
        vy1(k) = vy2(k)
        vz1(k) = vz2(k)
        vphi = -vx1(k)*dsin(phi1(k))+vy1(k)*dcos(phi1(k))
c
c Test to determine whether particle is lost
c
	psipa = terp2(r1(k),z1(k),nr,rgrid,nz,zgrid,psi,pxx,pyy,pxxyy,
     1      idm,ixd,iyd,isrch,isav,jsav)
	pphi = aa*rm*r2(k)*vphi-za*e*psipa/tpi
	vsq = vx1(k)**2+vy1(k)**2+vz1(k)**2
	write(8,*)t,psipa,vsq
        if (psipa > psip) nlost = nlost+1
        if (psipa > psip) nflag(k) = 1
c        i1=idint(1.d0/fi+0.5d0)
c        i2=(i/i1)*i1
c        if(i2.eq.i)then
c	   write(9,fmt='(1p,2d11.3)')r2(k),z2(k)
c	end if 
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






























