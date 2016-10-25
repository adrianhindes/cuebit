c
c cubspl     from nssl                version   8          05/19/77
c
c package cubspl
c
c
c latest revision        september 1976
c
c purpose                to perform one and two-dimensional cubic spline
c                        interpolation with choice of boundary
c                        conditions.  the function and selected
c                        derivatives may be evaluated at any point where
c                        interpolation is required.  in the
c                        two-dimensional case, the given data points
c                        must be on a rectangular grid, which need not
c                        be equally spaced.  the package cubspl contains
c                        seven routines.
c
c                        subroutine coeff1
c                          computes the coefficients for one-dimensional
c                          cubic spline interpolation using one of two
c                          boundary conditions at each end of the range.
c                            . second derivative given at boundary.
c                            . first derivative given at boundary.
c                            . periodic boundary condition.
c                            . first derivative determined by fitting a
c                              cubic to the four points nearest to the
c                              boundary.
c
c                        subroutine terp1
c                          using the coefficients computed by coeff1,
c                          this routine evaluates the function and/or
c                          first and second derivatives at any point
c                          where interpolation is required.
c
c                        subroutine coeff2
c                          computes the coefficients for two-dimensional
c                          bicubic spline interpolation with the same
c                          choice of boundary conditions as for coeff1.
c
c                        function terp2
c                          using the coefficients produced by coeff2,
c                          this subroutine evaluates the function or a
c                          selected derivative at any point where
c                          two-dimensional interpolation is required.
c
c                        subroutine trip
c                          a simple periodic, tridiagonal linear
c                          equation solver used by coeff1.
c
c                        subroutine search
c                          performs a binary search in a one-dimensional
c                          floating point table arranged in ascending
c                          order.  this routine is called by terp1 and
c                          terp2.
c
c                        subroutine intrp
c                          given coefficients provided by coeff1 and the
c                          position of the interpolation point in the
c                          independent variable table, this routine
c                          performs one-dimensional interpolation for
c                          the function value, first and second
c                          derivative, as desired.  this routine is
c                          called by terp1 and terp2.
c
c access cards           *fortran,s=ulib,n=cubspl
c                        *cosy
c                          these cards access all seven routines of the
c                          package, cubspl.
c
c space required         total space for all seven routines is
c                        2243 (octal) = 1187 (decimal).
c
c entry points           coeff1, terp1, coeff2, terp2, trip, search,
c                        intrp
c
c special conditions     tables of independent variables must be
c                        arranged in ascending order.  for
c                        two-dimensional interpolation, the data points
c                        must lie on a rectangular mesh.
c
c common blocks          none
c
c i/o                    none
c
c precision              single
c
c required ulib          none
c routines
c
c specialist             cicely ridley, ncar, boulder, colorado  80303
c
c language               fortran
c
c history                this package is based on the routines
c                            la e 102a, spl1d1
c                            la e 103a, spl1d2
c                            la e 104a, spl2d1
c                            la e 105a, spl2d2
c                        of the los alamos cubic spline package written
c                        by thomas j. jordan and bertha fagan, 1968.
c                        the routines have been streamlined and
c                        standardized.  the algorithm for
c                        two-dimensional interpolation is considerably
c                        modified.
c
c
      subroutine coeff1 (n,x,lenf,f,lenw,w,iop,int1,wk)
c      implicit double precision (a-h,o-z)
c
c
c dimension of           x(n),f(int*(n-1)+1),w(int*(n-1)+1),iop(2),
c arguments              wk(3*n+1)
c
c latest revision        february 1974
c
c usage                  call coeff1 (n,x,f,w,iop,int1,wk)
c
c arguments
c
c on input               n
c                          the number of data points.  n must be at
c                          least 4.
c
c                        x
c                          table of n independent variable values in
c                          ascending order.  dimension of x in calling
c                          program must be at least n.
c
c                        f
c                          table of n corresponding dependent variable
c                          values.  the values are separated by interval
c                          int.  this is usually unity for
c                          one-dimensional interpolation.  other values
c                          are useful when coeff1 is incorporated in a
c                          two-dimensional interpolation scheme (see
c                          coeff2).  dimension of f in the calling
c                          program must be at least (int*(n-1)+1).
c
c                        iop
c                          two element integer array defining boundary
c                          conditions at x(1) and x(n) according to the
c                          following code.
c
c                          for iop(1)
c                          = 1  second derivative given at x(1).  place
c                               value of second derivative in w(1)
c                               before call to coeff1.
c                          = 2  first derivative given at x(1).  place
c                               value of first derivative in w(1) before
c                               call.
c                          = 3  periodic boundary condition.  x(1) and
c                               x(n) are equivalent points.  f(1) and
c                               f(int*(n-1)+1) are equal.
c                          = 4  the first derivative at x(1) is
c                               calculated by fitting a cubic to points
c                               x(1) through x(4).
c                          similarly, iop(2) defines the boundary
c                          condition at x(n).  when iop(2) = 1 (or 2),
c                          the value of the second (or first) derivative
c                          must be placed in w(int*(n-1)+1).  note that
c                          if iop(1) = 3, consistency demands that
c                          iop(2) = 3 also.
c
c                        int1
c                          spacing in f and w tables.  for
c                          one-dimensional interpolation this will
c                          usually be unity.
c
c                        wk
c                          work area of dimension at least (3*n+1).
c   100
c on output              w
c                          table of second derivatives corresponding to
c                          given x and f values.  the separation of
c                          tabular entries is int (see above).
c                          dimension of w in the calling program must be
c                          at least (int*(n-1)+1).
c
c                          the arrays x, f, w are used as input for the
c                          routine terp1, which performs interpolation
c                          at a given value of the independent variable.
c
c space required         654 (octal) = 428 (decimal)
c
c timing                 the timing is linearly proportional to n, the
c                        number of data points.  for 21 data points, the
c                        time on the ncar cdc 7600 was approximately
c                        .24 milliseconds.
c
c
      integer iop(2),ii
      integer i,i1,i2,in,index1,int1,j0,j1,j2,j3,j4,jm
      integer lenf,lenw,lenz,mk,ml,n,nn
      double precision x(n),f(lenf),w(lenw),wk(n,4)
      double precision a12,a13,a14,a23,a24,a34,b2,g1,g2,g3,g4
      double precision xxx,xxy,y2,terp
c
c arithmetic statement function used to locate entries in f and w arrays
c
      ii(index1)=(index1-1)*int1+1
c
c arithmetic statement function for interpolation
c
      terp(g1,g2,g3,g4) = (1.d0/a12+1.d0/a13+1.d0/a14)*g1-
     1          a13*a14/(a12*a23*a24)*g2+a12*a14/(a13*a23*a34)*g3-
     2          a12*a13/(a14*a24*a34)*g4
c the following call is for gathering statistics on library use at ncar
c
c start to set up tridiagonal system
c
      j0 = 1
      do 101 i=2,n
	 jm = j0
	 j0 = j0+int1
	 wk(i,1) = x(i)-x(i-1)
	 wk(i,2) = (f(j0)-f(jm))/wk(i,1)
	 wk(i,3) = wk(i,1)/6.d0
	 wk(i,1) = wk(i,1)/3.d0
  101 continue
      nn = n
      mk = iop(1)
      ml = iop(2)
c
c apply boundary conditions at boundary 1
c
      go to (102,103,104,105),mk
c
c second derivative given at boundary 1
c
  102 continue
      wk(2,2) = wk(3,2)-wk(2,2)-wk(2,3)*w(1)
      wk(2,3) = 0.d0
      wk(2,1) = wk(2,1)+wk(3,1)
      i1 = 2
      nn = nn-1
      go to 106
c
c first derivative given at boundary 1
c
  103 continue
      xxx=w(1)
      wk(1,2) = wk(2,2)-w(1)
      wk(2,2) = wk(3,2)-wk(2,2)
      wk(1,3) = 0.d0
      wk(1,1) = wk(2,1)
      wk(2,1) = wk(2,1)+wk(3,1)
      i1 = 1
      go to 106
c
c periodic boundary condition
c
  104 continue
      y2 = wk(2,2)
      b2 = wk(2,1)
      wk(2,2) = wk(3,2)-wk(2,2)
      wk(2,1) = wk(3,1)+wk(2,1)
      i1 = 2
      nn = nn-1
      go to 106
c
c first derivative at boundary 1 from 4 point interpolation.
c
  105 continue
      a12 = x(1)-x(2)
      a13 = x(1)-x(3)
      a14 = x(1)-x(4)
      a23 = x(2)-x(3)
      a24 = x(2)-x(4)
      a34 = x(3)-x(4)
      j1 = 1
      j2 = j1+int1
      j3 = j2+int1
      j4 = j3+int1
      w(1) = terp(f(j1),f(j2),f(j3),f(j4))
      go to 103
c
c compute tridiagonal arrays
c
  106 continue
      i2 = n-2
      do 107 i=3,i2
	 wk(i,2) = wk(i+1,2)-wk(i,2)
	 wk(i,1) = wk(i+1,1)+wk(i,1)
  107 continue
c
c apply boundary conditions at boundary 2.
c
      in = ii(n)
      go to (108,109,110,111),ml
c
c second derivative given at boundary 2.
c
  108 continue
      wk(n-1,2) = wk(n,2)-wk(n-1,2)-wk(n,3)*w(in)
      wk(n,3) = 0.d0
      wk(n-1,1) = wk(n-1,1)+wk(n,1)
      nn = nn-1
      go to 112
c
c first derivative given at boundary 2.
c
  109 continue
      xxy=w(in)
      wk(n-1,2) = wk(n,2)-wk(n-1,2)
      wk(n,2) = -wk(n,2)+w(in)
      wk(n-1,1) = wk(n-1,1)+wk(n,1)
      wk(1,4) = 0.d0
      go to 112
c
c periodic boundary condition
c
  110 continue
      wk(n-1,2) = wk(n,2)-wk(n-1,2)
      wk(n,2) = y2-wk(n,2)
      wk(n-1,1) = wk(n-1,1)+wk(n,1)
      wk(n,1) = wk(n,1)+b2
      wk(1,4) = wk(2,3)
      go to 112
c
c first derivative at boundary 2 from 4 point interpolation.
c
  111 continue
      a12 = x(n)-x(n-1)
      a13 = x(n)-x(n-2)
      a14 = x(n)-x(n-3)
      a23 = x(n-1)-x(n-2)
      a24 = x(n-1)-x(n-3)
      a34 = x(n-2)-x(n-3)
      j1 = in
      j2 = j1-int1
      j3 = j2-int1
      j4 = j3-int1
      w(in) = terp(f(j1),f(j2),f(j3),f(j4))
      go to 109
  112 continue
c
c solve tridiagonal system
c
      lenz=(nn-1)*int1+1
      call trip (nn,wk(i1,3),wk(i1,1),wk(i1+1,3),wk(i1,2),lenz,
     1           w(ii(i1)),int1)
      go to (114,114,113,114),mk
  113 continue
      w(1) = w(in)
  114 continue
      return
      end
      subroutine trip (n,a,b,c,y,lenz,z,int1)
c      implicit double precision (a-h,o-z)
      integer ik,index1,in,int1,j,k,lenz,n,nm1,nm2,ii
      double precision  a(n),b(n),c(n),y(n),z(lenz)
      double precision bn,den,v,yn,terp
c
c arithmetic statement function used to locate entries in array z.
c
      ii(index1)=(index1-1)*int1+1
c
c gaussian elimination
c
      bn = b(n)
      yn = y(n)
      v = c(n)
      y(1) = y(1)/b(1)
      a(1) = a(1)/b(1)
      b(1) = c(1)/b(1)
      nm2 = n-2
      do 101 j=2,nm2
	 den = b(j)-a(j)*b(j-1)
	 b(j) = c(j)/den
	 y(j) = (y(j)-a(j)*y(j-1))/den
	 a(j) = -a(j)*a(j-1)/den
	 bn = bn-v*a(j-1)
	 yn = yn-v*y(j-1)
	 v = -v*b(j-1)
  101 continue
      den = b(n-1)-a(n-1)*b(n-2)
      b(n-1) = (c(n-1)-a(n-1)*a(n-2))/den
      y(n-1) = (y(n-1)-a(n-1)*y(n-2))/den
      bn = bn-v*a(n-2)
      yn = yn-v*y(n-2)
      v = a(n)-v*b(n-2)
c
c back substitution
c
      z(ii(n)) = (yn-v*y(n-1))/(bn-v*b(n-1))
      z(ii(n-1)) = y(n-1)-b(n-1)*z(ii(n))
      nm1 = n-1
      in = ii(n)
      do 102 j=2,nm1
	 k = n-j
	 ik = ii(k)
	 z(ik) = y(k)-b(k)*z(ik+int1)-a(k)*z(in)
  102 continue
      return
      end
      subroutine search (xbar,x,n,i)
c      implicit double precision (a-h,o-z)
      integer i,k,m,n,nm1
      double precision x(n),b,xbar
      data b/.69314718d0/
c
c if xbar is outside range of x table extrapolate
c
      if (xbar .gt. x(2)) go to 101
      i = 1
      return
  101 continue
      if (xbar .lt. x(n-1)) go to 102
      i = n-1
      return
  102 continue
c
c find maximum power of two less than n
c
      m = int((log(dble(n)))/b)
      i = 2**m
      if (i .ge. n) i = i/2
      k = i
      nm1 = n-1
c
c conduct binary search.
c
  103 continue
      k = k/2
      if (xbar .ge. x(i)) go to 104
      i = i-k
      go to 103
  104 continue
      if (xbar .le. x(i+1)) return
      i = min(i+k,nm1)
      go to 103
      end
      subroutine intrp (n,x,lenf,f,lenw,w,y,i,int1,tab,itab)
c      implicit double precision (a-h,o-z)
      integer i,i0,index1,int1,ip,lenf,lenw,n,ii
      double precision x(n),f(lenf),w(lenw),tab(3)
      double precision a,b,c,fl0,fl0sq,flk,flki,flp,flpsq
      double precision sixth,y
c     logical         itab(3)
      integer         itab(3)
      data sixth /0.1666666666667d0/
      save sixth
c
c arithmetic statement function used to locate entries in f and w arrays
c
      ii(index1)=(index1-1)*int1+1
c
c perform interpolation or extrapolation
c
      flk = x(i+1)-x(i)
      flki = 1.d0/flk
      flp = x(i+1)-y
      flpsq = flp*flp
      fl0 = y-x(i)
      fl0sq = fl0*fl0
      i0 = ii(i)
      ip = i0+int1
      if (itab(1) .ne. 1) go to 102
c
c calculate f(y)
c
      a = (w(i0)*flpsq*flp+w(ip)*fl0sq*fl0)*sixth*flki
      b = (f(ip)*flki-w(ip)*flk*sixth)*fl0
      c = (f(i0)*flki-w(i0)*flk*sixth)*flp
      tab(1) = a+b+c
  102 continue
      if (itab(2) .ne. 1) go to 104
c
c calculate first derivative at y
c
      a = (w(ip)*fl0sq-w(i0)*flpsq)*0.5d0*flki
      b = (f(ip)-f(i0))*flki
      c = (w(i0)-w(ip))*flk*sixth
      tab(2) = a+b+c
  104 continue
      if (itab(3) .ne. 1) go to 106
c
c calculate second derivative at y
c
      tab(3) = (w(i0)*flp+w(ip)*fl0)*flki
  106 continue
      return
      end
      subroutine terp1 (n,x,lenf,f,lenw,w,y,int1,tab,itab,isrch)
c      implicit double precision (a-h,o-z)
      save
c
c
c dimension of           x(n),f(int*(n-1)+1),w(int*(n-1)+1),tab(3),
c arguments              itab(3)
c
c latest revision        february 1974
c
c usage                  call terp1 (n,x,f,w,y,int1,tab,itab)
c
c arguments
c
c on input               n
c                          the number of data points.  n must be at
c                          least 4.
c
c                        x
c                          table of n independent variable values in
c                          ascending order.  dimension of x in the
c                          calling program must be at least n.
c
c                        f
c                          table of n corresponding dependent variable
c                          values separated by interval int, usually
c                          unity for one-dimensional interpolation.
c                          dimension of f in the calling program must be
c                          at least (int*(n-1)+1).
c
c                        w
c                          table of second derivatives computed by
c                          coeff1.  the separation of tabular entries is
c                          int.  dimension of w in the calling program
c                          must be at least (int*(n-1)+1).
c
c                        y
c                          value of the independent variable at which
c                          interpolation is required.  if y lies outside
c                          the range of the table, extrapolation takes
c                          place.
c
c                        int1
c                          spacing of tabular entries in f and w arrays.
c                          this is usually unity for one-dimensional
c                          interpolation.
c
c
c                        itab
c                          three element integer array defining
c                          interpolation to be performed at y.
c                            if itab(1) = 1, the function value is
c                                            returned in tab(1).
c                            if itab(2) = 1, the first derivative is
c                                            returned in tab(2).
c                            if itab(3) = 1, the second derivative is
c                                            returned in tab(3).
c                          if itab(i) = 0 for i = 1, 2 or 3, the
c                          corresponding function value or derivative is
c                          not computed and tab(i) is not referenced.
c
c (Added by M O'Brien    isrch = 1 if the interpolation position needs
c  25/1/89)                to be calculated (by search) otherwise it is
c                          assumed it has been calculated in a previous
c                          call
c
c on output              tab
c                          three element array in which interpolated
c                          function value, first and second derivatives
c                          are returned as dictated by itab (see above).
c
c space required         47 (octal) = 39 (decimal)
c
c timing                 this procedure is fast.  the maximum time for
c                        the binary search is proportional to log(n).
c                        the time for function and derivative evaluation
c                        is independent of n.  for 21 data points the
c                        time taken on the ncar 7600 is about
c                        80 microseconds.
c
c
      integer int1,isav,isrch,lenf,lenf2,lenw,lenw2,n
      double precision x(n),f(lenf),w(lenw),tab(3),y
      integer itab(3)
c the following call is for gathering statistics on library use at ncar
c
c perform search
c
      if(isrch.eq.1) call search (y,x,n,isav)
c
c carry out interpolation (or extrapolation)
c
      lenf2=(n-1)*int1+1
      lenw2=(n-1)*int1+1
      call intrp (n,x,lenf2,f,lenw2,w,y,isav,int1,tab,itab)
      return
      end
c
c
      subroutine coeff2 (nx,x,ny,y,f,fxx,fyy,fxxyy,idm,ibd,lenwk,wk)
c      implicit double precision (a-h,o-z)
c
c
c dimension of           x(nx),y(ny),f(idm,ny),fxx(idm,ny),fyy(idm,ny),
c arguments              fxxyy(idm,ny),ibd(4),wk(3*max0(nx,ny)+1)
c                        (idm must be .ge. nx)
c
c latest revision        february 1974
c
c usage                  call coeff2 (nx,x,ny,y,f,fxx,fyy,fxxyy,idm,ibd,
c                                     wk)
c
c arguments
c
c on input               nx
c                          number of grid points in the x-direction.  nx
c                          must be at least 4.
c
c                        x
c                          table of nx values of the first independent
c                          variable arranged in ascending order.
c                          dimension of x in the calling program must be
c                          at least nx.
c
c                        ny
c                          number of grid points in the y-direction.  ny
c                          must be at least 4.
c
c                        y
c                          table of ny values of the second independent
c                          variable arranged in ascending order.
c                          dimension of y in the calling program must be
c                          at least ny.
c
c                        f
c                          two dimensional array of function values at
c                          the grid points defined by the arrays x and
c                          y.  dimension of f in the calling program is
c                          (idm, nyy) where
c                              idm .ge. nx
c                              nyy .ge. ny
c
c                        idm
c                          first dimension in the calling program of
c                          arrays f, fxx, fyy, fxxyy.  idm must be at
c                          least nx.
c
c                        ibd
c                          four element integer array defining boundary
c                          conditions according to the following code.
c                          for ibd(1)
c                          = 1  the second derivative of f with respect
c                               to x is given at (x(1),y(j)) for
c                               j = 1,ny,1.  values of this second
c                               derivative must be placed in fxx(1,j)
c                               for j = 1,ny,1, before calling coeff2.
c                          = 2  the first derivative of f with respect
c                               to x is given at (x(1),y(j)) for
c                               j = 1,ny,1.  values of the derivative
c                               must be placed in fxx(1,j) for
c                               j = 1,ny,1 before calling coeff2.
c                          = 3  periodic boundary condition in the
c                               x-direction.  (x(1),y(j)) and
c                               and (x(nx),y(j)) are equivalent points
c                               for j = 1,ny,1.  f(1,j) and f(nx,j) are
c                               equal.
c                          = 4  the first derivative of f with respect
c                               to x at (x(1),y(j)) is computed by
c                               fitting a cubic to f(1,j) through f(4,j)
c                               for j = 1,ny,1.
c
c                          similarly, ibd(2) defines the boundary
c                          condition at (x(nx),y(j)) for j = 1,ny,1.
c                          when ibd(2) = 1 (or 2) the values of the
c                          second (or first) derivative of f with
c                          respect to x are placed in fxx(nx,j) for
c                          j = 1,ny,1.
c                            note that if ibd(1) = 3, consistency
c                            requires that ibd(2) = 3 also.
c                          for ibd(3)
c                          = 1  the second derivative of f with respect
c                               to y is given at (x(i),y(1)).  place
c                               values of the derivative in fyy(i,1) for
c                               i = 1,nx,1 before calling coeff2.
c                          = 2  the first derivative of f with respect
c                               to y is given at (x(i),y(1)).  values of
c                               this derivative must be placed in
c                               fyy(i,1) for i = 1,nx,1 before calling
c                               coeff2.
c                          = 3  periodic boundary condition in the
c                               y-direction.  (x(i),y(1)) and
c                               (x(i),y(ny)) are equivalent points.
c                               f(i,1) and f(i,ny) are equal.
c                          = 4  the first derivative of f with respect
c                               to y at (x(i),y(1)) is computed by
c                               fitting a cubic to f(i,1) through f(i,4)
c                               for i = 1,nx,1.
c
c                          similary, ibd(4) defines the boundary
c                          condition at (x(i),y(ny)) for i = 1,nx,1 and
c                          given derivative values are placed in
c                          fyy(i,ny).
c                            note that consistency demands that if
c                            ibd(3) = 3, then ibd(4) = 3 also.
c
c                        wk
c                          work area of dimension lenwk at least
c                          (3*max0(nx,ny)+1)
c
c on output              fxx
c                          array of second derivatives of f with respect
c                          to x computed by coeff2.  fxx(i,j) is
c                          derivative at (x(i),y(j)).  as for f,
c                          dimension of fxx in the calling program is
c                          (idm,nyy).
c
c                        fyy
c                          array of second derivatives of f with respect
c                          to y computed by coeff2.  dimension of fyy in
c                          the calling program is (idm,nyy).
c
c                        fxxyy
c                          array of fourth derivatives
c                          (d/dx)**2*(d/dy)**2*f, computed by coeff2.
c                          dimension of fxxyy in the calling program is
c                          (idm,nyy).
c
c                        the arrays x, y, f, fxx, fyy, fxxyy are used as
c                        input for the routine terp2 which performs
c                        interpolation at required values of the two
c                        independent variables.
c
c space required         370 (octal) = 248 (decimal)
c
c timing                 the timing is proportional to nx*ny.  for a
c                        21 x 21 grid of data points, the time taken on
c                        the ncar 7600 was 19.5 milliseconds.
c
c
      integer i,idm,j,lenw,lenwk,lenf,nx,ny
      double precision x(nx),y(ny),f(idm,ny),fxx(idm,ny)
      double precision fyy(idm,ny),fxxyy(idm,ny),wk(lenwk)
      integer ibd(4),iloc(2),jloc(2)
      data iloc(1),iloc(2),jloc(1),jloc(2)/1,1,4,4/
c the following call is for gathering statistics on library use at ncar
c
c compute fxx
c
      lenf=nx
      lenw=nx
      do 101 j=1,ny
	 call coeff1 (nx,x,lenf,f(1,j),lenw,fxx(1,j),ibd(1),1,wk)
  101 continue
c
c compute fyy
c
      lenf=idm*(ny-1)+1
      lenw=idm*(ny-1)+1
      do 102 i=1,nx
	 call coeff1 (ny,y,lenf,f(i,1),lenw,fyy(i,1),ibd(3),idm,wk)
  102 continue
c
c check for periodic boundary condition in both directions
c
      if (ibd(1) .eq. 3) go to 103
      if (ibd(3) .eq. 3) go to 105
c
c calculate fxxyy along left and right boundaries
c
      lenf=idm*(ny-1)+1
      lenw=idm*(ny-1)+1
      call coeff1 (ny,y,lenf,fxx(1,1),lenw,fxxyy(1,1),jloc,idm,wk)
      call coeff1 (ny,y,lenf,fxx(nx,1),lenw,fxxyy(nx,1),jloc,idm,wk)
      go to 106
  103 continue
c
c periodic in x direction . calculate fxxyy along lower and upper
c boundaries.
c
      lenf=nx
      lenw=nx
      call coeff1 (nx,x,lenf,fyy(1,1),lenw,fxxyy(1,1),ibd(1),1,wk)
      call coeff1 (nx,x,lenf,fyy(1,ny),lenw,fxxyy(1,ny),ibd(1),1,wk)
c
c calculate remaining fxxyy
c
      do 104 i=1,nx
	 lenf=(idm-1)*ny+1
	 lenw=(idm-1)*ny+1
	 call coeff1 (ny,y,lenf,fxx(i,1),lenw,fxxyy(i,1),iloc,idm,wk)
  104 continue
      go to 108
  105 continue
c
c periodic in y direction. calculate fxxyy along left and right
c boundaries.
c
      lenf=(idm-1)*ny+1
      lenw=(idm-1)*ny+1
      call coeff1 (ny,y,lenf,fxx(1,1),lenw,fxxyy(1,1),ibd(3),idm,wk)
      call coeff1 (ny,y,lenf,fxx(nx,1),lenw,fxxyy(nx,1),ibd(3),idm,wk)
  106 continue
c
c calculate remaining fxxyy
c
      do 107 j=1,ny
	 lenf=nx
	 lenw=nx
	 call coeff1 (nx,x,lenf,fyy(1,j),lenw,fxxyy(1,j),iloc,1,wk)
  107 continue
  108 continue
      return
      end
      double precision function terp2 (xb,yb,nx,x,ny,y,f,fxx,fyy,fxxyy,
     1                                 idm,ixd,iyd,isrch,isav,jsav)
c      implicit double precision (a-h,o-z)
      save
c
c Modified terp2 by adding isrch:  Binary search of grid for
c                    interpolation point is carried out only if
c                    isrch=1.  Otherwise, it is assumed that values
c                    have been generated by a previous call or
c                    supplied by the calling routine.
c
c
c dimension of           x(nx),y(ny),f(idm,ny),fxx(idm,ny),fyy(idm,ny),
c arguments              fxxyy(idm,ny))
c                        (idm must be .ge. nx)
c
c latest revision        february 1974
c
c usage                  r = terp2 (xb,yb,nx,x,ny,y,f,fxx,fyy,fxxyy,idm,
c                                   ixd,iyd)
c
c arguments
c
c on input               xb, yb
c                          values of the independent variables, x and y,
c                          at which interpolation is required.
c
c                        nx
c                          number of grid points in the x-direction.  nx
c                          must be at least 4.
c
c                        x
c                          table of nx values of the independent
c                          variable, x, arranged in ascending order.
c                          dimension of x in the calling program must be
c                          at least nx.
c
c                        ny
c                          number of grid points in the y-direction.  ny
c                          must be at least 4.
c
c                        y
c                          table of ny values of the independent
c                          variable, y, arranged in ascending order.
c                          dimension of y in the calling program must be
c                          at least ny.
c
c                        f
c                          two-dimensional array of function values at
c                          grid points defined by the arrays x and y.
c                          dimension of f in the calling program is
c                          (idm,nyy), where
c                              idm .ge. nx
c                              nyy .ge. ny
c
c                        fxx
c                          array of second derivatives of f with respect
c                          to x computed by coeff2.  dimension of fxx in
c                          the calling program is (idm,nyy).  see under
c                          f above.
c
c                        fyy
c                          array of second derivatives of f with respect
c                          to y computed by coeff2.  dimension of fyy in
c                          the calling program is (idm,nyy).
c
c                        fxxyy
c                          array of fourth derivatives,
c                          (d/dx)**2*(d/dy)**2*f, computed by coeff2.
c                          dimension of fxxyy in the calling program is
c                          (idm,nyy).
c
c                        idm
c                          first dimension in the calling program of
c                          arrays f, fxx, fyy and fxxyy,
c                              idm .ge. nx
c
c                        ixd, iyd
c                          define derivative to be returned by the
c                          function terp2.  ixd, iyd may each take the
c                          the values 0, 1, 2.  the derivative returned
c                          is (d/dx)**ixd*(d/dy)**iyd*f.
c                            note that if ixd = iyd = 0, the function
c                            value itself is returned.
c
c space required         220 (octal) = 144 (decimal)
c
c timing                 this procedure is fast.  the maximum time for
c                        the binary search is proportional to
c                        log(nx*ny).  the time for function evaluation
c                        is independent of n.  for a 21 x 21 grid of
c                        data points, an average time for an
c                        interpolation on the ncar cdc 7600 is about
c                        .29 milliseconds.
c
c
      integer i,i1,idm,isav,isrch,ixd,iyd,j,j1,jj,jsav,nx,ny
      double precision xb,yb
      double precision x(nx),y(ny),f(idm,ny),fxx(idm,ny)
      double precision fyy(idm,ny),fxxyy(idm,ny),ff(2)
      double precision ww(2),tab(3)
      integer itab(3)
c
c search in x and y arrays.
c
      if(isrch.eq.1)  then
	call search (xb,x,nx,i)
	call search (yb,y,ny,j)
	isav=i
	jsav=j
      endif
c
c interpolate in x direction
c
      itab(1) = 0
      itab(2) = 0
      itab(3) = 0
      i1 = ixd+1
      itab(i1) = 1
      do 102 j1=1,2
	 jj = jsav+j1-1
	 call intrp (nx,x,nx,f(1,jj),nx,fxx(1,jj),xb,isav,1,tab,itab)
	 ff(j1) = tab(i1)
	 call intrp (nx,x,nx,fyy(1,jj),nx,fxxyy(1,jj),
     1               xb,isav,1,tab,itab)
	 ww(j1) = tab(i1)
  102 continue
c
c interpolate in y direction
c
      itab(i1) = 0
      j1 = iyd+1
      itab(j1) = 1
      call intrp (2,y(jsav),2,ff,2,ww,yb,1,1,tab,itab)
      terp2 = tab(j1)
      return
      end
c cubspl     from testlib             version   3          05/19/77
