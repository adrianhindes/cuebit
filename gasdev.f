      FUNCTION gasdev(idum)
      INTEGER idum
      REAL gasdev
c 
c This function uses intrinsic function rand. It returns a 
c normally-distributed set of random numbers with zero mean and 
c unit variance, using rand as the source of uniform deviates.
c
      INTEGER iset
      REAL fac,gset,rsq,v1,v2,ran1
      SAVE iset,gset
      DATA iset/0/
      if (iset.eq.0) then
  1      v1 = 2.0*rand(idum)-1.0
         v2 = 2.0*rand(idum)-1.0
         rsq = v1**2+v2**2
         if(rsq.ge.1.0.or.rsq.eq.0.0)goto 1
         fac = sqrt(-2.0*log(rsq)/rsq)
         gset = v1*fac
         gasdev = v2*fac
         iset = 1
      else
         gasdev = gset
         iset = 0
      endif
      return
      END



