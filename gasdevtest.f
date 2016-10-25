      program gasdevtest
      real gasdev,mean,variance
      real rr(5000)
      integer j,iseed
      external gasdev
      iseed = 56
      call srand(iseed)
      idum = 0
      mean = 0.0
      variance = 0.0
      do j = 1,5000
         rr(j) = gasdev(idum)
         mean = mean+rr(j)
      enddo
      mean = mean/5000.0
      variance = 0.0
      do j = 1,5000
         variance = variance+(rr(j)-mean)**2
      enddo
      variance = variance/5000.0
      write(12,fmt='()')
      write(6,fmt='("mean     =",f6.2)')mean
      write(12,fmt='()')
      write(6,fmt='("variance =",f6.2)')variance
      END
