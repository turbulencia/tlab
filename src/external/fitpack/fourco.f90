      subroutine fourco(t,n,c,alfa,m,ress,resc,wrk1,wrk2,ier)
!  subroutine fourco calculates the integrals
!            /t(n-3)
!ress(i) =      !        s(x)*sin(alfa(i)*x) dx    and
!      t(4)/
!            /t(n-3)
!resc(i) =      !        s(x)*cos(alfa(i)*x) dx, i=1,...,m,
!      t(4)/
!  where s(x) denotes a cubic spline which is given in its
!  b-spline representation.
!
!  calling sequence:
! call fourco(t,n,c,alfa,m,ress,resc,wrk1,wrk2,ier)
!
!  input parameters:
!t    : real array,length n, containing the knots of s(x).
!n    : integer, containing the total number of knots. n>=10.
!c    : real array,length n, containing the b-spline coefficients.
!alfa : real array,length m, containing the parameters alfa(i).
!m    : integer, specifying the number of integrals to be computed.
!wrk1 : real array,length n. used as working space
!wrk2 : real array,length n. used as working space
!
!  output parameters:
!ress : real array,length m, containing the integrals ress(i).
!resc : real array,length m, containing the integrals resc(i).
!ier  : error flag:
!  ier=0 : normal return.
!  ier=10: invalid input data (see restrictions).
!
!  restrictions:
!n >= 10
!t(4) < t(5) < ... < t(n-4) < t(n-3).
!t(1) <= t(2) <= t(3) <= t(4).
!t(n-3) <= t(n-2) <= t(n-1) <= t(n).
!
!  other subroutines required: fpbfou,fpcsin
!
!  references :
!dierckx p. : calculation of fouriercoefficients of discrete
!         functions using cubic splines. j. computational
!         and applied mathematics 3 (1977) 207-209.
!dierckx p. : curve and surface fitting with splines, monographs on
!         numerical analysis, oxford university press, 1993.
!
!  author :
!p.dierckx
!dept. computer science, k.u.leuven
!celestijnenlaan 200a, b-3001 heverlee, belgium.
!e-mail : Paul.Dierckx@cs.kuleuven.ac.be
!
!  latest update : march 1987
!
!  ..scalar arguments..
      integer n,m,ier
!  ..array arguments..
      real t(n),c(n),wrk1(n),wrk2(n),alfa(m),ress(m),resc(m)
!  ..local scalars..
      integer i,j,n4
      real rs,rc
!  ..
      n4 = n-4
!  before starting computations a data check is made. in the input data
!  are invalid, control is immediately repassed to the calling program.
      ier = 10
      if(n.lt.10) go to 50
      j = n
      do 10 i=1,3
        if(t(i).gt.t(i+1)) go to 50
        if(t(j).lt.t(j-1)) go to 50
        j = j-1
  10  continue
      do 20 i=4,n4
        if(t(i).ge.t(i+1)) go to 50
  20  continue
      ier = 0
!  main loop for the different alfa(i).
      do 40 i=1,m
!  calculate the integrals
!wrk1(j) = integral(nj,4(x)*sin(alfa*x))    and
!wrk2(j) = integral(nj,4(x)*cos(alfa*x)),  j=1,2,...,n-4,
!  where nj,4(x) denotes the normalised cubic b-spline defined on the
!  knots t(j),t(j+1),...,t(j+4).
         call fpbfou(t,n,alfa(i),wrk1,wrk2)
!  calculate the integrals ress(i) and resc(i).
         rs = 0.
         rc = 0.
         do 30 j=1,n4
            rs = rs+c(j)*wrk1(j)
            rc = rc+c(j)*wrk2(j)
  30     continue
         ress(i) = rs
         resc(i) = rc
  40  continue
  50  return
      end
