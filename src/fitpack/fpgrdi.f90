      subroutine fpgrdi(ifsu,ifsv,ifbu,ifbv,iback,u,mu,v,mv,z,mz,dz,&
       iop0,iop1,tu,nu,tv,nv,p,c,nc,sq,fp,fpu,fpv,mm,mvnu,spu,spv,&
       right,q,au,av1,av2,bu,bv,aa,bb,cc,cosi,nru,nrv)
!  ..

      IMPLICIT NONE

#include "types.h"

!  ..scalar arguments..
      TREAL p,sq,fp
      TINTEGER ifsu,ifsv,ifbu,ifbv,iback,mu,mv,mz,iop0,iop1,nu,&
           nv,nc,mm,mvnu
!  ..array arguments..
      TREAL u(mu),v(mv),z(mz),dz(3),tu(nu),tv(nv),c(nc),fpu(nu),&
           fpv(nv),spu(mu,4),spv(mv,4),right(mm),q(mvnu),au(nu,5),&
           av1(nv,6),av2(nv,4),aa(2,mv),bb(2,nv),cc(nv),cosi(2,nv),&
           bu(nu,5),bv(nv,5)
      TINTEGER nru(mu),nrv(mv)
!  ..local scalars..
      TREAL arg,co,dz1,dz2,dz3,fac,fac0,pinv,piv,si,term,one,three,&
           half
      TINTEGER i,ic,ii,ij,ik,iq,irot,it,iz,i0,i1,i2,i3,j,jj,jk,jper,&
       j0,j1,k,k1,k2,l,l0,l1,l2,mvv,ncof,nrold,nroldu,nroldv,number,&
       numu,numu1,numv,numv1,nuu,nu4,nu7,nu8,nu9,nv11,nv4,nv7,nv8,n1
!  ..local arrays..
      TREAL h(5),h1(5),h2(4)
!  ..function references..
      TINTEGER min0
!  ..subroutine references..
!fpback,fpbspl,fpgivs,fpcyt1,fpcyt2,fpdisc,fpbacp,fprota
!  ..
!  let
!       |   (spu)    |            |   (spv)    |
!(au) = | ---------- |     (av) = | ---------- |
!       | (1/p) (bu) |            | (1/p) (bv) |
!
!                        | z  1 0 |
!                    q = | ------ |
!                        | 0  1 0 |
!
!  with c      : the (nu-4) x (nv-4) matrix which contains the b-spline
!        coefficients.
!   z      : the mu x mv matrix which contains the function values.
!   spu,spv: the mu x (nu-4), resp. mv x (nv-4) observation matrices
!        according to the least-squares problems in the u-,resp.
!        v-direction.
!   bu,bv  : the (nu-7) x (nu-4),resp. (nv-7) x (nv-4) matrices
!        containing the discontinuity jumps of the derivatives
!        of the b-splines in the u-,resp.v-variable at the knots
!  the b-spline coefficients of the smoothing spline are then calculated
!  as the least-squares solution of the following over-determined linear
!  system of equations
!
!(1)  (av) c (au)1 = q
!
!  subject to the constraints
!
!(2)  c(i,nv-3+j) = c(i,j), j=1,2,3 ; i=1,2,...,nu-4
!
!(3)  if iop0 = 0  c(1,j) = dz(1)
!    iop0 = 1  c(1,j) = dz(1)
!              c(2,j) = dz(1)+(dz(2)*cosi(1,j)+dz(3)*cosi(2,j))*
!                       tu(5)/3. = cc(j) , j=1,2,...nv-4
!
!(4)  if iop1 = 1  c(nu-4,j) = 0, j=1,2,...,nv-4.
!
!  set constants
      one = C_1_R
      three = C_3_R
      half = C_05_R
!  initialization
      nu4 = nu-4
      nu7 = nu-7
      nu8 = nu-8
      nu9 = nu-9
      nv4 = nv-4
      nv7 = nv-7
      nv8 = nv-8
      nv11 = nv-11
      nuu = nu4-iop0-iop1-1
      if(p.gt.C_0_R) pinv = one/p
!  it depends on the value of the flags ifsu,ifsv,ifbu,ifbv and iop0 and
!  on the value of p whether the matrices (spu), (spv), (bu), (bv) and
!  (cosi) still must be determined.
      if(ifsu.ne.0) go to 30
!  calculate the non-zero elements of the matrix (spu) which is the ob-
!  servation matrix according to the least-squares spline approximation
!  problem in the u-direction.
      l = 4
      l1 = 5
      number = 0
      do 25 it=1,mu
        arg = u(it)
  10    if(arg.lt.tu(l1) .or. l.eq.nu4) go to 15
        l = l1
        l1 = l+1
        number = number+1
        go to 10
  15    call fpbspl(tu,nu,3,arg,l,h)
        do 20 i=1,4
          spu(it,i) = h(i)
  20    continue
        nru(it) = number
  25  continue
      ifsu = 1
!  calculate the non-zero elements of the matrix (spv) which is the ob-
!  servation matrix according to the least-squares spline approximation
!  problem in the v-direction.
  30  if(ifsv.ne.0) go to 85
      l = 4
      l1 = 5
      number = 0
      do 50 it=1,mv
        arg = v(it)
  35    if(arg.lt.tv(l1) .or. l.eq.nv4) go to 40
        l = l1
        l1 = l+1
        number = number+1
        go to 35
  40    call fpbspl(tv,nv,3,arg,l,h)
        do 45 i=1,4
          spv(it,i) = h(i)
  45    continue
        nrv(it) = number
  50  continue
      ifsv = 1
      if(iop0.eq.0) go to 85
!  calculate the coefficients of the interpolating splines for cos(v)
!  and sin(v).
      do 55 i=1,nv4
         cosi(1,i) = C_0_R
         cosi(2,i) = C_0_R
  55  continue
      if(nv7.lt.4) go to 85
      do 65 i=1,nv7
         l = i+3
         arg = tv(l)
         call fpbspl(tv,nv,3,arg,l,h)
         do 60 j=1,3
            av1(i,j) = h(j)
  60     continue
         cosi(1,i) = cos(arg)
         cosi(2,i) = sin(arg)
  65  continue
      call fpcyt1(av1,nv7,nv)
      do 80 j=1,2
         do 70 i=1,nv7
            right(i) = cosi(j,i)
  70     continue
         call fpcyt2(av1,nv7,right,right,nv)
         do 75 i=1,nv7
            cosi(j,i+1) = right(i)
  75     continue
         cosi(j,1) = cosi(j,nv7+1)
         cosi(j,nv7+2) = cosi(j,2)
         cosi(j,nv4) = cosi(j,3)
  80  continue
  85  if(p.le.C_0_R) go to  150
!  calculate the non-zero elements of the matrix (bu).
      if(ifbu.ne.0 .or. nu8.eq.0) go to 90
      call fpdisc(tu,nu,5,bu,nu)
      ifbu = 1
!  calculate the non-zero elements of the matrix (bv).
  90  if(ifbv.ne.0 .or. nv8.eq.0) go to 150
      call fpdisc(tv,nv,5,bv,nv)
      ifbv = 1
!  substituting (2),(3) and (4) into (1), we obtain the overdetermined
!  system
! (5)  (avv) (cr) (auu)1 = (qq)
!  from which the nuu*nv7 remaining coefficients
! c(i,j) , i=2+iop0,3+iop0,...,nu-4-iop1 ; j=1,2,...,nv-7 ,
!  the elements of (cr), are then determined in the least-squares sense.
!  simultaneously, we compute the resulting sum of squared residuals sq.
 150  dz1 = dz(1)
      do 155 i=1,mv
         aa(1,i) = dz1
 155  continue
      if(nv8.eq.0 .or. p.le.C_0_R) go to 165
      do 160 i=1,nv8
         bb(1,i) = C_0_R
 160  continue
 165  mvv = mv
      if(iop0.eq.0) go to 220
      fac = tu(5)/three
      dz2 = dz(2)*fac
      dz3 = dz(3)*fac
      do 170 i=1,nv4
         cc(i) = dz1+dz2*cosi(1,i)+dz3*cosi(2,i)
 170  continue
      do 190 i=1,mv
         number = nrv(i)
         fac = C_0_R
         do 180 j=1,4
            number = number+1
            fac = fac+cc(number)*spv(i,j)
 180     continue
         aa(2,i) = fac
 190  continue
      if(nv8.eq.0 .or. p.le.C_0_R) go to 220
      do 210 i=1,nv8
         number = i
         fac = C_0_R
         do 200 j=1,5
            fac = fac+cc(number)*bv(i,j)
            number = number+1
 200     continue
         bb(2,i) = fac*pinv
 210  continue
      mvv = mvv+nv8
!  we first determine the matrices (auu) and (qq). then we reduce the
!  matrix (auu) to upper triangular form (ru) using givens rotations.
!  we apply the same transformations to the rows of matrix qq to obtain
!  the (mv+nv8) x nuu matrix g.
!  we store matrix (ru) into au and g into q.
 220  l = mvv*nuu
!  initialization.
      sq = C_0_R
      do 230 i=1,l
        q(i) = C_0_R
 230  continue
      do 240 i=1,nuu
        do 240 j=1,5
          au(i,j) = C_0_R
 240  continue
      l = 0
      nrold = 0
      n1 = nrold+1
      do 420 it=1,mu
        number = nru(it)
!  find the appropriate column of q.
 250    do 260 j=1,mvv
           right(j) = C_0_R
 260    continue
        if(nrold.eq.number) go to 280
        if(p.le.C_0_R) go to 410
!  fetch a new row of matrix (bu).
        do 270 j=1,5
          h(j) = bu(n1,j)*pinv
 270    continue
        i0 = 1
        i1 = 5
        go to 310
!  fetch a new row of matrix (spu).
 280    do 290 j=1,4
          h(j) = spu(it,j)
 290    continue
!  find the appropriate column of q.
        do 300 j=1,mv
          l = l+1
          right(j) = z(l)
 300    continue
        i0 = 1
        i1 = 4
 310    if(nu7-number .eq. iop1) i1 = i1-1
        j0 = n1
!  take into account that we eliminate the constraints (3)
 320     if(j0-1.gt.iop0) go to 360
         fac0 = h(i0)
         do 330 j=1,mv
            right(j) = right(j)-fac0*aa(j0,j)
 330     continue
         if(mv.eq.mvv) go to 350
         j = mv
         do 340 jj=1,nv8
            j = j+1
            right(j) = right(j)-fac0*bb(j0,jj)
 340     continue
 350     j0 = j0+1
         i0 = i0+1
         go to 320
 360     irot = nrold-iop0-1
         if(irot.lt.0) irot = 0
!  rotate the new row of matrix (auu) into triangle.
        do 390 i=i0,i1
          irot = irot+1
          piv = h(i)
          if(piv.eq.C_0_R) go to 390
!  calculate the parameters of the givens transformation.
          call fpgivs(piv,au(irot,1),co,si)
!  apply that transformation to the rows of matrix (qq).
          iq = (irot-1)*mvv
          do 370 j=1,mvv
            iq = iq+1
            call fprota(co,si,right(j),q(iq))
 370      continue
!  apply that transformation to the columns of (auu).
          if(i.eq.i1) go to 390
          i2 = 1
          i3 = i+1
          do 380 j=i3,i1
            i2 = i2+1
            call fprota(co,si,h(j),au(irot,i2))
 380      continue
 390    continue
! we update the sum of squared residuals
        do 395 j=1,mvv
          sq = sq+right(j)**2
 395    continue
 400    if(nrold.eq.number) go to 420
 410    nrold = n1
        n1 = n1+1
        go to 250
 420  continue
!  we determine the matrix (avv) and then we reduce her to
!  upper triangular form (rv) using givens rotations.
!  we apply the same transformations to the columns of matrix
!  g to obtain the (nv-7) x (nu-5-iop0-iop1) matrix h.
!  we store matrix (rv) into av1 and av2, h into c.
!  the nv7 x nv7 upper triangular matrix (rv) has the form
!      | av1 1     |
!   (rv) = |     1 av2 |
!      |  0  1     |
!  with (av2) a nv7 x 4 matrix and (av1) a nv11 x nv11 upper
!  triangular matrix of bandwidth 5.
      ncof = nuu*nv7
!  initialization.
      do 430 i=1,ncof
        c(i) = C_0_R
 430  continue
      do 440 i=1,nv4
        av1(i,5) = C_0_R
        do 440 j=1,4
          av1(i,j) = C_0_R
          av2(i,j) = C_0_R
 440  continue
      jper = 0
      nrold = 0
      do 770 it=1,mv
        number = nrv(it)
 450    if(nrold.eq.number) go to 480
        if(p.le.C_0_R) go to 760
!  fetch a new row of matrix (bv).
        n1 = nrold+1
        do 460 j=1,5
          h(j) = bv(n1,j)*pinv
 460    continue
!  find the appropiate row of g.
        do 465 j=1,nuu
          right(j) = C_0_R
 465    continue
        if(mv.eq.mvv) go to 510
        l = mv+n1
        do 470 j=1,nuu
          right(j) = q(l)
          l = l+mvv
 470    continue
        go to 510
!  fetch a new row of matrix (spv)
 480    h(5) = C_0_R
        do 490 j=1,4
          h(j) = spv(it,j)
 490    continue
!  find the appropiate row of g.
        l = it
        do 500 j=1,nuu
          right(j) = q(l)
          l = l+mvv
 500    continue
!  test whether there are non-zero values in the new row of (avv)
!  corresponding to the b-splines n(j,v),j=nv7+1,...,nv4.
 510     if(nrold.lt.nv11) go to 710
         if(jper.ne.0) go to 550
!  initialize the matrix (av2).
         jk = nv11+1
         do 540 i=1,4
            ik = jk
            do 520 j=1,5
               if(ik.le.0) go to 530
               av2(ik,i) = av1(ik,j)
               ik = ik-1
 520        continue
 530        jk = jk+1
 540     continue
         jper = 1
!  if one of the non-zero elements of the new row corresponds to one of
!  the b-splines n(j;v),j=nv7+1,...,nv4, we take account of condition
!  (2) for setting up this row of (avv). the row is stored in h1( the
!  part with respect to av1) and h2 (the part with respect to av2).
 550     do 560 i=1,4
            h1(i) = C_0_R
            h2(i) = C_0_R
 560     continue
         h1(5) = C_0_R
         j = nrold-nv11
         do 600 i=1,5
            j = j+1
            l0 = j
 570        l1 = l0-4
            if(l1.le.0) go to 590
            if(l1.le.nv11) go to 580
            l0 = l1-nv11
            go to 570
 580        h1(l1) = h(i)
            go to 600
 590        h2(l0) = h2(l0) + h(i)
 600     continue
!  rotate the new row of (avv) into triangle.
         if(nv11.le.0) go to 670
!  rotations with the rows 1,2,...,nv11 of (avv).
         do 660 j=1,nv11
            piv = h1(1)
            i2 = min0(nv11-j,4)
            if(piv.eq.C_0_R) go to 640
!  calculate the parameters of the givens transformation.
            call fpgivs(piv,av1(j,1),co,si)
!  apply that transformation to the columns of matrix g.
            ic = j
            do 610 i=1,nuu
               call fprota(co,si,right(i),c(ic))
               ic = ic+nv7
 610        continue
!  apply that transformation to the rows of (avv) with respect to av2.
            do 620 i=1,4
               call fprota(co,si,h2(i),av2(j,i))
 620        continue
!  apply that transformation to the rows of (avv) with respect to av1.
            if(i2.eq.0) go to 670
            do 630 i=1,i2
               i1 = i+1
               call fprota(co,si,h1(i1),av1(j,i1))
 630        continue
 640        do 650 i=1,i2
               h1(i) = h1(i+1)
 650        continue
            h1(i2+1) = C_0_R
 660     continue
!  rotations with the rows nv11+1,...,nv7 of avv.
 670     do 700 j=1,4
            ij = nv11+j
            if(ij.le.0) go to 700
            piv = h2(j)
            if(piv.eq.C_0_R) go to 700
!  calculate the parameters of the givens transformation.
            call fpgivs(piv,av2(ij,j),co,si)
!  apply that transformation to the columns of matrix g.
            ic = ij
            do 680 i=1,nuu
               call fprota(co,si,right(i),c(ic))
               ic = ic+nv7
 680        continue
            if(j.eq.4) go to 700
!  apply that transformation to the rows of (avv) with respect to av2.
            j1 = j+1
            do 690 i=j1,4
               call fprota(co,si,h2(i),av2(ij,i))
 690        continue
 700     continue
! we update the sum of squared residuals
         do 705 i=1,nuu
           sq = sq+right(i)**2
 705     continue
         go to 750
!  rotation into triangle of the new row of (avv), in case the elements
!  corresponding to the b-splines n(j;v),j=nv7+1,...,nv4 are all zero.
 710     irot =nrold
         do 740 i=1,5
            irot = irot+1
            piv = h(i)
            if(piv.eq.C_0_R) go to 740
!  calculate the parameters of the givens transformation.
            call fpgivs(piv,av1(irot,1),co,si)
!  apply that transformation to the columns of matrix g.
            ic = irot
            do 720 j=1,nuu
               call fprota(co,si,right(j),c(ic))
               ic = ic+nv7
 720        continue
!  apply that transformation to the rows of (avv).
            if(i.eq.5) go to 740
            i2 = 1
            i3 = i+1
            do 730 j=i3,5
               i2 = i2+1
               call fprota(co,si,h(j),av1(irot,i2))
 730        continue
 740     continue
! we update the sum of squared residuals
         do 745 i=1,nuu
           sq = sq+right(i)**2
 745     continue
 750     if(nrold.eq.number) go to 770
 760     nrold = nrold+1
         go to 450
 770  continue
!  test whether the b-spline coefficients must be determined.
      if(iback.ne.0) return
!  backward substitution to obtain the b-spline coefficients as the
!  solution of the linear system    (rv) (cr) (ru)1 = h.
!  first step: solve the system  (rv) (c1) = h.
      k = 1
      do 780 i=1,nuu
         call fpbacp(av1,av2,c(k),nv7,4,c(k),5,nv)
         k = k+nv7
 780  continue
!  second step: solve the system  (cr) (ru)1 = (c1).
      k = 0
      do 800 j=1,nv7
        k = k+1
        l = k
        do 790 i=1,nuu
          right(i) = c(l)
          l = l+nv7
 790    continue
        call fpback(au,right,nuu,5,right,nu)
        l = k
        do 795 i=1,nuu
           c(l) = right(i)
           l = l+nv7
 795    continue
 800  continue
!  calculate from the conditions (2)-(3)-(4), the remaining b-spline
!  coefficients.
      ncof = nu4*nv4
      i = nv4
      j = 0
      do 805 l=1,nv4
         q(l) = dz1
 805  continue
      if(iop0.eq.0) go to 815
      do 810 l=1,nv4
         i = i+1
         q(i) = cc(l)
 810  continue
 815  if(nuu.eq.0) go to 850
      do 840 l=1,nuu
         ii = i
         do 820 k=1,nv7
            i = i+1
            j = j+1
            q(i) = c(j)
 820     continue
         do 830 k=1,3
            ii = ii+1
            i = i+1
            q(i) = q(ii)
 830     continue
 840  continue
 850  if(iop1.eq.0) go to 870
      do 860 l=1,nv4
         i = i+1
         q(i) = C_0_R
 860  continue
 870  do 880 i=1,ncof
         c(i) = q(i)
 880  continue
!  calculate the quantities
!res(i,j) = (z(i,j) - s(u(i),v(j)))**2 , i=1,2,..,mu;j=1,2,..,mv
!fp = sumi=1,mu(sumj=1,mv(res(i,j)))
!fpu(r) = sum11i(sumj=1,mv(res(i,j))) , r=1,2,...,nu-7
!          tu(r+3) <= u(i) <= tu(r+4)
!fpv(r) = sumi=1,mu(sum11j(res(i,j))) , r=1,2,...,nv-7
!          tv(r+3) <= v(j) <= tv(r+4)
      fp = C_0_R
      do 890 i=1,nu
        fpu(i) = C_0_R
 890  continue
      do 900 i=1,nv
        fpv(i) = C_0_R
 900  continue
      iz = 0
      nroldu = 0
!  main loop for the different grid points.
      do 950 i1=1,mu
        numu = nru(i1)
        numu1 = numu+1
        nroldv = 0
        do 940 i2=1,mv
          numv = nrv(i2)
          numv1 = numv+1
          iz = iz+1
!  evaluate s(u,v) at the current grid point by making the sum of the
!  cross products of the non-zero b-splines at (u,v), multiplied with
!  the appropiate b-spline coefficients.
          term = C_0_R
          k1 = numu*nv4+numv
          do 920 l1=1,4
            k2 = k1
            fac = spu(i1,l1)
            do 910 l2=1,4
              k2 = k2+1
              term = term+fac*spv(i2,l2)*c(k2)
 910        continue
            k1 = k1+nv4
 920      continue
!  calculate the squared residual at the current grid point.
          term = (z(iz)-term)**2
!  adjust the different parameters.
          fp = fp+term
          fpu(numu1) = fpu(numu1)+term
          fpv(numv1) = fpv(numv1)+term
          fac = term*half
          if(numv.eq.nroldv) go to 930
          fpv(numv1) = fpv(numv1)-fac
          fpv(numv) = fpv(numv)+fac
 930      nroldv = numv
          if(numu.eq.nroldu) go to 940
          fpu(numu1) = fpu(numu1)-fac
          fpu(numu) = fpu(numu)+fac
 940    continue
        nroldu = numu
 950  continue
      return
      end
