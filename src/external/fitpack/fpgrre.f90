      subroutine fpgrre(ifsx,ifsy,ifbx,ifby,x,mx,y,my,z,mz,kx,ky,tx,nx,&
       ty,ny,p,c,nc,fp,fpx,fpy,mm,mynx,kx1,kx2,ky1,ky2,spx,spy,right,q,&
       ax,ay,bx,by,nrx,nry)
!  ..

      IMPLICIT NONE

#include "types.h"

!  ..scalar arguments..
      TREAL p,fp
      TINTEGER ifsx,ifsy,ifbx,ifby,mx,my,mz,kx,ky,nx,ny,nc,mm,mynx,&
       kx1,kx2,ky1,ky2
!  ..array arguments..
      TREAL x(mx),y(my),z(mx,my),tx(nx),ty(ny),c(nc),spx(mx,kx1),&
           spy(my,ky1), right(mm),q(mynx),ax(nx,kx2),bx(nx,kx2),&
           ay(ny,ky2),by(ny,ky2), fpx(nx),fpy(ny)
      TINTEGER nrx(mx),nry(my)
!  ..local scalars..
      TREAL arg,cos,fac,pinv,piv,sin,term,one,half
      TINTEGER i,ibandx,ibandy,ic,iq,irot,it,i1,i2,i3,j,k,k1,k2,l,&
       l1,l2,ncof,nk1x,nk1y,nrold,nroldx,nroldy,number,numx,numx1,&
       numy,numy1,n1
!  ..local arrays..
      TINTEGER li, lj
      TREAL h(7), hh(6), r1, f
      TREAL stor1,stor2
!  ..subroutine references..
!fpback,fpbspl,fpgivs,fpdisc,fprota
!  ..
!  the b-spline coefficients of the smoothing spline are calculated as
!  the least-squares solution of the over-determined linear system of
!  equations  (ay) c (ax) = q       where
!
!       |   (spx)    |            |   (spy)    |
!(ax) = | ---------- |     (ay) = | ---------- |
!       | (1/p) (bx) |            | (1/p) (by) |
!
!                        | z   0 |
!                    q = | ------ |
!                        | 0   0 |
!
!  with c      : the (ny-ky-1) x (nx-kx-1) matrix which contains the
!        b-spline coefficients.
!   z      : the my x mx matrix which contains the function values.
!   spx,spy: the mx x (nx-kx-1) and  my x (ny-ky-1) observation
!        matrices according to the least-squares problems in
!        the x- and y-direction.
!   bx,by  : the (nx-2*kx-1) x (nx-kx-1) and (ny-2*ky-1) x (ny-ky-1)
!        matrices which contain the discontinuity jumps of the
!        derivatives of the b-splines in the x- and y-direction.
      one = C_1_R
      r1 = C_1_R
      half = C_05_R
      nk1x = nx-kx1
      nk1y = ny-ky1
      if(p.gt.0.) pinv = one/p
!  it depends on the value of the flags ifsx,ifsy,ifbx and ifby and on
!  the value of p whether the matrices (spx),(spy),(bx) and (by) still
!  must be determined.
      if(ifsx.ne.0) go to 50
!  calculate the non-zero elements of the matrix (spx) which is the
!  observation matrix according to the least-squares spline approximat-
!  ion problem in the x-direction.
      l = kx1
      l1 = kx2
      number = 0
      do 40 it=1,mx
        arg = x(it)
!  10    if(arg.lt.tx(l1) .or. l.eq.nk1x) go to 20
!l = l1
!l1 = l+1
!number = number+1
!go to 10
        
        DO l=l,nk1x-1
           IF ( arg .LT. tx(l+1) ) GOTO 20
           number = number + 1
        ENDDO

!  20    call fpbspl(tx,nx,kx,arg,l,h)

 20     h(1) = r1
!CDIR LOOPCNT=6
        DO j=1,kx
!CDIR LOOPCNT=6
           DO i=1,j
              hh(i) = h(i)
           ENDDO
           h(1) = C_0_R
!CDIR LOOPCNT=6
           DO i=1,j
              li = l+i
              lj = li-j
              f = hh(i)/(tx(li)-tx(lj))
              h(i) = h(i)+f*(tx(li)-arg)
              h(i+1) = f*(arg-tx(lj))
           ENDDO
        ENDDO

!CDIR LOOPCNT=6
        do 30 i=1,kx1
          spx(it,i) = h(i)
  30    continue
        nrx(it) = number
  40  continue
      ifsx = 1
  50  if(ifsy.ne.0) go to 100
!  calculate the non-zero elements of the matrix (spy) which is the
!  observation matrix according to the least-squares spline approximat-
!  ion problem in the y-direction.
      l = ky1
      l1 = ky2
      number = 0
      do 90 it=1,my
        arg = y(it)
!  60    if(arg.lt.ty(l1) .or. l.eq.nk1y) go to 70
!l = l1
!l1 = l+1
!number = number+1
!go to 60

        DO l=l,nk1y-1
           IF ( arg .LT. ty(l+1) ) GOTO 70
           number = number + 1
        ENDDO

!  70    call fpbspl(ty,ny,ky,arg,l,h)
        
 70     h(1) = r1
!CDIR LOOPCNT=6
        DO j=1,ky
!CDIR LOOPCNT=6
           DO i=1,j
              hh(i) = h(i)
           ENDDO
           h(1) = C_0_R
!CDIR LOOPCNT=6
           DO i=1,j
              li = l+i
              lj = li-j
              f = hh(i)/(ty(li)-ty(lj))
              h(i) = h(i)+f*(ty(li)-arg)
              h(i+1) = f*(arg-ty(lj))
           ENDDO
        ENDDO        

!CDIR LOOPCNT=6
        do 80 i=1,ky1
          spy(it,i) = h(i)
  80    continue
        nry(it) = number
  90  continue
      ifsy = 1
 100  if(p.le.0.) go to 120
!  calculate the non-zero elements of the matrix (bx).
      if(ifbx.ne.0 .or. nx.eq.2*kx1) go to 110
      call fpdisc(tx,nx,kx2,bx,nx)
      ifbx = 1
!  calculate the non-zero elements of the matrix (by).
 110  if(ifby.ne.0 .or. ny.eq.2*ky1) go to 120
      call fpdisc(ty,ny,ky2,by,ny)
      ifby = 1
!  reduce the matrix (ax) to upper triangular form (rx) using givens
!  rotations. apply the same transformations to the rows of matrix q
!  to obtain the my x (nx-kx-1) matrix g.
!  store matrix (rx) into (ax) and g into q.
 120  l = my*nk1x
!  initialization.
      do 130 i=1,l
        q(i) = 0.
 130  continue
      do 140 i=1,nk1x
!CDIR LOOPCNT=7
        do 140 j=1,kx2
          ax(i,j) = 0.
 140  continue
      nrold = 0
!  ibandx denotes the bandwidth of the matrices (ax) and (rx).
      ibandx = kx1
      do 270 it=1,mx
        number = nrx(it)
 150    if(nrold.eq.number) go to 180
        if(p.le.0.) go to 260
        ibandx = kx2
!  fetch a new row of matrix (bx).
        n1 = nrold+1
!CDIR LOOPCNT=7
        do 160 j=1,kx2
          h(j) = bx(n1,j)*pinv
 160    continue
!  find the appropriate column of q.
        do 170 j=1,my
          right(j) = 0.
 170    continue
        irot = nrold
        go to 210
!  fetch a new row of matrix (spx).
 180    h(ibandx) = 0.
!CDIR LOOPCNT=6
        do 190 j=1,kx1
          h(j) = spx(it,j)
 190    continue
!  find the appropriate column of q.
        do 200 j=1,my
          right(j) = z(it,j)
 200    continue
        irot = number
!  rotate the new row of matrix (ax) into triangle.
 210    do 240 i=1,ibandx
          irot = irot+1
          piv = h(i)
          if(piv.eq.0.) go to 240
!  calculate the parameters of the givens transformation.
          call fpgivs(piv,ax(irot,1),cos,sin)
!  apply that transformation to the rows of matrix q.
          iq = (irot-1)*my
          do 220 j=1,my
            iq = iq+1
!    call fprota(cos,sin,right(j),q(iq))
            stor1 = right(j)
            stor2 = q(iq)
            q(iq) = cos*stor2+sin*stor1
            right(j) = cos*stor1-sin*stor2
 220      continue
!  apply that transformation to the columns of (ax).
          if(i.eq.ibandx) go to 250
          i2 = 1
          i3 = i+1
          do 230 j=i3,ibandx
            i2 = i2+1
!    call fprota(cos,sin,h(j),ax(irot,i2))
            stor1 = h(j)
            stor2 = ax(irot,i2)
            ax(irot,i2) = cos*stor2+sin*stor1
            h(j) = cos*stor1-sin*stor2
 230      continue
 240    continue
 250    if(nrold.eq.number) go to 270
 260    nrold = nrold+1
        go to 150
 270  continue
!  reduce the matrix (ay) to upper triangular form (ry) using givens
!  rotations. apply the same transformations to the columns of matrix g
!  to obtain the (ny-ky-1) x (nx-kx-1) matrix h.
!  store matrix (ry) into (ay) and h into c.
      ncof = nk1x*nk1y
!  initialization.
      do 280 i=1,ncof
        c(i) = 0.
 280  continue
      do 290 i=1,nk1y
!CDIR LOOPCNT=7
        do 290 j=1,ky2
          ay(i,j) = 0.
 290  continue
      nrold = 0
!  ibandy denotes the bandwidth of the matrices (ay) and (ry).
      ibandy = ky1
      do 420 it=1,my
        number = nry(it)
 300    if(nrold.eq.number) go to 330
        if(p.le.0.) go to 410
        ibandy = ky2
!  fetch a new row of matrix (by).
        n1 = nrold+1
!CDIR LOOPCNT=7
        do 310 j=1,ky2
          h(j) = by(n1,j)*pinv
 310    continue
!  find the appropiate row of g.
        do 320 j=1,nk1x
          right(j) = 0.
 320    continue
        irot = nrold
        go to 360
!  fetch a new row of matrix (spy)
 330    h(ibandy) = 0.
!CDIR LOOPCNT=6
        do 340 j=1,ky1
          h(j) = spy(it,j)
 340    continue
!  find the appropiate row of g.
        l = it
        do 350 j=1,nk1x
          right(j) = q(l)
          l = l+my
 350    continue
        irot = number
!  rotate the new row of matrix (ay) into triangle.
 360    do 390 i=1,ibandy
          irot = irot+1
          piv = h(i)
          if(piv.eq.0.) go to 390
!  calculate the parameters of the givens transformation.
          call fpgivs(piv,ay(irot,1),cos,sin)
!  apply that transformation to the colums of matrix g.
          ic = irot
          do 370 j=1,nk1x
!    call fprota(cos,sin,right(j),c(ic))
             stor1 = right(j)
             stor2 = c(ic)
             c(ic) = cos*stor2+sin*stor1
             right(j) = cos*stor1-sin*stor2
            ic = ic+nk1y
 370      continue
!  apply that transformation to the columns of matrix (ay).
          if(i.eq.ibandy) go to 400
          i2 = 1
          i3 = i+1
          do 380 j=i3,ibandy
            i2 = i2+1
!    call fprota(cos,sin,h(j),ay(irot,i2))
            stor1 = h(j)
            stor2 = ay(irot,i2)
            ay(irot,i2) = cos*stor2+sin*stor1
            h(j) = cos*stor1-sin*stor2
 380      continue
 390    continue
 400    if(nrold.eq.number) go to 420
 410    nrold = nrold+1
        go to 300
 420  continue
!  backward substitution to obtain the b-spline coefficients as the
!  solution of the linear system    (ry) c (rx) = h.
!  first step: solve the system  (ry) (c1) = h.
      k = 1
      do 450 i=1,nk1x
        call fpback(ay,c(k),nk1y,ibandy,c(k),ny)
        k = k+nk1y
 450  continue
!  second step: solve the system  c (rx) = (c1).
      k = 0
      do 480 j=1,nk1y
        k = k+1
        l = k
        do 460 i=1,nk1x
          right(i) = c(l)
          l = l+nk1y
 460    continue
        call fpback(ax,right,nk1x,ibandx,right,nx)
        l = k
        do 470 i=1,nk1x
          c(l) = right(i)
          l = l+nk1y
 470    continue
 480  continue
!  calculate the quantities
!res(i,j) = (z(i,j) - s(x(i),y(j)))**2 , i=1,2,..,mx;j=1,2,..,my
!fp = sumi=1,mx(sumj=1,my(res(i,j)))
!fpx(r) = sumi(sumj=1,my(res(i,j))) , r=1,2,...,nx-2*kx-1
!          tx(r+kx) <= x(i) <= tx(r+kx+1)
!fpy(r) = sumi=1,mx(sumj(res(i,j))) , r=1,2,...,ny-2*ky-1
!          ty(r+ky) <= y(j) <= ty(r+ky+1)
      fp = 0.
      do 490 i=1,nx
        fpx(i) = 0.
 490  continue
      do 500 i=1,ny
        fpy(i) = 0.
 500  continue
      nk1y = ny-ky1
      nroldx = 0
!  main loop for the different grid points.
      do 550 i1=1,mx
        numx = nrx(i1)
        numx1 = numx+1
        nroldy = 0
        do 540 i2=1,my
          numy = nry(i2)
          numy1 = numy+1
!  evaluate s(x,y) at the current grid point by making the sum of the
!  cross products of the non-zero b-splines at (x,y), multiplied with
!  the appropiate b-spline coefficients.
          term = 0.
          k1 = numx*nk1y+numy
!CDIR LOOPCNT=6
          do 520 l1=1,kx1
            k2 = k1
            fac = spx(i1,l1)
!CDIR LOOPCNT=6
            do 510 l2=1,ky1
              k2 = k2+1
              term = term+fac*spy(i2,l2)*c(k2)
 510        continue
            k1 = k1+nk1y
 520      continue
!  calculate the squared residual at the current grid point.
          term = (z(i1,i2)-term)**2
!  adjust the different parameters.
          fp = fp+term
          fpx(numx1) = fpx(numx1)+term
          fpy(numy1) = fpy(numy1)+term
          fac = term*half
          if(numy.eq.nroldy) go to 530
          fpy(numy1) = fpy(numy1)-fac
          fpy(numy) = fpy(numy)+fac
 530      nroldy = numy
          if(numx.eq.nroldx) go to 540
          fpx(numx1) = fpx(numx1)-fac
          fpx(numx) = fpx(numx)+fac
 540    continue
        nroldx = numx
 550  continue
      return
      end

