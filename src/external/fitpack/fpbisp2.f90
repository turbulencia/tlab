      subroutine fpbisp2(tx,nx,ty,ny,c,kx,ky,xy,mxy,z)

      IMPLICIT NONE

#include "types.h"

!  ..scalar arguments..
      TINTEGER nx,ny,kx,ky,mxy
!  ..array arguments..
      TREAL tx(nx),ty(ny),c((nx-kx-1)*(ny-ky-1))
      TREAL xy(2,mxy), z(mxy)
!  ..local scalars..
      TINTEGER kx1,ky1,l,l1,l2,nkx1,nky1,ij,i1,j1
      TREAL argx,argy,sp,tbx,tex,tby,tey
!  ..local arrays..
      TREAL h(6,2), r1, hh(5), f
      TINTEGER lx, ly, li, lj, i, j
!  ..subroutine references..
!fpbspl
! ..
      kx1 = kx+1
      nkx1 = nx-kx1
      ky1 = ky+1
      nky1 = ny-ky1
      tbx = tx(kx1)
      tex = tx(nkx1+1)
      tby = ty(ky1)
      tey = ty(nky1+1)
      r1 = C_1_R

      do 40 ij=1,mxy
         l = kx1
         l1 = l+1
         argx = xy(1,ij)
         argx = MAX(tbx, argx)
! if(argx.lt.tbx) argx = tbx
         argx = MIN(tex, argx)
! if(argx.gt.tex) argx = tex
! 10      if(argx.lt.tx(l1) .or. l.eq.nkx1) go to 20
! l = l1
! l1 = l+1
! go to 10
         
         DO l=kx1, nkx1-1
            IF ( argx.lt.tx(l+1) ) GOTO 20
         ENDDO
! 20      call fpbspl(tx,nx,kx,argx,l,h(1,1))

 20      h(1,1) = r1
!CDIR LOOPCNT=6
         DO j=1,kx
!CDIR LOOPCNT=6
            DO i=1,j
               hh(i) = h(i,1)
            ENDDO
            h(1,1) = C_0_R
!CDIR LOOPCNT=6
            DO i=1,j
               li = l+i
               lj = li-j
               f = hh(i)/(tx(li)-tx(lj))
               h(i,1) = h(i,1)+f*(tx(li)-argx)
               h(i+1,1) = f*(argx-tx(lj))
            ENDDO
         ENDDO

         lx = l-kx1

         l = ky1
         l1 = l+1
         argy = xy(2,ij)
         argy = MAX(tby, argy)
! if(argy.lt.tby) argy = tby
         argy = MIN(tey, argy)
! if(argy.gt.tey) argy = tey
! 50      if(argy.lt.ty(l1) .or. l.eq.nky1) go to 60
! l = l1
! l1 = l+1
! go to 50
         
         DO l=ky1, nky1-1
            IF ( argy.lt.ty(l+1) ) GOTO 60
         ENDDO
! 60      call fpbspl(ty,ny,ky,argy,l,h(1,2))
         
 60      h(1,2) = r1
!CDIR LOOPCNT=6
         DO j=1,ky
!CDIR LOOPCNT=6
            DO i=1,j
               hh(i) = h(i,2)
            ENDDO
            h(1,2) = C_0_R
!CDIR LOOPCNT=6
            DO i=1,j
               li = l+i
               lj = li-j
               f = hh(i)/(ty(li)-ty(lj))
               h(i,2) = h(i,2)+f*(ty(li)-argy)
               h(i+1,2) = f*(argy-ty(lj))
            ENDDO
         ENDDO


         ly = l-ky1

         l1 = lx*nky1+ly
         sp = C_0_R
         
!CDIR LOOPCNT=6
         do i1=1,kx1
            l2 = l1
!CDIR LOOPCNT=6
            do j1=1,ky1
               l2 = l2+1
               sp = sp+c(l2)*h(i1,1)*h(j1,2)
            ENDDO
            l1 = l1+nky1
         ENDDO
         z(ij) = sp

! PRINT *, ij, xy(1,ij), xy(2,ij), z(ij)
 40   continue
      
      return
      end

