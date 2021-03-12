      subroutine fpadno(maxtr,up,left,right,info,count,merk,jbind,&
       n1,ier)
!  subroutine fpadno adds a branch of length n1 to the triply linked
!  tree,the information of which is kept in the arrays up,left,right
!  and info. the information field of the nodes of this new branch is
!  given in the array jbind. in linking the new branch fpadno takes
!  account of the property of the tree that
!info(k) < info(right(k)) ; info(k) < info(left(k))
!  if necessary the subroutine calls subroutine fpfrno to collect the
!  free nodes of the tree. if no computer words are available at that
!  moment, the error parameter ier is set to 1.
!  ..
!  ..scalar arguments..
      integer maxtr,count,merk,n1,ier
!  ..array arguments..
      integer up(maxtr),left(maxtr),right(maxtr),info(maxtr),jbind(n1)
!  ..local scalars..
      integer k,niveau,point
      logical bool
!  ..subroutine references..
!fpfrno
!  ..
      point = 1
      niveau = 1
  10  k = left(point)
      bool = .true.
  20  if(k.eq.0) go to 50
      if(info(k)-jbind(niveau)) 30,40,50
  30  point = k
      k = right(point)
      bool = .false.
      go to 20
  40  point = k
      niveau = niveau+1
      go to 10
  50  if(niveau.gt.n1) go to 90
      count = count+1
      if(count.le.maxtr) go to 60
      call fpfrno(maxtr,up,left,right,info,point,merk,n1,count,ier)
      if(ier.ne.0) go to 100
  60  info(count) = jbind(niveau)
      left(count) = 0
      right(count) = k
      if(bool) go to 70
      bool = .true.
      right(point) = count
      up(count) = up(point)
      go to 80
  70  up(count) = point
      left(point) = count
  80  point = count
      niveau = niveau+1
      k = 0
      go to 50
  90  ier = 0
 100  return
      end
