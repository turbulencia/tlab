!# Compile with
!# mpxlf90_r -d -qsuffix=cpp=f90 -WF,-DUSE_MPI -O3 -qhot -qstrict -qsave -q64 -bdatapsize:64K -bstackpsize:64K -btextpsize:64K -o vmpi.x vmpi.f90

#define TREAL      REAL(8)
#define TINTEGER   INTEGER(4)
#define C_0_R      0.0d0

!########################################################################
!# Tool/Library
!#
!########################################################################
!# HISTORY
!#
!########################################################################
!# DESCRIPTION
!#
!########################################################################
!# ARGUMENTS
!#
!########################################################################
MODULE DNS_MPI
  IMPLICIT NONE
  SAVE

  INTEGER ims_pro, ims_npro, ims_err, ims_tag

  TINTEGER, PARAMETER :: MAX_MPI_PROC  = 2048
  TINTEGER, PARAMETER :: MAX_MPI_TYPES =   10

  TINTEGER nloc(MAX_MPI_TYPES)

  TINTEGER ndisp(MAX_MPI_PROC, MAX_MPI_TYPES)
  TINTEGER mdisp(MAX_MPI_PROC, MAX_MPI_TYPES)

  TINTEGER nip(MAX_MPI_PROC, MAX_MPI_TYPES)
  TINTEGER mip(MAX_MPI_PROC)

  INTEGER tfsend(MAX_MPI_PROC, MAX_MPI_TYPES)
  INTEGER tfrecv(MAX_MPI_PROC, MAX_MPI_TYPES)

END MODULE DNS_MPI

!########################################################################
!# Tool/Library
!#
!########################################################################
!# HISTORY
!#
!########################################################################
!# DESCRIPTION
!#
!########################################################################
!# ARGUMENTS
!#
!########################################################################
PROGRAM VMPI

  USE TLAB_MPI_VARS

  IMPLICIT NONE

#ifdef USE_MPI
#include "mpif.h"
#endif

  TINTEGER imax, jmax, kmax, kmax_total
  TINTEGER inb_flow, inb_flow_array, inb_scal, inb_scal_array, n_txc_size
  TINTEGER ip, it, is, it_max, is_max
  TINTEGER isize, i1, ierr

  TREAL, DIMENSION(:,:), ALLOCATABLE, SAVE :: q, z1, h, zh1
  TREAL, DIMENSION(:),   ALLOCATABLE, SAVE :: wrk1d, wrk2d, wrk3d, txc

  TARGET q

! Pointers to existing allocated sspace
  TREAL, DIMENSION(:),   POINTER           :: u, v, w

! ###################################################################
  i1 = 1

  imax       = 2048
  jmax       = 1536
  kmax       = 1
  kmax_total = 2048

  inb_flow       = 3
  inb_flow_array = inb_flow
  inb_scal        = 1
  inb_scal_array  = inb_scal
  n_txc_size         = imax*jmax*kmax*6

  ALLOCATE(wrk1d(MAX(imax,MAX(jmax,kmax_total))         *15))
  ALLOCATE(wrk2d(MAX(imax*jmax,MAX(imax*kmax,jmax*kmax))* 6))
  ALLOCATE(wrk3d(imax*jmax*kmax))

! big arrays
  ALLOCATE(q(imax*jmax*kmax,  inb_flow_array),STAT=ierr)
  ALLOCATE(h(imax*jmax*kmax,  inb_flow),      STAT=ierr)
  ALLOCATE(z1(imax*jmax*kmax, inb_scal_array),STAT=ierr)
  ALLOCATE(zh1(imax*jmax*kmax,inb_scal),      STAT=ierr)
  ALLOCATE(txc(n_txc_size),                   STAT=ierr)

! pointers
  u   => q(:,1)
  v   => q(:,2)
  w   => q(:,3)

! ###################################################################
! from DNS_START
  call MPI_INIT(ims_err)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,ims_npro,ims_err)
  call MPI_COMM_RANK(MPI_COMM_WORLD,ims_pro, ims_err)

! ###################################################################
! from DNS_READ_GLOBAL
  DO ip = 1,ims_npro
     mip(ip) = kmax
  ENDDO

! ###################################################################
! from TLAB_MPI_INIT
  isize = imax*jmax
  CALL TLAB_MPI_TYPE(ims_npro, ims_pro, kmax_total, isize,              &
       i1, i1, i1, i1, mip, nip(1,1), nloc(1), ndisp(1,1), mdisp(1,1), &
       tfsend(1,1), tfrecv(1,1))

  CALL TLAB_MPI_TAGRESET

! ###################################################################
  q(:,:)  = C_0_R; h(:,:)   = C_0_R
  z1(:,:) = C_0_R; zh1(:,:) = C_0_R

  it_max = 1000
  is_max = 5

  DO it = 0,it_max
     IF ( MOD(it,10) .EQ. 0 ) THEN
        IF ( ims_pro .EQ. 0 ) THEN
           OPEN(UNIT=22, FILE='dns.out', STATUS='unknown',POSITION='APPEND')
           WRITE(22,*) 'Iteration ...:', it, ims_tag
           CLOSE(22)
        ENDIF
     ENDIF

     DO is = 1,is_max
        CALL TLAB_MPI_TRPF_K(u, wrk3d, ndisp(1,1), mdisp(1,1), tfsend(1,1), tfrecv(1,1))
        CALL TLAB_MPI_TRPB_K(wrk3d, u, ndisp(1,1), mdisp(1,1), tfsend(1,1), tfrecv(1,1))
        CALL TLAB_MPI_TRPF_K(u, wrk3d, ndisp(1,1), mdisp(1,1), tfsend(1,1), tfrecv(1,1))
        CALL TLAB_MPI_TRPB_K(wrk3d, u, ndisp(1,1), mdisp(1,1), tfsend(1,1), tfrecv(1,1))

        CALL TLAB_MPI_TRPF_K(v, wrk3d, ndisp(1,1), mdisp(1,1), tfsend(1,1), tfrecv(1,1))
        CALL TLAB_MPI_TRPB_K(wrk3d, v, ndisp(1,1), mdisp(1,1), tfsend(1,1), tfrecv(1,1))
        CALL TLAB_MPI_TRPF_K(v, wrk3d, ndisp(1,1), mdisp(1,1), tfsend(1,1), tfrecv(1,1))
        CALL TLAB_MPI_TRPB_K(wrk3d, v, ndisp(1,1), mdisp(1,1), tfsend(1,1), tfrecv(1,1))

        CALL TLAB_MPI_TRPF_K(w, wrk3d, ndisp(1,1), mdisp(1,1), tfsend(1,1), tfrecv(1,1))
        CALL TLAB_MPI_TRPB_K(wrk3d, w, ndisp(1,1), mdisp(1,1), tfsend(1,1), tfrecv(1,1))
        CALL TLAB_MPI_TRPF_K(w, wrk3d, ndisp(1,1), mdisp(1,1), tfsend(1,1), tfrecv(1,1))
        CALL TLAB_MPI_TRPB_K(wrk3d, w, ndisp(1,1), mdisp(1,1), tfsend(1,1), tfrecv(1,1))

        CALL TLAB_MPI_TRPF_K(txc, wrk3d, ndisp(1,1), mdisp(1,1), tfsend(1,1), tfrecv(1,1))
        CALL TLAB_MPI_TRPB_K(wrk3d, txc, ndisp(1,1), mdisp(1,1), tfsend(1,1), tfrecv(1,1))

        CALL TLAB_MPI_TRPF_K(z1, wrk3d, ndisp(1,1), mdisp(1,1), tfsend(1,1), tfrecv(1,1))
        CALL TLAB_MPI_TRPB_K(wrk3d, z1, ndisp(1,1), mdisp(1,1), tfsend(1,1), tfrecv(1,1))
        CALL TLAB_MPI_TRPF_K(z1, wrk3d, ndisp(1,1), mdisp(1,1), tfsend(1,1), tfrecv(1,1))
        CALL TLAB_MPI_TRPB_K(wrk3d, z1, ndisp(1,1), mdisp(1,1), tfsend(1,1), tfrecv(1,1))

     ENDDO

  ENDDO

  CALL MPI_FINALIZE(ims_err)

END PROGRAM VMPI

!########################################################################
!# Tool/Library
!#
!########################################################################
!# HISTORY
!#
!########################################################################
!# DESCRIPTION
!#
!########################################################################
!# ARGUMENTS
!#
!########################################################################
SUBROUTINE TLAB_MPI_TYPE(ims_npro, ims_pro, kmax_total, nblock, &
     nd, md, n1, n2, mip, nip, nloc, ndisp, mdisp, tfsend, tfrecv)

  IMPLICIT NONE

#ifdef USE_MPI
#include "mpif.h"
#endif

  INTEGER ims_npro, ims_pro
  TINTEGER nblock, kmax_total
  TINTEGER nd, md, n1, n2, nloc
  TINTEGER mip(*), nip(*), ndisp(*), mdisp(*)
  INTEGER tfsend(*), tfrecv(*)

! -----------------------------------------------------------------------
  TINTEGER i, kmax

#ifdef USE_MPI
  INTEGER ims_ss, ims_rs, ims_err
  INTEGER ims_tmp1, ims_tmp2, ims_tmp3
#endif

! #######################################################################
  DO i=1,ims_npro
     nip(i) = ((nblock*mip(i))/kmax_total)
  ENDDO

  kmax = mip(ims_pro+1)
  nloc = nip(ims_pro+1)

! Calculate Displacements in Forward Send/Receive
  ndisp(1) = 0
  mdisp(1) = 0
  DO i=2,ims_npro
     ndisp(i) = ndisp(i-1) + nd * nip(i-1)
     mdisp(i) = mdisp(i-1) + nloc * md * mip(i-1)
  ENDDO

! #######################################################################
! nloc/1/kmax_total
#ifdef USE_MPI
  DO i=1, ims_npro

     ims_tmp1 = kmax*n1
     ims_tmp2 = nip(i)*n2
     ims_tmp3 = nblock*n2
     CALL MPI_TYPE_VECTOR(ims_tmp1, ims_tmp2, ims_tmp3, MPI_REAL8, tfsend(i), ims_err)
     CALL MPI_TYPE_COMMIT(tfsend(i), ims_err)

     ims_tmp1 = mip(i)*n1
     ims_tmp2 = nloc*n2
     CALL MPI_TYPE_VECTOR(ims_tmp1, ims_tmp2, ims_tmp2, MPI_REAL8, tfrecv(i), ims_err)
     CALL MPI_TYPE_COMMIT(tfrecv(i), ims_err)

     CALL MPI_TYPE_SIZE(tfsend(i), ims_ss, ims_err)
     CALL MPI_TYPE_SIZE(tfrecv(i), ims_rs, ims_err)

     IF ( ims_ss .NE. ims_rs ) THEN
        PRINT *, 'Message   : ', i, ' size is wrong'
        PRINT *, 'Send size : ', ims_ss
        PRINT *, 'Recv size : ', ims_rs
     ENDIF
  ENDDO
#endif

  RETURN
END SUBROUTINE TLAB_MPI_TYPE

!########################################################################
!# Tool/Library
!#
!########################################################################
!# HISTORY
!#
!#
!########################################################################
!# DESCRIPTION
!#
!########################################################################
!# ARGUMENTS
!#
!########################################################################
SUBROUTINE TLAB_MPI_TRPF_K(a, b, ndsp, mdsp, tsend, trecv)

  USE TLAB_MPI_VARS

  IMPLICIT NONE

#ifdef USE_MPI
#include "mpif.h"
#endif

  TREAL a(*)
  TREAL b(*)
  TINTEGER ndsp(ims_npro)
  TINTEGER mdsp(ims_npro)
  INTEGER tsend(ims_npro)
  INTEGER trecv(ims_npro)

! -----------------------------------------------------------------------
#ifdef USE_MPI
  TINTEGER i, l
  INTEGER status(MPI_STATUS_SIZE,2*MAX_MPI_PROC)
  INTEGER mpireq(2*MAX_MPI_PROC)
  INTEGER ip
#endif

! #######################################################################

#ifdef USE_MPI

! #######################################################################
! Same processor
! #######################################################################
!  PRINT *, ims_pro, ims_pro, ndsp(ims_pro+1), tsend(ims_pro+1), ims_tag
  CALL MPI_ISEND(a(ndsp(ims_pro+1)+1), 1, tsend(ims_pro+1), ims_pro,&
       ims_tag, MPI_COMM_WORLD, mpireq(1), ims_err)

!  PRINT *, ims_pro, ims_pro, mdsp(ims_pro+1), trecv(ims_pro+1), ims_tag
  CALL MPI_IRECV(b(mdsp(ims_pro+1)+1), 1, trecv(ims_pro+1), ims_pro,&
       ims_tag, MPI_COMM_WORLD, mpireq(2), ims_err)

  CALL MPI_WAITALL(2, mpireq, status, ims_err)

! #######################################################################
! Different processors
! #######################################################################
  l = 3
  DO i=1, ims_npro
     ip = i-1
     IF ( ip .NE. ims_pro ) THEN
!        PRINT *, ims_pro, ip, ndsp(i), tsend(i), ims_tag
        CALL MPI_ISEND(a(ndsp(i)+1), 1, tsend(i), ip,&
             ims_tag, MPI_COMM_WORLD, mpireq(l), ims_err)
        l = l + 1

!        PRINT *, ims_pro, ip, mdsp(i), trecv(i), ims_tag
        CALL MPI_IRECV(b(mdsp(i)+1), 1, trecv(i), ip,&
             ims_tag, MPI_COMM_WORLD, mpireq(l), ims_err)
        l = l + 1

     ENDIF
  ENDDO

  CALL MPI_WAITALL(ims_npro*2-2, mpireq(3), status(1,3), ims_err)

  CALL TLAB_MPI_TAGUPDT

#endif

  RETURN
END SUBROUTINE TLAB_MPI_TRPF_K

!########################################################################
!# Tool/Library
!#
!########################################################################
!# HISTORY
!#
!########################################################################
!# DESCRIPTION
!#
!########################################################################
!# ARGUMENTS
!#
!########################################################################
SUBROUTINE TLAB_MPI_TRPB_K(b, a, ndsp, mdsp, tsend, trecv)

  USE TLAB_MPI_VARS

  IMPLICIT NONE

#ifdef USE_MPI
#include "mpif.h"
#endif

  TREAL a(*)
  TREAL b(*)
  TINTEGER ndsp(ims_npro)
  TINTEGER mdsp(ims_npro)
  INTEGER tsend(ims_npro)
  INTEGER trecv(ims_npro)

#ifdef USE_MPI
  INTEGER status(MPI_STATUS_SIZE,2*MAX_MPI_PROC)
  INTEGER mpireq(2*MAX_MPI_PROC)
  INTEGER ip
  TINTEGER i, l
#endif

! #######################################################################
#ifdef USE_MPI

! #######################################################################
! Same processor
! #######################################################################
  CALL MPI_ISEND(b(mdsp(ims_pro+1)+1), 1, trecv(ims_pro+1), &
       ims_pro, ims_tag, MPI_COMM_WORLD, mpireq(1), ims_err)

  CALL MPI_IRECV(a(ndsp(ims_pro+1)+1), 1, tsend(ims_pro+1), &
       ims_pro, ims_tag, MPI_COMM_WORLD, mpireq(2), ims_err)

  CALL MPI_WAITALL(2, mpireq, status, ims_err)

! #######################################################################
! Different processors
! #######################################################################
  l = 3
  DO i=1, ims_npro
     ip = i-1
     IF ( ip .NE. ims_pro ) THEN
        CALL MPI_ISEND(b(mdsp(i)+1), 1, trecv(i), ip,&
             ims_tag, MPI_COMM_WORLD, mpireq(l), ims_err)
        l = l + 1

        CALL MPI_IRECV(a(ndsp(i)+1), 1, tsend(i), ip,&
             ims_tag, MPI_COMM_WORLD, mpireq(l), ims_err)
        l = l + 1

     ENDIF
  ENDDO

  CALL MPI_WAITALL(ims_npro*2-2, mpireq(3), status(1,3), ims_err)

  CALL TLAB_MPI_TAGUPDT

#endif

  RETURN
END SUBROUTINE TLAB_MPI_TRPB_K

!########################################################################
!# Tool/Library
!#
!########################################################################
!# HISTORY
!#
!########################################################################
!# DESCRIPTION
!#
!########################################################################
!# ARGUMENTS
!#
!########################################################################
SUBROUTINE TLAB_MPI_TAGUPDT

  USE TLAB_MPI_VARS, ONLY : ims_tag

  IMPLICIT NONE

  ims_tag = ims_tag + 1

  IF ( ims_tag .GT. 32000 ) THEN
     CALL TLAB_MPI_TAGRESET
  ENDIF

  RETURN
END SUBROUTINE TLAB_MPI_TAGUPDT

!########################################################################
!# Tool/Library
!#
!########################################################################
!# HISTORY
!#
!########################################################################
!# DESCRIPTION
!#
!########################################################################
!# ARGUMENTS
!#
!########################################################################
SUBROUTINE TLAB_MPI_TAGRESET

  USE TLAB_MPI_VARS, ONLY : ims_tag

  IMPLICIT NONE

  ims_tag = 0

  RETURN
END SUBROUTINE TLAB_MPI_TAGRESET
