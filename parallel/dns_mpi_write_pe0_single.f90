#include "types.h"
#include "dns_error.h"
#include "dns_const_mpi.h"

!########################################################################
!# Tool/Library DNS
!#
!########################################################################
!# HISTORY
!#
!# 2008/09/09 - J.P. Mellado
!#              Created
!# 2012/12/01 - J.P. Mellado
!#              Modified for 2D decomposition
!#
!########################################################################
!# DESCRIPTION
!#
!########################################################################
SUBROUTINE DNS_MPI_WRITE_PE0_SINGLE(iunit, nx,ny,nz, subdomain, u, tmp1, tmp2)

  USE DNS_CONSTANTS, ONLY : efile
  USE DNS_MPI

  IMPLICIT NONE

#include "mpif.h"

  TINTEGER iunit, nx,ny,nz, subdomain(6)
  TREAL, DIMENSION(nx*ny*nz),        TARGET :: u, tmp1
  TREAL, DIMENSION(nx*ims_npro_i,*), TARGET :: tmp2

! -------------------------------------------------------------------
  TINTEGER nx_total, ny_total, nz_total
  TINTEGER nx_min, nx_max, ny_min, ny_max, nz_min, nz_max
  TINTEGER nyz

  TINTEGER ip_i, ip_k, joffset_loc, koffset_loc, id
  TINTEGER i, jk, j_loc, k_loc
  INTEGER mpio_size, mpio_ip
  INTEGER status(MPI_STATUS_SIZE)

  TREAL, DIMENSION(:), POINTER :: p_org

! ###################################################################
  nx_total = nx*ims_npro_i
  ny_total = ny
  nz_total = nz*ims_npro_k

  nx_min = subdomain(1); nx_max = subdomain(2)
  ny_min = subdomain(3); ny_max = subdomain(4)
  nz_min = subdomain(5); nz_max = subdomain(6)

  koffset_loc = 0
  joffset_loc = 0

  id = DNS_MPI_I_PARTIAL

! -------------------------------------------------------------------
! Transposing along Ox
! -------------------------------------------------------------------
  IF ( ims_npro_i .GT. 1 ) THEN
     CALL DNS_MPI_TRPF_I(u, tmp1, ims_ds_i(1,id), ims_dr_i(1,id), ims_ts_i(1,id), ims_tr_i(1,id))
     p_org => tmp1
     nyz = ims_size_i(id)
  ELSE
     p_org => u
     nyz = ny*nz
  ENDIF
  mpio_size = nyz*nx_total

  CALL MPI_BARRIER(MPI_COMM_WORLD, ims_err)

! -------------------------------------------------------------------
! Passing all data through PE#0
! -------------------------------------------------------------------
  IF ( ims_pro .EQ. 0 ) THEN

     DO ip_k = 1,ims_npro_k
        koffset_loc = nz *(ip_k-1)

        DO ip_i = 1,ims_npro_i
           joffset_loc = nyz *(ip_i-1) ! Remember that data is Ox-transposed

           mpio_ip = ims_npro_i *(ip_k-1) + ip_i-1
           IF ( mpio_ip .EQ. 0 ) THEN
              tmp2(1:mpio_size,1) = p_org(1:mpio_size)
           ELSE
              CALL MPI_RECV(tmp2, mpio_size, MPI_REAL8, mpio_ip, ims_tag, MPI_COMM_WORLD, status, ims_err)
           ENDIF

           DO jk = 1,nyz
              j_loc = MOD( (jk-1+joffset_loc) , ny_total ) + 1
              k_loc =    ( (jk-1+joffset_loc) / ny_total ) + 1 + koffset_loc

              IF ( ( j_loc .GE. ny_min ) .AND. ( j_loc .LE. ny_max )  .AND. &
                   ( k_loc .GE. nz_min ) .AND. ( k_loc .LE. nz_max ) ) THEN
                 WRITE(iunit) (SNGL(tmp2(i,jk)),i=nx_min,nx_max)
              ENDIF
              
           ENDDO

        ENDDO
     ENDDO

  ELSE
     CALL MPI_SEND(p_org, mpio_size, MPI_REAL8, 0, ims_tag, MPI_COMM_WORLD, ims_err)
  ENDIF

  RETURN
END SUBROUTINE DNS_MPI_WRITE_PE0_SINGLE
