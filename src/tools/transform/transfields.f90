#include "types.h"
#include "dns_const.h"
#include "dns_error.h"
#ifdef USE_MPI
#include "dns_const_mpi.h"
#endif

#define C_FILE_LOC "TRANSFORM"

PROGRAM TRANSFIELDS

  USE DNS_TYPES,  ONLY : filter_dt, grid_dt
  USE DNS_CONSTANTS
  USE DNS_GLOBAL
  USE TLAB_ARRAYS
#ifdef USE_MPI
  USE DNS_MPI
#endif

  IMPLICIT NONE

#include "integers.h"
#ifdef USE_MPI
#include "mpif.h"
#endif

! Parameter definitions
  TINTEGER, PARAMETER :: itime_size_max = 3000
  TINTEGER, PARAMETER :: iopt_size_max  = 512

  ! -------------------------------------------------------------------
  ! Additional local arrays
  TREAL, DIMENSION(:),   ALLOCATABLE, SAVE         :: x_dst,y_dst,z_dst
  TREAL, DIMENSION(:,:), ALLOCATABLE, SAVE         :: q_dst, s_dst

  TREAL, DIMENSION(:),     ALLOCATABLE, SAVE :: y_aux
  TREAL, DIMENSION(:,:,:), ALLOCATABLE, SAVE :: txc_aux

! -------------------------------------------------------------------
! Local variables
! -------------------------------------------------------------------
  TINTEGER opt_main, opt_function
  TINTEGER iq, is, ig, ip, j,k
  TINTEGER idummy, iread_flow, iread_scal, ierr
  CHARACTER*32 bakfile, flow_file, scal_file
  CHARACTER*64 str, line
  CHARACTER*512 sRes
  TINTEGER subdomain(6)

  TYPE(grid_dt), DIMENSION(3) :: g_dst
  TINTEGER imax_dst,jmax_dst,kmax_dst

  LOGICAL flag_crop, flag_extend
  TINTEGER jmax_aux, inb_scal_dst
  TREAL dummy

  TINTEGER itime_size, it
  TINTEGER itime_vec(itime_size_max)

  TINTEGER iopt_size
  TREAL opt_vec(iopt_size_max)

! ###################################################################
  bakfile = TRIM(ADJUSTL(ifile))//'.bak'

  CALL DNS_START

  CALL DNS_READ_GLOBAL(ifile)

#ifdef USE_MPI
  CALL DNS_MPI_INITIALIZE
#endif

! -------------------------------------------------------------------
! Allocating memory space
! -------------------------------------------------------------------

! -------------------------------------------------------------------
! File names
! -------------------------------------------------------------------
#include "dns_read_times.h"

! -------------------------------------------------------------------
! Read local options
! -------------------------------------------------------------------
  opt_main     =-1 ! default values
  opt_function = 0

  CALL SCANINICHAR(bakfile, ifile, 'PostProcessing', 'ParamTransform', '-1', sRes)
  iopt_size = iopt_size_max
  CALL LIST_REAL(sRes, iopt_size, opt_vec)

  IF ( sRes .EQ. '-1' ) THEN
#ifdef USE_MPI
#else
     WRITE(*,'(A)') 'Option ?'
     WRITE(*,'(A)') '1. Crop fields'
     WRITE(*,'(A)') '2. Extend fields in Ox and Oy'
     WRITE(*,'(A)') '3. Remesh fields'
     WRITE(*,'(A)') '4. Linear combination of fields'
     WRITE(*,'(A)') '5. Filter fields'
     WRITE(*,'(A)') '6. Transform scalar fields'
     WRITE(*,'(A)') '7. Blend fields'
     WRITE(*,'(A)') '8. Add mean profiles'
     READ(*,*) opt_main
#endif
  ELSE
     opt_main = INT(opt_vec(1))
  ENDIF

! -------------------------------------------------------------------
  IF ( opt_main .EQ. 4 ) THEN
     IF ( sRes .EQ. '-1' ) THEN
#ifdef USE_MPI
#else
        WRITE(*,*) 'Coefficients ?'
        READ(*,'(A512)') sRes
        iopt_size = iopt_size_max-1
        CALL LIST_REAL(sRes, iopt_size, opt_vec(2))
#endif
     ELSE
        iopt_size = iopt_size-1
     ENDIF
     IF ( iopt_size .EQ. 0 ) THEN
        CALL IO_WRITE_ASCII(lfile,'TRANSFORM. Performing arithmetic mean of fields.')
        iopt_size = itime_size
        opt_vec(2:) = C_1_R /M_REAL(itime_size)
     ENDIF
     IF ( iopt_size .NE. itime_size ) THEN
        CALL IO_WRITE_ASCII(efile,'TRANSFORM. Number of coefficient incorrect.')
        CALL DNS_STOP(DNS_ERROR_UNDEVELOP)
     ENDIF
     IF ( iopt_size .GT. iopt_size_max ) THEN
        CALL IO_WRITE_ASCII(efile,'TRANSFORM. Array opt_vec too small.')
        CALL DNS_STOP(DNS_ERROR_UNDEVELOP)
     ENDIF
  ENDIF

! -------------------------------------------------------------------
  IF ( opt_main .EQ. 5 .AND. FilterDomain(1)%type .EQ. DNS_FILTER_NONE ) THEN
     CALL IO_WRITE_ASCII(efile,'TRANSFORM. Filter information needs to be provided in block [Filter].')
     CALL DNS_STOP(DNS_ERROR_OPTION)
  ENDIF

! -------------------------------------------------------------------
  IF ( opt_main .EQ. 6 ) THEN
     IF ( sRes .EQ. '-1' ) THEN
#ifdef USE_MPI
#else
        WRITE(*,'(A)') 'Function type ?'
        WRITE(*,'(A)') '1. Vapor Saturation Humidity'
        WRITE(*,'(A)') '2. Linear'
        READ(*,*) opt_function
#endif
     ELSE
        opt_function = INT(opt_vec(2))
     ENDIF

     IF ( sRes .EQ. '-1' ) THEN
#ifdef USE_MPI
#else
        WRITE(*,*) 'Coefficients ?'
        READ(*,'(A512)') sRes
        iopt_size = iopt_size_max-2
        CALL LIST_REAL(sRes, iopt_size, opt_vec(3))
#endif
     ELSE
        iopt_size = iopt_size-2
     ENDIF
  ENDIF

! -------------------------------------------------------------------
  IF ( opt_main .EQ. 7 ) THEN ! 2nd and 3rd entries in opt_vec contain coeffs.
     IF ( sRes .EQ. '-1' ) THEN
        WRITE(*,*) 'Coefficients ?'
        READ(*,'(A512)') sRes
        iopt_size = 2
        CALL LIST_REAL(sRes, iopt_size, opt_vec(2))
     ELSE
        iopt_size = iopt_size-1
     ENDIF
     IF ( iopt_size .NE. 2 ) THEN
        CALL IO_WRITE_ASCII(efile,'TRANSFORM. Number of blend coefficient incorrect.')
        CALL DNS_STOP(DNS_ERROR_UNDEVELOP)
     ENDIF
  ENDIF

  IF ( opt_main .LT. 0 ) THEN ! Check
     CALL IO_WRITE_ASCII(efile, 'TRANSFORM. Missing input [ParamTransform] in dns.ini.')
     CALL DNS_STOP(DNS_ERROR_INVALOPT)
  ENDIF

  IF ( opt_main .EQ. 6 ) THEN; icalc_flow = 0; ENDIF ! Force not to process the flow fields

! -------------------------------------------------------------------
  CALL SCANINICHAR(bakfile, ifile, 'PostProcessing', 'Subdomain', '-1', sRes)

  IF ( sRes .EQ. '-1' ) THEN
#ifdef USE_MPI
#else
     WRITE(*,*) 'Subdomain limits ?'
     READ(*,'(A64)') sRes
#endif
  ENDIF
  idummy = 6
  CALL LIST_INTEGER(sRes, idummy, subdomain)

  IF ( idummy .LT. 6 ) THEN ! default
     subdomain(1) = 1; subdomain(2) = g(1)%size
     subdomain(3) = 1; subdomain(4) = g(2)%size
     subdomain(5) = 1; subdomain(6) = g(3)%size
  ENDIF

! -------------------------------------------------------------------
  IF      ( opt_main .EQ. 1 .OR. &      ! Crop
            opt_main .EQ. 3      ) THEN ! Remesh
     g_dst(1)%size = subdomain(2)-subdomain(1)+1
     g_dst(2)%size = subdomain(4)-subdomain(3)+1
     g_dst(3)%size = subdomain(6)-subdomain(5)+1

  ELSE IF ( opt_main .EQ. 2      ) THEN ! Extend
     g_dst(1)%size = g(1)%size + subdomain(2) + subdomain(1)
     g_dst(2)%size = g(2)%size + subdomain(4) + subdomain(3)
     g_dst(3)%size = g(3)%size

  ELSE
     g_dst(1)%size = g(1)%size
     g_dst(2)%size = g(2)%size
     g_dst(3)%size = g(3)%size
  ENDIF

#ifdef USE_MPI
  imax_dst = g_dst(1)%size/ims_npro_i
  jmax_dst = g_dst(2)%size
  kmax_dst = g_dst(3)%size/ims_npro_k
#else
  imax_dst = g_dst(1)%size
  jmax_dst = g_dst(2)%size
  kmax_dst = g_dst(3)%size
#endif

! -------------------------------------------------------------------
! Further allocation of memory space
! -------------------------------------------------------------------
  inb_txc = 0

  inb_scal_dst = inb_scal

  iread_flow = icalc_flow
  iread_scal = icalc_scal

  IF      ( opt_main .EQ. 3 ) THEN ! Remesh
     isize_txc_field = MAX(isize_txc_field,imax_dst*jmax_dst*kmax_dst)
     inb_txc         = 5
  ELSE IF ( opt_main .EQ. 5 ) THEN ! Filter
     inb_txc         = 4
  ELSE IF ( opt_main .EQ. 6 ) THEN
     inb_txc         = 5
     inb_scal_dst    = 1
  ENDIF

  IF ( ifourier .EQ. 1 ) inb_txc = MAX(inb_txc,1)

! -------------------------------------------------------------------
  isize_wrk3d = MAX(isize_txc_field,imax_dst*jmax_dst*kmax_dst)

! -------------------------------------------------------------------
  IF ( icalc_flow .EQ. 1 ) ALLOCATE(q_dst(imax_dst*jmax_dst*kmax_dst,inb_flow))
  IF ( icalc_scal .EQ. 1 ) ALLOCATE(s_dst(imax_dst*jmax_dst*kmax_dst,inb_scal_dst))

  IF ( opt_main .EQ. 3 ) THEN
     ALLOCATE(x_dst(g_dst(1)%size))
     ALLOCATE(y_dst(g_dst(2)%size))
     ALLOCATE(z_dst(g_dst(3)%size))
  ENDIF

  CALL TLAB_ALLOCATE(C_FILE_LOC)

! -------------------------------------------------------------------
! Read the grid
! -------------------------------------------------------------------
#include "dns_read_grid.h"

! -------------------------------------------------------------------
! Initialize filters
! -------------------------------------------------------------------
  IF ( opt_main .EQ. 5 ) THEN
     DO ig = 1,3
        CALL OPR_FILTER_INITIALIZE( g(ig), FilterDomain(ig), wrk1d )
     END DO
  ENDIF

! -------------------------------------------------------------------
! Initialize Poisson solver
! -------------------------------------------------------------------
  IF ( ifourier .EQ. 1 ) CALL OPR_FOURIER_INITIALIZE(txc, wrk1d,wrk2d,wrk3d)

  IF ( inb_txc .GE. 3 .AND. icalc_flow .GT. 0 ) CALL OPR_CHECK(imax,jmax,kmax, q, txc, wrk2d,wrk3d)

! -------------------------------------------------------------------
! Initialize cumulative field
! -------------------------------------------------------------------
  IF ( opt_main .EQ. 4 .OR. opt_main .EQ. 7 ) THEN
     IF ( icalc_flow .EQ. 1 ) q_dst = C_0_R
     IF ( icalc_scal .EQ. 1 ) s_dst = C_0_R
  ENDIF

! -------------------------------------------------------------------
! Initialize thermodynamic quantities
! -------------------------------------------------------------------
  CALL FI_PROFILES_INITIALIZE(wrk1d)

! -------------------------------------------------------------------
! Initialize remeshing
! -------------------------------------------------------------------
  IF ( opt_main .EQ. 3 ) THEN
     CALL IO_READ_GRID('grid.trn', g_dst(1)%size,g_dst(2)%size,g_dst(3)%size, &
          g_dst(1)%scale,g_dst(2)%scale,g_dst(3)%scale, x_dst,y_dst,z_dst)

! Check grids; Ox and Oz directions are assumed to be periodic
     dummy = (g_dst(1)%scale-g(1)%scale) / (x(g(1)%size,1)-x(g(1)%size-1,1))
     IF ( ABS(dummy) .GT. C_1EM3_R ) THEN
        CALL IO_WRITE_ASCII(efile, 'TRANSFORM. Ox scales are not equal at the end.')
        CALL DNS_STOP(DNS_ERROR_GRID_SCALE)
     ENDIF
     wrk1d(1:g(1)%size,1) = x(1:g(1)%size,1) ! we need extra space

     dummy = (g_dst(3)%scale-g(3)%scale) / (z(g(3)%size,1)-z(g(3)%size-1,1))
     IF ( ABS(dummy) .GT. C_1EM3_R ) THEN
        CALL IO_WRITE_ASCII(efile, 'TRANSFORM. Oz scales are not equal')
        CALL DNS_STOP(DNS_ERROR_GRID_SCALE)
     ENDIF
     wrk1d(1:g(3)%size,3) = z(1:g(3)%size,1) ! we need extra space

! In the Oy direction, we allow to have a different box
     jmax_aux = g(2)%size; subdomain = 0

     dummy = (y_dst(g_dst(2)%size)-y(g(2)%size,1)) / (y(g(2)%size,1)-y(g(2)%size-1,1))
     IF      ( dummy .GT.  C_1EM3_R ) THEN ! Extend
        flag_extend  = .TRUE.
        subdomain(4) = INT(dummy) +1       ! # planes to add at the top
        jmax_aux     = jmax_aux +subdomain(4)
     ELSE IF ( dummy .LT. -C_1EM3_R ) THEN ! Crop
        flag_crop = .TRUE.
        DO j = jmax-1,1,-1
           IF ( ( g(2)%nodes(j)-y_dst(g_dst(2)%size) )*( g(2)%nodes(j+1)-y_dst(g_dst(2)%size) ) .LT. C_0_R ) EXIT
        ENDDO
        subdomain(4) = j +1                ! top plane of cropped region
        jmax_aux     = subdomain(4)
        subdomain(3) = 1
     ENDIF

     dummy = (y_dst(1)-y(1,1)) / (y(2,1)-y(1,1))
     IF      ( dummy .LT. -C_1EM3_R ) THEN ! Extend
        flag_extend = .TRUE.
        subdomain(3) = INT(ABS(dummy)) +1       ! # planes to add at the bottom
        jmax_aux     = jmax_aux +subdomain(3)
     ELSE IF ( dummy .GT.  C_1EM3_R ) THEN ! Crop
        flag_crop = .TRUE.
        DO j = 1,jmax-1,1
           IF ( ( g(2)%nodes(j)-y_dst(1) )*( g(2)%nodes(j+1)-y_dst(1) ) .LT. C_0_R ) EXIT
        ENDDO
        subdomain(3) = j                   ! bottom plane of cropped region
        jmax_aux     = jmax_aux -subdomain(3) +1
     ENDIF

     IF ( flag_extend .AND. flag_crop ) THEN
        CALL IO_WRITE_ASCII(efile, 'TRANSFORM. Simultaneous extend and crop is undeveloped.')
        CALL DNS_STOP(DNS_ERROR_UNDEVELOP)
     ENDIF

! Reallocating memory space
     idummy      = MAX(jmax_aux,MAX(g(1)%size,g(3)%size))
     isize_wrk1d = MAX(isize_wrk1d,idummy)
     isize_wrk1d = isize_wrk1d + 1

     inb_txc         = inb_txc -1    ! Creating txc_aux
     idummy          = MAX(imax,imax_dst) *MAX(jmax_aux,MAX(jmax,jmax_dst)) *MAX(kmax,kmax_dst)
     isize_txc_field = MAX(isize_txc_field,idummy)
#ifdef USE_MPI
     idummy = kmax *jmax_aux
     IF ( MOD(idummy,ims_npro_i) .NE. 0 ) THEN ! add space for MPI transposition
        idummy = idummy      /ims_npro_i
        idummy =(idummy +1 ) *ims_npro_i
     ENDIF
     idummy = idummy *MAX(imax,imax_dst)
     isize_txc_field = MAX(isize_txc_field,idummy)
#endif
     isize_wrk3d     = isize_txc_field

     idummy = isize_wrk1d*7 + (isize_wrk1d+10)*36
     isize_wrk3d = MAX(isize_wrk3d,idummy)

     DEALLOCATE(txc,wrk3d)

     ALLOCATE(txc(isize_txc_field,inb_txc))
     ALLOCATE(wrk3d(isize_wrk3d))
     ALLOCATE(y_aux(isize_wrk1d))
     ALLOCATE(txc_aux(imax,jmax_aux,kmax))

! Creating grid
     IF ( flag_crop ) THEN
        WRITE(str,'(I3)') subdomain(4)
        CALL IO_WRITE_ASCII(lfile, 'Croping above '//TRIM(ADJUSTL(str))//' for remeshing...')
        WRITE(str,'(I3)') subdomain(3)
        CALL IO_WRITE_ASCII(lfile, 'Croping below '//TRIM(ADJUSTL(str))//' for remeshing...')
        CALL TRANS_CROP(i1,jmax,1, subdomain, g(2)%nodes, y_aux)

        y_aux(1)        = y_dst(1)             ! Using min and max of new grid
        y_aux(jmax_aux) = y_dst(g_dst(2)%size)

     ELSE
        y_aux(1+subdomain(3):g(2)%size+subdomain(3)) = y(1:g(2)%size,1) ! we need extra space

        IF ( subdomain(4) .GT. 0 ) THEN
           WRITE(str,'(I3)') subdomain(4)
           CALL IO_WRITE_ASCII(lfile, 'Adding '//TRIM(ADJUSTL(str))//' planes at the top for remeshing...')
           dummy = (y_dst(g_dst(2)%size)-y(g(2)%size,1)) / REAL(subdomain(4)) ! distributing the points uniformly
           DO ip = g(2)%size+subdomain(3)+1,g(2)%size+subdomain(3)+subdomain(4)
              y_aux(ip) = y_aux(ip-1) + dummy
           ENDDO
        ENDIF

        IF ( subdomain(3) .GT. 0 ) THEN
           WRITE(str,'(I3)') subdomain(3)
           CALL IO_WRITE_ASCII(lfile, 'Adding '//TRIM(ADJUSTL(str))//' planes at the bottom for remeshing...')
           dummy = (y_dst(1)-y(1,1)) / REAL(subdomain(3))
           DO ip = subdomain(3),1,-1
              y_aux(ip) = y_aux(ip+1) + dummy ! dummy is negative
           ENDDO
        ENDIF

     ENDIF

     g(2)%scale = g_dst(2)%scale  ! watch out, overwriting grid information
     g(2)%size  = jmax_aux

  ENDIF

! ###################################################################
! Postprocess given list of files
! ###################################################################
  DO it = 1,itime_size
     itime = itime_vec(it)

     WRITE(sRes,*) itime; sRes = 'Processing iteration It'//TRIM(ADJUSTL(sRes))
     CALL IO_WRITE_ASCII(lfile,sRes)

     IF ( iread_flow .EQ. 1 ) THEN ! Flow variables
        WRITE(flow_file,*) itime; flow_file = TRIM(ADJUSTL(tag_flow))//TRIM(ADJUSTL(flow_file))
        CALL DNS_READ_FIELDS(flow_file, i2, imax,jmax,kmax, inb_flow, i0, isize_wrk3d, q, wrk3d)
     ENDIF

     IF ( iread_scal .EQ. 1 ) THEN ! Scalar variables
        WRITE(scal_file,*) itime; scal_file = TRIM(ADJUSTL(tag_scal))//TRIM(ADJUSTL(scal_file))
        CALL DNS_READ_FIELDS(scal_file, i1, imax,jmax,kmax, inb_scal, i0, isize_wrk3d, s, wrk3d)
     ENDIF

! ###################################################################
! Cropping
! ###################################################################
     IF ( opt_main .EQ. 1 ) THEN
        IF ( subdomain(5) .NE. 1 .OR. subdomain(6) .NE. g(3)%size) THEN
           CALL IO_WRITE_ASCII(efile,'TRANSFORM. Cropping only in Oy.')
           CALL DNS_STOP(DNS_ERROR_UNDEVELOP)
        ENDIF
        IF ( subdomain(1) .NE. 1 .OR. subdomain(2) .NE. g(1)%size) THEN
           CALL IO_WRITE_ASCII(efile,'TRANSFORM. Cropping only in Oy.')
           CALL DNS_STOP(DNS_ERROR_UNDEVELOP)
        ENDIF
        IF ( subdomain(3) .LT. 1 .OR. subdomain(4) .GT. g(2)%size) THEN
           CALL IO_WRITE_ASCII(efile,'TRANSFORM. Cropping out of bounds in Oy.')
           CALL DNS_STOP(DNS_ERROR_UNDEVELOP)
        ENDIF

        IF ( icalc_flow .GT. 0 ) THEN
           DO iq = 1,inb_flow
              CALL IO_WRITE_ASCII(lfile,'Transfering data to new array...')
              CALL TRANS_CROP(imax,jmax,kmax, subdomain, q(1,iq), q_dst(1,iq))
           ENDDO
        ENDIF

        IF ( icalc_scal .GT. 0 ) THEN
           DO is = 1,inb_scal
              CALL IO_WRITE_ASCII(lfile,'Transfering data to new array...')
              CALL TRANS_CROP(imax,jmax,kmax, subdomain, s(1,is), s_dst(1,is))
           ENDDO
        ENDIF

! ###################################################################
! Extension
! ###################################################################
     ELSE IF ( opt_main .EQ. 2 ) THEN
        IF ( icalc_flow .GT. 0 ) THEN
           DO iq = 1,inb_flow
              CALL IO_WRITE_ASCII(lfile,'Transfering data to new array...')
              CALL TRANS_EXTEND(imax,jmax,kmax, subdomain, q(1,iq), q_dst(1,iq))
           ENDDO
        ENDIF

        IF ( icalc_scal .GT. 0 ) THEN
           DO is = 1,inb_scal
              CALL IO_WRITE_ASCII(lfile,'Transfering data to new array...')
              CALL TRANS_EXTEND(imax,jmax,kmax, subdomain, s(1,is), s_dst(1,is))
           ENDDO
        ENDIF

! ###################################################################
! Change grid
! ###################################################################
     ELSE IF ( opt_main .EQ. 3 ) THEN

        IF ( icalc_flow .GT. 0 ) THEN
           DO iq = 1,inb_flow
              CALL IO_WRITE_ASCII(lfile,'Transfering data to new array...')
              IF ( flag_crop ) THEN
                 CALL TRANS_CROP(imax,jmax,kmax, subdomain, q(:,iq), txc_aux)
                 DO k = 1,kmax
                    txc_aux(:,1,       k) = txc_aux(:,1,         k) &
                         + (y_aux(1)-y(subdomain(3),1)) *(txc_aux(:,2,k)-txc_aux(:,1,k)) /(y(subdomain(3)+1,1)-y(subdomain(3),1))
                    txc_aux(:,jmax_aux,k) = txc_aux(:,jmax_aux-1,k) &
                         + (y_aux(jmax_aux)-y(subdomain(4)-1,1)) *(txc_aux(:,jmax_aux,k)-txc_aux(:,jmax_aux-1,k)) /(y(subdomain(4),1)-y(subdomain(4)-1,1))
                 ENDDO
              ELSE
                 CALL TRANS_EXTEND(imax,jmax,kmax, subdomain, q(:,iq), txc_aux)
              ENDIF
              CALL OPR_INTERPOLATE(imax,jmax_aux,kmax, imax_dst,jmax_dst,kmax_dst, &
                   g, wrk1d(:,1),y_aux,wrk1d(:,3), x_dst,y_dst,z_dst, &
                   txc_aux,q_dst(:,iq), txc, isize_wrk3d, wrk3d)
           ENDDO
        ENDIF

        IF ( icalc_scal .GT. 0 ) THEN
           DO is = 1,inb_scal
              CALL IO_WRITE_ASCII(lfile,'Transfering data to new array...')
              IF ( flag_crop ) THEN
                 CALL TRANS_CROP(imax,jmax,kmax, subdomain, s(:,is), txc_aux)
                 DO k = 1,kmax
                    txc_aux(:,1,       k) = txc_aux(:,1,         k) &
                         + (y_aux(1)-y(subdomain(3),1)) *(txc_aux(:,2,k)-txc_aux(:,1,k)) /(y(subdomain(3)+1,1)-y(subdomain(3),1))
                    txc_aux(:,jmax_aux,k) = txc_aux(:,jmax_aux-1,k) &
                         + (y_aux(jmax_aux)-y(subdomain(4)-1,1)) *(txc_aux(:,jmax_aux,k)-txc_aux(:,jmax_aux-1,k)) /(y(subdomain(4),1)-y(subdomain(4)-1,1))
                 ENDDO
              ELSE
                 CALL TRANS_EXTEND(imax,jmax,kmax, subdomain, s(:,is), txc_aux)
              ENDIF
              CALL OPR_INTERPOLATE(imax,jmax_aux,kmax, imax_dst,jmax_dst,kmax_dst, &
                   g, wrk1d(:,1),y_aux,wrk1d(:,3), x_dst,y_dst,z_dst, &
                   txc_aux,s_dst(:,is), txc, isize_wrk3d, wrk3d)
           ENDDO
        ENDIF

! ###################################################################
! Linear combination of fields
! ###################################################################
     ELSE IF ( opt_main .EQ. 4 ) THEN
        IF ( icalc_flow .GT. 0 ) THEN
           q_dst = q_dst + q *opt_vec(it+1)
        ENDIF

        IF ( icalc_scal .GT. 0 ) THEN
           s_dst = s_dst + s *opt_vec(it+1)
        ENDIF

! ###################################################################
! Filter
! ###################################################################
     ELSE IF ( opt_main .EQ. 5 ) THEN
        IF ( icalc_flow .GT. 0 ) THEN
           DO iq = 1,inb_flow
              CALL IO_WRITE_ASCII(lfile,'Filtering...')
              q_dst(:,iq) = q(:,iq) ! in-place operation
              IF ( FilterDomain(1)%type .EQ. DNS_FILTER_HELMHOLTZ ) &  ! Bcs depending on field
                   FilterDomain(2)%BcsMin = FilterDomainBcsFlow(iq)
              CALL OPR_FILTER(imax,jmax,kmax, FilterDomain, q_dst(1,iq), wrk1d,wrk2d,txc)
           ENDDO
        ENDIF

        IF ( icalc_scal .GT. 0 ) THEN
           DO is = 1,inb_scal
              CALL IO_WRITE_ASCII(lfile,'Filtering...')
              s_dst(:,is) = s(:,is) ! in-place operation
              IF ( FilterDomain(1)%type .EQ. DNS_FILTER_HELMHOLTZ ) & ! Bcs depending on field
                   FilterDomain(2)%BcsMin = FilterDomainBcsScal(is)
              CALL OPR_FILTER(imax,jmax,kmax, FilterDomain, s_dst(1,is), wrk1d,wrk2d,txc)
           ENDDO
        ENDIF

! ###################################################################
! Transformation
! ###################################################################
     ELSE IF ( opt_main .EQ. 6 ) THEN
        IF      ( opt_function .EQ. 1 ) THEN
           CALL TRANS_FUNCTION(imax,jmax,kmax, s,s_dst, txc)

        ELSE IF ( opt_function .EQ. 2 ) THEN
           s_dst(:,1) = C_0_R
           DO is = 1,MIN(inb_scal,iopt_size)
              s_dst(:,1) = s_dst(:,1) + opt_vec(2+is) *s(:,is)
           ENDDO

        ENDIF

! ###################################################################
! Blend
! ###################################################################
     ELSE IF ( opt_main .EQ. 7 ) THEN
        IF ( it .EQ. 1 ) opt_vec(2) = y(1,1) + opt_vec(2) *g(2)%scale
        WRITE(sRes,*) opt_vec(2),opt_vec(3); sRes = 'Blending with '//TRIM(ADJUSTL(sRes))
        CALL IO_WRITE_ASCII(lfile,sRes)

        IF ( icalc_scal .GT. 0 ) THEN
           DO is = 1,inb_scal
              CALL TRANS_BLEND(imax,jmax,kmax, opt_vec(2),y, s(1,is),s_dst(1,is))
           ENDDO
        ENDIF

        IF ( icalc_flow .GT. 0 ) THEN ! Blended fields have rtime from last velocity field
           DO iq = 1,inb_flow
              CALL TRANS_BLEND(imax,jmax,kmax, opt_vec(2),y, q(1,iq),q_dst(1,iq))
           ENDDO
        ENDIF

        IF ( it .EQ. 1 ) opt_vec(3) = -opt_vec(3) ! flipping blending shape

! ###################################################################
! Adding mean profiles
! ###################################################################
     ELSE IF ( opt_main .EQ. 8 ) THEN
        IF ( icalc_flow .GT. 0 ) THEN
           DO iq = 1,inb_flow
              CALL IO_WRITE_ASCII(lfile,'Adding mean flow profiles...')
              CALL TRANS_ADD_MEAN(i0, iq, imax,jmax,kmax, y, q(1,iq), q_dst(1,iq))
           ENDDO
        ENDIF

        IF ( icalc_scal .GT. 0 ) THEN
           DO is = 1,inb_scal
              CALL IO_WRITE_ASCII(lfile,'Adding mean scal profiles...')
              CALL TRANS_ADD_MEAN(i1, is, imax,jmax,kmax, y, s(1,is), s_dst(1,is))
           ENDDO
        ENDIF

     ENDIF

! ###################################################################
! Writing transform fields
! ###################################################################
     IF ( opt_main .NE. 4 .AND. opt_main .NE. 7 ) THEN
        IF ( icalc_flow .GT. 0 ) THEN
           flow_file=TRIM(ADJUSTL(flow_file))//'.trn'
           CALL DNS_WRITE_FIELDS(flow_file, i2, imax_dst,jmax_dst,kmax_dst, inb_flow, isize_wrk3d, q_dst,wrk3d)
        ENDIF
        IF ( icalc_scal .GT. 0 ) THEN
           scal_file=TRIM(ADJUSTL(scal_file))//'.trn'
           CALL DNS_WRITE_FIELDS(scal_file, i1, imax_dst,jmax_dst,kmax_dst, inb_scal_dst, isize_wrk3d, s_dst,wrk3d)
        ENDIF
     ENDIF

  ENDDO

! ###################################################################
! Final operations
! ###################################################################
  IF ( opt_main .EQ. 4 .OR. opt_main .EQ. 7 ) THEN
     IF ( icalc_flow .GT. 0 ) THEN
        flow_file='flow.trn'
        CALL DNS_WRITE_FIELDS(flow_file, i2, imax_dst,jmax_dst,kmax_dst, inb_flow, isize_wrk3d, q_dst,wrk3d)
     ENDIF
     IF ( icalc_scal .GT. 0 ) THEN
        scal_file='scal.trn'
        CALL DNS_WRITE_FIELDS(scal_file, i1, imax_dst,jmax_dst,kmax_dst, inb_scal_dst, isize_wrk3d, s_dst,wrk3d)
     ENDIF
  ENDIF

  CALL DNS_STOP(0)
END PROGRAM TRANSFIELDS

!########################################################################
!# DESCRIPTION
!#
!# Crop array a into array b in the first two indices
!#
!########################################################################
SUBROUTINE TRANS_CROP(nx,ny,nz, subdomain, a, b)

  IMPLICIT NONE

  TINTEGER nx,ny,nz, subdomain(6)
  TREAL, DIMENSION(nx,ny,nz)                          :: a
  TREAL, DIMENSION(nx,subdomain(4)-subdomain(3)+1,nz) :: b

! -----------------------------------------------------------------------
  TINTEGER j, k

! #######################################################################
  DO k = 1,nz
     DO j = subdomain(3),subdomain(4)
        b(:,j-subdomain(3)+1,k) = a(:,j,k)
     ENDDO
  ENDDO

  RETURN
END SUBROUTINE TRANS_CROP

!########################################################################
!# DESCRIPTION
!#
!# Extend array a into array b in the first two indices
!#
!########################################################################
SUBROUTINE TRANS_EXTEND(nx,ny,nz, planes, a, b)

  IMPLICIT NONE

  TINTEGER nx, ny, nz, planes(6)
  TREAL, DIMENSION(nx,ny,nz)                                         :: a
  TREAL, DIMENSION(planes(1)+nx+planes(2),planes(3)+ny+planes(4),nz) :: b

! -----------------------------------------------------------------------
  TINTEGER j, k

! #######################################################################
  DO k = 1,nz
     b(1+planes(1):nx+planes(1),1+planes(3):ny+planes(3),k) = a(1:nx,1:ny,k)

! extension in i
     DO j = 1,ny
        b(             1:   planes(1)          ,j,k) = b( 1+planes(1),j,k)
        b(nx+planes(1)+1:nx+planes(1)+planes(2),j,k) = b(nx+planes(1),j,k)
     ENDDO

! extension in j; corners are now written
     DO j = 1,planes(3)
        b(:,j,k) = b(:, 1+planes(3),k)
     ENDDO

     DO j = ny+planes(3)+1,ny+planes(3)+planes(4)
        b(:,j,k) = b(:,ny+planes(3),k)
     ENDDO

  ENDDO

  RETURN
END SUBROUTINE TRANS_EXTEND

!########################################################################
!# DESCRIPTION
!#
!########################################################################
SUBROUTINE TRANS_ADD_MEAN(flag_mode, is, nx,ny,nz, y, a,b)

  USE DNS_CONSTANTS, ONLY : efile
  USE DNS_GLOBAL, ONLY : g, sbg, qbg

  IMPLICIT NONE

  TINTEGER flag_mode, is, nx,ny,nz
  TREAL, DIMENSION(*),        INTENT(IN)  :: y
  TREAL, DIMENSION(nx,ny,nz), INTENT(IN)  :: a
  TREAL, DIMENSION(nx,ny,nz), INTENT(OUT) :: b

  ! -----------------------------------------------------------------------
  TINTEGER j
  TREAL PROFILES, ycenter, dummy
  EXTERNAL PROFILES

  ! #######################################################################
  IF ( flag_mode .EQ. 0 ) THEN ! Velocity
    IF ( is .EQ. 1 ) THEN ! Only the mean velocity
      ycenter = y(1) + g(2)%scale *qbg(1)%ymean
      DO j = 1,ny
        dummy =  PROFILES&
        (qbg(1)%type, qbg(1)%thick, qbg(1)%delta, qbg(1)%mean, ycenter, qbg(1)%parameters, y(j))
        b(:,j,:) = dummy + a(:,j,:)
      ENDDO
    ELSE
      b = a
    ENDIF

  ELSE                         ! Scalars
    ycenter = y(1) + g(2)%scale *sbg(is)%ymean
    DO j = 1,ny
      dummy =  PROFILES&
      (sbg(is)%type, sbg(is)%thick, sbg(is)%delta, sbg(is)%mean, ycenter, sbg(is)%parameters, y(j))
      b(:,j,:) = dummy + a(:,j,:)
    ENDDO

  ENDIF

  RETURN
END SUBROUTINE TRANS_ADD_MEAN

!########################################################################
!# DESCRIPTION
!#
!# Calculate b = f(a)
!#
!########################################################################
SUBROUTINE TRANS_FUNCTION(nx,ny,nz, a,b, txc)

  USE DNS_GLOBAL, ONLY : inb_scal, epbackground
  USE THERMO_GLOBAL, ONLY : imixture, MRATIO, GRATIO, dsmooth
  USE THERMO_GLOBAL, ONLY : THERMO_AI, WGHT_INV

  IMPLICIT NONE

  TINTEGER nx,ny,nz
  TREAL, DIMENSION(nx*ny*nz)   :: a, b
  TREAL, DIMENSION(nx*ny*nz,*) :: txc

! -----------------------------------------------------------------------
  TREAL qt_0,qt_1, h_0,h_1, p
  TREAL LATENT_HEAT

! #######################################################################
  imixture = MIXT_TYPE_AIRWATER
  CALL THERMO_INITIALIZE
  MRATIO  = C_1_R
  dsmooth = C_0_R
  inb_scal= 1

  LATENT_HEAT = THERMO_AI(6,1,1)-THERMO_AI(6,1,3)

  qt_0 = 9.0d-3;     qt_1 = 1.5d-3
  h_0  = 0.955376d0; h_1  = 0.981965d0
  p    = 0.940d0

  txc(:,1) = h_0  + a(:)*(h_1 -h_0 ) ! total enthalpy
  txc(:,2) = qt_0 + a(:)*(qt_1-qt_0) ! total water, space for q_l
  txc(:,3) = C_0_R
  txc(:,4) = p                       ! pressure

  CALL THERMO_AIRWATER_PH(nx,ny,nz, txc(1,2), txc(1,1), epbackground,p)        ! Calculate q_l
  CALL THERMO_ANELASTIC_TEMPERATURE(nx,ny,nz, txc(1,1), epbackground, txc(1,5))

! Calculate saturation specific humidity
  CALL THERMO_POLYNOMIAL_PSAT(nx,ny,nz, txc(1,5), txc(1,1))
  txc(:,1) = C_1_R/(MRATIO*txc(:,4)/txc(:,1)-C_1_R) *WGHT_INV(2) /WGHT_INV(1)
  txc(:,1) = txc(:,1)/(C_1_R+txc(:,1))

! Calculate parameter \beta (assuming c_p = c_p,d)
  txc(:,3) = WGHT_INV(2)/WGHT_INV(1)/GRATIO*LATENT_HEAT*LATENT_HEAT / ( txc(:,5)*txc(:,5) )

! Calculate s
  b(:) = txc(:,2) - txc(:,1) * ( C_1_R + txc(:,3)*txc(:,2) ) / ( C_1_R + txc(:,3)*txc(:,1) )

  RETURN
END SUBROUTINE TRANS_FUNCTION

!########################################################################
!# DESCRIPTION
!#
!# b <- b + f(y)*a
!#
!########################################################################
SUBROUTINE TRANS_BLEND(nx,ny,nz, params, y, a, b)

  IMPLICIT NONE

  TINTEGER nx,ny,nz
  TREAL, DIMENSION(*)        :: params
  TREAL, DIMENSION(ny)       :: y
  TREAL, DIMENSION(nx,ny,nz) :: a, b

! -----------------------------------------------------------------------
  TINTEGER j
  TREAL shape, xi

! #######################################################################
  DO j = 1,ny
     xi = ( y(j) - params(1) ) / params(2)
     shape = C_05_R*( C_1_R + TANH(-C_05_R*xi) )
     b(:,j,:) = b(:,j,:) + shape* a(:,j,:)
  ENDDO

RETURN

END SUBROUTINE TRANS_BLEND
