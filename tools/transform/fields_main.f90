#include "types.h"
#include "dns_const.h"
#include "dns_error.h"
#ifdef USE_MPI
#include "dns_const_mpi.h"
#endif

#define C_FILE_LOC "TRANSFORM"

PROGRAM TRANSFIELDS

  USE DNS_TYPES,  ONLY : filter_structure
  USE DNS_CONSTANTS
  USE DNS_GLOBAL
#ifdef USE_MPI
  USE DNS_MPI
#endif

  IMPLICIT NONE

#include "integers.h"
#ifdef USE_MPI
#include "mpif.h"
#endif

! Parameter definitions
  TINTEGER, PARAMETER :: itime_size_max = 512
  TINTEGER, PARAMETER :: iopt_size_max  = 512

! -------------------------------------------------------------------
! Grid and associated arrays
  TREAL, DIMENSION(:,:), ALLOCATABLE, SAVE, TARGET :: x,y,z
  TREAL, DIMENSION(:),   ALLOCATABLE, SAVE         :: x_dst,y_dst,z_dst

! Fields
  TREAL, DIMENSION(:,:), ALLOCATABLE, SAVE, TARGET :: q,     s,    txc
  TREAL, DIMENSION(:,:), ALLOCATABLE, SAVE         :: q_dst, s_dst

! Work arrays
  TREAL, DIMENSION(:,:), ALLOCATABLE, SAVE :: wrk1d, wrk2d
  TREAL, DIMENSION(:),   ALLOCATABLE, SAVE :: wrk3d

! Filters
  TYPE(filter_structure), DIMENSION(3)                              :: FilterDomain
  TREAL,                  DIMENSION(:,:), ALLOCATABLE, SAVE, TARGET :: filter_x, filter_y, filter_z

! -------------------------------------------------------------------
! Local variables
! -------------------------------------------------------------------
  TINTEGER opt_main, opt_filter, opt_function
  TINTEGER iq, is, k, ip, ip_b, ip_t
  TINTEGER isize_wrk3d, idummy, iread_flow, iread_scal, ierr
  CHARACTER*32 inifile, bakfile, flow_file, scal_file
  CHARACTER*64 str, line
  CHARACTER*512 sRes
  TINTEGER subdomain(6)

  TINTEGER imax_dst,jmax_dst,kmax_dst, imax_total_dst,jmax_total_dst,kmax_total_dst, jmax_aux, inb_scal_dst
  TREAL scalex_dst, scaley_dst, scalez_dst
  TREAL dummy

! Filter variables
  TREAL alpha
  TINTEGER delta
  TREAL, DIMENSION(:), POINTER :: p_bcs

  TINTEGER itime_size, it
  TINTEGER itime_vec(itime_size_max)

  TINTEGER iopt_size
  TREAL opt_vec(iopt_size_max)

#ifdef USE_MPI
  TINTEGER npage, id
  INTEGER icount
#endif

! ###################################################################
  inifile = 'dns.ini'
  bakfile = TRIM(ADJUSTL(inifile))//'.bak'

  CALL DNS_INITIALIZE

  CALL DNS_READ_GLOBAL(inifile)

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
  opt_filter   = 0
  opt_function = 0
  
  CALL SCANINICHAR(bakfile, inifile, 'PostProcessing', 'ParamTransform', '-1', sRes)
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
     opt_main = DINT(opt_vec(1))
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
  IF ( opt_main .EQ. 5 ) THEN
     IF ( sRes .EQ. '-1' ) THEN
#ifdef USE_MPI
#else
        WRITE(*,'(A)') 'Filter type ?'
        WRITE(*,'(A)') '1. Compact Fourth-Order'
        WRITE(*,'(A)') '2. Explicit Sixth-Order'
        WRITE(*,'(A)') '3. Explicit Fourth-Order'
        WRITE(*,'(A)') '4. Deconvolution'
        WRITE(*,'(A)') '5. Alpha'
        WRITE(*,'(A)') '6. Spectral cut-off' 
        WRITE(*,'(A)') '7. Gaussian Filter' 
        WRITE(*,'(A)') '8. Tophat' 
        READ(*,*) opt_filter
#endif
     ELSE
        opt_filter = DINT(opt_vec(2))
     ENDIF
     
     alpha = C_0_R
     delta = 2
     IF ( sRes .EQ. '-1' ) THEN
#ifdef USE_MPI
#else
        IF      ( opt_filter .EQ. DNS_FILTER_COMPACT ) THEN
           WRITE(*,*) 'Alpha Coeffcient ?'
           READ(*,*) alpha
        ELSE IF ( opt_filter .EQ. DNS_FILTER_ALPHA   ) THEN
           WRITE(*,*) 'Alpha (multiple of space step) ?'
           READ(*,*) alpha
        ELSE IF ( opt_filter .EQ. DNS_FILTER_CUTOFF  ) THEN
           WRITE(*,*) 'Frequency interval ?'
           READ(*,*) opt_vec(3), opt_vec(4)
        ELSE IF ( opt_filter .EQ. DNS_FILTER_ERF     ) THEN 
           WRITE(*,*) 'Transition wavenumber (in physical units)' 
           WRITE(*,*) '   >0: high-pass filter' 
           WRITE(*,*) '   <0: low-pass  filter'
           WRITE(*,*) 'and Characteristic Width--in log units (relative to domain size)?'   
           WRITE(*,*) '   positive real number.'
           READ(*,*) opt_vec(3), opt_vec(4) 
        ELSE IF ( opt_filter .EQ. DNS_FILTER_TOPHAT ) THEN
           WRITE(*,*) 'Delta (in multiples of grid step)?'
           READ(*,*) delta
        ENDIF
#endif
     ELSE
        alpha = opt_vec(3)
        delta = DINT(opt_vec(3))
     ENDIF
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
        opt_function = DINT(opt_vec(2))
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
  CALL SCANINICHAR(bakfile, inifile, 'PostProcessing', 'Subdomain', '-1', sRes)

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
     subdomain(1) = 1; subdomain(2) = imax_total
     subdomain(3) = 1; subdomain(4) = jmax_total
     subdomain(5) = 1; subdomain(6) = kmax_total
  ENDIF
  
! -------------------------------------------------------------------
  IF      ( opt_main .EQ. 1 .OR. &      ! Crop
            opt_main .EQ. 3      ) THEN ! Remesh
     imax_total_dst = subdomain(2)-subdomain(1)+1
     jmax_total_dst = subdomain(4)-subdomain(3)+1
     kmax_total_dst = subdomain(6)-subdomain(5)+1

  ELSE IF ( opt_main .EQ. 2      ) THEN ! Extend 
     imax_total_dst = imax_total + subdomain(2) + subdomain(1)
     jmax_total_dst = jmax_total + subdomain(4) + subdomain(3)
     kmax_total_dst = kmax_total

  ELSE
     imax_total_dst = imax_total
     jmax_total_dst = jmax_total
     kmax_total_dst = kmax_total
  ENDIF

#ifdef USE_MPI
  imax_dst = imax_total_dst/ims_npro_i
  jmax_dst = jmax_total_dst
  kmax_dst = kmax_total_dst/ims_npro_k
#else
  imax_dst = imax_total_dst
  jmax_dst = jmax_total_dst
  kmax_dst = kmax_total_dst
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
     inb_txc         = 3
  ELSE IF ( opt_main .EQ. 6 ) THEN
     inb_txc         = 5
     inb_scal_dst    = 1
  ENDIF

  IF ( opt_filter .EQ. DNS_FILTER_ALPHA ) inb_txc = MAX(inb_txc,4) ! needs more space

  IF ( ifourier .EQ. 1 ) inb_txc = MAX(inb_txc,1)
 
! -------------------------------------------------------------------
  isize_wrk3d = MAX(isize_txc_field,imax_dst*jmax_dst*kmax_dst)

  IF ( opt_main .EQ. 3 ) THEN ! space for splines interpolation
     idummy      = MAX(imax_total_dst,MAX(jmax_total_dst,kmax_total_dst))
     isize_wrk1d = MAX(isize_wrk1d,idummy)
     isize_wrk1d = isize_wrk1d + 1

     idummy = isize_wrk1d*7 + (isize_wrk1d+10)*36
     isize_wrk3d = MAX(isize_wrk3d,idummy)
  ENDIF

! -------------------------------------------------------------------
  IF ( icalc_flow .EQ. 1 ) ALLOCATE(q_dst(imax_dst*jmax_dst*kmax_dst,inb_flow))
  IF ( icalc_scal .EQ. 1 ) ALLOCATE(s_dst(imax_dst*jmax_dst*kmax_dst,inb_scal_dst))

  ALLOCATE(wrk1d(isize_wrk1d,inb_wrk1d))
  ALLOCATE(wrk2d(isize_wrk2d,inb_wrk2d))

  IF ( opt_main .EQ. 3 ) THEN
     ALLOCATE(x_dst(imax_total_dst))
     ALLOCATE(y_dst(jmax_total_dst))
     ALLOCATE(z_dst(kmax_total_dst))
  ENDIF

#include "dns_alloc_arrays.h"

! -------------------------------------------------------------------
! Read the grid 
! -------------------------------------------------------------------
#include "dns_read_grid.h"

! -------------------------------------------------------------------
! Initialize filters
! -------------------------------------------------------------------
  IF ( opt_main .EQ. 5 ) THEN ! Filters
     FilterDomain(:)%type       = opt_filter
     FilterDomain(:)%delta      = delta
     FilterDomain(:)%alpha      = alpha
     FilterDomain(:)%size       = g(:)%size
     FilterDomain(:)%periodic   = g(:)%periodic
     FilterDomain(:)%uniform    = g(:)%uniform
     FilterDomain(:)%inb_filter = FilterDomain(:)%delta +1 ! Default
     FilterDomain(:)%bcs_min    = 0
     FilterDomain(:)%bcs_max    = 0
     
     DO is = 1,3
        IF ( FilterDomain(is)%size .EQ. 1 ) FilterDomain(is)%type = DNS_FILTER_NONE
        
        SELECT CASE( FilterDomain(is)%type )
           
        CASE( DNS_FILTER_4E, DNS_FILTER_ADM )
           FilterDomain(is)%inb_filter = 5
           
        CASE( DNS_FILTER_COMPACT )
           FilterDomain(is)%inb_filter = 6
        
        CASE( DNS_FILTER_TOPHAT )
           IF ( MOD(FilterDomain(is)%delta,2) .NE. 0 ) THEN
              CALL IO_WRITE_ASCII(efile, 'DNS_MAIN. Filter delta is not even.')
              CALL DNS_STOP(DNS_ERROR_PARAMETER)
           ENDIF
           
        END SELECT
     ENDDO
     
     ALLOCATE( filter_x( FilterDomain(1)%size, FilterDomain(1)%inb_filter))
     FilterDomain(1)%coeffs => filter_x
     ALLOCATE( filter_y( FilterDomain(2)%size, FilterDomain(2)%inb_filter))
     FilterDomain(2)%coeffs => filter_y
     ALLOCATE( filter_z( FilterDomain(3)%size, FilterDomain(3)%inb_filter))
     FilterDomain(3)%coeffs => filter_z
     
     DO is = 1,3     
        SELECT CASE( FilterDomain(is)%type )
           
        CASE( DNS_FILTER_4E, DNS_FILTER_ADM )
           CALL FLT_E4_INI(g(is)%scale, g(is)%nodes, FilterDomain(is))
           
        CASE( DNS_FILTER_COMPACT )
           CALL FLT_C4_INI(             g(is)%jac,   FilterDomain(is))
        
        CASE( DNS_FILTER_TOPHAT )
           CALL FLT_T1_INI(g(is)%scale, g(is)%nodes, FilterDomain(is), wrk1d)
           
        END SELECT
     END DO
     
#ifdef USE_MPI
     FilterDomain(1)%mpitype = DNS_MPI_I_PARTIAL 
     FilterDomain(3)%mpitype = DNS_MPI_K_PARTIAL
#endif

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
        IF ( subdomain(5) .NE. 1 .OR. subdomain(6) .NE. kmax_total) THEN
           CALL IO_WRITE_ASCII(efile,'TRANSFORM. Cropping only in Oy.')
           CALL DNS_STOP(DNS_ERROR_UNDEVELOP)           
        ENDIF
        IF ( subdomain(1) .NE. 1 .OR. subdomain(2) .NE. imax_total) THEN
           CALL IO_WRITE_ASCII(efile,'TRANSFORM. Cropping only in Oy.')
           CALL DNS_STOP(DNS_ERROR_UNDEVELOP)           
        ENDIF
        IF ( subdomain(3) .LT. 1 .OR. subdomain(4) .GT. jmax_total) THEN
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
        CALL IO_READ_GRID('grid.trn', imax_total_dst,jmax_total_dst,kmax_total_dst, &
             scalex_dst,scaley_dst,scalez_dst, x_dst,y_dst,z_dst)

! Check grids. In the Oy direction, we allow to have a larger box
        jmax_aux = jmax_total; subdomain = 0

! To be updated: we need to distinguish between periodic and nonperiodic cases.
!        dummy = (x_dst(imax_total_dst)-x(imax_total)) / (x(imax_total)-x(imax_total-1))
        dummy = (scalex_dst-scalex) / (x(imax_total,1)-x(imax_total-1,1))
        IF ( ABS(dummy) .GT. C_1EM3_R ) THEN
           CALL IO_WRITE_ASCII(efile, 'TRANSFORM. Ox scales are not equal at the end.')
           CALL DNS_STOP(DNS_ERROR_GRID_SCALE)
        ENDIF
        wrk1d(1:imax_total,1) = x(1:imax_total,1) ! we need extra space
        
        dummy = (y_dst(jmax_total_dst)-y(jmax_total,1)) / (y(jmax_total,1)-y(jmax_total-1,1))
!        dummy = (scaley_dst-scaley) / (y(jmax_total)-y(jmax_total-1))
        IF ( ABS(dummy) .GT. C_1EM3_R ) THEN
           IF ( dummy .GT. C_0_R ) THEN
              subdomain(4) = ABS(jmax_total_dst - jmax_total) ! additional planes at the top
              jmax_aux = jmax_aux + subdomain(4)
              dummy = (y_dst(jmax_total_dst)-y(jmax_total,1)) / INT(subdomain(4))
!              dummy = (scaley_dst-scaley) / INT(subdomain(4))
           ELSE
              CALL IO_WRITE_ASCII(efile, 'TRANSFORM. Oy scales are not equal at the end.')
              CALL DNS_STOP(DNS_ERROR_GRID_SCALE)
           ENDIF
        ENDIF
        wrk1d(1+subdomain(3):jmax_total+subdomain(3),2) = y(1:jmax_total,1) ! we need extra space
        DO ip = jmax_total+1,jmax_aux
           wrk1d(ip,2) = wrk1d(ip-1,2) + dummy
        ENDDO

        dummy = (y_dst(1)-y(1,1)) / (y(2,1)-y(1,1))
        IF ( ABS(dummy) .GT. C_1EM3_R ) THEN
           IF ( dummy .LT. C_0_R ) THEN
              subdomain(3) = ABS(jmax_total_dst - jmax_total) ! additional planes at the beginning
              jmax_aux = jmax_aux + subdomain(3)
              dummy = (y_dst(1)-y(1,1)) / INT(subdomain(3))
           ELSE
              CALL IO_WRITE_ASCII(efile, 'TRANSFORM. Oy scales are not equal at the beginning.')
              CALL DNS_STOP(DNS_ERROR_GRID_SCALE)
           ENDIF
        ENDIF
        wrk1d(1+subdomain(3):jmax_total+subdomain(3),2) = y(1:jmax_total,1) ! we need extra space
        DO ip = subdomain(3),1,-1
           wrk1d(ip,2) = wrk1d(ip+1,2) + dummy ! dummy is negative
        ENDDO
           
!        dummy = (z_dst(kmax_total_dst)-z(kmax_total)) / (z(kmax_total)-z(kmax_total-1))
        dummy = (scalez_dst-scalez) / (z(kmax_total,1)-z(kmax_total-1,1))
        IF ( ABS(dummy) .GT. C_1EM3_R ) THEN
           CALL IO_WRITE_ASCII(efile, 'TRANSFORM. Oz scales are not equal')
           CALL DNS_STOP(DNS_ERROR_GRID_SCALE)
        ENDIF
        wrk1d(1:kmax_total,3) = z(1:kmax_total,1) ! we need extra space

! Remeshing the fields
        IF ( icalc_flow .GT. 0 ) THEN
           DO iq = 1,inb_flow
              CALL IO_WRITE_ASCII(lfile,'Transfering data to new array...')
              CALL TRANS_EXTEND(imax,jmax,kmax, subdomain, q(:,iq), txc(:,1))
              CALL OPR_INTERPOLATE(imax,jmax_aux,kmax, imax_dst,jmax_dst,kmax_dst, i1bc,j1bc,k1bc, &
                   scalex,scaley_dst,scalez, wrk1d(:,1),wrk1d(:,2),wrk1d(:,3), x_dst,y_dst,z_dst, &
                   txc(:,1),q_dst(:,iq), txc(:,2), isize_wrk3d, wrk3d)
           ENDDO
        ENDIF
              
        IF ( icalc_scal .GT. 0 ) THEN
           DO is = 1,inb_scal
              CALL IO_WRITE_ASCII(lfile,'Transfering data to new array...')
              CALL TRANS_EXTEND(imax,jmax,kmax, subdomain, s(:,is), txc(:,1))
              CALL OPR_INTERPOLATE(imax,jmax_aux,kmax, imax_dst,jmax_dst,kmax_dst, i1bc,j1bc,k1bc, &
                   scalex,scaley_dst,scalez, wrk1d(:,1),wrk1d(:,2),wrk1d(:,3), x_dst,y_dst,z_dst, &
                   txc(:,1),s_dst(:,is), txc(:,2), isize_wrk3d, wrk3d)
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
        IF ( opt_filter .EQ. DNS_FILTER_ALPHA  ) dummy =-C_1_R /(alpha*g(1)%jac(1,1))**2
        IF ( opt_filter .EQ. DNS_FILTER_CUTOFF ) dummy = C_1_R /M_REAL(imax_total*kmax_total)
        IF ( opt_filter .EQ. DNS_FILTER_ERF    ) dummy = C_1_R /M_REAL(imax_total*kmax_total)

        IF ( icalc_flow .GT. 0 ) THEN
           DO iq = 1,inb_flow
              CALL IO_WRITE_ASCII(lfile,'Filtering...')
              
              IF      ( opt_filter .EQ. DNS_FILTER_CUTOFF .OR. opt_filter .EQ. DNS_FILTER_ERF ) THEN
                 txc(1:isize_field,1) = q(1:isize_field,iq)
                 CALL OPR_FOURIER_F(i2, imax,jmax,kmax, txc(:,1),txc(:,2), txc(:,3),wrk2d,wrk3d) 
                 IF      ( opt_filter .EQ. DNS_FILTER_CUTOFF ) THEN 
                    CALL TRANS_CUTOFF_2D(imax,jmax,kmax, opt_vec(3), txc(:,2)) 
                 ELSE IF ( opt_filter .EQ. DNS_FILTER_ERF    ) THEN 
                    CALL TRANS_ERF_2D(imax,jmax,kmax,opt_vec(3),txc(:,2))  
                 ENDIF
                 CALL OPR_FOURIER_B(i2, imax,jmax,kmax, txc(:,2), txc(:,1), wrk3d)
                 q_dst(1:isize_field,iq) = txc(1:isize_field,1) *dummy

              ELSE IF ( opt_filter .EQ. DNS_FILTER_ALPHA  ) THEN
                 ip   =                 1
                 ip_b =                 1
                 ip_t = imax*(jmax-1) + 1
                 DO k = 1,kmax
                    p_bcs => q(ip_b:,iq); wrk2d(ip:ip+imax-1,1) = p_bcs(1:imax); ip_b = ip_b + imax*jmax ! bottom
                    p_bcs => q(ip_t:,iq); wrk2d(ip:ip+imax-1,2) = p_bcs(1:imax); ip_t = ip_t + imax*jmax ! top
                    ip = ip + imax
                 ENDDO
                 q_dst(:,iq) = q(:,iq)*dummy ! forcing term
                 CALL OPR_HELMHOLTZ_FXZ(imax,jmax,kmax, g, i0, dummy,&
                      q_dst(1,iq), txc(1,2),txc(1,3), &
                      wrk2d(1,1),wrk2d(1,2), wrk1d,wrk1d(1,5),wrk3d)
                 
              ELSE ! Rest; ADM needs two arrays in txc, rest just 1
                 q_dst(:,iq) = q(:,iq) 
                 CALL OPR_FILTER(imax,jmax,kmax, FilterDomain, q_dst(1,iq), wrk1d,wrk2d,txc)
              ENDIF
           ENDDO
           
        ENDIF

        IF ( icalc_scal .GT. 0 ) THEN
           DO is = 1,inb_scal
              CALL IO_WRITE_ASCII(lfile,'Filtering...')
              
              IF      ( opt_filter .EQ. DNS_FILTER_CUTOFF .OR. opt_filter .EQ. DNS_FILTER_ERF ) THEN
                 txc(1:isize_field,1) = s(1:isize_field,is)
                 CALL OPR_FOURIER_F(i2, imax,jmax,kmax, txc(:,1),txc(:,2), txc(:,3),wrk2d,wrk3d) 
                 IF      ( opt_filter .EQ. DNS_FILTER_CUTOFF ) THEN 
                    CALL TRANS_CUTOFF_2D(imax,jmax,kmax, opt_vec(3), txc(:,2)) 
                 ELSE IF ( opt_filter .EQ. DNS_FILTER_ERF    ) THEN 
                    CALL TRANS_ERF_2D(imax,jmax,kmax, opt_vec(3), txc(:,2)) 
                 ENDIF
                 CALL OPR_FOURIER_B(i2, imax,jmax,kmax, txc(1,2), txc(1,1), wrk3d)
                 s_dst(1:isize_field,is) = txc(1:isize_field,1) *dummy

              ELSE IF ( opt_filter .EQ. DNS_FILTER_ALPHA ) THEN
                 ! Not yet included

              ELSE ! Rest; ADM needs two arrays in txc, rest just 1
                 s_dst(:,is) = s(:,is) 
                 CALL OPR_FILTER(imax,jmax,kmax, FilterDomain, s_dst(1,is), wrk1d,wrk2d,txc)
                 
              ENDIF
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
        IF ( it .EQ. 1 ) opt_vec(2) = y(1,1) + opt_vec(2)*scaley
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

  CALL DNS_END(0)

  STOP

END PROGRAM TRANSFIELDS

!########################################################################
!# Tool/Library
!#
!########################################################################
!# HISTORY
!#
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
!# Tool/Library
!#
!########################################################################
!# HISTORY
!#
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
!# Tool/Library
!#
!########################################################################
!# HISTORY
!#
!########################################################################
!# DESCRIPTION
!#
!########################################################################
SUBROUTINE TRANS_ADD_MEAN(flag_mode, is, nx,ny,nz, y, a,b)

  USE DNS_CONSTANTS, ONLY : efile
  USE DNS_GLOBAL, ONLY : imode_flow, scaley
  USE DNS_GLOBAL, ONLY : iprof_u, thick_u, delta_u, mean_u, prof_u, ycoor_u
  USE DNS_GLOBAL, ONLY : sbg

  IMPLICIT NONE

  TINTEGER flag_mode, is, nx,ny,nz
  TREAL, DIMENSION(*),        INTENT(IN)  :: y
  TREAL, DIMENSION(nx,ny,nz), INTENT(IN)  :: a
  TREAL, DIMENSION(nx,ny,nz), INTENT(OUT) :: b
  
! -----------------------------------------------------------------------
  TINTEGER j
  TREAL FLOW_SHEAR_TEMPORAL, ycenter, dummy
  EXTERNAL FLOW_SHEAR_TEMPORAL

! #######################################################################
  IF ( imode_flow .EQ. DNS_FLOW_SHEAR     ) THEN
     IF ( flag_mode .EQ. 0 ) THEN ! Velocity
        IF ( is .EQ. 1 ) THEN ! Only the mean velocity
        ycenter = y(1) + scaley*ycoor_u
        DO j = 1,ny
           dummy =  FLOW_SHEAR_TEMPORAL&
                (iprof_u, thick_u, delta_u, mean_u, ycenter, prof_u, y(j))
            b(:,j,:) = dummy + a(:,j,:)
        ENDDO
        ELSE
           b = a
        ENDIF

     ELSE                         ! Scalars
        ycenter = y(1) + scaley *sbg(is)%ymean
        DO j = 1,ny
           dummy =  FLOW_SHEAR_TEMPORAL&
                (sbg(is)%type, sbg(is)%thick, sbg(is)%delta, sbg(is)%mean, ycenter, sbg(is)%parameters, y(j))
           b(:,j,:) = dummy + a(:,j,:)
        ENDDO

     ENDIF
     
  ELSE 
     CALL IO_WRITE_ASCII(efile, 'TRANS_ADD_MEAN. Only shear mode developed.')
     CALL DNS_STOP(DNS_ERROR_UNDEVELOP)
     
  ENDIF

  RETURN
END SUBROUTINE TRANS_ADD_MEAN

!########################################################################
!# Tool/Library
!#
!########################################################################
!# HISTORY
!#
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

  CALL THERMO_AIRWATER_PH(nx,ny,nz, txc(1,2), txc(1,1), epbackground,p)
  CALL THERMO_ANELASTIC_TEMPERATURE(nx,ny,nz, txc(1,1), epbackground, txc(1,5))

! Calculate liquid water temperature from temperature (assuming c_p = c_p,d)
  txc(:,5) = txc(:,5) - LATENT_HEAT*txc(:,3)

! Calculate saturation pressure based on the liquid water temperature
  CALL THERMO_POLYNOMIAL_PSAT(nx,ny,nz, txc(1,5), txc(1,1))

! Calculate saturation specific humidity
  txc(:,1) = C_1_R/(MRATIO*txc(:,4)/txc(:,1)-C_1_R)*WGHT_INV(2)/WGHT_INV(1)
  txc(:,1) = txc(:,1)/(C_1_R+txc(:,1))
 
! Calculate parameter \beta (assuming c_p = c_p,d)
  txc(:,3) = WGHT_INV(2)/WGHT_INV(1)/GRATIO*LATENT_HEAT*LATENT_HEAT / ( txc(:,5)*txc(:,5) )

! Calculate s
  b(:) = txc(:,2) - txc(:,1) * ( C_1_R + txc(:,3)*txc(:,2) ) / ( C_1_R + txc(:,3)*txc(:,1) )

  RETURN
END SUBROUTINE TRANS_FUNCTION

!########################################################################
!# Tool/Library
!#
!########################################################################
!# HISTORY
!#
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

!########################################################################
!# Tool/Library
!#
!########################################################################
!# HISTORY
!#
!########################################################################
!# DESCRIPTION
!#
!# b <- b + f(y)*a
!#
!########################################################################
SUBROUTINE TRANS_CUTOFF_2D(nx,ny,nz, spc_param, a)

  USE DNS_GLOBAL, ONLY : kmax_total, isize_txc_dimz
  USE DNS_GLOBAL, ONLY : scalex, scalez
#ifdef USE_MPI
  USE DNS_MPI,    ONLY : ims_offset_i, ims_offset_k
#endif

  IMPLICIT NONE

  TINTEGER nx,ny,nz
  TREAL,    DIMENSION(*)                                  :: spc_param
  TCOMPLEX, DIMENSION(isize_txc_dimz/2,nz), INTENT(INOUT) :: a
  
! -----------------------------------------------------------------------
  TINTEGER i,j,k, iglobal,kglobal, ip
  TREAL fi,fk,f

! #######################################################################
  DO k = 1,nz
#ifdef USE_MPI
     kglobal = k + ims_offset_k
#else
     kglobal = k
#endif
     IF ( kglobal .LE. kmax_total/2+1 ) THEN; fk = M_REAL(kglobal-1)/scalez
     ELSE;                                    fk =-M_REAL(kmax_total+1-kglobal)/scalez; ENDIF

     DO i = 1,nx/2+1
#ifdef USE_MPI
        iglobal = i + ims_offset_i/2
#else
        iglobal = i
#endif
        fi = M_REAL(iglobal-1)/scalex
        
        f = SQRT(fi**2 + fk**2)
        
! apply spectral cutoff
        DO j = 1,ny
           ip = (j-1)*(nx/2+1) + i 
           IF ( (f-spc_param(1))*(spc_param(2)-f) .LT. C_0_R ) a(ip,k) = C_0_R
        ENDDO
        
     ENDDO
  ENDDO

  RETURN
END SUBROUTINE TRANS_CUTOFF_2D

!########################################################################
!# Tool/Library
!#
!########################################################################
!# HISTORY
!#
!########################################################################
!# DESCRIPTION
!#
!# Spectral filter with smooth (error-function) transition. 
!# The error function filter-response function is imposed 
!# in logarithmic wavenumber space. 
!# 
!#
!########################################################################
!# ARGUMENTS 
!#
!#    spc_param(1) physical frequency for transition 
!#                 > 0: High-pass 
!#                 < 0: Low-pass
!#    spc_param(2) width of transition in logarithmic wavenumber space
!# 
!########################################################################
SUBROUTINE TRANS_ERF_2D(nx,ny,nz, spc_param, a)

  USE DNS_GLOBAL, ONLY : kmax_total, isize_txc_dimz
  USE DNS_GLOBAL, ONLY : scalex, scalez
#ifdef USE_MPI
  USE DNS_MPI,    ONLY : ims_offset_i, ims_offset_k
#endif

  IMPLICIT NONE

  TINTEGER nx,ny,nz
  TREAL,    DIMENSION(*)                                  :: spc_param
  TCOMPLEX, DIMENSION(isize_txc_dimz/2,nz), INTENT(INOUT) :: a
  
! -----------------------------------------------------------------------
  TINTEGER i,j,k, iglobal,kglobal, ip 
  TINTEGER sign_pass, off_pass 
  TREAL fi,fk,f,fcut_log,damp

! ####################################################################### 
  IF ( spc_param(1) .GT. 0 ) THEN; 
     sign_pass = 1.   ! HIGHPASS
     off_pass  = 0.
  ELSE                ! spc_param(1) <= 0 
     sign_pass =-1.   ! LOWPASS  
     off_pass  = 1.
  ENDIF

  fcut_log = LOG(spc_param(1)) 
  DO k = 1,nz
#ifdef USE_MPI
     kglobal = k + ims_offset_k
#else
     kglobal = k
#endif
     IF ( kglobal .LE. kmax_total/2+1 ) THEN; fk = M_REAL(kglobal-1)/scalez
     ELSE;                                    fk =-M_REAL(kmax_total+1-kglobal)/scalez; ENDIF

     DO i = 1,nx/2+1
#ifdef USE_MPI
        iglobal = i + ims_offset_i/2
#else
        iglobal = i
#endif
        fi = M_REAL(iglobal-1)/scalex
        
        f = SQRT(fi**2 + fk**2) 
        IF ( f .GT. 0 )  THEN;  damp = (ERF((LOG(f) - fcut_log)/spc_param(2)) + 1.)/2.
        ELSE ;                  damp = C_0_R; 
        ENDIF

        ! Set to high- or low-pass
        ! high-pass: damp = 0.0 + damp 
        ! low-pass:  damp = 1.0 - damp 
        damp = off_pass + sign_pass*damp  

        ! apply filter
        DO j = 1,ny
           ip = (j-1)*(nx/2+1) + i 
           a(ip,k) = damp*a(ip,k)
        ENDDO
     ENDDO
  ENDDO

  RETURN
END SUBROUTINE TRANS_ERF_2D
