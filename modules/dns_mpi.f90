#include "types.h"
#include "dns_const_mpi.h"

MODULE DNS_MPI
  USE DNS_TYPES, ONLY : subarray_dt
  IMPLICIT NONE
  SAVE

  INTEGER  :: ims_comm_xz,     ims_comm_x,     ims_comm_z
  INTEGER  :: ims_comm_xz_aux, ims_comm_x_aux, ims_comm_z_aux

  INTEGER  :: ims_pro,  ims_pro_i,  ims_pro_j,  ims_pro_k
  INTEGER  :: ims_npro, ims_npro_i, ims_npro_j, ims_npro_k
  INTEGER  :: ims_offset_i, ims_offset_j, ims_offset_k
  INTEGER  :: ims_err, ims_tag
  TREAL    :: ims_time_min, ims_time_max, ims_time_trans
  TINTEGER :: ims_bcs_imax, ims_bcs_jmax

  INTEGER,  DIMENSION(:  ), ALLOCATABLE :: ims_map_i
  TINTEGER, DIMENSION(  :), ALLOCATABLE :: ims_size_i
  TINTEGER, DIMENSION(:,:), ALLOCATABLE :: ims_ds_i, ims_dr_i
  INTEGER,  DIMENSION(:,:), ALLOCATABLE :: ims_ts_i, ims_tr_i
  INTEGER(KIND=4), DIMENSION(:), ALLOCATABLE :: ims_plan_trps_i, ims_plan_trpr_i 


!  INTEGER,  DIMENSION(:  ), ALLOCATABLE :: ims_map_j
!  TINTEGER, DIMENSION(  :), ALLOCATABLE :: ims_size_j
!  TINTEGER, DIMENSION(:,:), ALLOCATABLE :: ims_ds_j, ims_dr_j
!  INTEGER,  DIMENSION(:,:), ALLOCATABLE :: ims_ts_j, ims_tr_j

  INTEGER,  DIMENSION(:  ), ALLOCATABLE :: ims_map_k
  TINTEGER, DIMENSION(  :), ALLOCATABLE :: ims_size_k
  TINTEGER, DIMENSION(:,:), ALLOCATABLE :: ims_ds_k, ims_dr_k
  INTEGER,  DIMENSION(:,:), ALLOCATABLE :: ims_ts_k, ims_tr_k
  INTEGER(KIND=4), DIMENSION(:), ALLOCATABLE :: ims_plan_trps_k, ims_plan_trpr_k

  TINTEGER, DIMENSION(:),   ALLOCATABLE :: ims_size_p ! Particle data  

  TYPE(subarray_dt), DIMENSION(MPIO_SUBARRAY_SIZE) :: mpio_aux

#ifdef USE_PSFFT  
  INTEGER :: ims_nb_thrsupp_provided 
  INTEGER,DIMENSION(2) :: ims_nb_proc_grid
  INTEGER,DIMENSION(3) :: ims_nb_msize
  INTEGER,DIMENSION(3) :: ims_nb_xsrt,ims_nb_xend,ims_nb_xsiz
  INTEGER,DIMENSION(3) :: ims_nb_ysrt,ims_nb_yend,ims_nb_ysiz
  INTEGER,DIMENSION(3) :: ims_nb_zsrt,ims_nb_zend,ims_nb_zsiz
#endif 

END MODULE DNS_MPI
