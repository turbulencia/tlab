#include "types.h"
#include "dns_const.h"
#include "dns_const_mpi.h"

module TLAB_MPI_VARS
    implicit none
    save

    integer :: ims_comm_xz, ims_comm_x, ims_comm_z
    integer :: ims_comm_xz_aux, ims_comm_x_aux, ims_comm_z_aux

    integer :: ims_pro, ims_pro_i, ims_pro_j, ims_pro_k
    integer :: ims_npro, ims_npro_i, ims_npro_j, ims_npro_k
    integer :: ims_offset_i, ims_offset_j, ims_offset_k
    integer :: ims_err
    TREAL :: ims_time_min, ims_time_max, ims_time_trans
    TINTEGER :: ims_bcs_imax, ims_bcs_jmax

    integer, dimension(:), allocatable :: ims_map_i
    TINTEGER :: ims_sizBlock_i
    TINTEGER, dimension(:), allocatable :: ims_size_i
    TINTEGER, dimension(:, :), allocatable :: ims_ds_i, ims_dr_i
    integer, dimension(:, :), allocatable :: ims_ts_i, ims_tr_i
    TINTEGER, dimension(:), allocatable :: ims_plan_trps_i, ims_plan_trpr_i
    TINTEGER :: ims_trp_mode_i

!  INTEGER,  DIMENSION(:  ), ALLOCATABLE :: ims_map_j
!  TINTEGER, DIMENSION(  :), ALLOCATABLE :: ims_size_j
!  TINTEGER, DIMENSION(:,:), ALLOCATABLE :: ims_ds_j, ims_dr_j
!  INTEGER,  DIMENSION(:,:), ALLOCATABLE :: ims_ts_j, ims_tr_j

    integer, dimension(:), allocatable :: ims_map_k
    TINTEGER :: ims_sizBlock_k
    TINTEGER, dimension(:), allocatable :: ims_size_k
    TINTEGER, dimension(:, :), allocatable :: ims_ds_k, ims_dr_k
    integer, dimension(:, :), allocatable :: ims_ts_k, ims_tr_k
    TINTEGER, dimension(:), allocatable :: ims_plan_trps_k, ims_plan_trpr_k
    TINTEGER :: ims_trp_mode_k

    TINTEGER, dimension(:), allocatable :: ims_size_p ! Particle data

    integer, dimension(:, :), allocatable :: ims_status
    integer, dimension(:), allocatable :: ims_request

#ifdef USE_PSFFT
    integer :: ims_nb_thrsupp_provided
    integer, dimension(2) :: ims_nb_proc_grid
    integer, dimension(3) :: ims_nb_msize
    integer, dimension(3) :: ims_nb_xsrt, ims_nb_xend, ims_nb_xsiz
    integer, dimension(3) :: ims_nb_ysrt, ims_nb_yend, ims_nb_ysiz
    integer, dimension(3) :: ims_nb_zsrt, ims_nb_zend, ims_nb_zsiz
#endif

end module TLAB_MPI_VARS
