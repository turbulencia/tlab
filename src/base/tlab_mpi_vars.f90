module TLabMPI_VARS
    use TLab_Constants, only: wp, wi
    implicit none
    save

    integer :: ims_comm_xz, ims_comm_x, ims_comm_z
    integer :: ims_comm_xz_aux, ims_comm_x_aux, ims_comm_z_aux

    integer :: ims_pro, ims_pro_i, ims_pro_j, ims_pro_k
    integer :: ims_npro, ims_npro_i, ims_npro_j, ims_npro_k
    integer :: ims_offset_i, ims_offset_j, ims_offset_k

    integer :: ims_err
    integer, allocatable :: ims_status(:, :)
    integer, allocatable :: ims_request(:)

    real(wp) :: ims_time_min, ims_time_max, ims_time_trans
    integer(wi) :: ims_bcs_imax, ims_bcs_jmax

    integer(wi), allocatable :: ims_size_i(:), ims_size_k(:)

end module TLabMPI_VARS
