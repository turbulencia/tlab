module TLabMPI_VARS
    use TLab_Constants, only: wp, wi
    implicit none
    save

    integer :: ims_comm_xz                                      ! Plane communicators
    integer :: ims_comm_x, ims_comm_z                           ! Directional communicators
    integer :: ims_comm_xz_aux
    integer :: ims_comm_x_aux, ims_comm_z_aux

    integer :: ims_pro                                          ! Local PE number in global communicator
    integer :: ims_pro_i, ims_pro_j, ims_pro_k                  ! Local PE number in directional communicators
    integer :: ims_npro                                         ! Total number of PE
    integer :: ims_npro_i, ims_npro_j, ims_npro_k               ! Number of PEs in directional communicators
    integer :: ims_offset_i, ims_offset_j, ims_offset_k         ! Grid offsets 

    integer :: ims_err

    real(wp) :: ims_time_min, ims_time_max, ims_time_trans      ! Profiling
    integer(wi) :: ims_bcs_imax, ims_bcs_jmax

    integer(wi), allocatable :: ims_size_i(:), ims_size_k(:)    ! Maybe public in module for TLabMPI_Transpose

end module TLabMPI_VARS
