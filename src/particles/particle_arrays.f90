module PARTICLE_ARRAYS
    use TLab_Types, only: wp, wi
    use TLab_Pointers_3D, only: pointers3d_dt
    use PARTICLE_TYPES
    implicit none
    save

    type(particle_dt) :: l_g                    ! particle tags and Oy-node information in local processor, changes in time
#ifdef USE_MPI
    integer(wi), allocatable :: ims_np_all(:)   ! vector with all # of particles per processor
#endif

    real(wp), allocatable :: l_q(:, :)          ! Lagrangian fields, flow vartiables
    real(wp), allocatable :: l_txc(:, :)        ! Temporary space for Lagrnagian fields

    real(wp), allocatable :: l_work(:)          ! halo space for field-particle interpolations
    real(wp), pointer :: halo_k(:, :, :, :), halo_i(:, :, :, :), halo_ik(:, :, :, :)
    type(pointers3d_dt), allocatable :: p_halo_i(:), p_halo_k(:), p_halo_ik(:)
#ifdef USE_MPI
    real(wp), pointer :: halo_mpi_send_i(:, :, :), halo_mpi_recv_i(:, :, :)
    real(wp), pointer :: halo_mpi_send_k(:, :, :), halo_mpi_recv_k(:, :, :)
    real(wp), pointer :: p_buffer_1(:), p_buffer_2(:)
#endif

    target l_txc, l_work, l_q

end module PARTICLE_ARRAYS
