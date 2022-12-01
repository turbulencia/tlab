module PARTICLE_ARRAYS
    use TLAB_TYPES, only: wp, wi
    use PARTICLE_TYPES
    implicit none
    save

    type(particle_dt) :: l_g                    ! particle tags and Oy-node information in local processor, changes in time
#ifdef USE_MPI
    integer(wi), allocatable :: ims_np_all(:)   ! vector with all # of particles per processor
#endif

    real(wp), allocatable :: l_q(:, :)          ! Lagrangian fields, flow vartiables
    real(wp), allocatable :: l_txc(:, :)        ! Temporary space for Lagrnagian fields
    real(wp), allocatable :: l_comm(:)          ! halo space for field-particle interpolations 

    target l_txc, l_comm, l_q
    
end module PARTICLE_ARRAYS
