module PARTICLE_ARRAYS
    use TLAB_TYPES, only: wp
    use PARTICLE_TYPES
    implicit none
    save

    real(wp), allocatable :: l_q(:, :)      ! Lagrangian fields, flow vartiables
    real(wp), allocatable :: l_txc(:, :)    ! Temporary space for Lagrnagian fields

end module PARTICLE_ARRAYS
