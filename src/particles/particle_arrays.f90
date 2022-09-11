module PARTICLE_ARRAYS
    use TLAB_TYPES, only: cp
    use PARTICLE_TYPES
    implicit none
    save

    real(cp), allocatable :: l_q(:, :)      ! Lagrangian fields, flow vartiables
    real(cp), allocatable :: l_txc(:, :)    ! Temporary space for Lagrnagian fields

end module PARTICLE_ARRAYS
