module PARTICLE_TYPES
    use TLAB_CONSTANTS, only: wp, wi, longi
    implicit none
    save

    type particle_dt
        sequence
        integer(wi) np                                      ! number of particles in local processor
        integer(longi), dimension(:), allocatable :: tags   ! tags of particles in local processor
        integer(wi), dimension(:), allocatable :: nodes     ! closest node below in Y direction
    end type particle_dt

end module PARTICLE_TYPES