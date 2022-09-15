module PARTICLE_TYPES
    use TLAB_TYPES, only: cp, ci, longi
    implicit none
    save

    type particle_dt
        sequence
        integer(ci) np                                      ! number of particles in local processor
        integer(longi), dimension(:), allocatable :: tags   ! tags of particles in local processor
        integer(ci), dimension(:), allocatable :: nodes     ! closest node below in Y direction
    end type particle_dt

end module PARTICLE_TYPES