module Avg_Vars
    use TLab_Constants, only: wi, wp
    use Avg_Types, only: phaseavg_dt
    implicit none
    save

    type(phaseavg_dt) :: PhAvg
    real(wp), dimension(:), allocatable, target :: avg_flow, avg_stress, avg_p, avg_scal
    integer(wi) :: nxy, nxz, nyz, nz_total
    integer(wi) :: avg_planes
    character(len=32), parameter :: avgu_name = 'avg_flow'
    character(len=32), parameter :: avgstr_name = 'avg_stress'
    character(len=32), parameter :: avgp_name = 'avg_p'
    character(len=32), parameter :: avgs_name = 'avg_scal'

    public :: avg_flow, avg_p, avg_scal, avg_stress, avg_planes
    public :: PhAvg

end module Avg_Vars