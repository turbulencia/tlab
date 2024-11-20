module Avg_Types
    use TLab_Constants, only: wi
    implicit none
    save

    type phaseavg_dt
        sequence
        logical :: active
        integer(wi) :: stride
        character(32) :: type
    end type phaseavg_dt

end module Avg_Types
