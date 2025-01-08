module TLab_Types
    use TLab_Constants, only: wp, wi, MAX_VARS, MAX_PARS
    implicit none
    save

    type term_dt
        sequence
        integer type
        integer scalar(MAX_VARS)                ! fields defining this term
        logical active(MAX_VARS), lpadding(3)   ! fields affected by this term
        real(wp) parameters(MAX_PARS)
        real(wp) auxiliar(MAX_PARS)
        real(wp) vector(3)
    end type term_dt

end module TLab_Types
