module FLT_Base
    implicit none
    private

    integer, parameter, public :: DNS_FILTER_BCS_PERIODIC = 0
    integer, parameter, public :: DNS_FILTER_BCS_BIASED = 1
    integer, parameter, public :: DNS_FILTER_BCS_FREE = 2
    integer, parameter, public :: DNS_FILTER_BCS_SOLID = 3
    integer, parameter, public :: DNS_FILTER_BCS_DIRICHLET = 4
    integer, parameter, public :: DNS_FILTER_BCS_NEUMANN = 5
    integer, parameter, public :: DNS_FILTER_BCS_ZERO = 6

end module FLT_Base
