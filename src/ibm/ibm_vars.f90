!########################################################################
!# HISTORY / AUTHORS
!#
!# 2022/04/01 - J. Kostelecky
!#              Created
!#
!########################################################################
!# DESCRIPTION OF MODLE
!#   IBM based on:
!#     "A simple and scalable immersed boundary method for high-ﬁdelity 
!#      simulations of ﬁxed and moving objects on a Cartesian mesh"
!#      https://doi.org/10.1016/j.apm.2021.06.026
!#                    
!#
!########################################################################

module IBM_VARS

  use TLAB_CONSTANTS, only : MAX_NSP, wp, wi
  use IBM_TYPES

  implicit none

  save 

  ! ibm_scalar mode
  integer(wi)                                    :: imode_ibm_scal

  ! descriptive geometry fields (saved)
  real(wp),    dimension(:), allocatable, target :: eps, epsp                   ! eps indicator field
  integer(wi), dimension(:), allocatable         :: nobi,    nobj,   nobk       ! number of objects in i/j/k 
  integer(wi), dimension(:), allocatable         :: nobi_b,  nobj_b, nobk_b     ! beginn of objects in i/j/k 
  integer(wi), dimension(:), allocatable         :: nobi_e,  nobj_e, nobk_e     ! end    of objects in i/j/k
  integer(wi)                                    :: nobi_max, nobj_max, nobk_max, nob_max
  real(wp)                                       :: max_height_objlo, max_height_objup
  integer(wi), dimension(:), allocatable         :: ibm_case_x, ibm_case_y, ibm_case_z  ! Store IBM case (won't be valid for moving objects)
  logical                                        :: IBM_ini_case_x, IBM_ini_case_y, IBM_ini_case_z

  ! modified field
  real(wp),    dimension(:), allocatable, target :: fld_ibm                     ! with splines in solid regions

  ! boundary values of scalar fields 
  real(wp),    dimension(MAX_NSP)                :: ibmscaljmin, ibmscaljmax 

  ! work array for splines
  real(wp),    dimension(:), allocatable         :: xa, xb, ya, yb

  ! gammas for conditional averages & scalar boundary values applied in solids
  real(wp),    dimension(:),   allocatable       :: dy, facu, facl
  real(wp),    dimension(:),   allocatable       :: gamma_0, gamma_1
  real(wp),    dimension(:,:), allocatable       :: scal_bcs

  ! flag (decides which fdm calls are with modified fields)
  logical                      :: ibm_burgers, ibm_partial, ibm_objup

  ! read_local
  logical                      :: ibm_restart
  integer(wi)                  :: nflu                        ! number of fluid points used for Splines 
  integer(wi)                  :: ibm_io                      ! IBM IO Type 
  
  ! array sizes
  integer(wi)                  :: isize_nobi,    isize_nobj,    isize_nobk
  integer(wi)                  :: isize_nobi_be, isize_nobj_be, isize_nobk_be
  integer(wi)                  :: nspl
  integer(wi)                  :: isize_wrk1d_ibm

  ! check IBM procs (active/idle)
  logical                      :: ims_pro_ibm_x, ims_pro_ibm_y, ims_pro_ibm_z
  
  ! ibm_dt geometry type 
  type(ibm_geo_dt)             :: ibm_geo                   ! create new geometry here

  ! name of io eps
  character(len=32), parameter :: eps_name       = 'eps0.1'
  character(len=32), parameter :: epsp_name      = 'epsp0.1'
  character(len=32), parameter :: eps_name_real  = eps_name(1:4)
  character(len=32), parameter :: epsp_name_real = epsp_name(1:5)

  ! io types
  integer(wi),       parameter :: IBM_IO_REAL  = 1
  integer(wi),       parameter :: IBM_IO_INT   = 2
  integer(wi),       parameter :: IBM_IO_BIT   = 3

end module IBM_VARS

!########################################################################