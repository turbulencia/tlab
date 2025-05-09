#include "dns_error.h"
#include "dns_const.h"
#ifdef USE_MPI

#endif

#define C_FILE_LOC "PARTICLE_PDF"

!########################################################################
!########################################################################
program PARTICLE_BUILD_PDF

    use TLab_Constants, only: wp, wi
    use TLAB_VARS
    use TLab_Arrays
    use TLab_WorkFlow, only: TLab_Write_ASCII
    use IO_Fields
#ifdef USE_MPI
    use TLabMPI_PROCS, only: TLabMPI_Initialize
    use TLabMPI_Transpose, only: TLabMPI_Trp_Initialize
#endif
    use Thermodynamics
    use NavierStokes, only: NavierStokes_Initialize_Parameters
    use THERMO_AIRWATER
    use PARTICLE_VARS
    use PARTICLE_ARRAYS
    use PARTICLE_PROCS
    use TLab_Grid

    implicit none

! -------------------------------------------------------------------

    integer(wi) nitera_first, nitera_last, nitera_save
    integer(wi) i

    character*64 fname
    character*32 bakfile
    real(wp) params(0)
    
    bakfile = trim(adjustl(ifile))//'.bak'

    call TLab_Start

    call TLab_Initialize_Parameters(ifile)
#ifdef USE_MPI
    call TLabMPI_Initialize(ifile)
    call TLabMPI_Trp_Initialize(ifile)
#endif
    call NavierStokes_Initialize_Parameters(ifile)
    call Thermodynamics_Initialize_Parameters(ifile)
    call Particle_Initialize_Parameters(ifile)

!  CALL DNS_READ_LOCAL(ifile) !for nitera stuff

! Get the local information from the tlab.ini
    call ScanFile_Int(bakfile, ifile, 'Iteration', 'Start', '0', nitera_first)
    call ScanFile_Int(bakfile, ifile, 'Iteration', 'End', '0', nitera_last)
    call ScanFile_Int(bakfile, ifile, 'Iteration', 'Restart', '50', nitera_save)

! -------------------------------------------------------------------
! Allocating memory space
! -------------------------------------------------------------------
    allocate (wrk1d(isize_wrk1d, inb_wrk1d))
    allocate (wrk2d(isize_wrk2d, 1))
    allocate (wrk3d(isize_wrk3d))

    allocate (s(isize_field, inb_scal_array))

    if (part%type == PART_TYPE_BIL_CLOUD_3 .or. part%type == PART_TYPE_BIL_CLOUD_4) then !Allocte memory to read fields
        allocate (txc(isize_field, 3))
    end if

    inb_part_txc = 1

    call Particle_Initialize_Memory(C_FILE_LOC)

    isize_wrk2d = max(isize_wrk2d, jmax*inb_part_interp)

    ! IF (jmax_part .EQ. 1) THEN
    !    jmax_part   = jmax ! 1 by default
    ! ENDIF

! -------------------------------------------------------------------
! Read the grid
! -------------------------------------------------------------------
    call TLab_Grid_Read(gfile, x, y, z)
    call FDM_Initialize(ifile)

!#######################################################################
!Setup the file properties
!#######################################################################
    do i = nitera_first, nitera_last, nitera_save

!#######################################################################
!READ ALL FILES
!#######################################################################
        write (fname, *) i; fname = trim(adjustl(tag_scal))//trim(adjustl(fname))
        call IO_Read_Fields(fname, imax, jmax, kmax, itime, inb_scal, 0, s, params)
        call THERMO_AIRWATER_LINEAR(imax*jmax*kmax, s, s(1, inb_scal_array))

        write (fname, *) i; fname = trim(adjustl(tag_part))//trim(adjustl(fname))
        call IO_READ_PARTICLE(fname, l_g, l_q)

! ######################################################################
! Save particle pathlines for particle_pdf
! ######################################################################
        write (fname, *) i; fname = "particle_pdf."//trim(adjustl(fname))
        call PARTICLE_PDF(fname, s, l_g, l_q, l_txc)

    end do

    call TLab_Stop(0)
end program
