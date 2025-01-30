#include "dns_const.h"
#include "dns_error.h"

#define C_FILE_LOC "INIRAND"

program INIRAND
    use TLab_Constants, only: wp, wi
    use TLab_Constants, only: ifile, gfile, lfile
    use TLAB_VARS
    use TLab_Arrays
    use TLab_WorkFlow, only: TLab_Write_ASCII, TLab_Stop, TLab_Start
    use TLab_Memory, only: TLab_Initialize_Memory
    use FDM, only: g, FDM_Initialize
#ifdef USE_MPI
    use TLabMPI_PROCS, only: TLabMPI_Initialize, ims_pro
    use TLabMPI_Transpose, only: TLabMPI_Transpose_Initialize
#endif
    use Thermodynamics, only: Thermodynamics_Initialize_Parameters
    use NavierStokes, only: NavierStokes_Initialize_Parameters
    use TLab_Background, only: TLab_Initialize_Background
    use IO_FIELDS
    use OPR_FOURIER
    use RAND_LOCAL

    implicit none

    ! -------------------------------------------------------------------
    integer(wi) iq, is

    ! ###################################################################
    call TLab_Start()

    call TLab_Initialize_Parameters(ifile)
#ifdef USE_MPI
    call TLabMPI_Initialize(ifile)
    call TLabMPI_Transpose_Initialize(ifile)
#endif

    call NavierStokes_Initialize_Parameters(ifile)
    call Thermodynamics_Initialize_Parameters(ifile)

    call TLab_Consistency_Check()

    call RAND_READ_LOCAL(ifile)

    inb_txc = 3

    call TLab_Initialize_Memory(C_FILE_LOC)

    call IO_READ_GRID(gfile, g(1)%size, g(2)%size, g(3)%size, g(1)%scale, g(2)%scale, g(3)%scale, wrk1d(:, 1), wrk1d(:, 2), wrk1d(:, 3))
    call FDM_Initialize(x, g(1), wrk1d(:, 1), wrk1d(:, 4))
    call FDM_Initialize(y, g(2), wrk1d(:, 2), wrk1d(:, 4))
    call FDM_Initialize(z, g(3), wrk1d(:, 3), wrk1d(:, 4))

    call TLab_Initialize_Background(ifile)

    ! ###################################################################
    call TLab_Write_ASCII(lfile, 'Initializing random fiels.')

#ifdef USE_MPI
    seed = seed + ims_pro         ! seed for random generator
#endif
    seed = -abs(seed)

    if (fourier_on) then
        call OPR_FOURIER_INITIALIZE()
    end if

    call OPR_CHECK()

    itime = 0; rtime = 0.0_wp

    do iq = 1, inb_flow
        call RAND_FIELD(ucov(iq), q(1, iq), txc(1, 1), txc(1, 2), txc(1, 3))
    end do
    if (ipdf == 2) then ! Gaussian PDF
        call RAND_COVARIANCE(ucov, q(:, 1), q(:, 2), q(:, 3))
    end if
    call IO_WRITE_FIELDS('flow.rand', imax, jmax, kmax, itime, inb_flow, q, io_header_q(1:1))

    do is = 1, inb_scal
        call RAND_FIELD(ucov(is), s(1, is), txc(1, 1), txc(1, 2), txc(1, 3))
    end do
    call IO_WRITE_FIELDS('scal.rand', imax, jmax, kmax, itime, inb_scal, s, io_header_s(1:inb_scal))

    call TLab_Stop(0)
end program INIRAND
