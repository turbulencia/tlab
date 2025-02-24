#include "dns_const.h"
#include "dns_error.h"

#define C_FILE_LOC "INIRAND"

program INIRAND
    use TLab_Constants, only: wp, wi
    use TLab_Constants, only: ifile, gfile, lfile
    use TLab_Memory, only: inb_txc
    use TLab_Time, only: itime, rtime
    use TLab_Memory, only: inb_flow, inb_scal
    use TLab_Arrays
    use TLab_WorkFlow, only: TLab_Write_ASCII, TLab_Stop, TLab_Start, fourier_on
    use TLab_Memory, only: TLab_Initialize_Memory
    use FDM, only: g, FDM_Initialize
#ifdef USE_MPI
    use TLabMPI_PROCS, only: TLabMPI_Initialize
    use TLabMPI_VARS, only: ims_pro
    use TLabMPI_Transpose, only: TLabMPI_Transpose_Initialize
#endif
    use Thermodynamics, only: Thermodynamics_Initialize_Parameters
    use NavierStokes, only: NavierStokes_Initialize_Parameters
    use TLab_Background, only: TLab_Initialize_Background
    use TLab_Grid
    use IO_Fields
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

    call TLab_Grid_Read(gfile, x, y, z, [g(1)%size, g(2)%size, g(3)%size])
    call FDM_Initialize(g(1), x)
    call FDM_Initialize(g(2), y)
    call FDM_Initialize(g(3), z)

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
    call IO_Write_Fields('flow.rand', imax, jmax, kmax, itime, inb_flow, q, io_header_q(1:1))

    do is = 1, inb_scal
        call RAND_FIELD(ucov(is), s(1, is), txc(1, 1), txc(1, 2), txc(1, 3))
    end do
    call IO_Write_Fields('scal.rand', imax, jmax, kmax, itime, inb_scal, s, io_header_s(1:inb_scal))

    call TLab_Stop(0)
end program INIRAND
