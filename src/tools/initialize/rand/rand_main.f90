#include "dns_const.h"
#include "dns_error.h"

#define C_FILE_LOC "INIRAND"

program INIRAND
    use TLab_Constants, only: wp, wi
    use TLab_Constants, only: ifile, gfile, lfile
    use TLab_Time, only: itime, rtime
    use TLab_Memory, only: inb_flow, inb_scal
    use TLab_Arrays
    use TLab_WorkFlow, only: TLab_Write_ASCII, TLab_Stop, TLab_Start, fourier_on
    use TLab_Memory, only: TLab_Initialize_Memory
#ifdef USE_MPI
    use TLabMPI_PROCS, only: TLabMPI_Initialize
    use TLabMPI_Transpose, only: TLabMPI_Trp_Initialize
#endif
    use IO_Fields
    use TLab_Grid
    use FDM, only: FDM_Initialize
    use Thermodynamics, only: Thermodynamics_Initialize_Parameters
    use NavierStokes, only: NavierStokes_Initialize_Parameters
    use TLab_Background, only: TLab_Initialize_Background
    use OPR_Fourier, only: OPR_Fourier_Initialize
    use RAND_LOCAL

    implicit none

    ! -------------------------------------------------------------------
    integer(wi) iq, is

    ! ###################################################################
    call TLab_Start()

    call TLab_Initialize_Parameters(ifile)
#ifdef USE_MPI
    call TLabMPI_Initialize(ifile)
    call TLabMPI_Trp_Initialize(ifile)
#endif

    call TLab_Grid_Read(gfile, x, y, z)
    call FDM_Initialize(ifile)

    call NavierStokes_Initialize_Parameters(ifile)
    call Thermodynamics_Initialize_Parameters(ifile)

    call TLab_Consistency_Check()

    call Inirand_Initialize_Parameters(ifile)

    ! ###################################################################
    call TLab_Initialize_Memory(C_FILE_LOC)

    call TLab_Initialize_Background(ifile)

    if (fourier_on) then
        call OPR_Fourier_Initialize()
    end if

    call OPR_CHECK()

    ! ###################################################################
    itime = 0; rtime = 0.0_wp

    ! ###################################################################
    call TLab_Write_ASCII(lfile, 'Initializing random fields.')

    ! -------------------------------------------------------------------
    do iq = 1, inb_flow
        call RAND_FIELD(ucov(iq), q(1, iq), txc(1, 1), txc(1, 2), txc(1, 3))
    end do
    if (pdf%type == TYPE_DF_GAUSSIAN) then ! Gaussian PDF
        call RAND_COVARIANCE(ucov, q(:, 1), q(:, 2), q(:, 3))
    end if

    io_header_q(1)%params(1) = rtime
    call IO_Write_Fields('flow.rand', imax, jmax, kmax, itime, inb_flow, q, io_header_q(1:1))

    ! -------------------------------------------------------------------
    do is = 1, inb_scal
        call RAND_FIELD(ucov(is), s(1, is), txc(1, 1), txc(1, 2), txc(1, 3))
    end do

    io_header_s(:)%params(1) = rtime
    call IO_Write_Fields('scal.rand', imax, jmax, kmax, itime, inb_scal, s, io_header_s(1:inb_scal))

    call TLab_Stop(0)
end program INIRAND
