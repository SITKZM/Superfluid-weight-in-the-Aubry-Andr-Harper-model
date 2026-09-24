! ifx -o b.out Dxx_AAH.f90 -qmkl -qopenmp -mcmodel=medium -shared-intel -heap-arrays
!
! Superfluid weight (conventional + geometric contributions) of a 1D chain
! superconductor with a quasi-periodic (Aubry-Andre-Harper) on-site potential.
!
! The order parameter Delta is assumed to be uniform over the sites and is set
! here as a program parameter (delta_mag below), instead of being read from the
! output of the self-consistent calculation.
program main
    use,intrinsic :: iso_fortran_env
    implicit none

    ! --- Constants ---
    double precision, parameter :: pi = 4.0d0*atan(1.0d0)
    double precision, parameter :: tau = (1 + sqrt(5.0d0)) / 2
    double precision :: tau_a
    complex(kind(0d0)), parameter :: iunit = (0.0d0, 1.0d0)

    ! =========================================================
    ! Parameters set inside the program
    ! =========================================================
    ! Uniform order parameter |Delta|
    double precision, parameter :: delta_mag = 1.0d-2

    ! Number of k points in the (reduced) Brillouin zone
    integer, parameter :: Nk = 100

    ! --- Parameters read from the external file (params.in) ---
    double precision :: temperature ! fixed to absolute zero in this program
    integer :: N
    double precision :: U, t, mu, lambda
    double precision :: a, b
    integer :: n_c
    double precision :: criterion
    integer :: Nstop

    ! U, a, b, n_c, criterion, Nstop and temperature are not used here; they are kept so that the same params.in can be shared with the other programs.
    namelist /input_params/ temperature, N, U, t, mu, lambda, a, b, n_c, criterion, Nstop

    ! --- Internal variables ---
    integer :: M, ik
    double precision :: k_val, dk
    character(256) :: Dxx_file
    integer :: ios, un

    ! Decomposed components of the superfluid weight
    double precision :: Dxx_tot, Dxx_conv, Dxx_geom

    ! Dynamic arrays (all in double precision)
    complex(kind(0d0)), allocatable :: H0(:, :)
    complex(kind(0d0)), allocatable :: Jx(:, :)

    ! Workspace for ZHEEVD
    double precision, allocatable :: epsilons(:)
    integer :: lwork, lrwork, liwork, info
    complex(kind(0d0)), allocatable :: work(:)
    double precision, allocatable :: rwork(:)
    integer, allocatable :: iwork(:)

    integer(int64) :: time_begin_c, time_end_c, CountPerSec, CountMax

    ! =========================================================
    ! 1. Read and set the parameters
    ! =========================================================
    open(newunit=un, file="params.in", status="old", iostat=ios)
    if (ios /= 0) stop "Error: params.in not found."
    read(un, nml=input_params)
    close(un)

    ! =========================================================
    ! 2. Dynamic memory allocation
    ! =========================================================
    allocate(H0(N, N), Jx(N, N))

    ! The output file is labelled by Delta (a parameter here) instead of U.
    write(Dxx_file, '("./Dxx_AAH/Delta=", F6.4, "_mu=", F4.2, "_lambda=", F4.2, ".txt")') &
        delta_mag, mu, lambda

    ! =========================================================
    ! 3. Calculation of the superfluid weight tensor and file output
    ! =========================================================
    call system_clock(time_begin_c, CountPerSec, CountMax)

    ! Explicit evaluation of the workspace sizes
    lwork = N**2 + 2*N
    lrwork = 2*N**2 + 5*N + 1
    liwork = 5*N + 3
    allocate(epsilons(N))
    allocate(work(lwork), rwork(lrwork), iwork(liwork))

    ! Number of periods M and period tau_a of the rational approximant
    M = nint(dble(N) / tau)
    tau_a = dble(N) / dble(M)

    ! Discretization of the wave number and start of the loop
    dk = (2.0d0 * pi / dble(N)) / dble(Nk)
    Dxx_tot = 0.0d0; Dxx_conv = 0.0d0; Dxx_geom = 0.0d0
    do ik = 0, Nk - 1
        k_val = -pi / dble(N) + dk*(dble(ik) + 0.5d0)
        ! Build the current matrix in the band basis from the normal-state Hamiltonian
        call make_normal_hamiltonian(k_val)
        call ZHEEVD('V', 'U', N, H0, N, epsilons, work, lwork, rwork, lrwork, iwork, liwork, info)
        call make_current_operator(k_val)
        Jx = matmul(conjg(transpose(H0)), matmul(Jx, H0)) ! transform to the band basis

        ! Calculate the superfluid weight
        call calculate_superfluid_weight()
        print *, ik + 1, "completed."
    end do

    call write_files()

    call system_clock(time_end_c)
    print *, "Run time=", dble(time_end_c - time_begin_c) / CountPerSec, "sec"

contains
    subroutine make_normal_hamiltonian(k)
        double precision, intent(in) :: k
        integer :: i
        complex(kind(0d0)) :: phase

        H0 = (0.0d0, 0.0d0)
        do i = 1, N
            H0(i, i) = -mu + lambda*cos(2.0d0*pi*dble(i - 1) / tau_a)
        end do

        phase = exp(iunit*k)
        do i = 1, N - 1
            H0(i, i + 1) = -t*phase
            H0(i + 1, i) = -t*conjg(phase)
        end do
        H0(N, 1) = -t*phase
        H0(1, N) = -t*conjg(phase)
    end subroutine make_normal_hamiltonian

    subroutine make_current_operator(k)
        double precision, intent(in) :: k
        integer :: i
        complex(kind(0d0)) :: phase

        Jx = (0.0d0, 0.0d0)
        phase = exp(iunit*k)
        do i = 1, N - 1
            Jx(i, i + 1) = -iunit*t*phase
            Jx(i + 1, i) = conjg(-iunit*t*phase)
        end do
        Jx(N, 1) = -iunit*t*phase
        Jx(1, N) = conjg(-iunit*t*phase)
    end subroutine make_current_operator

    subroutine calculate_superfluid_weight()
        double precision :: Dxx_tot_k, Dxx_conv_k, Dxx_geom_k
        integer :: i, j
        double precision :: En, Em, weight, current_term
        double precision, parameter :: tol = 1.0d-10

        Dxx_tot_k = 0.0d0
        Dxx_conv_k = 0.0d0
        Dxx_geom_k = 0.0d0

        ! Double loop over the band indices i, j
        do i = 1, N
            En = sqrt(epsilons(i)**2 + delta_mag**2)
            do j = 1, N
                ! Evaluate E_j
                Em = sqrt(epsilons(j)**2 + delta_mag**2)
                current_term = abs(Jx(i, j))**2

                ! Weight of each term of the superfluid weight
                weight = (delta_mag**2 / dble(N) / dble(Nk))*(1.0d0 / (En*Em*(En + Em)))*2.0d0*current_term

                ! Add to the total (Dxx_tot_k)
                Dxx_tot_k = Dxx_tot_k + weight

                ! Dxx_conv_k collects the terms with \epsilon_n = \epsilon_m
                ! (judged within the tolerance tol to account for numerical error)
                if (abs(epsilons(i) - epsilons(j)) < tol) then
                    Dxx_conv_k = Dxx_conv_k + weight
                end if
            end do
        end do

        ! Dxx_geom_k is the total contribution minus the conventional one
        Dxx_geom_k = Dxx_tot_k - Dxx_conv_k

        ! Accumulate the result at wave number k over all wave numbers
        Dxx_tot = Dxx_tot + Dxx_tot_k
        Dxx_conv = Dxx_conv + Dxx_conv_k
        Dxx_geom = Dxx_geom + Dxx_geom_k
    end subroutine calculate_superfluid_weight

    subroutine write_files()
        integer :: uwrite
        open(newunit=uwrite, file=trim(Dxx_file))
        write(uwrite, '(A)') "# Dxx_tot        Dxx_conv        Dxx_geom"
        write(uwrite, '(3(ES15.6, 1x))') Dxx_tot, Dxx_conv, Dxx_geom
        close(uwrite)
        print *, "Dxx components successfully saved to: ", trim(Dxx_file)
    end subroutine write_files
end program main
