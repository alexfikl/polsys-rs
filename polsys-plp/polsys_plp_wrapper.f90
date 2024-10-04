! SPDX-FileCopyrightText: 2024 Alexandru Fikl <alexfikl@gmail.com>
! SPDX-License-Identifier: MIT

module polsys_plp_wrapper
    use :: POLSYS, only:POLYNOMIAL, PARTITION, PARTITION_SIZES
    use, intrinsic :: iso_c_binding, only: c_int, c_double_complex, c_double

    implicit none

    interface
        subroutine BEZOUT_PLP(N, MAXT, TOL, BPLP)
            use REAL_PRECISION, ONLY: R8

            INTEGER, INTENT(IN):: N, MAXT
            REAL(KIND=R8), INTENT(IN OUT):: TOL
            INTEGER, INTENT(OUT):: BPLP
        end subroutine
    end interface

contains

    ! {{{ Polynomial

    subroutine is_polynomial_allocated(flag) bind(c)
        integer(c_int), intent(out) :: flag

        if (allocated(POLYNOMIAL)) then
            flag = 1
        else
            flag = 0
        end if
    end subroutine

    !> @brief Deallocates the global polynomial value.
    !!
    !! This is automatically called by `init_polynomial`.
    subroutine deallocate_polynomial(ierr) bind(c)
        ! routine arguments
        integer(c_int), intent(out) :: ierr

        ! local variables
        integer :: i, j

        ierr = 0
        if (allocated(POLYNOMIAL)) then
            do i = 1, size(POLYNOMIAL)
                if (ASSOCIATED(POLYNOMIAL(i)%TERM)) then
                    do j = 1, POLYNOMIAL(i)%NUM_TERMS
                        if (ASSOCIATED(POLYNOMIAL(i)%TERM(j)%DEG)) then
                            deallocate (POLYNOMIAL(i)%TERM(j)%DEG, stat=ierr)
                            if (ierr .ne. 0) return
                            nullify (POLYNOMIAL(i)%TERM(j)%DEG)
                        end if
                    end do

                    deallocate (POLYNOMIAL(i)%TERM, stat=ierr)
                    if (ierr .ne. 0) return
                    nullify (POLYNOMIAL(i)%TERM)
                end if
            end do

            deallocate (POLYNOMIAL, stat=ierr)
        end if

    end subroutine

    !> @brief Initialize a polynomial from its coefficients.
    !!
    !! This function must be called before calling any solve method from the
    !! POLSYS library. It initializes an internal type that contains the description
    !! of the system of polynomial equations.
    !!
    !! Consider the following two-dimensional example
    !!
    !! \f[
    !!      \begin{aligned}
    !!      3 x^2 + y^2 - 1 & = 0, \\
    !!      2 x y + y^2 - 3 & = 0.
    !! \f]
    !!
    !! This function is not thread-safe as it initializes a global variable
    !! in the POPLSYS library.
    !!
    !! @param[in] n                 number of equations and number of variables.
    !! @param[in] m                 total number of coefficients
    !! @param[in] n_coeffs_per_eq   number of coefficients for each of the equations
    !!                              in the system. For the example above, we have
    !!                              a 3 coefficients in each equations, so this
    !!                              would be [3, 3].
    !! @param[in] coefficients      actual coefficients for each equations for
    !!                              a total of `m = sum(n_coeffs_per_eq)`. For
    !!                              example, for the system above, this would be
    !!                              [3, 1, -1, 2, 1, -3]
    !! @param[in] degrees           the degrees of each variable for each coefficient
    !!                              for a total of `n * m`.
    !!                              For example, for the system above, this would be
    !!                              two pairs [2, 0, 0, 2, 0, 0, 1, 1, 0, 2, 0, 0],
    !!                              i.e. for coefficient k and variable j, the
    !!                              degree is `degrees(n * (k - 1) + j)` (starting
    !!                              from k = j = 1).
    !! @param[out] ierr             Error code.
    subroutine init_polynomial(n, m, n_coeffs_per_eq, &
                               coefficients, degrees, ierr) bind(c)
        ! routine arguments
        integer(c_int), intent(in), value :: n
        integer(c_int), intent(in), value :: m
        integer(c_int), dimension(n), intent(in) :: n_coeffs_per_eq
        complex(c_double_complex), dimension(m), intent(in) :: coefficients
        integer(c_int), dimension(n*m), intent(in) :: degrees
        integer(c_int), intent(out) :: ierr

        ! local variables
        integer :: i, j, k, n_i

        if (m .ne. sum(n_coeffs_per_eq)) then
            ierr = 20
            return
        end if

        call deallocate_polynomial(ierr)
        if (ierr .ne. 0) then
            ! NOTE: bumping the error codes so we distinguish deallocation form
            ! allocation errors in the calling code
            ierr = ierr + 4
            return
        end if

        allocate (POLYNOMIAL(n), stat=ierr)
        if (ierr .ne. 0) return

        k = 1
        do i = 1, n
            n_i = n_coeffs_per_eq(i)
            POLYNOMIAL(i)%NUM_TERMS = n_i
            allocate (POLYNOMIAL(i)%TERM(n_i), stat=ierr)
            if (ierr .ne. 0) return

            do j = 1, n_i
                allocate (POLYNOMIAL(i)%TERM(j)%DEG(n + 1), stat=ierr)
                if (ierr .ne. 0) return

                POLYNOMIAL(i)%TERM(j)%COEF = coefficients(k)
                POLYNOMIAL(i)%TERM(j)%DEG(1:n) = degrees(n * (k - 1) + 1:n * k)
                k = k + 1
            end do
        end do

        ierr = 0
    end subroutine

    ! }}}

    ! {{{ Partition

    subroutine is_partition_allocated(flag) bind(c)
        ! routine arguments
        integer(c_int), intent(out) :: flag

        if (allocated(PARTITION) .and. allocated(PARTITION_SIZES)) then
            flag = 1
        else
            flag = 0
        end if

    end subroutine

    !> @brief Deallocate the global partition.
    !!
    !! This is automatically called by `init_partition`.
    subroutine deallocate_partition(ierr) bind(c)
        ! routine arguments
        integer(c_int), intent(out) :: ierr

        ! local variables
        integer :: i, j

        ierr = 0
        if (allocated(PARTITION)) then
            do i = 1, size(PARTITION)
                if (associated(PARTITION(i)%SET)) then
                    do j = 1, PARTITION_SIZES(i)
                        if (associated(PARTITION(i)%SET(j)%INDEX)) then
                            deallocate (PARTITION(i)%SET(j)%INDEX, stat=ierr)
                            if (ierr .ne. 0) return
                            nullify (PARTITION(i)%SET(j)%INDEX)
                        end if

                        if (associated(PARTITION(i)%SET(j)%START_COEF)) then
                            deallocate (PARTITION(i)%SET(j)%START_COEF, stat=ierr)
                            if (ierr .ne. 0) return
                            nullify (PARTITION(i)%SET(j)%START_COEF)
                        end if
                    end do

                    deallocate (PARTITION(i)%SET, stat=ierr)
                    if (ierr .ne. 0) return
                    nullify (PARTITION(i)%SET)
                end if
            end do

            deallocate (PARTITION, stat=ierr)
            if (ierr .ne. 0) return
        end if

        if (allocated(PARTITION_SIZES)) then
            deallocate (PARTITION_SIZES, stat=ierr)
            if (ierr .ne. 0) return
        end if
    end subroutine

    !> @brief Initialize the global partition.
    !!
    !! @param[in] n         number of equations and variables
    !! @param[in] m         number of sets per equation
    !! @param[in] p         number of indices per equation
    subroutine init_partition(n, m, p, n_sets_per_partition, n_indices_per_set, &
                              indices, ierr) bind(c)
        ! routine arguments
        integer(c_int), intent(in), value :: n
        integer(c_int), intent(in), value :: m
        integer(c_int), intent(in), value :: p
        integer(c_int), dimension(n), intent(in) :: n_sets_per_partition
        integer(c_int), dimension(m), intent(in) :: n_indices_per_set
        integer(c_int), dimension(p), intent(in) :: indices
        integer(c_int), intent(out) :: ierr

        ! local variables
        integer :: i, j, g_i, g_j, n_i, n_j

        call deallocate_partition(ierr)
        if (ierr .ne. 0) then
            ! NOTE: bumping the error codes so we distinguish deallocation form
            ! allocation errors in the calling code
            ierr = ierr + 4
            return
        end if

        allocate (PARTITION_SIZES(n), stat=ierr)
        if (ierr .ne. 0) return

        allocate (PARTITION(n), stat=ierr)
        if (ierr .ne. 0) return

        ! global indices into n_indices_per_set and indices
        g_i = 1
        g_j = 1

        PARTITION_SIZES(1:n) = n_sets_per_partition(1:n)
        do i = 1, n             ! for each partition
            n_i = n_sets_per_partition(i)

            allocate (PARTITION(i)%SET(n_i), stat=ierr)
            if (ierr .ne. 0) return

            do j = 1, n_i       ! for each set in the partition
                n_j = n_indices_per_set(g_i + j - 1)

                allocate (PARTITION(i)%SET(j)%INDEX(n_j), stat=ierr)
                if (ierr .ne. 0) return

                PARTITION(i)%SET(j)%NUM_INDICES = n_j
                PARTITION(i)%SET(j)%INDEX(1:n_j) = indices(g_j:g_j + n_j - 1)

                ! NOTE: these are filled in by POLSYS_PLP during the solve
                ! PARTITION(i)%SET(j)%SET_DEG = 0
                ! PARTITION(i)%SET(j)%START_COEF = 0

                g_j = g_j + n_j
            end do

            g_i = g_i + n_i
        end do

        ierr = 0
    end subroutine

    ! }}}

    ! {{{ Bezout

    subroutine bezout_plp_wrapper(n, maxt, tol, bplp, ierr) bind(c)
        ! routine arguments
        integer(c_int), intent(in), value :: n
        integer(c_int), intent(in), value :: maxt
        real(c_double), intent(in), value:: tol
        integer(c_int), intent(out) :: bplp
        integer(c_int), intent(out) :: ierr

        ! local variables
        real(c_double) :: tol_plp

        ierr = 0
        bplp = 0

        if (tol .le. 0.0_c_double) then
            ierr = 1
            return
        end if

        if (.not. allocated(POLYNOMIAL)) then
            ierr = 2
            return
        end if

        if (.not. allocated(PARTITION_SIZES) .or. .not. allocated(PARTITION)) then
            ierr = 3
            return
        end if

        tol_plp = tol
        call BEZOUT_PLP(n, maxt, tol_plp, bplp)

    end subroutine

    ! }}}
end module
