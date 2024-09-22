! SPDX-FileCopyrightText: 2024 Alexandru Fikl <alexfikl@gmail.com>
! SPDX-License-Identifier: MIT

module polsys_plp_wrapper
    use POLSYS
    use iso_c_binding

    implicit none
contains

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
    !! @param[in] n_coeffs_per_eq   number of coefficients for each of the equations
    !!                              in the system. For the example above, we have
    !!                              a 3 coefficients in each equations, so this
    !!                              would be [3, 3].
    !! @param[in] coefficients      actual coefficients for each equations for
    !!                              a total of `m = sum(n_coeffs_per_eq)`. For
    !!                              example, for the system above, this would be
    !!                              [3, 1, -1, 2, 1, -3]
    !! @param[in] degrees           the degrees of each variable for each coefficient.
    !!                              For example, for the system above, this would be
    !!                              two pairs
    !!                              [[2, 0], [0, 2], [0, 0], [1, 1], [0, 2], [0, 0]].
    !! @param[out] ierr             Error code.
    subroutine init_polynomial(n, n_coeffs_per_eq, coefficients, degrees, ierr) bind(c)
        ! routine arguments
        integer(c_int), intent(in) :: n
        integer(c_int), dimension(n), intent(in) :: n_coeffs_per_eq
        complex(c_double_complex), dimension(:), intent(in) :: coefficients
        integer(c_int), dimension(:, :), intent(in) :: degrees
        integer(c_int), intent(out) :: ierr

        ! local variables
        integer :: i, j, k, n_i

        write(0, *) "In 'init_polynomial'"

        call CLEANUP_POL()
        allocate(POLYNOMIAL(n))

        k = 0
        do i = 1, n
            n_i = n_coeffs_per_eq(i)
            POLYNOMIAL(i)%NUM_TERMS = n_i
            allocate(POLYNOMIAL(i)%TERM(n_i))

            do j = 1, n_i
                allocate(POLYNOMIAL(i)%TERM(j)%DEG(n + 1))
                POLYNOMIAL(i)%TERM(j)%COEF = coefficients(k)
                POLYNOMIAL(i)%TERM(j)%DEG(1:n) = degrees(k, 1:n)
                k = k + 1
            end do
        end do

        ierr = 0
    end subroutine

    subroutine init_partition(n, num_sets, num_indices, set_index, ierr) bind(c)
        ! routine arguments
        integer(c_int), intent(in) :: n
        integer(c_int), dimension(:), intent(in) :: num_sets
        integer(c_int), dimension(:, :), intent(in) :: num_indices
        integer(c_int), dimension(:, :, :), intent(in) :: set_index
        integer(c_int), intent(out) :: ierr

        ! local variables
        integer :: i, j, n_i

        write(0, *) "In 'init_partition'"

        call CLEANUP_PAR()

        allocate(PARTITION_SIZES(n))
        PARTITION_SIZES(1:n) = num_sets(1:n)

        allocate(PARTITION(n))
        do i = 1, n
            allocate(PARTITION(i)%SET(NUM_SETS(i)))

            do j = 1, NUM_SETS(i)
                n_i = num_indices(i, j)

                allocate(PARTITION(i)%SET(j)%INDEX(n_i))
                PARTITION(i)%SET(j)%NUM_INDICES = n_i
                PARTITION(i)%SET(j)%INDEX(1:n_i) = set_index(i, j, 1:n_i)

                ! NOTE: these are filled in by POLSYS_PLP during the solve
                ! PARTITION(i)%SET(j)%SET_DEG = 0
                ! PARTITION(i)%SET(j)%START_COEF = 0
            end do
        end do

        ierr = 0
    end subroutine
end module
