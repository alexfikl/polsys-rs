! SPDX-FileCopyrightText: 2024 Alexandru Fikl <alexfikl@gmail.com>
! SPDX-License-Identifier: MIT

module polsys_plp_wrapper
    use :: POLSYS, only:POLYNOMIAL, PARTITION, PARTITION_SIZES, BEZOUT_PLP, POLSYS_PLP
    use, intrinsic :: iso_c_binding, only: c_int32_t, c_double_complex, c_double, c_bool

    implicit none

    private

    public :: is_polynomial_allocated
    public :: deallocate_polynomial
    public :: init_polynomial

    public :: is_partition_allocated
    public :: deallocate_partition
    public :: init_partition

    public :: bezout_plp_wrap
    public :: polsys_plp_wrap

contains

    ! {{{ Polynomial

    subroutine is_polynomial_allocated(flag) bind(c)
        integer(c_int32_t), intent(out) :: flag

        if (allocated(POLYNOMIAL)) then
            flag = 1
        else
            flag = 0
        end if
    end subroutine is_polynomial_allocated

    !> @brief Deallocates the global polynomial value.
    !!
    !! This is automatically called by `init_polynomial`.
    subroutine deallocate_polynomial(ierr) bind(c)
        ! routine arguments
        integer(c_int32_t), intent(out) :: ierr

        ! local variables
        integer :: i, j

        ierr = 0
        if (allocated(POLYNOMIAL)) then
            do i = 1, size(POLYNOMIAL)
                if (ASSOCIATED(POLYNOMIAL(i)%TERM)) then
                    do j = 1, POLYNOMIAL(i)%NUM_TERMS
                        if (ASSOCIATED(POLYNOMIAL(i)%TERM(j)%DEG)) then
                            deallocate (POLYNOMIAL(i)%TERM(j)%DEG, stat=ierr)
                            if (ierr /= 0) then
                                ierr = ierr + 3
                                return
                            end if
                            nullify (POLYNOMIAL(i)%TERM(j)%DEG)
                        end if
                    end do

                    deallocate (POLYNOMIAL(i)%TERM, stat=ierr)
                    if (ierr /= 0) then
                        ierr = ierr + 3
                        return
                    end if
                    nullify (POLYNOMIAL(i)%TERM)
                end if
            end do

            deallocate (POLYNOMIAL, stat=ierr)
            if (ierr /= 0) then
                ierr = ierr + 3
            end if
        end if
    end subroutine deallocate_polynomial

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
    !!      3 x^2 + y^2 - 1 & = 0,
    !!      2 x y + y^2 - 3 & = 0.
    !!      \end{aligned}
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
        integer(c_int32_t), value, intent(in) :: n
        integer(c_int32_t), value, intent(in) :: m
        integer(c_int32_t), dimension(n), intent(in) :: n_coeffs_per_eq
        complex(c_double_complex), dimension(m), intent(in) :: coefficients
        integer(c_int32_t), dimension(n*m), intent(in) :: degrees
        integer(c_int32_t), intent(out) :: ierr

        ! local variables
        integer :: i, j, k, n_i

        if (m /= sum(n_coeffs_per_eq)) then
            ierr = 105
            return
        end if

        call deallocate_polynomial(ierr)
        if (ierr /= 0) then
            return
        end if

        allocate (POLYNOMIAL(n), stat=ierr)
        if (ierr /= 0) return

        k = 1
        do i = 1, n
            n_i = n_coeffs_per_eq(i)
            POLYNOMIAL(i)%NUM_TERMS = n_i
            allocate (POLYNOMIAL(i)%TERM(n_i), stat=ierr)
            if (ierr /= 0) return

            do j = 1, n_i
                allocate (POLYNOMIAL(i)%TERM(j)%DEG(n + 1), stat=ierr)
                if (ierr /= 0) return

                POLYNOMIAL(i)%TERM(j)%COEF = coefficients(k)
                POLYNOMIAL(i)%TERM(j)%DEG(1:n) = degrees(n * (k - 1) + 1:n * k)
                k = k + 1
            end do
        end do
    end subroutine init_polynomial

    ! }}}

    ! {{{ Partition

    subroutine is_partition_allocated(flag) bind(c)
        ! routine arguments
        integer(c_int32_t), intent(out) :: flag

        if (allocated(PARTITION) .and. allocated(PARTITION_SIZES)) then
            flag = 1
        else
            flag = 0
        end if

    end subroutine is_partition_allocated

    !> @brief Deallocate the global partition.
    !!
    !! This is automatically called by `init_partition`.
    subroutine deallocate_partition(ierr) bind(c)
        ! routine arguments
        integer(c_int32_t), intent(out) :: ierr

        ! local variables
        integer :: i, j

        ierr = 0
        if (allocated(PARTITION)) then
            do i = 1, size(PARTITION)
                if (associated(PARTITION(i)%SET)) then
                    do j = 1, PARTITION_SIZES(i)
                        if (associated(PARTITION(i)%SET(j)%INDEX)) then
                            deallocate (PARTITION(i)%SET(j)%INDEX, stat=ierr)
                            if (ierr /= 0) then
                                ierr = ierr + 3
                                return
                            end if
                            nullify (PARTITION(i)%SET(j)%INDEX)
                        end if

                        if (associated(PARTITION(i)%SET(j)%START_COEF)) then
                            deallocate (PARTITION(i)%SET(j)%START_COEF, stat=ierr)
                            if (ierr /= 0) then
                                ierr = ierr + 3
                                return
                            end if
                            nullify (PARTITION(i)%SET(j)%START_COEF)
                        end if
                    end do

                    deallocate (PARTITION(i)%SET, stat=ierr)
                    if (ierr /= 0) then
                        ierr = ierr + 3
                        return
                    end if
                    nullify (PARTITION(i)%SET)
                end if
            end do

            deallocate (PARTITION, stat=ierr)
            if (ierr /= 0) then
                ierr = ierr + 3
                return
            end if
        end if

        if (allocated(PARTITION_SIZES)) then
            deallocate (PARTITION_SIZES, stat=ierr)
            if (ierr /= 0) then
                ierr = ierr + 3
                return
            end if
        end if
    end subroutine deallocate_partition

    !> @brief Initialize the global partition.
    !!
    !! @param[in] n         number of equations and variables
    !! @param[in] m         number of sets per equation
    !! @param[in] p         number of indices per equation
    subroutine init_partition(n, m, p, n_sets_per_partition, n_indices_per_set, &
                              indices, ierr) bind(c)
        ! routine arguments
        integer(c_int32_t), value, intent(in) :: n
        integer(c_int32_t), value, intent(in) :: m
        integer(c_int32_t), value, intent(in) :: p
        integer(c_int32_t), dimension(n), intent(in) :: n_sets_per_partition
        integer(c_int32_t), dimension(m), intent(in) :: n_indices_per_set
        integer(c_int32_t), dimension(p), intent(in) :: indices
        integer(c_int32_t), intent(out) :: ierr

        ! local variables
        integer :: i, j, g_i, g_j, n_i, n_j

        call deallocate_partition(ierr)
        if (ierr /= 0) then
            return
        end if

        allocate (PARTITION_SIZES(n), stat=ierr)
        if (ierr /= 0) return

        allocate (PARTITION(n), stat=ierr)
        if (ierr /= 0) return

        ! global indices into n_indices_per_set and indices
        g_i = 1
        g_j = 1

        PARTITION_SIZES(1:n) = n_sets_per_partition(1:n)
        do i = 1, n             ! for each partition
            n_i = n_sets_per_partition(i)

            allocate (PARTITION(i)%SET(n_i), stat=ierr)
            if (ierr /= 0) return

            do j = 1, n_i       ! for each set in the partition
                n_j = n_indices_per_set(g_i + j - 1)

                allocate (PARTITION(i)%SET(j)%INDEX(n_j), stat=ierr)
                if (ierr /= 0) return

                PARTITION(i)%SET(j)%NUM_INDICES = n_j
                PARTITION(i)%SET(j)%INDEX(1:n_j) = indices(g_j:g_j + n_j - 1)

                ! NOTE: these are filled in by POLSYS_PLP during the solve
                PARTITION(i)%SET(j)%SET_DEG = 0
                nullify (PARTITION(i)%SET(j)%START_COEF)

                g_j = g_j + n_j
            end do

            g_i = g_i + n_i
        end do
    end subroutine init_partition

    ! }}}

    ! {{{ Bezout

    subroutine bezout_plp_wrap(n, maxt, tol, bplp, ierr) bind(c)
        ! routine arguments
        integer(c_int32_t), value, intent(in) :: n
        integer(c_int32_t), value, intent(in) :: maxt
        real(c_double), value, intent(in):: tol
        integer(c_int32_t), intent(out) :: bplp
        integer(c_int32_t), intent(out) :: ierr

        ! local variables
        real(c_double) :: tol_plp

        ierr = 0
        bplp = 0

        if (tol <= 0.0_c_double) then
            ierr = 106
            return
        end if

        if (.not. allocated(POLYNOMIAL)) then
            ierr = 101
            return
        end if

        if (.not. allocated(PARTITION_SIZES) .or. .not. allocated(PARTITION)) then
            ierr = 102
            return
        end if

        tol_plp = tol
        call BEZOUT_PLP(n, maxt, tol_plp, bplp)

    end subroutine bezout_plp_wrap

    ! }}}

    ! {{{ POLSYS_PLP

    subroutine polsys_plp_wrap(n, tracktol, finaltol, singtol, &
                               sspar, bplp, iflag1, iflag2, &
                               arclen, lambda, roots, nfe, scale_factors, &
                               numrr, recall, no_scaling, user_f_df) bind(c)
        ! routine arguments
        integer(c_int32_t), intent(in), value :: n
        real(c_double), intent(in), value:: tracktol
        real(c_double), intent(in), value:: finaltol
        real(c_double), intent(in), value:: singtol
        real(c_double), dimension(8), intent(inout) :: sspar
        integer(c_int32_t), intent(in) :: bplp
        integer(c_int32_t), intent(out) :: iflag1
        integer(c_int32_t), dimension(bplp), intent(out) :: iflag2
        real(c_double), dimension(bplp), intent(out) :: arclen
        real(c_double), dimension(bplp), intent(out) :: lambda
        complex(c_double_complex), dimension((n + 1)*bplp), intent(out) :: roots
        integer(c_int32_t), dimension(bplp), intent(out) :: nfe
        real(c_double), dimension(n), intent(out) :: scale_factors
        integer(c_int32_t), intent(in), value :: numrr
        integer(c_int32_t), intent(in), value :: recall
        integer(c_int32_t), intent(in), value :: no_scaling
        integer(c_int32_t), intent(in), value :: user_f_df

        ! local variables
        integer(c_int32_t) :: bplp_f
        real(c_double) :: singtol_f
        integer(c_int32_t), dimension(:), pointer :: iflag2_f
        real(c_double), dimension(:), pointer :: arclen_f
        real(c_double), dimension(:), pointer :: lambda_f
        integer(c_int32_t), dimension(:), pointer :: nfe_f
        complex(c_double_complex), dimension(:, :), pointer :: roots_f

        iflag1 = 0
        if (tracktol <= 0.0_c_double &
            .or. finaltol <= 0.0_c_double &
            .or. singtol <= 0.0_c_double) then
            iflag1 = 106
            return
        end if

        if (.not. allocated(POLYNOMIAL)) then
            iflag1 = 101
            return
        end if

        if (.not. allocated(PARTITION_SIZES) .or. .not. allocated(PARTITION)) then
            iflag1 = 102
            return
        end if

        singtol_f = singtol
        bplp_f = bplp

        call POLSYS_PLP(n, tracktol, finaltol, singtol_f, &
                        sspar, bplp_f, iflag1, iflag2_f, &
                        arclen_f, lambda_f, roots_f, nfe_f, scale_factors, &
                        numrr, &
                        recall == 1, &
                        no_scaling == 1, &
                        user_f_df == 1)

        if (iflag1 == 0) then
            iflag2 = iflag2_f
            arclen = arclen_f
            lambda = lambda_f
            roots = reshape(roots_f, [(n + 1) * bplp])
            nfe = nfe_f
        end if
    end subroutine polsys_plp_wrap

    ! }}}

end module polsys_plp_wrapper

subroutine TARGET_SYSTEM_USER(N, PROJ_COEF, XC, F, DF)
    use:: REAL_PRECISION, only: R8

    integer, intent(in):: N
    complex(kind=R8), intent(in), dimension(N + 1):: PROJ_COEF, XC
    complex(kind=R8), intent(out):: F(N), DF(N, N + 1)

    stop
end subroutine TARGET_SYSTEM_USER
