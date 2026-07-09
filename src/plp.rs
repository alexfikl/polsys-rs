// SPDX-FileCopyrightText: 2024 Alexandru Fikl <alexfikl@gmail.com>
// SPDX-License-Identifier: MIT

use num::complex::Complex64;
use std::array;
use std::iter;

use crate::bindings;

// {{{ Errors

#[derive(Eq, PartialEq, Debug, thiserror::Error)]
pub enum PolsysError {
    /// Error flag returned from the Fortran library when there is no error. This
    /// should never be returned by the Rust wrappers in this library.
    #[error("Finished successfully")]
    NoError = 0,
    /// Unknown error
    #[error("Unknown error occurred")]
    UnknownError = 1024,

    /// Error in system routine attempting to do allocation.
    #[error("Error in system routine attempting to do allocation")]
    AllocateSystemFailed = 1,
    /// An invalid data object has been specified for allocation.
    #[error("An invalid data object has been specified for allocation")]
    AllocateInvalidObject = 2,
    /// Both system and object errors in allocation.
    #[error("Failed to allocate")]
    AllocateFailed = 3,
    /// Error in system routine attempting to do deallocation.
    #[error("Error in system routine attempting to do deallocation")]
    DeallocateSystemFailed = 4,
    /// An invalid data object has been specified for deallocation.
    #[error("An invalid data object has been specified for deallocation")]
    DeallocateInvalidObject = 5,
    /// Both system and object errors in deallocation.
    #[error("Failed to deallocate")]
    DeallocateFailed = 6,

    /// The polynomial or partition data structure is malformed (wrong
    /// array sizes, missing arrays, mismatched dimensions).
    #[error("Polynomial or partition is improperly allocated")]
    PolynomialInvalid = -1,
    /// Terms with negative powers found in the polynomial.
    #[error("Input polynomial has negative powers")]
    NegativePower = -2,
    /// One of the equations in the system is a constant.
    #[error("Constant equation (all degrees are zero)")]
    ConstantEquation = -3,
    /// Partition sizes do not add up to the number of variables.
    #[error("Partition indices do not cover all variables")]
    InconsistentPartitionSize = -4,
    /// The sets in a partition component do not cover variables 1..N
    /// exactly once (some variable index is missing or duplicated).
    #[error("Partition does not cover variables exactly once")]
    RepeatedPartition = -5,
    /// The Bezout number or the path-tracking arrays from a previous
    /// solve are inconsistent with the current call.
    #[error("Recall data is inconsistent with previous solve")]
    InconsistentRecall = -6,
    /// Number of scale factors does not match number of equations.
    #[error("Less scale factors specified than system size")]
    InsufficientScaleFactors = -7,

    /// Input dimensions do not match.
    #[error("Input dimensions do not match")]
    DimensionMismatch = 100,
    /// Polynomial is not allocated on calls.
    #[error("Polynomial has not been initialized (call Polynomial::init)")]
    PolynomialNotAllocated = 101,
    /// Partition was not allocated
    #[error("Partition has not been initialized (call Partition::init)")]
    PartitionNotAllocated = 102,
    /// The indices given in a partition does not match the number of variables.
    #[error("Partition index count does not match system size")]
    PartitionInvalidIndexCount = 103,
    /// The indices given in a partition have incorrect values (not in 1 <= i <= n).
    #[error("Partition contains incorrect indices (not in 1 <= i <= n)")]
    PartitionIncorrectIndex = 104,
    /// Number of coefficients does not match given coefficients per equation.
    #[error("Number of coefficients does not match number of coefficients per equation")]
    PolynomialInvalidCoefficientCount = 105,
    /// Invalid tolerance given (e.g. < 0).
    #[error("Invalid tolerance given (e.g. < 0)")]
    InvalidTolerance = 106,
    /// Invalid Bezout number (negative).
    #[error("Bezout number is negative or otherwise invalid")]
    InvalidBezoutNumber = 107,
    /// Polynomial is already initialized.
    #[error("Polynomial has already been initialized (cannot call init again)")]
    PolynomialDoubleInitialize = 108,
    /// Partition is already initialized.
    #[error("Partition has already been initialized (cannot call init again)")]
    PartitionDoubleInitialize = 109,
}

impl From<i32> for PolsysError {
    fn from(flag: i32) -> Self {
        match flag {
            0 => PolsysError::NoError,
            // initialization errors (returned by `init_polynomial`)
            1 => PolsysError::AllocateSystemFailed,
            2 => PolsysError::AllocateInvalidObject,
            3 => PolsysError::AllocateFailed,
            4 => PolsysError::DeallocateSystemFailed,
            5 => PolsysError::DeallocateInvalidObject,
            6 => PolsysError::DeallocateFailed,
            // solver errors (returned by `polsys_plp`)
            -1 => PolsysError::PolynomialInvalid,
            -2 => PolsysError::NegativePower,
            -3 => PolsysError::ConstantEquation,
            -4 => PolsysError::InconsistentPartitionSize,
            -5 => PolsysError::RepeatedPartition,
            -6 => PolsysError::InconsistentRecall,
            -7 => PolsysError::InsufficientScaleFactors,
            // custom
            100 => PolsysError::DimensionMismatch,
            101 => PolsysError::PolynomialNotAllocated,
            102 => PolsysError::PartitionNotAllocated,
            103 => PolsysError::PartitionInvalidIndexCount,
            104 => PolsysError::PartitionIncorrectIndex,
            105 => PolsysError::PolynomialInvalidCoefficientCount,
            106 => PolsysError::InvalidTolerance,
            107 => PolsysError::InvalidBezoutNumber,
            108 => PolsysError::PolynomialDoubleInitialize,
            109 => PolsysError::PartitionDoubleInitialize,
            _ => PolsysError::UnknownError,
        }
    }
}

#[derive(Eq, PartialEq, Debug, thiserror::Error)]
pub enum PathTrackingError {
    /// An unknown flag returned by the routine.
    #[error("An unknown error has occurred")]
    UnknownError = -1,

    /// Tracking tolerance was not met (should re-run with a higher tolerance).
    #[error("The specified error tolerance could not be met (increase TRACKTOL)")]
    TrackingToleranceFailed = 2,
    /// Maximum number of steps allowed was exceeded.
    #[error("The maximum number of steps allowed was exceeded (increase NUMRR)")]
    MaximumStepsExceeded = 3,
    /// Jacobian does not have full rank.
    #[error(
        "The Jacobian matrix does not have full rank (the zero curve of the homotopy map cannot be followed any further)"
    )]
    BadJacobian = 4,
    /// The tracking algorithm has lost the zero curve of the homotopy map and
    /// is not making progress.
    #[error(
        "The zero curve of the homotopy has been lost and no progress can be made (TRACKTOL and FINALTOL may be too lenient)"
    )]
    ZeroCurveLost = 5,
    /// Normal flow Newton iteration failed to converge.
    #[error(
        "The normal flow Newton iteration failed to converge (TRACKTOL or FINALTOL may be too stringent)"
    )]
    NewtonConvergenceFailed = 6,
    /// Failed to find a root.
    #[error("Failed to find a root in 10*NUMRR iterations")]
    RootSearchFailed = 7,
}

/// Outcome of tracking a single homotopy path.
///
/// An `Ok(n)` means the path converged after `n` end-game retries.
/// An `Err(e)` reports the specific failure mode (see
/// [`PathTrackingError`]).
#[derive(Debug)]
pub struct PathTrackingResult(Result<u32, PathTrackingError>);

impl PathTrackingResult {
    /// Returns `true` if the path completed successfully.
    pub fn is_ok(&self) -> bool {
        self.0.is_ok()
    }
}

impl From<i32> for PathTrackingResult {
    fn from(flag: i32) -> Self {
        let result = match flag {
            2 => Err(PathTrackingError::TrackingToleranceFailed),
            3 => Err(PathTrackingError::MaximumStepsExceeded),
            4 => Err(PathTrackingError::BadJacobian),
            5 => Err(PathTrackingError::ZeroCurveLost),
            6 => Err(PathTrackingError::NewtonConvergenceFailed),
            7 => Err(PathTrackingError::RootSearchFailed),
            n => {
                if n > 10 && (n - 1) % 10 == 0 {
                    Ok(((n - 1) / 10) as u32)
                } else {
                    Err(PathTrackingError::UnknownError)
                }
            }
        };

        PathTrackingResult(result)
    }
}

// }}}

// {{{ Polynomial

/// A system of `N` polynomial equations in `N` variables.
///
/// Each equation is a list of terms, where each term pairs a degree tuple
/// (exponents for the `N` variables) with a complex coefficient. Constant
/// terms use an all-zeros degree tuple.
///
/// Use [`term`] to build term tuples and [`Polynomial::new`] to construct
/// the system.
#[derive(Clone, Debug)]
pub struct Polynomial<const N: usize> {
    /// Number of coefficients per equation. This can be used to slice into the
    /// coefficients and degrees arrays.
    n_coeffs_per_eq: [i32; N],
    /// A flat list of complex coefficients for the whole system, i.e. the first
    /// `n_coeffs_per_eq[0]` coefficients belong to the first equation.
    coefficients: Vec<Complex64>,
    /// A flat list of degrees for the whole system, i.e. the first
    /// `n * n_coeffs_per_eq[0]` degrees belong to the first equation.
    degrees: Vec<i32>,
}

fn is_polynomial_allocated() -> bool {
    let mut flag: i32 = 0;

    unsafe { bindings::is_polynomial_allocated(&mut flag) }

    flag != 0
}

fn deallocate_polynomial() -> i32 {
    let mut ierr: i32 = 0;
    unsafe { bindings::deallocate_polynomial(&mut ierr) };

    ierr
}

/// Builds a polynomial term from a degree tuple and a coefficient.
pub fn term<const N: usize>(deg: [i32; N], coeff: impl Into<Complex64>) -> ([i32; N], Complex64) {
    (deg, coeff.into())
}

impl<const N: usize> Polynomial<N> {
    /// Creates a polynomial system from a list of equations.
    ///
    /// `system` must have exactly `N` entries (one per equation). Each entry
    /// is a `Vec` of `(degree_tuple, coefficient)` pairs as returned by, e.g.,
    /// [`term`].
    ///
    /// Returns [`PolsysError::DimensionMismatch`] if `system.len() != N`.
    pub fn new(system: Vec<Vec<([i32; N], Complex64)>>) -> Result<Self, PolsysError> {
        if system.len() != N {
            return Err(PolsysError::DimensionMismatch);
        }

        let n_coeffs_per_eq: [i32; N] = array::from_fn(|i| system[i].len() as i32);

        let m = n_coeffs_per_eq.iter().sum::<i32>();
        let mut coefficients: Vec<Complex64> = Vec::with_capacity(m as usize);
        let mut degrees: Vec<i32> = Vec::with_capacity(N * m as usize);

        for p in system.iter() {
            for (degree, c) in p.iter() {
                coefficients.push(*c);
                degrees.extend(degree.iter());
            }
        }

        Ok(Polynomial {
            n_coeffs_per_eq,
            coefficients,
            degrees,
        })
    }

    fn deallocate(&mut self) -> Result<&mut Self, PolsysError> {
        let ierr = deallocate_polynomial();

        if ierr == 0 {
            Ok(self)
        } else {
            Err(PolsysError::from(ierr))
        }
    }

    /// Initialises the underlying Fortran library with this polynomial.
    ///
    /// Must be called once (and only once) before solving. Most users do not
    /// need to call this directly -- [`bezout`] and [`PolsysSolver::solve`]
    /// call it automatically.
    ///
    /// Returns [`PolsysError::PolynomialDoubleInitialize`] if a polynomial
    /// is already initialised.
    pub fn init(&mut self) -> Result<&mut Self, PolsysError> {
        if is_polynomial_allocated() {
            return Err(PolsysError::PolynomialDoubleInitialize);
        }

        let mut ierr: i32 = 0;

        unsafe {
            bindings::init_polynomial(
                N as i32,
                self.coefficients.len() as i32,
                self.n_coeffs_per_eq.as_ptr(),
                self.coefficients.as_ptr(),
                self.degrees.as_ptr(),
                &mut ierr,
            )
        }

        if ierr == 0 {
            Ok(self)
        } else {
            Err(PolsysError::from(ierr))
        }
    }

    /// Returns the degree of each equation in the system.
    ///
    /// The degree of an equation is the maximum total degree among its terms
    /// (sum of the exponents in the degree tuple).
    pub fn degrees(&self) -> [i32; N] {
        let mut start = 0usize;
        array::from_fn(|i| {
            let m = self.n_coeffs_per_eq[i] as usize;
            let from = N * start;
            let to = N * (start + m);
            let degree = self.degrees[from..to]
                .chunks(N)
                .map(|chunk| chunk.iter().sum())
                .max()
                .unwrap_or(0);
            start += m;
            degree
        })
    }

    /// Total degree of the system -- the product of the per-equation
    /// degrees. For a 1-homogeneous partition this equals the Bezout
    /// number (the number of paths tracked).
    pub fn total_degree(&self) -> i32 {
        self.degrees().iter().product()
    }

    /// Number of equations and variables (`N`).
    pub fn len(&self) -> usize {
        N
    }

    /// Returns `true` if `N == 0`.
    pub fn is_empty(&self) -> bool {
        N == 0
    }
}

impl<const N: usize> Drop for Polynomial<N> {
    fn drop(&mut self) {
        let _ = self.deallocate();
    }
}

// }}}

// {{{ Partition

/// A variable partition for the PLP homotopy formulation.
///
/// The partition groups variables into sets per equation.  Use
/// [`make_homogeneous_partition`], [`make_m_homogeneous_partition`], or
/// [`make_plp_homogeneous_partition`] to construct one.
#[derive(Clone, Debug)]
pub struct Partition {
    pub n_sets_per_partition: Vec<i32>,
    pub n_indices_per_set: Vec<i32>,
    pub indices: Vec<i32>,
}

fn is_partition_allocated() -> bool {
    let mut flag: i32 = 0;

    unsafe { bindings::is_partition_allocated(&mut flag) }

    flag != 0
}

fn deallocate_partition() -> i32 {
    let mut ierr: i32 = 0;
    unsafe { bindings::deallocate_partition(&mut ierr) };

    ierr
}

impl Partition {
    /// Initialises the underlying Fortran library with this partition.
    ///
    /// Like [`Polynomial::init`], this is usually called automatically by
    /// [`bezout`] and [`PolsysSolver::solve_with_partition`].
    ///
    /// Returns [`PolsysError::PartitionDoubleInitialize`] if a partition
    /// is already initialised.
    pub fn init(&mut self) -> Result<&mut Self, PolsysError> {
        if is_partition_allocated() {
            return Err(PolsysError::PartitionDoubleInitialize);
        }

        let mut ierr: i32 = 0;

        unsafe {
            bindings::init_partition(
                self.n_sets_per_partition.len() as i32,
                self.n_indices_per_set.len() as i32,
                self.indices.len() as i32,
                self.n_sets_per_partition.as_ptr(),
                self.n_indices_per_set.as_ptr(),
                self.indices.as_ptr(),
                &mut ierr,
            )
        }

        if ierr == 0 {
            Ok(self)
        } else {
            Err(PolsysError::from(ierr))
        }
    }

    fn deallocate(&mut self) -> Result<&mut Self, PolsysError> {
        let ierr = deallocate_partition();

        if ierr == 0 {
            Ok(self)
        } else {
            Err(PolsysError::from(ierr))
        }
    }

    /// Number of partition components (should equal `N`).
    pub fn len(&self) -> usize {
        self.n_sets_per_partition.len()
    }

    /// Returns `true` if the partition has no components.
    pub fn is_empty(&self) -> bool {
        self.n_sets_per_partition.is_empty()
    }
}

/// Builds an m-homogeneous partition.
///
/// `n` is the number of variables and `part` is a list of sets (each set
/// is a list of variable indices, 1-based).  The same sets are replicated
/// for every equation in the system.
///
/// Returns [`PolsysError::PartitionInvalidIndexCount`] if the total
/// number of indices across all sets does not equal `n`, or
/// [`PolsysError::PartitionIncorrectIndex`] if any index is not in
/// `1..=n`.
pub fn make_m_homogeneous_partition(
    n: usize,
    part: Vec<Vec<u32>>,
) -> Result<Partition, PolsysError> {
    let flat_count: usize = part.iter().map(|x| x.len()).sum();
    if flat_count != n {
        return Err(PolsysError::PartitionInvalidIndexCount);
    }

    if part.iter().flatten().any(|&i| i < 1 || i > n as u32) {
        return Err(PolsysError::PartitionIncorrectIndex);
    }

    let m = part.len();
    let n_sets_per_partition = iter::repeat_n(m as i32, n).collect();
    let n_indices_per_set_eq: Vec<i32> = part.iter().map(|x| x.len() as i32).collect();
    let n_indices_per_set = (0..n)
        .flat_map(|_| n_indices_per_set_eq.iter().copied())
        .collect();
    let indices = (0..n)
        .flat_map(|_| part.iter().flatten().copied())
        .map(|x| x as i32)
        .collect();

    Ok(Partition {
        n_sets_per_partition,
        n_indices_per_set,
        indices,
    })
}

impl Drop for Partition {
    fn drop(&mut self) {
        let _ = self.deallocate();
    }
}

/// Builds a 1-homogeneous partition (all variables in a single set).
///
/// Equivalent to `make_m_homogeneous_partition(n, vec![1..n])`.
pub fn make_homogeneous_partition(n: usize) -> Result<Partition, PolsysError> {
    make_m_homogeneous_partition(n, vec![(1..(n as u32) + 1).collect()])
}

/// Builds a PLP (partitioned linear product) partition.
///
/// `n` is the number of variables and `part` has one entry per equation.
/// Each entry is itself a list of sets, and each set is a list of
/// variable indices (1-based).
///
/// Unlike [`make_m_homogeneous_partition`], this allows each equation to
/// have a different partition structure.
///
/// Returns [`PolsysError::PartitionInvalidIndexCount`] if
/// `part.len() != n`.
pub fn make_plp_homogeneous_partition(
    n: usize,
    part: Vec<Vec<Vec<u32>>>,
) -> Result<Partition, PolsysError> {
    if part.len() != n {
        return Err(PolsysError::PartitionInvalidIndexCount);
    }

    let n_sets_per_partition = part.iter().map(|p| p.len() as i32).collect();
    let n_indices_per_set = part
        .iter()
        .flat_map(|p| p.iter().map(|pp| pp.len() as i32))
        .collect();
    let indices = part.iter().flatten().flatten().map(|&x| x as i32).collect();

    Ok(Partition {
        n_sets_per_partition,
        n_indices_per_set,
        indices,
    })
}

// }}}

// {{{ Bezout

/// Computes the generalised Bezout number for the given partition.
///
/// The Bezout number equals the number of homotopy paths that will be
/// tracked during the solve.  For a 1-homogeneous partition it is the
/// product of the per-equation degrees (see [`Polynomial::total_degree`]).
///
/// Calls [`Polynomial::init`] and [`Partition::init`] as needed.
///
/// Returns [`PolsysError::InvalidBezoutNumber`] if the result is negative.
pub fn bezout<const N: usize>(
    poly: &mut Polynomial<N>,
    part: &mut Partition,
    tol: f64,
) -> Result<usize, PolsysError> {
    poly.init()?;
    part.init()?;

    let mut bplp: i32 = 0;
    let mut ierr: i32 = 0;

    unsafe {
        bindings::bezout_plp_wrap(
            poly.len() as i32,
            *poly.n_coeffs_per_eq.iter().max().unwrap(),
            tol,
            &mut bplp,
            &mut ierr,
        )
    }

    if bplp < 0 {
        return Err(PolsysError::InvalidBezoutNumber);
    }

    if ierr == 0 {
        Ok(bplp as usize)
    } else {
        Err(PolsysError::from(ierr))
    }
}

// }}}

// {{{ Solve

/// Outcome of a [`PolsysSolver::solve`].
///
/// It contains the roots found, the per-path number of function evaluations,
/// and the per-path tracking status.
#[derive(Debug)]
pub struct SolveResult {
    /// Number of affine variables in the solved system. Each root occupies
    /// `n_vars + 1` contiguous entries in [`SolveResult::roots`].
    pub n_vars: usize,
    /// All roots, laid out contiguously: root `i` occupies the entries
    /// `roots[i * (n_vars + 1) .. (i + 1) * (n_vars + 1)]`. The final entry of
    /// each root is the homogenizing coordinate. Use
    /// [`SolveResult::affine_roots`] for per-root access.
    pub roots: Vec<Complex64>,
    /// Number of function evaluations along each tracked path.
    pub nfe: Vec<i32>,
    /// Tracking status for each path.
    pub path_status: Vec<PathTrackingResult>,
}

impl SolveResult {
    /// Number of roots found.
    pub fn n_roots(&self) -> usize {
        self.roots.len() / (self.n_vars + 1)
    }

    /// Iterates over the affine part of each root: the first `n_vars`
    /// coordinates of every root, dropping the trailing homogenizing
    /// coordinate. Each slice borrows from [`SolveResult::roots`] directly, so
    /// no allocation is performed.
    pub fn affine_roots(&self) -> impl Iterator<Item = &[Complex64]> {
        self.roots
            .chunks(self.n_vars + 1)
            .map(move |r| &r[..self.n_vars])
    }
}

/// Configurable driver for solving a polynomial system by globally convergent
/// homotopy continuation.
///
/// Build with [`PolsysSolver::new`] (or [`PolsysSolver::default`]), tune via
/// the builder methods, then call [`PolsysSolver::solve`].
#[derive(Clone, Copy, Debug)]
pub struct PolsysSolver {
    /// Path-tracking tolerance for the homotopy-continuation phase.
    pub tracktol: f64,
    /// Tolerance for accepting a root at the end of a tracked path.
    pub finaltol: f64,
    /// Tolerance for detecting singular or ill-conditioned roots.
    pub singtol: f64,
    /// Number of restarts attempted when a path fails to converge (the Fortran
    /// routine's `NUMRR` argument).
    pub n_path_steps: i32,
    /// Reuse the path data produced by a previous solve.
    pub recall: bool,
    /// Disable automatic scaling of the polynomial equations.
    pub no_scaling: bool,
    /// Seed for the Fortran random-number generator. `None` (the default)
    /// leaves the RNG at its implementation-defined default; `Some(v)` seeds
    /// with `v` (broadcast across the seed array) before every solve.
    pub seed: Option<i32>,
}

impl Default for PolsysSolver {
    fn default() -> Self {
        Self {
            tracktol: 1.0e-8,
            finaltol: 1.0e-8,
            singtol: 1.0e-8,
            n_path_steps: 1,
            recall: false,
            no_scaling: false,
            seed: None,
        }
    }
}

impl PolsysSolver {
    /// Creates a new solver with default settings.
    pub fn new() -> Self {
        Self::default()
    }

    /// Sets the path-tracking tolerance.
    pub fn with_tracktol(mut self, tracktol: f64) -> Self {
        self.tracktol = tracktol;
        self
    }

    /// Sets the root-acceptance tolerance.
    pub fn with_finaltol(mut self, finaltol: f64) -> Self {
        self.finaltol = finaltol;
        self
    }

    /// Sets the singularity-detection tolerance.
    pub fn with_singtol(mut self, singtol: f64) -> Self {
        self.singtol = singtol;
        self
    }

    /// Sets the number of restarts attempted on path failure.
    pub fn with_n_path_steps(mut self, n_path_steps: i32) -> Self {
        self.n_path_steps = n_path_steps;
        self
    }

    /// Sets whether to reuse path data from a previous solve.
    pub fn with_recall(mut self, recall: bool) -> Self {
        self.recall = recall;
        self
    }

    /// Sets whether to disable automatic equation scaling.
    pub fn with_no_scaling(mut self, no_scaling: bool) -> Self {
        self.no_scaling = no_scaling;
        self
    }

    /// Sets the RNG seed.
    ///
    /// `None` leaves the Fortran default and `Some(v)` broadcasts `v` across
    /// the seed array. This is in general not needed, but if you see the tracking
    /// fail over and over, you may want to remove any non-determinism.
    pub fn with_seed(mut self, seed: Option<i32>) -> Self {
        self.seed = seed;
        self
    }

    /// Solves the system using the default 1-homogeneous partition.
    pub fn solve<const N: usize>(
        &self,
        poly: &mut Polynomial<N>,
    ) -> Result<SolveResult, PolsysError> {
        let mut part = make_homogeneous_partition(poly.len())?;
        self.solve_with_partition(poly, &mut part)
    }

    /// Solves the system using an explicit variable partition.
    pub fn solve_with_partition<const N: usize>(
        &self,
        poly: &mut Polynomial<N>,
        part: &mut Partition,
    ) -> Result<SolveResult, PolsysError> {
        let bplp = bezout(poly, part, self.tracktol)?;

        let n = poly.len();
        let mut iflag1: i32 = 0;
        let mut sspar = [0.0f64; 8];

        let mut iflag2 = vec![0i32; bplp];
        let mut arclen = vec![0.0f64; bplp];
        let mut lambda = vec![0.0f64; bplp];
        let mut roots = vec![Complex64::new(0.0, 0.0); (n + 1) * bplp];
        let mut nfe = vec![0i32; bplp];
        let mut scale_factors = [0.0f64; N];

        unsafe {
            bindings::polsys_plp_wrap(
                n as i32,
                self.tracktol,
                self.finaltol,
                self.singtol,
                sspar.as_mut_ptr(),
                bplp as i32,
                &mut iflag1,
                iflag2.as_mut_ptr(),
                arclen.as_mut_ptr(),
                lambda.as_mut_ptr(),
                roots.as_mut_ptr(),
                nfe.as_mut_ptr(),
                scale_factors.as_mut_ptr(),
                self.n_path_steps,
                self.recall as i32,
                self.no_scaling as i32,
                0,
                self.seed.unwrap_or(-1),
            )
        }

        if iflag1 == 0 {
            let path_status = iflag2
                .iter()
                .map(|iflag| PathTrackingResult::from(*iflag))
                .collect();

            Ok(SolveResult {
                n_vars: n,
                roots,
                nfe,
                path_status,
            })
        } else {
            Err(PolsysError::from(iflag1))
        }
    }
}

// }}}

// {{{ Tests

#[cfg(test)]
mod tests {
    use super::*;
    use num::complex::c64;

    fn norm(a: &[Complex64], b: &[Complex64]) -> f64 {
        a.iter()
            .zip(b.iter())
            .map(|(x, y)| (x - y).norm_sqr())
            .sum::<f64>()
            .sqrt()
    }

    fn assert_roots_contain(found: &[&[Complex64]], reference: &[Vec<Complex64>], tol: f64) {
        for ref_root in reference {
            let nearest = found
                .iter()
                .map(|r| norm(r, ref_root))
                .fold(f64::INFINITY, f64::min);
            assert!(
                nearest < tol,
                "reference root {ref_root:?}: no match within tol={tol} (nearest={nearest})",
            );
        }
    }

    #[test]
    fn test_polysys_error_from_i32() {
        let err = PolsysError::from(0);
        assert_eq!(err, PolsysError::NoError);

        let err = PolsysError::from(-7);
        assert_eq!(err, PolsysError::InsufficientScaleFactors);

        let err = PolsysError::from(999);
        assert_eq!(err, PolsysError::UnknownError);
    }

    #[test]
    fn test_polsys_error_display() {
        let err = PolsysError::NoError;
        assert_eq!(err.to_string(), "Finished successfully");
    }

    #[test]
    fn test_path_tracking_result_from_i32() {
        let result = PathTrackingResult::from(2);
        match result.0 {
            Ok(_) => unreachable!(),
            Err(e) => assert_eq!(e, PathTrackingError::TrackingToleranceFailed),
        }

        let result = PathTrackingResult::from(1 + 10 * 2);
        match result.0 {
            Ok(n) => assert_eq!(n, 2),
            Err(_) => unreachable!(),
        }

        let result = PathTrackingResult::from(2 + 10 * 2);
        match result.0 {
            Ok(_) => unreachable!(),
            Err(e) => assert_eq!(e, PathTrackingError::UnknownError),
        }
    }

    #[test]
    fn test_polynomial_init() {
        let mut ierr: i32 = 0;
        let n = 2;
        let n_coeffs_per_eq = [3, 3];
        let coefficients = [
            c64(3.0, 1.3),
            c64(1.0, 0.1),
            c64(-1.0, 0.2),
            c64(2.0, 0.5),
            c64(1.0, 0.5),
            c64(-3.0, 1.0),
        ];
        let degrees = [[2, 0], [0, 2], [0, 0], [1, 1], [0, 2], [0, 0]].as_flattened();

        // NOTE: ensure the polynomial is deallocated
        let _ = deallocate_polynomial();
        assert!(!is_polynomial_allocated());

        unsafe {
            bindings::init_polynomial(
                n,
                coefficients.len() as i32,
                n_coeffs_per_eq.as_ptr(),
                coefficients.as_ptr(),
                degrees.as_ptr(),
                &mut ierr,
            );
        }

        assert_eq!(ierr, 0);
        assert!(is_polynomial_allocated());
        let _ = deallocate_polynomial();
    }

    #[test]
    fn test_polynomial_wrapper() {
        let mut poly = Polynomial::new(vec![
            vec![
                ([2, 0], c64(3.0, 1.3)),
                ([0, 2], c64(1.0, 0.1)),
                ([0, 0], c64(-1.0, 0.2)),
            ],
            vec![
                ([1, 1], c64(2.0, 0.5)),
                ([0, 2], c64(1.0, 0.5)),
                ([0, 0], c64(-3.0, 1.0)),
            ],
        ])
        .unwrap();

        poly.init().unwrap();
        assert!(is_polynomial_allocated());

        assert_eq!(poly.len(), 2);
        assert_eq!(poly.degrees(), [2, 2]);
        assert_eq!(poly.total_degree(), 4);

        poly.deallocate().unwrap();
    }

    #[test]
    fn test_partition_init() {
        let mut ierr: i32 = 0;

        let n_sets_per_partition = [2, 2, 2, 2, 2, 2, 2, 2];
        let n_indices_per_set = [
            [4, 4],
            [4, 4],
            [4, 4],
            [4, 4],
            [4, 4],
            [4, 4],
            [4, 4],
            [4, 4],
        ]
        .as_flattened();
        let indices = [
            [[1, 2, 5, 6], [3, 4, 7, 8]],
            [[1, 2, 5, 6], [3, 4, 7, 8]],
            [[1, 2, 5, 6], [3, 4, 7, 8]],
            [[1, 2, 5, 6], [3, 4, 7, 8]],
            [[1, 2, 5, 6], [3, 4, 7, 8]],
            [[1, 2, 5, 6], [3, 4, 7, 8]],
            [[1, 2, 5, 6], [3, 4, 7, 8]],
            [[1, 2, 5, 6], [3, 4, 7, 8]],
        ]
        .as_flattened()
        .as_flattened();

        unsafe {
            bindings::init_partition(
                n_sets_per_partition.len() as i32,
                n_indices_per_set.len() as i32,
                indices.len() as i32,
                n_sets_per_partition.as_ptr(),
                n_indices_per_set.as_ptr(),
                indices.as_ptr(),
                &mut ierr,
            );
        }

        assert_eq!(ierr, 0);
        let _ = deallocate_partition();
    }

    #[test]
    fn test_partition_m_homogeneous() {
        let mut part = make_m_homogeneous_partition(3, vec![vec![1, 2], vec![3]]).unwrap();

        part.init().unwrap();
        assert!(is_partition_allocated());

        assert_eq!(part.n_sets_per_partition, [2, 2, 2]);
        assert_eq!(part.n_indices_per_set, [2, 1, 2, 1, 2, 1]);
        assert_eq!(part.indices, [1, 2, 3, 1, 2, 3, 1, 2, 3]);

        part.deallocate().unwrap();
    }

    #[test]
    fn test_partition_plp_homogeneous() {
        let mut part = make_plp_homogeneous_partition(
            3,
            vec![
                vec![vec![1, 2], vec![3]],
                vec![vec![1], vec![2, 3]],
                vec![vec![1], vec![2], vec![3]],
            ],
        )
        .unwrap();

        part.init().unwrap();
        assert!(is_partition_allocated());

        assert_eq!(part.n_sets_per_partition, [2, 2, 3]);
        assert_eq!(part.n_indices_per_set, [2, 1, 1, 2, 1, 1, 1]);
        assert_eq!(part.indices, [1, 2, 3, 1, 2, 3, 1, 2, 3]);

        part.deallocate().unwrap();
    }

    #[test]
    fn test_bezout_plp() {
        let mut poly = Polynomial::new(vec![
            vec![
                ([2, 0], c64(3.0, 1.3)),
                ([0, 2], c64(1.0, 0.1)),
                ([0, 0], c64(-1.0, 0.2)),
            ],
            vec![
                ([1, 1], c64(2.0, 0.5)),
                ([0, 2], c64(1.0, 0.5)),
                ([0, 0], c64(-3.0, 1.0)),
            ],
        ])
        .unwrap();
        let mut part = make_homogeneous_partition(poly.len()).unwrap();

        let bplp = bezout(&mut poly, &mut part, 1.0e-8).unwrap();
        println!("BPLP: {}", bplp);
        assert_eq!(bplp, 4);

        poly.deallocate().unwrap();
        part.deallocate().unwrap();

        let bplp = bezout(&mut poly, &mut part, 1.0e-8).unwrap();
        println!("BPLP: {}", bplp);
        assert_eq!(bplp, 4);

        poly.deallocate().unwrap();
        part.deallocate().unwrap();
    }

    #[test]
    fn test_polsys_plp_2d_system() {
        let mut poly = Polynomial::new(vec![
            vec![
                ([2, 0], c64(3.0, 1.3)),
                ([0, 2], c64(1.0, 0.1)),
                ([0, 0], c64(-1.0, 0.2)),
            ],
            vec![
                ([1, 1], c64(2.0, 0.5)),
                ([0, 2], c64(1.0, 0.5)),
                ([0, 0], c64(-3.0, 1.0)),
            ],
        ])
        .unwrap();
        let mut part = make_homogeneous_partition(poly.len()).unwrap();

        let result = PolsysSolver::new()
            .solve_with_partition(&mut poly, &mut part)
            .unwrap();

        for (i, root) in result.affine_roots().enumerate() {
            println!("Root {i}: {root:?}");
        }

        assert_eq!(result.n_roots(), 4);
        assert_eq!(result.path_status.len(), 4);
        assert_eq!(result.nfe.len(), 4);
    }

    #[test]
    fn test_solve_roots_dimensions() {
        let mut poly = Polynomial::new(vec![
            vec![
                ([2, 0], c64(3.0, 1.3)),
                ([0, 2], c64(1.0, 0.1)),
                ([0, 0], c64(-1.0, 0.2)),
            ],
            vec![
                ([1, 1], c64(2.0, 0.5)),
                ([0, 2], c64(1.0, 0.5)),
                ([0, 0], c64(-3.0, 1.0)),
            ],
        ])
        .unwrap();
        let mut part = make_homogeneous_partition(poly.len()).unwrap();

        let result = PolsysSolver::new()
            .solve_with_partition(&mut poly, &mut part)
            .unwrap();

        assert_eq!(result.n_roots(), 4);
        assert_eq!(result.roots.len(), result.n_roots() * (result.n_vars + 1));
        for root in result.roots.chunks(result.n_vars + 1) {
            assert_eq!(root.len(), result.n_vars + 1);
        }

        assert_eq!(result.n_vars, 2);
        let affine: Vec<_> = result.affine_roots().collect();
        assert_eq!(affine.len(), 4);
        for root in &affine {
            assert_eq!(root.len(), 2);
        }
    }

    /// Solves the 3-variable system:
    /// ```text
    ///     x^2 + y + z - 1 = 0
    ///     x + y^2 + z - 1 = 0
    ///     x + y + z^2 - 1 = 0
    /// ```
    /// which has five finite solutions in C:
    ///     (1, 0, 0),
    ///     (0, 1, 0),
    ///     (0, 0, 1),
    ///     (-1+sqrt(2), -1+sqrt(2), -1+sqrt(2))
    ///     (-1-sqrt(2), -1-sqrt(2), -1-sqrt(2)).
    ///
    /// The roots (1, 0, 0), (0, 1, 0), and (0, 0, 1) are singular (Jacobian
    /// determinant is zero), so homotopy paths converge poorly there.
    ///
    /// Reference: D. Cox, J. Little, and Donal O'Shea. Ideals, Varieties,
    /// and Algorithms: An Introduction to Computational Algebraic Geometry and
    /// Commutative Algebra. Springer Science & Business Media, 2013, p. 122.
    #[test]
    fn test_polsys_plp_cox_little_oshea_3vars() {
        let mut poly = Polynomial::<3>::new(vec![
            vec![
                term([2, 0, 0], 1.0),
                term([0, 1, 0], 1.0),
                term([0, 0, 1], 1.0),
                term([0, 0, 0], -1.0),
            ],
            vec![
                term([1, 0, 0], 1.0),
                term([0, 2, 0], 1.0),
                term([0, 0, 1], 1.0),
                term([0, 0, 0], -1.0),
            ],
            vec![
                term([1, 0, 0], 1.0),
                term([0, 1, 0], 1.0),
                term([0, 0, 2], 1.0),
                term([0, 0, 0], -1.0),
            ],
        ])
        .unwrap();

        // 1-homogeneous partition: Bezout number equals the total degree 2^3 = 8.
        assert_eq!(poly.total_degree(), 8);

        let mut part = make_homogeneous_partition(poly.len()).unwrap();

        let result = PolsysSolver::new()
            .with_tracktol(1.0e-8)
            .with_finaltol(1.0e-14)
            .with_singtol(0.0)
            .solve_with_partition(&mut poly, &mut part)
            .unwrap();

        // Eight paths are tracked.
        assert_eq!(result.n_roots(), 8);

        let sqrt2 = 2.0f64.sqrt();
        let found: Vec<&[Complex64]> = result.affine_roots().collect();

        let reference: [Vec<Complex64>; 5] = [
            vec![c64(1.0, 0.0), c64(0.0, 0.0), c64(0.0, 0.0)],
            vec![c64(0.0, 0.0), c64(1.0, 0.0), c64(0.0, 0.0)],
            vec![c64(0.0, 0.0), c64(0.0, 0.0), c64(1.0, 0.0)],
            vec![c64(-1.0 + sqrt2, 0.0); 3],
            vec![c64(-1.0 - sqrt2, 0.0); 3],
        ];

        assert_roots_contain(&found, &reference, 1.0e-10);
    }

    #[test]
    fn test_polsys_plp_univariate() {
        let mut poly =
            Polynomial::<1>::new(vec![vec![term([2], 1.0), term([1], -3.0), term([0], 2.0)]])
                .unwrap();

        assert_eq!(poly.total_degree(), 2);

        let mut part = make_homogeneous_partition(poly.len()).unwrap();

        let result = PolsysSolver::new()
            .solve_with_partition(&mut poly, &mut part)
            .unwrap();

        assert_eq!(result.n_roots(), 2);

        let found: Vec<&[Complex64]> = result.affine_roots().collect();
        let reference = [vec![c64(1.0, 0.0)], vec![c64(2.0, 0.0)]];

        assert_roots_contain(&found, &reference, 1.0e-6);
    }
}

// }}}
