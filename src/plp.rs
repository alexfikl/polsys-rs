// SPDX-FileCopyrightText: 2024 Alexandru Fikl <alexfikl@gmail.com>
// SPDX-License-Identifier: MIT

use num::complex::Complex64;
use std::collections::HashMap;
use std::fmt;
use std::iter;

use crate::bindings;

// {{{ Errors

#[derive(Eq, PartialEq, Debug)]
pub enum PolsysError {
    /// Error flag returned from the Fortran library when there is no error. This
    /// should never be returned by the Rust wrappers in this library.
    NoError = 0,
    /// Unknown error
    UnknownError = 1024,

    /// Error in system routine attempting to do allocation.
    AllocateSystemFailed = 1,
    /// An invalid data object has been specified for allocation.
    AllocateInvalidObject = 2,
    /// Both system and object errors in allocation.
    AllocateFailed = 3,
    /// Error in system routine attempting to do deallocation.
    DeallocateSystemFailed = 4,
    /// An invalid data object has been specified for deallocation.
    DeallocateInvalidObject = 5,
    /// Both system and object errors in deallocation.
    DeallocateFailed = 6,

    /// Dimensions or setup of input polynomial and partition do not match.
    PolynomialInvalid = -1,
    /// Terms with negative powers found in the polynomial.
    NegativePower = -2,
    /// One of the equations in the system is a constant.
    ConstantEquation = -3,
    /// Partition sizes do not add up to the number of equations.
    InconsistentPartitionSize = -4,
    /// A partition is defined more than once.
    RepeatedPartition = -5,
    /// Array sizes used for recall do not have valid sizes.
    InconsistentRecall = -6,
    /// Number of scale factors does not match number of equations.
    InsufficientScaleFactors = -7,

    /// The indices given in a partition does not match the number of variables.
    PartitionInvalidIndexCount = -8,
    /// The indices given in a partition have incorrect values (not in 1 <= i <= n).
    PartitionIncorrectIndex = -9,
}

impl From<i32> for PolsysError {
    fn from(flag: i32) -> Self {
        match flag {
            // initialization errors (returned by `init_polynomial`)
            0 => PolsysError::NoError,
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
            -8 => PolsysError::PartitionInvalidIndexCount,
            -9 => PolsysError::PartitionIncorrectIndex,
            _ => PolsysError::UnknownError,
        }
    }
}

impl PolsysError {
    fn as_str(&self) -> &'static str {
        match *self {
            PolsysError::NoError => "Finished successfully",
            PolsysError::UnknownError => "Unknown error ocured",
            // allocation
            PolsysError::AllocateSystemFailed => {
                "Error in system routine attempting to do allocation"
            }
            PolsysError::AllocateInvalidObject => {
                "An invalid data object has been specified for allocation"
            }
            PolsysError::AllocateFailed => {
                "Failed to allocate invalid object"
            }
            PolsysError::DeallocateSystemFailed => {
                "Error in system routine attempting to do deallocation"
            }
            PolsysError::DeallocateInvalidObject => {
                "An invalid data object has been specified for deallocation"
            }
            PolsysError::DeallocateFailed => {
                "Failed to allocate invalid object"
            }
            // polsys_plp
            PolsysError::PolynomialInvalid => {
                "Input polynomial and partitions dimensions or coefficient counts do not match"
            }
            PolsysError::NegativePower => {
                "Input polynomial has negative powers"
            }
            PolsysError::ConstantEquation => {
                "Constant equation (all degrees are zero)"
            }
            PolsysError::InconsistentPartitionSize => {
                "Partition sizes do not match up to system size"
            }
            PolsysError::RepeatedPartition => {
                "Repeated terms present in the partition"
            }
            PolsysError::InconsistentRecall => {
                "Recall is not consistent (BPLP differs)"
            }
            PolsysError::InsufficientScaleFactors => {
                "Less scale factors specified than system size"
            }
            // partition
            PolsysError::PartitionInvalidIndexCount => {
                "Partition index count does not match system size"
            }
            PolsysError::PartitionIncorrectIndex => {
                "Partition contains incorrect indices (not in 1 <= i <= n)"
            }
        }
    }
}

impl fmt::Display for PolsysError {
    fn fmt(&self, fmt: &mut fmt::Formatter<'_>) -> fmt::Result {
        fmt.write_str(self.as_str())
    }
}

#[derive(Eq, PartialEq, Debug)]
pub enum PathTrackingError {
    /// An unnknown flag returned by the routine.
    UnknownError = -1,

    /// Tracking tolerance was not met (should re-run with a higher tolerance).
    TrackingToleranceFailed = 2,
    /// Maximum number of steps allowed was exceeded.
    MaximumStepsExceeded = 3,
    /// Jacobian does not have full rank.
    BadJacobian = 4,
    /// The tracking algorithm has lost the zero curve of the homotopy map and
    /// is not making progress.
    ZeroCurveLost = 5,
    /// Normal flow Newton iteration failed to converge.
    NewtonConvergenceFailed = 6,
    /// Failed to find a root.
    RootSearchFailed = 7,
}

impl PathTrackingError {
    fn as_str(&self) -> &'static str {
        match *self {
            PathTrackingError::UnknownError => "An unknown error has occurred",
            PathTrackingError::TrackingToleranceFailed => "The specified error tolerance could not be met (increase TRACKTOL)",
            PathTrackingError::MaximumStepsExceeded => "The maximum number of steps allowed was exceeded (increase NUMRR)",
            PathTrackingError::BadJacobian => "The Jacobian matrix does not have full rank (the zero curve of the homotopy map cannot be followed any further)",
            PathTrackingError::ZeroCurveLost => "The zero curve of the homotopy has been lost and no progress can be made (TRACKTOL and FINALTOL may be too lenient)",
            PathTrackingError::NewtonConvergenceFailed => "The normal flow Newton iteration failed to converge (TRACKTOL or FINALTOL may be too stringent)",
            PathTrackingError::RootSearchFailed => "Failed to find a root in 10*NUMRR iterations"
        }
    }
}

impl fmt::Display for PathTrackingError {
    fn fmt(&self, fmt: &mut fmt::Formatter<'_>) -> fmt::Result {
        fmt.write_str(self.as_str())
    }
}

#[derive(Debug)]
#[allow(dead_code)]
pub struct PathTrackingResult(Result<u32, PathTrackingError>);

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

/// }}}

// {{{ Polynomial

#[derive(Clone, Debug)]
pub struct Polynomial<const N: usize> {
    /// A full description of the polynomial system as a vector of mappings from
    /// degree tuples (i_0, i_1, ..., i_n) to complex coefficients.
    system: Vec<HashMap<[i32; N], Complex64>>,

    /// Number of coefficients per equation. This can be used to slice into the
    /// coefficients and degrees arrays.
    n_coeffs_per_eq: Vec<i32>,
    /// A flat list of complex coefficients for the whole system, i.e. the first
    /// `n_coeffs_per_eq[0]` coefficients belong to the first equation.
    coefficients: Vec<Complex64>,
    /// A flat list of degrees for the whole system, i.e. the first
    /// `n * n_coeffs_per_eq[0]` degrees belong to the first equation.
    degrees: Vec<i32>,

    /// A flag to denote that the polynomial was correctly initialized.
    is_initialized: bool,
}

pub fn is_polynomial_allocated() -> bool {
    let mut flag: i32 = 0;

    unsafe { bindings::is_polynomial_allocated(&mut flag) }

    flag != 0
}

pub fn deallocate_polynomial() -> i32 {
    let mut ierr: i32 = 0;
    unsafe { bindings::deallocate_polynomial(&mut ierr) };

    ierr
}

impl<const N: usize> Polynomial<N> {
    pub fn new(system: Vec<HashMap<[i32; N], Complex64>>) -> Self {
        let n_coeffs_per_eq: Vec<i32> = system.iter().map(|c| c.len() as i32).collect();

        let n = system.len() as i32;
        let m = n_coeffs_per_eq.iter().sum::<i32>();
        let mut coefficients: Vec<Complex64> = Vec::with_capacity(m as usize);
        let mut degrees: Vec<i32> = Vec::with_capacity((n * m) as usize);

        for p in system.iter() {
            for (degree, c) in p.iter() {
                coefficients.push(*c);
                degrees.extend(degree.iter());
            }
        }

        Polynomial {
            system,
            n_coeffs_per_eq,
            coefficients,
            degrees,
            is_initialized: false,
        }
    }

    #[allow(dead_code)]
    fn deallocate(&mut self) -> Result<&mut Self, PolsysError> {
        let ierr = deallocate_polynomial();

        if ierr == 0 {
            self.is_initialized = false;
            Ok(self)
        } else {
            Err(PolsysError::from(ierr))
        }
    }

    pub fn init(&mut self) -> Result<&mut Self, PolsysError> {
        if self.is_initialized {
            return Ok(self);
        }

        let mut ierr: i32 = 0;

        unsafe {
            bindings::init_polynomial(
                self.system.len() as i32,
                self.coefficients.len() as i32,
                self.n_coeffs_per_eq.as_ptr(),
                self.coefficients.as_ptr(),
                self.degrees.as_ptr(),
                &mut ierr,
            )
        }

        if ierr == 0 {
            self.is_initialized = true;
            Ok(self)
        } else {
            Err(PolsysError::from(ierr))
        }
    }

    pub fn degrees(&self) -> Vec<i32> {
        let n = self.system.len();

        self.n_coeffs_per_eq
            .iter()
            .scan(0, |start, &m| {
                let from = *start as usize;
                let to = (*start + m) as usize;
                *start = to as i32;

                Some(
                    self.degrees[(n * from)..(n * to)]
                        .chunks(n)
                        .map(|chunk| chunk.iter().sum())
                        .max(),
                )
            })
            .map(|d| d.unwrap())
            .collect()
    }

    pub fn total_degree(&self) -> i32 {
        self.degrees().iter().product()
    }

    pub fn len(&self) -> usize {
        self.system.len()
    }

    pub fn is_empty(&self) -> bool {
        self.system.is_empty()
    }
}

// }}}

// {{{ Partition

pub struct Partition {
    pub n_sets_per_partition: Vec<i32>,
    pub n_indices_per_set: Vec<i32>,
    pub indices: Vec<i32>,

    pub is_initialized: bool,
}

pub fn is_partition_allocated() -> bool {
    let mut flag: i32 = 0;

    unsafe { bindings::is_partition_allocated(&mut flag) }

    flag != 0
}

pub fn deallocate_partition() -> i32 {
    let mut ierr: i32 = 0;
    unsafe { bindings::deallocate_partition(&mut ierr) };

    ierr
}

impl Partition {
    pub fn init(&mut self) -> Result<&mut Self, PolsysError> {
        if self.is_initialized {
            return Ok(self);
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
            self.is_initialized = true;
            Ok(self)
        } else {
            Err(PolsysError::from(ierr))
        }
    }

    pub fn len(&self) -> usize {
        self.n_sets_per_partition.len()
    }

    pub fn is_empty(&self) -> bool {
        self.n_sets_per_partition.is_empty()
    }
}

pub fn make_m_homogeneous_partition(
    n: usize,
    part: Vec<Vec<u32>>,
) -> Result<Partition, PolsysError> {
    let flat_part: Vec<u32> = part.concat();
    if flat_part.len() != n {
        return Err(PolsysError::PartitionInvalidIndexCount);
    }

    if flat_part.iter().any(|&i| i < 1 || i > (n as u32)) {
        return Err(PolsysError::PartitionIncorrectIndex);
    }

    let m = part.len();
    let n_sets_per_partition = iter::repeat(m as i32).take(n).collect();
    let n_indices_per_set_eq: Vec<i32> = part.iter().map(|x| x.len() as i32).collect();
    let n_indices_per_set = (0..n).flat_map(|_| n_indices_per_set_eq.clone()).collect();
    let indices = (0..n)
        .flat_map(|_| flat_part.iter())
        .map(|x| *x as i32)
        .collect();

    println!("n_sets_per_partition {:?}", n_sets_per_partition);
    println!("n_indices_per_set {:?}", n_indices_per_set);
    println!("indices {:?}", indices);

    Ok(Partition {
        n_sets_per_partition,
        n_indices_per_set,
        indices,
        is_initialized: false,
    })
}

pub fn make_homogeneous_partition(n: usize) -> Result<Partition, PolsysError> {
    make_m_homogeneous_partition(n, vec![(1..(n as u32) + 1).collect()])
}

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
    let indices = part.concat().concat().iter().map(|x| *x as i32).collect();

    Ok(Partition {
        n_sets_per_partition,
        n_indices_per_set,
        indices,
        is_initialized: false,
    })
}

// }}}

// {{{ Bezout

pub fn bezout<const N: usize>(
    poly: &mut Polynomial<N>,
    part: &mut Partition,
    tol: f64,
) -> Result<i32, PolsysError> {
    let _ = deallocate_polynomial();
    poly.init()?;
    let _ = deallocate_partition();
    part.init()?;

    let mut bplp: i32 = 0;
    let mut ierr: i32 = 0;

    unsafe {
        bindings::bezout_plp_wrapper(
            poly.len() as i32,
            *poly.n_coeffs_per_eq.iter().max().unwrap(),
            tol,
            &mut bplp,
            &mut ierr,
        )
    }

    if ierr == 0 {
        Ok(bplp)
    } else {
        Err(PolsysError::from(ierr))
    }
}

// }}}

// {{{ Tests

#[cfg(test)]
mod tests {
    use super::*;
    use num::complex::c64;
    use std::collections::HashMap;

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
        assert_eq!(err.as_str(), "Finished successfully");
        assert_eq!(format!("{err}"), "Finished successfully");
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

        assert!(is_polynomial_allocated());
    }

    #[test]
    fn test_polynomial_wrapper() {
        let mut poly = Polynomial::new(vec![
            HashMap::from([
                ([2, 0], c64(3.0, 1.3)),
                ([0, 2], c64(1.0, 0.1)),
                ([0, 0], c64(-1.0, 0.2)),
            ]),
            HashMap::from([
                ([1, 1], c64(2.0, 0.5)),
                ([0, 2], c64(1.0, 0.5)),
                ([0, 0], c64(-3.0, 1.0)),
            ]),
        ]);

        // NOTE: ensure the polynomial is deallocated
        let _ = deallocate_polynomial();
        assert!(!is_polynomial_allocated());

        poly.init().unwrap();
        assert!(is_polynomial_allocated());

        assert_eq!(poly.len(), 2);
        assert_eq!(poly.degrees(), [2, 2]);
        assert_eq!(poly.total_degree(), 4);
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
                2,
                4,
                n_sets_per_partition.as_ptr(),
                n_indices_per_set.as_ptr(),
                indices.as_ptr(),
                &mut ierr,
            );
        }

        assert_eq!(ierr, 0);
    }

    #[test]
    fn test_partition_m_homogeneous() {
        let mut part =
            make_m_homogeneous_partition(3, vec![vec![1, 2], vec![3]]).unwrap();

        let _ = deallocate_partition();
        assert!(!is_partition_allocated());

        part.init().unwrap();
        assert!(is_partition_allocated());

        assert_eq!(part.n_sets_per_partition, [2, 2, 2]);
        assert_eq!(part.n_indices_per_set, [2, 1, 2, 1, 2, 1]);
        assert_eq!(part.indices, [1, 2, 3, 1, 2, 3, 1, 2, 3]);
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

        let _ = deallocate_partition();
        assert!(!is_partition_allocated());

        part.init().unwrap();
        assert!(is_partition_allocated());

        assert_eq!(part.n_sets_per_partition, [2, 2, 3]);
        assert_eq!(part.n_indices_per_set, [2, 1, 1, 2, 1, 1, 1]);
        assert_eq!(part.indices, [1, 2, 3, 1, 2, 3, 1, 2, 3]);
    }

    #[test]
    fn test_bezout_plp() {
        let mut poly = Polynomial::new(vec![
            HashMap::from([
                ([2, 0], c64(3.0, 1.3)),
                ([0, 2], c64(1.0, 0.1)),
                ([0, 0], c64(-1.0, 0.2)),
            ]),
            HashMap::from([
                ([1, 1], c64(2.0, 0.5)),
                ([0, 2], c64(1.0, 0.5)),
                ([0, 0], c64(-3.0, 1.0)),
            ]),
        ]);
        let mut part = make_homogeneous_partition(2).unwrap();

        let bplp = bezout(&mut poly, &mut part, 1.0e-8);
    }
}

// }}}
