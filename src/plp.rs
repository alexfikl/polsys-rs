// SPDX-FileCopyrightText: 2024 Alexandru Fikl <alexfikl@gmail.com>
// SPDX-License-Identifier: MIT

use num::complex::Complex64;
use std::collections::HashMap;
use std::fmt;
use std::iter;

use crate::bindings;

/// {{{ Errors

#[derive(Debug)]
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
}

impl From<i32> for PolsysError {
    fn from(flag: i32) -> Self {
        match flag {
            // initialization errors (returned by `init_polynomial`)
            0 => PolsysError::NoError,
            1 => PolsysError::AllocateSystemFailed,
            2 => PolsysError::AllocateInvalidObject,
            3 => PolsysError::AllocateFailed,
            11 => PolsysError::DeallocateSystemFailed,
            12 => PolsysError::DeallocateInvalidObject,
            13 => PolsysError::DeallocateFailed,
            // solver errors (returned by `polsys_plp`)
            -1 => PolsysError::PolynomialInvalid,
            -2 => PolsysError::NegativePower,
            -3 => PolsysError::ConstantEquation,
            -4 => PolsysError::InconsistentPartitionSize,
            -5 => PolsysError::RepeatedPartition,
            -6 => PolsysError::InconsistentRecall,
            -7 => PolsysError::InsufficientScaleFactors,
            _ => PolsysError::UnknownError,
        }
    }
}

impl PolsysError {
    fn as_str(&self) -> &'static str {
        match *self {
            PolsysError::NoError => "Finished successfully",
            PolsysError::UnknownError => "Unknown error ocured",
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
        }
    }
}

impl fmt::Display for PolsysError {
    fn fmt(&self, fmt: &mut fmt::Formatter<'_>) -> fmt::Result {
        fmt.write_str(self.as_str())
    }
}

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
                if n > 10 {
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
}

// }}}

// {{{ Partition

pub struct Partition {
    /// Number of sets per partition.
    pub m: usize,
    /// Number of indices per set.
    pub p: usize,
    pub n_sets_per_partition: Vec<i32>,
    pub n_indices_per_set: Vec<i32>,
    pub indices: Vec<i32>,

    pub is_initialized: bool,
}

impl Partition {
    pub fn init(&mut self) -> Result<&mut Self, PolsysError> {
        if self.is_initialized {
            return Ok(self);
        }

        let mut ierr: i32 = 0;

        unsafe {
            bindings::init_partition(
                self.len() as i32,
                self.m as i32,
                self.p as i32,
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
}

pub fn make_m_homogeneous_partition(n: usize, part: Vec<Vec<u32>>) -> Partition {
    let m = part.len();
    let p = 10;
    let n_sets_per_partition = iter::repeat(m as i32).take(n).collect();
    let n_indices_per_set_eq: Vec<i32> = part.iter().map(|x| x.len() as i32).collect();
    let n_indices_per_set = (0..n).flat_map(|_| n_indices_per_set_eq.clone()).collect();
    let indices = (0..n)
        .flat_map(|_| part.iter().flatten())
        .map(|x| *x as i32)
        .collect();

    Partition {
        m,
        p,
        n_sets_per_partition,
        n_indices_per_set,
        indices,
        is_initialized: false,
    }
}

pub fn make_homogeneous_partition(n: usize) -> Partition {
    make_m_homogeneous_partition(n, vec![(1..(n as u32) + 1).collect()])
}

// pub fn make_plp_homogeneous_partition() {}

// }}}
