// SPDX-FileCopyrightText: 2024 Alexandru Fikl <alexfikl@gmail.com>
// SPDX-License-Identifier: MIT

use num::complex::Complex64;
use std::collections::HashMap;

use crate::bindings;

// {{{ InitializeError

#[derive(Debug)]
pub enum InitializeError {
    /// Error in system routine attempting to do allocation.
    SystemAllocateFailed = 1,
    /// An invalid data object has been specified for allocation.
    InvalidAllocateObject = 2,
    /// Both system and object errors in allocation.
    AllocateFailed = 3,
    /// Error in system routine attempting to do deallocation.
    SystemDeallocateFailed = 11,
    /// An invalid data object has been specified for deallocation.
    InvalidDeallocateObject = 12,
    /// Both system and object errors in deallocation.
    DeallocateFailed = 13,
}

impl From<i32> for InitializeError {
    fn from(flag: i32) -> Self {
        match flag {
            1 => InitializeError::SystemAllocateFailed,
            2 => InitializeError::InvalidAllocateObject,
            3 => InitializeError::AllocateFailed,
            11 => InitializeError::SystemDeallocateFailed,
            12 => InitializeError::InvalidDeallocateObject,
            13 => InitializeError::DeallocateFailed,
            _ => panic!("Unknown InitializeError value: {}", flag),
        }
    }
}

// }}}

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

    pub fn init(&mut self) -> Result<&mut Self, InitializeError> {
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
            Err(InitializeError::from(ierr))
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

pub struct Partition {}

pub fn make_homogeneous_partition() {}

pub fn make_m_homogeneous_partition() {}

pub fn make_plp_homogeneous_partition() {}

// }}}

// {{{ SolveError

pub enum SolverError {
    /// Dimensions of inputs do not match.
    DimensionMismatch = -1,
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

impl From<i32> for SolverError {
    fn from(flag: i32) -> Self {
        match flag {
            -1 => SolverError::DimensionMismatch,
            -2 => SolverError::NegativePower,
            -3 => SolverError::ConstantEquation,
            -4 => SolverError::InconsistentPartitionSize,
            -5 => SolverError::RepeatedPartition,
            -6 => SolverError::InconsistentRecall,
            -7 => SolverError::InsufficientScaleFactors,
            _ => panic!("Unknown SolverError value: {}", flag),
        }
    }
}

pub enum PathTrackingError {
    /// The i-th homotopy path should be re-tracked.
    PathRetracked = -2,
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

impl From<i32> for PathTrackingError {
    fn from(flag: i32) -> Self {
        match flag {
            -2 => PathTrackingError::PathRetracked,
            2 => PathTrackingError::TrackingToleranceFailed,
            3 => PathTrackingError::MaximumStepsExceeded,
            4 => PathTrackingError::BadJacobian,
            5 => PathTrackingError::ZeroCurveLost,
            6 => PathTrackingError::NewtonConvergenceFailed,
            7 => PathTrackingError::RootSearchFailed,
            _ => panic!("Unknown PathTrackingStatus value: {}", flag),
        }
    }
}

pub type PathTrackingResult = Result<u32, PathTrackingError>;

// }}}

// {{{ Solve

// }}}
