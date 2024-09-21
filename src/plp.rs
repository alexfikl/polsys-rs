// SPDX-FileCopyrightText: 2024 Alexandru Fikl <alexfikl@gmail.com>
// SPDX-License-Identifier: MIT

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
            _ => panic!("Unknown ErrorFlag value: {}", flag),
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
