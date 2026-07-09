// SPDX-FileCopyrightText: 2024 Alexandru Fikl <alexfikl@gmail.com>
// SPDX-License-Identifier: MIT

//! Rust wrapper around the POLSYS_PLP Fortran library, a globally
//! convergent homotopy continuation solver for polynomial systems.
//!
//! # Example
//!
//! Solve `x^2 - 3x + 2 = 0` (roots 1 and 2):
//!
//! ```
//! use polsys::{PolsysSolver, Polynomial, term};
//!
//! let mut poly = Polynomial::<1>::new(vec![vec![
//!     term([2], 1.0),   //   x^2
//!     term([1], -3.0),  // -3 x
//!     term([0], 2.0),   // +2
//! ]])?;
//!
//! let result = PolsysSolver::new().solve(&mut poly)?;
//! assert_eq!(result.n_roots(), 2);
//! # Ok::<(), Box<dyn std::error::Error>>(())
//! ```

mod bindings;
mod plp;

#[doc(inline)]
pub use crate::plp::*;
