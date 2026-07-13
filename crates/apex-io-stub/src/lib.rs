// Copyright The SfM Tool Authors
// SPDX-License-Identifier: Apache-2.0

//! Evaluation stub for `apex-io` (see Cargo.toml for why it exists).
//!
//! Provides exactly the names `apex-solver` 1.3.0's library code touches:
//! `IoError` (wrapped by `ApexSolverError::Io`) and the `BalLoader` /
//! `G2oLoader` / `ToroLoader` / `Graph` re-exports from its lib.rs. None of
//! the dataset-loading functionality is present; the sfmtool evaluation
//! builds its optimization problems directly from arrays.

use std::fmt;

/// Stand-in for `apex_io::IoError`; never constructed by the evaluation.
#[derive(Debug)]
pub struct IoError(pub String);

impl fmt::Display for IoError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "apex-io stub error: {}", self.0)
    }
}

impl std::error::Error for IoError {}

/// Name-only stand-in for `apex_io::BalLoader` (unused by the evaluation).
pub struct BalLoader;

/// Name-only stand-in for `apex_io::G2oLoader` (unused by the evaluation).
pub struct G2oLoader;

/// Name-only stand-in for `apex_io::ToroLoader` (unused by the evaluation).
pub struct ToroLoader;

/// Name-only stand-in for `apex_io::Graph` (unused by the evaluation).
pub struct Graph;
