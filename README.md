![slsqp](media/logo.png)
============

Modern Fortran Edition of the SLSQP Optimizer

### Status

[![Language](https://img.shields.io/badge/-Fortran-734f96?logo=fortran&logoColor=white)](https://github.com/topics/fortran)
[![GitHub release](https://img.shields.io/github/release/jacobwilliams/slsqp.svg)](https://github.com/jacobwilliams/slsqp/releases/latest)
[![Build Status](https://github.com/jacobwilliams/slsqp/actions/workflows/CI.yml/badge.svg)](https://github.com/jacobwilliams/slsqp/actions)
[![codecov](https://codecov.io/gh/jacobwilliams/slsqp/branch/master/graph/badge.svg)](https://codecov.io/gh/jacobwilliams/slsqp)
[![last-commit](https://img.shields.io/github/last-commit/jacobwilliams/slsqp)](https://github.com/jacobwilliams/slsqp/commits/master)

### Description

This is an updated version of the SLSQP nonlinear constrained optimization code. It can be used to solve nonlinear programming problems that seek to minimize a scalar performance index subject to nonlinear equality and inequality constraints as well as bounds on the variables.

Updates to the original code include:

* It has been translated into free-form source.
* It is now thread safe. The original version was not thread safe due to the use of saved variables in one of the subroutines.
* It no longer uses obsolescent and non-standard Fortran features. It should now be 100% standard compliant (Fortran 2008).
* It now has an easy-to-use object-oriented interface. The `slsqp_class` is used for all interactions with the solver. Methods include `initialize()`, `optimize()`, and `destroy()`.
* It includes updated versions of some of the third-party routines used in the original code (BLAS, LINPACK, and NNLS).
* Some new features were added to support printing error  messages and reporting iterations to the user.
* The user can now specify the max and min `alpha` to use during the line search.
* The user can supply a routine to compute the gradients of the objective function and constriants, or allow the code to estimate them using finite differences (backward, forward, or central).
* The documentation strings in the code have been converted to [FORD](https://github.com/Fortran-FOSS-Programmers/ford) format, allowing for [nicely formatted documentation](https://jacobwilliams.github.io/slsqp/) to be auto-generated.
* A couple of bug fixes noted elsewhere have been applied.

### License

  * The original sourcecode and the modifications are released under a [permissive BSD-style license](https://github.com/jacobwilliams/slsqp/blob/master/LICENSE).

### Building SLSQP

#### **Fortran Package Manager**

The library can be built with the [Fortran Package Manager](https://github.com/fortran-lang/fpm) using the provided `fpm.toml` file like so:

```bash
fpm build --release
```

By default, the library is built with double precision (`real64`) real values. Explicitly specifying the real kind can be done using the following processor flags:

Preprocessor flag | Kind  | Number of bytes
----------------- | ----- | ---------------
`REAL32`  | `real(kind=real32)`  | 4
`REAL64`  | `real(kind=real64)`  | 8
`REAL128` | `real(kind=real128)` | 16

For example, to build a single precision version of the library, use:

```
fpm build --profile release --flag "-DREAL32"
```

To use SLSQP within your fpm project, add the following to your `fpm.toml` file:

```toml
[dependencies]
slsqp = { git="https://github.com/jacobwilliams/slsqp.git" }
```

or, to use a specific version:
```toml
[dependencies]
slsqp = { git="https://github.com/jacobwilliams/slsqp.git", tag = "1.3.0" }
```

### Development

  * Development continues on [GitHub](https://github.com/jacobwilliams/slsqp).

### Documentation

  The latest API documentation can be found [here](https://jacobwilliams.github.io/slsqp/). This was generated from the source code using [FORD](https://github.com/Fortran-FOSS-Programmers/ford).

### References

* [Original sourcecode at NETLIB](http://www.netlib.org/toms/733)
* D. Kraft, "[A software package for sequential quadratic programming](https://degenerateconic.com/uploads/2018/03/DFVLR_FB_88_28.pdf)",
  Technical Report DFVLR-FB 88-28, Institut für Dynamik der Flugsysteme,
  Oberpfaffenhofen, July 1988.
* D. Kraft, "[Algorithm 733: TOMP--Fortran modules for optimal control calculations](http://dl.acm.org/citation.cfm?id=192124),"
  ACM Transactions on Mathematical Software, Vol. 20, No. 3, p. 262-281 (1994).
* C. L. Lawson, R. J. Hanson, "Solving Least Squares Problems", Prentice Hall, 1974. (Revised 1995 edition)

