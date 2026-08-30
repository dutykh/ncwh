# Quasinormal modes of noncommutative geometry-inspired wormholes

[![Class. Quantum Grav.](https://img.shields.io/badge/Class.%20Quantum%20Grav.-42%20%282025%29%20085003-005a9c)](https://iopscience.iop.org/article/10.1088/1361-6382/adc2c6)
[![DOI](https://img.shields.io/badge/DOI-10.1088%2F1361--6382%2Fadc2c6-1a7f37)](https://doi.org/10.1088/1361-6382/adc2c6)
[![arXiv](https://img.shields.io/badge/arXiv-2504.02370-b31b1b)](https://arxiv.org/abs/2504.02370)
[![License: LGPL v2.1](https://img.shields.io/badge/License-LGPL%20v2.1-4c1)](LICENSE)
[![Maple](https://img.shields.io/badge/Maple-18%2B-ff6f00)](https://www.maplesoft.com/products/Maple/)
[![MATLAB](https://img.shields.io/badge/MATLAB-Advanpix%20MCT-0076a8)](https://www.advanpix.com/)
[![QNMs Hall of Fame](https://img.shields.io/badge/QNMs-Hall%20of%20Fame-8a2be2)](https://qnms.denys-dutykh.com/)
[![GitHub](https://img.shields.io/badge/GitHub-dutykh%2Fncwh-181717?logo=github)](https://github.com/dutykh/ncwh)
[![Last commit](https://img.shields.io/github/last-commit/dutykh/ncwh)](https://github.com/dutykh/ncwh/commits/main)

**Authors**: Dr. Davide Batic, Dr. Denys Dutykh, and Ms. Joudy Jamal Beek
(Khalifa University of Science and Technology, Abu Dhabi, UAE)

Computational supplemental material for the study of quasinormal modes (QNMs) of
noncommutative geometry-inspired wormholes, computed with a Chebyshev collocation
spectral method coupled to a multiprecision polynomial eigenvalue solver.

## Contents

- [Overview](#overview)
- [The paper](#the-paper)
- [QNMs Hall of Fame](#qnms-hall-of-fame)
- [The method in a nutshell](#the-method-in-a-nutshell)
- [What is in this repository](#what-is-in-this-repository)
- [Requirements](#requirements)
- [Quick start](#quick-start)
- [Reproducing the classical limit](#reproducing-the-classical-limit)
- [Practical notes](#practical-notes)
- [How to cite](#how-to-cite)
- [License](#license)
- [Authors and contact](#authors-and-contact)

## Overview

Noncommutative geometry replaces the point-like mass distribution of the classical
Schwarzschild solution by a minimal-width Gaussian, which regularises the central
singularity and, for a suitable range of the rescaled mass, opens a traversable
throat. The resulting spacetime supports a discrete set of damped oscillations, the
quasinormal modes, whose complex frequencies `ω` carry the fingerprint of the
underlying geometry.

This repository contains everything needed to reproduce the frequencies reported in
the accompanying paper, for scalar, electromagnetic and vector-type gravitational
perturbations:

1. the Maple worksheets that carry out the symbolic reduction of the perturbation
   equation and the boundary analysis at the throat and at spatial infinity;
2. the Maple procedures that assemble the spectral matrices at a prescribed
   resolution;
3. the Matlab driver that solves the resulting polynomial eigenvalue problem in
   multiprecision arithmetic and displays the spectrum.

The classical Schwarzschild wormhole is treated by the same pipeline and serves as
the reference case: as the rescaled mass parameter grows, the noncommutative
frequencies converge to the classical ones, which is the accuracy test reported in
the paper.

## The paper

> Davide Batic, Denys Dutykh and Joudy Jamal Beek,
> *A spectral approach for quasinormal frequencies of noncommutative
> geometry-inspired wormholes*,
> **Classical and Quantum Gravity**, Volume 42, Number 8, 085003 (2025).

| Resource | Link |
|:--|:--|
| Published version | [iopscience.iop.org/article/10.1088/1361-6382/adc2c6](https://iopscience.iop.org/article/10.1088/1361-6382/adc2c6) |
| DOI | [10.1088/1361-6382/adc2c6](https://doi.org/10.1088/1361-6382/adc2c6) |
| Preprint | [arXiv:2504.02370](https://arxiv.org/abs/2504.02370) [gr-qc] |
| Code | [github.com/dutykh/ncwh](https://github.com/dutykh/ncwh) |

Citation line: Davide Batic *et al* 2025 *Class. Quantum Grav.* **42** 085003.
Published 2 April 2025.

## QNMs Hall of Fame

This repository is one instalment of a longer series. All of it lives at

**[qnms.denys-dutykh.com](https://qnms.denys-dutykh.com/)**

the *Quasi-Normal Modes Hall of Fame*, which gathers, under the banner of
*Celestial Music*, the spectral-method QNM studies carried out at the Mathematics
Department of Khalifa University. Every entry there links its journal version, its
arXiv preprint and its companion code repository, so the site is the natural place to
start if you want the wider picture rather than this single spacetime. The wormhole
paper accompanying the present codes is listed among them.

Sister studies published in the same programme cover Gauss–Bonnet black holes, the
Kazakov–Solodukhin quantum-corrected black hole, Bonanno–Reuter black holes, Planck
stars, noncommutative dirty black holes, noncommutative Schwarzschild black holes,
Lee-Wick black holes, massive static phantom wormholes and Morris–Thorne wormholes.
Several of them reuse, almost verbatim, the pipeline documented below.

## The method in a nutshell

The radial perturbation equation is a second-order ODE on a half-line. Four steps turn
it into a finite matrix problem.

**1. Fix the boundary behaviour.** An ansatz factors out the oscillatory and
divergent behaviour at the throat and at spatial infinity, leaving a function that is
regular at both ends. Two admissible factorisations are used, labelled BC1 and BC2;
they differ by a `√(1 - 1/x)` prefactor, which shifts the indicial exponents and
therefore leads to two different, mutually confirming, discretisations.

**2. Compactify.** The Möbius change of variable

```text
z = 2 / (1 - y)
```

maps the physical domain `z ∈ [1, ∞)` onto `y ∈ [-1, 1]`, with `y = -1` the
throat and `y = +1` spatial infinity.

**3. Collocate.** The unknown is expanded in Chebyshev polynomials of the first kind,

```text
         n-1
P(y)  =   Σ   a_j T_j(y)
         j=0
```

and the transformed ODE is imposed at the `n` Chebyshev roots

```text
y_i = cos( (2i - 1) π / (2n) ),    i = 1, ..., n
```

Because the coefficients of the ODE are polynomials of degree two in the frequency
`ω`, the collocation conditions collect into three `n × n` matrices `M0`, `M1`, `M2`.

**4. Solve.** The frequencies are the roots of the quadratic polynomial eigenvalue
problem

```text
( M0 + i ω M1 + ω² M2 ) a = 0
```

solved with `polyeig` in extended precision. Spectra obtained at several resolutions
are superposed: the eigenvalues that persist are the physical QNMs, the rest are
numerical artefacts of the truncation.

## What is in this repository

The layout is flat. Nine files, in three groups matching the three stages above.

### Naming convention

Two labels are combined in every filename.

- **Model.** Files without an `s` treat the **noncommutative** wormhole (`NCWH`),
  identified by the rescaled mass `mu`, the error function and the Gaussian smearing
  term. Files with the `s` prefix (`sbc1`, `sbc2`) and the worksheets marked
  `FROM THROAT` treat the **classical Schwarzschild** wormhole, whose metric function
  reduces to the closed form `f(x) = 1 - 1/x` with the throat fixed at `x = 1`.
- **Boundary condition.** `BC1` and `BC2` are the two factorisations of step 1 above.
  For the classical worksheets this distinction is carried by the `01-` and `02-`
  numeric prefix instead.

Only the non-extreme case is covered here, as recorded in the `NONEXTREME` filenames
and in the header comments of the assembly procedures.

### Symbolic derivations (Maple worksheets)

Each worksheet reproduces, cell by cell, the reduction of the perturbation equation
to the collocated form, together with the one-sided limits at `y = ±1` that
classify the endpoints and populate the classification tables of the paper. They are
saved as executed worksheets, so the outputs are visible without rerunning them.

| File | Case | Content |
|:--|:--|:--|
| `01-NONEXTREME CASE_NCWH_BC1-SUPPLEMENTAL MATERIAL.mw` | Noncommutative, BC1 | Multiplication of the radial equation by `f(z)`, the BC1 ansatz, the transformed ODE, the definitions of `f(y)` and `V(y)`, the order of the zeroes and poles of the coefficients `S2`, `S1`, `S0`, the auxiliary functions `N1`, `N0` and `C2`, `C1`, `C0`, a consistency check against the assembled operator, and the endpoint limits. |
| `01-NONEXTREME CASE_NCWH_BC2-SUPPLEMENTAL MATERIAL.mw` | Noncommutative, BC2 | The same programme for the ansatz `psi(z) = sqrt(1-1/z) * z^(i*Omega/2) * exp(i*x0*Omega*(z-1)) * Phi(z)`, ending with the definitions of `L00hat` to `L22hat` and their limits at both endpoints. |
| `01-FROM THROAT_SCHWARZSCHILD WORMHOLE SUPPL MATERIAL.mw` | Classical, BC1 | Normal form of the ODE, indicial analysis at the throat and the asymptotic analysis at infinity, the ansatz without the square-root prefactor, the compactification, and the limits of the coefficients. |
| `02-FROM THROAT_SCHWARZSCHILD WORMHOLE SUPPL MATERIAL.mw` | Classical, BC2 | Identical route with the ansatz carrying the extra `sqrt(1-1/x)` factor. |

### Matrix assembly (Maple procedures)

Each file defines a single procedure, `MatrixAssembler`, which performs step 3 and
writes the three matrices to `data/M0_<n>.mat`, `data/M1_<n>.mat` and
`data/M2_<n>.mat` in ASCII MATLAB format.

| File | Case | Signature |
|:--|:--|:--|
| `matrixassembler_bc1.mpl` | Noncommutative, BC1 | `MatrixAssembler(d, n, s, L, mu, p)` |
| `matrixassembler_bc2.mpl` | Noncommutative, BC2 | `MatrixAssembler(d, n, s, L, mu, p)` |
| `matrixassembler_sbc1.mpl` | Classical, BC1 | `MatrixAssembler(d, n, s, L, p)` |
| `matrixassembler_sbc2.mpl` | Classical, BC2 | `MatrixAssembler(d, n, s, L, p)` |

Arguments:

| Name | Meaning |
|:--|:--|
| `d` | number of digits used in the symbolic and floating-point computations |
| `n` | number of Chebyshev modes, equal to the matrix size |
| `s` | spin of the perturbation: `0` scalar, `1` electromagnetic, `2` vector-type gravitational |
| `L` | angular momentum, required to satisfy `L >= s` |
| `mu` | rescaled, dimensionless mass (noncommutative cases only) |
| `p` | absolute path of the directory containing `data/` |

In the noncommutative cases the procedure first locates the outer throat `x0` as the
root of

```text
1 - erf(μ x) / x + (2 μ / √π) exp(-μ² x²) = 0
```

by `fsolve`, reports it, and only then builds the potential and the nine ODE
coefficients. In the classical cases the coefficients are closed-form polynomials in
`y` and no root finding is needed.

### Eigenvalue solver (Matlab)

`chaseQNMs.m` performs step 4. It loops over the resolutions `n = 200, 150, 100`,
sets the working precision to `n` decimal digits, reads the three matrices produced
by Maple, solves `polyeig(M0, 1i*M1, M2)`, and overlays the three spectra in the
complex `ω` plane with distinct markers. Superposition is the convergence
diagnostic: converged frequencies coincide across resolutions, spurious ones migrate.

### Which worksheet goes with which assembler

Nothing in the filenames records this, so it is worth stating explicitly. The final
coefficient definitions of each worksheet are, term by term, the `L00` to `L22`
entries of the corresponding procedure.

| Symbolic derivation | Numerical implementation |
|:--|:--|
| `01-NONEXTREME CASE_NCWH_BC1-SUPPLEMENTAL MATERIAL.mw` | `matrixassembler_bc1.mpl` |
| `01-NONEXTREME CASE_NCWH_BC2-SUPPLEMENTAL MATERIAL.mw` | `matrixassembler_bc2.mpl` |
| `01-FROM THROAT_SCHWARZSCHILD WORMHOLE SUPPL MATERIAL.mw` | `matrixassembler_sbc1.mpl` |
| `02-FROM THROAT_SCHWARZSCHILD WORMHOLE SUPPL MATERIAL.mw` | `matrixassembler_sbc2.mpl` |

## Requirements

| Component | Notes |
|:--|:--|
| Maple | The worksheets are saved in the Maple 18 format, so Maple 18 or any later release will open them. The assembly procedures use only `LinearAlgebra` and the standard library. |
| Matlab | Any recent release. Only `polyeig` and base plotting are used. |
| [Advanpix Multiprecision Computing Toolbox](https://www.advanpix.com/) | Required by `chaseQNMs.m` for the `mp` class. The QNM spectrum of a spectral discretisation of this size is not resolvable in double precision, so this dependency is not optional. The toolbox is commercial. |

## Quick start

**1. Create the output directory.** It is not tracked in the repository:

```bash
git clone https://github.com/dutykh/ncwh.git
cd ncwh
mkdir -p data
```

**2. Assemble the matrices in Maple.** Read the procedure for the case of interest,
then call it once per resolution. For the noncommutative wormhole with the first
boundary condition, a scalar perturbation with `L = 2` and rescaled mass `mu = 1.5`:

```maple
read "matrixassembler_bc1.mpl":
MatrixAssembler(200, 200, 0, 2, 1.5, "/absolute/path/to/ncwh"):
MatrixAssembler(200, 150, 0, 2, 1.5, "/absolute/path/to/ncwh"):
MatrixAssembler(200, 100, 0, 2, 1.5, "/absolute/path/to/ncwh"):
```

The three resolutions match the `list` array at the top of `chaseQNMs.m`. Keep `d`
greater than or equal to the largest `n`, since the Matlab driver sets its working
precision to `n` digits.

For the classical Schwarzschild wormhole the call drops the mass argument:

```maple
read "matrixassembler_sbc1.mpl":
MatrixAssembler(200, 200, 0, 2, "/absolute/path/to/ncwh"):
```

**3. Solve and plot in Matlab.** Point the first `addpath` in `chaseQNMs.m` at your
Advanpix installation, then run:

```matlab
chaseQNMs
```

The figure shows the three spectra superposed. The frequencies that lie on top of one
another across all three resolutions are the converged quasinormal modes.

To sweep spin or angular momentum, change `s` and `L` in step 2 and rerun step 3. To
change the boundary condition or the model, restart Maple and read the corresponding
`.mpl` file (see the note below).

## Reproducing the classical limit

The convergence of the noncommutative frequencies to the classical Schwarzschild ones
at large rescaled mass is the accuracy check reported in the paper, and it is
reproduced entirely with the files here:

1. assemble and solve with `matrixassembler_sbc1.mpl` (or `sbc2`) for a chosen `s`
   and `L`, and record the resulting frequencies;
2. repeat with `matrixassembler_bc1.mpl` (or `bc2`, matching the boundary condition)
   for the same `s` and `L`, with `mu` taken progressively larger;
3. compare. The noncommutative spectrum should approach the classical one as `mu`
   grows.

Running BC1 and BC2 for the same physical configuration gives a second, independent
check: the two discretisations are built from different factorisations and must agree
on the physical frequencies.

## Practical notes

- **The `data/` directory must exist before the first Maple call.** `ExportMatrix`
  will not create it.
- **All four procedures share the name `MatrixAssembler`**, so only one of them can
  be loaded in a Maple session at a time. Restart Maple before switching case.
  The noncommutative procedures also leave the computed throat `x0` in the global
  session, which is a further reason to start from a clean state.
- **The working precision grows with the resolution.** `chaseQNMs.m` calls
  `mp.Digits(n)`, so a run at `n = 200` operates at 200 decimal digits. Assembly in
  Maple must be carried out at least as accurately.
- **The marker list must be long enough.** `chaseQNMs.m` stops with an error if
  `Marks` has fewer entries than `list`.

## How to cite

If these codes contribute to your work, please cite the paper:

> D. Batic, D. Dutykh and J. Jamal Beek, *A spectral approach for quasinormal
> frequencies of noncommutative geometry-inspired wormholes*, Class. Quantum Grav.
> **42** (2025) 085003.

```bibtex
@article{Batic2025ncwh,
  author        = {Batic, Davide and Dutykh, Denys and Jamal Beek, Joudy},
  title         = {A spectral approach for quasinormal frequencies of
                   noncommutative geometry-inspired wormholes},
  journal       = {Classical and Quantum Gravity},
  volume        = {42},
  number        = {8},
  pages         = {085003},
  year          = {2025},
  month         = apr,
  publisher     = {IOP Publishing},
  doi           = {10.1088/1361-6382/adc2c6},
  eprint        = {2504.02370},
  archivePrefix = {arXiv},
  primaryClass  = {gr-qc},
  url           = {https://doi.org/10.1088/1361-6382/adc2c6}
}
```

## License

Released under the GNU Lesser General Public License, version 2.1. The full text is
in [`LICENSE`](LICENSE).

## Authors and contact

| Author | Affiliation | Contact |
|:--|:--|:--|
| Dr. Davide Batic | Mathematics Department, Khalifa University of Science and Technology, Abu Dhabi, UAE | [davide.batic@ku.ac.ae](mailto:davide.batic@ku.ac.ae) |
| Dr. Denys Dutykh | Mathematics Department, Khalifa University of Science and Technology, Abu Dhabi, UAE | [denys.dutykh@ku.ac.ae](mailto:denys.dutykh@ku.ac.ae) |
| Ms. Joudy Feras Jamal Beek | Physics Department, Khalifa University of Science and Technology, Abu Dhabi, UAE | [KU100061683@alumni.ku.ac.ae](mailto:KU100061683@alumni.ku.ac.ae) |

Questions, corrections and reports of successful reuse are welcome. For the broader
series of spectral QNM studies, see the
[QNMs Hall of Fame](https://qnms.denys-dutykh.com/).
