# HP-MSF Synchronizability

# UNDER CONSTRUCTION

This repository is currently under active development.

Results, figures, numerical implementations, and documentation will be added
progressively as the different stages of the project are completed. The
organization and content of the repository may therefore change during the
development of the work.

---

## About this repository

This repository contains the numerical implementation and results of a study
of **synchronizability in networks of Hastings–Powell ecological systems**
within the **Master Stability Function (MSF)** framework.

The work starts from the dynamics of the isolated Hastings–Powell system and
then progressively introduces different ecological coupling mechanisms.

The general organization of the project is

```math
\text{Isolated HP system}
\longrightarrow
\text{Diffusive migration}
\longrightarrow
\text{Nonlinear pairwise couplings}
\longrightarrow
\text{Higher-order interactions}.
```

The current development is focused on the first two stages.

---

## 1. Isolated Hastings–Powell system

The first stage is devoted to the isolated Hastings–Powell model.

The objective is to characterize the local dynamics of the system and establish
a reliable dynamical reference before introducing interactions between
ecological patches.

This includes the analysis of the dynamical behavior of the model across its
parameter space and the identification of parameter values located within a
well-characterized chaotic region.

The selected chaotic regime defines the reference trajectory
$\mathbf{x}_s(t)$ on the synchronization manifold used in the subsequent
Master Stability Function analysis.

---

## 2. Diffusive migration

The second stage considers networks of identical Hastings–Powell systems
coupled through simple diffusive migration.

For the local state vector

```math
\mathbf{x}=(x,y,z)^{\mathsf T},
```

the seven possible non-empty combinations of migrating populations are
considered:

- $x$
- $y$
- $z$
- $xy$
- $xz$
- $yz$
- $xyz$

For diffusive coupling, the transverse variational dynamics can be written as

```math
\dot{\boldsymbol{\xi}}
=
\left[
D\mathbf{f}(\mathbf{x}_s(t))
-
rD\mathbf{g}
\right]
\boldsymbol{\xi},
```

where

```math
r=\sigma\lambda
```

is the normalized coupling parameter.

The Master Stability Function is defined by the largest Lyapunov exponent of
the transverse variational system,

```math
\Lambda=\Lambda(r).
```

Regions satisfying

```math
\Lambda(r)<0
```

correspond to transverse stability of the synchronized solution.

For a given network, synchronizability requires all nontrivial transverse modes

```math
r_i=\sigma\lambda_i,
\qquad i=2,\ldots,N,
```

to lie within a region where

```math
\Lambda(r_i)<0.
```

The seven diffusive migration schemes provide the reference cases for the
subsequent analysis of more general ecological coupling mechanisms.

---

## Repository structure

```text
hp-msf-synchronizability/
│
├── 01_isolated_hp/
│   ├── README.md
│   ├── notebooks/
│   ├── data/
│   └── figures/
│
├── 02_diffusive_migration/
│   ├── README.md
│   ├── notebooks/
│   ├── data/
│   └── figures/
│
├── src/
│
├── README.md
└── .gitignore
```

### `01_isolated_hp`

Contains the analysis and numerical characterization of the isolated
Hastings–Powell system.

### `02_diffusive_migration`

Contains the MSF analysis for the seven simple diffusive migration schemes.

### `src`

Contains the common numerical routines used throughout the project, including
the Hastings–Powell dynamics, Jacobians, numerical integration, Lyapunov
exponent calculations, coupling definitions, and MSF routines.

---

## Numerical organization

The repository separates reusable numerical routines from the individual
numerical experiments.

Jupyter notebooks are used for the analysis, visualization, and organization
of numerical results, while common routines are maintained in reusable Python
modules.

Numerical results are stored together with the corresponding data and figures
so that the calculations can be reproduced and traced back to the parameters
and numerical settings used to generate them.

---

## Next stages

After the analysis of simple diffusive migration, the framework will be applied
to more general pairwise ecological coupling mechanisms.

These interactions will include nonlinear and ecologically motivated coupling
protocols associated with processes such as population movement, trophic
interactions, and density-dependent effects.

The subsequent stage will extend the analysis to **higher-order interactions**,
where the dynamics of a focal ecological patch may depend simultaneously on
the states of multiple interacting patches.

---

## Synchronizability and synchronization

**Synchronizability** and **synchronization** are treated as distinct concepts
throughout this repository.

The Master Stability Function characterizes the transverse stability of the
synchronized solution and therefore provides conditions for the
**synchronizability** of the network.

Direct synchronization in the complete coupled dynamical system is a separate
dynamical question and is studied through direct network simulations, which
can also provide an independent validation of the MSF predictions.

---

## Current status

**UNDER CONSTRUCTION**

The repository is being populated progressively as new numerical results are
obtained and validated.

Results, figures, numerical methods, and documentation may be updated during
the development of the project.
