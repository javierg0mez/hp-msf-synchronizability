# HP-MSF Synchronizability

# UNDER CONSTRUCTION

This repository is currently under active development.

Results, figures, numerical implementations, and documentation are added
progressively as the different stages of the project are completed. The
organization and content of the repository may therefore change during the
development of the work.

---

## About this repository

This repository contains the numerical implementation and results of a study
of **synchronizability in networks of Hastings-Powell (HP) ecological systems**
within the **Master Stability Function (MSF)** framework.

The work begins with a detailed analysis of the isolated HP
system. This first stage establishes the local chaotic dynamics that will later
be used as the reference dynamics in the study of coupled ecological systems.

The project is organized progressively as

```math
\text{Isolated HP system}
\longrightarrow
\text{Diffusive migration}
\longrightarrow
\text{Nonlinear pairwise couplings}
\longrightarrow
\text{Higher-order interactions}.
```

The current repository contains the numerical analysis of the isolated
HP system.

---

## 1. Isolated HP system

The first stage is devoted to establishing and characterizing the chaotic
dynamics of the isolated HP model before introducing interactions
between ecological patches.

It is divided into two parts.

### 1.1 Chaos characterization

The first part studies the dynamical behavior of the model as the mortality
parameter $d_1$ is varied.

The analysis combines bifurcation information and Lyapunov exponents to
identify parameter values located inside chaotic regions and to distinguish
them from nearby periodic windows.

The numerical work includes global and fine parameter scans, Lyapunov-spectrum
calculations, continuation tests, and additional checks at the selected working
points.

The corresponding notebooks, numerical data, and figures are stored in

```text
01_isolated_hp/01_chaos_characterization/
```

### 1.2 Isolated dynamics

Once suitable chaotic working points have been identified, the second part
examines the dynamics of the isolated system in greater detail.

This includes the evolution of the three populations, phase-space
representations, Lyapunov-exponent behavior, sensitivity to initial
conditions, and comparisons between the selected chaotic regimes.

The corresponding notebook, processed data, and figures are stored in

```text
01_isolated_hp/02_isolated_dynamics/
```

These results provide the dynamical reference for the subsequent MSF analysis.

---

## Repository structure

```text
hp-msf-synchronizability/
│
├── 01_isolated_hp/
│   │
│   ├── README.md
│   │
│   ├── 01_chaos_characterization/
│   │   ├── hp_core.py
│   │   ├── hp_global_d1_scan.ipynb
│   │   ├── hp_fine_d1_scans.ipynb
│   │   ├── hp_working_points_checks.ipynb
│   │   ├── data/
│   │   └── figures/
│   │
│   └── 02_isolated_dynamics/
│       ├── hp_core.py
│       ├── hp_isolated_dynamics.ipynb
│       ├── data/
│       └── figures/
│
├── README.md
└── .gitignore
```

---

## Numerical organization

Jupyter notebooks are used to organize the numerical experiments, analysis,
and visualization of the results.

Numerical outputs are stored separately from the notebooks in `data/`, while
the corresponding graphical results are stored in `figures/`.

The numerical implementation is organized so that the simulations can be
traced back to the model parameters and numerical settings used to generate
the results.

Large intermediate trajectory files that can be regenerated directly from the
numerical implementation are not included in the repository.

---

## Next stages

The next stage will introduce **simple diffusive migration** between identical
HP systems and apply the classical MSF
formalism to the seven possible non-empty combinations of migrating
populations:

- $x$
- $y$
- $z$
- $xy$
- $xz$
- $yz$
- $xyz$

These cases will provide the reference MSF results before moving to more
general nonlinear and ecologically motivated pairwise coupling mechanisms.

The subsequent stage will extend the analysis to **higher-order interactions**.

---

## Current status

**UNDER CONSTRUCTION**

The isolated-system analysis is currently available in the repository.

Additional numerical results, coupling mechanisms, figures, and documentation
will be incorporated progressively as the project develops.
