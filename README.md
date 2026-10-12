# Indirect Low-Thrust Solar Electric Propulsion Trajectory Optimization

M.Sc. thesis research project on indirect optimal control of low-thrust Solar Electric Propulsion trajectories for deep-space missions.

## Overview

This repository contains three numerical validation case studies developed as part of a Master of Science research thesis on indirect optimal control for low-thrust mission design.

The case studies are used to generate numerical evidence supporting the proposed analytical formulation, with a focus on low-thrust asteroid rendezvous and Earth-Moon transfer problems.

The work is organized into three main case studies:

1. 2BP asteroid rendezvous with variable SEP parameters
2. 2BP asteroid rendezvous with simplified constant propulsion parameters
3. CR3BP Earth-Moon transfer with simplified constant propulsion parameters

## Results and Validation

Key results:

* Achieved convergence of time-, energy- and fuel-optimal trajectories with final heliocentric errors of < 1e-5 km, < 1e-11 km/s. Results are validated through Hamiltonian testing and comparison with DART Group proprietary trajectory optimization tools.
* The continuation scheme reliably characterizes the solution space. Energy-to-fuel optimal continuation succeeds in more than 95% of cases. In the remaining instances, the continuation parameter epsilon is ~ 0, yielding a negligible error.
* Proved that under specified analytical conditions, optimal costate paths of energy- and fuel-optimal solutions with minimum time-of-flight can be computed through scalar scaling from the associated time-optimal costate path solution.

Representative outputs:

* Time- and fuel-optimal solutions employing low-thrust SEP propulsion
* Low-thrust mission analysis porkchop plots
* Convergence tests of energy- and fuel-optimal solutions with minimum time-of-flight computed through scalar costate path scaling
* Comparison of switching functions of energy- and fuel-optimal solutions obtained through scalar costate path scaling

Representative figures are available in 'results/' (PNG format).
See 'results/results.txt' for figure-by-figure explanations.
The full methodology and results are documented in 'docs/thesis\_report.pdf'.

## Repository structure

The repository is organized by case study and propulsion modeling assumptions:

* 'src/2BP\_VP' - Two-body problem case study with variable SEP parameters

  * 'TO\_dyn' - Indirect time-optimal dynamics
  * 'FO\_dyn' - Energy-to-fuel optimal homotopy dynamics
  * 'post\_processing' - Data analysis and porkchop plot generation tools

* 'src/2BP\_CP' - Two-body problem case study with constant propulsion parameters

  * 'TO\_dyn'
  * 'FO\_dyn'
  * 'post\_processing'

* 'src/CR3BP\_CP' - CR3BP Earth-Moon transfer case study

  * 'CR3BP\_dyn' - Indirect time-optimal dynamics

* 'docs' - Thesis report and post-processing thesis-analyses scripts

  * 'thesis\_report.pdf'
  * 'reproducibility.txt' - Detailed instructions for code execution
  * 'docs/2BP\_VP/post\_processing' - Data analysis and porkchop plot generation tools
  * 'docs/2BP\_CP/post\_processing'

* 'results/' - Key result figures (PNG)
* 'figures/' - Source figures

## Development notes

This repository was uploaded after project completion.
Commit history does not reflect the original development timeline.

## Reproducibility and external dependencies

Detailed instructions for correct code execution are provided in 'docs/reproducibility.txt'.

This project relies on external SPICE kernels for ephemerides and constants, not included in this repository.
To run the simulations, the following NAIF kernels are required:

* 'naif0012.tls'
* 'de440.bsp'
* 'gm\_de440.tpc'
* '3054374\_2000SG344.bsp'
* '3550232\_2010UE51.bsp'
* '3568303\_2011MD.bsp'
* '20478784\_2012UV136.bsp'
* '3702319\_2014YD.bsp'

The code calls the kernels through meta-kernel named 'metaker.tm', not included in this repository.

