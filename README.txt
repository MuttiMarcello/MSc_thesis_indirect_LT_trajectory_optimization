# Indirect Low-Thrust Solar Electric Propulsion Trajectory Optimization

M.Sc. research thesis project on indirect optimal control of low-thrust Solar Electric Propulsion trajectories for deep-space missions.

## Overview

This repository contains three numerical validation case studies developed as part of a Master of Science research thesis on indirect optimal control for low-thrust mission design.

The case studies are used to generate numerical evidence supporting the proposed analytical formulation, with a focus on low-thrust asteroid rendezvous and Earth-Moon transfer problems.

The work is organized into three main case studies:
1. 2BP asteroid rendezvous with variable SEP parameters
2. 2BP asteroid rendezvous with simplified constant propulsion parameters
3. CR3BP Earth-Moon transfer with simplified constant propulsion parameters

The full methodology and results, including figures and performance maps, are documented in 'docs/thesis_report.pdf'.

## Repository structure

The repository is organized by case study and propulsion modeling assumptions:

- 'src/2BP_VP' - Two-body problem case study with variable SEP parameters
  - 'TO_dyn' - Indirect time-optimal dynamics
  - 'FO_dyn' - Energy-to-fuel optimal homotopy dynamics
  - 'post_processing' - Data analysis and porkchop plot generation tools

- 'src/2BP_CP' - Two-body problem case study with constant propulsion parameters
  - 'TO_dyn'
  - 'FO_dyn'
  - 'post_processing'

- 'src/CR3BP_CP' - CR3BP Earth-Moon transfer case study
  - 'CR3BP_dyn' - Indirect time-optimal dynamics

- 'docs' - Thesis report and post-processing thesis-analyses scripts
  - 'thesis_report.pdf'
  - 'reproducibility.txt' - Detailed instructions for code execution
  - 'docs/2BP_VP/post_processing' - Data analysis and porkchop plot generation tools
  - 'docs/2BP_CP/post_processing'

## Reproducibility and external dependencies

Detailed instructions for correct code execution are provided in 'docs/reproducibility.txt'.

This project relies on external SPICE kernels for ephemerides and constants, not included in this repository.
To run the simulations, the following NAIF kernels are required:
- 'naif0012.tls'
- 'de440.bsp'
- 'gm_de440.tpc'
- '3054374_2000SG344.bsp'
- '3550232_2010UE51.bsp'
- '3568303_2011MD.bsp'
- '20478784_2012UV136.bsp'
- '3702319_2014YD.bsp'

The code calls the kernels through meta-kernel named 'metaker.tm', not included in this repository.
