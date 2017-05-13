# WAVEFORM - a softWAre enVironmEnt For nOnlinear inveRse probleMs

Waveform is a flexible and modular approach to solving PDE-constrained inverse problems that aims to
- accurately reflect the underlying mathematics of the problem
- produce *readable* code that reflects the underlying mathematics
- scales easily from small 2D problems to large-scale 3D problems with *minimal* code modifications
- allow users to integrate their own preconditioners, stencils, linear solvers, etc. easily in to the entire inversion framework

For more details, as well as the full design of the software, see (https://arxiv.org/abs/1703.09268). All functions are written by Curt Da Silva, unless indicated otherwise in the function documentation.

# Requirements
Waveform requires Matlab (tested on R2015b), the Parallel Toolbox (for using parallelism), and the following packages:

- SPOT: A linear-operator toolbox for Matlab, available at https://github.com/mpf/spot
- pSPOT: Parallel extensions to SPOT, available at https://github.com/slimgroup/pSPOT

Optional packages, for running certain examples
- Curvelab, for the Curvelet transform, available at http://www.curvelet.org/

# Installation

At the terminal, type
```
git clone git@github.com:slimgroup/WAVEFORM.git
matlab
```

In the Matlab console, once you have the prerequisite packages installed in your path, type
```
addpath(genpath('WAVEFORM'));
compile_waveform_mex;
```

# Examples
The examples from the attendant paper, as well as some other scripts demonstrating how to use the framework, can be found in the /examples/ subdirectory.

# Citing
If you find this software useful in your own work, please cite

"A Unified 2D/3D Large Scale Software Environment for Nonlinear Inverse Problems", Curt Da Silva and Felix J. Herrmann, 2017.
