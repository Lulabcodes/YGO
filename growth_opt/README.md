## Coarse-Grained _S. cerevisiae_ Model: Growth Optimization Application

This repository contains the source code and visualization scripts for a coarse-grained, systemic, and dynamic model of _Saccharomyces cerevisiae_. It also includes scripts for reproducing the results presented in our paper:

- Viviana Nguyen, Pu Xue, Yifei Li, Huimin Zhao, and Ting Lu, "Controlling Circuitry Dictates the Growth Optimization of Saccharomyces cerevisiae".

### Files

#### model
Contains the model infrastructure.

#### utils
Includes utilities that specify colors, labels, and other elements.

#### read_parametersv2.m
Model parameters.

#### shared_setup.m
Contains shared simulation setups, such as initial conditions.

#### Figures
Comprises multiple subfolders. Each subfolder contains a MATLAB file that calls the appropriate files from the `Model`, `Utils`, `read_parametersv2.m`, and `Shared_Setup` . These scripts allow you to run simulations and reproduce the results presented in our paper:
- Viviana Nguyen, Pu Xue, Yifei Li, Huimin Zhao, and Ting Lu, "Controlling Circuitry Dictates the Growth Optimization of Saccharomyces cerevisiae".
