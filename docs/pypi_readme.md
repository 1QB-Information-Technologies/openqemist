![logo](http://1qbit.com/wp-content/uploads/2019/05/1qbitlogo.png "1QBit is awesome!")

# OpenQEMIST
[![Build Status](https://travis-ci.com/1QB-Information-Technologies/openqemist.svg?token=zt4rNJ8MTUGcpVsToGyy&branch=master)](https://travis-ci.com/1QB-Information-Technologies/openqemist)

Harnessing the combined power of emerging quantum computing technologies and
state-of-the-art classical techniques, the Quantum-Enabled Molecular ab Initio
Simulation Toolkit, or QEMIST, is 1QBit’s innovative solution to a fundamental
and intractable problem in chemistry: ab initio simulation of molecules.

QEMIST is designed to enable the accurate calculation of molecular properties by
leveraging advanced problem decomposition (PD) techniques and quantum computing.
The variety of PD techniques implemented in QEMIST enables massively parallel
simulations by breaking down a computational chemistry task into smaller,
independent subproblems. These subproblems can use a combination of interfaces
to various classical and quantum solvers to achieve a higher level of accuracy
for large-scale, practical molecular simulations.

OpenQEMIST provides access to a portion of the functionalities of QEMIST as
open source software under an Apache 2.0 license. For more information about the
full functionality of QEMIST and to obtain additional information, please
consult our main [product page](https://1qbit.com/qemist).

## Getting started

To get started, install the package, then see the [Jupyter notebooks](https://github.com/1QB-Information-Technologies/openqemist/tree/master/examples)
for example usage.

To use the Rigetti integration, see the [instructions here](http://docs.rigetti.com/en/stable/start.html)
to install and configure the QVM and Quilc compiler. Alternatively, use the
provided docker-compose file to set up the environment.

To use the Microsoft QDK integration follow the [instructions here](https://docs.microsoft.com/en-us/quantum/install-guide/?view=qsharp-preview)
to install the required dependencies.

## Architechture of OpenQEMIST
OpenQEMIST is organized into three layers: problem decomposition, electronic
structure solvers, and hardware backends. The problem decomposition layer is
responsible for splitting the input molecule into smaller subproblems and
treating these using one particular eigenvalue solver (conceivably, fragments
could be treated using multiple solvers), then processing these results into an
overall output energy. Some electronic structure solvers use classical methods,
while others use wrappers over quantum algorithms running on quantum computing
emulators and simulators from quantum platform providers. The quantum solver
backend layer implements a common interface over libraries, emulators, and
simulators of quantum hardware.

As OpenQEMIST includes only a portion of the functionalities incorporated in QEMIST, DMET
is the only problem decomposition technique open sourced in our initial release.
This release includes the Full CI and coupled-cluster with single and double
excitations (CCSD) electronic structure solvers, as well as a quantum electronic
structure solver based on the Variational Quantum Eigensolver (VQE) algorithm.

On the hardware backend, this initial release is integrated with the Microsoft
Quantum Development Kit for running the VQE algorithm.

## Contributing
We welcome contributions to OpenQEMIST! Please open an issue or submit a pull request on GitHub to start the process.

## Citing
If you use OpenQEMIST in your research, please cite

Takeshi Yamazaki, Shunji Matsuura, Ali Narimani, Anushervon Saidmuradov, and Arman Zaribafiyan "Towards the Practical Application of Near-Term Quantum Computers in Quantum Chemistry Simulations: A Problem Decomposition Approach" Published on [arXiv](https://arxiv.org/abs/1806.01305) on Jun 4, 2018.


_Copyright 1QBit 2019. This software is released under the Apache Software License version 2.0._
