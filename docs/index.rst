Welcome to the 1QBit OpenQEMIST documentation!
==============================================

.. toctree::
   :maxdepth: 4
   :caption: Contents:

   installation
   refman/modules
   dmet_microsoft_qsharp
   vqe_microsoft_qsharp
   dmet_multisolver
   rigetti_example
   release_notes


Harnessing the combined power of emerging quantum computing technologies and
state-of-the-art classical techniques, QEMIST is 1QBit’s innovative solution to
a fundamental and intractable problem in chemistry: ab initio simulation of molecules. 
 
Ab initio simulation starts from first-principle quantum mechanical theory,
describing the electronic structure of a molecule without making any empirical
assumptions. The accurate prediction of the electronic structure of a molecule
is key to the design of new materials, such as drug compounds and catalyst
molecules, by helping to anticipate a material’s properties before its synthesis
in the lab. However, obtaining this information using classical computers is
computationally intensive, and the resources required for an exact solution
scale exponentially with the size of the problem. Attempts to provide approximate
approaches to this problem on classical computers have been to date either limited
to small-sized systems or compromising on the accuracy of the simulation.

QEMIST is designed to enable the accurate calculation of molecular properties by
leveraging advanced problem decomposition (PD) techniques and quantum computing.
The variety of PD techniques implemented in QEMIST enables massively parallel
simulations by breaking down a computational chemistry task into smaller,
independent subproblems. These subproblems can use a combination of interfaces
to various classical and quantum solvers to achieve a higher level of accuracy
for large-scale, practical molecular simulations. QEMIST is developed primarily
in Python, and its API facilitates the use of quantum development environments
at a lower-level.

1QBit has released 1QBit has released OpenQEMIST_, offering the open source
community an entry point to quantum computing and quantum chemistry simulation.
It provides access to a portion of the functionalities of QEMIST, as open source
software under an Apache 2.0 license, for the benefit of the rapidly growing
community of quantum computing researchers and developers.

Both QEMIST and OpenQEMIST are by design agnostic with respect to quantum
development platforms, in alignment with 1QBit’s unique hardware-agnostic approach.
This enables users, developers, and industry experts to harness the power of the
most-advanced computing resources and algorithms, without the need to learn the
intricacies of each individual hardware platform or to maintain a complex infrastructure. 

In addition to hybrid quantum–classical computing frameworks, OpenQEMIST is
equipped with classical electronic structure solvers, such as the coupled-cluster
and configuration interaction methods. Any of these classical solvers can be easily
utilized as a highly accurate solver in PD techniques, such as the density matrix
embedding theory (DMET) framework implemented in QEMIST.

In addition to hybrid quantum–classical computing frameworks, OpenQEMIST is equipped
with classical electronic structure solvers, such as the coupled-cluster and
configuration interaction methods. Any of these classical solvers can be easily
utilized as a highly accurate solver in PD techniques, such as the density
matrix embedding theory (DMET) framework implemented in QEMIST.

With a growing community of quantum computing researchers and industry partners,
1QBit aims to pave the road toward practical quantum-enabled simulation to
further accelerate innovation in materials science and drug discovery.

Unlock the full potential of QEMIST for your particular problem, contact_ us.

.. _OpenQEMIST: https://1qbit.com/qemist

.. _contact: https://1qbit.com/contact-us/

Highlights of the First Version of OpenQEMIST
==============================================

Further details are available in the OpenQEMIST reference manual. Source code
and installation instructions are available on GitHub_.

1. **Density matrix embedding theory (DMET) based VQE:** 1QBit has been exploring
PD techniques, which have the potential to scale up the size of molecules that
can be simulated by reducing the required quantum resources while maintaining
the accuracy of the electronic structure calculation. DMET, a promising PD
technique, divides a molecule into fragments and determines the electronic
structure of each subsystem using a highly accurate calculation method. Examples
of the DMET and DMET-VQE workflows are available in the "DMET Example" section of
this documentation. The DMET module of OpenQEMIST is available on GitHub_ along with an interactive Jupyter notebook_.

2. **Variational Quantum Eigensolver (VQE) simulations:**  VQE is a hybrid
quantum–classical algorithm for performing quantum chemistry simulations that
requires a shallower circuit than conventional algorithms. An example of the VQE
workflow is available in the “VQE Example” section of this documentation. The
VQE module of OpenQEMIST is available on GitHub_ along with an interactive
Jupyter notebook_.

3. **Integration of the Microsoft QDK:**  The above VQE sample demonstrates the
integration of Microsoft Quantum Development Kit (QDK) with QEMIST. For more
information on 1QBit’s collaboration with Microsoft and the release of
OpenQEMIST, please see the Microsoft blog_.

4. **Integration with PySCF:** OpenQEMIST is equipped with an interface to PySCF
a Python-based electronic structure calculation package. OpenQEMIST utilizes
PySCF to generate the second quantized molecular Hamiltonian for quantum-enabled
simulations. The interface to PySCF allows users to access different classical
solvers in QEMIST.

.. _Github: https://github.com/1QB-Information-Technologies/openqemist

.. _notebook: https://github.com/1QB-Information-Technologies/openqemist/tree/master/examples

.. _blog: https://cloudblogs.microsoft.com/quantum/2019/05/06/new-plans-to-open-source-more-of-the-quantum-development-kit/

Contributing to OpenQEMIST
==========================
We welcome contributions to OpenQEMIST! Please open an issue or submit a pull request on GitHub to start the process.

Citing OpenQEMIST
=================
If you use OpenQEMIST in your research, please cite

Takeshi Yamazaki, Shunji Matsuura, Ali Narimani, Anushervon Saidmuradov, and Arman Zaribafiyan "Towards the Practical Application of Near-Term Quantum Computers in Quantum Chemistry Simulations: A Problem Decomposition Approach" Published on arXiv_ on Jun 4, 2018.

.. _arXiv: https://arxiv.org/abs/1806.01305

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
