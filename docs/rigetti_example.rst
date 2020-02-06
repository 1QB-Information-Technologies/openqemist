OpenQEMIST and Rigetti example
==============================

This notebook shows how OpenQEMIST can be combined with the Rigetti
stack to use the Variational Quantum Eigensolver (VQE) as an electronic
structure solver, and combine it with a problem-decomposition technique
such as Density Matrix Embedding Theory (DMET).

VQE Example
-----------

This tutorial assumes that the user has correctly set up and configured
the OpenQEMIST package. The Variational Quantum Eigensolver
(VQE):math:`^{1,2}` is a hybrid quantum-classical algorithm for
simulating quantum systems. We here focus on VQE within the context of
solving the molecular electronic structure problem for the ground-state
energy of a molecular system. In VQE, we first prepare the trial
wavefunction (quantum state)
:math:`\vert \Psi(\vec{\theta}) \rangle = U(\vec{\theta}) \vert 0 \rangle`
based on an ansatz that depends on :math:`m` parameters defining
:math:`\vec{\theta}=(\theta_1, \theta_2, \ldots, \theta_m)`. The
expectation value of the Hamiltonian (:math:`\hat{H}`),
:math:`\langle \Psi(\vec{\theta}) \vert \hat{H} \vert \Psi(\vec{\theta}) \rangle`,
will then be simulated.

The expectation value can be minimized based on the variational
principle,

.. raw:: latex
`\begin{equation}
E = \min_{\vec{\theta}} \frac{\langle \Psi(\vec{\theta}) \vert \hat{H} \vert \Psi(\vec{\theta}) \rangle}{\langle \Psi(\vec{\theta}) \vert \Psi(\vec{\theta}) \rangle} \geq E_{\text{gs}}\nonumber
\end{equation}`

which ensures that the energy computed will be an upper bound to the
true ground-state energy :math:`E_{\text{gs}}`. This allows us using
classical minimizers to find optimal parameters :math:`\vec{\theta}` for
the ground-state energy :math:`E_{\text{gs}}`.

VQE can be performed using OpenQEMIST in conjuction with the Rigetti
stack for calculating the ground state energy of a molecular system. The
unitary coupled-cluster ansatz can be used to prepare the trial
wavefunction :math:`\vert \Psi(\vec{\theta}) \rangle`. In this notebook,
we will show you an example using a small molecule, the hydrogen
molecule (H:math:`_\text{2}`), for a simulation using VQE.

.. code:: ipython3

    from openqemist.electronic_structure_solvers import VQESolver, FCISolver
    from openqemist.quantum_solvers import RigettiParametricSolver
    
    from pyscf import gto
    
    # Build the molecule
    H2 = [['H', [0.0,   0.0,   0.0]], ['H', [0.0,   0.0,   0.74137727]]]
    mol = gto.Mole()
    mol.atom = H2
    mol.basis = "sto-3g"
    mol.charge = 0
    mol.spin = 0
    mol.build()
    
    # Configure the solver object
    vqe_solver = VQESolver()
    vqe_solver.hardware_backend_type = RigettiParametricSolver
    vqe_solver.ansatz_type = RigettiParametricSolver.Ansatze.UCCSD


We can now simulate the molecule and get its energy.

.. code:: ipython3

    energy_fci = FCISolver().simulate(mol)
    energy_vqe = vqe_solver.simulate(mol)
    
    print("\nFCI energy = ", energy_fci)
    print("VQE energy = ", energy_vqe)


.. parsed-literal::

    FCI energy =  -1.1372704220924401
    VQE energy =  -1.1372704138513532


It is possible to use different initial parameters for the optimization:

.. code:: ipython3

    # Using custom initial parameters
    # Getting the dimension of the initial parameters vector
    num_var_params = vqe_solver.hardware_backend.amplitude_dimension
    # Set the intial parameters for the solver
    vqe_solver.initial_var_params = [0.01 for i in range(num_var_params)]
    
    vqe_solver.simulate(mol)


.. parsed-literal::

    VQE : initial variational parameters: 
     [0.01, 0.01] 
    
.. parsed-literal::

    -1.1372665775495083



Using the QVM shot-based simulator
----------------------------------

To use the QVM, we can use the ``backend_parameters`` attribute of the
``VQESolver`` object. The VQE object then configures the hardware
backend automatically. Because the QVM is slower than the default
wavefunction simulator backend, we specify an optimizer function that
returns after a few iterations, in the interest of showing the usage of
the solver in a reasonable time. See the documentation for more details
about using custom optimizers. This interface is what would also be used
to target a QPU backend in the future.

.. code:: ipython3

    def quick_optimizer(backend, amplitudes):
            from scipy.optimize import minimize
    
            print("Running using custom optimizer.")
            
            # We use a force the optimizer to return after 2 iterations.
            result = minimize(backend, amplitudes, method='COBYLA',
                    options={'disp':True, 'maxiter':2})
    
            return result.fun
    
    vqe = VQESolver()
    vqe.optimizer = quick_optimizer

To use the QVM, we can use the ``backend_parameters`` attribute of the
``VQESolver`` object. The VQE object then configures the hardware
backend automatically. We can then run the simulation with the object.
The number of shots can also be set with this parameter.

Note that because we restricted the optimizer to 2 iterations and
reduced the number of shots, the resulting energy will not be accurate.

.. code:: ipython3

    vqe.hardware_backend_type = RigettiParametricSolver
    vqe.ansatz_type = RigettiParametricSolver.Ansatze.UCCSD
    vqe.backend_parameters = {'backend': '4q-qvm', 'n_shots': 10}
    
    energy = vqe.simulate(mol)
    print("Unconverged QMV energy: ", energy)

.. parsed-literal::
    Unconverged QMV energy:  -1.146711378495224


DMET Example
------------

At the current early stage of quantum hardware, the available
computational resource is yet very limited. Thus, it is still
challenging to perform accurate electronic structure calculations on
actual quantum hardware. Simulation on classical computer requires large
computational cost as well. Therefore, we need to reduce the problem
size while maintaining the accuracy of electronic structure calculation
to solve a problem for small sized molecules to perform quantum
simulations.

Density Matrix Embedding Theory (DMET):math:`^{3,4}` is a powerful
problem decomposition technique to reduce the problem size, while
maintaining the accuracy of the electronic structure calculation. The
DMET method decomposes a molecule into fragments, and each fragment is
treated as an open quantum system that is entangled with each of the
other fragments, all taken together to be that fragment’s surrounding
environment (or “bath”). VQE algorithm can be used with DMET using
OpenQEMIST in conjuction with the Rigetti stack.

In this notebook, we will show you an example of H\ :math:`_\text{4}`
molecule for DMET simulation using VQE as an electronic structure
solver.

.. code:: ipython3

    from openqemist.problem_decomposition import DMETProblemDecomposition
    from openqemist.problem_decomposition.electron_localization import meta_lowdin_localization
    
    H4 = [['H', [0.7071067811865476,   0.0,                 0.0]],
          ['H', [0.0,                  0.7071067811865476,  0.0]],
          ['H', [-1.0071067811865476,  0.0,                 0.0]],
          ['H', [0.0,                 -1.0071067811865476,  0.0]]]
    
    mol = gto.Mole()
    mol.atom = H4
    mol.basis = "minao"
    mol.charge = 0
    mol.spin = 0
    mol.build()
    
    dmet = DMETProblemDecomposition()
    dmet.verbose = True
    dmet.electron_localization_method = meta_lowdin_localization
    # Set the DMET object to use the solver that we configured above
    dmet.electronic_structure_solver = vqe_solver
    energy_vqe = dmet.simulate(mol, [1,1,1,1])
    
    print("The DMET energy is: ", energy_vqe)


.. parsed-literal::

    The DMET energy is:  -1.9916111022655567


References
----------

1. Alberto Peruzzo, Jarrod McClean, Peter Shadbolt, Man-Hong Yung,
   Xiao-Qi Zhou, Peter J. Love, Alán Aspuru-Guzik, and Jeremy L.
   O’Brien, “A variational eigenvalue solver on a photonic quantum
   processor”, Nat. Commun., 5, 4213 (2014).
2. Jarrod R. McClean, Jonathan Romero, Ryan Babbush, and Alán
   Aspuru-Guzik, “The theory of variational hybrid quantum-classical
   algorithms”, New J. Phys., 18, 023023 (2016).
3. Gerald Knizia and Garnet K.-L. Chan, “Density Matrix Embedding: A
   Simple Alternative to Dynamical Mean-Field Theory”, Phys. Rev. Lett.,
   109, 186404 (2012).
4. Sebastian Wouters, Carlos A. Jiménez-Hoyos, Qiming Sun, and Garnet
   K.-L. Chan, “A Practical Guide to Density Matrix Embedding Theory in
   Quantum Chemistry”, J. Chem. Theory Comput., 12, pp. 2706–2719
   (2016).
