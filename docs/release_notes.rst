Release Notes
==============

0.1.1
++++++

- Fix to the RDM generation in the MicrosoftQSharpParametericSolver.

0.2.0
+++++

- DMETProblemDecomposition can now use a different ElectronicStructureSolver
  for each fragment. See the documentation or examples for more details.

- Added integration for IBM Qiskit as a new ParametricQuantumSolver.

- Moved the handling of initial parameters for VQE out of the `vqe` module into
  the concrete implementations of the `ParametricQuantumSolver` interface. The
  interface now requires that a `default_initial_var_parameters` function be
  implemented for each concrete implementation.

- The `initial_amplitude_function` functional argument of the `VQESolver` has
  been removed. The solver will use the `default_initial_var_parameters`
  function of the hardware backend if the `initial_var_parameters` are not
  specified. The only way to specify custom values for the intial variational
  parameters is through the `initial_var_parameters` attribute.

- `ParametricQuantumSolver` objects constructors now take an optional dictionary
  mapping strings to values. This can be passed from the VQE caller to further
  configure the hardware backend.

- Added `RigettiParametricSolver` to provide an integration with the Rigetti
  stack.
