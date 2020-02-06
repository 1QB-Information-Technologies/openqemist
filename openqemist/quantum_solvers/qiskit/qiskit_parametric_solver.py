#   Copyright 2019 1QBit
#
#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0
#
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.

"""
    Functions to compute the energy of a system with a parametric circuit
    and extract the RDM values
    with IBM Qiskit package
"""

from enum import Enum
import numpy as np
from pyscf import ao2mo, mp, scf, gto
from itertools import product

# Import qiskit simulator
from qiskit import Aer

# Lib from Qiskit Aqua and Chemistry
from qiskit.aqua import Operator, QuantumInstance
from qiskit.aqua.algorithms import VQE
from qiskit.chemistry.qmolecule import QMolecule
from qiskit.chemistry import FermionicOperator
from qiskit.chemistry.aqua_extensions.components.variational_forms import UCCSD
from qiskit.chemistry.aqua_extensions.components.initial_states import HartreeFock

from ..parametric_quantum_solver import ParametricQuantumSolver

class QiskitParametricSolver(ParametricQuantumSolver):
    """
        Performs an energy estimation for a molecule with a supprted parametric circuit.
        Can be initialized with PySCF objects (molecule and mean-field)
    """

    class Ansatze(Enum):
        """ Enumeration of the ansatz circuits that are supported."""
        UCCSD = 0

    def __init__(self, ansatz, molecule, mean_field=None, backend_options=None):
        """Initialize the settings for simulation.
        If the mean field is not provided it is automatically calculated.

        Args:
            ansatz (QiskitParametricSolver.Ansatze): Ansatz for the quantum solver.
            molecule (pyscf.gto.Mole): The molecule to simulate.
            mean_field (pyscf.scf.RHF): The mean field of the molecule.
            backend_options (dict): Extra parameters that control the behaviour
                of the solver.
        """

        # Check the ansatz
        assert(isinstance(ansatz, QiskitParametricSolver.Ansatze))
        self.ansatz = ansatz

        # Calculate the mean field if the user has not already done it.
        if not mean_field:
            mean_field = scf.RHF(molecule)
            mean_field.verbose = 0
            mean_field.scf()

            if (mean_field.converged == False):
                orb_temp = mean_field.mo_coeff
                occ_temp = mean_field.mo_occ
                nr = scf.newton(mean_field)
                energy = nr.kernel(orb_temp, occ_temp)
                mean_field = nr

        # Check the convergence of the mean field
        if not mean_field.converged:
            warnings.warn("QiskitParametricSolver simulating with mean field not converged.",
                    RuntimeWarning)

        # Initialize molecule quantities
        # ------------------------------
        self.num_orbitals = mean_field.mo_coeff.shape[0]
        self.num_spin_orbitals = self.num_orbitals * 2
        self.num_particles = molecule.nelec[0] + molecule.nelec[1]
        self.nuclear_repulsion_energy = gto.mole.energy_nuc(molecule)

        # Build one body integrals in Qiskit format
        one_body_integrals_pyscf = mean_field.mo_coeff.T @ mean_field.get_hcore() @ mean_field.mo_coeff
        self.one_body_integrals = QMolecule.onee_to_spin(one_body_integrals_pyscf)

        # Build two body integrals in Qiskit format
        eri = ao2mo.incore.full(mean_field._eri, mean_field.mo_coeff, compact=False)
        eri2 = eri.reshape(tuple([self.num_orbitals]*4))
        self.two_body_integrals = QMolecule.twoe_to_spin(eri2)

        # Build Hamiltonians
        # ------------------
        # Set the fermionic Hamiltonian
        self.f_hamiltonian = FermionicOperator(h1=self.one_body_integrals, h2=self.two_body_integrals)

        # Transform the fermionic Hamiltonian into qubit Hamiltonian
        self.map_type = 'jordan_wigner'
        self.qubit_hamiltonian = self.f_hamiltonian.mapping(map_type=self.map_type, threshold=1e-8)
        self.qubit_hamiltonian.chop(10 ** -10)
        self.n_qubits = self.qubit_hamiltonian.num_qubits

        # Qubits, mapping, backend
        backend_opt = ('statevector_simulator', 'matrix')
        # backend_opt = ('qasm_simulator', 'paulis')
        self.backend = Aer.get_backend(backend_opt[0])
        self.backend_opt = backend_opt

        # Set ansatz and initial parameters
        # ---------------------------------

        # Define initial state
        self.init_state = HartreeFock(self.qubit_hamiltonian.num_qubits, self.num_spin_orbitals, self.num_particles,
                                 self.map_type, False)  # No qubit reduction

        # Select ansatz, set the dimension of the amplitudes
        self.var_form = UCCSD(self.qubit_hamiltonian.num_qubits, depth=1,
                         num_orbitals=self.num_spin_orbitals, num_particles=self.num_particles,
                         initial_state=self.init_state, qubit_mapping=self.map_type,
                         two_qubit_reduction=False, num_time_slices=1)

        self.amplitude_dimension = self.var_form._num_parameters
        self.optimized_amplitudes = []

    def simulate(self, var_params):
        """ Evaluate the parameterized circuit for the input amplitudes.

        Args:
            var_params (list): The initial amplitudes (float64).
        Returns:
            float64: The total energy (energy).
        Raise:
            ValueError: If the dimension of the amplitude list is incorrect.
        """
        if len(var_params) != self.amplitude_dimension:
            raise ValueError("Incorrect dimension for amplitude list.")

        # Use the Qiskit VQE class to perform a single energy evaluation
        from qiskit.aqua.components.optimizers import COBYLA
        cobyla = COBYLA(maxiter=0)
        vqe = VQE(self.qubit_hamiltonian, self.var_form, cobyla, initial_point=var_params)
        quantum_instance = QuantumInstance(backend=self.backend, shots=1000000)
        results = vqe.run(quantum_instance)

        energy = results['eigvals'][0] + self.nuclear_repulsion_energy

        # Save the amplitudes so we have the optimal ones for RDM calculation
        self.optimized_amplitudes = var_params

        return energy

    def get_rdm(self):
        """
            Obtain the 1- and 2-RDM matrices for given variational parameters.
            This makes sense for problem decomposition methods if these amplitudes are the ones
            that minimizes energy.

        Returns:
            (numpy.array, numpy.array): One & two-particle RDMs (rdm1_np & rdm2_np, float64).
        Raises:
            RuntimeError: If no simulation has been run before calling this method.
        """

        if len(self.optimized_amplitudes) == 0:
            raise RuntimeError('Cannot retrieve RDM because method "Simulate" needs to run first.')

        # Initialize RDM matrices and other work arrays
        n_orbital = self.num_orbitals
        one_rdm = np.zeros(tuple([n_orbital]*2))
        two_rdm = np.zeros(tuple([n_orbital]*4))

        tmp_h1 = np.zeros(self.one_body_integrals.shape)
        tmp_h2 = np.zeros(self.two_body_integrals.shape)

        # h1 and h2 are the one- and two-body integrals for the whole system
        # They are in spin-orbital basis, seemingly with all alpha orbitals first and then all beta second
        # eg lines and columns refer to 1a, 2a, 3a ... , then 1b, 2b, 3b ....

        # Compute values for 1-RDM matrix
        # -------------------------------
        for mol_coords, _ in np.ndenumerate(one_rdm):

            rdm_contribution = 0.
            one_rdm[mol_coords] = 0.0

            # Find all entries of one-RDM that contributes to the computation of that entry of one-rdm
            for spin_coords in product([mol_coords[0], mol_coords[0] + self.num_spin_orbitals // 2],
                                       [mol_coords[1], mol_coords[1] + self.num_spin_orbitals // 2]):

                # Skip values too close to zero
                if abs(self.one_body_integrals[spin_coords]) < 1e-10:
                    continue

                # Ignore all Fermionic Hamiltonian term except one
                tmp_h1[spin_coords] = 1.
                coeff = -1. if (spin_coords[0] // n_orbital != spin_coords[1] // n_orbital) else 1.

                # Accumulate contribution of the term to RDM value
                tmp_ferOp = FermionicOperator(h1=tmp_h1, h2=tmp_h2)
                tmp_qubitOp = tmp_ferOp.mapping(map_type=self.map_type, threshold=1e-8)
                tmp_qubitOp.chop(10 ** -10)
                if tmp_qubitOp.num_qubits == 0:
                    continue
                self.qubit_hamiltonian = tmp_qubitOp
                ene_temp = self.simulate(self.optimized_amplitudes) - self.nuclear_repulsion_energy
                rdm_contribution += coeff * ene_temp

                # Reset entries of tmp_h1
                tmp_h1[spin_coords] = 0.

            # Write the value to the 1-RDM matrix
            one_rdm[mol_coords] = rdm_contribution

        # Compute values for 2-RDM matrix
        # -------------------------------
        for mol_coords, _ in np.ndenumerate(two_rdm):

            rdm_contribution = 0.
            two_rdm[mol_coords] = 0.0

            # Find all entries of h1 that contributes to the computation of that entry of one-rdm
            for spin_coords in product([mol_coords[0], mol_coords[0] + self.num_spin_orbitals // 2],
                                       [mol_coords[1], mol_coords[1] + self.num_spin_orbitals // 2],
                                       [mol_coords[2], mol_coords[2] + self.num_spin_orbitals // 2],
                                       [mol_coords[3], mol_coords[3] + self.num_spin_orbitals // 2]):

                # Skip values too close to zero
                if abs(self.two_body_integrals[spin_coords]) < 1e-10:
                    continue

                # Set entries to the right coefficient for tmp_h1
                tmp_h2[spin_coords] = 1.

                # Count alphas and betas. If odd, coefficient is -1, else its 1.
                # Or maybe its a quadrant thin? Check with Yukio ?
                n_betas_total = sum([spin_orb // n_orbital for spin_orb in spin_coords])
                if (n_betas_total == 0) or (n_betas_total == 4):
                    coeff = 2.0
                elif n_betas_total == 2:
                    coeff = -1.0 if (spin_coords[0] // n_orbital != spin_coords[1] // n_orbital) else 1.0

                # Accumulate contribution of the term to RDM value
                tmp_ferOp = FermionicOperator(h1=tmp_h1, h2=tmp_h2)
                tmp_qubitOp = tmp_ferOp.mapping(map_type=self.map_type, threshold=1e-8)
                tmp_qubitOp.chop(10 ** -10)
                if tmp_qubitOp.num_qubits == 0:
                    continue

                self.qubit_hamiltonian = tmp_qubitOp
                ene_temp = self.simulate(self.optimized_amplitudes) - self.nuclear_repulsion_energy
                rdm_contribution += coeff * ene_temp

                # Reset entries of tmp_h2
                tmp_h2[spin_coords] = 0.

            # Write the value to the 1-RDM matrix
            two_rdm[mol_coords] = rdm_contribution

        return one_rdm, two_rdm

    def default_initial_var_parameters(self):
        """ Returns initial variational parameters for a VQE simulation.

        Returns initial variational parameters for the circuit that is generated
        for a given ansatz.

        Returns:
            list: Initial parameters.
        """
        if self.ansatz == self.__class__.Ansatze.UCCSD:
            return self.var_form.preferred_init_points
        else:
            raise RuntimeError("Unsupported ansatz for automatic parameter generation")
