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

from enum import Enum

from ..parametric_quantum_solver import ParametricQuantumSolver

import os
import tempfile
import warnings

import numpy as np

# Import pyscf and functions making use of it
from pyscf import gto, scf
from .integrals_pyscf import compute_integrals_fragment
from .generate_uccsd_operators import count_amplitudes, compute_cluster_operator
from .broombridge_dummy import _dummy_0_2_yaml

class MicrosoftQSharpParametricSolver(ParametricQuantumSolver):
    """Performs an energy estimation for a molecule with a parametric circuit.

    Performs energy estimations for a given molecule and a choice of ansatz
    circuit that is supported.

    Attributes:
        n_samples (int): The number of samples to take from the hardware emulator.
        optimized_amplitudes (list): The optimized amplitudes.
        verbose(bool): Toggles the printing of debug statements.
    """

    class Ansatze(Enum):
        """ Enumeration of the ansatz circuits that are supported."""
        UCCSD = 0

    def __init__(self, ansatz, molecule, mean_field=None, backend_options=None):
        """Initialize the settings for simulation.

        If the mean field is not provided it is automatically calculated.

        Args:
            ansatz (OpenFermionParametricSolver.Ansatze): Ansatz for the quantum solver.
            molecule (pyscf.gto.Mole): The molecule to simulate.
            mean_field (pyscf.scf.RHF): The mean field of the molecule.
            backend_options (dict): Extra parameters that control the behaviour
                of the solver.
        """
        # Import python packages for Microsoft Python interops
        import qsharp
        import qsharp.chemistry as qsharpchem

        assert(isinstance(ansatz, MicrosoftQSharpParametricSolver.Ansatze))
        self.ansatz = ansatz

        self.verbose = False

        # Initialize the number of samples to be used by the MicrosoftQSharp backend
        self.n_samples = 1e18

        # Initialize the amplitudes (parameters to be optimized)
        self.optimized_amplitudes = []

        # Obtain fragment info with PySCF
        # -----------------------------------------

        # Compute mean-field if not provided. Check that it has converged
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

        if not mean_field.converged:
            warnings.warn("MicrosoftQSharpParametricSolver simulating with mean field not converged.",
                    RuntimeWarning)

        self.molecule = molecule
        self.mean_field = mean_field

        self.n_orbitals = len(mean_field.mo_energy)
        self.n_spin_orbitals = 2 * self.n_orbitals
        self.n_electrons = molecule.nelectron
        nuclear_repulsion = mean_field.energy_nuc()

        # Compute and set values of electronic integrals
        # ----------------------------------------------

        # Get data-structure to store problem description
        fd, path = tempfile.mkstemp(suffix='.yaml')
        try:
            # Write the dummp_0.2.yaml file to a temporary file
            with os.fdopen(fd, 'w') as tmp:
                tmp.write(_dummy_0_2_yaml)

            molecular_data = qsharpchem.load_broombridge(path)

            # Compute one and two-electron integrals, store them in the Microsoft data-structure
            integrals_one, integrals_two = compute_integrals_fragment(molecule, mean_field)
            molecular_data.problem_description[0].hamiltonian['OneElectronIntegrals']['Values'] = integrals_one
            molecular_data.problem_description[0].hamiltonian['TwoElectronIntegrals']['Values'] = integrals_two
            molecular_data.problem_description[0].coulomb_repulsion['Value'] = nuclear_repulsion

            # Compute and set values of UCCSD operators
            # -----------------------------------------

            # Generate UCCSD one- and two-body operators
            n_amplitudes = count_amplitudes(self.n_spin_orbitals, self.n_electrons)
            self.amplitude_dimension = n_amplitudes
            amplitudes = 0.01 * np.ones((n_amplitudes), dtype=np.float64)
            self.preferred_var_params = amplitudes
            ref,t = compute_cluster_operator(self.n_spin_orbitals, self.n_electrons, amplitudes)

            # Load a dummy inputstate object from the dummy Broombridge file, and set its values
            self.inputstate = qsharpchem.load_input_state(path, "UCCSD |G>")

            if self.verbose:
                print("inputstate energy :\n", self.inputstate.Energy)
                print("inputstate mcfdata :\n", self.inputstate.MCFData)
                print("inputstate method :\n", self.inputstate.Method)
                print("inputstate scfdata :\n", self.inputstate.SCFData)
                print("inputstate uccdata :\n", self.inputstate.UCCData, "\n\n\n")

            self.inputstate.UCCData['Reference'] = ref
            self.inputstate.UCCData['Excitations'] = t

            if self.verbose:
                print("inputstate :\n", self.inputstate.UCCData)
                print("------------\n")

            # Generate Fermionic and then qubit Hamiltonians
            # ----------------------------------------------

            # C# Chemistry library : Compute fermionic Hamiltonian
            self.ferm_hamiltonian = molecular_data.problem_description[0].load_fermion_hamiltonian()
            if self.verbose:
                print("ferm_hamiltonian:\n", self.ferm_hamiltonian.terms)
                print("------------\n")

            # C# Chemistry library : Compute the Pauli Hamiltonian using the Jordan-Wigner transform
            self.jw_hamiltonian = qsharpchem.encode(self.ferm_hamiltonian, self.inputstate)
            if self.verbose:
                print("jw_hamiltonian ::", self.jw_hamiltonian)
                print("------------\n")

            # Retrieve energy offset and number of qubits
            self.n_qubits = self.jw_hamiltonian[0]
            self.energy_offset = self.jw_hamiltonian[3]
        finally:
            # Cleanup the temp file
            os.remove(path)


    def simulate(self, amplitudes):
        """Perform the simulation for the molecule.

        If the mean field is not provided it is automatically calculated.

        Args:
            amplitudes (list): The initial amplitudes (float64).

        Returns:
            float64: The total energy (energy).
        """
        import qsharp
        import qsharp.chemistry as qsharpchem
        # Import the "EstimateEnergy" Q# operation from the QDK Chemistry library
        estimate_energy = qsharp.QSharpCallable("Microsoft.Quantum.Chemistry.JordanWigner.VQE.EstimateEnergy", "")

        # Test if right number of amplitudes have been passed
        if len(amplitudes) != self.amplitude_dimension:
            raise ValueError("Incorrect dimension for amplitude list.")

        amplitudes = list(amplitudes)
        self.jw_hamiltonian = self._set_amplitudes(amplitudes, self.jw_hamiltonian)

        # Compute energy
        energy = estimate_energy.simulate(jwHamiltonian=self.jw_hamiltonian, nSamples=self.n_samples)

        # Update optimal amplitudes
        self.optimized_amplitudes = amplitudes

        return energy


    def get_rdm(self):
        """Obtain the RDMs from the optimized amplitudes.

        Obtain the RDMs from the optimized amplitudes by using the
        same function for energy evaluation.
        The RDMs are computed by using each fermionic Hamiltonian term,
        transforming them and computing the elements one-by-one.
        Note that the Hamiltonian coefficients will not be multiplied
        as in the energy evaluation.
        The first element of the Hamiltonian is the nuclear repulsion
        energy term, not the Hamiltonian term.

        Returns:
            (numpy.array, numpy.array): One & two-particle RDMs (rdm1_np & rdm2_np, float64).
        """
        import qsharp
        import qsharp.chemistry as qsharpchem

        amplitudes = self.optimized_amplitudes
        one_rdm = np.zeros((self.n_orbitals, self.n_orbitals))
        two_rdm = np.zeros((self.n_orbitals, self.n_orbitals, self.n_orbitals, self.n_orbitals))

        # Loop over all single fermionic hamiltonian term to get RDM values
        all_terms = self.ferm_hamiltonian.terms
        import copy
        fh_copy = copy.deepcopy(self.ferm_hamiltonian)

        for ii in all_terms:
            for jj in ii[1]:
                # Only use a single fermionic term, set its coefficient to 1.
                term_type = ii[0]
                jj = (jj[0], 1.0)
                single_fh = (term_type, [jj])

                fh_copy.terms = [single_fh]
                # Compute qubit Hamiltonian (C# Chemistry library)
                self.jw_hamiltonian = qsharpchem.encode(fh_copy, self.inputstate)

                # Compute RDM value
                RDM_value = self.simulate(amplitudes)

                # Update RDM matrices
                ferm_ops = single_fh[1][0][0][0]
                indices = [ferm_op[1] for ferm_op in ferm_ops]

                # 1-RDM matrix
                if (len(term_type) == 2):
                    i, j = indices[0]//2, indices[1]//2
                    if (i == j):
                        one_rdm[i, j] += RDM_value
                    else:
                        one_rdm[i, j] += RDM_value
                        one_rdm[j, i] += RDM_value

                # 2-RDM matrix (works with Microsoft Chemistry library sign convention)
                elif (len(term_type) == 4):
                    i, j, k, l = indices[0]//2, indices[1]//2, indices[2]//2, indices[3]//2

                    if((indices[0]==indices[3]) and (indices[1]==indices[2])):
                        if((indices[0]%2 == indices[2]%2) and (indices[1]%2 == indices[3]%2)):
                            two_rdm[i,l,j,k] += RDM_value
                            two_rdm[j,k,i,l] += RDM_value
                            two_rdm[i,k,j,l] -= RDM_value
                            two_rdm[j,l,i,k] -= RDM_value
                        else:
                            two_rdm[i,l,j,k] += RDM_value
                            two_rdm[j,k,i,l] += RDM_value
                    else:
                        if((indices[0]%2 == indices[3]%2) and (indices[1]%2 == indices[2]%2)):
                            two_rdm[i,l,j,k] += RDM_value
                            two_rdm[j,k,i,l] += RDM_value
                            two_rdm[l,i,k,j] += RDM_value
                            two_rdm[k,j,l,i] += RDM_value
                            if((indices[0]%2 == indices[2]%2) and (indices[1]%2 == indices[3]%2)):
                                two_rdm[i,k,j,l] -= RDM_value
                                two_rdm[j,l,i,k] -= RDM_value
                                two_rdm[k,i,l,j] -= RDM_value
                                two_rdm[l,j,k,i] -= RDM_value
                        else:
                            two_rdm[i,k,j,l] -= RDM_value
                            two_rdm[j,l,i,k] -= RDM_value
                            two_rdm[k,i,l,j] -= RDM_value
                            two_rdm[l,j,k,i] -= RDM_value

        return (one_rdm, two_rdm)

    def default_initial_var_parameters(self):
        """ Returns initial variational parameters for a VQE simulation.

        Returns initial variational parameters for the circuit that is generated
        for a given ansatz.

        Returns:
            list: Initial parameters.
        """
        if self.ansatz == self.__class__.Ansatze.UCCSD:
            from .._variational_parameters import mp2_variational_parameters
            return mp2_variational_parameters(self.molecule, self.mean_field)
        else:
            raise RuntimeError("Unsupported ansatz for automatic parameter generation")

    def _set_amplitudes(self, amplitudes, jw_hamiltonian):
        """ Update variational parameters stored in the Q# JW data-structure """
        # Unpack data-structure
        a1, a2, input_state, a3 = jw_hamiltonian
        b1, operator = input_state

        # Update the cluster operator with the new variational parameters
        ref, new_operator = compute_cluster_operator(self.n_qubits, self.n_electrons, amplitudes, True, operator)

	# Re-pack the data-structure
        input_state = (b1, new_operator)
        jw_hamiltonian = (a1, a2, input_state, a3)
        return jw_hamiltonian
