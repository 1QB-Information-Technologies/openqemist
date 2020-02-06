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
Perform quantum simulation using Rigetti's stack, using packages such as pyquil and forestopenfermion.
"""

from enum import Enum
import numpy as np

from pyscf import ao2mo, mp, scf
from openfermion.transforms import jordan_wigner
from openfermion.hamiltonians import MolecularData
from openfermionprojectq import uccsd_trotter_engine, uccsd_singlet_evolution
from openfermion.utils import uccsd_singlet_paramsize, uccsd_singlet_generator

from pyquil.quil import Program, Pragma
from pyquil.paulis import exponentiate, commuting_sets, sX, sY, sZ, sI
from pyquil.gates import *
from pyquil.api import WavefunctionSimulator, get_qc
from pyquil.experiment import ExperimentSetting, TomographyExperiment, TensorProductState
from forestopenfermion import qubitop_to_pyquilpauli

from ..parametric_quantum_solver import ParametricQuantumSolver


class RigettiParametricSolver(ParametricQuantumSolver):
    """Performs an energy estimation for a molecule with a parametric circuit.

    Performs energy estimations for a given molecule and a choice of ansatz
    circuit that is supported.

    Uses the CCSD method to solve the electronic structure problem.
    PySCF program will be utilized.
    Users can also provide a function that takes a `pyscf.gto.Mole`
    as its first argument and `pyscf.scf.RHF` as its second.

    Attributes:
        optimized_amplitudes (list): The optimized UCCSD amplitudes.
        of_mole (openfermion.hamiltonian.MolecularData): Molecular Data in Openfermion.
        f_hamiltonian (openfermion.ops.InteractionOperator): Fermionic Hamiltonian.
        qubit_hamiltonian (openfermion.transforms.jordan_wigner): Qubit Hamiltonian.
        n_qubits (int): Number of qubits.
    """


    class Ansatze(Enum):
        """ Enumeration of the ansatz circuits that are supported."""
        UCCSD = 0


    def __init__(self, ansatz, molecule, mean_field = None,
                 backend_options = {"backend":"wavefunction_simulator"}):
        """Initialize the settings for simulation.

        If the mean field is not provided it is automatically calculated.

        Args:
            ansatz (OpenFermionParametricSolver.Ansatze): Ansatz for the quantum solver.
            molecule (pyscf.gto.Mole): The molecule to simulate.
            mean_field (pyscf.scf.RHF): The mean field of the molecule.
        """

        # Check the ansatz
        assert(isinstance(ansatz, RigettiParametricSolver.Ansatze))
        self.ansatz = ansatz

        # Calculate the mean field if the user has not already done it.
        if not mean_field:
            mean_field = scf.RHF(molecule)
            mean_field.verbose = 0
            mean_field.scf()

            if not mean_field.converged:
                orb_temp = mean_field.mo_coeff
                occ_temp = mean_field.mo_occ
                nr = scf.newton(mean_field)
                energy = nr.kernel(orb_temp, occ_temp)
                mean_field = nr

        # Check the convergence of the mean field
        if not mean_field.converged:
            warnings.warn("RigettiParametricSolver simulating with mean field not converged.",
                    RuntimeWarning)

        self.molecule = molecule
        self.mean_field = mean_field

        # Initialize the amplitudes (parameters to be optimized)
        self.optimized_amplitudes = []

        # Set the parameters for Openfermion
        self.of_mole = self._build_of_molecule(molecule, mean_field)

        # Set the fermionic Hamiltonian
        self.f_hamiltonian = self.of_mole.get_molecular_hamiltonian()

        # Transform the fermionic Hamiltonian into qubit Hamiltonian
        self.qubit_hamiltonian = jordan_wigner(self.f_hamiltonian)
        self.qubit_hamiltonian.compress()
        # Also stores the Rigetti/Forest qubit Hamiltonian
        self.forest_qubit_hamiltonian = qubitop_to_pyquilpauli(self.qubit_hamiltonian)

        # Set the dimension of the amplitudes
        if ansatz == RigettiParametricSolver.Ansatze.UCCSD:
            no_occupied = int(np.ceil(molecule.nelectron / 2))
            no_virtual = len(mean_field.mo_energy) - no_occupied
            no_t1 = no_occupied * no_virtual
            no_t2 = no_t1 * (no_t1 + 1) / 2
            self.amplitude_dimension = int(no_t1 + no_t2)

        # Set the number of qubits
        self.n_qubits = self.of_mole.n_qubits

        # Instantiate backend for computation
        self.backend_options = dict()
        if ("backend" in backend_options) and (backend_options["backend"] != "wavefunction_simulator"):
            self.backend_options["backend"] = get_qc(backend_options["backend"])
        else:
            self.backend_options["backend"] = WavefunctionSimulator()
        self.backend_options["n_shots"] = backend_options["n_shots"] if ("n_shots" in backend_options) else 1000


    def simulate(self, amplitudes):
        """Perform the simulation for the molecule.

        Args:
            amplitudes (list): The initial amplitudes (float64).
        Returns:
            float64: The total energy (energy).
        Raise:
            ValueError: If the dimension of the amplitude list is incorrect.
        """
        if len(amplitudes) != self.amplitude_dimension:
            raise ValueError("Incorrect dimension for amplitude list.")

        #Anti-hermitian operator and its qubit form
        generator = uccsd_singlet_generator(amplitudes, self.of_mole.n_qubits, self.of_mole.n_electrons)
        jw_generator = jordan_wigner(generator)
        pyquil_generator = qubitop_to_pyquilpauli(jw_generator)

        p = Program(Pragma('INITIAL_REWIRING', ['"GREEDY"']))
        # Set initial wavefunction (Hartree-Fock)
        for i in range(self.of_mole.n_electrons):
            p.inst(X(i))
        # Trotterization (unitary for UCCSD state preparation)
        for term in pyquil_generator.terms :
            term.coefficient = np.imag(term.coefficient)
            p += exponentiate(term)
        p.wrap_in_numshots_loop(self.backend_options["n_shots"])

        #  Do not simulate if no operator was passed
        if len(self.qubit_hamiltonian.terms) == 0:
            return 0.
        else:
            # Run computation using the right backend
            if isinstance(self.backend_options["backend"], WavefunctionSimulator):
                energy = self.backend_options["backend"].expectation(prep_prog=p, pauli_terms=self.forest_qubit_hamiltonian)
            else:
                # Set up experiment, each setting corresponds to a particular measurement basis
                settings = [ExperimentSetting(in_state=TensorProductState(), out_operator=forest_term) for forest_term in self.forest_qubit_hamiltonian.terms]
                experiment = TomographyExperiment(settings=settings, program=p)
                print(experiment, "\n")
                results = self.backend_options["backend"].experiment(experiment)

                energy = 0.
                coefficients = [forest_term.coefficient for forest_term in self.forest_qubit_hamiltonian.terms]
                for i in range(len(results)):
                    energy += results[i].expectation * coefficients[i]

            energy = np.real(energy)

        # Save the amplitudes so we have the optimal ones for RDM calculation
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

        Raises:
            RuntimeError: If no simulation has been run.
        """
        if len(self.optimized_amplitudes) == 0:
            raise RuntimeError("Cannot retrieve RDM because no simulation has been run.")

        # Save our accurate hamiltonian
        tmp_hamiltonian = self.qubit_hamiltonian

        # Initialize the RDM arrays
        rdm1_np=np.zeros((self.of_mole.n_orbitals,)*2)
        rdm2_np=np.zeros((self.of_mole.n_orbitals,)*4)

        # Loop over each element of Hamiltonian (non-zero value)
        for ikey,key in enumerate(self.f_hamiltonian):
            length=len(key)
            # Treat one-body and two-body term accordingly
            if(length==2):
                pele, qele = int(key[0][0]), int(key[1][0])
                iele, jele = pele//2, qele//2
            if(length==4):
                pele, qele, rele, sele = int(key[0][0]), int(key[1][0]), int(key[2][0]), int(key[3][0])
                iele, jele, kele, lele = pele//2, qele//2, rele//2, sele//2

            # Select the Hamiltonian element (Set coefficient to one)
            hamiltonian_temp = self.of_mole.get_molecular_hamiltonian()
            for jkey,key2 in enumerate(hamiltonian_temp):
                hamiltonian_temp[key2] = 1. if (key == key2 and ikey !=0) else 0.

            # Qubitize the element
            qubit_hamiltonian2 = jordan_wigner(hamiltonian_temp)
            qubit_hamiltonian2.compress()

            # Overwrite with the temp hamiltonian
            self.qubit_hamiltonian = qubit_hamiltonian2
            self.forest_qubit_hamiltonian = qubitop_to_pyquilpauli(self.qubit_hamiltonian)

            # Calculate the energy with the temp hamiltonian
            opt_energy2 = self.simulate(self.optimized_amplitudes)

            # Put the values in np arrays (differentiate 1- and 2-RDM)
            if(length==2):
                rdm1_np[iele,jele] = rdm1_np[iele,jele] + opt_energy2
            elif(length==4):
                if((iele!=lele) or (jele!=kele)):
                    rdm2_np[lele,iele,kele,jele] = rdm2_np[lele,iele,kele,jele] + 0.5 * opt_energy2
                    rdm2_np[iele,lele,jele,kele] = rdm2_np[iele,lele,jele,kele] + 0.5 * opt_energy2
                else:
                    rdm2_np[iele,lele,jele,kele] = rdm2_np[iele,lele,jele,kele] + opt_energy2

        # Restore the accurate hamiltonian
        self.qubit_hamiltonian = tmp_hamiltonian

        return rdm1_np, rdm2_np

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

    def _build_of_molecule(self, molecule, mean_field):
        """Initialize the instance of Openfermion MolecularData class.

        Interface the pyscf and Openfermion data.
        `pyscf.ao2mo` is used to transform the AO integrals into
        the MO integrals.

        Args:
            molecule (pyscf.gto.Mole): The molecule to simulate.
            mean_field (pyscf.scf.RHF): The mean field of the molecule.

        Returns:
            openfermion.hamiltonian.MolecularData: Molecular Data in Openfermion (of_mole).
        """
        of_mole = MolecularData(geometry=molecule.atom, basis=molecule.basis,
                multiplicity=molecule.spin + 1)

        of_mole.mf = mean_field
        of_mole.mol = molecule
        of_mole.n_atoms = molecule.natm
        of_mole.atoms = [row[0] for row in molecule.atom],
        of_mole.protons = 0
        of_mole.nuclear_repulsion = molecule.energy_nuc()
        of_mole.charge = molecule.charge
        of_mole.n_electrons = molecule.nelectron
        of_mole.n_orbitals = len(mean_field.mo_energy)
        of_mole.n_qubits = 2 * of_mole.n_orbitals
        of_mole.hf_energy = mean_field.e_tot
        of_mole.orbital_energies = mean_field.mo_energy
        of_mole.mp2_energy = None
        of_mole.cisd_energy = None
        of_mole.fci_energy = None
        of_mole.ccsd_energy = None
        of_mole.general_calculations = {}
        of_mole._canonical_orbitals = mean_field.mo_coeff
        of_mole._overlap_integrals = mean_field.get_ovlp()
        of_mole.h_core = mean_field.get_hcore()
        of_mole._one_body_integrals = of_mole._canonical_orbitals.T @ of_mole.h_core @ of_mole._canonical_orbitals
        twoint = mean_field._eri
        eri = ao2mo.restore(8, twoint, of_mole.n_orbitals)
        eri = ao2mo.incore.full(eri, of_mole._canonical_orbitals)
        eri = ao2mo.restore(1, eri, of_mole.n_orbitals)
        of_mole._two_body_integrals = np.asarray(eri.transpose(0,2,3,1), order='C')

        return of_mole
