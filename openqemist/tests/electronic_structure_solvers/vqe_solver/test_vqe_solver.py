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

import unittest
from enum import Enum

from pyscf import gto, scf

from openqemist.electronic_structure_solvers import VQESolver
from openqemist.quantum_solvers.parametric_quantum_solver import ParametricQuantumSolver

H2 = """
   H 0.00 0.00 0.0
   H 0.00 0.00 0.74137727
   """

H4 = """
    H   0.7071067811865476   0.0                 0.0
    H   0.0                  0.7071067811865476  0.0
    H  -1.0071067811865476   0.0                 0.0
    H   0.0                 -1.0071067811865476  0.0
"""

class MockQuantumSolver(ParametricQuantumSolver):
    """ Class that mocks out the abstract ParametricQuantumSolver for testing.

    For the purpose of this test, all we need is that we return the same thing,
    causing the optimizer to converge. We assume that the optimizer is working
    correctly.
    """

    class Ansatze(Enum):
        UCCSD = 0

    def __init__(self, ansatz, molecule, mean_field=None, backend_parameters={}):
        self.initial_wavefunction = None
        self.n_qubits = None
        self.amplitude_dimension = None
        self.preferred_var_params = [0.]

    def simulate(self, amplitudes):
        return 3

    def get_rdm(self):
        pass

    def default_initial_var_parameters(self):
        pass

class VQESolverTest(unittest.TestCase):

    def test_mock(self):
        """ Test that the solver returns with the default optimizer from scipy.

        Since simlate is always returning the same number, this should force the
        optimizer to converge.
        """
        mol = gto.Mole()
        mol.atom = H2
        mol.basis = "sto-3g"
        mol.charge = 0
        mol.spin = 0
        mol.build()

        solver = VQESolver()
        solver.hardware_backend_type = MockQuantumSolver
        solver.ansatz_type = MockQuantumSolver.Ansatze.UCCSD

        energy = solver.simulate(mol)

        self.assertEqual(energy, 3)

    def test_get_rdm_no_sim(self):
        """Tests that exception is raised when calling `get_rdm` before `simulate`
        """
        solver = VQESolver()
        solver.hardware_backend_type = MockQuantumSolver
        solver.ansatz_type = MockQuantumSolver.Ansatze.UCCSD

        self.assertRaises(RuntimeError, solver.get_rdm)

if __name__ == "__main__":
    unittest.main()
