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
from pyscf import gto, scf

from openqemist.electronic_structure_solvers import VQESolver, FCISolver
from openqemist.quantum_solvers import QiskitParametricSolver

H2 = [['H', [0.0,   0.0,   0.0]], ['H', [0.0,   0.0,   0.74137727]]]
H4 = [['H', [0.7071067811865476,   0.0,                 0.0]],
      ['H', [0.0,                  0.7071067811865476,  0.0]],
      ['H', [-1.0071067811865476,  0.0,                 0.0]],
      ['H', [0.0,                 -1.0071067811865476,  0.0]]]


class VQESolverTest(unittest.TestCase):

    def test_h2_sto3g(self):
        """ Test the converged energy of VQE with no initial variational parameters
            provided by the user """

        mol = gto.Mole()
        mol.atom = H2
        mol.basis = "sto-3g"
        mol.charge = 0
        mol.spin = 0
        mol.build()

        solver = VQESolver()
        solver.hardware_backend_type = QiskitParametricSolver
        solver.ansatz_type = QiskitParametricSolver.Ansatze.UCCSD

        energy = solver.simulate(mol)
        self.assertAlmostEqual(energy, -1.1372704178510415, delta=1e-3)

    @unittest.skip("Pending faster implementation of get_rdm.")
    def test_h2_321g(self):
        """ Test the converged energy of VQE with initial variational parameters
            provided by the user """

        mol = gto.Mole()
        mol.atom = H2
        mol.basis = "3-21g"
        mol.charge = 0
        mol.spin = 0
        mol.build()

        solver = VQESolver()
        solver.hardware_backend_type = QiskitParametricSolver
        solver.ansatz_type = QiskitParametricSolver.Ansatze.UCCSD

        solver.initial_var_params = [-0.01039529, -0.04685435, -0.01858744, -0.01118045, -0.04674074,
                                     -0.01848484, -0.12702138, -0.0222594, 0.04799664, -0.02237422,
                                     -0.04972733, 0.01266251, 0.04764409, 0.01265669, -0.06169727]

        energy = solver.simulate(mol)
        self.assertAlmostEqual(energy, -1.1478300615818977, delta=1e-3)

    @unittest.skip("Pending faster implementation of get_rdm.")
    def test_h4_sto3g(self):
        """ Test the converged energy of VQE with initial variational parameters provided by the user """

        mol = gto.Mole()
        mol.atom = H4
        mol.basis = "sto3g"
        mol.charge = 0
        mol.spin = 0
        mol.build()

        solver = VQESolver()
        solver.hardware_backend_type = QiskitParametricSolver
        solver.ansatz_type = QiskitParametricSolver.Ansatze.UCCSD

        energy = solver.simulate(mol)
        self.assertAlmostEqual(energy, -1.97784, delta=1e-3)

if __name__ == "__main__":
    unittest.main()
