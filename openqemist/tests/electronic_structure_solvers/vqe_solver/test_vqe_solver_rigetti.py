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

import unittest, time
from enum import Enum

from pyscf import gto, scf

from openqemist.electronic_structure_solvers import VQESolver, FCISolver
from openqemist.quantum_solvers import RigettiParametricSolver

H2 = """
   H 0.00 0.00 0.0
   H 0.00 0.00 0.74137727
   """

H4 = [['H', [0.7071067811865476,   0.0,                 0.0]],
      ['H', [0.0,                  0.7071067811865476,  0.0]],
      ['H', [-1.0071067811865476,  0.0,                 0.0]],
      ['H', [0.0,                 -1.0071067811865476,  0.0]]]

LiH = [['H', [0.0,   0.0,   0.0]], ['Li', [0.0,   0.0,   1.0]]]


class VQESolverTest(unittest.TestCase):

    def setUp(self):
        self.tstart = time.time()

    def tearDown(self):
        self.tstop = time.time()
        print("Time elapsed in test ", self.id(), " :: ", "{:0.2f}".format(self.tstop-self.tstart), " seconds.")

    def test_h2_sto3g(self):
        """ Test the converged energye of VQE """

        mol = gto.Mole()
        mol.atom = H2
        mol.basis = "sto-3g"
        mol.charge = 0
        mol.spin = 0
        mol.build()

        solver = VQESolver()
        solver.hardware_backend_type = RigettiParametricSolver
        solver.ansatz_type = RigettiParametricSolver.Ansatze.UCCSD

        energy = solver.simulate(mol)
        self.assertAlmostEqual(energy, -1.1372704178510415, delta=1e-3)

    def test_h2_321g(self):
        """ Test the converged energye of VQE """

        mol = gto.Mole()
        mol.atom = H2
        mol.basis = "3-21g"
        mol.charge = 0
        mol.spin = 0
        mol.build()

        solver = VQESolver()
        solver.hardware_backend_type = RigettiParametricSolver
        solver.ansatz_type = RigettiParametricSolver.Ansatze.UCCSD

        solver.initial_var_params = [-0.003013364325590502, 0.005694780930328863, 0.00035386539723315964,
                                     0.04244748156400211, 0.02380828037656643, 0.018789486618923344,
                                     0.0007322518664444671, 0.04376425859874745, 8.503425743572844e-05]
        energy = solver.simulate(mol)
        self.assertAlmostEqual(energy, -1.1478300615818977, delta=1e-3)

    @unittest.skip("This test is too long")
    def test_lih_sto3g(self):
        """ Test the converged energie of VQE """

        mol = gto.Mole()
        mol.atom = LiH
        mol.basis = "sto-3g"
        mol.charge = 0
        mol.spin = 0
        mol.build()

        solver = VQESolver()
        solver.hardware_backend_type = RigettiParametricSolver
        solver.ansatz_type = RigettiParametricSolver.Ansatze.UCCSD

        energy = solver.simulate(mol)
        self.assertAlmostEqual(energy, -7.78445501, delta=1e-3)
        print("LiH optimal energy :: ", energy)


if __name__ == "__main__":
    unittest.main()
