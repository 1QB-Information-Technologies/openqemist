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
import numpy as np
from pyscf import gto, scf, mp, ao2mo

# Solve known issue with import triggered by Qiskit on MacOS
import matplotlib
matplotlib.use('TkAgg')

from openqemist.quantum_solvers import QiskitParametricSolver

# References for H2
H2 = [['H', [0.0,   0.0,   0.0]], ['H', [0.0,   0.0,   0.74137727]]]
var_params_H2 = [1.57080667, 1.57080677, 1.45773205]

# References for H4
H4 = [['H', [0.7071067811865476,   0.0,                 0.0]],
      ['H', [0.0,                  0.7071067811865476,  0.0]],
      ['H', [-1.0071067811865476,  0.0,                 0.0]],
      ['H', [0.0,                 -1.0071067811865476,  0.0]]]
var_params_H4 = [ 1.70365728e-07,  3.23361923e-07,  3.02419386e-07,  3.91014755e-07,
              -7.29790417e-08, -5.04036340e-08, -1.59990812e-07,  3.88548623e-07,
              -2.36729681e-02, -6.08235149e-03, -5.84984602e-02,  1.37519834e-01,
              -6.57438101e-03, -7.61316069e-02,  7.11402219e-02,  2.44463543e-02,
              -5.70748861e-02,  7.11323344e-02, -7.00812412e-01, -1.17355797e-02,
              1.68824018e-01,  3.26678809e-02, -1.16455633e-02, -3.04246090e-02,
              8.53925483e-02,  8.53971803e-02]

# References for LiH
LiH = [['H', [0.0,   0.0,   0.0]], ['Li', [0.0,   0.0,   1.0]]]
var_params_LiH = [3.80151865e-04, -1.47450109e-07, -9.94540842e-10, -4.38354455e-04,
              3.03110120e-02,  7.45678067e-08,  6.14755877e-07, -2.70233698e-04,
              3.79328596e-04, -2.09796676e-07,  1.42329546e-07, -4.37315843e-04,
              3.03067527e-02,  4.03294655e-07,  1.13872869e-06, -2.65637049e-04,
              -3.36974356e-03, -5.97450527e-08,  2.16450224e-07, -1.60094182e-03,
              -5.35353398e-04,  1.09121463e-07,  1.85075888e-07,  4.17976351e-03,
              -9.33627490e-08, -1.73397444e-03,  6.14535395e-07,  3.57545585e-07,
              2.61608302e-07, -2.98313516e-03,  1.09748417e-08,  1.45943621e-07,
              3.42396311e-08,  2.35864895e-07, -1.73139999e-03,  3.77089281e-07,
              1.62156379e-08,  1.65632057e-07, -2.98324120e-03, -4.87712296e-08,
              -1.60119321e-03,  2.77927947e-07,  5.14926451e-07, -7.89996341e-04,
              1.46883285e-03,  1.32483826e-07, -1.29768569e-07, -2.96418687e-03,
              -5.35587213e-04, -1.55756668e-07, -2.59200965e-08,  1.46848934e-03,
              -2.22725583e-02,  2.58345645e-07, -1.90302420e-07,  4.32667495e-02,
              -1.95755915e-08, -2.98522226e-03,  1.16721428e-07, -9.90547205e-08,
              -1.62125703e-07, -3.64300682e-02, -3.86611297e-07,  9.54364419e-08,
              5.61002199e-08,  9.10937568e-08, -2.98642982e-03, -1.51121518e-07,
              2.32808711e-07, -2.04874352e-07, -3.64597327e-02, -1.59731636e-07,
              4.19347730e-03, -2.40669131e-07,  2.03306367e-08, -2.96527197e-03,
              4.33626919e-02,  5.40661718e-08, -2.53819532e-08, -9.69127071e-02,
              8.89873288e-08,  1.79085116e-08,  2.74127564e-03,  7.38160055e-09,
              1.30979386e-07, -1.04477106e-07, -1.74214500e-07,  6.40586397e-08,
              2.74198904e-03, -2.53755031e-08,  2.59000867e-08,  9.23596599e-09 ]


def matricize_2rdm(two_rdm, n_electrons, n_orbitals):
    """ Turns the two_rdm tensor into a matrix for test purposes """

    l = 0
    sq = n_orbitals * n_orbitals
    jpqrs = np.zeros((n_orbitals, n_orbitals), dtype=np.int)
    for i in range(n_orbitals):
        for j in range(n_orbitals):
            jpqrs[i, j] = l
            l += 1

    rho = np.zeros((sq, sq))
    for i in range(n_orbitals):
        for j in range(n_orbitals):
            ij = jpqrs[i, j]
            for k in range(n_orbitals):
                for l in range(n_orbitals):
                    kl = jpqrs[k, l]
                    rho[ij, kl] += two_rdm[i, k, j, l]
    return rho


class QiskitParametricSolverParametricSolverTest(unittest.TestCase):

    def test_no_mf_H2(self):
        """ Tests number of variational parameters as well as simulate and get_RDM methods """

        mol = gto.Mole()
        mol.atom = H2
        mol.basis = "sto-3g"
        mol.charge = 0
        mol.spin = 0
        mol.build()

        # Compute mean-field or load pre-computed mean-field object
        # The parametric solver will compute the mean-field automatically if user passes nothing
        mean_field = scf.RHF(mol)
        mean_field.verbose = 0
        mean_field.scf()

        ansatz = QiskitParametricSolver.Ansatze.UCCSD
        solver = QiskitParametricSolver(ansatz, mol, mean_field)

        # Test that the number of variational parameters is as expected
        self.assertEqual(solver.amplitude_dimension, 3)

        # # Test "simulate"
        energy = solver.simulate(var_params_H2)
        self.assertAlmostEqual(energy, -1.13727, delta=1e-5)

        # Compute RDM matrices
        one_rdm, two_rdm = solver.get_rdm()

        # Test traces of matrices
        n_elec, n_orb = mol.nelectron, mol.nao_nr()

        self.assertAlmostEqual(np.trace(one_rdm), n_elec, msg='Trace of one_rdm does not match number of electrons',
                               delta=1e-6)

        rho = matricize_2rdm(two_rdm, n_elec, n_orb)
        self.assertAlmostEqual(np.trace(rho), n_elec * (n_elec - 1),
                               msg='Trace of two_rdm does not match n_elec * (n_elec-1)', delta=1e-6)

    #TODO: Check running time after QEM2-108 and move to benchmarks if needed
    @unittest.skip("This test can take 20 mins")
    def test_no_mf_H4(self):
        """ Tests number of var_params as well as simulate and get_RDM methods """

        mol = gto.Mole()
        mol.atom = H4
        mol.basis = "sto-3g"
        mol.charge = 0
        mol.spin = 0
        mol.build()

        ansatz = QiskitParametricSolver.Ansatze.UCCSD
        solver = QiskitParametricSolver(ansatz, mol)

        # Test that the number of variational parameters is as expected
        self.assertEqual(solver.amplitude_dimension, 26)

        # Test "simulate"
        energy = solver.simulate(var_params_H4)
        self.assertAlmostEqual(energy, -1.97784, delta=1e-5)

        # Compute RDM matrices
        one_rdm, two_rdm = solver.get_rdm()

        # Test traces of matrices
        n_elec, n_orb = mol.nelectron, mol.nao_nr()

        self.assertAlmostEqual(np.trace(one_rdm), n_elec, msg='Trace of one_rdm does not match number of electrons',
                               delta=1e-6)
        rho = matricize_2rdm(two_rdm, n_elec, n_orb)
        self.assertAlmostEqual(np.trace(rho), n_elec * (n_elec - 1),
                               msg='Trace of two_rdm does not match n_elec * (n_elec-1)', delta=1e-6)

    #TODO: Check running time after QEM2-108 and move to benchmarks if needed
    @unittest.skip("This test takes ~4h to perform")
    def test_no_mf_LiH(self):
        """ Tests get_RDM methods: assume energy is correct and reconstruct from RDM """

        mol = gto.Mole()
        mol.atom = LiH
        mol.basis = "sto-3g"
        mol.charge = 0
        mol.spin = 0
        mol.build()

        ansatz = QiskitParametricSolver.Ansatze.UCCSD
        solver = QiskitParametricSolver(ansatz, mol)

        # Test that the number of variational parameters is as expected
        self.assertEqual(solver.amplitude_dimension, 92)

        # Test "simulate"
        energy = solver.simulate(var_params_LiH)
        #self.assertAlmostEqual(energy, -7.78445501, delta=1e-5)
        # energy for FCI LiH: ~ -7.78445501 ?

        # Compute RDM matrices
        one_rdm, two_rdm = solver.get_rdm()

        # Test traces of matrices
        n_elec, n_orb = mol.nelectron, mol.nao_nr()

        self.assertAlmostEqual(np.trace(one_rdm), n_elec, msg='Trace of one_rdm does not match number of electrons',
                               delta=1e-6)
        rho = matricize_2rdm(two_rdm, n_elec, n_orb)
        self.assertAlmostEqual(np.trace(rho), n_elec * (n_elec - 1),
                               msg='Trace of two_rdm does not match n_elec * (n_elec-1)', delta=1e-6)

    def test_simulate_dimension_throw(self):
        """Tests that all the values are set correctly in the constructor."""

        mol = gto.Mole()
        mol.atom = H2
        mol.basis = "sto-3g"
        mol.charge = 0
        mol.spin = 0
        mol.build()

        mean_field = scf.RHF(mol)
        mean_field.verbose = 0
        mean_field.scf()

        ansatz = QiskitParametricSolver.Ansatze.UCCSD
        solver = QiskitParametricSolver(ansatz, mol, mean_field)

        # This should throw
        self.assertRaises(ValueError, solver.simulate, [0])


if __name__ == "__main__":
    unittest.main()
