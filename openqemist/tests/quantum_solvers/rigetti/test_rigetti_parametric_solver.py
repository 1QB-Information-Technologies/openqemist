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
import time
from pyscf import gto, scf, mp, ao2mo

from openqemist.quantum_solvers import RigettiParametricSolver

# References for H2
H2 = [['H', [0.0,   0.0,   0.0]], ['H', [0.0,   0.0,   0.74137727]]]
var_params_H2 = [5.87296965e-05, 5.66658421e-02]

# References for H4
H4 = [['H', [0.7071067811865476,   0.0,                 0.0]],
      ['H', [0.0,                  0.7071067811865476,  0.0]],
      ['H', [-1.0071067811865476,  0.0,                 0.0]],
      ['H', [0.0,                 -1.0071067811865476,  0.0]]]
var_params_H4 = [-3.96898484e-04,  4.59786847e-05,  3.95285013e-05,  1.85885610e-04,
                 1.05759154e-02,  3.47363359e-01,  3.42657596e-02,  1.45006203e-02,
                 7.43941871e-02,  7.57255601e-03, -1.83407761e-01, -1.03261491e-01,
                 1.34258277e-02, -3.78096407e-02]

# References for LiH
LiH = [['H', [0.0,   0.0,   0.0]], ['Li', [0.0,   0.0,   1.0]]]
var_params_LiH = [-4.04489267e-04, -3.23963592e-02, 5.90161572e-07, -1.35660710e-05,
                  5.89732960e-07, -1.35569584e-05, 3.22762229e-04, -6.12257984e-07,
                  1.68285018e-03, 1.11546116e-02, 8.66964813e-04, 1.84483160e-02,
                  8.61076185e-04, 1.84605771e-02, 3.77549747e-04, 4.82008394e-02,
                  5.74735823e-04, 5.21309956e-05, -2.04050098e-05, 5.21365670e-05,
                  -2.04011517e-05, 1.63227074e-03, -4.23918658e-03, -3.72788448e-06,
                  -2.88652083e-07, -3.73133954e-06, -2.90058581e-07, -1.45093094e-03,
                  -4.36941970e-02, 2.98175992e-03, 3.60631133e-05, -1.18229228e-05,
                  -1.04558703e-05, -1.21734343e-05, -1.17681016e-05, -3.44081538e-06,
                  -8.72733250e-06, -7.14705798e-06, 2.98193961e-03, -1.04616189e-05,
                  -1.21791692e-05, -8.72101499e-06, -7.14230627e-06, 2.55056180e-03]


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


class RigettiParametricSolverParametricSolverTest(unittest.TestCase):

    def setUp(self):
        self.tstart = time.time()


    def tearDown(self):
        self.tstop = time.time()
        print("Time elapsed in test ", self.id(), " :: ", "{:0.2f}".format(self.tstop-self.tstart), " seconds.")


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

        ansatz = RigettiParametricSolver.Ansatze.UCCSD
        backend_options = {"backend": "wavefunction_simulator"}
        solver = RigettiParametricSolver(ansatz, mol, mean_field, backend_options=backend_options)

        # Test that the number of variational parameters is as expected
        self.assertEqual(solver.amplitude_dimension, 2)

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


    def test_no_mf_H2_QVM_QPU(self):
        """ Test to run energy evaluation with the QVM / QPU backend accesses through Rigetti's get_qc """

        mol = gto.Mole()
        mol.atom = H2
        mol.basis = "sto-3g"
        mol.charge = 0
        mol.spin = 0
        mol.build()

        mean_field = scf.RHF(mol)
        mean_field.verbose = 0
        mean_field.scf()

        ansatz = RigettiParametricSolver.Ansatze.UCCSD
        backend_options = {"backend": "4q-qvm", "shots":1000}
        solver = RigettiParametricSolver(ansatz, mol, mean_field, backend_options=backend_options)

        energy = solver.simulate(var_params_H2)
        self.assertAlmostEqual(energy, -1.13727, delta=1e-2)


    @unittest.skip("This test is too long")
    def test_no_mf_H4(self):
        """ Tests number of var_params as well as simulate and get_RDM methods """

        mol = gto.Mole()
        mol.atom = H4
        mol.basis = "sto-3g"
        mol.charge = 0
        mol.spin = 0
        mol.build()

        ansatz = RigettiParametricSolver.Ansatze.UCCSD
        solver = RigettiParametricSolver(ansatz, mol)

        # Test that the number of variational parameters is as expected
        self.assertEqual(solver.amplitude_dimension, 14)

        # Test "simulate"
        energy = solver.simulate(var_params_H4)
        self.assertAlmostEqual(energy, -1.9778, delta=1e-4)

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

        ansatz = RigettiParametricSolver.Ansatze.UCCSD
        solver = RigettiParametricSolver(ansatz, mol, mean_field)

        # This should throw
        self.assertRaises(ValueError, solver.simulate, [0])


if __name__ == "__main__":
    unittest.main()
