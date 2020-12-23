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

from pyscf import gto

from openqemist.electronic_structure_solvers import FCISolver, CCSDSolver
from openqemist.problem_decomposition import DMETProblemDecomposition
from openqemist.problem_decomposition.electron_localization import iao_localization, meta_lowdin_localization


H4_RING = """
H   0.7071067811865476   0.0                 0.0
H   0.0                  0.7071067811865476  0.0
H  -1.0071067811865476   0.0                 0.0
H   0.0                 -1.0071067811865476  0.0
"""


class DMETProblemDecompositionMicrosoftQSharpTest(unittest.TestCase):

    def test_h4ring_vqe_uccsd_microsoft_qsharp_size1(self):
        """
        DMET on H4 ring with fragment size one, using VQE-UCCSD backend
        from Microsoft.
        """

        mol = gto.Mole()
        mol.atom = H4_RING
        mol.basis = "minao"
        mol.charge = 0
        mol.spin = 0
        mol.build()

        # Initialize VQE object with MicrosoftQSharp backend
        from openqemist.quantum_solvers import MicrosoftQSharpParametricSolver
        from openqemist.electronic_structure_solvers import VQESolver
        vqe = VQESolver()
        vqe.hardware_backend_type = MicrosoftQSharpParametricSolver
        vqe.ansatz_type = MicrosoftQSharpParametricSolver.Ansatze.UCCSD

        # Run DMET
        dmet = DMETProblemDecomposition()
        dmet.electron_localization_method = meta_lowdin_localization
        dmet.electronic_structure_solver = vqe
        energy_vqe = dmet.simulate(mol, [1,1,1,1])

        self.assertAlmostEqual(energy_vqe, -1.9916120594, delta=1e-3)
        
    @unittest.skip("DMET couldn't converge after 50 iterations. "
                   "Reminder: Microsoft simulator is not noiseless in practice, as we must draw shots.")
    def test_h4ring_vqe_uccsd_microsoft_qsharp_size2(self):
        """
        DMET on H4 ring with fragment size two, using VQE-UCCSD backend
        from Microsoft.                                                                                                                                                  
        """

        mol = gto.Mole()
        mol.atom = H4_RING
        mol.basis = "minao"
        mol.charge = 0
        mol.spin = 0
        mol.build()

        # Initialize VQE object with MicrosoftQSharp backend                                                                                                              
        from openqemist.quantum_solvers import MicrosoftQSharpParametricSolver
        from openqemist.electronic_structure_solvers import VQESolver
        vqe = VQESolver()
        vqe.hardware_backend_type = MicrosoftQSharpParametricSolver
        vqe.ansatz_type = MicrosoftQSharpParametricSolver.Ansatze.UCCSD

        # Run DMET                                                                                                                                                        
        dmet = DMETProblemDecomposition()
        dmet.electron_localization_method = meta_lowdin_localization
        dmet.electronic_structure_solver = vqe
        energy_vqe = dmet.simulate(mol, [2,2])

        self.assertAlmostEqual(energy_vqe, -1.9916120594, delta=1e-3)

        
if __name__ == "__main__":
    unittest.main()
