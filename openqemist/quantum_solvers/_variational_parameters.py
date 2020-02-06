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

"""Prepare the initial amplitudes for quantum simulation.

The electronic structure calculation employing the
second-order Møller–Plesset perturbation theory (MP2)
to prepare the initial T2 amplitudes is done here.

"""

def mp2_variational_parameters(molecule, mean_field):
    """ Computes and prepares the MP2 inital variational paramerters.

    Compute the inital variational parameters with PySCF MP2 calculation,
    and then reorders the elements into the QEMIST convention.

    MP2 only has doubles (T2) amplitudes, thus the single (T1) amplitudes
    are set to a small non-zero value and added.

    The ordering for QEMIST is single, double (diagonal), double (non-diagonal).

    Args:
        molecule (pyscf.gto.Mole): The molecule to simulate.
        mean_field (pyscf.scf.RHF): The mean field of the molecule.

    Returns:
        list: The initial variational parameters (float64).

    Note:
        This is only provided as a convenience function. Users are responsible
        for ensuring that the parameters returned are appropriate and of correct
        dimenstion.
    """
    from pyscf import mp
    import numpy as np
    import itertools
    mp2_fragment = mp.MP2(mean_field)
    mp2_fragment.verbose = 0
    mp2_correlation_energy, mp2_t2 = mp2_fragment.kernel()
    scf_energy = mean_field.e_tot
    mp2_total_energy = scf_energy + mp2_correlation_energy

    # Reordering the amplitudes into QEMIST convention
    n_spatial_orbitals = len(mean_field.mo_energy)
    n_occupied = int(np.ceil(molecule.nelectron / 2))
    n_virtual = n_spatial_orbitals - n_occupied

    singles = []
    doubles_1 = []
    doubles_2 = []

    icount = 0

    # Get singles and doubles amplitudes associated with one
    # spatial occupied-virtual pair
    for p, q in itertools.product(range(n_virtual), range(n_occupied)):
        virtual_spatial = p
        occupied_spatial = q

        # Get singles amplitude
        # Just get up amplitude, since down should be the same
        singles.append(2.0e-5)

        temp=-mp2_t2[q,q,p,p]/2.0
        if(abs(temp)<1.0e-15):
            doubles_1.append(0.0)
        else:
            doubles_1.append(-mp2_t2[q,q,p,p]/2.0)

    # Get doubles amplitudes associated with two spatial occupied-virtual pairs
    for (p, q), (r, s) in itertools.combinations(
            itertools.product(range(n_virtual), range(n_occupied)), 2):
        # Get indices of spatial orbitals
        virtual_spatial_1 = p
        occupied_spatial_1 = q
        virtual_spatial_2 = r
        occupied_spatial_2 = s

        # Get amplitude
        doubles_2.append(-mp2_t2[occupied_spatial_1,occupied_spatial_2,virtual_spatial_1,virtual_spatial_2])

    return singles + doubles_1 + doubles_2

