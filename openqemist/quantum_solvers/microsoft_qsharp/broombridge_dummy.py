_dummy_0_2_yaml = """

"$schema": https://raw.githubusercontent.com/Microsoft/Quantum/master/Chemistry/Schema/broombridge-0.2.schema.json

bibliography:
- {url: 'https://www.nwchem-sw.org'}
format: {version: '0.2'}
generator: {source: nwchem, version: '6.8'}
problem_description:
- basis_set: {name: sto-3g, type: gaussian}
  coulomb_repulsion: {units: hartree, value: 0.0}
  energy_offset: {units: hartree, value: 0.0}
  fci_energy: {lower: 0.0, units: hartree, upper: 0.0, value: 0.0}
  geometry:
    atoms:
    - coords: [0.0, 0.0, 0.0]
      name: X
    coordinate_system: cartesian
    symmetry: c1
    units: angstrom
  hamiltonian:
    one_electron_integrals:
      format: sparse
      units: hartree
      values:
      - [1, 1, 0.0]
    two_electron_integrals:
      format: sparse
      index_convention: mulliken
      units: hartree
      values:
      - [1,1,1,1,0.0]
  initial_state_suggestions:
  - label: "UCCSD |G>"
    method: unitary_coupled_cluster
    cluster_operator:             # Initial state that cluster operator is applied to.
        reference_state: [1.0, (1a)+, (1b)+, '|vacuum>']             # A one-body cluster term is t^{q}_{p} a^\dag_q a_p             # A one-body unitary cluster term is t^{q}_{p}(a^\dag_q a_p- a^\dag_p a_q)
        one_body_amplitudes:             # t^{q}_{p} p q 
            - [999.0, "(1a)+", "(1a)"]
        two_body_amplitudes: 			# t^{pq}_{rs} p q r s 			# If this is a PQQR term, the middle two indices must coincide.
            - [999.0, "(1a)+", "(1b)+", "(1a)", "(1b)"]
  metadata: {molecule_name: dummy}
  n_electrons: 0
  n_orbitals: 0
  scf_energy: {units: hartree, value: 0.0}
  scf_energy_offset: {units: hartree, value: 0.0}

"""