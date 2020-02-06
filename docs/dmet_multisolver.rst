
Using multiple types of solvers for DMET fragments
==================================================

This notebook demonstrates how to use a different electronic structure
solver for each framgent that is produced by the DMET problem
decomposition. For details and theory about DMET, see the other examples
provided with the documentation.

We set up the molecule that we use in this example.

.. code:: ipython3

    import py3Dmol
    
    H8_He2 = """
    HE         1.6180339887          0.0000000000          0.0000000000
    HE         1.3090169944          0.9510565163          0.0000000000
    H          0.5000000000          1.5388417686          0.0000000000
    H         -0.5000000000          1.5388417686          0.0000000000
    H         -1.3090169944          0.9510565163          0.0000000000
    H         -1.6180339887          0.0000000000          0.0000000000
    H         -1.3090169944         -0.9510565163          0.0000000000
    H         -0.5000000000         -1.5388417686          0.0000000000
    H          0.5000000000         -1.5388417686          0.0000000000
    H          1.3090169944         -0.9510565163          0.0000000000
    """
    view = py3Dmol.view(width=400,height=400)
    view.addModel("10\n" + H8_He2,'xyz',{'keepH': True})
    view.setStyle({'sphere':{}})
    view.setStyle({'model':0},{'sphere':{'colorscheme':'cyanCarbon','scale':'0.2'}})
    view.zoomTo()
    view.show()
    
    from pyscf import gto
    mol = gto.Mole() # Instantiate the molecule class in PySCF
    mol.atom = H8_He2   # The coordinates of the atoms of the 10-hydrogen-atom ring are defined above
    mol.basis = "3-21g" # Use "minao" as the basis set
    mol.charge = 0 # Assign the charge of the molecule 
    mol.spin = 0 # Assign the spin of the molecule
    mol.build() # Build the molecule object



.. raw:: html

    <div id="3dmolviewer_15682344872293563"  style="position: relative; width: 400px; height: 400px">
            <p id="3dmolwarning_15682344872293563" style="background-color:#ffcccc;color:black">You appear to be running in JupyterLab (or JavaScript failed to load for some other reason).  You need to install the 3dmol extension: <br>
            <tt>jupyter labextension install jupyterlab_3dmol</tt></p>
            </div>
    <script>
    
    var loadScriptAsync = function(uri){
      return new Promise((resolve, reject) => {
        var tag = document.createElement('script');
        tag.src = uri;
        tag.async = true;
        tag.onload = () => {
          resolve();
        };
      var firstScriptTag = document.getElementsByTagName('script')[0];
      firstScriptTag.parentNode.insertBefore(tag, firstScriptTag);
    });
    };
    
    if(typeof $3Dmolpromise === 'undefined') {
    $3Dmolpromise = null;
      $3Dmolpromise = loadScriptAsync('https://3dmol.csb.pitt.edu/build/3Dmol.js');
    }
    
    var viewer_15682344872293563 = null;
    var warn = document.getElementById("3dmolwarning_15682344872293563");
    if(warn) {
        warn.parentNode.removeChild(warn);
    }
    $3Dmolpromise.then(function() {
    viewer_15682344872293563 = $3Dmol.createViewer($("#3dmolviewer_15682344872293563"),{backgroundColor:"white"});
    	viewer_15682344872293563.addModel("10\n\nHE         1.6180339887          0.0000000000          0.0000000000\nHE         1.3090169944          0.9510565163          0.0000000000\nH          0.5000000000          1.5388417686          0.0000000000\nH         -0.5000000000          1.5388417686          0.0000000000\nH         -1.3090169944          0.9510565163          0.0000000000\nH         -1.6180339887          0.0000000000          0.0000000000\nH         -1.3090169944         -0.9510565163          0.0000000000\nH         -0.5000000000         -1.5388417686          0.0000000000\nH          0.5000000000         -1.5388417686          0.0000000000\nH          1.3090169944         -0.9510565163          0.0000000000\n","xyz",{"keepH": true});
    	viewer_15682344872293563.setStyle({"sphere": {}});
    	viewer_15682344872293563.setStyle({"model": 0},{"sphere": {"colorscheme": "cyanCarbon", "scale": "0.2"}});
    	viewer_15682344872293563.zoomTo();
    viewer_15682344872293563.render();
    });
    </script>




.. parsed-literal::

    <pyscf.gto.mole.Mole at 0x7fcaf7d8d198>



As usual, we begin with importing the modules from OpenQEMIST.

.. code:: ipython3

    from openqemist import electronic_structure_solvers as ess
    from openqemist import problem_decomposition as pd

The DMET object can be used in two ways. The default behaviour of the
solver is to use an instance of an electronic structure solver to solve
all the fragments that are produced by the decomposition. This is shown
below:

.. code:: ipython3

    dmet = pd.DMETProblemDecomposition()
    dmet.electronic_structure_solver = ess.CCSDSolver()
    
    energy = dmet.simulate(mol, [2,2,2,2,2])
    
    print("DMET energy with CCSD is ", energy)


.. parsed-literal::

    DMET energy with CCSD is  -9.759969890236498


The DMET object can also use a different electronic structure solver to
solve each fragment. This is done with the optional ``fragment_solvers``
parameter. The value passed here should be a list of
``ElectronicStructureSolver`` instances that has as many elements as
there are fragments. This is shown below.

.. code:: ipython3

    # Create instances of the solvers that we want to use
    fci = ess.FCISolver()
    from openqemist import quantum_solvers as qs
    vqe = ess.VQESolver()
    vqe.hardware_backend_type = qs.MicrosoftQSharpParametricSolver
    vqe.ansatz_type = qs.MicrosoftQSharpParametricSolver.Ansatze.UCCSD
    
    # Use the VQE sovler to solve two helium fragments and the FCI solver for the hydrogen
    solvers = [vqe, vqe] + [fci for i in range(8)]
    energy = dmet.simulate(mol, [1,1,1,1,1,1,1,1,1,1], fragment_solvers=solvers)
    
    print("Mixed solver energy is ", energy)


.. parsed-literal::

    VQE : initial amplitudes
     [2e-05, 2e-05, 0.008336404499183074, 0.01731333728608542, -0.015990523785525974] 
    
    
    
    		Optimal UCCSD Singlet Energy: -3.1643666897843596
    		Optimal UCCSD Singlet Amplitudes: [-0.00430023 -0.00066218  0.00965449  0.02171837 -0.02179854]
    		Number of Function Evaluations :  88
    VQE : initial amplitudes
     [2e-05, 2e-05, 0.008336404498097656, 0.017313337286374863, -0.015990523784726086] 
    
    
    
    		Optimal UCCSD Singlet Energy: -3.1643666904553904
    		Optimal UCCSD Singlet Amplitudes: [-0.00427418 -0.00065326  0.00964003  0.02173165 -0.02180745]
    		Number of Function Evaluations :  85
    VQE : initial amplitudes
     [2e-05, 2e-05, 0.00833478236891564, 0.017313593028937303, -0.01598928337740948] 
    
    
    
    		Optimal UCCSD Singlet Energy: -3.1645570854987706
    		Optimal UCCSD Singlet Amplitudes: [-0.00433078 -0.00063696  0.00964952  0.0217214  -0.02181834]
    		Number of Function Evaluations :  79
    VQE : initial amplitudes
     [2e-05, 2e-05, 0.008334782367830496, 0.017313593029226616, -0.015989283376609684] 
    
    
    
    		Optimal UCCSD Singlet Energy: -3.164557084483473
    		Optimal UCCSD Singlet Amplitudes: [-0.00435145 -0.00065564  0.00964747  0.02172714 -0.02181109]
    		Number of Function Evaluations :  87
    VQE : initial amplitudes
     [2e-05, 2e-05, 0.008331805035645306, 0.01731406248067778, -0.015987006134357218] 
    
    
    
    		Optimal UCCSD Singlet Energy: -3.1649040972985785
    		Optimal UCCSD Singlet Amplitudes: [-0.00432664 -0.0006375   0.00964491  0.02173269 -0.02180754]
    		Number of Function Evaluations :  78
    VQE : initial amplitudes
     [2e-05, 2e-05, 0.008331805034560592, 0.017314062480967073, -0.01598700613355758] 
    
    
    
    		Optimal UCCSD Singlet Energy: -3.16490409774925
    		Optimal UCCSD Singlet Amplitudes: [-0.00429697 -0.00064323  0.00963988  0.0217343  -0.0218027 ]
    		Number of Function Evaluations :  84
    VQE : initial amplitudes
     [2e-05, 2e-05, 0.008332267250790905, 0.017313989596596194, -0.01598735971010735] 
    
    
    
    		Optimal UCCSD Singlet Energy: -3.1648492991218578
    		Optimal UCCSD Singlet Amplitudes: [-0.00429425 -0.00065077  0.00965251  0.02172731 -0.02180487]
    		Number of Function Evaluations :  91
    VQE : initial amplitudes
     [2e-05, 2e-05, 0.008332267249706118, 0.01731398959688549, -0.015987359709307668] 
    
    
    
    		Optimal UCCSD Singlet Energy: -3.164849297647699
    		Optimal UCCSD Singlet Amplitudes: [-0.00433949 -0.000644    0.00965045  0.02172857 -0.02180613]
    		Number of Function Evaluations :  76
    Mixed solver energy is  -9.76142175612128


