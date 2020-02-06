# Minimal makefile for Sphinx documentation
#

# You can set these variables from the command line.
SPHINXOPTS    =
SPHINXBUILD   = sphinx-build
SOURCEDIR     = .
BUILDDIR      = _build

# Put it first so that "make" without argument is like "make help".
help:
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

.PHONY: help Makefile

# Catch-all target: route all unknown targets to Sphinx using the new
# "make mode" option.  $(O) is meant as a shortcut for $(SPHINXOPTS).
%: Makefile
	@sphinx-apidoc  -o refman/ .. ../setup.py ../openqemist/tests/ ../benchmarks
	@$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

multisolver_notebook:
	@jupyter nbconvert --ExecutePreprocessor.timeout=6000 --execute --to rst ../examples/dmet_multisolver.ipynb
	@mv ../examples/dmet_multisolver.rst .

dmet_notebook:
	@jupyter nbconvert --ExecutePreprocessor.timeout=6000 --execute --to rst ../examples/dmet_microsoft_qsharp.ipynb
	@mv ../examples/dmet_microsoft_qsharp.rst .

vqe_notebook:
	@jupyter nbconvert --ExecutePreprocessor.timeout=6000 --execute --to rst ../examples/vqe_microsoft_qsharp.ipynb
	@mv ../examples/vqe_microsoft_qsharp.rst .

rigetti_notebook:
	@jupyter nbconvert --ExecutePreprocessor.timeout=6000 --execute --to rst ../examples/rigetti_example.ipynb
	@mv ../examples/rigetti_example.rst .

notebooks: multisolver_notebook dmet_notebook vqe_notebook qiskit_notebook

