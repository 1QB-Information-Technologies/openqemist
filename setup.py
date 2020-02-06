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

import setuptools

with open("docs/pypi_readme.md", "r") as readme:
    long_description = readme.read() # heh

# This requires the a line of the module init to be something like:
# __version__ = "0.0.1"
with open("openqemist/__init__.py", "r") as init:
    line = init.readline()
    while line:
        if "__version__" in line:
            break
        line = init.readline()

    package_version = line.split("\"")[1]

setuptools.setup(
    name="openqemist",
    version=package_version,
    author="Rudi Plesch",
    author_email="rudi.plesch@1qbit.com",
    description="Quantum chemistry simulation library.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/1QB-Information-Technologies/openqemist",
    packages=setuptools.find_packages(),
    install_requires=['pyscf==1.6', 'numpy', 'scipy', 'qsharp', 'qiskit==0.13',
        'pyquil', 'forestopenfermion'],
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
        "License :: OSI Approved :: Apache Software License"
    ],
)
