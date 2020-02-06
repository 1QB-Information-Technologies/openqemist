# OpenQEMIST Docker image
This is a fully configured Docker image for running OpenQEMIST, and building the
documentation. OpenQEMIST can run in a container built from the included
Dockerfile. A docker-compose file is also provided to configure the Rigetti
quilc and qvm images.

## Getting started
On some platforms, limiting the number of threads that OpenMP uses to the number
of available processors improves performance. To do this, run the docker container
with `--env OMP_NUM_THREADS=2` to limit the number of threads. We have observed
this number to work well, but some experimentation should be done on each new
platform.

## Running Examples
If running Jupyter notebooks inside the docker container, run the docker
container with the port that Jupyter is serving on exposed. The default is 8888,
so use the `-p 8888:8888` argument.
Then run the Jypiter notebook with
```bash
$jupyter notebook --no-browser --allow-root --ip=0.0.0.0
```
## Rigetti configuration
Using the docker-compose file configures the environment to use the Rigetti
