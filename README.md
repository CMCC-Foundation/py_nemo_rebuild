Pythonic NEMO Rebuild Tool
==========================

Rebuild NEMO/XIOS multiple output/restart files in a single file.

Getting Started
---------------

This is a python script which rebuilds NEMO/XIOS output/restart files
split over multiple subdomains in a single file.
While the script can be run sequentially (i.e. single task), it is designed
to run in parallel (SPMD) leveraging the MPI library, in order to speed up
the reconstruction.

### Prerequisites

* MPI library & run-time environment (mpirun, etc...) [optional]
* NetCDF4/HDF5 libraries [parallel version optional]
* Python
* numpy
* netCDF4
* mpi4py [optional]

### Setup an anaconda environment

The rebuild script can be run sequentially:

```
python nemo_rebuild.py <options> <files>
```

In order to take advantage of the parallelization, an MPI enabled environment
is required:

* mpi4py python module
* parallel NetCDF4/HDF5 libraries
* netCDF4 python module
* MPI run time environment (mpirun, etc...)

The simplest way to set up such an environment is through Anaconda:

```
conda create -n mpipy -c conda-forge python=3 numpy netcdf4=*=mpi_mpich* mpi4py

conda activate mpipy
```

BEWARE: when updating the environment conda keeps trying to switch to the OpenMPI-based packages.
Please avoid to update the environment for now, a fix is under research.

### Run the script in parallel

This script should be invoked as:

```
mpirun -n N python -m mpy4py nemo_rebuild.py ...
```
where N is the number of MPI tasks required.

### Command line options

```
python nemo_rebuild.py -h
usage: mpirun -n N python -m mpi4py nemo_rebuild.py [-h] -i IN_FILE [-o OUT_FILE]
                                                    [-n NUMDOM] [-f FILL] [-v VARIABLES]
                                                    [-r] [--verbose]

NEMO output/restart file rebuilder

optional arguments:
  -h, --help            show this help message and exit
  -i IN_FILE, --input IN_FILE
                        input file name
  -o OUT_FILE, --output OUT_FILE
                        output file name
  -n NUMDOM, --numdom NUMDOM
                        Number of domains (input files)
  -f FILL, --fill FILL  Fill missing domains with fill if no _FillValue is present (default 0)
  -v VARIABLES, --variable VARIABLES
                        Variable(s) to rebuild
  -r, --remove_halo     Remove global domain halo (rebuilt restart file won't work)
  --verbose             Verbose mode

Rebuild NEMO/XIOS multiple output/restart files in a single file
```

IN_FILE can be one of:
 * Full name of a subdomain file, e.g. filebase_0000.nc
 * Base name of a subdomain file, e.g. filebase

Authors
-------

* Pier Giuseppe Fogli ()


