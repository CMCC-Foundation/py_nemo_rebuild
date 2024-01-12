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

### Setup Anaconda environment

The working environment can be setup using the provided anaconda3 `environment.yml` file:
```
conda env create -f environment.yml

conda activate nemo_rebuild
```
> [!IMPORTANT]
> In order to avoid issues during later updates, it is necessary to pin the netCDF4 and MPICH python modules in the anaconda environment :
> 
> ```
> echo 'netcdf4 * *mpich*' >> ${CONDA_PREFIX}/conda-meta/pinned
> echo 'mpich * *external*' >> ${CONDA_PREFIX}/conda-meta/pinned
> ```

> [!IMPORTANT]
> The choice of the external build of the MPICH package enable the use of the MPI library installed on the HPC platform (Intel MPI, etc...), provided that the environment is properly set up (e.g. module load impi-2021.6.0/2021.6.0).

### Usage of nemo\_rebuild.py

First of all the correct environment needs to be set up with:
```
module load impi-2021.6.0/2021.6.0     # or whatever MPI library is available on the HPC platoform
conda activate nemo_rebuild
```
Once the environment is properly set up the script nemo_rebuild.py can be run:

```
python nemo_rebuild.py -h
usage: mpirun -n N python nemo_rebuild.py [-h] -i IN_FILE [-o OUT_FILE] [-n NUMDOM] [-f FILL] [-v VARIABLES] [-r] [-c]
                                          [--verbose] [-V]

NEMO output/restart file rebuilder

options:
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
  -c, --classic         use NetCDF4 Classic format (default: False)
  --verbose             Verbose mode
  -V, --version         show program's version number and exit

Rebuild NEMO/XIOS multiple output/restart files in a single file
```

**IN\_FILE** can be one of:
 * Full name of a subdomain file, e.g. filebase\_0000.nc
 * Base name of a subdomain file, e.g. filebase

Examples
--------

The rebuild script can be run either sequentially or using MPI parallelization as described below.

### Sequential run

We want to rebuild a set of NEMO restart files from ORCA1 configuration. The minimal command line to perform a simple **sequential** reconstruction is
```
python nemo_rebuild.py -i ORCA1_01305240_restart
```
and the script will automatically detect the number of subdomain files.
In this case it's not mandatory to set up the MPI environment (module load ...).

### Parallel run

Let's rebuild NEMO output files from ORCA025 configuration by using the **MPI** interface to speed up things ... 
```
mpirun -n N python nemo_rebuild.py -i ORCA025_1m_20000101_20000131_grid_T
```
where N is the number of MPI tasks required.
The number N of MPI tasks cannot exceed the number of subdomains and does not necessarily have to be an integer divisor of the number of subdomains.

Install as a package
--------------------

When installed as a package, an executable named nemo\_rebuild\_py will be added to your path, so one can call it from anywhere, provided that the conda and MPI environments are properly set up. The package can be installed by launching:
```
pip install .
```
from the root folder of the repository.

Authors
-------

* [Pier Giuseppe Fogli](https://github.com/pgf) 
* [Tomas Lovato](https://github.com/tomaslovato)
* [Momme Butensch&#246;n](https://github.com/mommebutenschoen)

