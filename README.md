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

In order to avoid issues during later updates, It is necessary to pin the netCDF4 python module in the anaconda environment :

```
echo 'netcdf4 *mpich*' >> ${CONDA_PREFIX}/conda-meta/pinned
```

### Usage of nemo\_rebuild.py

```
python nemo_rebuild.py -h
usage: mpirun -n N python nemo_rebuild.py [-h] -i IN_FILE [-o OUT_FILE]
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

**IN\_FILE** can be one of:
 * Full name of a subdomain file, e.g. filebase\_0000.nc
 * Base name of a subdomain file, e.g. filebase

Examples
--------

#### The rebuild script can be run either sequentially or using MPI parallelization as described below.

We want to rebuild a set of nemo restart files from ORCA1 configuration. The minimal command line to perform a simple **sequential** reconstruction is
```
python nemo_rebuild.py -i ORCA1_01305240_restart
```
and the script will automatically detect the number of subdomain files.

Let's rebuild nemo output files from ORCA025 configuration by using the **MPI**  interface to speed up things ... 
```
mpirun -n N python nemo_rebuild.py -i ORCA025_1m_20000101_20000131_grid_T
```
where N is the number of MPI tasks required.

Install package
---------------

When installed as package, an executable names nemo_rebuild_py will be added to your path, so you can call it from wherever which executees the main python function of the script. The package can be installed by launching
```
pip install .
```
from the root folder of the repository.

Authors
-------

* Pier Giuseppe Fogli 
* Tomas Lovato
* Momme Butensch#&X00F6;n

