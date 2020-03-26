#!/usr/bin/env python

from __future__ import print_function
import sys
import os
import re
import datetime as dt
import numpy as np
from netCDF4 import Dataset

####################################################################

# SemVer version
_major_version = 0
_minor_version = 3
_patch = 0
_release = "beta"

_version = '{0:d}.{1:d}.{2:d}-{3:s}'.format(_major_version, _minor_version,
                                            _patch, _release)
_date = '25-03-2020'


#
def nemo_rebuild(in_file=None,
                 out_file=None,
                 numdom=None,
                 variables=None,
                 nohalo=False,
                 fill=0,
                 verbose=False):
    """
    Rebuild NEMO/XIOS multiple output/restart files in a single file.

    Parameters
    ----------
    in_file  : (string) input file name template. May be one of:
                - Full name of a domain file, e.g. filebase_0000.nc
                - Base name of a domain file, e.g. filebase
    out_file : (string) rebuilt output file name. Default: filebase.nc
    numdom   : (integer) number of domains. Default: DOMAIN_number_total global attribute.
    variables: (string) rebuild only selected variables (comma separated list),
               e.g. variables='thetao,so,zos'
    nohalo   : (bool) Remove global domain halo (rebuilt restart file won't work) (default False)
    fill     : (numeric) Fill missing domains with this value if no _FillValue is present (default 0)
    verbose  : (bool) Verbose mode (default False)
    """
    #
    # This script should be invoked as: mpirun -n N python -m mpy4py nemo_rebuild.py ...
    # in order to avoid possible MPI deadlocks.
    # Load the mpi4py module if not provided on the command line (-m mpi4py)
    have_mpi4py = True
    if 'mpi4py' not in sys.modules.keys():
        print(
            'This script should be invoked as: mpirun -n N python -m mpy4py nemo_rebuild.py ...'
        )
        print(
            '-m mpi4py not provided on the command line. Loading the mpi4py module...'
        )
        try:
            import mpi4py
        except:
            print('Unable to load the mpi4py module. Executing sequentially.')
            have_mpi4py = False
    #
    # MPI related initialization
    rank = 0
    size = 1
    if (have_mpi4py):
        mpi4py = sys.modules['mpi4py']
        MPI = mpi4py.MPI
        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()
        size = comm.Get_size()
    #
    master = (rank == 0)
    parallel = (size > 1)
    #
    verbose = (verbose and master)
    #
    # Check input file name (path)
    if in_file is None:
        raise ValueError('Missing input file.')
    #
    # Check domains number
    if (numdom != None):
        if (numdom < 0):
            raise ValueError('NUMDOM must be nonnegative.')
        elif (numdom == 1):
            print('NUMDOM=1, nothing to be done. Exiting.')
            return
    #
    # TODO: Rebuild of a subset of variables needs more logic to identify
    # and handle auxiliary variables (coordinates, bounds, scalars, etc...)
    if (variables != None):
        raise RuntimeError(
            'Rebuild of a subset of variables not implemented yet!')
    #
    # Remove .nc file extension if any
    mtch = re.search('\.nc', in_file)
    if (mtch):
        in_file = in_file[:mtch.start()]
    #
    # Match & remove process number at the end of the file name if any (_[0-9]*)
    mtch = re.search('_\d+$', in_file)
    if (mtch):
        nume = in_file[mtch.start() + 1:mtch.endpos]
        in_file = in_file[:mtch.start()]
        pnd = len(nume)
        del nume
    else:
        # Try to figure out the number of digits in the process number if not given in in_file
        for pnd in range(6, -1, -1):
            test_file = in_file + '_' + '0' * pnd + '.nc'
            if os.path.isfile(test_file):
                break
        if (pnd == 0):
            raise RuntimeError('No input file found! Check in_file argument.')
        del test_file
    #
    del mtch
    #
    if (verbose):
        print('Input files name template: ',
              in_file + '_' + 'x' * pnd + '.nc',
              flush=True)
    #
    # Output file name (path)
    if out_file is None:
        out_file = in_file + '.nc'
    #
    if (verbose):
        print('Output file name: ', out_file, flush=True)
    #
    # Creation of the output file
    #
    excluded_gattrs = ('ibegin', 'jbegin', 'ni', 'nj', 'DOMAIN_number',
                       'DOMAIN_dimensions_ids', 'DOMAIN_size_local',
                       'DOMAIN_position_first', 'DOMAIN_position_last',
                       'DOMAIN_halo_size_start', 'DOMAIN_halo_size_end')
    #
    # open input first file
    in_file0 = in_file + '_' + '{0:0{width}d}'.format(rank, width=pnd) + '.nc'
    incid = Dataset(in_file0, 'r')
    #
    # Get number of domains from global attributes
    if (numdom != incid.DOMAIN_number_total):
        if (numdom != None):
            print(
                'WARNING: NUMDOM overwritten with DOMAIN_number_total global attribute!',
                flush=True)
        numdom = incid.DOMAIN_number_total
    #
    if (verbose):
        print(rank, 'Total domains number= ', numdom, flush=True)
    #
    # Get global dimensions
    gnx = incid.DOMAIN_size_global[0]
    gny = incid.DOMAIN_size_global[1]
    # Global dimensions with the global domain halo removed
    if (nohalo):
        gnx -= 2
        gny -= 1
    #
    # Create output file
    oncid = Dataset(out_file, 'w', format='NETCDF4_CLASSIC', parallel=parallel)
    #
    # Copy/update global attributes, except excluded ones
    gattrs = {
        k: v
        for k, v in incid.__dict__.items() if k not in excluded_gattrs
    }
    #
    if (nohalo):
        orig = gattrs['DOMAIN_size_global']
        gattrs['original_DOMAIN_size_global'] = orig
        gattrs['DOMAIN_size_global'] = np.array([gnx, gny], dtype=orig.dtype)
        orig = gattrs.pop('comment', None)
        if (orig == None):
            gattrs['comment'] = 'Global domain halo removed'
        else:
            gattrs['comment'] = orig + '\nGlobal domain halo removed'
    #
    orig = gattrs.pop('history', None)
    history = dt.datetime.strftime(dt.datetime.now(), '%c')
    if (have_mpi4py):
        history = history + ': mpirun -n {0:d} python -m mpy4py '.format(
            size) + ' '.join(sys.argv[:])
    else:
        history = history + ': python ' + ' '.join(sys.argv[:])
    #
    if (orig == None):
        gattrs['history'] = history
    else:
        gattrs['history'] = history + '\n' + orig
    #
    gattrs['nemo_rebuild_version'] = _version + ' (' + _date + ')'
    #
    oncid.setncatts(gattrs)
    #
    del gattrs
    del excluded_gattrs
    del orig
    del history
    #
    # Define dimensions & get global dims to be rebuilt
    gdimids = incid.DOMAIN_dimensions_ids
    gdims = ['', '']
    for name, dim in incid.dimensions.items():
        if (verbose):
            print(dim)
        if ((dim._dimid + 1) == gdimids[0]):
            oncid.createDimension(name, gnx)
            gdims[0] = name
        elif ((dim._dimid + 1) == gdimids[1]):
            oncid.createDimension(name, gny)
            gdims[1] = name
        else:
            oncid.createDimension(name,
                                  len(dim) if not dim.isunlimited() else None)
    #
    if (verbose):
        print('Global dimensions to be rebuilt: ', gdims)
    del dim
    del gdimids
    #
    # Define variables and their attributes
    # if _FillValue is not defined, define a temporary fill value (default 0)
    # to fill the missing land subdomains,
    # to be deleted at the end of the rebuilt process
    fv_del = {}
    for name, var in incid.variables.items():
        if (verbose):
            print('\nVar Def:\n', rank, name, var, flush=True)
        vattrs = var.__dict__
        fv = vattrs.pop('_FillValue', None)
        mv = vattrs.pop('missing_value', None)
        if (fv == None):
            if (mv == None):
                fv = fill
                fv_del[name] = True
            else:
                fv = mv
        elif (mv != None):
            if (fv != mv):
                mv = fv
            vattrs['missing_value'] = mv
        #
        ovid = oncid.createVariable(var.name,
                                    var.dtype,
                                    var.dimensions,
                                    fill_value=fv)
        ovid.setncatts(vattrs)
    #
    del vattrs
    del var
    #
    del in_file0
    #
    # Close files
    incid.close()
    #
    if (numdom == None):
        oncid.close()
        raise RuntimeError('NUMDOM=None.')
    #
    # TODO: handle the case where (numdom%size)!=0
    if (numdom % size != 0):
        oncid.close()
        raise RuntimeError(
            'The number of MPI tasks must be a divisor of the number of domains! '
        )
    #
    # Main parallel loop
    for niter in range(numdom // size):
        #
        # Open input files
        lin_file = in_file + '_' + '{0:0{width}d}'.format(niter * size + rank,
                                                          width=pnd) + '.nc'
        if (verbose):
            print('\n',
                  rank,
                  niter,
                  'Input file name: ',
                  lin_file,
                  '\n',
                  flush=True)
        incid = Dataset(lin_file, 'r')
        #
        # NetCDF IDs of the dimensions to be rebuilt
        gdimids = incid.DOMAIN_dimensions_ids
        #
        # Local domain size & index
        lnx = incid.DOMAIN_size_local[0]
        lny = incid.DOMAIN_size_local[1]
        li1 = 1
        lj1 = 1
        li2 = lnx
        lj2 = lny
        #
        # Local domain halo
        hi1 = incid.DOMAIN_halo_size_start[0]
        hj1 = incid.DOMAIN_halo_size_start[1]
        hi2 = incid.DOMAIN_halo_size_end[0]
        hj2 = incid.DOMAIN_halo_size_end[1]
        #
        # Take halo into account
        li1 += hi1
        lj1 += hj1
        li2 -= hi2
        lj2 -= hj2
        #
        # Global domain position index
        gi1 = incid.DOMAIN_position_first[0] + li1 - 1
        gj1 = incid.DOMAIN_position_first[1] + lj1 - 1
        gi2 = incid.DOMAIN_position_first[0] + li2 - 1
        gj2 = incid.DOMAIN_position_first[1] + lj2 - 1
        #
        # Convert position index from Fortran to Python (0 based)
        li1 -= 1
        lj1 -= 1
        gi1 -= 1
        gj1 -= 1
        #
        #bnd=False
        #if (gi1==0 or gi2>=incid.DOMAIN_size_global[0] or gj2>=incid.DOMAIN_size_global[1]):
        #    bnd=True
        #    print(rank, niter, 'PRE', gi1,gi2,gj1,gj2,li1,li2,lj1,lj2, lnx, lny, flush=True)
        #
        # Remove global domain halo if requested
        if (nohalo):
            if (gi1 == 0):  # First column
                li1 += 1
                gi2 -= 1
            else:  #
                gi1 -= 1
                gi2 -= 1
            if ((gi2 + 1) == incid.DOMAIN_size_global[0]):  # Last column
                li2 -= 1
                gi2 -= 1
            if (gj2 == incid.DOMAIN_size_global[1]):  # Last row
                lj2 -= 1
                gj2 -= 1
        #
        if ((gj2 - gj1) != (lj2 - lj1) or (gi2 - gi1) != (li2 - li1)):
            print(rank, niter, '\nError with index:', flush=True)
            print(rank,
                  niter,
                  gi1,
                  gi2,
                  gj1,
                  gj2,
                  li1,
                  li2,
                  lj1,
                  lj2,
                  flush=True)
            oncid.close()
            incid.close()
            raise RuntimeError('Error with index!')
        #
        #if (bnd):
        #    print(rank, niter, 'POST', gi1,gi2,gj1,gj2,li1,li2,lj1,lj2, lnx, lny, flush=True)
        #
        # Loop over input variables
        for name, var in incid.variables.items():
            #
            if (verbose):
                print('\n', rank, niter, name, var, flush=True)
            #
            # Set collective I/O mode
            ovid = oncid.variables[name]
            if (parallel):
                ovid.set_collective(True)
            #
            # Variables to be rebuilt
            if (gdims[0] in var.dimensions and gdims[1] in var.dimensions):
                #
                if (var.ndim == 2):
                    ovid[gj1:gj2, gi1:gi2] = var[lj1:lj2, li1:li2]
                elif (var.ndim == 3):
                    if (var._nunlimdim > 0):  # Record variable
                        ovid[:, gj1:gj2, gi1:gi2] = var[:, lj1:lj2, li1:li2]
                    else:  # Coordinate bounds
                        ovid[gj1:gj2, gi1:gi2, :] = var[lj1:lj2, li1:li2, :]
                elif (var.ndim == 4):
                    ovid[:, :, gj1:gj2, gi1:gi2] = var[:, :, lj1:lj2, li1:li2]
                #
                # Variables NOT to be rebuilt
            elif (niter == 0):
                ovid = oncid.variables[name]
                ovid[:] = var[:]
        #
        incid.close()
    #
    del ovid
    del incid
    #
    # Remove temporary _FillValue
    if (len(fv_del) > 0):
        for name, var in oncid.variables.items():
            if (fv_del.pop(name, False)):
                var.delncattr('_FillValue')
    #
    oncid.close()
    #
    del fv_del
    del oncid


if __name__ == "__main__":
    import argparse
    #
    parser = argparse.ArgumentParser(
        prog='mpirun -n N python -m mpi4py ' + sys.argv[0],
        description='NEMO output/restart file rebuilder',
        epilog=
        'Rebuild NEMO/XIOS multiple output/restart files in a single file')
    parser.add_argument("-i",
                        "--input",
                        action="store",
                        dest="in_file",
                        default=None,
                        required=True,
                        help="input file name")
    parser.add_argument("-o",
                        "--output",
                        action="store",
                        dest="out_file",
                        default=None,
                        required=False,
                        help="output file name")
    parser.add_argument("-n",
                        "--numdom",
                        action="store",
                        dest="numdom",
                        default=None,
                        type=int,
                        required=False,
                        help="Number of domains (input files)")
    parser.add_argument(
        "-f",
        "--fill",
        action="store",
        dest="fill",
        default=0,
        required=False,
        help=
        "Fill missing domains with %(dest)s if no _FillValue is present (default 0)"
    )
    parser.add_argument("-v",
                        "--variable",
                        action="store",
                        dest="variables",
                        default=None,
                        required=False,
                        help="Variable(s) to rebuild")
    parser.add_argument(
        "-r",
        "--remove_halo",
        action="store_true",
        dest="nohalo",
        default=False,
        required=False,
        help="Remove global domain halo (rebuilt restart file won't work)")
    parser.add_argument("--verbose",
                        action="store_true",
                        dest="verbose",
                        default=False,
                        required=False,
                        help="Verbose mode")
    #
    args = parser.parse_args()
    #
    nemo_rebuild(in_file=args.in_file,
                 out_file=args.out_file,
                 numdom=args.numdom,
                 variables=args.variables,
                 nohalo=args.nohalo,
                 fill=args.fill,
                 verbose=args.verbose)


