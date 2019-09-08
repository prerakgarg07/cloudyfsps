#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import (division, print_function, absolute_import,
                        unicode_literals)
import os
import sys
import numpy as np
import fsps
from past.utils import old_div
from cloudyfsps.ASCIItools import (writeASCII, compileASCII, checkCompiled, compiledExists)
from cloudyfsps.cloudyInputTools import *
from cloudyfsps.generalTools import calcForLogQ
from cloudyfsps.cloudyOutputTools import *
from cloudyfsps.outputFormatting import *
#import ipdb
#runMake, formatAllOutput, writeFormattedOutput)

# this code snippet goes through every step needed
# to integrate FSPS into Cloudy.
# This example uses stellar pops with a constant SFH
# as the input ionizing source.
# 1. Write an ascii file in Cloudy format with grid
#    of FSPS spectra in all available ages and
#    metallicities
# 2. Compile asii file into binary format required
#    for Cloudy use. Assumes $CLOUDY_EXE is set to
#    your /path/to/cloudy.exe
# 3. Writes Cloudy input files for a subset of grid
#    parameters.
# 4. Runs Cloudy on the *.in files
# 5. Formats the various output files

zsun = 0.0142
exec_write_ascii = False
exec_write_input = False
exec_run_cloudy = False
exec_write_output = True
exec_gen_FSPS_grid = True

# Function to write the ascii file.
# This is where you set the properties of the
# ionizing spectrum (SSP/CSFH, IMF, FBHB, etc)

def hden_ascii(fileout, **kwargs):
    # change these parameters to modify the ionizing source grid
    # default mode is to produce an ascii grid in age and Z,
    # though different variables and more dimensions are possible.
    sp_dict = dict(zcontinuous=1,
                   imf_type=2,
                   sfh=0,
                   const=0.0,
                   sf_start=0.0)
    sp = fsps.StellarPopulation(**sp_dict)
    # all ages and Zs
    ages = 10.**sp.log_age
    logZs = np.log10(old_div(sp.zlegend,zsun))
    modpars = [(age, logZ) for age in ages for logZ in logZs]
    lam = sp.wavelengths
    all_fluxs = []
    for logZ in logZs:
        sp.params['logzsol'] = logZ
        all_fluxs.append(sp.get_spectrum()[1]) #lsun per hz
    nmod = len(modpars)
    # flatten flux for writing
    flat_flux = np.array([all_fluxs[j][i]
                          for i in range(len(ages))
                          for j in range(len(logZs))])
    # this function is flexible, ndim can be 3/4/n.
    # in this example, however, ndim is 2 (age, logz).
    writeASCII(fileout, lam, flat_flux, modpars,
               nx=len(lam), ndim=2, npar=2, nmod=nmod)
    return
#---------------------------------------------------------------------
# ASCII FILE: WRITE AND COMPILE
#---------------------------------------------------------------------
# assumes you have $CLOUDY_EXE and $CLOUDY_DATA_PATH set as sys vars.

# name of ascii file
ascii_file = "FSPS_MIST_C3K.ascii"

# or if there is an already-compiled one you want to use, specify here
compiled_ascii = "{}.mod".format(ascii_file.split(".")[0])

if exec_write_ascii:
    print("Executing write ascii sequence...")
    if not compiledExists(ascii_file):
        print("No compiled model exists...Writing.")
        hden_ascii(ascii_file)
        print("Compiling {} with Cloudy".format(ascii_file))
        compileASCII(ascii_file)
        print("Checking to see if compilation was successful...")
        if checkCompiled(ascii_file):
            print("Your model {} is ready to run.".format(compiled_ascii))
        else:
            sys.exit()
    else:
        print("{} already exists.".format(compiled_ascii))

#---------------------------------------------------------------------
# WRITE CLOUDY INPUT
#---------------------------------------------------------------------
# local folder to read and write *.in, *.out files
mod_dir = '/ufrc/narayanan/prerakgarg/cloudy_runs/c17_c3k_mist/'
#mod_dir = "/home/prerakgarg/redo/c17/"
mod_prefix = 'ZAU'

# GRID PARAMETERS FOR CLOUDY RUN
#--------------
# ages between 1 and 7 Myr
#ages = np.linspace(1., 7., 7)*1.e6i


ages = np.array([0.5, 1, 2, 3, 4, 5, 6, 7, 10, 20])*1.e6
# stellar metallicities
logZs = np.array([-1.98, -1.5, -1.0, -0.6, -0.4, -0.3, -0.2, -0.1, 0.0, 0.1, 0.198])
# ionization parameters between -4 and -1
logUs = np.array([-4.0, -3.5, -3.0, -2.5, -2.0, -1.5, -1.0])

'''
ages = np.array([1, 2,])*1.e6
# stellar metallicities
logZs = np.array([-1.5])
# ionization parameters between -4 and -1
logUs = np.array([-4.0])
#logUs = [-1.]
# Hydrogen density between 30 and 400
#nhs = np.arange(50., 450., 200) # density of gas, cm-3
'''
nhs = [100]

# Other default parameters based off of Byler+2017
Rinners =  np.array([19.]) # inner radius of HII region, 3pc
efrac = -1.0 # calculation is stopped when H is 10^efrac % neutral
set_name='dopita' # abundances from Dopita+2001
dust=False # don't include dust in nebula
extra_output=False # include lots of outputs
#-----------------------------------------------------------------

# iterate through all of the above parameters
# calcForLogQ just calculates Q = U*4*pi*Ri^2*nH

pars = np.array([(Z, a, U, R, calcForLogQ(logU=U, Rinner=10.0**R, nh=n), n, efrac)
                 for Z in logZs
                 for a in ages
                 for U in logUs
                 for R in Rinners
                 for n in nhs])

if exec_write_input:
    print('Writing input files...')
    writeParamFiles(dir_=mod_dir,
                    model_prefix=mod_prefix,
                    cloudy_mod=compiled_ascii,
                    run_cloudy=False, # don't run yet
                    ages=ages,
                    logZs=logZs,
                    logUs=logUs,
                    r_inners=Rinners,
                    nhs=nhs,
                    use_Q=True,
                    # if False, will use logU;
                    # does not matter in this case,
                    # since Q is calculated at
                    # each specified logU.
                    verbose=False, # don't print output to screen
                    set_name=set_name,
                    dust=dust,
                    extra_output=extra_output)
    print('Wrote {} param files'.format(len(pars)))
else:
    print('Skipping input writing.')


#---------------------------------------------------------------------
# RUN CLOUDY ON ALL INPUT FILES
#---------------------------------------------------------------------
if exec_run_cloudy:
    print("Running Cloudy....")
    runMake(dir_=mod_dir, n_proc=4, model_name=mod_prefix)
    print("Cloudy finished.")
else:
    print("Not running Cloudy. Skipping to formatting output.")


#---------------------------------------------------------------------
# FORMAT OUTPUT
#---------------------------------------------------------------------
if exec_write_output:
    print("Formatting output files...\n")
    formatAllOutput(mod_dir, mod_prefix, write_line_lum=False)
else:
    print("\n\nNot formatting output. DONE.")

if exec_gen_FSPS_grid:
    print("Creating FSPS input grids...")
    writeFormattedOutput(mod_dir, mod_prefix, "")
else:
    print("\n\nNot formatting FSPS output. DONE.")
