#!/usr/bin/env python

import pynbody
import numpy as np
import h5py
import os
from argparse import ArgumentParser

# Define arguments ###
parser = ArgumentParser()
parser.add_argument("-i", "--infile", action="store", dest="inputFile", type=str, default=None,
                    help="The path to the simulation file to be read")
parser.add_argument("-o", "--outfile", action="store", dest="outputFile", type=str, default=None,
                    help="The path to the SimSpin file to be produced")

# Parse arguments ###
args = parser.parse_args()

# Validate inputs ###

if args.inputFile is None:
    print("\nERROR: No --infile argument given! \n For help, use \"python create_SimSpinFile.py -h\".\n --- ABORTING ---\n")
    exit()  # Check that an input file has been specified

if args.outputFile is None:
    print("\nERROR: No --outfile argument given! \n For help, use \"python create_SimSpinFile.py -h\".\n --- ABORTING ---\n")
    exit()  # Check that an output file has been specified

if os.path.exists(args.inputFile):
    inputFile = args.inputFile
else:
    print("\nERROR: --infile does not exist!\n --- ABORTING ---\n")
    exit()  # Check that inputFile path leads to a file


print("Input file: ", args.inputFile)

ds = pynbody.load(args.inputFile)  # open supplied simulation file

try:
    ptype = ds.header.npart
    print("This input file is a Gadget binary file with ", ptype[0], "gas particles, ", ptype[1], "DM particles, ",
          ptype[2], "disc particles and ", ptype[3], "bulge particles.")

except AttributeError:
    ptype = ds.properties["NumPart_ThisFile"]
    print("This input file is a Gadget HDF5 file with ", ptype[0], "gas particles, ", ptype[1], "DM particles, ",
          ptype[2], "disc particles and ", ptype[3], "bulge particles.")
    # determining if the inputFile is binary or HDF5

hf = h5py.File(args.outputFile, "w")    # creating HDF5 file for SimSpin
hf.create_dataset("x", data = ds["x"])
hf.create_dataset("y", data = ds["y"])
hf.create_dataset("z", data = ds["z"])
hf.create_dataset("vx", data = ds["vx"])
hf.create_dataset("vy", data = ds["vy"])
hf.create_dataset("vz", data = ds["vz"])
hf.create_dataset("Mass", data = ds["mass"])
hf.create_dataset("Part_Type", data = np.array([0] * ptype[0] + [1] * ptype[1] + [2] * ptype[2] + [3] * ptype[3]))
hf.close()
print("New SimSpin file written at: ", args.outputFile)