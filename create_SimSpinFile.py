#!/usr/bin/env python

import pynbody
import numpy as np
import h5py
import os
from argparse import ArgumentParser

# Define arguments ###
parser = ArgumentParser()
parser.add_argument("-i", "--infile", action="store", dest="inputFile", type=str, default=None,
                    help="The path to the simulation file you wish to process with SimSpin")
parser.add_argument("-o", "--outfile", action="store", dest="outputFile", type=str, default=None,
                    help="The path to the SimSpin HDF5 file produced")

# Parse arguments ###
args = parser.parse_args()

# Validate inputs ###

if args.inputFile is None:
    print("\nERROR: No --infile argument given! \n "
          "For help, use \"python create_SimSpinFile.py -h\".\n"
          " --- ABORTING ---\n")
    exit()  # Check that an input file has been specified

if args.outputFile is None:
    print("\nERROR: No --outfile argument given! \n "
          "For help, use \"python create_SimSpinFile.py -h\".\n"
          " --- ABORTING ---\n")
    exit()  # Check that an output file has been specified

if os.path.exists(args.inputFile):
    inputFile = args.inputFile
else:
    print("\nERROR: --infile does not exist!\n "
          "--- ABORTING ---\n")
    exit()  # Check that inputFile path leads to a file


print("Input file: ", args.inputFile)

ds = pynbody.load(args.inputFile)  # open supplied simulation file
ds.physical_units()

try:
    ptype = ds.header.npart
    print("This input file is a Gadget binary file with ", ptype[0], "gas particles, ", ptype[1],
          "DM particles, ", ptype[2], "disc particles", ptype[3], "bulge particles,", ptype[4],
          "stellar particles and", ptype[5], "boundary particles.")

except AttributeError:
    ptype = ds.properties["NumPart_ThisFile"]
    print("This input file is a Gadget HDF5 file with ", ptype[0], "gas particles, ", ptype[1],
          "DM particles, ", ptype[2], "disc particles and ", ptype[3], "bulge particles,", ptype[4],
          "stellar particles and", ptype[5], "boundary particles.")
    # determining if the inputFile is binary or HDF5

hf = h5py.File(args.outputFile, "w")    # creating HDF5 file for SimSpin
pcount = np.cumsum(ptype) # calculating markers for the positions that each ptype starts and ends

if ptype[0] != 0:  # if there are ptype[0] = gas particles in the file
    part0 = hf.create_group("PartType0")
    x0 = part0.create_dataset("x", data=ds["x"][0:pcount[0]])
    x0.attrs["Units"] = str(ds["x"].units)
    x0.attrs["PartNum"] = ptype[0]
    y0 = part0.create_dataset("y", data=ds["y"][0:pcount[0]])
    y0.attrs["Units"] = str(ds["y"].units)
    y0.attrs["PartNum"] = ptype[0]
    z0 = part0.create_dataset("z", data=ds["z"][0:pcount[0]])
    z0.attrs["Units"] = str(ds["z"].units)
    z0.attrs["PartNum"] = ptype[0]
    vx0 = part0.create_dataset("vx", data=ds["vx"][0:pcount[0]])
    vx0.attrs["Units"] = str(ds["vx"].units)
    vx0.attrs["PartNum"] = ptype[0]
    vy0 = part0.create_dataset("vy", data=ds["vy"][0:pcount[0]])
    vy0.attrs["Units"] = str(ds["vy"].units)
    vy0.attrs["PartNum"] = ptype[0]
    vz0 = part0.create_dataset("vz", data=ds["vz"][0:pcount[0]])
    vz0.attrs["Units"] = str(ds["vz"].units)
    vz0.attrs["PartNum"] = ptype[0]
    mass0 = part0.create_dataset("Mass", data=ds["mass"].in_units('1e10 Msol')[0:pcount[0]])
    mass0.attrs["Units"] = str(ds["mass"].in_units('1e10 Msol').units)
    mass0.attrs["PartNum"] = ptype[0]
    
if ptype[1] != 0:  # if there are ptype[1] = DM particles in the file
    part1 = hf.create_group("PartType1")
    x1 = part1.create_dataset("x", data=ds["x"][pcount[0]:pcount[1]])
    x1.attrs["Units"] = str(ds["x"].units)
    x1.attrs["PartNum"] = ptype[1]
    y1 = part1.create_dataset("y", data=ds["y"][pcount[0]:pcount[1]])
    y1.attrs["Units"] = str(ds["y"].units)
    y1.attrs["PartNum"] = ptype[1]
    z1 = part1.create_dataset("z", data=ds["z"][pcount[0]:pcount[1]])
    z1.attrs["Units"] = str(ds["z"].units)
    z1.attrs["PartNum"] = ptype[1]
    vx1 = part1.create_dataset("vx", data=ds["vx"][pcount[0]:pcount[1]])
    vx1.attrs["Units"] = str(ds["vx"].units)
    vx1.attrs["PartNum"] = ptype[1]
    vy1 = part1.create_dataset("vy", data=ds["vy"][pcount[0]:pcount[1]])
    vy1.attrs["Units"] = str(ds["vy"].units)
    vy1.attrs["PartNum"] = ptype[1]
    vz1 = part1.create_dataset("vz", data=ds["vz"][pcount[0]:pcount[1]])
    vz1.attrs["Units"] = str(ds["vz"].units)
    vz1.attrs["PartNum"] = ptype[1]
    mass1 = part1.create_dataset("Mass", data=ds["mass"].in_units('1e10 Msol')[pcount[0]:pcount[1]])
    mass1.attrs["Units"] = str(ds["mass"].in_units('1e10 Msol').units)
    mass1.attrs["PartNum"] = ptype[1]
    
if ptype[2] != 0:  # if there are ptype[2] = Disc particles in the file
    part2 = hf.create_group("PartType2")
    x2 = part2.create_dataset("x", data=ds["x"][pcount[1]:pcount[2]])
    x2.attrs["Units"] = str(ds["x"].units)
    x2.attrs["PartNum"] = ptype[2]
    y2 = part2.create_dataset("y", data=ds["y"][pcount[1]:pcount[2]])
    y2.attrs["Units"] = str(ds["y"].units)
    y2.attrs["PartNum"] = ptype[2]
    z2 = part2.create_dataset("z", data=ds["z"][pcount[1]:pcount[2]])
    z2.attrs["Units"] = str(ds["z"].units)
    z2.attrs["PartNum"] = ptype[2]
    vx2 = part2.create_dataset("vx", data=ds["vx"][pcount[1]:pcount[2]])
    vx2.attrs["Units"] = str(ds["vx"].units)
    vx2.attrs["PartNum"] = ptype[2]
    vy2 = part2.create_dataset("vy", data=ds["vy"][pcount[1]:pcount[2]])
    vy2.attrs["Units"] = str(ds["vy"].units)
    vy2.attrs["PartNum"] = ptype[2]
    vz2 = part2.create_dataset("vz", data=ds["vz"][pcount[1]:pcount[2]])
    vz2.attrs["Units"] = str(ds["vz"].units)
    vz2.attrs["PartNum"] = ptype[2]
    mass2 = part2.create_dataset("Mass", data=ds["mass"].in_units('1e10 Msol')[pcount[1]:pcount[2]])
    mass2.attrs["Units"] = str(ds["mass"].in_units('1e10 Msol').units)
    mass2.attrs["PartNum"] = ptype[2]
    
if ptype[3] != 0:  # if there are ptype[3] = Bulge particles in the file
    part3 = hf.create_group("PartType3")
    x3 = part3.create_dataset("x", data=ds["x"][pcount[2]:pcount[3]])
    x3.attrs["Units"] = str(ds["x"].units)
    x3.attrs["PartNum"] = ptype[3]
    y3 = part3.create_dataset("y", data=ds["y"][pcount[2]:pcount[3]])
    y3.attrs["Units"] = str(ds["y"].units)
    y3.attrs["PartNum"] = ptype[3]
    z3 = part3.create_dataset("z", data=ds["z"][pcount[2]:pcount[3]])
    z3.attrs["Units"] = str(ds["z"].units)
    z3.attrs["PartNum"] = ptype[3]
    vx3 = part3.create_dataset("vx", data=ds["vx"][pcount[2]:pcount[3]])
    vx3.attrs["Units"] = str(ds["vx"].units)
    vx3.attrs["PartNum"] = ptype[3]
    vy3 = part3.create_dataset("vy", data=ds["vy"][pcount[2]:pcount[3]])
    vy3.attrs["Units"] = str(ds["vy"].units)
    vy3.attrs["PartNum"] = ptype[3]
    vz3 = part3.create_dataset("vz", data=ds["vz"][pcount[2]:pcount[3]])
    vz3.attrs["Units"] = str(ds["vz"].units)
    vz3.attrs["PartNum"] = ptype[3]
    mass3 = part3.create_dataset("Mass", data=ds["mass"].in_units('1e10 Msol')[pcount[2]:pcount[3]])
    mass3.attrs["Units"] = str(ds["mass"].in_units('1e10 Msol').units)
    mass3.attrs["PartNum"] = ptype[3]
    
if ptype[4] != 0:  # if there are ptype[4] = Stellar particles in the file
    part4 = hf.create_group("PartType4")
    x4 = part4.create_dataset("x", data=ds["x"][pcount[3]:pcount[4]])
    x4.attrs["Units"] = str(ds["x"].units)
    x4.attrs["PartNum"] = ptype[4]
    y4 = part4.create_dataset("y", data=ds["y"][pcount[3]:pcount[4]])
    y4.attrs["Units"] = str(ds["y"].units)
    y4.attrs["PartNum"] = ptype[4]
    z4 = part4.create_dataset("z", data=ds["z"][pcount[3]:pcount[4]])
    z4.attrs["Units"] = str(ds["z"].units)
    z4.attrs["PartNum"] = ptype[4]
    vx4 = part4.create_dataset("vx", data=ds["vx"][pcount[3]:pcount[4]])
    vx4.attrs["Units"] = str(ds["vx"].units)
    vx4.attrs["PartNum"] = ptype[4]
    vy4 = part4.create_dataset("vy", data=ds["vy"][pcount[3]:pcount[4]])
    vy4.attrs["Units"] = str(ds["vy"].units)
    vy4.attrs["PartNum"] = ptype[4]
    vz4 = part4.create_dataset("vz", data=ds["vz"][pcount[3]:pcount[4]])
    vz4.attrs["Units"] = str(ds["vz"].units)
    vz4.attrs["PartNum"] = ptype[4]
    mass4 = part4.create_dataset("Mass", data=ds["mass"].in_units('1e10 Msol')[pcount[3]:pcount[4]])
    mass4.attrs["Units"] = str(ds["mass"].in_units('1e10 Msol').units)
    mass4.attrs["PartNum"] = ptype[4]
    imass4 = part4.create_dataset("InitialMass", data=ds["InitialMass"].in_units('1e10 Msol')[pcount[3]:pcount[4]])
    imass4.attrs["Units"] = str(ds["InitialMass"].in_units('1e10 Msol').units)
    imass4.attrs["PartNum"] = ptype[4]
    age = part4.create_dataset("Age", data=ds.star["aform"])
    age.attrs["Units"] = "Expansion factor, a, where star particle forms"
    age.attrs["PartNum"] = ptype[4]
    Z = part4.create_dataset("Metallicity", data=ds.star["metals"])
    Z.attrs["Units"] = "Mass fraction of elements heavier than Helium"
    Z.attrs["PartNum"] = ptype[4]
    
hf.close()

print("New SimSpin file written at: ", args.outputFile)
