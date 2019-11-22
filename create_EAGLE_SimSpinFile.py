#!/usr/bin/env python

import numpy as np
import h5py
import os
import csv
from read_eagle import EagleSnapshot
from argparse import ArgumentParser

parser = ArgumentParser()

parser.add_argument("-i", "--infile", action="store", dest="inputEagleFile", type=str, default=None,
                    help="The path to the directory that contains the EAGLE particle data that you wish to process.") 

parser.add_argument("-n", "--nfiles", action="store", dest="inputNFiles", type=int, default=None,
                    help="The number of hdf5 files that make up the data at that snapshot.") 
                    
parser.add_argument("-d", "--indata", action="store", dest="inputData", type=str, default=None,
                    help="The path to the CSV file of particle locations within EAGLE.") 

parser.add_argument("-o", "--outfile", action="store", dest="outputFile", type=str, default=None,
                    help="The path to the SimSpin HDF5 file produced.")

# Parse arguments ###
args = parser.parse_args()


# Validate inputs ###

if args.inputEagleFile is None:
    print("\nERROR: No --infile argument given! \n "
          "For help, use \"python create_SimSpinFile.py -h\".\n"
          " --- ABORTING ---\n")
    exit()  # Check that an input file has been specified

if args.inputNFiles is None:
    print("\nERROR: No --nfiles argument given! \n "
          "For help, use \"python create_SimSpinFile.py -h\".\n"
          " --- ABORTING ---\n")
    exit()  # Check that the number of input particle files has been specified

if args.inputData is None:
    print("\nERROR: No --indata argument given! \n "
          "For help, use \"python create_SimSpinFile.py -h\".\n"
          " --- ABORTING ---\n")
    exit()  # Check that an input data file has been specified

if args.outputFile is None:
    print("\nERROR: No --outfile argument given! \n "
          "For help, use \"python create_SimSpinFile.py -h\".\n"
          " --- ABORTING ---\n")
    exit()  # Check that an output file has been specified

if os.path.exists(f"{args.inputEagleFile}.0.hdf5"):
    inputEagleFile = args.inputEagleFile
else:
    print("\nERROR: --infile does not exist!\n "
          "--- ABORTING ---\n")
    exit()  # Check that inputFile path leads to a real location

if os.path.exists(args.inputData):
    inputData = args.inputData
else:
    print("\nERROR: --indata does not exist!\n "
          "--- ABORTING ---\n")
    exit()  # Check that inputData path leads to a real file

# BEGIN ANALYSIS ---------------------------------------------------------------------------------

class SimSpinFile_ReadEagle:

    def __init__(self, inputEagleFile, nfiles, gn, sgn, centre, load_region_length):

        self.a, self.h, self.boxsize = self.read_header(inputEagleFile)

        # Load data
        self.data = self.read_galaxy(inputEagleFile, nfiles, gn, sgn, centre, load_region_length)

    def read_header(self, inputEagleFile):
        # Read various attributes from the header group 
        f = h5py.File(f"{inputEagleFile}.0.hdf5", 'r')
        self.a = f['Header'].attrs.get('Time') # Scale factor, a
        self.h = f['Header'].attrs.get('HubbleParam') # h
        self.boxsize = f["Header"].attrs.get("BoxSize")
        f.close()
        return self.a, self.h, self.boxsize
    
    def read_galaxy(self, inputEagleFile, nfiles, gn, sgn, centre, load_region_length):
        # For a given galaxy (defined by its GroupNumber and SubGroupNumber)
        # extract the coordinatesm velocitiesm masses, age and metallicity
        # using the read_eagle routine. Conversion factors are still loaded directly.
        # from the hdf5 files.
        self.data = {}

        # Put centre into units of cMpc/h
        centre = centre * self.h

        # Select region to load
        region = np.array([
            (centre[0]-0.5*load_region_length), (centre[0]+0.5*load_region_length),
            (centre[1]-0.5*load_region_length), (centre[1]+0.5*load_region_length),
            (centre[2]-0.5*load_region_length), (centre[2]+0.5*load_region_length) 
        ])

        # Initialise the read_eagle module
        for i in range(nfiles+1):
            eagle_data = EagleSnapshot(f"{inputEagleFile}.{i}.hdf5")
            eagle_data.select_region(*region)

            for att in ["GroupNumber", "SubGroupNumber", "Coordinates", "Velocity", "Mass", "StellarFormationTime", "SmoothedMetallicity"]:
                tmp = eagle_data.read_dataset(4, att)
                if i == 0:
                    self.data[att] = tmp
                else:
                    self.data[att] = np.append(self.data[att], tmp, axis = 0)
            
            eagle_data.close()

            print(f"Done reading EAGLE snapshot file {i}.")

        # Mask to selected GroupNumber and SubGroupNumber.
        mask = np.logical_and(self.data['GroupNumber'] == gn, self.data['SubGroupNumber'] == sgn)
        for att in self.data.keys():
            self.data[att] = self.data[att][mask]
            
        # Periodic wrap coordiantes around centre
        boxsize = self.boxsize/self.h
        self.data["Coordinates"] = np.mod(self.data['Coordinates'] - centre+0.5*boxsize, boxsize) + centre-0.5*boxsize

        # Put back to simulation units
        self.data["Coordinates"] = np.multiply(self.data["Coordinates"], ((self.a**1) * (self.h**-1) * 1000), dtype="f8") # putting coordinates back into kpc
        self.data["Velocity"] = np.multiply(self.data["Velocity"], (self.a**0.5), dtype="f8") # putting velocity into km/s
        self.data["Mass"] = np.multiply(self.data["Mass"], (self.h**-1), dtype="f8") # putting mass into 1e10 Msolar

        return self.data

print("Input file location: ", args.inputEagleFile)
print("Number of files to loop over: ", (args.inputNFiles+1))
print("Input data: ", args.inputData)

galaxyID, group_number, subgroup_number, cop_x, cop_y, cop_z, cut_out_region = ([] for i in range(7))

with open(args.inputData) as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=",")
    line_count = 0
    for row in csv_reader:
        if line_count == 0:
            print(f'Input galaxy data supplied: {", ".join(row)}')
            line_count += 1
        else:
            galaxyID.append(row[0])
            group_number.append(int(row[1]))
            subgroup_number.append(int(row[2]))
            cop_x.append(float(row[3]))
            cop_y.append(float(row[4]))
            cop_z.append(float(row[5]))
            cut_out_region.append(float(row[9]))


no_galaxies = len(galaxyID)
print(f"There are {no_galaxies} galaxies contained in the request file. Beginning to read EAGLE snapshots...")

for j in range(no_galaxies):

    stars = SimSpinFile_ReadEagle(inputEagleFile, args.inputNFiles, gn=group_number[j], sgn=subgroup_number[j], centre=np.array([cop_x[j], cop_y[j], cop_z[j]]), load_region_length=cut_out_region[j])

    hf = h5py.File(f"{args.outputFile}/SimSpin_galaxyID{galaxyID[j]}.hdf5", "w")    # creating HDF5 file for SimSpin
    part4 = hf.create_group("PartType4")
    x4 = part4.create_dataset("x", data=stars.data["Coordinates"][:,0])
    x4.attrs["Units"] = "kpc"
    y4 = part4.create_dataset("y", data=stars.data["Coordinates"][:,1])
    y4.attrs["Units"] = "kpc"
    z4 = part4.create_dataset("z", data=stars.data["Coordinates"][:,2])
    z4.attrs["Units"] = "kpc"
    vx4 = part4.create_dataset("vx", data=stars.data["Velocity"][:,0])
    vx4.attrs["Units"] = "km s-1"
    vy4 = part4.create_dataset("vy", data=stars.data["Velocity"][:,1])
    vy4.attrs["Units"] = "km s-1"
    vz4 = part4.create_dataset("vz", data=stars.data["Velocity"][:,2])
    vz4.attrs["Units"] = "km s-1"
    mass4 = part4.create_dataset("Mass", data=stars.data["Mass"])
    mass4.attrs["Units"] = "1e10 M_solar"
    sft4 = part4.create_dataset("StellarFormationTime", data=stars.data["StellarFormationTime"])
    sft4.attrs["Units"] = "Expansion factor, a, at time of star birth."
    metal4 = part4.create_dataset("Metallicity", data=stars.data["SmoothedMetallicity"])
    metal4.attrs["Units"] = "Mass fraction of elements heavier than Helium"

    hf.close()
    print(f"Done reading galaxy {j+1} of {no_galaxies}. New SimSpin File written: SimSpin_galaxyID{galaxyID[j]}.hdf5")

print(f"New SimSpin files written at: {args.outputFile}")


        



