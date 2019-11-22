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

if os.path.exists(f"{args.inputEagleFile}"):
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

class CreateFile_ReadEagle:

    def __init__(self, inputEagleFile, gn, sgn, centre, load_region_length):

        self.a, self.h, self.boxsize = self.read_header(inputEagleFile)

        # Load data
        self.data, self.region, self.numpart_total = self.read_galaxy(inputEagleFile, gn, sgn, centre, load_region_length)

    def read_header(self, inputEagleFile):
        # Read various attributes from the header group 
        f = h5py.File(f"{inputEagleFile}", 'r')
        self.a = f['Header'].attrs.get('Time') # Scale factor, a
        self.h = f['Header'].attrs.get('HubbleParam') # h
        self.boxsize = f["Header"].attrs.get("BoxSize")
        f.close()
        return self.a, self.h, self.boxsize
    
    def read_galaxy(self, inputEagleFile, gn, sgn, centre, load_region_length):
        # For a given galaxy (defined by its GroupNumber and SubGroupNumber)
        # extract the coordinatesm velocitiesm masses, age and metallicity
        # using the read_eagle routine. Conversion factors are still loaded directly.
        # from the hdf5 files.
        self.data = {}

        # Put centre into units of cMpc/h
        centre = centre * self.h
        boxsize = self.boxsize/self.h

        # Select region to load
        self.region = np.array([
            (centre[0]-0.5*load_region_length), (centre[0]+0.5*load_region_length),
            (centre[1]-0.5*load_region_length), (centre[1]+0.5*load_region_length),
            (centre[2]-0.5*load_region_length), (centre[2]+0.5*load_region_length) 
        ])

        self.numpart_total = np.zeros(6, dtype="uint64")

        # Initialise the read_eagle module
        
        eagle_data = EagleSnapshot(f"{inputEagleFile}")
        eagle_data.select_region(*self.region)

        # Loop over particle types to process
        for itype in range(6):
            # Get number of particles to read
            nop = eagle_data.count_particles(itype) 

            if nop > 0: # if there are particles in this particle type within the region
                mask = np.logical_and(eagle_data.read_dataset(itype, '/GroupNumber') == gn, eagle_data.read_dataset(itype, '/SubGroupNumber') == sgn)
                nopig = np.sum(mask)
                if nopig > 0: # if there are particles of this type within the galaxy
                    print(f"Particle type {itype} - keeping {nopig} particles in galaxy from EAGLE snapshot file")

                    self.data[f"PartType{itype}"] = {}

                    read_datasets = eagle_data.datasets(itype)

                    for dset_name in read_datasets:
                        tmp = eagle_data.read_dataset(itype, dset_name)[mask]
                        self.data[f"PartType{itype}"][dset_name] = tmp

                    
        eagle_data.close()    
        

        f = h5py.File(f"{inputEagleFile}", "r")
        for ptype in self.data.keys():
            # Periodic wrap coordiantes around centre
            self.data[f"{ptype}"]["/Coordinates"] = np.mod(self.data[f"{ptype}"]['/Coordinates'] - centre+0.5*boxsize, boxsize) + centre-0.5*boxsize

            # Put back to simulation units
            #for att in self.data[f"{ptype}"].keys():
            #    cgs = f[f"{ptype}{att}"].attrs.get('CGSConversionFactor')
            #    aexp = f[f"{ptype}{att}"].attrs.get("aexp-scale-exponent")
            #    hexp = f[f"{ptype}{att}"].attrs.get("h-scale-exponent")
            #    self.data[f"{ptype}"][att] = np.multiply(self.data[f"{ptype}"][att], cgs * self.a**aexp * self.h**hexp, dtype="f8")
            
            self.numpart_total[int(ptype[8])] = len(self.data[f"{ptype}"]["/ParticleIDs"])

        f.close()

        return self.data, self.region, self.numpart_total 


print("Input file location: ", inputEagleFile)
print("Input data: ", inputData)

galaxyID, group_number, subgroup_number, cop_x, cop_y, cop_z, cut_out_region = ([] for i in range(7))

with open(inputData) as csv_file:
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

    gal_data = CreateFile_ReadEagle(inputEagleFile, gn=group_number[j], sgn=subgroup_number[j], centre=np.array([cop_x[j], cop_y[j], cop_z[j]]), load_region_length=cut_out_region[j])
    types = (1,1,1,1,1,1)
    
    f = h5py.File(f"{inputEagleFile}", 'r')
    hf = h5py.File(f"{args.outputFile}/EAGLE_galaxyID{galaxyID[j]}.hdf5", "w")    # creating HDF5 file for SimSpin

    header = hf.create_group("Header")
    for (name,val) in f["Header"].attrs.items():
        header.attrs[name] = val

    # Update particle numbers in header
    nptot    = np.zeros(6, dtype="uint32")
    nptot_hw = np.zeros(6, dtype="uint32")
    nptot_hw[:] = gal_data.numpart_total >> 32
    nptot[:]    = gal_data.numpart_total - (nptot_hw << 32)
    header.attrs["NumPart_Total"] = nptot
    header.attrs["NumPart_Total_HighWord"] = nptot_hw
    header.attrs["NumPart_ThisFile"] = nptot

    # Now only have a single file
    header.attrs["NumFilesPerSnapshot"] = 1

    # Copy other groups with run information
    for group_name in ("Config",
                       "Constants",
                       "Parameters",
                       "Parameters/ChemicalElements",
                       "RuntimePars",
                       "Units"):
        group = hf.create_group(group_name)
        for (name,val) in f[group_name].attrs.items():
            group.attrs[name] = val

    header.attrs["ExtractedFromSnapshot"] = 28

    # Add region spec and type flags
    header.attrs["RegionExtracted"] = np.array(gal_data.region, dtype="float64")
    header.attrs["TypesExtracted"]  = np.array(types, dtype="int32")
    header.attrs["SamplingRate"] = 1.0

    for ptype in gal_data.data.keys():
        if gal_data.numpart_total[int(ptype[8])] > 0: # if there are any ptype particles associated with the galaxy
            for att in gal_data.data[ptype].keys():
                chunks = [s for s in gal_data.data[ptype][att].shape]
                chunks[0] = min((8192, chunks[0]))
                hf.create_dataset(f"{ptype}/{att}", 
                                data=gal_data.data[ptype][att],
                                chunks=tuple(chunks),
                                shuffle=True,
                                compression="gzip",
                                compression_opts=6)
                for (name,val) in f[f"{ptype}/{att}"].attrs.items():
                    hf[f"{ptype}/{att}"].attrs[name] = val


    f.close()
    hf.close()
    print(f"Done reading galaxy {j+1} of {no_galaxies}. New EAGLE particle file written: EAGLE_galaxyID{galaxyID[j]}.hdf5")
    print()
    
print()
print(f"New SimSpin files written at: {args.outputFile}")
print()


        



