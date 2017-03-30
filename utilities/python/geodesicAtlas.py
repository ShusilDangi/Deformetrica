import functions as func
import numpy as np
import os
import csv

CP_DIR = '/Volumes/WDElements/InternshipPartII/ProofOfConcept9M/Atlas30'
CP_FILENAME = 'CP_final.txt'
MOM_FILENAME = 'Mom_final.txt'
PARAMS_FILENAME = 'paramDiffeos.xml'
ATLAS_FILENAME = 'atlasFilesOrder.txt'

# Read the file order in which the moments are stored in the atlas
AtlasFileOrder = []
with open(os.path.join(CP_DIR,ATLAS_FILENAME), 'rb') as csvfile:
    filereader = csv.reader(csvfile, delimiter=',', quotechar='|')
    for row in filereader:
        for fileName in row:
            AtlasFileOrder.append(fileName)

# Read Control Points, Kernel-width, and compute Gaussian Kernel
CP = func.readCP(os.path.join(CP_DIR,CP_FILENAME))
kernelWidth = func.readKernelWidth(os.path.join(CP_DIR,PARAMS_FILENAME))
K = func.gaussianKernel(CP,kernelWidth)

# Read moments for all normal patients in Atlas from the mean shape
MOM = func.readMOMAtlas(os.path.join(CP_DIR,MOM_FILENAME))
dist = []
for mom in MOM:
    dist.append(func.dotProduct(mom,mom,K))     # Geodesic Distance

# Sort the files according to the Geodesic Distance from the mean shape
dist = np.asarray(dist)
sortInd = np.argsort(dist)
for ind in sortInd:
    print(AtlasFileOrder[ind],dist[ind])
