import functions as func
import os

# CP_DIR = '/Volumes/WDElements/InternshipPartII/ProofOfConcept9M/CranioPatients/DataWeight100/MeanToCranio'
CP_DIR = '/Volumes/WDElements/InternshipPartII/ProofOfConcept9M/CranioPatients/DataWeight100/MeanToCranio'
CP_FILENAME = 'CP_final.txt'
MOM_FILENAME = 'Mom_final.txt'
PARAMS_FILENAME = 'paramDiffeos.xml'

# Read Momentum file paths in the CP_DIR
MOM_FILEPATH = func.pathToFile(CP_DIR,MOM_FILENAME)

# Read Control Points, Kernel-width, and compute Gaussian Kernel
CP = func.readCP(os.path.join(CP_DIR,CP_FILENAME))
kernelWidth = func.readKernelWidth(os.path.join(CP_DIR,PARAMS_FILENAME))
K = func.gaussianKernel(CP,kernelWidth)

# Read Momentum files and compute Geodesic Distances
for file in MOM_FILEPATH:
    MOM = func.readMOM(file)
    for mom in MOM:
        dist = func.dotProduct(mom,mom,K)
        print(file, dist)
