import functions as func
import os

# CP_DIR = '/Volumes/WDElements/InternshipPartII/ProofOfConcept9M/CranioPatients/MeanToCranio_Data1'
CP_DIR = '/Volumes/WDElements/InternshipPartII/ProofOfConcept9M/CranioPatients/MeanToCranio_Data001'
CP_FILENAME = 'CP_final.txt'
MOM_FILENAME = 'Mom_final.txt'
PARAMS_FILENAME = 'paramDiffeos.xml'

# Read Momentum file paths in the CP_DIR
MOM_DIR = func.pathToFile(CP_DIR,MOM_FILENAME)

# Read Control Points, Kernel-width, and compute Gaussian Kernel
CP = func.readCP(os.path.join(CP_DIR,CP_FILENAME))
kernelWidth = func.readKernelWidth(os.path.join(CP_DIR,PARAMS_FILENAME))
K = func.gaussianKernel(CP,kernelWidth)

# Read Momentum files and compute Geodesic Distances
for d in MOM_DIR:
    MOM = func.readMOM(os.path.join(d,MOM_FILENAME))
    for mom in MOM:
        dist = func.dotProduct(mom,mom,K)		# Geodesic Distance
        print(MOM_FILENAME,dist)
