# Deformetrica
Diffeomorphic Registration using Deformetrica Library with some additional functionality and python utility files

Refer to the Deformetrica website for Installation details, User Manuals, Tutorials and more here:
http://www.deformetrica.org/

Added Utility Functions:
-------------------------------------------------------------------------
geodesicDist
-------------------------------------------------------------------------
This command computes the geodesic distance between the two subjects. User should provide the XML file with diffeomorphic 
parameters (paramDiffeos.xml), final position of the control points (CP.txt), and the computed moments (MOM.txt). 
The syntax is:

  `geodesicDist paramDiffeos.xml CP.txt MOM.txt`

The geodesic distance is given by
geoDist = sum(mom\*K\*mom)

mom - Momenta, K - Gaussian Kernel

---------------------------------------------------------------------------
geodesicDistAtlas
---------------------------------------------------------------------------
This command computes the geodesic distance between the mean subject and the subjects used to construct the atlas. User should
provide the XML file with diffeomorphic parameters (paramDiffeos.xml), final position of the control points (CP.txt), and the 
computed moments for all subjects in a single file (MOM.txt). The syntax is:

`geodesicDistAtlas paramDiffeos.xml CP.txt MOM.txt`

----------------------------------------------------------------------------
Regularization weights for the Atlas Computation
----------------------------------------------------------------------------
Sum of geodesic distances from the mean shape to all the subjects in the atlas is used to regularize the Atlas construction.
By default, the geodesic distances are weighted equally, forcing the mean shape to be equi-distant from all the subjects. 
The regularization weights will force the mean shape to be closer to the subjects with higher weights. 
Typical use of this weighting scheme could be to regress a shape using Weighted Nearest Neighbor scheme.

The weights should to be provided through the "paramDiffeo.xml" file under the field "regularity-weights"

e.g. `<regularity-weights>0.1,0.2,0.3,0.2,0.1</regularity-weights>`

If the total number of regularity-weights does not match the number of subjects used to construct the Atlas, default uniform
weighting scheme will be used.

-----------------------------------------------------------------------------
Python Utility Files
-----------------------------------------------------------------------------
Functions.py
-----------------------------------------------------------------------------
Contains several utility functions:
- Read Control Point file
- Read Moment File
- Read Atlas Moments File
- Read Kernel Width from the paramDiffeos.xml file
- Add fields to the XML file (add <regularity-weights> to 'paramDiffeo.xml')
- Compute the Gaussian Kernel from Control Points and Kernel Width
- Convert Moments to Velocity field
- Convert Velocity field to Moments
- Compute the Geodesic Distance using Gaussian Kernel and Moments
- Compute the dot product in Kernel space
- Read the directory path to specified FileName within the RootDirectory
- Find closest Normal subjects in Atlas from the outlier patient based on approximate geodesic distances

------------------------------------------------------------------------------
geodesicCranio.py, geodesicAtlas.py
------------------------------------------------------------------------------
File demonstrating the computation of geodesic distances.

Between two subjects -> geodesicCranio.py

Mean to Atlas subjects -> geodesicAtlas.py

