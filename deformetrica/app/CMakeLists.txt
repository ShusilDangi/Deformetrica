CMake_Minimum_Required(VERSION 2.8.7)

Project(Deformetrica)
Set (DEFORMETRICA_VERSION_MAJOR 2)
Set (DEFORMETRICA_VERSION_MINOR 1)

If (COMMAND cmake_policy)
    cmake_policy(SET CMP0003 NEW)
Endif (COMMAND cmake_policy)



# ITK :
Find_Package(ITK)
If (ITK_FOUND)
    Include (${ITK_USE_FILE})
Else (ITK_FOUND)
    Message (FATAL_ERROR "ITK not found. Please set ITK_DIR.")
Endif (ITK_FOUND)

# VTK :
Find_Package(VTK)
If (VTK_FOUND)
    If ( "${VTK_MAJOR_VERSION}" LESS 6 )
        Message ("VTK ${VTK_MAJOR_VERSION}.${VTK_MINOR_VERSION} is used.")
        Include (${USE_VTK_FILE})
    Else ( "${VTK_MAJOR_VERSION}" LESS 6 )
        Find_Package(VTK COMPONENTS vtkCommonCore vtkCommonDataModel vtkWrappingTools vtkIOLegacy vtkFiltersCore NO_MODULE)
        If (VTK_FOUND)
            Message ("VTK ${VTK_MAJOR_VERSION}.${VTK_MINOR_VERSION} is used.")
            Include (${VTK_USE_FILE})
        Else (VTK_FOUND)
            Message (FATAL_ERROR "Some components of VTK have not been found.")
        Endif (VTK_FOUND)
    Endif ( "${VTK_MAJOR_VERSION}" LESS 6 )
Else (VTK_FOUND)
    Message (FATAL_ERROR "VTK not found. Please set VTK_DIR.")
Endif (VTK_FOUND)

# Cuda :
Find_Package(CUDA)
If (CUDA_FOUND)
    Option (USE_CUDA "Use cuda implementation of convolutions" ON) 
Else (CUDA_FOUND)
    Message ("CUDA was not found.")
Endif (CUDA_FOUND)

Include_Directories(
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/../src/algo
    ${CMAKE_CURRENT_SOURCE_DIR}/../src/algo/deformable_objects
    ${CMAKE_CURRENT_SOURCE_DIR}/../src/algo/deformations
    ${CMAKE_CURRENT_SOURCE_DIR}/../src/algo/energy_minimization
    ${CMAKE_CURRENT_SOURCE_DIR}/../src/algo/integrators
    ${CMAKE_CURRENT_SOURCE_DIR}/../src/algo/kernels
    ${CMAKE_CURRENT_SOURCE_DIR}/../src/algo/utilities
    ${CMAKE_CURRENT_SOURCE_DIR}/../src/io
    ${PROJECT_BINARY_DIR}
)

# Release version by default :
If (NOT DEFINED CMAKE_BUILD_TYPE:STRING)
    Set (CMAKE_BUILD_TYPE:STRING Release)
Endif (NOT DEFINED CMAKE_BUILD_TYPE:STRING)

# Extra flags for the C++ compiler :
If (CMAKE_BUILD_TYPE MATCHES Release)
    Set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2")
Endif (CMAKE_BUILD_TYPE MATCHES Release)


# Set computations on double precision by default :
If (NOT DEFINED SET_DOUBLE_PRECISION)
    Set (SET_DOUBLE_PRECISION YES)
Endif (NOT DEFINED SET_DOUBLE_PRECISION)
If (SET_DOUBLE_PRECISION)
    Option (USE_DOUBLE_PRECISION "Deformetrica will make computations in double precision" ON) 
Else (SET_DOUBLE_PRECISION)
    Message ("Computations will be in simple precision")
Endif (SET_DOUBLE_PRECISION)


# Configure a header file to pass some of the CMake settings to the source code :
Configure_File (
    "${PROJECT_SOURCE_DIR}/DeformetricaConfig.h.in"
    "${PROJECT_BINARY_DIR}/DeformetricaConfig.h"
)


If (USE_CUDA)

    Message("Deformetrica will be compiled with Cuda")

    Set( MY_CUDA_FLAGS -gencode arch=compute_11,code=sm_11;-gencode arch=compute_12,code=sm_12;-gencode arch=compute_13,code=sm_13;-gencode arch=compute_20,code=sm_20;-gencode arch=compute_30,code=sm_30;-gencode arch=compute_35,code=sm_35)
    List(APPEND CUDA_NVCC_FLAGS ${MY_CUDA_FLAGS})

    # Library containing the Cuda operations
    Set (LIB_TYPE STATIC)
    Cuda_Add_Library (
    CudaConvolution ${LIB_TYPE}
        ${CMAKE_CURRENT_SOURCE_DIR}/../lib/cuda_convolutions/GaussFunction.h
        ${CMAKE_CURRENT_SOURCE_DIR}/../lib/cuda_convolutions/RadialFunction.h
        ${CMAKE_CURRENT_SOURCE_DIR}/../lib/cuda_convolutions/ScalarRadialKernel.h
        ${CMAKE_CURRENT_SOURCE_DIR}/../lib/cuda_convolutions/GpuConv1D.cu
    )
    Set (CUDA_LIBS ${CUDA_LIBS} CudaConvolution)

Cuda_Add_Executable(sparseMatching2
  ../src/io/SparseDiffeoParametersXMLFile.cxx
  ../src/io/DeformableObjectParametersXMLFile.cxx
  ../src/io/SparseDiffeoParameters.cxx
  ../src/io/DeformableObjectParameters.cxx
  ../src/algo/utilities/SimpleTimer.cxx
  ../src/algo/utilities/myvtkPolyDataNormals.cxx
  matching/sparseMatching2.cxx
)

Cuda_Add_Executable(sparseMatching3
  ../src/io/SparseDiffeoParametersXMLFile.cxx
  ../src/io/DeformableObjectParametersXMLFile.cxx
  ../src/io/SparseDiffeoParameters.cxx
  ../src/io/DeformableObjectParameters.cxx
  ../src/algo/utilities/SimpleTimer.cxx
  ../src/algo/utilities/myvtkPolyDataNormals.cxx
  matching/sparseMatching3.cxx
)

Cuda_Add_Executable(sparseAtlas2
  ../src/io/SparseDiffeoParametersXMLFile.cxx
  ../src/io/DeformableObjectParametersXMLFile.cxx
  ../src/io/SparseDiffeoParameters.cxx
  ../src/io/DeformableObjectParameters.cxx
  ../src/algo/utilities/SimpleTimer.cxx
  ../src/algo/utilities/myvtkPolyDataNormals.cxx
  atlas_building/sparseAtlas2.cxx
)

Cuda_Add_Executable(sparseAtlas3
  ../src/io/SparseDiffeoParametersXMLFile.cxx
  ../src/io/DeformableObjectParametersXMLFile.cxx
  ../src/io/SparseDiffeoParameters.cxx
  ../src/io/DeformableObjectParameters.cxx
  ../src/algo/utilities/SimpleTimer.cxx
  ../src/algo/utilities/myvtkPolyDataNormals.cxx
  atlas_building/sparseAtlas3.cxx
)

Cuda_Add_Executable(ShootAndFlow2
  ../src/io/SparseDiffeoParametersXMLFile.cxx
  ../src/io/DeformableObjectParametersXMLFile.cxx
  ../src/io/SparseDiffeoParameters.cxx
  ../src/io/DeformableObjectParameters.cxx
  ../src/algo/utilities/myvtkPolyDataNormals.cxx
  shooting/ShootAndFlow2.cxx
)

Cuda_Add_Executable(ShootAndFlow3
  ../src/io/SparseDiffeoParametersXMLFile.cxx
  ../src/io/DeformableObjectParametersXMLFile.cxx
  ../src/io/SparseDiffeoParameters.cxx
  ../src/io/DeformableObjectParameters.cxx
  ../src/algo/utilities/myvtkPolyDataNormals.cxx
  shooting/ShootAndFlow3.cxx
)

Cuda_Add_Executable(geodesicDist
  ../src/io/SparseDiffeoParametersXMLFile.cxx
  ../src/io/DeformableObjectParametersXMLFile.cxx
  ../src/io/SparseDiffeoParameters.cxx
  ../src/io/DeformableObjectParameters.cxx
  ../src/algo/utilities/myvtkPolyDataNormals.cxx
  geodesicDist/geodesicDist.cxx
)

Cuda_Add_Executable(geodesicDistMulti
  ../src/io/SparseDiffeoParametersXMLFile.cxx
  ../src/io/DeformableObjectParametersXMLFile.cxx
  ../src/io/SparseDiffeoParameters.cxx
  ../src/io/DeformableObjectParameters.cxx
  ../src/algo/utilities/myvtkPolyDataNormals.cxx
  geodesicDist/geodesicDistMulti.cxx
)

else(USE_CUDA)

    Message("Deformetrica will be compiled without Cuda")

Add_Executable(sparseMatching2
  ../src/io/SparseDiffeoParametersXMLFile.cxx
  ../src/io/DeformableObjectParametersXMLFile.cxx
  ../src/io/SparseDiffeoParameters.cxx
  ../src/io/DeformableObjectParameters.cxx
  ../src/algo/utilities/SimpleTimer.cxx
  ../src/algo/utilities/myvtkPolyDataNormals.cxx
  matching/sparseMatching2.cxx
)

Add_Executable(sparseMatching3
  ../src/io/SparseDiffeoParametersXMLFile.cxx
  ../src/io/DeformableObjectParametersXMLFile.cxx
  ../src/io/SparseDiffeoParameters.cxx
  ../src/io/DeformableObjectParameters.cxx
  ../src/algo/utilities/SimpleTimer.cxx
  ../src/algo/utilities/myvtkPolyDataNormals.cxx
  matching/sparseMatching3.cxx
)

Add_Executable(sparseAtlas2
  ../src/io/SparseDiffeoParametersXMLFile.cxx
  ../src/io/DeformableObjectParametersXMLFile.cxx
  ../src/io/SparseDiffeoParameters.cxx
  ../src/io/DeformableObjectParameters.cxx
  ../src/algo/utilities/SimpleTimer.cxx
  ../src/algo/utilities/myvtkPolyDataNormals.cxx
  atlas_building/sparseAtlas2.cxx
)

Add_Executable(sparseAtlas3
  ../src/io/SparseDiffeoParametersXMLFile.cxx
  ../src/io/DeformableObjectParametersXMLFile.cxx
  ../src/io/SparseDiffeoParameters.cxx
  ../src/io/DeformableObjectParameters.cxx
  ../src/algo/utilities/SimpleTimer.cxx
  ../src/algo/utilities/myvtkPolyDataNormals.cxx
  atlas_building/sparseAtlas3.cxx
)

Add_Executable(ShootAndFlow2
    ../src/io/SparseDiffeoParametersXMLFile.cxx
    ../src/io/DeformableObjectParametersXMLFile.cxx
    ../src/io/SparseDiffeoParameters.cxx
    ../src/io/DeformableObjectParameters.cxx
    ../src/algo/utilities/myvtkPolyDataNormals.cxx
    shooting/ShootAndFlow2.cxx
)

Add_Executable(ShootAndFlow3
  ../src/io/SparseDiffeoParametersXMLFile.cxx
  ../src/io/DeformableObjectParametersXMLFile.cxx
  ../src/io/SparseDiffeoParameters.cxx
  ../src/io/DeformableObjectParameters.cxx
  ../src/algo/utilities/myvtkPolyDataNormals.cxx
  shooting/ShootAndFlow3.cxx
)

Add_Executable(geodesicDist
  ../src/io/SparseDiffeoParametersXMLFile.cxx
  ../src/io/DeformableObjectParametersXMLFile.cxx
  ../src/io/SparseDiffeoParameters.cxx
  ../src/io/DeformableObjectParameters.cxx
  ../src/algo/utilities/myvtkPolyDataNormals.cxx
  geodesicDist/geodesicDist.cxx
)

Add_Executable(geodesicDistMulti
  ../src/io/SparseDiffeoParametersXMLFile.cxx
  ../src/io/DeformableObjectParametersXMLFile.cxx
  ../src/io/SparseDiffeoParameters.cxx
  ../src/io/DeformableObjectParameters.cxx
  ../src/algo/utilities/myvtkPolyDataNormals.cxx
  geodesicDist/geodesicDistMulti.cxx
)

endif(USE_CUDA)



# Target link :
Set ( VTK_DEPENDANCIES "" )
If ( "${VTK_MAJOR_VERSION}" LESS 6 )
    List (APPEND VTK_DEPENDANCIES
        vtkCommon
        vtkGraphics
        vtkIO
        vtkRendering
        vtkFiltering
    )
Else ( "${VTK_MAJOR_VERSION}" LESS 6 )
    List (APPEND VTK_DEPENDANCIES
        ${VTK_LIBRARIES}
    )
Endif ( "${VTK_MAJOR_VERSION}" LESS 6 )


Target_Link_Libraries(sparseMatching2
  ${ITK_LIBRARIES}
  ${VTK_DEPENDANCIES}
  ${CUDA_LIBS}
)

Target_Link_Libraries(sparseMatching3
  ${ITK_LIBRARIES}
  ${VTK_DEPENDANCIES}
  ${CUDA_LIBS}
)

Target_Link_Libraries(sparseAtlas2
  ${ITK_LIBRARIES}
  ${VTK_DEPENDANCIES}
  ${CUDA_LIBS}
)

Target_Link_Libraries(sparseAtlas3
  ${ITK_LIBRARIES}
  ${VTK_DEPENDANCIES}
  ${CUDA_LIBS}
)

Target_Link_Libraries(ShootAndFlow2
  ${ITK_LIBRARIES}
  ${VTK_DEPENDANCIES}
  ${CUDA_LIBS}
)

Target_Link_Libraries(ShootAndFlow3
  ${ITK_LIBRARIES}
  ${VTK_DEPENDANCIES}
  ${CUDA_LIBS}
)

Target_Link_Libraries(geodesicDist
  ${ITK_LIBRARIES}
  ${VTK_DEPENDANCIES}
  ${CUDA_LIBS}
)

Target_Link_Libraries(geodesicDistMulti
  ${ITK_LIBRARIES}
  ${VTK_DEPENDANCIES}
  ${CUDA_LIBS}
)


# Dependencies (only for VTK < 6) :
If ( "${VTK_MAJOR_VERSION}" LESS 6 )
Add_Dependencies (sparseMatching2
  ${ITK_LIBRARIES}
  ${VTK_DEPENDANCIES}
)

Add_Dependencies (sparseMatching3
  ${ITK_LIBRARIES}
  ${VTK_DEPENDANCIES}
)

Add_Dependencies (sparseAtlas2
  ${ITK_LIBRARIES}
  ${VTK_DEPENDANCIES}
)

Add_Dependencies (sparseAtlas3
  ${ITK_LIBRARIES}
  ${VTK_DEPENDANCIES}
)

Add_Dependencies (ShootAndFlow2
  ${ITK_LIBRARIES}
  ${VTK_DEPENDANCIES}
)

Add_Dependencies (ShootAndFlow3
  ${ITK_LIBRARIES}
  ${VTK_DEPENDANCIES}
)

Add_Dependencies (geodesicDist
  ${ITK_LIBRARIES}
  ${VTK_DEPENDANCIES}
)

Add_Dependencies (geodesicDistMulti
  ${ITK_LIBRARIES}
  ${VTK_DEPENDANCIES}
)
Endif ( "${VTK_MAJOR_VERSION}" LESS 6 )



# Testing With CTest :
Enable_Testing()

Add_Test(
  NAME surface_matching
  WORKING_DIRECTORY ../../tests/surface_matching/
  COMMAND sparseMatching3 paramDiffeos.xml paramSurface.xml sourceSurfaceMesh.vtk targetSurfaceMesh.vtk
)

Add_Test(
  NAME atlas_building
  WORKING_DIRECTORY ../../tests/atlas_building/
  COMMAND sparseAtlas3 paramDiffeos.xml 1 paramSurface.xml template.vtk banana_1.vtk banana_2.vtk banana_3.vtk
)





