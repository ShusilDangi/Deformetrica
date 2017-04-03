CMake_Minimum_Required(VERSION 2.8.7)

Project(SurfUtil)

If (COMMAND cmake_policy)
    cmake_policy(SET CMP0003 NEW)
Endif (COMMAND cmake_policy)

#SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/CMake")

# ITK :
Find_Package(ITK)
If (ITK_FOUND)
    Include (${USE_ITK_FILE})
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

# Release version by default :
If (NOT DEFINED CMAKE_BUILD_TYPE:STRING)
    Set (CMAKE_BUILD_TYPE:STRING Release)
Endif (NOT DEFINED CMAKE_BUILD_TYPE:STRING)

# Extra flags for the C++ compiler :
If (CMAKE_BUILD_TYPE MATCHES Release)
    Set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2")
Endif (CMAKE_BUILD_TYPE MATCHES Release)





# Executable :
Add_Executable(MapsEllipsoidWithSource MapsEllipsoidWithSource.cxx)

# Target link library :
Target_Link_Libraries(MapsEllipsoidWithSource
    ${ITK_LIBRARIES}
    ${VTK_LIBRARIES}
)


