#!/bin/bash

# =========================================================================== #
#                    Automatic installation of Deformetrica                   #
# Description :                                                               #
# -------------                                                               #
#   This script will automatically install the CMake software together with   #
# the ITK & VTK libraries (if needed) and compile the project.                #
# =========================================================================== #


# Path to download the several softwares :
path_to_download_cmake="http://www.cmake.org/files/v2.8/cmake-2.8.12.2.tar.gz"
path_to_download_itk="http://sourceforge.net/projects/itk/files/itk/4.5/InsightToolkit-4.5.2.tar.gz/download"
path_to_download_vtk="http://www.vtk.org/files/release/6.1/VTK-6.1.0.tar.gz"



# OS name & architecture :
OSLOWER=`uname -s 2>/dev/null | tr "[:lower:]" "[:upper:]"` # (Linux or Darwin)
OS_ARCH=`uname -m | sed -e "s/i386/i686/"`

if [ $OS_ARCH != "x86_64" ] ; then 
    echo "Warning : Deformetrica has been tested only on 64 bit architecture (your architecture : $OS_ARCH)"
fi

if [ $OSLOWER != "LINUX" ] && [ $OSLOWER = "DARWIN" ]
then
     echo "Unknown operating system : the script will close."
     exit 1
fi



# CMake :
cmake_installed_or_uptodate=true
path_to_cmake=`command -v cmake`
if [ $? -eq 0 ] ; then # 'cmake' command found
     cmake_version=`cmake --version`
     cmake_version=${cmake_version#"cmake version "}
     cmake_version=${cmake_version:0:6}
     cmake_version=${cmake_version//./}
     if [ ${cmake_version} -lt 288 ] ; then
          echo "Version of CMake is too old"
          cmake_installed_or_uptodate=false
     fi
else
     cmake_installed_or_uptodate=false
fi

if [ $cmake_installed_or_uptodate == false ] ; then
     if [ $OSLOWER = "LINUX" ]
     then
          wget "$path_to_download_cmake" -O "cmake-latest.tar.gz"
     elif [ $OSLOWER = "DARWIN" ]
     then
          curl "$path_to_download_cmake" -o "cmake-latest.tar.gz"
     fi
     echo "Extracting CMake..."
     tar -zxf cmake-latest.tar.gz
     rm -vf cmake-latest.tar.gz
     folder_source_cmake=`ls | grep cmake-[0-9]*`
     if [ -d $folder_source_cmake ]; then
          echo "Installation of $folder_source_cmake"
          mkdir cmake
          path_to_cmake=`pwd`"/cmake/bin/cmake"
          cd $folder_source_cmake
          ./bootstrap --prefix=../cmake && make && make install
          # Cleaning :
          cd ../
          rm -rf $folder_source_cmake
     else
          echo "Error while extracting CMake"
          exit 1
     fi
fi



echo "##################################################"
echo "# Installation of ITK :"
echo "##################################################"
read -p 'Enter the path to ITK (if not installed, press ENTER)  : ' path_to_itk
if [ -z $path_to_itk ]; then
     if [ $OSLOWER = "LINUX" ]
     then
          wget "$path_to_download_itk" -O "InsightToolkit-latest.tar.gz"
     elif [ $OSLOWER = "DARWIN" ]
     then
          curl "$path_to_download_itk" -o "InsightToolkit-latest.tar.gz"
     fi
     if [ $? != 0 ]; then
          echo "Error when downloading ITK"
          exit 1
     else
          # Extracting ITK :
          echo "Extracting ITK..."
          tar -zxf InsightToolkit-latest.tar.gz
          rm -vf InsightToolkit-latest.tar.gz
          folder_source_itk=`ls | grep InsightToolkit-[0-9]*`
          if [ -d $folder_source_itk ]; then
               # Installation of ITK :
               mkdir ITKb;
               cd ITKb/;
               $path_to_cmake -D USE_FFTWD=ON -D USE_FFTWF=ON ../$folder_source_itk/
               make -j 2
               path_to_itk=`pwd`
               cd ../
          else
               echo "Error while extracting ITK."
               exit 1
          fi
     fi
fi



echo "##################################################"
echo "# Installation of VTK :"
echo "##################################################"
read -p 'Enter the path to VTK (if not installed, press ENTER)  : ' path_to_vtk
if [ -z $path_to_vtk ]; then
     if [ $OSLOWER = "LINUX" ]
     then
          wget "$path_to_download_vtk" -O "vtk-latest.tar.gz"
     elif [ $OSLOWER = "DARWIN" ]
     then
          curl "$path_to_download_vtk" -o "vtk-latest.tar.gz"
     fi
     if [ $? != 0 ]; then
          echo "Error when downloading VTK."
          exit 1
     else
          echo "Extracting VTK..."
          tar -zxf vtk-latest.tar.gz
          rm -vf vtk-latest.tar.gz
          folder_source_vtk=`ls | grep VTK*[0-9]*`
          if [ -d $folder_source_vtk ]; then
               # Installation of VTK :
               mkdir VTKb
               cd VTKb/
               $path_to_cmake -D Module_vtkWrappingTools=ON ../$folder_source_vtk/
               make -j 2
               path_to_vtk=`pwd`
               cd ..
          else
               echo "Error while extracting VTK."
               exit 1
          fi
     fi
fi
mkdir -p ../bin
cd ../bin/



echo "##################################################"
echo "# Compilation of deformetrica :"
echo "##################################################"
$path_to_cmake -D CMAKE_BUILD_TYPE=Release -D ITK_DIR=$path_to_itk -D VTK_DIR=$path_to_vtk ../app
make

