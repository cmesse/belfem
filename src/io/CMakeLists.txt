set( LIBNAME io )

set( SOURCES
        filetools.cpp
        HDF5_Tools.cpp
        cl_HDF5.cpp
        cl_Ascii.cpp
        cl_InputFile.cpp
        cl_Input_Section.cpp
        cl_XML.cpp )

include( ${BELFEM_CONFIG_DIR}/scripts/Add_Library.cmake )