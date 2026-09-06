set( LIBNAME io )

set( SOURCES
        filetools.cpp
        hdf5_tools.cpp
        cl_HDF5.cpp
        cl_Ascii.cpp
        cl_CsvFile.cpp
        cl_InputFile.cpp
        cl_Input_Settings.cpp
        cl_Input_Section.cpp
        cl_XML.cpp )

include( ${BELFEM_CONFIG_DIR}/scripts/Add_Library.cmake )