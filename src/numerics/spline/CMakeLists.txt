set( LIBNAME spline )

set( SOURCES
        cl_Spline.cpp
        fn_Create_Glue_Poly.cpp )

include_directories( ${BELFEM_SOURCE_DIR}/numerics/spline )

include( ${BELFEM_CONFIG_DIR}/scripts/Add_Library.cmake )