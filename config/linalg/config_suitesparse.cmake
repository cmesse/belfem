if( USE_SUITESPARSE )

    message(WARNING
        "USE_SUITESPARSE=ON links against SuiteSparse/UMFPACK, which is partly "
        "licensed under GPL-2.0-or-later. Any binary built with this option "
        "becomes a combined work governed by the GPL, not by BELFEM's BSD "
        "license, and must NOT be redistributed under BELFEM's terms. This "
        "option is intended for non-distributed local builds only. For "
        "redistributable builds, use the SuperLU solver (BSD-3-Clause).")

    belfem_find_package(
            SUITESPARSE
            REQUIRED
            HEADERS umfpack.h
            LIBRARIES umfpack klu cholmod ccolamd colamd camd btf amd suitesparseconfig
            INCLUDE_SUFFIXES suitesparse . )

    list( APPEND BELFEM_DEFS "BELFEM_SUITESPARSE" )
    list( APPEND BELFEM_INCLUDES "${SUITESPARSE_INCLUDE_DIR}" )
    list( APPEND BELFEM_RPATH "${SUITESPARSE_LIB_DIR}" )
    # not part of the scls graph; needs blas (1) and metis (3)
    belfem_link_libraries( 4 ${SUITESPARSE_LIBRARIES} )
endif()
