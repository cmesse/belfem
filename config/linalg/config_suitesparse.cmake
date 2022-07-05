if( USE_SUITESPARSE )
    list( APPEND BELFEM_DEFS "BELFEM_SUITESPARSE" )
    set( BELFEM_MATRIX_LIBS "-lumfpack -lklu -lcholmod -lmetis -lccolamd -lcolamd -lcamd -lbtf -lamd -lsuitesparseconfig ${BELFEM_MATRIX_LIBS}")

    list( APPEND BELFEM_INCLUDES "$ENV{TPLS}/include/suitesparse" )
endif()