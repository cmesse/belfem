add_subdirectory( core )
add_subdirectory( comm )
add_subdirectory( containers )
add_subdirectory( io )

add_subdirectory( math )
add_subdirectory( sparse )
add_subdirectory( numerics )
add_subdirectory( physics )
add_subdirectory( mesh )
add_subdirectory( fem )
add_subdirectory( homology )


if( USE_MAXWELL )
    add_subdirectory( lbl )
endif()