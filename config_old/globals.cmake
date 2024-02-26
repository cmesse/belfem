# signing is needed for macOS
set( BELFEM_SIGNKEY "Christian Messe" )


set( LIBPREFIX "belfem" )
set( LIBSUFFIX "a" )
set( LIBDIR "lib" )
set( BINDIR "bin" )

# basic libraries
set( BELFEM_LIBLIST_BASE
        core
        comm
        containers
        io
        graph
        sparse
        spline )
