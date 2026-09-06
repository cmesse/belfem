# signing is needed for macOS
set( BELFEM_SIGNKEY "Christian Messe" )


set( LIBPREFIX "belfem" )

# CMake target of the project library ( Add_BelfemLibrary.cmake ). Distinct
# from LIBPREFIX so the plain name `belfem` stays free for the solver
# executable; the archive on disk is still lib${LIBPREFIX}.a
set( LIBTARGET "belfem_lib" )

# executables deployed by `make install` ( Add_Executable.cmake ). gas only
# exists under USE_GASMODELS and msh2exo only under USE_EXODUS; a name with no
# target is ignored.
#
# banner is deliberately absent: its whole body is print_banner(), which
# `belfem --version` now does. The target is kept — it is the cheapest
# whole-library link check in the tree ( `make banner` ) — but it is not
# something an installation needs to carry.
set( BELFEM_INSTALL_EXECUTABLES belfem material gas db2exo msh2exo )
set( LIBDIR "lib" )
set( BINDIR "bin" )

# basic libraries
set( BELFEM_LIBLIST_BASE
        core
        containers
        comm
        io
        graph
        sparse
        spline )
