# Read the SCLS toolchain flavor and derive the build defaults from it.
#
# This MUST stay above the option() calls in CMakeLists.txt. option() sets a
# default only when the cache entry does not exist yet, so a default computed
# after the declaration has no effect - and no visible effect either, because the
# first configure simply caches the un-preset value and every later one reads it
# back. The full SCLS setup ( link directories, rpath, includes ) still runs
# later, in find_scls.cmake; only the environment read happens here.
#
# What this file does NOT do is search the host or check anything. $SCLS_FLAVOR
# is a statement by the user about which toolchain they sourced, and everything
# below is a *default*: an explicit -DUSE_MKL=... / -DUSE_PARDISO=... /
# -DUSE_DEBUG=... on the command line or in ccmake always wins, and a value
# already in the cache is never touched.
#
#   flavor              USE_MKL   USE_PARDISO   USE_DEBUG
#   name contains mkl   ON        ON            OFF
#   debug               OFF       OFF           ON
#   anything else       OFF       OFF           OFF
#   no $SCLS            OFF       OFF           OFF
#
# USE_PARDISO defaults in parity with USE_MKL: PARDISO is part of MKL, so the
# flavor that provides one provides the other. Any vendor other than MKL is
# assumed to behave like reference BLAS/LAPACK, so no other flavor name means
# anything here.

set( BELFEM_SCLS_FLAVOR "" )
set( BELFEM_DEFAULT_USE_MKL     OFF )
set( BELFEM_DEFAULT_USE_PARDISO OFF )
set( BELFEM_DEFAULT_USE_DEBUG   OFF )

if( DEFINED ENV{SCLS} )

    # the flavors live at /opt/scls/<flavor>, so the basename carries the same
    # statement as $SCLS_FLAVOR and keeps a shell that exports only $SCLS working
    get_filename_component( BELFEM_SCLS_BASENAME "$ENV{SCLS}" NAME )

    if( DEFINED ENV{SCLS_FLAVOR} )
        set( BELFEM_SCLS_FLAVOR "$ENV{SCLS_FLAVOR}" )

        # a half-sourced module leaves the two disagreeing. warn, then trust the
        # explicit variable: it is the more specific statement
        if( NOT BELFEM_SCLS_FLAVOR STREQUAL BELFEM_SCLS_BASENAME )
            message( WARNING
                     "SCLS_FLAVOR is '${BELFEM_SCLS_FLAVOR}' but SCLS points at '$ENV{SCLS}' "
                     "( flavor '${BELFEM_SCLS_BASENAME}' ). Using SCLS_FLAVOR. "
                     "Re-source the toolchain if that is not what you want." )
        endif()
    else()
        set( BELFEM_SCLS_FLAVOR "${BELFEM_SCLS_BASENAME}" )
    endif()

    if( BELFEM_SCLS_FLAVOR MATCHES "mkl" )
        set( BELFEM_DEFAULT_USE_MKL     ON )
        set( BELFEM_DEFAULT_USE_PARDISO ON )
    elseif( BELFEM_SCLS_FLAVOR STREQUAL "debug" )
        set( BELFEM_DEFAULT_USE_DEBUG ON )
    endif()

endif()
