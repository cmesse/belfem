# set the def telling that we use the gcc
set(CMAKE_CXX_STANDARD 17)

if( BELFEM_USE_CLANG )
    list( APPEND BELFEM_DEFS "BELFEM_CLANG" )
else()
    list( APPEND BELFEM_DEFS "BELFEM_GCC" )
endif()

#-------------------------------------------------------------------------------
# Compiler Optimization flags
#-------------------------------------------------------------------------------

# test if debug flags are used
if( USE_DEBUG )
    list( APPEND BELFEM_DEFS "DEBUG" )
    set( BELFEM_CXXFLAGS "-Og -g" )
    set( BELFEM_CFLAGS ${BELFEM_CXXFLAGS} )
    set( BELFEM_FCFLAGS  "-O0  -fcheck=bounds -fbacktrace -fallow-argument-mismatch" )
else()
    list( APPEND BELFEM_DEFS "NDEBUG" )

    # get native flag
    if( APPLE )
        set( BELFEM_FC_NATIVE "-mtune=native")
    else()
        set( BELFEM_FC_NATIVE "-march=native")
    endif()

    set( BELFEM_CFLAGS "-O2" )
    set( BELFEM_CXXFLAGS ${BELFEM_CFLAGS} )
    set( BELFEM_FCFLAGS  "-O2 ${BELFEM_FC_NATIVE} -fallow-argument-mismatch" )
endif()

if( USE_OPENMP AND NOT USE_MKL )
    set( BELFEM_CFLAGS "${BELFEM_CFLAGS} -fopenmp" )
    set( BELFEM_CXXFLAGS "${BELFEM_CXXFLAGS} -fopenmp" )
    set( BELFEM_FCFLAGS "${BELFEM_FCFLAGS} -fopenmp" )
    set( BELFEM_OPENMPLIBS "-lgomp" )
    list( APPEND BELFEM_DEFS "OMP" )
else()
    set( BELFEM_OPENMPLIBS "" )
endif()

#-------------------------------------------------------------------------------
# Warnings
#-------------------------------------------------------------------------------
if( USE_WARNINGS )
    # Add some strict compiler checks.
    # -pedantic-errors     Make all pedantic warnings into errors. -pedantic issues
    #                      all the warnings demanded by strict ISO C and ISO C++.
    # -Wall                This enables all the warnings about constructions
    #                      that some users consider questionable, and
    #                      that are easy to avoid (or modify to prevent the warning),
    #                      even in conjunction with macros.
    #                      However, it does not enable all warnings available.
    # -Werror              Make all warnings into errors.
    # -Wconversion         Give warning if conversion between data types occurs.
    # -Wno-long-long       Do not issue a warning if a long long variable type is used.
    # -Wno-error=maybe-uninitialized
    #                      Do not issue an error if a maybe-uninitialized warning is thrown
    # -fno-strict-aliasing Do not enforce strict aliasing.
    #                      Strict aliasing means that pointer arguments in a function are assumed to not alias.
    #                      For example, the following code would not compile: foo * a; bar * b; b = (foo *) a;
    #                      Because the pointers point to fundamentally different types.
    if ( BELFEM_USE_CLANG )
        set( BELFEM_CXXFLAGS "${BELFEM_CXXFLAGS} -Wall -Werror=uninitialized" )
    else()
        set( BELFEM_CXXFLAGS "${BELFEM_CXXFLAGS} -Wall -Werror -Wno-long-long -pedantic-errors -Wno-error=maybe-uninitialized" )
    endif()
endif()
#-------------------------------------------------------------------------------
# Special Stuff
#-------------------------------------------------------------------------------

# Build 64-bit binaries.
set( BELFEM_CXXFLAGS "${BELFEM_CXXFLAGS} -m64")

# Add UTF-8 support for identifier names.
set( BELFEM_CXXFLAGS "${BELFEM_CXXFLAGS} -fextended-identifiers")


# Add support to use shared libraries as input files.
# -rdynamic Pass the flag -export-dynamic to the ELF linker,
#           on targets that support it. This instructs the linker
#           to add all symbols, not only used ones,
#           to the dynamic symbol table. This option is needed for
#           some uses of dlopen or to allow obtaining backtraces
#           from within a program.
if ( BELFEM_USE_CLANG )
    set(BELFEM_LIBS "${BELFEM_LIBS} -Wl,-export_dynamic")
    # @todo temporary hack
    set( BELFEM_COMPILER "Apple Clang" )
    set(BELFEM_CXXFLAGS "${BELFEM_CXXFLAGS}  -std=17 -stdlib=libc++")
else()
    set(BELFEM_CXXFLAGS "${BELFEM_CXXFLAGS} -rdynamic")
endif()

set(BELFEM_CXXFLAGS "${BELFEM_CXXFLAGS} -fPIC")

#-------------------------------------------------------------------------------
# Fortran Specific Stuff
#-------------------------------------------------------------------------------

# preprocessor for fortran
set( BELFEM_FCFLAGS   "${BELFEM_FCFLAGS} -cpp" )

# fortran library for gcc
set( BELFEM_FORTRANLIBS "-lgfortran" )

if ( BELFEM_USE_CLANG AND APPLE)
    # @todo temporary hack
    set( BELFEM_FORTRANLIBS "-stdlib=libc++ /opt/gcc/latest/lib/libgfortran.a /opt/gcc/latest/lib/libquadmath.a /opt/gcc/latest/lib/gcc/x86_64-apple-darwin21.5.0/12.1.0/libgcc.a")
    set( CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,-no_compact_unwind" )
else()
    set( BELFEM_FORTRANLIBS "-lgfortran" )
    set( CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,-rpath,$ENV{TPLS}/lib -L$ENV{TPLS}/lib " )
endif()