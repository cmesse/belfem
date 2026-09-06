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

    # get native flag. GCC's native CPU detection on Apple arm64 is not
    # reliable, so the tuning flag is only set on x86
    if( NOT CMAKE_SYSTEM_PROCESSOR MATCHES "^(x86_64|AMD64|i[3-6]86)$" )
        set( BELFEM_FC_NATIVE "" )
    elseif( APPLE )
        set( BELFEM_FC_NATIVE "-mtune=native")
    else()
        set( BELFEM_FC_NATIVE "-march=native")
    endif()

    set( BELFEM_CFLAGS "-O2" )
    set( BELFEM_CXXFLAGS ${BELFEM_CFLAGS} )
    set( BELFEM_FCFLAGS  "-O2 ${BELFEM_FC_NATIVE} -fallow-argument-mismatch" )
endif()

if( USE_OPENMP )
    set( BELFEM_CFLAGS "${BELFEM_CFLAGS} -fopenmp" )
    set( BELFEM_CXXFLAGS "${BELFEM_CXXFLAGS} -fopenmp" )
    set( BELFEM_FCFLAGS "${BELFEM_FCFLAGS} -fopenmp" )
    # -fopenmp is already in BELFEM_CFLAGS / CXXFLAGS / FCFLAGS above, and CMake passes
    # those to the link step, where the GCC driver expands it to -lgomp. Naming the
    # library again here duplicates it; see the note on BELFEM_FORTRANLIBS below.
    set( BELFEM_OPENMPLIBS "" )
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

# Build 64-bit binaries. -m64 is an x86 option; aarch64 GCC rejects it.
if( CMAKE_SYSTEM_PROCESSOR MATCHES "^(x86_64|AMD64|i[3-6]86)$" )
    set( BELFEM_CXXFLAGS "${BELFEM_CXXFLAGS} -m64")
endif()

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
    set( CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,-export_dynamic" )
    # @todo temporary hack
    set( BELFEM_COMPILER "Apple Clang" )
    set(BELFEM_CXXFLAGS "${BELFEM_CXXFLAGS}  -std=17 -stdlib=libc++")
else()
    set(BELFEM_CXXFLAGS "${BELFEM_CXXFLAGS} -rdynamic")
endif()

# -fPIC is not set here: CMAKE_POSITION_INDEPENDENT_CODE in the top-level
# CMakeLists.txt covers C, C++ and Fortran uniformly, which a per-language
# flag does not (the Fortran objects in libbelfem_sparse were not PIC).

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
elseif( "gfortran" IN_LIST CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES )
    # CMake already appends the Fortran runtime for us, because enable_language(Fortran)
    # detected it: CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES carries gfortran and quadmath and
    # the generator emits them at the end of the link line. Adding -lgfortran here as well
    # puts it on twice, which Apple's linker reports as
    #     ld: warning: ignoring duplicate libraries: '-lgfortran'
    set( BELFEM_FORTRANLIBS "" )
else()
    set( BELFEM_FORTRANLIBS "-lgfortran" )
endif()
