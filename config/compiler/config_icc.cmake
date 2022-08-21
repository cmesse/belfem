# set the def telling that we use the gcc

list( APPEND BELFEM_DEFS "BELFEM_INTEL" )

#-------------------------------------------------------------------------------
# Linking and Include Flags
#-------------------------------------------------------------------------------

# use the intel archive program
set( CMAKE_AR xiar )

#-------------------------------------------------------------------------------
# Compiler Optimization Flags
#-------------------------------------------------------------------------------

# test if debug flags are used
if( USE_DEBUG )
    list( APPEND BELFEM_DEFS "DEBUG" )
    set( BELFEM_CXXFLAGS "-O0 -g -debug -Wno-deprecated" )
    set( BELFEM_CXXFLAGS ${BELFEM_CFLAGS} )
    set( BELFEM_FCFLAGS  "-O0 -g -warn all -check all -fpe0 -traceback -debug extended -fstack-security-check -fbounds-check" )
else()
    list( APPEND BELFEM_DEFS "NDEBUG" )

    set( BELFEM_CFLAGS "-O1 -xHost" )
    set( BELFEM_CXXFLAGS ${BELFEM_CFLAGS} )
    set( BELFEM_FCFLAGS  "-O1" )
endif()


#-------------------------------------------------------------------------------
# Warnings
#-------------------------------------------------------------------------------
if( USE_WARNINGS )
    # Add some strict compiler checks.
    set( BELFEM_CXXFLAGS "${BELFEM_CXXFLAGS} -Wall -Werror -Wno-uninitialized -wd161 -wd10014 -wd11012" )
endif()

#-------------------------------------------------------------------------------
# Special Stuff
#-------------------------------------------------------------------------------

# Build 64-bit binaries.
set( BELFEM_CXXFLAGS "${BELFEM_CXXFLAGS} -m64 -std=gnu++14")

# Add UTF-8 support for identifier names.
#set( BELFEM_CXXFLAGS "${BELFEM_CXXFLAGS} -fextended-identifiers")


# Add support to use shared libraries as input files.
# -rdynamic Pass the flag -export-dynamic to the ELF linker,
#           on targets that support it. This instructs the linker
#           to add all symbols, not only used ones,
#           to the dynamic symbol table. This option is needed for
#           some uses of dlopen or to allow obtaining backtraces
#           from within a program.
set(BELFEM_CXXFLAGS "${BELFEM_CXXFLAGS} -rdynamic")
set(BELFEM_CXXFLAGS "${BELFEM_CXXFLAGS} -fPIC")

#set( CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS}" )

#-------------------------------------------------------------------------------
# Fortran Specific Stuff
#-------------------------------------------------------------------------------

# preprocessor for fortran
set( BELFEM_FCFLAGS   "${BELFEM_FCFLAGS} -fpp" )

# fortran library for intel
set( BELFEM_FORTRANLIBS "-lbfp754 -lchkpwrap_h -lchkpwrap_h_w -ldecimal -lifcoremt -lifport -limf -liomp5 -liompstubs5 -lipgo -lirc -lirc_s -lirng -listrconv" )
#set( BELFEM_FORTRANLIBS "-lifcore" )
set( CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -L${INTELROOT}/lib -Wl,-rpath,${INTELROOT}/lib" )

#set( BELFEM_FORTRANLIBS "-lbfp754 -lchkpwrap_h -lchkpwrap_h_w -ldecimal -lifcore -lifcoremt -lifport -limf -liomp5 -liompstubs5 -lipgo -lirc -lirc_s -lirng -listrconv -lsvml" )

if( USE_OPENMP )
    set( BELFEM_CFLAGS "${BELFEM_CFLAGS} -qopenmp" )
    set( BELFEM_CXXFLAGS "${BELFEM_CXXFLAGS} -qopenmp" )
    set( BELFEM_FCFLAGS "${BELFEM_FCFLAGS} -qopenmp" )
    set( BELFEM_OPENMPLIBS "-liomp5" )
    list( APPEND BELFEM_DEFS "OMP" )
else()
    set( BELFEM_OPENMPLIBS "" )
endif()
