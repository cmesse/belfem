# Remove the canonical system library directories from a list of rpath entries.
#
# /usr/lib64, /usr/lib, /lib64 and /lib are searched by the linker and the loader
# by default. Recording one of them in -rpath adds nothing at run time, but at
# link time GNU ld resolves a shared library's own DT_NEEDED entries through the
# -rpath directories first, in order, ahead of LD_LIBRARY_PATH, ahead of the
# library's own RUNPATH and ahead of the default directories. A system copy of a
# library the toolchain prefix also ships ( PMIx for Open MPI is the observed case )
# then shadows the prefix copy on every executable link.
#
# On Debian-family hosts the multiarch directory ( lib/x86_64-linux-gnu ) is a
# system libdir too, and the package finder searches it.
#
# Entries are compared by real path on both sides, because /lib64 is a symlink
# into /usr on most distributions, and trailing slashes are stripped. Other
# entries are kept as written and in order.
function( belfem_prune_system_libdirs aListVar )
    set( tSystemLibdirs /usr/lib64 /usr/lib /lib64 /lib )
    if( CMAKE_LIBRARY_ARCHITECTURE )
        list( APPEND tSystemLibdirs
              "/usr/lib/${CMAKE_LIBRARY_ARCHITECTURE}"
              "/lib/${CMAKE_LIBRARY_ARCHITECTURE}" )
    endif()
    set( tSystemReal )
    foreach( tDir ${tSystemLibdirs} )
        get_filename_component( tDir "${tDir}" REALPATH )
        list( APPEND tSystemReal "${tDir}" )
    endforeach()
    set( tPruned )
    foreach( tItem ${${aListVar}} )
        string( REGEX REPLACE "/+$" "" tItem "${tItem}" )
        get_filename_component( tReal "${tItem}" REALPATH )
        if( NOT tReal IN_LIST tSystemReal )
            list( APPEND tPruned "${tItem}" )
        endif()
    endforeach()
    set( ${aListVar} "${tPruned}" PARENT_SCOPE )
endfunction()
