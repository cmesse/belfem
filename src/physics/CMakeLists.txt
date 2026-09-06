add_subdirectory( database )
add_subdirectory( materials )

# Gas tables and gas models. Equilibrium chemistry, equations of state and
# transport properties: general thermophysics built on a published method
# (NASA RP-1311, unlimited distribution). Finite-rate combustion is a separate
# module and stays in the nonfree tree.
#
# The reference tables live in share/fluid and are found at runtime through
# gBelfemDataPath, set from $BELFEM_DATA; see fn_GT_data_path.
if( USE_GASMODELS )
    add_subdirectory( gastables )
    add_subdirectory( gasmodels )
endif()