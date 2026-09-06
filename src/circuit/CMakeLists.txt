set( LIBNAME circuit )

set( SOURCES
        cl_ElectricNode.cpp
        cl_Component.cpp
        cl_TwoTerminals.cpp
        cl_FEMTwoTerminals.cpp
        cl_Resistor.cpp
        cl_ElectricalCircuit.cpp
        cl_ElectricalCircuitFactory.cpp
        cl_VoltageSource.cpp
        cl_CurrentSource.cpp
        cl_Switch.cpp
        cl_Inductor.cpp
        cl_Capacitor.cpp
        cl_Diode.cpp
        cl_Superconductor.cpp
        Component_Enums.cpp
        fn_spice_number.cpp
        cl_NetlistParser.cpp
        cl_NgspiceCircuitFactory.cpp
)

# trimmed 2026-08-15 ( Cut A ): the base dirs ( containers, linalg, io,
# core, comm, sparse, math/graph, ... ) are injected by Add_Library.cmake;
# the rows below are the ones circuit TUs actually include directly.
# fem/kernel stays until Cut B relocates the PhysicalBoundaryCondition
# construction out of ElectricalCircuitFactory
include_directories( ${BELFEM_SOURCE_DIR}/numerics/ode )
include_directories( ${BELFEM_SOURCE_DIR}/numerics/sources )
include_directories( ${BELFEM_SOURCE_DIR}/fem/kernel )


include( ${BELFEM_CONFIG_DIR}/scripts/Add_Library.cmake )
