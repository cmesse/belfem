set( LIBNAME graph )

set( SOURCES
        cl_Graph_Vertex.cpp
        fn_Graph_BFS.cpp
        fn_Graph_find_pseudo_peripheral_node.cpp
        fn_Graph_symrcm.cpp
        fn_Graph_clear.cpp
        fn_Graph_sort.cpp
        )

include( ${BELFEM_CONFIG_DIR}/scripts/Add_Library.cmake )