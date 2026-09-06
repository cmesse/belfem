set( LIBNAME graph )

set( SOURCES
        cl_Graph_Vertex.cpp
        fn_Graph_bfs.cpp
        fn_Graph_dfs.cpp
        fn_Graph_spfa.cpp
        fn_Graph_clear.cpp
        fn_Graph_sort.cpp
        fn_Graph_find_connected_partitions.cpp
        fn_Graph_find_pseudo_peripheral_vertex.cpp
        fn_Graph_find_pseudo_peripheral_node.cpp
        fn_Graph_multibfs.cpp
        fn_Graph_symrcm.cpp
        fn_Graph_METIS.cpp
        fn_Graph_ParMETIS.cpp
        fn_Graph_SCOTCH.cpp
        fn_Graph_PTSCOTCH.cpp
        )

include( ${BELFEM_CONFIG_DIR}/scripts/Add_Library.cmake )
