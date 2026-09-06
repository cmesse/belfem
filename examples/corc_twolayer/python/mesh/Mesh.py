
import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import reverse_cuthill_mckee

from mesh.Node import Node
from mesh.Element import Element
from mesh.AsciiFile import AsciiFile
from mesh.Physical import Physical

#===============================================================================
## An abstract element class
class Mesh:
#-------------------------------------------------------------------------------
    # number of dimensions per GMSH element type
    DimensionPerElement = [-1, 1, 2, 2, 3, 3, 3, -1, 1, 2, 2, 3, 3, 3, -1, 0, 2, 3, 3, -1, -1, 2, -1, -1, -1, -1, 1, -1,
                       -1, 3]

    def __init__(self):

        # number of used nodes
        self.numNodes = 0

        # List of nodes
        self.nodes = []

        # List of elements
        self.elements = []

        # Number of Elements
        self.numElements = 0

        # scaling of this mesh
        self.Scale = 1

        # Number of Geometry groups
        self.NumGroups = 0

        self.PhysicalGroups = []

        # List of PeriodicLink — preserved through read/compact/write
        self.periodic = []

# -------------------------------------------------------------------------------

    def flag_all_nodes( self ):
        for N in self.nodes:
            N.flag = True

# -------------------------------------------------------------------------------

    def unflag_all_nodes(self):
        for N in self.nodes:
            N.flag = False

# -------------------------------------------------------------------------------

    ## create consecutive numbers for the nodes
    def _update_node_ids(self):

        # unflag all nodes
        self.unflag_all_nodes()

        # loop over all elements
        for E in self.elements:
            # loop over all nodes
            for N in E.nodes :
                N.flag = True

        # ID counter
        self.numNodes = 0

        for N in self.nodes :
            if N.flag :
                # increment the counter
                self.numNodes += 1

                # Store the ID in the node
                N.id = self.numNodes
# -------------------------------------------------------------------------------

    def add_physical(self, dimension: int, phystag: int, label: str ):
        self.PhysicalGroups.append( Physical( dimension, phystag, label ) )

# -------------------------------------------------------------------------------

    ## create consecutive numbers for the elements
    def _update_element_ids(self):

        # ID counter
        Count = 1
        for E in self.elements :
            # Store the ID in the node
            E.id = Count

            # increment the counter
            Count += 1

# -------------------------------------------------------------------------------

    ## Print the header
    def _save_header(self, Buffer ):

        # write version
        Buffer.append('$MeshFormat')
        Buffer.append('2.2 0 8')
        Buffer.append('$EndMeshFormat')
# -------------------------------------------------------------------------------

    ## Print the physical groups if they exist
    def _save_physical( self, Buffer ):
        if len( self.PhysicalGroups ) > 0 :
            Buffer.append('$PhysicalNames')
            Buffer.append('{:d}'.format( len( self.PhysicalGroups ) ) )
            for p in self.PhysicalGroups :
                Buffer.append( p.gmsh() )
            Buffer.append('$EndPhysicalNames')
# -------------------------------------------------------------------------------

    ## Print the nodes
    def _save_nodes(self, Buffer ):
        # Update Node IDs
        self._update_node_ids()

        Buffer.append('$Nodes')
        Buffer.append(str( self.numNodes ) )

        for N in self.nodes:
            # test if node is used
            if N.flag :
                # Write line into Buffer
                Buffer.append('{:d} {:21.12f} {:21.12f} {:21.12f}'.format(
                    N.id,
                    N.x * self.Scale,
                    N.y * self.Scale,
                    N.z * self.Scale ) )
        Buffer.append('$EndNodes')

# -------------------------------------------------------------------------------

    def remove_unused_nodes( self ):
        # copy nodes in temporary container
        TempNodes = [None] * ( len( self.nodes ) )

        # flag nodes that are connected
        self._update_node_ids()

        AllCount = 0
        UsedCount = 0

        # copy nodes into temporary array
        for Node in self.nodes:
            TempNodes[AllCount] = Node
            AllCount += 1
            if Node.flag :
                UsedCount += 1

        # copy used nodes back into mesh
        self.nodes.clear()
        self.nodes = [ None ] * UsedCount

        UsedCount = 0
        for Node in TempNodes :
            if Node.flag :
                self.nodes[ UsedCount ] = Node
                UsedCount += 1

        self.numNodes = UsedCount

# -------------------------------------------------------------------------------

    ## Print the elements
    def _save_elements(self, Buffer ):
        # Update Node IDs
        self._update_element_ids()

        Buffer.append('$Elements')
        Buffer.append(str( len( self.elements )))

        for E in self.elements:
            # Write line into Buffer
            Buffer.append( E.gmsh() )

        Buffer.append('$EndElements')

# -------------------------------------------------------------------------------

    ## Save the file to a gmsh mesh
    def save(self, Path: str):

        print('creating ' + Path)

        # Update Element IDs
        self._update_element_ids()

        # create a new and empty Ascii file
        F = AsciiFile()

        # save the header
        self._save_header( F.Buffer )

        # write the physical groups if they exist
        self._save_physical( F.Buffer )

        # save the nodes
        self._save_nodes( F.Buffer )

        # save the elements
        self._save_elements( F.Buffer )

        # save the periodic links if they exist
        self._save_periodic( F.Buffer )

        # tidy up buffer
        F.clean_buffer()

        # save file to disk
        F.save( Path )

# -------------------------------------------------------------------------------

    ## Print the periodic links
    def _save_periodic(self, Buffer ):
        if not self.periodic :
            return

        Buffer.append('$Periodic')
        Buffer.append('{:d}'.format( len( self.periodic ) ) )
        for link in self.periodic :
            Buffer.append('{:d} {:d} {:d}'.format(
                link.dim, link.slave_tag, link.master_tag ) )
            if link.affine :
                Buffer.append( link.affine )
            Buffer.append('{:d}'.format( len( link.node_pairs ) ) )
            for s, m in link.node_pairs :
                Buffer.append('{:d} {:d}'.format( s, m ) )
        Buffer.append('$EndPeriodic')

# -------------------------------------------------------------------------------

    ## optimize bandwidth for better sovler performance
    def optimize_adjacency_bandwidth( self ):

        # number of nodes on mesh
        numNodes = len(self.nodes)

        # set node indices
        Count = 0
        for N in self.nodes:
            N.index = Count
            Count = Count + 1

        # set element indices
        Count = 0
        for E in self.elements:
            E.index = Count
            Count = Count + 1

        # counts how many elements are connected with each node
        elemCount = np.zeros(numNodes, dtype=int)

        # loop over all elements
        for E in self.elements:
            # loop over all nodes of this element
            for N in E.Nodes:
                elemCount[N.index] = elemCount[N.index] + 1

        # Allocate Element Counter
        maxElemCount = max( elemCount )

        elemsPerNode = np.zeros( ( numNodes, maxElemCount ), dtype=int )

        # reset counter
        elemCount = np.zeros(numNodes, dtype=int)

        # loop over all elements
        for E in self.elements:
            # loop over all nodes of this element
            for N in E.Nodes:
                # get element counter
                e = elemCount[N.index]

                # add element to list
                elemsPerNode[N.index][e] = E.index

                # increment counter
                elemCount[N.index] = e + 1

        # collect nodes connected with nodes
        nodeCount = np.zeros((numNodes, 1), dtype=int)

        Adjency = numNodes * [None]

        # loop over all nodes
        for n in range(numNodes):
            # reset counter
            Count = 0

            # loop over all elements connected to this node
            for e in range(elemCount[n]):
                # get element
                E = self.elements[elemsPerNode[n][e]]

                # add number of nodes to counter
                Count = Count + len(E.Nodes)

            # create vector
            NodeIndices = np.zeros(Count, dtype=int)

            # reset counter
            Count = 0
            for e in range(elemCount[n]):
                # get element
                E = self.elements[elemsPerNode[n][e]]

                # loop over all nodes of element
                for N in E.Nodes:
                    # add node index to list
                    NodeIndices[Count] = N.index

                    # increment counter
                    Count = Count + 1

            # add row to adjency
            Adjency[n] = np.unique(NodeIndices)

        # pointers array
        pointers = np.zeros(numNodes + 1, dtype=int)

        # reset counter
        Count = 0

        for n in range(numNodes):
            Count = Count + len(Adjency[n])
            pointers[n + 1] = Count

        # indices array
        indices = np.zeros(Count, dtype=int)

        # reset counter
        Count = 0

        # loop over all nodes
        for i in range(numNodes):
            # loop over all connected nodes
            for k in Adjency[i]:
                indices[Count] = k
                Count += 1

        # create csr matrix
        Graph = csr_matrix((np.ones(Count), indices, pointers), shape=(numNodes, numNodes))

        NewOrder = reverse_cuthill_mckee(Graph)

        # copy nodes into temporary array
        Nodes = numNodes * [None]

        Count = 0
        for N in self.nodes:
            Nodes[Count] = N
            Count += 1

        self.nodes = numNodes * [None]

        for k in range(numNodes):
            self.nodes[k] = Nodes[NewOrder[k]]

        Count = 0

        for N in self.nodes:
            N.index = Count
            Count += 1
            N.id = Count

    # -------------------------------------------------------------------------------

    def collect_elements(self, dimension: int, geotag: int ):
        # counter
        c = 0

        # loop over all elements and count
        for e in self.elements :
            if self.DimensionPerElement[ e.GMSH ] == dimension and e.geometryTag == geotag :
                c+=1

        # allocate container
        E = [ None ] * c

        c = 0
        for e in self.elements :
            if self.DimensionPerElement[ e.GMSH ] == dimension and e.geometryTag == geotag :
                E[ c ] = e
                c+=1
        return E