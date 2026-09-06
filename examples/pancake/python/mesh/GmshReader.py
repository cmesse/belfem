import numpy as np

from mesh.AsciiFile import AsciiFile

from mesh.Node import Node
from mesh.Element import Element
from mesh.Vertex import Vertex
from mesh.Line2 import Line2
from mesh.Line3 import Line3

from mesh.Tri3 import Tri3
from mesh.Tri6 import Tri6

from mesh.Quad4 import Quad4
from mesh.Quad8 import Quad8
from mesh.Quad9 import Quad9

from mesh.Tet4 import Tet4
from mesh.Tet10 import Tet10

from mesh.Pyra5 import Pyra5
from mesh.Pyra14 import Pyra14

from mesh.Penta6 import Penta6
from mesh.Penta15 import Penta15
from mesh.Penta18 import Penta18

from mesh.Hex8 import Hex8
from mesh.Hex20 import Hex20
from mesh.Hex27 import Hex27

from mesh.Mesh import Mesh
from mesh.PeriodicLink import PeriodicLink


class GmshReader( AsciiFile ) :

    def __init__(self, path: str ):
        AsciiFile.__init__( self, path )
        self.clean_buffer()

        self.mesh = Mesh()

        self._nodeOffset = 0
        self._elemOffset = 0
        self._periodicOffset = -1

        self._read_flags()
        self._read_nodes()
        self._read_elements()
        self._read_periodic()

    def _read_flags( self ) :
        n = len(self.Buffer)
        l = 0
        for line in self.Buffer :
            if line == "$Nodes" :
                break
            l += 1
        self.mesh.numNodes = int( self.Buffer[l+1] )
        self._nodeOffset = l + 2

        for l in range( self._nodeOffset + self.mesh.numNodes, n ) :
            if self.Buffer[l] == "$Elements" :
                self.mesh.numElements = int( self.Buffer[l+1] )
                self._elemOffset = l + 2
                break

        for l in range( self._elemOffset + self.mesh.numElements, n ) :
            if self.Buffer[l] == "$Periodic" :
                self._periodicOffset = l + 1
                break

    def _read_nodes(self):
        self.mesh.nodes= self.mesh.numNodes * [ None ]

        # the msh2 format allows non-dense, unordered node ids, so elements
        # must look nodes up by id, not by position in the file
        self._nodeMap = {}

        for k in range( self.mesh.numNodes ) :
            data = self.Buffer[ self._nodeOffset + k ].split()

            id = int(data[0])
            x = float(data[1])
            y = float( data[2])
            z = float( data[3])

            n = Node( x, y, z )
            n.id = id
            n.index = k
            self.mesh.nodes[ k ] = n
            self._nodeMap[ id ] = n

    def _read_elements(self):
        self.mesh.elements = self.mesh.numElements * [ None ]

        for e in range( self.mesh.numElements ) :

            # get data from buffer
            sdata = self.Buffer[ self._elemOffset + e ].split()
            idata = np.zeros( len( sdata ), dtype=int )
            for k in range( len(sdata)) :
                idata[k] = int(sdata[k])
    
            E = self._create_element( idata[1] )
            E.id = idata[0]
            E.physicalTag = idata[3]
            E.geometryTag = idata[4]
            n = len( E.nodes )
            off = idata[2] + 3
            for k in range( n ) :
                E.nodes[k] = self._nodeMap[ idata[off+k] ]
            self.mesh.elements[e] = E


    def _read_periodic(self):
        if self._periodicOffset < 0 :
            return

        l = self._periodicOffset
        numLinks = int( self.Buffer[ l ] )
        l += 1

        links = []
        for _ in range( numLinks ) :
            header = self.Buffer[ l ].split()
            dim = int( header[0] )
            slave_tag = int( header[1] )
            master_tag = int( header[2] )
            l += 1

            affine = None
            if self.Buffer[ l ].startswith("Affine") :
                affine = self.Buffer[ l ]
                l += 1

            numPairs = int( self.Buffer[ l ] )
            l += 1

            pairs = [ None ] * numPairs
            for k in range( numPairs ) :
                data = self.Buffer[ l ].split()
                pairs[ k ] = ( int( data[0] ), int( data[1] ) )
                l += 1

            links.append( PeriodicLink( dim, slave_tag, master_tag, affine, pairs ) )

        self.mesh.periodic = links

    def _create_element(self, type: int ):
        if type == 1 :
            return Line2()
        elif type == 2 :
            return Tri3()
        elif type == 3 :
            return Quad4()
        elif type == 4 :
            return Tet4()
        elif type == 5 :
            return Hex8()
        elif type == 6 :
            return Penta6()
        elif type == 7 :
            return Pyra5()
        elif type == 8 :
            return Line3()
        elif type == 9 :
            return Tri6()
        elif type == 10 :
            return Quad9()
        elif type == 11 :
            return Tet10()
        elif type == 12 :
            return Hex27()
        elif type == 13 :
            return Penta18()
        elif type == 14 :
            return Pyra14()
        elif type == 15 :
            return Vertex()
        elif type == 16 :
            return Quad8()
        elif type == 17 :
            return Hex20()
        elif type == 18 :
            return Penta15()
        else:
            raise Exception("Unknown element type")


