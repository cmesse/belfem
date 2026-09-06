clear all;
clc;
close all;

mesh = gmsh.load( "penta18.msh" );



nSides = 5;
nEdges = 9;

tHEX = [0, 1, 2, 3, 4, 5, 6, 7, 8, 11, 12, 9, 13, 10, 14, 15, 16, 19, 17, 18, 21, 25, 23, 24, 26, 22, 20];
tPENTA= [1, 2, 3, 4, 5, 6, 7, 9, 10, 8, 11, 12, 13, 15, 14, 16, 18, 17];

tFile = fopen( 'outfile.txt', 'w');
fprintf( tFile, '// FACES\n' );
for s = 1:nSides
    fprintf( tFile, '                case( %i ):\n', s-1 );
    NODES = mesh.ELEMENTS( s ).NODES;
    n = length( NODES );
    for k = 1:n
        fprintf( tFile, '                    aNodes( %i ) = mNodes[ %i ];\n', k-1, tPENTA(NODES(k))-1);
    end
    fprintf( tFile, '                    break;\n');
end


fprintf( tFile, '// EDGES\n' );
for e = 1:nEdges
    fprintf( tFile, '                case( %i ):\n', e-1 );
    NODES = mesh.ELEMENTS( e+nSides ).NODES;
    n = length( NODES );
    for k = 1:n
        fprintf( tFile, '                    aNodes( %i ) = mNodes[ %i ];\n', k-1, tPENTA(NODES(k))-1);
    end
    fprintf( tFile, '                    break;\n');
end
fclose( tFile );
%ELEMENT_TYPE = mesh.ELEMENTS( mesh.number_of_elements ).type;
%ELEMENT = mesh.ELEMENTS( mesh.number_of_elements ).NODES';

%nNodes = mesh.ELEMENTS( mesh.number_of_elements ).number_of_nodes;

%NODES = [ (1:nNodes)' mesh.NODE_COORDS( ELEMENT, : ) ]

