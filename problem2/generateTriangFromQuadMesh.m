function [nodes,elemT]=generateTriangFromQuadMesh(nodes, elemQ)
%
% Generate a triangular mesh from a quadragular one:
% (it is not appropiate for big meshes)
% Input: quadrilateral mesh defined by nodes and elemQ
% Output: triangular mesh defined by same nodes and new elemT
%
numElemQ=size(elemQ);
elemT=[];
for i=1:numElemQ
    elemT=[elemT; elemQ(i,1),elemQ(i,2),elemQ(i,3)];
    elemT=[elemT; elemQ(i,1),elemQ(i,3),elemQ(i,4)];
end
