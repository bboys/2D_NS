function [feMesh] = createMesh(basisOrder, nX, nY, uniformMesh)
% if only basisOrder is given, ask user for input
if nargin < 3
	if basisOrder == 1
		sizeQuery = 'Number of 2x2 patches in each direction';
	else
		sizeQuery = 'Number of elements in each direction';
	end
	% nr of elements/patches in each direction (for Q1P0 defines the 2x2 patches)
	nX = default(sizeQuery, 2^4);
	nX = ceil(nX/2)*2; % make sure nX is even for nonuniform mesh
	nY = nX; 
end

if nargin < 4
	uniformMesh = default({'Which mesh spacing to use', 'Uniform',...
	 'Non-uniform'}, 2, {'', 'more nodes at the boundary'});
end

% for now only rectangular elements
[ feMesh ] = createRectMesh(basisOrder, nX, nY, uniformMesh);

end