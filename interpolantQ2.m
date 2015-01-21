% function interpVal = interpolantQ2(feMesh, gridPos)
% outputs the value of the interpolant at positions given by gridPos,
% a 2 x ... array. The interpolant is locally (per element) defined by
% the 9 basis functions as in Q2. gridPos is rectilinear.
nrNodes = feMesh.problemSize(3)*feMesh.problemSize(4);
nrGridPoints = size(gridPos,2);
interpVal = zeros(nrGridPoints,1);

% first consider any points defined by gridEval that coincide with gridpoints
[meshNodes, idxMesh, idxGrid] = intersect(feMesh.node',...
	gridPos', 'rows');

% non-mesh grid idx
nonMeshIdx = setdiff(1:nrGridPoints, idxGrid);
nonMesh = gridPos(:,nonMeshIdx);

% nearest mesh point
distArray = bsxfun(@minus, nonMesh(1,:)', feMesh.node(1,:)).^2 +...
	bsxfun(@minus, nonMesh(2,:)', feMesh.node(2,:)).^2
[minDist, minIdx] = min(distArray,[],2);

% determine in which element it belongs








nodeIndices = [];