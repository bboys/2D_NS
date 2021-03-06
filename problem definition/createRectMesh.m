function [ feMesh ] = createRectMesh(basisOrder, nX, nY, uniformMesh)
% creates a quadrilateralation of nX by nY elements (or 2x2 patches for order = 1)
% note that an element is defined by 9 nodes if we consider Q2, the 
% midpoints are added, also in the case of 'Q1P0' for stabilisation
% assuming Crouzeix-Raviart quadrilaterals (number of pressure basisF)


nrPBasisF = 1/2*(basisOrder)*(basisOrder + 1); % dim(P_k-1)
meshType = 'quad';

if uniformMesh == 1
	xVal = linspace(0, 1, nX + 1);
	yVal = linspace(0, 1, nY + 1);
elseif uniformMesh == 2
	xVal = linspace(0,1,nX/2 + 1).^2; % more nodes at boundary
	xVal = [xVal - 1, 1 - xVal(end-1:-1:1)];

	% make sure box is 1x1 (size scales the effective Re!)
	xVal = 1/2 + xVal/2;

	yVal = linspace(0,1,nY/2 + 1).^2; % more nodes at boundary
	yVal = [yVal - 1, 1 - yVal(end-1:-1:1)];

	% make sure box is 1x1 (size scales the effective Re!)
	yVal = 1/2 + yVal/2;
end
xVal = xVal(:);
yVal = yVal(:);

xDif = xVal(2:end) - xVal(1:end-1);
yDif = yVal(2:end) - yVal(1:end-1);
if basisOrder >= 2
	% interpolate to ensure affine equivalence (add basisOrder - 1 points
	% between each two points)

	temp = bsxfun(@times,[0:basisOrder - 1] , xDif(:)/basisOrder);
	tempEndp = xVal(end);

	xVal = bsxfun(@plus, xVal(1:end-1), temp)';
	xVal = [xVal(:); tempEndp];

	temp = bsxfun(@times,[0:basisOrder - 1] , yDif(:)/basisOrder);
	tempEndp = yVal(end);

	yVal = bsxfun(@plus, yVal(1:end-1), temp)';
	yVal = [yVal(:); tempEndp];
elseif basisOrder == 1
	patchSize = bsxfun(@times, yDif, xDif');

	% interpolate to allow stabilisation over 2x2 element patches 5.3.2 Silvester
	tempX = xVal;
	xVal(1:2:2*nX + 1) = xVal;
	xVal(2:2:2*nX) = (tempX(2:end) + tempX(1:end-1))/2;

	tempY = yVal;
	yVal(1:2:2*nY + 1) = yVal;
	yVal(2:2:2*nY) = (tempY(2:end) + tempY(1:end-1))/2;
	
	xDif = xVal(2:end) - xVal(1:end-1);
	yDif = yVal(2:end) - yVal(1:end-1);

	nX = nX*2;
	nY = nY*2;
end

% create an array that gives the surface area per element
areaList = bsxfun(@times, yDif, xDif');
areaList = areaList(:);

% create an array such that the ith column gives the coordinates of the ith
% node. nodes are numbered from bottom to top, left to right.
[X,Y] = meshgrid(xVal,yVal);
nodeList = [X(:)';Y(:)'];

% size per element
hXList = repmat(xDif(:)',nY,1); hXList = hXList(:);
hYList = repmat(yDif(:),nX,1);

% create array that links elements to nodes, that is, the ith column gives
% the (order+1)^2 node indices which define the ith element. elements are numbered
% similar to the nodes, nodes within elements as well.
elementList = zeros((basisOrder + 1)^2,nX*nY);
tempRange = bsxfun(@plus,([1:nX]' - 1)*nY, [1:nY]);

tempOut = bsxfun(@plus,(basisOrder + basisOrder^2*nY)*([1:nX]' - 1) + 1,...
	basisOrder*[0:nY-1]);
elementList(1,tempRange(:)) = tempOut(:);

elementList([2:basisOrder + 1],:) = bsxfun(@plus, [1:basisOrder]', elementList(1,:));

for i = 2:basisOrder + 1;
	elementList((i-1)*(basisOrder + 1) + 1:i*(basisOrder + 1),:) =...
		elementList((i-2)*(basisOrder + 1) + 1:(i-1)*(basisOrder + 1),:) +...
		basisOrder*nY + 1;
end

% create array such that the ith column gives the node indices of the ith 
% edge.
edge = [1:basisOrder*nY];
edge1 = [edge; edge + 1]; % left
edge = [(basisOrder*nY + 1):(basisOrder*nY + 1):(basisOrder*nY+1)*(basisOrder*nX)];
edge2 = [edge; edge + (basisOrder*nY + 1)]; % top
edge = [(basisOrder*nY + 1)*(basisOrder*nX + 1):-1:((basisOrder*nY + 1)*basisOrder*nX + 2)]; 
edge3 = [edge; edge - 1]; % right
edge = [((basisOrder*nY + 1)*basisOrder*nX + 1):(-(basisOrder*nY + 1)):(basisOrder*nY + 2)];
edge4 = [edge; edge - (basisOrder*nY + 1)]; % bottom

% types of boundary conditions (see page 43, Segal)
% 1 : dirichlet, u = g1,
% 2 : normal velocity u_n = g2, sigma_nt = g3 = 0, 
% 3 : tangential velocity u_t = g4, sigma_nn = g5 = 0, 
% 4 : sigma_nt = 0 = sigma_nn 

boundary(1).edges = [edge1];
boundary(2).edges = [edge2];
boundary(3).edges = [edge3];
boundary(4).edges = [edge4];

% store sizes
problemSize = [nX,nY,(basisOrder*nX+1),(basisOrder*nY+1)]; % nr elements in x,y dir, nr of nodes

feMesh = struct('node',nodeList,'elt',elementList,'boundary',boundary,...
    'area',areaList,'problemSize',problemSize,'eltSize',[hXList';hYList'],...
    'gridX',xVal,'gridY',yVal, 'basisOrder', basisOrder...
    ,'meshType', meshType);

end
