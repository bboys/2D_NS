basisOrder = 'Q1P0';
Re = 50; % reynolds number

[ localMatrix ] = createRectBasis(basisOrder );
nrPBasisF = size(localMatrix.pressure.x,1);


% create mesh
nX = 2^1; % nr of "original nodes" in each direction 
nY = nX; % 2^2; 

% x = linspace(0,1,nX/2 + 1).^3; % more nodes at boundary
% x = [x - 1, 1 - x(end-1:-1:1)];

x = linspace(0,1,nX);
y = x; 

feMesh = createRectMesh(x, y, basisOrder);
nrNodes = feMesh.problemSize(3)*feMesh.problemSize(4);
nrElts = feMesh.problemSize(1)*feMesh.problemSize(2);

fprintf('%7.0d elements and %7.0d nodes \n',[nrElts,nrNodes])

% create initial velocity field (& solve pressure?)
velocity = zeros(2*nrNodes,1); % should match bdy conditions

% assemble matrices
M = massAssembly( feMesh, localMatrix.mass); % mass matrix
DS = diffusionAssembly( feMesh, localMatrix.stiff); % diffusive matrix
D = laplaceAssembly( feMesh, localMatrix.stiff);
L = massPAssembly( feMesh, localMatrix.pressure); % "pressure mass" matrix
globalMatrix = struct('M',M,'D',D,'L',L);

% determine boundaries
inflowNodes = unique(feMesh.boundary.gamma1(:));
outflowNodes =  unique(feMesh.boundary.gamma3(:));

% homogeneous Dirichlet Nodes
homDNodes = unique([feMesh.boundary.gamma2(:);feMesh.boundary.gamma4(:)]);

fixedVel = unique([inflowNodes; inflowNodes + nrNodes; homDNodes;...
	homDNodes + nrNodes; outflowNodes]);
freeVel = setdiff(1:2*nrNodes, fixedVel);

freeSol = [freeVel, 2*nrNodes + [1:nrPBasisF*nrElts]]; % include pressure DOF
fixedSol = setdiff(1:2*nrNodes + nrPBasisF*nrElts, freeSol);

% test convective matrix with u = 3*x + 1, v = 1 => N1*[u v]' = [3*u 0]'
velTemp = 3*(linspace(0,1,feMesh.problemSize(3))) + 1;

velocity(nrNodes + 1:2*nrNodes) = reshape(repmat(velTemp,feMesh.problemSize(4),1),nrNodes,1);
velocity(1:nrNodes) = 1;

% D(freeVel,1:2*nrNodes)*velocity % (should be zero for laplace)
v = velocity;
calcVorticity(feMesh, v)

% system of ODEs
% 	d/dt M * u + N(u) * u - L' * p + 1/Re * D * u = f
% 										    L * u = 0

test = rand(size(v));
N1 = nonLinearAssembly1(feMesh, localMatrix.nonlin, test );
N1s = nonLinearAssembly1Standard(feMesh, localMatrix.nonlin, test );

t = N1 - N1s;
norm(t(:))
