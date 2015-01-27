% create local matrices
basisOrder = 'Q1P0';
[ localMatrix ] = createRectBasis(basisOrder );

nrPBasisF = size(localMatrix.pressure.x,1);


% create mesh
nX = 2^4; % nr of "original nodes" in each direction 
nY = 2^4; 

x = linspace(-1,1,nX); y = linspace(-1,1,nY);
feMesh = createRectMesh(x, y, basisOrder);
nrNodes = feMesh.problemSize(3)*feMesh.problemSize(4);
nrElts = feMesh.problemSize(1)*feMesh.problemSize(2);

fprintf('%7.0d elements and %7.0d nodes \n',[nrElts,nrNodes])


% assemble matrices
% M = massAssembly( feMesh, localMatrix.mass); % mass matrix
% D = diffusionAssembly( feMesh, localMatrix.stiff); % diffusive matrix
D = laplaceAssembly( feMesh, localMatrix.stiff); % alternative (not using Sij)
L = massPAssembly( feMesh, localMatrix.pressure); % "pressure mass" matrix



% relabel to improve structure
relabel = 1:nrNodes; relabel = [relabel; relabel + nrNodes]; 
relabel = relabel(:);

% solve [D,  -L'; -L, 0] [u;p] = 0, under bdy conditions
% u = [1 - y^2; 0 ] on gamma1 (left, inflow)
% u = 0 on gamma2,4 (top & bottom, no slip wall)
% v_n = 0 on gamma3 & p = 0 (right, outflow)

% determine boundaries
inflowNodes = unique(feMesh.boundary.gamma1(:));
outflowNodes =  unique(feMesh.boundary.gamma3(:));

% homogeneous Dirichlet Nodes
homDNodes = unique([feMesh.boundary.gamma2(:);feMesh.boundary.gamma4(:)]);

fixedVel = unique([inflowNodes; inflowNodes + nrNodes; homDNodes;...
	homDNodes + nrNodes]);
freeVel = setdiff(1:2*nrNodes, fixedVel);

freeSol = [freeVel, 2*nrNodes + [1:nrPBasisF*nrElts]]; % include pressure DOF
fixedSol = setdiff(1:2*nrNodes + nrPBasisF*nrElts, freeSol);

% fill in the boundary conditions
solVec = zeros(2*nrNodes + nrPBasisF*nrElts,1);
solVec(inflowNodes) = 1 - feMesh.node(2, inflowNodes).^2; % gamma1

% determine "free" part of matrices and part that multiplies with "fixed DOF"
Dfree = D(freeVel, freeVel); Dfixed = D(freeVel, fixedVel);
Lfree = L(:, freeVel); Lfixed = L(:, fixedVel);
Afree = [Dfree, -Lfree';Lfree, sparse(nrPBasisF*nrElts,...
	nrPBasisF*nrElts)];


% construct RHS vector containing bdy conditions
rhsVec = [-Dfixed*solVec(fixedVel); -Lfixed*solVec(fixedVel)];

% solve the linear system 
solVec(freeSol) = Afree\rhsVec;

% plot solution
plotSol(feMesh, solVec, ' ');


stokesArray = reshape(solVec(1:nrNodes),feMesh.problemSize(4),feMesh.problemSize(3));
centerStokes = stokesArray(:,round(feMesh.problemSize(3)/2));
yPoints = feMesh.node(2,1:feMesh.problemSize(4));

figure
plot(yPoints, centerStokes, 'kx-')
hold on
plot(yPoints, 1 - yPoints.^2, 'ro')
hold off
