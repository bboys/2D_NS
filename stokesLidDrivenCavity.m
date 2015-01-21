% create local matrices
basisOrder = 'Q2P1';
Re = 100; % reynolds number
lidVel = 1;
[ localMatrix ] = createRectBasis(basisOrder );

nrPBasisF = size(localMatrix.pressure.x,1);


% create mesh
nX = 2^5; % nr of "original nodes" in each direction 
nY = 2^5; 

x = linspace(-1,1,nX); y = linspace(-1,1,nY);
feMesh = createRectMesh(x, y, basisOrder);
nrNodes = feMesh.problemSize(3)*feMesh.problemSize(4);
nrElts = feMesh.problemSize(1)*feMesh.problemSize(2);

fprintf('%7.0d elements and %7.0d nodes \n',[nrElts,nrNodes])


% assemble matrices
M = massAssembly( feMesh, localMatrix.mass); % mass matrix
D = diffusionAssembly( feMesh, localMatrix.stiff); % diffusive matrix
L = massPAssembly( feMesh, localMatrix.pressure); % "pressure mass" matrix


% relabel to improve structure
relabel = 1:nrNodes; relabel = [relabel; relabel + nrNodes]; 
relabel = relabel(:);

% solve [D,  -L'; L, 0] [u;p] = 0, under bdy conditions


% determine free nodes (interior)
lidNodes = unique(feMesh.boundary.gamma2(:)); % top
% homogeneous Dirichlet Nodes
homDNodes = unique([feMesh.boundary.gamma1(:);feMesh.boundary.gamma3(:);...
	feMesh.boundary.gamma4(:)]);


fixedVel = [lidNodes; lidNodes + nrNodes; homDNodes; homDNodes + nrNodes;];
freeVel = setdiff(1:2*nrNodes, fixedVel);

% pressure is unique up to constant
freePressure = 2:nrPBasisF*nrElts;

freeSol = [freeVel, 2*nrNodes + freePressure]; % include pressure DOF

% fill in the boundary conditions

solVec = zeros(2*nrNodes + nrPBasisF*nrElts,1);
solVec(lidNodes) = lidVel*(1 - feMesh.node(1, lidNodes).^4);

Dfree = D(freeVel, freeVel); Dfixed = D(freeVel, fixedVel);
Lfree = L(freePressure, freeVel); Lfixed = L(freePressure, fixedVel);

Mfree = [Dfree, -Lfree';-Lfree, sparse(nrPBasisF*nrElts - 1,...
	nrPBasisF*nrElts - 1)];

rhsVec = [-Dfixed*solVec(fixedVel); Lfixed*solVec(fixedVel)];

tic
solVec(freeSol) = Mfree\rhsVec;
toc

plotSol(feMesh, solVec, ['Re = ',num2str(Re)]);

% plot centerline solution
figure
stokesArray = reshape(solVec(1:nrNodes),feMesh.problemSize(4),feMesh.problemSize(3));
centerStokes = stokesArray(:,round(feMesh.problemSize(3)/2));
plot(linspace(-1,1,length(centerStokes)),centerStokes/lidVel,'k')
title(['Centerline velocity, Re = ',num2str(Re)])



