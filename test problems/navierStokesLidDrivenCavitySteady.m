% system 
% 	N(u) * u - L' * p + 1/Re * D * u = 0
% 							   L * u = 0

cd('../') % go up one level
addToPath();

lidVel = 1; % velocity of moving lid
announce(sprintf(['Finite element implementation of the steady Navier Stokes equations',... 
	' for the lid-driven cavity problem on a [0,1]^2 box with a lid velocity of %d.',...
	' Please enter parameters.'], lidVel), 1, 1)

% create local matrices
[ localMatrix, basisOrder ] = createBasis();
nrPBasisF = size(localMatrix.pdivv.x,1);

% create mesh 
feMesh = createMesh(basisOrder);
globalMatrix.stabC = createStabC(feMesh, localMatrix);

Re = default('Reynolds number', 1000); % reynolds number

nrElts = feMesh.problemSize(1)*feMesh.problemSize(2);
nrNodes = feMesh.problemSize(3)*feMesh.problemSize(4);

% parameters for newton iteration
tol = default('Tolerance for nonlin iteration', 1e-8); 
maxIt = default('Maximum number of nonlin iterations', 15);

announce(['Problemsize: ', sprintf('%d elements and %d nodes',...
	[nrElts,nrNodes]),'. Assembling global matrices...'], 1, 0)

tic;
% assemble matrices
globalMatrix.L = PdivVAssembly( feMesh, localMatrix.pdivv); % "pressure mass" matrix
% D = diffusionAssembly( feMesh, localMatrix.stiff); % diffusive matrix
globalMatrix.D = laplaceAssembly( feMesh, localMatrix.stiff); % alternative (not using Sij)
time.assembly = toc;

% apply bdy conditions
basisType = 'Crouzeix-Raviart';

% define bdys
feMesh.boundary(1).type = 1;
feMesh.boundary(1).func = [0; 0];
feMesh.boundary(3).type = 1;
feMesh.boundary(3).func = [0; 0];
feMesh.boundary(4).type = 1;
feMesh.boundary(4).func = [0; 0];

feMesh.boundary(2).type = 1;
feMesh.boundary(2).func = str2func('@(x,y) cavityLidDirichlet(x,y)');

% apply bdy conditions
[nodeType, solVec] = applyBdyCond(feMesh, basisType);

tic;
% solve corresponding stokes problem as initial guess
M = [globalMatrix.D, -globalMatrix.L'; -globalMatrix.L -globalMatrix.stabC];
rhsVec = -M(:, nodeType.fixedVel)*solVec(nodeType.fixedVel);
solVec(nodeType.freeSol) = matrixSolve(M(nodeType.freeSol,...
	nodeType.freeSol), rhsVec(nodeType.freeSol),1);
time.stokes = toc;

announce(sprintf(['Corresponding Stokes problem solved in %4.2f seconds. ',...
 'Starting nonlinear solve...'],...
	time.stokes), 0, 1)

stokesSol = solVec;

tic;
% solve nonlinear equation
solVec = nonLinSolve(feMesh, globalMatrix, localMatrix, solVec, nodeType,...
	Re, tol, maxIt);
time.nonlin = toc;

announce(sprintf(['Nonlinear problem solved in %4.2f seconds. '],...
	time.nonlin), 1, 1, [0.5 0 0])

plotYN = default({'Plot centerline u velocity', 'Yes', 'No'}, 1);
announce('', 0, 1)

if plotYN == 1
	cavityCenterPlot(solVec, feMesh, stokesSol, Re, lidVel)
end