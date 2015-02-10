% system 
% 	N(u) * u - L' * p + 1/Re * D * u = 0
% 							   L * u = 0

cd('../') % go up one level
addToPath();


% clear all
% close all
% clc


lidVel = 1; % velocity of moving lid
announce(sprintf(['Finite element implementation of the steady Navier Stokes equations',... 
	' for the lid-driven cavity problem on a [0,1]^2 box with a lid velocity of %d.',...
	' Please enter parameters.'], lidVel), 1, 1)

% create local matrices
[ localMatrix, basisOrder ] = createRectBasis();
nrPBasisF = size(localMatrix.pdivv.x,1);

% create mesh 
feMesh = createRectMesh(basisOrder);

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
% M = vmassAssembly( feMesh, localMatrix.vmass); % mass matrix
globalMatrix.L = PdivVAssembly( feMesh, localMatrix.pdivv); % "pressure mass" matrix
% D = diffusionAssembly( feMesh, localMatrix.stiff); % diffusive matrix
globalMatrix.D = laplaceAssembly( feMesh, localMatrix.stiff); % alternative (not using Sij)
time.assembly = toc;

% determine free nodes (interior)
% lidNodes = unique(feMesh.boundary.gamma2(:)); % requires regularisation
lidNodes = feMesh.boundary.gamma2(:, 2:end-1);
lidNodes = unique(lidNodes(:));

% homogeneous Dirichlet Nodes
homDNodes = unique([feMesh.boundary.gamma1(:);feMesh.boundary.gamma3(:);...
	feMesh.boundary.gamma4(:)]);

fixedVel = [lidNodes; lidNodes + nrNodes; homDNodes; homDNodes + nrNodes;];
freeVel = setdiff(1:2*nrNodes, fixedVel);

fixedPressure = [1];
freePressure = setdiff(1:nrPBasisF*nrElts, fixedPressure);

freeSol = [freeVel, 2*nrNodes + freePressure]; % include pressure DOF
fixedSol = [fixedVel', 2*nrNodes + fixedPressure];

nodeType = struct('freeVel', freeVel, 'freePressure', freePressure,...
	'freeSol', freeSol, 'fixedVel', fixedVel, 'fixedPressure', fixedPressure);

% fill in the boundary conditions
solVec = zeros(2*nrNodes + nrPBasisF*nrElts,1);
% solVec(lidNodes) = lidVel*(1 - feMesh.node(1, lidNodes).^4); % regularised
solVec(lidNodes) = lidVel;

tic;
% solve corresponding stokes problem as initial guess
M = [globalMatrix.D, -globalMatrix.L'; -globalMatrix.L -feMesh.stabC];
rhsVec = -M(:, fixedVel)*solVec(fixedVel);
solVec(freeSol) = matrixSolve(M(freeSol, freeSol), rhsVec(freeSol),1);
time.stokes = toc;

%  clear some variables
clear -regexp ^free ^fixed homDNodes lidNodes

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