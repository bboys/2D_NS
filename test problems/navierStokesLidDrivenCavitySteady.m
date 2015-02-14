% system 
% 	N(u) * u - L' * p + 1/Re * D * u = 0
% 							   L * u = 0


cd('../') % go up one level
addToPath();
% initialize setup
setup = [];

lidVel = 1; % velocity of moving lid
announce(sprintf(['Finite element implementation of the steady Navier Stokes equations',... 
	' for the lid-driven cavity problem on a [0,1]^2 box with a lid velocity of %d.',...
	' Please enter parameters.'], lidVel), 1, 1)

% create local matrices
localMatrix = createBasis(2);
nrPBasisF = size(localMatrix.pdivv.x,1);

% create mesh 
feMesh = createMesh(localMatrix.basisOrder, 16, 16, 2);
globalMatrix.stabC = createStabC(feMesh, localMatrix);

Re = 500; % default('Reynolds number', 500); % reynolds number

nrElts = feMesh.problemSize(1)*feMesh.problemSize(2);
nrNodes = feMesh.problemSize(3)*feMesh.problemSize(4);

announce(['Problemsize: ', sprintf('%d elements and %d nodes',...
	[nrElts,nrNodes]),'. Assembling global matrices...'], 0, 1)

% parameters
[setup] = setupLinSolve(setup);
announce('',1,0)

[setup] = setupNonLin(setup);

if setup.linsolve.precon == 1 | setup.nonlin.precon == 1
	announce('',1,0)
	[setup] = setupAMG(setup);
end



tic;
% assemble matrices
globalMatrix.L = PdivVAssembly( feMesh, localMatrix.pdivv); % "pressure mass" matrix
% D = diffusionAssembly( feMesh, localMatrix.stiff); % diffusive matrix
globalMatrix.D = laplaceAssembly( feMesh, localMatrix.stiff); % alternative (not using Sij)
globalMatrix.Q = pmassAssembly(feMesh, localMatrix.pmass);
time.assembly = toc;

% define bdys
feMesh.boundary(1).type = 1; feMesh.boundary(1).func = [0; 0];
feMesh.boundary(3).type = 1; feMesh.boundary(3).func = [0; 0];
feMesh.boundary(4).type = 1; feMesh.boundary(4).func = [0; 0];
feMesh.boundary(2).type = 1;
feMesh.boundary(2).func = str2func('@(x,y) cavityLidDirichlet(x,y)');

% apply bdy conditions
[nodeType, solVec] = applyBdyCond(feMesh, localMatrix.basisType);

tic;
% solve corresponding stokes problem as initial guess
M = [globalMatrix.D, -globalMatrix.L'; -globalMatrix.L -globalMatrix.stabC];
rhsVec = -M(nodeType.freeSol, nodeType.fixedVel)*solVec(nodeType.fixedVel);
[solVec(nodeType.freeSol), relres, resminres] =...
	matrixSolve(M(nodeType.freeVel, nodeType.freeVel),...
	globalMatrix.L(nodeType.freePressure, nodeType.freeVel),...
	globalMatrix.stabC(nodeType.freePressure, nodeType.freePressure),...
	rhsVec, setup,...
	globalMatrix.Q(nodeType.freePressure, nodeType.freePressure), 1);
time.stokes = toc;

announce(sprintf(['Corresponding Stokes problem solved in %4.2f seconds,',...
 ' with a relative residual of %4.2e in %d MVs. ',...
 'Starting nonlinear solve...'],...
	time.stokes, relres, length(resminres)), 0, 1)

stokesSol = solVec;

tic;
% solve nonlinear equation
solVec = nonLinSolve(feMesh, globalMatrix, localMatrix, solVec, nodeType,...
	Re, setup);
time.nonlin = toc;

announce(sprintf(['Nonlinear problem solved in %4.2f seconds. '],...
	time.nonlin), 1, 1, [0.5 0 0])

plotYN = default({'Plot centerline u velocity and velocity magnitude',...
	'Yes', 'No'}, 1);
announce('', 0, 1)

if plotYN == 1
	figure
	cavityCenterPlot(solVec, feMesh, stokesSol, Re, lidVel)
	figure
	plotVel(feMesh, solVec, ['Navier-Stokes, Re = ', num2str(Re)],...
		'interp')
end