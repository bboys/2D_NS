% system 
% 	N(u) * u - L' * p + 1/Re * D * u = 0
% 							   L * u = 0
cd('../') % go up one level
addToPath();


clear all
close all
clc


lidVel = 1; % velocity of moving lid
announce(sprintf(['Finite element implementation of the steady Navier Stokes equations',... 
	' for the lid-driven cavity problem on a [0,1]^2 box with a lid velocity of %d.',...
	' Please enter parameters.'], lidVel), 1, 1)

% create local matrices
[ localMatrix, basisOrder ] = createRectBasis();
nrPBasisF = size(localMatrix.pressure.x,1);

% create mesh 
feMesh = createRectMesh(basisOrder);

Re = default('Reynolds number', 1000); % reynolds number

nrElts = feMesh.problemSize(1)*feMesh.problemSize(2);
nrNodes = feMesh.problemSize(3)*feMesh.problemSize(4);

% parameters for newton iteration
tol = default('Tolerance for nonlin iteration', 1e-12); 
maxIt = default('Maximum number of nonlin iterations', 15);


announce(['Problemsize: ', sprintf('%d elements and %d nodes',...
	[nrElts,nrNodes]),'. Assembling global matrices...'], 1, 0)



tic;
% assemble matrices
% M = massAssembly( feMesh, localMatrix.mass); % mass matrix
L = massPAssembly( feMesh, localMatrix.pressure); % "pressure mass" matrix
% D = diffusionAssembly( feMesh, localMatrix.stiff); % diffusive matrix
D = laplaceAssembly( feMesh, localMatrix.stiff); % alternative (not using Sij)
globalMatrix = struct('D', D, 'L', L);
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

typesNodes = struct('freeVel', freeVel, 'freePressure', freePressure,...
	'freeSol', freeSol, 'fixedVel', fixedVel, 'fixedPressure', fixedPressure);

% fill in the boundary conditions
solVec = zeros(2*nrNodes + nrPBasisF*nrElts,1);
% solVec(lidNodes) = lidVel*(1 - feMesh.node(1, lidNodes).^4); % regularised
solVec(lidNodes) = lidVel;

tic;
% solve corresponding stokes problem as initial guess
M = [D, -L'; -L -feMesh.stabC];
rhsVec = -M(:, fixedVel)*solVec(fixedVel);
solVec(freeSol) = matrixSolve(M(freeSol, freeSol), rhsVec(freeSol),1);
time.stokes = toc;

announce(sprintf(['Corresponding Stokes problem solved in %4.2f seconds. ',...
 'Starting nonlinear solve...'],...
	time.stokes), 0, 1)

stokesSol = solVec;

tic;
% solve nonlinear equation
solVec = nonLinSolve(feMesh, globalMatrix, localMatrix, solVec, typesNodes,...
	Re, tol, maxIt);
time.nonlin = toc;

announce(sprintf(['Nonlinear problem solved in %4.2f seconds. ',...
 'Plotting results...'],...
	time.nonlin), 1, 1, [0.5 0 0 ])

% plot stuff
% plotData = default();

solArray = reshape(solVec(1:nrNodes),feMesh.problemSize(4),feMesh.problemSize(3));
centerSol = solArray(:,round(feMesh.problemSize(3)/2));

stokesArray = reshape(stokesSol(1:nrNodes),feMesh.problemSize(4),feMesh.problemSize(3));
centerStokes = stokesArray(:,round(feMesh.problemSize(3)/2));

% plotSol(feMesh, solVec, ['Re = ',num2str(Re)]);

yPoints = feMesh.node(2,1:feMesh.problemSize(4));

figure 
plot(yPoints,centerSol/lidVel,'r')
hold on
plot(yPoints,centerStokes/lidVel,'k')

% reference solution
load referenceSols
refSolStokes(:,2) = -refSolStokes(:,2)*2 + 1;

plot(refSolYVal, refSol400, 'rx' )
plot(1/2 + comsolNSRe1000(:,1)/2,  comsolNSRe1000(:,2), 'r+' )

plot(1/2 + comsolCavityNSRe1(:,1)/2, comsolCavityNSRe1(:,2), 'kx')

title(['Centerline velocity, Re = ',num2str(Re)])
legend('NS','Stokes', 'Ghia 400 NS', 'Comsol NS 1000', 'Stokes Comsol')

hold off
% grid on