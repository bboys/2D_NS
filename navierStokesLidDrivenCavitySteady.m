% system 
% 	N(u) * u - L' * p + 1/Re * D * u = 0
% 							   L * u = 0

% clear all
close all

% create local matrices
basisOrder = 'Q2P1';
Re = 400; % reynolds number
lidVel = 1; % velocity of moving lid

% parameters for newton iteration
tol = 1e-12; 
maxIt = 1;

% reference solution
load referenceSols
refSolStokes(:,2) = -refSolStokes(:,2)*2 + 1;

% initialise
[ localMatrix ] = createRectBasis(basisOrder );
nrPBasisF = size(localMatrix.pressure.x,1);

% create mesh
nX = 2^7; % nr of "original nodes" in each direction (for Q1P0 defines the 2x2 patches)
nY = nX; 

x = linspace(0,1,nX/2 + 1).^2; % more nodes at boundary
x = [x - 1, 1 - x(end-1:-1:1)];

% x = linspace(-1,1,nX);
y = x; 

% make sure box is 1x1 (size scales the effective Re!)
x = x/2;
y = y/2;

feMesh = createRectMesh(x, y, basisOrder);
nrNodes = feMesh.problemSize(3)*feMesh.problemSize(4);
nrElts = feMesh.problemSize(1)*feMesh.problemSize(2);

fprintf('%7.0d elements and %7.0d nodes \n',[nrElts,nrNodes])

% assemble matrices
% M = massAssembly( feMesh, localMatrix.mass); % mass matrix
L = massPAssembly( feMesh, localMatrix.pressure); % "pressure mass" matrix
% D = diffusionAssembly( feMesh, localMatrix.stiff); % diffusive matrix
D = laplaceAssembly( feMesh, localMatrix.stiff); % alternative (not using Sij)
globalMatrix = struct('D', D, 'L', L);

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
% solVec(lidNodes) = lidVel*(1 - feMesh.node(1, lidNodes).^6); % regularised
solVec(lidNodes) = lidVel;

% solve corresponding stokes problem as initial guess
M = [D, -L'; -L -feMesh.stabC];
rhsVec = -M(:, fixedVel)*solVec(fixedVel);

solVec(freeSol) = matrixSolve(M(freeSol, freeSol), rhsVec(freeSol),1);

stokesSol = solVec;

% solve nonlinear equation
% solVec = newtonSolve(feMesh, globalMatrix, localMatrix, solVec, typesNodes,...
% 	Re, tol, maxIt);

% plot stuff
solArray = reshape(solVec(1:nrNodes),feMesh.problemSize(4),feMesh.problemSize(3));
centerSol = solArray(:,round(feMesh.problemSize(3)/2));

stokesArray = reshape(stokesSol(1:nrNodes),feMesh.problemSize(4),feMesh.problemSize(3));
centerStokes = stokesArray(:,round(feMesh.problemSize(3)/2));

% plotSol(feMesh, solVec, ['Re = ',num2str(Re)]);

yPoints = feMesh.node(2,1:feMesh.problemSize(4));

figure 
plot(2*yPoints,centerSol/lidVel,'r')
hold on
plot(2*yPoints,centerStokes/lidVel,'k')

% reference sols
plot(2*refSolYVal - 1, refSol400, 'rx' )
plot(comsolNSRe1000(:,1),  comsolNSRe1000(:,2), 'r+' )

plot(comsolCavityNSRe1(:,1), comsolCavityNSRe1(:,2), 'kx')

title(['Centerline velocity, Re = ',num2str(Re)])
legend('NS','Stokes', 'Ghia 400 NS', 'Comsol NS 1000', 'Stokes Comsol')

hold off
% grid on