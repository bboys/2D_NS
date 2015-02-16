warning off

% create local matrices
[ localMatrix, basisOrder, basisType ] = createBasis(2);

% create mesh 
nX = 2^4;
feMesh = createMesh(localMatrix.basisOrder, nX, nX, 2);
globalMatrix.stabC = createStabC(feMesh, localMatrix);

% parameters for linsolve
setup.linsolve.solver = 1; 
setup.linsolve.precon = 1; % 1 = amg, 2 = ichol
setup.linsolve.tol  = 1e-9;
setup.linsolve.maxIt = 200;

% parameters for amg
setup.amg.levels = 5;
setup.amg.maxIt = 2;

setup.amg.nrPreSmooth = 1;
setup.amg.nrPostSmooth = 1;
setup.amg.smoothType = 1;

setup.amg.coarseMethod = 2;
setup.amg.interpMethod = 2;
setup.amg.theta = 0.6;
setup.amg.tol = 1e-11;

% assemble matrices
globalMatrix.L = PdivVAssembly( feMesh, localMatrix.pdivv); % "pressure mass" matrix
globalMatrix.D = laplaceAssembly( feMesh, localMatrix.stiff); % alternative (not using Sij)
globalMatrix.Q = pmassAssembly(feMesh, localMatrix.pmass);

% define bdys (lid driven cavity)
feMesh.boundary(1).type = 1; feMesh.boundary(1).func = [0; 0];
feMesh.boundary(3).type = 1; feMesh.boundary(3).func = [0; 0];
feMesh.boundary(4).type = 1; feMesh.boundary(4).func = [0; 0];
feMesh.boundary(2).type = 1;
feMesh.boundary(2).func = str2func('@(x,y) cavityLidDirichlet(x,y)');

% apply bdy conditions
[nodeType, solVec] = applyBdyCond(feMesh, localMatrix.basisType);

M = [globalMatrix.D, -globalMatrix.L'; -globalMatrix.L -globalMatrix.stabC];
rhsVec = -M(nodeType.freeSol, nodeType.fixedVel)*solVec(nodeType.fixedVel);

preconVals = [2 2 2 2]; % 1 = amg, 2 = ichol
preconName = {'amg', 'ichol', 'none'};
solverVals = [2 2 2 2]; % 1 = minres, 2 = gmres
gmresRestarts = [16 20 24 28];
solverName = {'minres','gmres'};
nrTests = length(preconVals);
timeArray = zeros(nrTests,1);
colors = {'b', 'r', 'k', 'g', 'c'};
legendStr = {};

for test = 1:nrTests
	setup.linsolve.solver = solverVals(test);
	setup.linsolve.precon = preconVals(test);
	setup.linsolve.gmresRestart = gmresRestarts(test);
	tic
	[solVec(nodeType.freeSol), relres, resVecs{test}, precon{test}] =...
		matrixSolve(M(nodeType.freeVel, nodeType.freeVel),...
		globalMatrix.L(nodeType.freePressure, nodeType.freeVel),...
		globalMatrix.stabC(nodeType.freePressure, nodeType.freePressure),...
		rhsVec, setup,...
		globalMatrix.Q(nodeType.freePressure, nodeType.freePressure), 1);
	timeArray(test) = toc;
	legendStr = [legendStr, [preconName{preconVals(test)},', ' ,...
		solverName{solverVals(test)},sprintf(': %4.2f seconds', timeArray(test))]];

end

% figure
% for level = 1:setup.amg.levels + 1
% 	subplot(2,3,level)
% 	spy(precon{2}.level(level).matrix)
% 	title(sprintf('Theta = %3.2f, fill = %4.2f %%', [setup.amg.theta,...
% 		100*nnz(precon{2}.level(level).matrix)/...
% 		prod(size(precon{2}.level(level).matrix))]))
% end


figure
semilogy(resVecs{1}, colors{1})
hold on
for test = 2:nrTests
	semilogy(resVecs{test}/norm(rhsVec), colors{test})
end
hold off
legend(legendStr)
title(['nX = nY = ', num2str(nX),', relative residual'])