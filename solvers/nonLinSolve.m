function [solVec] = nonLinSolve(feMesh, globalMatrix, localMatrix, solVec ,...
	typesNodes, Re, setup)
% hybrid method using both picard and newton (alternating)

tol = setup.nonlin.tol;
maxIt = setup.nonlin.maxIt;

setup.linsolve.solver = 0; % backslash (option is changed only inside this function)
% setup.linsolve.precon = 1; % amg

nrNodes = feMesh.problemSize(3)*feMesh.problemSize(4);
nrElts = feMesh.problemSize(1)*feMesh.problemSize(2);
nrPBasisF = size(localMatrix.pdivv.x,1);

freeSol = typesNodes.freeSol;
freePressure = typesNodes.freePressure;
freeVel = typesNodes.freeVel;
fixedVel = typesNodes.fixedVel;
fixedPressure = typesNodes.fixedPressure;

stepsD = 2; % nr of diverging picard steps before breaking
stepsC = 3; % nr of converging picard steps before switching to newton 
initPicard = 3;
errorEst = zeros(maxIt + 1,1);
deltaSol = 0;
iter = 0;
method = 'picard'; % 1 = newton, 0 = picard (start)
convList = ones(maxIt + 1, 1); % if 0, then errotEst(iter) < errorEst(iter + 1)

% start with picard 
% if 1 step divergence > switch to picard
% if picard diverges for stepsD steps then break
% if picard converges for stepsC steps, then switch to newton again

colWidth = [4, 10, 10, 8];
printTableSep(colWidth); % print horizontal separator
printTable({'i', 'residual', 'update', 'method'},...
 	repmat({'string'}, 1, 4), colWidth);
printTableSep(colWidth);


while iter <= maxIt
	nonLin2 = nonLinearAssembly2(feMesh, localMatrix.nonlin, solVec);

	% construct matrix
	if strcmp(method, 'newton')
		% newton
		nonLin1 = nonLinearAssembly1(feMesh, localMatrix.nonlin, solVec);
		velMatrix = nonLin1 + nonLin2 + 1/Re*globalMatrix.D;
	else
		% picard
		velMatrix =  nonLin2 + 1/Re*globalMatrix.D;
	end
	% A = [velMatrix, -globalMatrix.L'; -globalMatrix.L, ...
	% 	-globalMatrix.stabC];

	% calculate residuals
	momentumRes = nonLin2*solVec(1:2*nrNodes) -...
		globalMatrix.L'*solVec(2*nrNodes + 1:end) +...
		1/Re*globalMatrix.D*solVec(1:2*nrNodes);

	contRes = - globalMatrix.L*solVec(1:2*nrNodes) -...
		globalMatrix.stabC*solVec(2*nrNodes + 1:end);

	% the rhs
	rhsVec = [-momentumRes; contRes];

	% error of current solution
	errorEst(iter + 1) = norm(rhsVec(freeSol));
	if (iter > 0) & (errorEst(iter) < errorEst(iter + 1))
		convList(iter + 1) = 0;
		% method is diverging (newton)
		method = 'picard';
	end

	if strcmp(method, 'picard') & (all(convList(...
		max(1, iter - (stepsC - 2)):iter + 1))) & (iter > initPicard - 1)
		% picard is converging stepsC steps, so switch to newton
		method = 'newton';
	end

	% print convergence info
	printTable({iter, errorEst(iter + 1), norm(deltaSol), method},...
	 	{'%d', '%4.2e', '%4.2e', 'string'},colWidth);

	% check if tol is reached or stepsD steps of divergence
	if (errorEst(iter + 1) < tol) |...
		(all(~convList(max(1, iter - (stepsD - 1)):iter + 1)))
		break
	end

	% while statement is useless .. but this is easier
	if iter <= maxIt
		
		% solve the updates
		deltaSol = matrixSolve(velMatrix(freeVel, freeVel),...
			globalMatrix.L(freePressure, freeVel),... 
			globalMatrix.stabC(freePressure, freePressure),...
			rhsVec(freeSol), setup,...
			globalMatrix.Q(freePressure, freePressure),2);

		% update solution
		solVec(freeSol) = solVec(freeSol) + deltaSol;
	end

	iter = iter + 1;
end
printTableSep(colWidth);
end