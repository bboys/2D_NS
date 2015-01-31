function [solVec] = nonLinSolve(feMesh, globalMatrix, localMatrix, solVec ,...
	typesNodes, Re, tol, maxIt)
% hybrid method using both picard and newton (alternating)


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
errorEst = zeros(maxIt + 1,1);
errorEst(1) = Inf;
iter = 0;
method = 0; % 1 = newton, 0 = picard (start)
convList = ones(maxIt + 1, 1); % if 0, then errotEst(iter) < errorEst(iter + 1)
initPicard = 3;
% start with picard 
% if 1 step divergence > switch to picard
% if picard diverges for stepsD steps then break
% if picard converges for stepsC steps, then switch to newton again

colWidth = [4, 9, 9];
printTableSep(colWidth); % print horizontal separator
printTable({'i', 'residual', 'update'},...
 	repmat({'string'}, 1, 3), colWidth);
printTableSep(colWidth);

while errorEst(iter + 1) > tol
	nonLin2 = nonLinearAssembly2(feMesh, localMatrix.nonlin, solVec);

	% construct matrix
	if method == 1
		% newton
		nonLin1 = nonLinearAssembly1(feMesh, localMatrix.nonlin, solVec);
		velMatrix = nonLin1 + nonLin2 + 1/Re*globalMatrix.D;
	else
		% picard
		velMatrix =  nonLin2 + 1/Re*globalMatrix.D;
	end
	A = [velMatrix, -globalMatrix.L'; -globalMatrix.L, ...
		-feMesh.stabC];

	% calculate residuals
	momentumRes = nonLin2*solVec(1:2*nrNodes) -...
		globalMatrix.L'*solVec(2*nrNodes + 1:end) +...
		1/Re*globalMatrix.D*solVec(1:2*nrNodes);

	contRes = - globalMatrix.L*solVec(1:2*nrNodes) -...
		feMesh.stabC*solVec(2*nrNodes + 1:end);

	% the rhs
	rhsVec = [-momentumRes; contRes];

	% solve the updates
	deltaSol = matrixSolve(A(freeSol, freeSol),...
		rhsVec(freeSol),0);

	% update solution
	solVec(freeSol) = solVec(freeSol) + deltaSol;
	
	% increase counter & check maxit
	iter = iter + 1;

	errorEst(iter + 1) = norm(rhsVec(freeSol));
	if (errorEst(iter) < errorEst(iter + 1))

		convList(iter + 1) = 0;
		% method is diverging (newton)
		method = 0;
	end

	
	if (method == 0) & (all(convList(max(1, iter - (stepsC - 2)):iter + 1))) &...
	 iter > initPicard - 1
		% picard is converging stepsC steps, so switch to newton
		method = 1;
	end

	% check if maxIt is reached or stepsD steps of divergence
	if ((iter + 1 > maxIt) & (errorEst(iter + 1) > tol)) |...
		(all(~convList(max(1, iter - (stepsD - 1)):iter + 1)))
		% fprintf(['Iteration stopped without reaching tolerance of %2.1e', ...
		% 	', stopped after %3.0f iterations with a residual of %2.1e \n'],...
		% 	[tol,iter,errorEst(iter + 1)])
		break
	end
	
	% print convergence info
	% fprintf('%2.0f, residual of equations %4.3e, norm of update %4.3e, %d\n',...
	% 	[iter, norm(rhsVec(freeSol)), norm(deltaSol), method])
	printTable({iter, norm(rhsVec(freeSol)), norm(deltaSol)},...
	 	{'%d', '%4.2e', '%4.2e'},colWidth);
	
end
printTableSep(colWidth);
end