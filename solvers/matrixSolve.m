function [x] = matrixSolve(A, b, symm)

% solve stokes system (upper left block is symmetric)
if symm
	x = A\b;
	
	
	% using ILUPack
	
	% options = AMGinit(A);

	% tic
	% [PREC,options] = AMGfactor(A, options);
	% toc

	% tic
	% [x, options] = AMGsolver(A, PREC, options, b);
	% toc

	% PREC=AMGdelete(PREC);

else
	% solve system arrising from nonlinear solve (upper left block is unsymmetric)
	x = A\b;

	% save('nonlinsystem', 'A', 'b')  % store for testing solve methods
end
end