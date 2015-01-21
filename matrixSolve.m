function [x] = matrixSolve(A, b, symm)
if symm
	% x = A\b;
% 	using ILUPack
	options = AMGinit(A);

	tic
	[PREC,options] = AMGfactor(A, options);
	toc

	tic
	[x, options] = AMGsolver(A, PREC, options, b);
	toc

	PREC=AMGdelete(PREC);
else
	x = A\b;
end
end