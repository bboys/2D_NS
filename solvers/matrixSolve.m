function [x] = matrixSolve(A, b, typeA)
% solves system of linear equations of the form
%
% [A, B*; B, -C] [u, p] = [f, 0]
%
% where C may be zero (if no stabilization is needed). Moreover when solving
% the stokes problem, A is symmetric. 

if strcmp(typeA, 'stokes')
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
	% temporary: use standard solve
	x = A\b;

end
end