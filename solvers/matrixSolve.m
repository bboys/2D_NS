function [x, relres, resminres] = matrixSolve(A, L, C, f, iterPar, typeA, Q)
% 
% solves system of linear equations of the form
%
% [A, -L'; -L, -C] [u, p] = [f, 0]
%
% where C may be zero (if no stabilization is needed). Moreover when solving
% the stokes problem, A is symmetric. 
%
% Q is the pressure mass matrix which is used as a preconditioner

[nrp, nrv] = size(L);
M = [A, -L'; -L -C];

% check input
if nargin < 7
	Q = [];
end
if nargin < 6
	typeA = 'backslash';
end
if (nargin < 5) | (isempty(iterPar) == 1)
	tol = 1e-11;
	maxIt = floor(nrv/10);
else
	tol = iterPar.Ltol;
	maxIt = iterPar.LmaxIt;
end

if strcmp(typeA, 'stokes')

	% build preconditioner
	setup.type = 'ict'; 
	setup.droptol = 1e-4;

	P1 = ichol(A, setup);
	P2 = sqrt(spdiags(diag(Q),0,nrp,nrp));
	P = [P1 sparse(nrv, nrp); sparse(nrp, nrv) P2];

	[x, ~, relres, ~, resminres] = minres(M, f, tol, maxIt, P, P');

elseif strcmp(typeA, 'NS')
	% nonlinear system, hence nonsymmetric matrix A
	x = M\f;
	relres = 0;
	resminres = 0;


else  % strcmp(typeA, 'backslash')

	% nonlinear system, hence nonsymmetric matrix A
	x = M\f;
	relres = 0;
	resminres = 0;
end
end

% using ILUPack

% options = AMGinit(M);

% tic
% [PREC,options] = AMGfactor(M, options);
% toc

% tic
% [x, options] = AMGsolver(M, PREC, options, f);
% toc

% PREC=AMGdelete(PREC);