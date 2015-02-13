function [x, relres, resminres] = matrixSolve(A, L, C, f, setup, typeA, Q)
% 
% solves system of linear equations of the form
%
% [A, -L'; -L, -C] [u, p] = [f, 0]
%
% where C may be zero (if no stabilization is needed). Moreover when solving
% the stokes problem, A is symmetric. 
%
% Q is the pressure mass matrix which is used as a (block-) preconditioner

[nrp, nrv] = size(L);
M = [A, -L'; -L -C];

% check input
if nargin < 7
	Q = [];
end


if strcmp(typeA, 'stokes')
	if setup.linsolve.stokesPrecon == 2	
		% build preconditioner
		icholsetup.type = 'ict'; 
		icholsetup.droptol = 1e-4;

		P1 = ichol(A, icholsetup);
		P2 = sqrt(spdiags(diag(Q),0,nrp,nrp));
		precon = [P1 sparse(nrv, nrp); sparse(nrp, nrv) P2];

	elseif setup.linsolve.stokesPrecon == 1
		% AMG preconditioner
		[precon, setup] = createAmgSystem(A, setup); % can be done more efficiently!
		precon.nrv = nrv;
		precon.nrp = nrp;
	end

	[x, ~, relres, ~, resminres] = myMinres(M, f, zeros(size(f)), setup, precon, Q);

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