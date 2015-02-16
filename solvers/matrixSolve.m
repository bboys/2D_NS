function [x, relres, resvec, precon] = matrixSolve(A, L, C, f, setup, Q, typeSolve)
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

if typeSolve == 2 % non symmetric system coming from nonlinear NS eqns

	if setup.nonlin.precon == 1 % amg
		% AMG preconditioner
		[precon, setup] = createAmgSystem(A, setup); % can be done more efficiently!
		precon.nrv = nrv;
		precon.nrp = nrp;
		precon.Q = Q;
		precon.type = setup.nonlin.precon;
	end

	if setup.nonlin.solver == 1 
		% backslash
		x = M\f;
	elseif setup.nonlin.solver == 2 
		% gmres
		[x, resvec] = gmresPrecon(M, f, zeros(size(f)), setup, precon);
		relres = resvec(end);
	end

elseif typeSolve == 1

	if setup.linsolve.precon == 2 
		% build preconditioner
		icholsetup.type = 'ict'; 
		icholsetup.droptol = 1e-4;

		P1 = ichol(A, icholsetup);
		P2 = sqrt(spdiags(diag(Q),0,nrp,nrp));
		precon.P = [P1 sparse(nrv, nrp); sparse(nrp, nrv) P2];
		precon.Q = Q;
		precon.type = setup.linsolve.precon;
	elseif setup.linsolve.precon == 1
		% AMG preconditioner
		[precon, setup] = createAmgSystem(A, setup); % can be done more efficiently!
		precon.nrv = nrv;
		precon.nrp = nrp;
		precon.Q = Q;
		precon.type = setup.linsolve.precon;
	else
		precon = [];
	end

	

	if setup.linsolve.solver == 1
		% minres
		 [x, ~, relres, ~, resvec] = minresPrecon(M, f, zeros(size(f)), setup,...
		  precon);
	elseif setup.linsolve.solver == 2
		% gmres
		[x, resvec] = gmresPrecon(M, f, zeros(size(f)), setup, precon);
		relres = resvec(end);

	end
end


end

