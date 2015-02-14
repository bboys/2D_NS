function [setup] = setupNonLin(setup, setupType)
if nargin < 1
	setup = [];
end
if nargin < 2
	setup.nonlin.tol = default('Tolerance for nonlin iteration', 1e-7); 
	setup.nonlin.maxIt = default('Maximum number of nonlin iterations', 15);
	setup.nonlin.solver = default({'Solver for nonlinsolve', 'backslash', 'gmres'},1);
	setup.nonlin.precon = default({'Preconditioner for nonlinsolve','amg','ilu', 'none'},3);
elseif strcmp(setupType, 'default')
	setup.nonlin.linsolver = 1;
	setup.nonlin.linprecon = 3;
	setup.nonlin.tol  = 1e-7;
	setup.nonlin.maxIt = 15;
end
