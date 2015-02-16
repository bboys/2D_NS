function [setup] = setupLinSolve(setup, setupType)
if nargin < 1
	setup = [];
end
if nargin < 2
	setup.linsolve.solver = default({'Solver for linsolve','minres','gmres'},1);
	if setup.linsolve.solver == 2
		setup.linsolve.gmresRestart = default('Restart parameter for GMRES', 20);
	end
	setup.linsolve.precon = default({'Preconditioner for linsolve','amg','ichol', 'none'},2);
	setup.linsolve.tol = default('Tolerance for linsolve', 1e-11); 
	setup.linsolve.maxIt = default('Maximum number of MVs for linsolve', 500);

elseif strcmp(setupType, 'default')
	setup.linsolve.solver = 1;
	setup.linsolve.precon = 2;
	setup.linsolve.tol  = 1e-11;
	setup.linsolve.maxIt = 500;
	setup.linsolve.gmresRestart = 20;
end
end
