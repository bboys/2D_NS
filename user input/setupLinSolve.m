function [setup] = setupLinSolve(setup, setupType)
if nargin < 1
	setup = [];
end
if nargin < 2
	setup.linsolve.stokesPrecon = default({'Preconditioner for Stokes','amg','ichol'},1);
	setup.linsolve.tol = default('Tolerance for linsolve', 1e-11); 
	setup.linsolve.maxIt = default('Maximum number of MVs for linsolve', 500);
elseif strcmp(setupType, 'default')
	setup.linsolve.stokesPrecon = 1;
	setup.linsolve.tol  = 1e-11;
	setup.linsolve.maxIt = 500;
end
end
