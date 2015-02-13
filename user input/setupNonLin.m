function [setup] = setupNonLin(setup, setupType)
if nargin < 1
	setup = [];
end
if nargin < 2
	setup.nonlin.tol = default('Tolerance for nonlin iteration', 1e-7); 
	setup.nonlin.maxIt = default('Maximum number of nonlin iterations', 15);
elseif strcmp(setupType, 'default')
	setup.nonlin.tol  = 1e-7;
	setup.nonlin.maxIt = 15;
end
