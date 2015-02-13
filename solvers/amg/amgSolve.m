function [u, relres] = amgSolve(A, b, u, setup, amgSystem)
% perform setup.amg.maxIt number of V-cycles on A*u = b, with initial guess u
% depth of V-cycle is determined by setup.amg.levels
%
% A should be symmetric and positive definite

if nargin < 4
	setup.amg.levels = 4;
	setup.amg.maxIt = 15;

	setup.amg.nrPreSmooth = 1;
	setup.amg.nrPostSmooth = 1;
	setup.amg.smoothType = 'GS'

	setup.amg.coarseMethod = 'RS';
	setup.amg.interpMethod = 'classical';
	setup.amg.theta = 0.8;
end

if nargin < 5
	% allows re-use of interpolation operator
	[amgSystem, setup] = createAmgSystem(A, setup);
end

if nargin < 3
	u = zeros(size(b));
end

if nargout > 1
	relres = zeros(setup.amg.maxIt + 1, 1);
	normb = norm(b);
	relres = norm(b - A*u)/normb;
end

% main loop
for iter = 1:setup.amg.maxIt
	u = amgVCycle(u, b, 1, setup.amg.levels, amgSystem, setup);
	if nargout > 1
		relres(iter + 1) = norm(b - A*u)/normb;
	end
end

end