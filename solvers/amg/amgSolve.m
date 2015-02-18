function [u, relres] = amgSolve(A, b, u, setup, amgSystem)
% perform setup.amg.maxIt number of V-cycles on A*u = b, with initial guess u
% depth of V-cycle is determined by setup.amg.levels
%
% A should be symmetric and positive definite

if nargin < 4 | isempty(setup)
	% load default values
	setup = setupAMG([], 'default');
end

if nargin < 5
	% allows re-use of interpolation operator & smoothing matrices (A = L + U + D)
	[amgSystem, setup] = createAmgSystem(A, setup);
elseif isempty(A)
	A = amgSystem.level(1).matrix;
end

if nargin < 3 | isempty(u)
	u = zeros(size(b));
end

relres = zeros(setup.amg.maxIt, 1);
normb = norm(b);

% main loop
for iter = 1:setup.amg.maxIt
	[u, res, amgSystem] = amgVCycle(u, b, 1, setup.amg.levels, amgSystem, setup);
	relres(iter) = norm(res)/normb;

	% check convergence
	if (relres(iter) < setup.amg.tol)
		relres(iter + 1) = norm(b - A*u)/normb;
		relres(iter + 2:end) = [];
		break
	end
end


end