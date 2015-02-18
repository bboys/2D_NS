function [amgSystem, setup] = createAmgSystem(A, setup)
% creates struct amgSystem containing the matrices 
% A^h, ..., A^(2^levels h), and corresponding interpolation
% operators
%
% amgSystem.level(i).matrix, and
% amgSystem.level(i).interp, such that
%
% I^{2h}_h = amgSystem.level(1).interp'
% amgSystem.level(1).matrix = A^h = A (original matrix, finest grid)
%
% levels indicates the maximum level (1 is sufficient for two cycle AMG)

% initialize
amgSystem.level(1).matrix = A;

% apply same setup on each level
for level = 1:setup.amg.levels
	amgSystem.level(level).interp = createInterpOp(amgSystem.level(level).matrix, setup);
	if size(amgSystem.level(level).interp,2) == 1
		setup.amg.levels = level - 1;
		tempStr = sprintf('Warning: AMG system stopped at level %d.',[level - 1]);
		announce(tempStr)
		break
	end
	amgSystem.level(level + 1).matrix =...
		amgSystem.level(level).interp'*amgSystem.level(level).matrix*...
		amgSystem.level(level).interp;

	% precompute smoothing matrices
	if setup.amg.smoothType == 1 % GS
		if ~isfield(amgSystem.level(level), 'Ls') | isempty(amgSystem.level(level).Ls)
			% make sure this is done only once
			amgSystem.level(level).Ls = tril(amgSystem.level(level).matrix);
			amgSystem.level(level).U = triu(amgSystem.level(level).matrix,1);
		end
	elseif setup.amg.smoothType == 2 % symmetric GS
		if ~isfield(amgSystem.level(level), 'D') | isempty(amgSystem.level(level).D)
			% make sure this is done only once
			n = size(amgSystem.level(level).matrix,1);
			amgSystem.level(level).D = spdiags(diag(...
				amgSystem.level(level).matrix), 0, n,n);
			amgSystem.level(level).U = triu(amgSystem.level(level).matrix,1);
			amgSystem.level(level).L = tril(amgSystem.level(level).matrix,-1);

			amgSystem.level(level).Us = amgSystem.level(level).U + amgSystem.level(level).D;
			amgSystem.level(level).Ls = amgSystem.level(level).L + amgSystem.level(level).D;
		end	
	end
end



end