function [u, amgSystem] = amgSmoother(u, level, f, setup, amgSystem)
if setup.amg.smoothType == 1 % GS

	if ~isfield(amgSystem.level(level), 'L') | isempty(amgSystem.level(level).L)
		% make sure this is done only once
		amgSystem.level(level).L = tril(amgSystem.level(level).matrix);
		amgSystem.level(level).U = triu(amgSystem.level(level).matrix,1);
	end

	u = amgSystem.level(level).L\(f - amgSystem.level(level).U*u);

end


end